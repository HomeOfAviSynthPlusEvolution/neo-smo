#!/usr/bin/env python3
"""Side-by-side VapourSynth differential verification for Phase 3: zsmooth.dll vs neo-smo.dll."""

import argparse
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from isolated import run, suppress_crash_ui
suppress_crash_ui()
from fp16_reference import Reference
import numpy as np
import vapoursynth as vs


def make_test_clip(core: vs.Core, fmt: int, width: int = 70, height: int = 42, length: int = 1, seed: int = 42) -> vs.VideoNode:
    base = core.std.BlankClip(format=fmt, width=width, height=height, length=length)
    f_info = base.format

    def populate(n: int, f: vs.VideoFrame) -> vs.VideoFrame:
        rng = np.random.default_rng(seed + n * 17)
        fout = f.copy()
        for p in range(f_info.num_planes):
            pw = fout.width if p == 0 else (fout.width >> f_info.subsampling_w)
            ph = fout.height if p == 0 else (fout.height >> f_info.subsampling_h)
            arr = np.asarray(fout[p])
            yy, xx = np.indices((ph, pw))
            is_chroma = (f_info.color_family == vs.YUV and p > 0)

            if f_info.sample_type == vs.INTEGER:
                peak = (1 << f_info.bits_per_sample) - 1
                noise = rng.integers(0, peak + 1, size=(ph, pw), dtype=np.int64)
                grad = ((xx * 7 + yy * 13 + n * 19) % (peak + 1)).astype(np.int64)
                vals = (noise * 3 + grad) // 4
                vals[0, 0] = peak
                vals[0, pw - 1] = 0
                vals[ph - 1, 0] = 0
                vals[ph - 1, pw - 1] = peak
                vals[ph // 2, pw // 2] = peak
                vals[ph // 3, pw // 3] = 0
                arr[:, :] = vals.astype(arr.dtype)
            else:
                raw = rng.random(size=(ph, pw), dtype=np.float32)
                grad = (((xx * 5 + yy * 11 + n * 7) % 64) / 63.0).astype(np.float32)
                vals = 0.75 * raw + 0.25 * grad
                if is_chroma:
                    vals = vals - np.float32(0.5)
                vals[0, 0] = np.float32(0.5 if is_chroma else 1.0)
                vals[0, pw - 1] = np.float32(-0.5 if is_chroma else 0.0)
                vals[ph // 2, pw // 2] = np.float32(0.5 if is_chroma else 1.0)
                vals[ph // 3, pw // 3] = np.float32(-0.5 if is_chroma else 0.0)
                arr[:, :] = vals.astype(arr.dtype)
        return fout

    return core.std.ModifyFrame(clip=base, clips=base, selector=populate)


def compare_frames(ref_frame: vs.VideoFrame, cand_frame: vs.VideoFrame, fmt: vs.VideoFormat, tol: float):
    max_abs_err = 0.0
    max_region_errs = {"interior": 0.0, "edge": 0.0, "corner": 0.0}
    worst_detail = ""

    for p in range(fmt.num_planes):
        ra = np.asarray(ref_frame[p])
        ca = np.asarray(cand_frame[p])
        if fmt.sample_type == vs.FLOAT and fmt.bits_per_sample == 16:
            ra = ra.view(np.float16).astype(np.float32)
            ca = ca.view(np.float16).astype(np.float32)
        else:
            ra = ra.astype(np.float64)
            ca = ca.astype(np.float64)

        if not np.all(np.isfinite(ra)) or not np.all(np.isfinite(ca)):
            return False, float("inf"), max_region_errs, f"plane={p}: non-finite sample"
        diff = np.abs(ra - ca)
        h, w = diff.shape
        yy, xx = np.indices((h, w))
        ey = (yy < 3) | (yy >= h - 3)
        ex = (xx < 3) | (xx >= w - 3)
        masks = {
            "interior": ~(ey | ex),
            "edge": (ey | ex) & ~(ey & ex),
            "corner": ey & ex,
        }
        for rname, m in masks.items():
            if np.any(m):
                rmax = float(np.max(diff[m]))
                if rmax > max_region_errs[rname]:
                    max_region_errs[rname] = rmax

        p_max = float(np.max(diff))
        if p_max > max_abs_err:
            max_abs_err = p_max
            idx = np.unravel_index(int(np.argmax(diff)), diff.shape)
            worst_detail = f"plane={p} y={idx[0]} x={idx[1]} ref={ra[idx]} cand={ca[idx]} err={p_max}"

    passed = max_abs_err <= tol
    return passed, max_abs_err, max_region_errs, worst_detail


def main():
    parser = argparse.ArgumentParser(description="Phase 3 Differential correctness test against zsmooth.dll")
    parser.add_argument("--plugin", type=Path, required=True, help="Path to neo-smo.dll")
    parser.add_argument("--reference", type=Path, required=True, help="Path to zsmooth.dll")
    parser.add_argument("--worker-format", help=argparse.SUPPRESS)
    args = parser.parse_args()

    if args.worker_format is None:
        failures = 0
        for name in ("YUV420P8", "YUV420P10", "YUV444P16", "YUV444PH", "YUV444PS", "RGB24"):
            print(f"[RUN Phase 3] {name} (isolated process)", flush=True)
            failures += bool(run([sys.executable, "-u", str(Path(__file__).resolve()),
                                  "--plugin", str(args.plugin.resolve()),
                                  "--reference", str(args.reference.resolve()),
                                  "--worker-format", name], timeout=60))
        sys.exit(1 if failures else 0)

    core = vs.core
    core.num_threads = 4
    core.std.LoadPlugin(str(args.reference.resolve()))
    core.std.LoadPlugin(str(args.plugin.resolve()))
    reference = Reference(core, args.plugin)

    formats = [
        ("YUV420P8", vs.YUV420P8, 0.0),
        ("YUV420P10", vs.YUV420P10, 0.0),
        ("YUV444P16", vs.YUV444P16, 0.0),
        ("YUV444PH", vs.YUV444PH, 1e-3),
        ("YUV444PS", vs.YUV444PS, 1e-6),
        ("RGB24", vs.RGB24, 0.0),
    ]

    total_cases = 0
    failed_cases = []

    for fmt_name, fmt_id, float_tol in formats:
        if fmt_name != args.worker_format:
            continue
        clip = make_test_clip(core, fmt_id, width=164, height=42, length=1, seed=300 + int(fmt_id))

        fmt_max_err = 0.0

        # 1. InterQuartileMean: radii 0..3 and per-plane array
        for r in [0, 1, 2, 3, [1, 2, 3]]:
            ref_out = reference.InterQuartileMean(clip, radius=r)
            cand_out = core.neo_smo.InterQuartileMean(clip, radius=r)
            tol = 0.0 if (r == 0 or fmt_id not in (vs.YUV444PH, vs.YUV444PS)) else float_tol
            ok, err, reg, detail = compare_frames(ref_out.get_frame(0), cand_out.get_frame(0), clip.format, tol=tol)
            fmt_max_err = max(fmt_max_err, err)
            total_cases += 1
            if not ok:
                failed_cases.append(f"[FAIL] InterQuartileMean({fmt_name}, radius={r}): {detail} regions={reg}")

        # InterQuartileMean planes subset
        ref_out = reference.InterQuartileMean(clip, radius=2, planes=[0, 2])
        cand_out = core.neo_smo.InterQuartileMean(clip, radius=2, planes=[0, 2])
        tol = 0.0 if fmt_id not in (vs.YUV444PH, vs.YUV444PS) else float_tol
        ok, err, reg, detail = compare_frames(ref_out.get_frame(0), cand_out.get_frame(0), clip.format, tol=tol)
        fmt_max_err = max(fmt_max_err, err)
        total_cases += 1
        if not ok:
            failed_cases.append(f"[FAIL] InterQuartileMean_planes({fmt_name}): {detail}")

        # 2. SmartMedian: default threshold for radii 0..3
        for r in [0, 1, 2, 3, [1, 2, 3]]:
            ref_out = reference.SmartMedian(clip, radius=r)
            cand_out = core.neo_smo.SmartMedian(clip, radius=r)
            tol = 0.0 if (r == 0 or fmt_id not in (vs.YUV444PH, vs.YUV444PS)) else float_tol
            ok, err, reg, detail = compare_frames(ref_out.get_frame(0), cand_out.get_frame(0), clip.format, tol=tol)
            fmt_max_err = max(fmt_max_err, err)
            total_cases += 1
            if not ok:
                ok = reference.smart_median_threshold_match(clip, r, cand_out.get_frame(0))
            if not ok:
                failed_cases.append(f"[FAIL] SmartMedian({fmt_name}, radius={r}): {detail} regions={reg}")

        # SmartMedian: explicit threshold & scalep
        for r in [1, 2]:
            for scalep in [False, True]:
                th = 40.0 if scalep else (0.15 if fmt_id in (vs.YUV444PH, vs.YUV444PS) else 40.0)
                ref_out = reference.SmartMedian(clip, radius=r, threshold=th, scalep=scalep)
                cand_out = core.neo_smo.SmartMedian(clip, radius=r, threshold=th, scalep=scalep)
                tol = 0.0 if fmt_id not in (vs.YUV444PH, vs.YUV444PS) else float_tol
                ok, err, reg, detail = compare_frames(ref_out.get_frame(0), cand_out.get_frame(0), clip.format, tol=tol)
                fmt_max_err = max(fmt_max_err, err)
                total_cases += 1
                if not ok:
                    failed_cases.append(f"[FAIL] SmartMedian_th({fmt_name}, r={r}, th={th}, scalep={scalep}): {detail}")

        # SmartMedian: per-plane threshold array & planes subset
        th_arr = [20.0, 50.0, 80.0] if clip.format.sample_type == vs.INTEGER else [0.1, 0.2, 0.3]
        ref_out = reference.SmartMedian(clip, radius=2, threshold=th_arr, planes=[0, 1])
        cand_out = core.neo_smo.SmartMedian(clip, radius=2, threshold=th_arr, planes=[0, 1])
        tol = 0.0 if fmt_id not in (vs.YUV444PH, vs.YUV444PS) else float_tol
        ok, err, reg, detail = compare_frames(ref_out.get_frame(0), cand_out.get_frame(0), clip.format, tol=tol)
        fmt_max_err = max(fmt_max_err, err)
        total_cases += 1
        if not ok:
            failed_cases.append(f"[FAIL] SmartMedian_planes({fmt_name}): {detail}")

        # Small width stress check
        small_clip = make_test_clip(core, fmt_id, width=16, height=16, length=1, seed=999 + int(fmt_id))
        _ = core.neo_smo.InterQuartileMean(small_clip, radius=3).get_frame(0)
        _ = core.neo_smo.SmartMedian(small_clip, radius=3).get_frame(0)

        print(f"[{'FAIL' if failed_cases else 'PASS'}] {fmt_name:<10} | Phase 3 cases verified | max_abs_err = {fmt_max_err:.3e}")


    # Threshold conversion and round-half-up must not change the filter decision.
    # Repeat the witness across SIMD interiors, tails and spatial borders.
    if args.worker_format == "YUV420P8":
        witnesses = [
            (vs.GRAY8, 17.75, [0, 0, 0, 0, 255, 2, 2, 2, 2], 255),
            (vs.GRAY16, 17.75, [0, 0, 0, 0, 65535, 2, 2, 2, 2], 65535),
            (vs.GRAY16, 2704.0, [724, 724, 724, 996, 65535, 1208, 1208, 1208, 1208], 65535),
            (vs.GRAYH, 0.9188, [0, 0, 0, 0, 1, 0.1, 0.1, 0.1, 0.1], float(np.float16(0.1)) if reference.native else 1.0),
        ]
        for fmt_id, threshold, patch, expected in witnesses:
            base = core.std.BlankClip(format=fmt_id, width=164, height=42)
            dtype = np.float16 if fmt_id == vs.GRAYH else np.uint8 if fmt_id == vs.GRAY8 else np.uint16
            data = np.tile(np.asarray(patch, dtype=dtype).reshape(3, 3), (14, 55))[:, :164]
            def fill_witness(n, f, data=data):
                out = f.copy()
                np.asarray(out[0])[:] = data
                return out
            clip = core.std.ModifyFrame(base, clips=base, selector=fill_witness)
            ref = reference.SmartMedian(clip, radius=1, threshold=threshold).get_frame(0)
            cand = core.neo_smo.SmartMedian(clip, radius=1, threshold=threshold).get_frame(0)
            a, b = np.asarray(ref[0]), np.asarray(cand[0])
            total_cases += 1
            if not np.array_equal(a, b) or float(b[7, 7]) != expected:
                failed_cases.append(f"[FAIL] SmartMedian threshold witness fmt={fmt_id} th={threshold}: "
                                    f"reference={a[7, 7]} candidate={b[7, 7]} expected={expected}")

    print(f"Phase 3 Differential Verification Summary: {total_cases - len(failed_cases)}/{total_cases} cases passed.")
    if failed_cases:
        for msg in failed_cases:
            print(msg, file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
