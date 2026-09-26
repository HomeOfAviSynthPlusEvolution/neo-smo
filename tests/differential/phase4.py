#!/usr/bin/env python3
"""Side-by-side VapourSynth differential verification for Phase 4: zsmooth.dll vs neo-smo.dll."""

import argparse
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from isolated import run, suppress_crash_ui
suppress_crash_ui()
import numpy as np
import vapoursynth as vs


def make_test_clip(core: vs.Core, fmt: int, width: int = 70, height: int = 42, length: int = 12, seed: int = 42, with_scenechange: bool = False) -> vs.VideoNode:
    base = core.std.BlankClip(format=fmt, width=width, height=height, length=length)
    f_info = base.format

    def populate(n: int, f: vs.VideoFrame) -> vs.VideoFrame:
        rng = np.random.default_rng(seed + n * 31)
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

        if with_scenechange:
            # Set scene change flags at frame 3 and frame 8
            props = fout.props
            props["_SceneChangePrev"] = 1 if n == 4 else 0
            props["_SceneChangeNext"] = 1 if n == 7 else 0

        return fout

    return core.std.ModifyFrame(clip=base, clips=base, selector=populate)


def compare_frames(ref_frame: vs.VideoFrame, cand_frame: vs.VideoFrame, fmt: vs.VideoFormat, tol: float):
    max_abs_err = 0.0
    max_region_errs = {"interior": 0.0, "edge": 0.0, "corner": 0.0}
    worst_detail = ""
    total_violations = 0

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
        violations = int(np.count_nonzero(diff > tol))
        total_violations += violations
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
    parser = argparse.ArgumentParser(description="Phase 4 Differential correctness test against zsmooth.dll")
    parser.add_argument("--plugin", type=Path, required=True, help="Path to neo-smo.dll")
    parser.add_argument("--reference", type=Path, required=True, help="Path to zsmooth.dll")
    parser.add_argument("--worker-format", help=argparse.SUPPRESS)
    args = parser.parse_args()

    if args.worker_format is None:
        failures = 0
        for name in ("YUV420P8", "YUV420P10", "YUV444P16", "YUV444PH", "YUV444PS", "RGB24"):
            print(f"[RUN Phase 4] {name} (isolated process)", flush=True)
            failures += bool(run([sys.executable, "-u", str(Path(__file__).resolve()),
                                  "--plugin", str(args.plugin.resolve()),
                                  "--reference", str(args.reference.resolve()),
                                  "--worker-format", name], timeout=90))
        sys.exit(1 if failures else 0)

    core = vs.core
    core.num_threads = 4
    core.std.LoadPlugin(str(args.reference.resolve()))
    core.std.LoadPlugin(str(args.plugin.resolve()))

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
        clip = make_test_clip(core, fmt_id, width=164, height=42, length=12, seed=400 + int(fmt_id), with_scenechange=False)
        clip_sc = make_test_clip(core, fmt_id, width=164, height=42, length=12, seed=400 + int(fmt_id), with_scenechange=True)
        rep_clip = make_test_clip(core, fmt_id, width=164, height=42, length=12, seed=600 + int(fmt_id), with_scenechange=False)

        fmt_max_err = 0.0

        # 1. TemporalMedian: radius 1, 2 (without & with scenechange)
        for r in [1, 2]:
            ref_out = core.zsmooth.TemporalMedian(clip, radius=r)
            cand_out = core.neo_smo.TemporalMedian(clip, radius=r)
            for fn in [0, r, r + 2, 11 - r, 11]:
                ok, err, reg, detail = compare_frames(ref_out.get_frame(fn), cand_out.get_frame(fn), clip.format, tol=0.0)
                fmt_max_err = max(fmt_max_err, err)
                total_cases += 1
                if not ok:
                    failed_cases.append(f"[FAIL] TemporalMedian({fmt_name}, r={r}, f={fn}): {detail}")

        # TemporalMedian with scenechange
        ref_out = core.zsmooth.TemporalMedian(clip_sc, radius=2, scenechange=True)
        cand_out = core.neo_smo.TemporalMedian(clip_sc, radius=2, scenechange=True)
        for fn in [2, 4, 7, 9]:
            ok, err, reg, detail = compare_frames(ref_out.get_frame(fn), cand_out.get_frame(fn), clip.format, tol=0.0)
            fmt_max_err = max(fmt_max_err, err)
            total_cases += 1
            if not ok:
                failed_cases.append(f"[FAIL] TemporalMedian_sc({fmt_name}, f={fn}): {detail}")

        # 2. TemporalSoften: radius 1, 2 (without & with scenechange)
        for r in [1, 2]:
            ref_out = core.zsmooth.TemporalSoften(clip, radius=r)
            cand_out = core.neo_smo.TemporalSoften(clip, radius=r)
            tol = 0.0 if fmt_id not in (vs.YUV444PH, vs.YUV444PS) else float_tol
            for fn in [r, r + 3]:
                ok, err, reg, detail = compare_frames(ref_out.get_frame(fn), cand_out.get_frame(fn), clip.format, tol=tol)
                fmt_max_err = max(fmt_max_err, err)
                total_cases += 1
                if not ok:
                    failed_cases.append(f"[FAIL] TemporalSoften({fmt_name}, r={r}, f={fn}): {detail}")

        # TemporalSoften with scenechange
        ref_out = core.zsmooth.TemporalSoften(clip_sc, radius=2, scenechange=-1)
        cand_out = core.neo_smo.TemporalSoften(clip_sc, radius=2, scenechange=-1)
        tol = 0.0 if fmt_id not in (vs.YUV444PH, vs.YUV444PS) else float_tol
        for fn in [3, 5, 8]:
            ok, err, reg, detail = compare_frames(ref_out.get_frame(fn), cand_out.get_frame(fn), clip.format, tol=tol)
            fmt_max_err = max(fmt_max_err, err)
            total_cases += 1
            if not ok:
                failed_cases.append(f"[FAIL] TemporalSoften_sc({fmt_name}, f={fn}): {detail}")

        # 3. TemporalRepair: modes 0..4
        for m in [0, 1, 2, 3, 4, [2, 2, 2]]:
            ref_out = core.zsmooth.TemporalRepair(clip, rep_clip, mode=m)
            cand_out = core.neo_smo.TemporalRepair(clip, rep_clip, mode=m)
            tol = 0.0 if fmt_id not in (vs.YUV444PH, vs.YUV444PS) else float_tol
            for fn in [0, 1, 5, 11]:
                ok, err, reg, detail = compare_frames(ref_out.get_frame(fn), cand_out.get_frame(fn), clip.format, tol=tol)
                fmt_max_err = max(fmt_max_err, err)
                total_cases += 1
                if not ok:
                    failed_cases.append(f"[FAIL] TemporalRepair({fmt_name}, m={m}, f={fn}): {detail}")

        # 4. DegrainMedian: modes 0..5, norow, interlaced
        for m in [0, 1, 2, 3, 4, 5]:
            ref_out = core.zsmooth.DegrainMedian(clip, mode=m)
            cand_out = core.neo_smo.DegrainMedian(clip, mode=m)
            tol = 0.0 if fmt_id not in (vs.YUV444PH, vs.YUV444PS) else float_tol
            ok, err, reg, detail = compare_frames(ref_out.get_frame(2), cand_out.get_frame(2), clip.format, tol=tol)
            fmt_max_err = max(fmt_max_err, err)
            total_cases += 1
            if not ok:
                failed_cases.append(f"[FAIL] DegrainMedian({fmt_name}, m={m}): {detail}")

        # DegrainMedian with norow & interlaced
        ref_out = core.zsmooth.DegrainMedian(clip, mode=1, norow=True, interlaced=True)
        cand_out = core.neo_smo.DegrainMedian(clip, mode=1, norow=True, interlaced=True)
        tol = 0.0 if fmt_id not in (vs.YUV444PH, vs.YUV444PS) else float_tol
        ok, err, reg, detail = compare_frames(ref_out.get_frame(3), cand_out.get_frame(3), clip.format, tol=tol)
        fmt_max_err = max(fmt_max_err, err)
        total_cases += 1
        if not ok:
            failed_cases.append(f"[FAIL] DegrainMedian_interlaced_norow({fmt_name}): {detail}")

        # 5. FluxSmoothT
        ref_out = core.zsmooth.FluxSmoothT(clip)
        cand_out = core.neo_smo.FluxSmoothT(clip)
        tol = 0.0 if fmt_id not in (vs.YUV444PH, vs.YUV444PS) else float_tol
        for fn in [0, 1, 4, 11]:
            ok, err, reg, detail = compare_frames(ref_out.get_frame(fn), cand_out.get_frame(fn), clip.format, tol=tol)
            fmt_max_err = max(fmt_max_err, err)
            total_cases += 1
            if not ok:
                failed_cases.append(f"[FAIL] FluxSmoothT({fmt_name}, f={fn}): {detail}")

        # 6. FluxSmoothST
        ref_out = core.zsmooth.FluxSmoothST(clip)
        cand_out = core.neo_smo.FluxSmoothST(clip)
        tol = 0.0 if fmt_id not in (vs.YUV444PH, vs.YUV444PS) else float_tol
        for fn in [0, 1, 4, 11]:
            ok, err, reg, detail = compare_frames(ref_out.get_frame(fn), cand_out.get_frame(fn), clip.format, tol=tol)
            fmt_max_err = max(fmt_max_err, err)
            total_cases += 1
            if not ok:
                failed_cases.append(f"[FAIL] FluxSmoothST({fmt_name}, f={fn}): {detail}")

        # Directed parameters, temporal endpoints, larger windows and optional properties.
        def check_case(name, source, kwargs, frames):
            nonlocal total_cases
            try:
                ref = getattr(core.zsmooth, name)(source, **kwargs)
                cand = getattr(core.neo_smo, name)(source, **kwargs)
                for frame in frames:
                    total_cases += 1
                    ok, err, _, detail = compare_frames(ref.get_frame(frame), cand.get_frame(frame), source.format, float_tol)
                    if not ok:
                        failed_cases.append(f"[FAIL] {name} {fmt_name} {kwargs} frame={frame}: {detail}")
            except vs.Error as error:
                total_cases += 1
                failed_cases.append(f"[FAIL] {name} {fmt_name} {kwargs}: {error}")

        long_clip = make_test_clip(core, fmt_id, width=164, height=10, length=23, seed=97)
        for radius in range(3, 11):
            check_case("TemporalMedian", long_clip, dict(radius=radius), [radius, 11, 22-radius])
        for radius in (1, 4, 10):
            for sc in (0, -1):
                check_case("TemporalSoften", clip, dict(radius=radius, scenechange=sc), [0, 1, 5, 10, 11])
        # A low-contrast source ensures thresholds admit samples, unlike full-range noise.
        def low_contrast(n, f):
            out = f.copy()
            for plane in range(out.format.num_planes):
                a = np.asarray(out[plane])
                yy, xx = np.indices(a.shape)
                jitter = (xx * 3 + yy * 5 + n * 7) % 17
                values = 80 + jitter if out.format.sample_type == vs.INTEGER else 0.1 + jitter / 1024
                a[:] = values.astype(a.dtype)
            return out
        smooth = core.std.ModifyFrame(clip=clip, clips=clip, selector=low_contrast)
        native = 3.75 if clip.format.sample_type == vs.INTEGER else 0.00375
        for scalep, value in ((False, native), (True, 3.75)):
            check_case("TemporalSoften", smooth, dict(radius=4, threshold=[value], scalep=scalep), [0, 4, 11])
            for mode in range(6):
                for interlaced, norow in ((False, False), (True, True)):
                    check_case("DegrainMedian", smooth, dict(mode=mode, limit=[value], scalep=scalep,
                               interlaced=interlaced, norow=norow), [3])
            check_case("FluxSmoothT", smooth, dict(temporal_threshold=[value], scalep=scalep), [3])
            for temporal, spatial in ((value, value), (-1, value), (value, -1), (-1, -1)):
                check_case("FluxSmoothST", smooth, dict(temporal_threshold=[temporal], spatial_threshold=[spatial], scalep=scalep), [3])
        check_case("FluxSmoothT", smooth, dict(temporal_threshold=[-1]), [3])
        for name, kwargs in (
            ("TemporalSoften", dict(threshold=[native, 0, native], planes=[0, 2])),
            ("DegrainMedian", dict(mode=[0, 3, 5], limit=[native, 0, native])),
            ("FluxSmoothT", dict(temporal_threshold=[native, -1, native], planes=[0, 2])),
            ("FluxSmoothST", dict(temporal_threshold=[native, -1, -1], spatial_threshold=[-1, native, -1])),
        ):
            check_case(name, smooth, kwargs, [3])
        gray = core.std.ShufflePlanes(smooth, planes=0, colorfamily=vs.GRAY)
        check_case("TemporalSoften", gray, dict(threshold=[native, native, native]), [0, 3, 11])
        # The repair clip may be shorter; VS repeats its last available frame.
        for mode in range(5):
            check_case("TemporalRepair", clip, dict(repairclip=rep_clip[:2], mode=mode), [1, 5, 10])
            check_case("TemporalRepair", clip, dict(repairclip=clip, mode=mode), [1, 5])

        # Upstream's array reader repeats mode[0]. Validate intended per-plane
        # behavior against separate scalar-mode calls instead of copying that bug.
        mixed = core.neo_smo.TemporalRepair(clip, rep_clip, mode=[0, 2, 4]).get_frame(5)
        for plane, mode in enumerate((0, 2, 4)):
            expected = core.zsmooth.TemporalRepair(clip, rep_clip, mode=mode).get_frame(5)
            aa = np.asarray(expected[plane]).astype(np.float64)
            bb = np.asarray(mixed[plane]).astype(np.float64)
            total_cases += 1
            if not np.all(np.abs(aa - bb) <= float_tol):
                failed_cases.append(f"[FAIL] TemporalRepair per-plane mode {fmt_name} plane={plane}")

        for length in (1, 2, 3):
            check_case("TemporalSoften", smooth[:length], dict(radius=10), range(length))

        # Small width stress check
        small_clip = make_test_clip(core, fmt_id, width=16, height=16, length=5, seed=999 + int(fmt_id))
        small_rep = make_test_clip(core, fmt_id, width=16, height=16, length=5, seed=888 + int(fmt_id))
        _ = core.neo_smo.TemporalMedian(small_clip, radius=1).get_frame(2)
        _ = core.neo_smo.TemporalSoften(small_clip, radius=1).get_frame(2)
        _ = core.neo_smo.TemporalRepair(small_clip, small_rep, mode=1).get_frame(2)
        _ = core.neo_smo.DegrainMedian(small_clip, mode=1).get_frame(2)
        _ = core.neo_smo.FluxSmoothT(small_clip).get_frame(2)
        _ = core.neo_smo.FluxSmoothST(small_clip).get_frame(2)

        print(f"[{'FAIL' if failed_cases else 'PASS'}] {fmt_name:<10} | Phase 4 cases verified | max_abs_err = {fmt_max_err:.3e}")

    print(f"Phase 4 Differential Verification Summary: {total_cases - len(failed_cases)}/{total_cases} cases passed.")
    if failed_cases:
        for msg in failed_cases:
            print(msg, file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
