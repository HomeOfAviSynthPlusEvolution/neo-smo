#!/usr/bin/env python3
"""Side-by-side VapourSynth differential verification: zsmooth.dll vs neo-smo.dll."""

import argparse
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from isolated import run, suppress_crash_ui
suppress_crash_ui()
import numpy as np
import vapoursynth as vs


def make_test_clip(core: vs.Core, fmt: int, width: int = 70, height: int = 42, seed: int = 42) -> vs.VideoNode:
    base = core.std.BlankClip(format=fmt, width=width, height=height, length=1)
    f_info = base.format

    def populate(n: int, f: vs.VideoFrame) -> vs.VideoFrame:
        rng = np.random.default_rng(seed + n)
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
                grad = ((xx * 7 + yy * 13) % (peak + 1)).astype(np.int64)
                vals = (noise * 3 + grad) // 4
                # Inject impulses and extreme values at edges and interior
                vals[0, 0] = peak
                vals[0, pw - 1] = 0
                vals[ph - 1, 0] = 0
                vals[ph - 1, pw - 1] = peak
                vals[ph // 2, pw // 2] = peak
                vals[ph // 3, pw // 3] = 0
                arr[:, :] = vals.astype(arr.dtype)
            else:
                raw = rng.random(size=(ph, pw), dtype=np.float32)
                grad = (((xx * 5 + yy * 11) % 64) / 63.0).astype(np.float32)
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
    parser = argparse.ArgumentParser(description="Differential correctness test against zsmooth.dll")
    parser.add_argument("--plugin", type=Path, required=True, help="Path to neo-smo.dll")
    parser.add_argument("--reference", type=Path, required=True, help="Path to zsmooth.dll")
    parser.add_argument("--worker-format", help=argparse.SUPPRESS)
    args = parser.parse_args()

    if args.worker_format is None:
        failures = 0
        for name in ("YUV420P8", "YUV420P10", "YUV444P16", "YUV444PH", "YUV444PS", "RGB24"):
            print(f"[RUN] {name} (isolated process)", flush=True)
            failures += bool(run([sys.executable, "-u", str(Path(__file__).resolve()),
                                  "--plugin", str(args.plugin.resolve()),
                                  "--reference", str(args.reference.resolve()),
                                  "--worker-format", name], timeout=45))
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
        clip = make_test_clip(core, fmt_id, width=164, height=42, seed=100 + int(fmt_id))

        fmt_max_err = 0.0

        # 1. Median: radii 0..3 and per-plane combinations
        for r in [0, 1, 2, 3, [1, 2, 3]]:
            ref_out = core.zsmooth.Median(clip, radius=r)
            cand_out = core.neo_smo.Median(clip, radius=r)
            # Median is a pure order-statistic filter: tolerance is strictly 0.0 even on FP16/FP32
            ok, err, reg, detail = compare_frames(ref_out.get_frame(0), cand_out.get_frame(0), clip.format, tol=0.0)
            fmt_max_err = max(fmt_max_err, err)
            total_cases += 1
            if not ok:
                failed_cases.append(f"[FAIL] Median({fmt_name}, radius={r}): {detail} regions={reg}")

        # Median with planes subset
        ref_out = core.zsmooth.Median(clip, radius=2, planes=[0, 2])
        cand_out = core.neo_smo.Median(clip, radius=2, planes=[0, 2])
        ok, err, reg, detail = compare_frames(ref_out.get_frame(0), cand_out.get_frame(0), clip.format, tol=0.0)
        fmt_max_err = max(fmt_max_err, err)
        total_cases += 1
        if not ok:
            failed_cases.append(f"[FAIL] Median({fmt_name}, radius=2, planes=[0,2]): {detail}")

        # 2. VerticalCleaner: modes 0, 1, 2, and [1, 2, 0]
        for m in [0, 1, 2, [1, 2, 0]]:
            ref_out = core.zsmooth.VerticalCleaner(clip, mode=m)
            cand_out = core.neo_smo.VerticalCleaner(clip, mode=m)
            tol = 0.0 if m in (0, 1) else float_tol
            ok, err, reg, detail = compare_frames(ref_out.get_frame(0), cand_out.get_frame(0), clip.format, tol=tol)
            fmt_max_err = max(fmt_max_err, err)
            total_cases += 1
            if not ok:
                failed_cases.append(f"[FAIL] VerticalCleaner({fmt_name}, mode={m}): {detail} regions={reg}")

        # 3. RemoveGrain: all modes 0..24 + multi-plane mode array
        for m in list(range(0, 25)) + [[4, 12, 20]]:
            ref_out = core.zsmooth.RemoveGrain(clip, mode=m)
            cand_out = core.neo_smo.RemoveGrain(clip, mode=m)
            tol = 0.0 if (m in (0, 1, 2, 3, 4, 17)) else float_tol
            ok, err, reg, detail = compare_frames(ref_out.get_frame(0), cand_out.get_frame(0), clip.format, tol=tol)
            fmt_max_err = max(fmt_max_err, err)
            total_cases += 1
            if not ok:
                failed_cases.append(f"[FAIL] RemoveGrain({fmt_name}, mode={m}): {detail} regions={reg}")

        # 4. Small-width stress check for neo_smo (width < vector_len + 2*radius)
        small_clip = make_test_clip(core, fmt_id, width=16, height=16, seed=999 + int(fmt_id))
        _ = core.neo_smo.Median(small_clip, radius=3).get_frame(0)
        _ = core.neo_smo.VerticalCleaner(small_clip, mode=2).get_frame(0)
        _ = core.neo_smo.RemoveGrain(small_clip, mode=20).get_frame(0)

        print(f"[{'FAIL' if failed_cases else 'PASS'}] {fmt_name:<10} | 36 filter/mode cases verified | max_abs_err = {fmt_max_err:.3e}")

    print(f"Differential Verification Summary: {total_cases - len(failed_cases)}/{total_cases} cases passed.")
    if failed_cases:
        for msg in failed_cases:
            print(msg, file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
