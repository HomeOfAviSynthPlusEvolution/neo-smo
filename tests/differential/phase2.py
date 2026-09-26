#!/usr/bin/env python3
"""Side-by-side VapourSynth differential verification for Phase 2: zsmooth.dll vs neo-smo.dll."""

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
    parser = argparse.ArgumentParser(description="Phase 2 Differential correctness test against zsmooth.dll")
    parser.add_argument("--plugin", type=Path, required=True, help="Path to neo-smo.dll")
    parser.add_argument("--reference", type=Path, required=True, help="Path to zsmooth.dll")
    parser.add_argument("--worker-format", help=argparse.SUPPRESS)
    args = parser.parse_args()

    if args.worker_format is None:
        failures = 0
        for name in ("YUV420P8", "YUV420P10", "YUV444P16", "YUV444PH", "YUV444PS", "RGB24"):
            print(f"[RUN Phase 2] {name} (isolated process)", flush=True)
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
        clip = make_test_clip(core, fmt_id, width=164, height=42, length=7, seed=200 + int(fmt_id))
        repair_clip = make_test_clip(core, fmt_id, width=164, height=42, length=7, seed=500 + int(fmt_id))

        fmt_max_err = 0.0

        # 1. Repair: all modes 0..24 + multi-plane mode array
        for m in list(range(0, 25)) + [[2, 10, 18], [1, 24, 0]]:
            ref_out = reference.Repair(clip, repair_clip, mode=m)
            cand_out = core.neo_smo.Repair(clip, repair_clip, mode=m)
            # Modes 0, 1, 2, 3, 4, 11, 12, 13, 14, 17 are pure sorting/min/max; others involve diffs/clamping
            tol = 0.0 if (m in (0, 1, 2, 3, 4, 11, 12, 13, 14, 17)) else float_tol
            ok, err, reg, detail = compare_frames(ref_out.get_frame(0), cand_out.get_frame(0), clip.format, tol=tol)
            fmt_max_err = max(fmt_max_err, err)
            total_cases += 1
            if not ok:
                failed_cases.append(f"[FAIL] Repair({fmt_name}, mode={m}): {detail} regions={reg}")

        # 2. Clense: standard (prev/next = clip) across multiple frames (frame 0, 1, 3, 6)
        ref_out = reference.Clense(clip)
        cand_out = core.neo_smo.Clense(clip)
        for frame_idx in (0, 1, 3, 6):
            ok, err, reg, detail = compare_frames(ref_out.get_frame(frame_idx), cand_out.get_frame(frame_idx), clip.format, tol=0.0)
            fmt_max_err = max(fmt_max_err, err)
            total_cases += 1
            if not ok:
                failed_cases.append(f"[FAIL] Clense({fmt_name}, frame={frame_idx}): {detail} regions={reg}")

        # Clense with explicit previous & next clips
        prev_clip = make_test_clip(core, fmt_id, width=164, height=42, length=7, seed=700 + int(fmt_id))
        next_clip = make_test_clip(core, fmt_id, width=164, height=42, length=7, seed=900 + int(fmt_id))
        ref_out = reference.Clense(clip, previous=prev_clip, next=next_clip, planes=[0, 1])
        cand_out = core.neo_smo.Clense(clip, previous=prev_clip, next=next_clip, planes=[0, 1])
        for frame_idx in (0, 2, 6):
            ok, err, reg, detail = compare_frames(ref_out.get_frame(frame_idx), cand_out.get_frame(frame_idx), clip.format, tol=0.0)
            fmt_max_err = max(fmt_max_err, err)
            total_cases += 1
            if not ok:
                failed_cases.append(f"[FAIL] Clense_explicit({fmt_name}, frame={frame_idx}): {detail} regions={reg}")

        # 3. ForwardClense: across frames (0, 1, 4, 5, 6)
        ref_out = reference.ForwardClense(clip)
        cand_out = core.neo_smo.ForwardClense(clip)
        for frame_idx in (0, 1, 4, 5, 6):
            ok, err, reg, detail = compare_frames(ref_out.get_frame(frame_idx), cand_out.get_frame(frame_idx), clip.format, tol=float_tol)
            fmt_max_err = max(fmt_max_err, err)
            total_cases += 1
            if not ok:
                failed_cases.append(f"[FAIL] ForwardClense({fmt_name}, frame={frame_idx}): {detail} regions={reg}")

        # ForwardClense with planes subset
        ref_out = reference.ForwardClense(clip, planes=[0])
        cand_out = core.neo_smo.ForwardClense(clip, planes=[0])
        ok, err, reg, detail = compare_frames(ref_out.get_frame(2), cand_out.get_frame(2), clip.format, tol=float_tol)
        fmt_max_err = max(fmt_max_err, err)
        total_cases += 1
        if not ok:
            failed_cases.append(f"[FAIL] ForwardClense_planes({fmt_name}): {detail} regions={reg}")

        # 4. BackwardClense: across frames (0, 1, 2, 5, 6)
        ref_out = reference.BackwardClense(clip)
        cand_out = core.neo_smo.BackwardClense(clip)
        for frame_idx in (0, 1, 2, 5, 6):
            ok, err, reg, detail = compare_frames(ref_out.get_frame(frame_idx), cand_out.get_frame(frame_idx), clip.format, tol=float_tol)
            fmt_max_err = max(fmt_max_err, err)
            total_cases += 1
            if not ok:
                failed_cases.append(f"[FAIL] BackwardClense({fmt_name}, frame={frame_idx}): {detail} regions={reg}")

        # BackwardClense with planes subset
        ref_out = reference.BackwardClense(clip, planes=[1])
        cand_out = core.neo_smo.BackwardClense(clip, planes=[1])
        ok, err, reg, detail = compare_frames(ref_out.get_frame(3), cand_out.get_frame(3), clip.format, tol=float_tol)
        fmt_max_err = max(fmt_max_err, err)
        total_cases += 1
        if not ok:
            failed_cases.append(f"[FAIL] BackwardClense_planes({fmt_name}): {detail} regions={reg}")

        # Small width stress check
        small_clip = make_test_clip(core, fmt_id, width=16, height=16, length=5, seed=999 + int(fmt_id))
        small_rep = make_test_clip(core, fmt_id, width=16, height=16, length=5, seed=888 + int(fmt_id))
        _ = core.neo_smo.Repair(small_clip, small_rep, mode=20).get_frame(0)
        _ = core.neo_smo.Clense(small_clip).get_frame(2)
        _ = core.neo_smo.ForwardClense(small_clip).get_frame(1)
        _ = core.neo_smo.BackwardClense(small_clip).get_frame(3)

        # VSHelper's format comparison deliberately ignores FPS and length.
        # Short references repeat their last frame; request frames out of order.
        for length in (1, 3, 9):
            other = make_test_clip(core, fmt_id, width=164, height=42, length=length, seed=1234)
            other = core.std.AssumeFPS(other, fpsnum=25, fpsden=1)
            for name, kwargs in [('Repair', {'repairclip': other, 'mode': 20}),
                                 ('Clense', {'previous': other, 'next': other}),
                                 ('Clense', {'previous': clip, 'next': other})]:
                ref_out = getattr(reference, name)(clip, **kwargs)
                cand_out = getattr(core.neo_smo, name)(clip, **kwargs)
                for frame_idx in (6, 2, 0, 5, 1, 3):
                    ok, err, reg, detail = compare_frames(ref_out.get_frame(frame_idx),
                        cand_out.get_frame(frame_idx), clip.format, tol=float_tol)
                    total_cases += 1
                    fmt_max_err = max(fmt_max_err, err)
                    if not ok:
                        failed_cases.append(f'[FAIL] {name} length={length} frame={frame_idx}: {detail}')

        # Very short temporal inputs must return boundary frames unchanged.
        for length in (1, 2):
            short = clip[:length]
            for name in ('Clense', 'ForwardClense', 'BackwardClense'):
                out = getattr(core.neo_smo, name)(short)
                for n in range(length):
                    ok, _, _, detail = compare_frames(short.get_frame(n), out.get_frame(n), short.format, 0.0)
                    total_cases += 1
                    if not ok:
                        failed_cases.append(f'[FAIL] {name} short boundary: {detail}')

        print(f"[{'FAIL' if failed_cases else 'PASS'}] {fmt_name:<10} | Phase 2 cases verified | max_abs_err = {fmt_max_err:.3e}")

    print(f"Phase 2 Differential Verification Summary: {total_cases - len(failed_cases)}/{total_cases} cases passed.")
    if failed_cases:
        for msg in failed_cases:
            print(msg, file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
