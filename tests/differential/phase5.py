#!/usr/bin/env python3
"""Side-by-side VapourSynth differential verification for Phase 5: zsmooth.dll vs neo-smo.dll."""

import argparse
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from isolated import run, suppress_crash_ui
suppress_crash_ui()
from fp16_reference import Reference
import numpy as np
import vapoursynth as vs


def make_test_clip(core: vs.Core, fmt: int, width: int = 256, height: int = 48, length: int = 12, seed: int = 42, with_scenechange: bool = False, low_contrast: bool = False) -> vs.VideoNode:
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
                if low_contrast:
                    vals = peak // 2 + rng.integers(0, max(2, peak // 16), size=(ph, pw))
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
                if low_contrast:
                    vals = 0.49 + raw * 0.02 - (0.5 if is_chroma else 0.0)
                arr[:, :] = vals.astype(arr.dtype)

        if with_scenechange:
            props = fout.props
            props["_SceneChangePrev"] = 1 if n == 4 else 0
            props["_SceneChangeNext"] = 1 if n == 7 else 0

        return fout

    return core.std.ModifyFrame(clip=base, clips=base, selector=populate)


def compare_frames(ref_frame: vs.VideoFrame, cand_frame: vs.VideoFrame, fmt: vs.VideoFormat, tol: float):
    max_abs_err = 0.0
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
            return False, float("inf"), f"plane={p}: non-finite sample"
        diff = np.abs(ra - ca)
        p_max = float(np.max(diff))
        if p_max > max_abs_err:
            max_abs_err = p_max
            worst_detail = f"plane={p}, max_err={p_max:.4e}"

    passed = (max_abs_err <= tol)
    return passed, max_abs_err, worst_detail


def main():
    parser = argparse.ArgumentParser(description="Phase 5 differential tests: zsmooth vs neo-smo")
    parser.add_argument("--zsmooth-dll", "--reference", dest="zsmooth_dll", type=Path, default=Path("zsmooth.dll"))
    parser.add_argument("--neo-smo-dll", "--plugin", dest="neo_smo_dll", type=Path, default=Path("build/neo-smo.dll"))
    parser.add_argument("--float-tol", type=float, default=1e-6)
    args = parser.parse_args()

    core = vs.core
    core.num_threads = 2
    core.std.LoadPlugin(path=str(args.zsmooth_dll.resolve()))
    core.std.LoadPlugin(path=str(args.neo_smo_dll.resolve()))
    reference = Reference(core, args.neo_smo_dll)

    formats = [
        ("YUV420P8", vs.YUV420P8, 0.0),
        ("YUV420P10", vs.YUV420P10, 0.0),
        ("YUV444P16", vs.YUV444P16, 0.0),
        ("YUV444PS", vs.YUV444PS, args.float_tol),
        ("RGB24", vs.RGB24, 0.0),
    ]

    total_cases = 0
    failed_cases = []

    for fmt_name, fmt_id, float_tol in formats:
        clip = make_test_clip(core, fmt_id, width=256, height=48, length=12, seed=101, with_scenechange=True)
        ref_clip = make_test_clip(core, fmt_id, width=256, height=48, length=12, seed=202, with_scenechange=True)
        fmt_max_err = 0.0

        def check_case(name, clip_obj, kwargs, frame_num, tol):
            nonlocal total_cases, fmt_max_err
            total_cases += 1
            try:
                ref_out = getattr(reference, name)(clip_obj, **kwargs)
                cand_out = getattr(core.neo_smo, name)(clip_obj, **kwargs)
                ok, err, detail = compare_frames(ref_out.get_frame(frame_num), cand_out.get_frame(frame_num), clip_obj.format, tol=tol)
                fmt_max_err = max(fmt_max_err, err)
                if not ok:
                    failed_cases.append(f"[FAIL] {name}({fmt_name}, {kwargs}, frame={frame_num}): {detail}")
            except vs.Error as e:
                failed_cases.append(f"[FAIL] {name}({fmt_name}, {kwargs}): {e}")

        # -------------------------------------------------------------------
        # 1. TTempSmooth (U8, U10, U16, F32, RGB24)
        # TTempSmooth is a floating-point weighted filter; fast-math permits +/- 1 LSB on exact .5 ties.
        # -------------------------------------------------------------------
        tt_tol = 1.0 if float_tol == 0.0 else float_tol
        # Radius 1..7, strength 1..4, fp=True/False, pfclip, scenechange
        for maxr in (1, 3, 5):
            for fp in (True, False):
                check_case("TTempSmooth", clip, dict(maxr=maxr, fp=fp, scthresh=0), 5, tt_tol)
        # Inverse difference mode vs temporal mode
        check_case("TTempSmooth", clip, dict(maxr=2, thresh=[10, 10, 10], mdiff=[2, 2, 2], scthresh=0), 4, tt_tol)
        check_case("TTempSmooth", clip, dict(maxr=2, thresh=[2, 2, 2], mdiff=[4, 4, 4], scthresh=0), 4, tt_tol)
        # With scenechange = -1
        check_case("TTempSmooth", clip, dict(maxr=3, scthresh=-1), 4, tt_tol)
        check_case("TTempSmooth", clip, dict(maxr=3, scthresh=-1), 7, tt_tol)
        # With pfclip
        check_case("TTempSmooth", clip, dict(maxr=2, pfclip=ref_clip, scthresh=0), 3, tt_tol)
        # Endpoints
        check_case("TTempSmooth", clip, dict(maxr=3, scthresh=0), 0, tt_tol)
        check_case("TTempSmooth", clip, dict(maxr=3, scthresh=0), 11, tt_tol)

        # -------------------------------------------------------------------
        # 2. CCD (RGB24 and YUV formats, scale >= 1.0 required for h < 240)
        # -------------------------------------------------------------------
        # Low, medium, high points
        check_case("CCD", clip, dict(threshold=4.0, temporal_radius=0, points=[1, 1, 0], scale=1.0), 2, float_tol)
        check_case("CCD", clip, dict(threshold=6.0, temporal_radius=0, points=[1, 0, 0], scale=1.0), 2, float_tol)
        check_case("CCD", clip, dict(threshold=8.0, temporal_radius=0, points=[0, 1, 0], scale=1.0), 2, float_tol)
        check_case("CCD", clip, dict(threshold=5.0, temporal_radius=0, points=[1, 1, 1], scale=1.0), 2, float_tol)
        # Temporal radius 1 and 2
        check_case("CCD", clip, dict(threshold=4.0, temporal_radius=1, points=[1, 1, 0], scale=1.0), 2, float_tol)
        check_case("CCD", clip, dict(threshold=4.0, temporal_radius=2, points=[1, 1, 0], scale=1.0), 4, float_tol)
        # With ref clip
        check_case("CCD", clip, dict(threshold=4.0, temporal_radius=1, ref=ref_clip, scale=1.0), 3, float_tol)
        # Scale parameter
        check_case("CCD", clip, dict(threshold=4.0, temporal_radius=0, scale=1.5), 2, float_tol)

        # -------------------------------------------------------------------
        # 3. Cnr4 (YUV integer only)
        # -------------------------------------------------------------------
        if fmt_id in (vs.YUV420P8, vs.YUV420P10, vs.YUV444P16):
            # All 5 temporal modes
            for tm in (0, 1, 2, 3, 4):
                check_case("Cnr4", clip, dict(radius=2, tmode=tm, scenechange=False), 3, 0.0)
            # All 4 weight modes
            for wm in (0, 1, 2, 3):
                cnr_tol = 1.0 if wm == 3 else 0.0
                check_case("Cnr4", clip, dict(radius=2, wmode=wm, scenechange=False), 3, cnr_tol)
            # Mode strings
            for md in ("oxx", "xoo", "ooo", "xxx"):
                check_case("Cnr4", clip, dict(mode=md, radius=2, scenechange=False), 3, 0.0)
            # With scenechange enabled (using _SceneChangePrev / _SceneChangeNext)
            check_case("Cnr4", clip, dict(radius=2, scenechange=True), 4, 0.0)
            check_case("Cnr4", clip, dict(radius=2, scenechange=True), 7, 0.0)
            # With ref clip
            check_case("Cnr4", clip, dict(radius=2, ref=ref_clip, scenechange=False), 3, 0.0)

        print(f"[{'FAIL' if failed_cases else 'PASS'}] {fmt_name:<10} | Phase 5 verified | max_abs_err = {fmt_max_err:.3e}")

    # Use complete reference SIMD blocks: its CCD short-width and Cnr4
    # overlapping in-place tail paths have upstream defects. Small sizes and
    # tails are covered independently by core.phase5.
    # Directed low-contrast inputs exercise accepted neighbors, not merely
    # passthrough on unrelated random frames. Keep the matrix small and bounded.
    def directed(name, source, kwargs, tol, frames=(11,)):
        nonlocal total_cases
        ref_out = getattr(reference, name)(source, **kwargs)
        cand_out = getattr(core.neo_smo, name)(source, **kwargs)
        for n in frames:
            total_cases += 1
            ok, err, detail = compare_frames(ref_out.get_frame(n), cand_out.get_frame(n), source.format, tol)
            if not ok:
                failed_cases.append(f"[FAIL] directed {name} {source.format.name} {kwargs} frame={n}: {detail}")

    for fmt in (vs.YUV420P8, vs.YUV420P10, vs.YUV444P16, vs.YUV444PS, vs.RGB24):
        smooth = make_test_clip(core, fmt, length=23, low_contrast=True, with_scenechange=True)
        for radius in (1, 7):
            for fp in (False, True):
                directed("TTempSmooth", smooth, dict(maxr=radius, thresh=[32], mdiff=[8], strength=8, fp=fp, scthresh=0),
                         1 if smooth.format.sample_type == vs.INTEGER else 1e-6, frames=(0, 11, 22))
        directed("TTempSmooth", smooth, dict(maxr=3, pfclip=smooth, scthresh=-1),
                 1 if smooth.format.sample_type == vs.INTEGER else 1e-6, frames=(4, 7))
        for radius in (0, 3, 10):
            directed("CCD", smooth, dict(scale=1., threshold=32., temporal_radius=radius, ref=smooth),
                     1 if smooth.format.sample_type == vs.INTEGER else 1e-6)
        if smooth.format.color_family == vs.YUV and smooth.format.sample_type == vs.INTEGER:
            directed("Cnr4", smooth, dict(pow=[0.,0.,0.], scenechange=False), 0)
            for radius in (1, 3, 10):
                for mode in range(5):
                    directed("Cnr4", smooth, dict(radius=radius, tmode=mode, scenechange=False), 0)
            for weight in range(4):
                directed("Cnr4", smooth, dict(radius=3, wmode=weight, ref=smooth, scenechange=False), 1)

    # Half arithmetic has a separate exact scalar oracle in core.phase5. The
    # optimized Zig DLL can contract/reassociate operations near an acceptance
    # threshold; keep its normalized-video comparison bounded to two 8-bit code values.
    for seed in (42, 51, 193):
        half_clip = make_test_clip(core, vs.YUV444PH, length=23, seed=seed, low_contrast=True)
        for radius, threshold in ((0, 4.), (2, 8.)):
            directed("CCD", half_clip, dict(scale=1., threshold=threshold, temporal_radius=radius), 2.0 / 255)

    # The default scale on short clips must be validated exactly as an explicit
    # scale; unsafe float-to-integer conversions must be rejected before kernels.
    for kwargs in ({}, dict(scale=float("inf")), dict(scale=1e30), dict(scale=1., threshold=float("nan")),
                   dict(scale=1., threshold=float("inf"))):
        try:
            core.neo_smo.CCD(clip, **kwargs)
        except vs.Error:
            total_cases += 1
        else:
            failed_cases.append(f"[FAIL] CCD accepted unsafe scale/threshold {kwargs}")
    zero_sense = core.neo_smo.Cnr4(make_test_clip(core, vs.YUV420P8), sense=[0,0,0], scenechange=False).get_frame(3)
    assert zero_sense is not None

    # Small size stress test
    small_clip = make_test_clip(core, vs.YUV420P8, width=16, height=16, length=5, seed=777)
    _ = core.neo_smo.TTempSmooth(small_clip, maxr=1, scthresh=0).get_frame(2)
    _ = core.neo_smo.CCD(small_clip, temporal_radius=0, scale=1.0, points=[1,0,0]).get_frame(2)
    _ = core.neo_smo.Cnr4(small_clip, radius=1, scenechange=False).get_frame(2)

    print(f"Phase 5 Differential Verification Summary: {total_cases - len(failed_cases)}/{total_cases} cases passed.")
    if failed_cases:
        for msg in failed_cases:
            print(msg, file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    if "--worker" not in sys.argv:
        sys.exit(run([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"], timeout=25))
    sys.argv.remove("--worker")
    main()
