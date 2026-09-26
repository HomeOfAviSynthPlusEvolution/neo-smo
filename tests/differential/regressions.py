"""Parameter admission and directed FP16 / tie-heavy Zig comparisons."""
import argparse
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from isolated import run, suppress_crash_ui
suppress_crash_ui()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plugin', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    if not args.worker:
        return bool(run([sys.executable, '-u', str(Path(__file__).resolve()),
                         '--plugin', str(args.plugin.resolve()), '--reference', str(args.reference.resolve()), '--worker']))
    import numpy as np
    import vapoursynth as vs
    core = vs.core
    core.num_threads = 4
    core.std.LoadPlugin(str(args.plugin.resolve()))
    core.std.LoadPlugin(str(args.reference.resolve()))
    base = core.std.BlankClip(format=vs.GRAY8, width=164, height=42)
    checks = 0
    for name, key, maximum in [('Median', 'radius', 3), ('RemoveGrain', 'mode', 24), ('VerticalCleaner', 'mode', 2)]:
        for value in [-2**63, -4294967295, -1, maximum + 1, 2**31, 4294967297, 2**63 - 1]:
            try:
                getattr(core.neo_smo, name)(base, **{key: [value]})
            except vs.Error:
                checks += 1
            else:
                raise AssertionError(f'{name} accepted {value}')
    for value in [-2**63, -4294967295, -1, 25, 2**31, 4294967297, 2**63 - 1]:
        try:
            core.neo_smo.Repair(base, base, mode=[value])
        except vs.Error:
            checks += 1
        else:
            raise AssertionError(f'Repair accepted {value}')
    for planes in [[], [-2**63], [4294967296], [0, 0], [1]]:
        try:
            core.neo_smo.Median(base, planes=planes)
        except vs.Error:
            checks += 1
        else:
            raise AssertionError(f'accepted planes={planes}')
    for name in ['RemoveGrain', 'VerticalCleaner']:
        try:
            getattr(core.neo_smo, name)(base, mode=[])
        except vs.Error:
            checks += 1
        else:
            raise AssertionError(f'{name} accepted empty mode')
    try:
        core.neo_smo.Repair(base, base, mode=[])
    except vs.Error:
        checks += 1
    else:
        raise AssertionError('Repair accepted empty mode')
    for name in ['InterQuartileMean', 'SmartMedian']:
        for value in [-2**63, -4294967295, -1, 4, 2**31, 4294967297, 2**63 - 1]:
            try:
                getattr(core.neo_smo, name)(base, radius=[value])
            except vs.Error:
                checks += 1
            else:
                raise AssertionError(f'{name} accepted {value}')
    for th in [-1.0, 256.0, 1000.0]:
        try:
            core.neo_smo.SmartMedian(base, threshold=[th])
        except vs.Error:
            checks += 1
        else:
            raise AssertionError(f'SmartMedian accepted invalid threshold {th}')
    for r in [-1, 0, 11, 100]:
        try:
            core.neo_smo.TemporalMedian(base, radius=r)
        except vs.Error:
            checks += 1
        else:
            raise AssertionError(f'TemporalMedian accepted radius {r}')
    for sc in [-2, 255, 1000]:
        try:
            core.neo_smo.TemporalSoften(base, scenechange=sc)
        except vs.Error:
            checks += 1
        else:
            raise AssertionError(f'TemporalSoften accepted invalid scenechange {sc}')
    for m in [-1, 5, 10]:
        try:
            core.neo_smo.TemporalRepair(base, base, mode=[m])
        except vs.Error:
            checks += 1
        else:
            raise AssertionError(f'TemporalRepair accepted invalid mode {m}')
    for m in [-1, 6, 10]:
        try:
            core.neo_smo.DegrainMedian(base, mode=[m])
        except vs.Error:
            checks += 1
        else:
            raise AssertionError(f'DegrainMedian accepted invalid mode {m}')
    try:
        core.neo_smo.DegrainMedian(base, limit=[0, 0, 0])
    except vs.Error:
        checks += 1
    else:
        raise AssertionError('DegrainMedian accepted all limits 0')
    for name, parameter in [("TemporalSoften", "threshold"), ("DegrainMedian", "limit"),
                            ("FluxSmoothT", "temporal_threshold"), ("FluxSmoothST", "spatial_threshold")]:
        for value in [float("nan"), float("inf"), -float("inf")]:
            try:
                getattr(core.neo_smo, name)(base, **{parameter: [value]})
            except vs.Error:
                checks += 1
            else:
                raise AssertionError(f"{name} accepted non-finite {parameter}")

    # Phase 5 Parameter Admission Checks
    yuv_base = core.std.BlankClip(format=vs.YUV420P8, width=164, height=42)
    for maxr in [-1, 0, 8, 10]:
        try:
            core.neo_smo.TTempSmooth(yuv_base, maxr=maxr)
        except vs.Error:
            checks += 1
        else:
            raise AssertionError(f"TTempSmooth accepted invalid maxr {maxr}")
    for th in [-1, 0, 257, 1000]:
        try:
            core.neo_smo.TTempSmooth(yuv_base, thresh=[th])
        except vs.Error:
            checks += 1
        else:
            raise AssertionError(f"TTempSmooth accepted invalid thresh {th}")
    for md in [-1, 256, 1000]:
        try:
            core.neo_smo.TTempSmooth(yuv_base, mdiff=[md])
        except vs.Error:
            checks += 1
        else:
            raise AssertionError(f"TTempSmooth accepted invalid mdiff {md}")
    for str_val in [-1, 0, 9, 100]:
        try:
            core.neo_smo.TTempSmooth(yuv_base, strength=str_val)
        except vs.Error:
            checks += 1
        else:
            raise AssertionError(f"TTempSmooth accepted invalid strength {str_val}")
    for sc in [-2.0, 101.0, 500.0]:
        try:
            core.neo_smo.TTempSmooth(yuv_base, scthresh=sc)
        except vs.Error:
            checks += 1
        else:
            raise AssertionError(f"TTempSmooth accepted invalid scthresh {sc}")

    # CCD
    for tr in [-1, 11, 20]:
        try:
            core.neo_smo.CCD(yuv_base, temporal_radius=tr)
        except vs.Error:
            checks += 1
        else:
            raise AssertionError(f"CCD accepted invalid temporal_radius {tr}")
    for sc in [-1.0, 0.0, 0.5, 0.99]:
        try:
            core.neo_smo.CCD(yuv_base, scale=sc)
        except vs.Error:
            checks += 1
        else:
            raise AssertionError(f"CCD accepted invalid scale {sc}")
    for pts in [[], [1], [1, 1], [1, 1, 1, 1], [0, 0, 0]]:
        try:
            core.neo_smo.CCD(yuv_base, points=pts, scale=1.0)
        except vs.Error:
            checks += 1
        else:
            raise AssertionError(f"CCD accepted invalid points {pts}")

    # Cnr4
    for r in [-1, 0, 11, 20]:
        try:
            core.neo_smo.Cnr4(yuv_base, radius=r)
        except vs.Error:
            checks += 1
        else:
            raise AssertionError(f"Cnr4 accepted invalid radius {r}")
    for tm in [-1, 5, 10]:
        try:
            core.neo_smo.Cnr4(yuv_base, tmode=tm)
        except vs.Error:
            checks += 1
        else:
            raise AssertionError(f"Cnr4 accepted invalid tmode {tm}")
    for wm in [-1, 4, 10]:
        try:
            core.neo_smo.Cnr4(yuv_base, wmode=wm)
        except vs.Error:
            checks += 1
        else:
            raise AssertionError(f"Cnr4 accepted invalid wmode {wm}")
    for md in ["", "o", "ox", "oxxx", "abc", "oyo"]:
        try:
            core.neo_smo.Cnr4(yuv_base, mode=md)
        except vs.Error:
            checks += 1
        else:
            raise AssertionError(f"Cnr4 accepted invalid mode {md}")
    for s_val in [-2, 256, 1000]:
        try:
            core.neo_smo.Cnr4(yuv_base, sense=[s_val, 47, 47])
        except vs.Error:
            checks += 1
        else:
            raise AssertionError(f"Cnr4 accepted invalid sense {s_val}")
    for st_val in [-2, 256, 1000]:
        try:
            core.neo_smo.Cnr4(yuv_base, str=[st_val, 255, 255])
        except vs.Error:
            checks += 1
        else:
            raise AssertionError(f"Cnr4 accepted invalid str {st_val}")
    for p_val in [-1.0, -0.5]:
        try:
            core.neo_smo.Cnr4(yuv_base, pow=[p_val, 1.0, 1.0])
        except vs.Error:
            checks += 1
        else:
            raise AssertionError(f"Cnr4 accepted invalid pow {p_val}")

    print(f'{checks} invalid parameter cases rejected')

    count = 0
    max_float_error = 0.0
    for fmt, kind in [(vs.GRAYH, 'subnormal'), (vs.GRAYH, 'signed_subnormal'),
                      (vs.GRAYH, 'normal_boundary'), (vs.GRAY8, 'ties'),
                      (vs.GRAY16, 'ties'), (vs.GRAYS, 'ties')]:
        for seed in [44, 81, 193]:
            base = core.std.BlankClip(format=fmt, width=164, height=42)
            rng = np.random.default_rng(seed)
            if kind == 'subnormal':
                data = rng.integers(0, 1024, (42, 164), dtype=np.uint16).view(np.float16)
            elif kind == 'signed_subnormal':
                data = (rng.integers(0, 1024, (42, 164), dtype=np.uint16) |
                        (rng.integers(0, 2, (42, 164), dtype=np.uint16) << 15)).view(np.float16)
            elif kind == 'normal_boundary':
                data = rng.integers(1000, 1048, (42, 164), dtype=np.uint16).view(np.float16)
            else:
                data = rng.integers(0, 4, (42, 164))
                if fmt == vs.GRAYS:
                    data = data.astype(np.float32) / 4
            def fill(n, f):
                out = f.copy()
                np.asarray(out[0])[:] = data
                return out
            clip = core.std.ModifyFrame(base, clips=base, selector=fill)
            for name, key, modes in [('Median', 'radius', range(4)), ('RemoveGrain', 'mode', range(25)), ('VerticalCleaner', 'mode', range(3))]:
                for mode in modes:
                    ref = getattr(core.zsmooth, name)(clip, **{key: mode}).get_frame(0)
                    cand = getattr(core.neo_smo, name)(clip, **{key: mode}).get_frame(0)
                    a, b = np.asarray(ref[0]), np.asarray(cand[0])
                    if clip.format.sample_type == vs.FLOAT:
                        aa, bb = a.astype(np.float64), b.astype(np.float64)
                        # Normalized-video tolerance, not bitwise equivalence.
                        atol, rtol = (2**-20, 2**-10) if fmt == vs.GRAYH else (1e-6, 1e-6)
                        error = np.abs(aa - bb)
                        ok = np.isfinite(aa).all() and np.isfinite(bb).all() and np.all(error <= atol + rtol * np.abs(aa))
                        max_float_error = max(max_float_error, float(error.max()))
                    else:
                        ok = np.array_equal(a, b)
                    if not ok:
                        raise AssertionError(f'{clip.format.name} {kind} seed={seed} {name}({mode}): {np.count_nonzero(a != b)} differing samples')
                    count += 1

            rep_base = core.std.BlankClip(format=fmt, width=164, height=42)
            def fill_rep(n, f):
                out = f.copy()
                # Preserve subnormals and dense ties in both inputs.
                np.asarray(out[0])[:] = (data + 1) if clip.format.sample_type != vs.FLOAT else np.roll(data, (5, 17), (0, 1))
                return out
            rep_clip = core.std.ModifyFrame(rep_base, clips=rep_base, selector=fill_rep)
            for mode in range(25):
                ref = core.zsmooth.Repair(clip, rep_clip, mode=mode).get_frame(0)
                cand = core.neo_smo.Repair(clip, rep_clip, mode=mode).get_frame(0)
                a, b = np.asarray(ref[0]), np.asarray(cand[0])
                if clip.format.sample_type == vs.FLOAT:
                    aa, bb = a.astype(np.float64), b.astype(np.float64)
                    atol, rtol = (2**-20, 2**-10) if fmt == vs.GRAYH else (1e-6, 1e-6)
                    error = np.abs(aa - bb)
                    ok = np.isfinite(aa).all() and np.isfinite(bb).all() and np.all(error <= atol + rtol * np.abs(aa))
                    max_float_error = max(max_float_error, float(error.max()))
                else:
                    ok = np.array_equal(a, b)
                if not ok:
                    raise AssertionError(f'{clip.format.name} {kind} seed={seed} Repair({mode}): {np.count_nonzero(a != b)} differing samples')
                count += 1
    print(f'{count} directed differential cases passed (integer exact, float bounded); max_float_abs={max_float_error}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
