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
