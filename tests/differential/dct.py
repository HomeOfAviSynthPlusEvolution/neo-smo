#!/usr/bin/env python3
"""Bounded DCTFilter API, edge extension and upstream comparison tests."""
import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from isolated import run, suppress_crash_ui

suppress_crash_ui()
import numpy as np
import vapoursynth as vs


def array(frame, plane):
    a = np.asarray(frame[plane])
    if frame.format.sample_type == vs.FLOAT and frame.format.bits_per_sample == 16:
        a = a.view(np.float16)
    return a


def make_clip(core, fmt, width, height):
    base = core.std.BlankClip(format=fmt, width=width, height=height, length=4)

    def fill(n, f):
        out = f.copy()
        rng = np.random.default_rng(731 + n)
        for p in range(out.format.num_planes):
            a = array(out, p)
            if out.format.sample_type == vs.INTEGER:
                a[:] = rng.integers(0, 1 << out.format.bits_per_sample, a.shape)
            else:
                a[:] = rng.random(a.shape) - (0.5 if out.format.color_family == vs.YUV and p else 0)
        out.props["DctTest"] = n + 71
        return out

    return core.std.ModifyFrame(base, clips=base, selector=fill)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--plugin", required=True)
    parser.add_argument("--reference", required=True)
    args = parser.parse_args()
    core = vs.core
    core.num_threads = 2
    core.std.LoadPlugin(path=str(Path(args.plugin).resolve()))
    core.std.LoadPlugin(path=str(Path(args.reference).resolve()))
    cases = 0
    maxima = {}
    for family in (vs.GRAY, vs.RGB, vs.YUV):
        for sample, bits in ((vs.INTEGER, 8), (vs.INTEGER, 10), (vs.INTEGER, 16), (vs.FLOAT, 16), (vs.FLOAT, 32)):
            fmt = core.query_video_format(family, sample, bits, 1 if family == vs.YUV else 0,
                                          1 if family == vs.YUV else 0)
            sizes = ((2, 2), (18, 22), (144, 32)) if family == vs.YUV else ((1, 1), (19, 11), (144, 32))
            for width, height in sizes:
                clip = make_clip(core, fmt.id, width, height)
                for factors, planes in (([1] * 8, None), ([1, .9, .7, .5, .3, .2, .1, 0], [0]),
                                        ([0] * 8, None), ([1, 0, 0, 0, 0, 0, 0, 0],
                                                         [1] if fmt.num_planes == 3 else None)):
                    kw = {} if planes is None else {"planes": planes}
                    candidate = core.neo_smo.DCTFilter(clip, factors=factors, **kw)
                    reference = core.zsmooth.DCTFilter(clip, factors=factors, **kw)
                    assert candidate.width == width and candidate.height == height
                    futures = [candidate.get_frame_async(n) for n in range(4)]
                    for n, future in enumerate(futures):
                        got, want, src = future.result(), reference.get_frame(n), clip.get_frame(n)
                        assert got.props["DctTest"] == n + 71
                        for p in range(fmt.num_planes):
                            a, b = array(got, p), array(want, p)
                            if planes is not None and p not in planes:
                                assert np.array_equal(a, array(src, p)), "unselected plane changed"
                            error = float(np.max(np.abs(a.astype(np.float64) - b.astype(np.float64))))
                            tolerance = 1 if sample == vs.INTEGER else (0.001 if bits == 16 else 2e-6)
                            assert np.all(np.isfinite(a)) and error <= tolerance, (
                                fmt.name, width, height, factors, planes, n, p, error)
                            maxima[fmt.name] = max(maxima.get(fmt.name, 0), error)
                    cases += 1
    clip = make_clip(core, vs.GRAY8, 8, 8)
    for factors in ([], [1] * 7, [1] * 9, [float("nan")] * 8, [float("inf")] * 8, [-.01] * 8, [1.01] * 8):
        try:
            core.neo_smo.DCTFilter(clip, factors=factors)
        except vs.Error:
            pass
        else:
            raise AssertionError("accepted invalid factors: " + str(factors))
    for planes in ([], [0, 0], [-1], [1]):
        try:
            core.neo_smo.DCTFilter(clip, factors=[1] * 8, planes=planes)
        except vs.Error:
            pass
        else:
            raise AssertionError("accepted invalid planes: " + str(planes))
    print(f"DCTFilter: {cases} cases x 4 frames passed; invalid parameters rejected")
    print("Maximum absolute difference by format:", maxima)


if __name__ == "__main__":
    if "--worker" not in sys.argv:
        sys.exit(run([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"], timeout=45))
    sys.argv.remove("--worker")
    main()
