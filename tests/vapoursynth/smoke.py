"""Bounded real-host coverage of all 19 filters, without an upstream plugin."""
import argparse
from pathlib import Path

import numpy as np
import vapoursynth as vs


def pixels(frame, plane):
    data = np.asarray(frame[plane])
    if frame.format.sample_type == vs.FLOAT and frame.format.bits_per_sample == 16:
        data = data.view(np.float16)
    return data


def source(core, format_id):
    base = core.std.BlankClip(width=96, height=64, length=7, format=format_id)

    def fill(n, f):
        out = f.copy()
        rng = np.random.default_rng(731 + n)
        for p in range(out.format.num_planes):
            a = pixels(out, p)
            if out.format.sample_type == vs.INTEGER:
                a[:] = rng.integers(0, 1 << out.format.bits_per_sample, a.shape)
            else:
                a[:] = rng.random(a.shape) - (0.5 if out.format.color_family == vs.YUV and p else 0)
        out.props['SmokeFrame'] = n
        out.props['_SceneChangePrev'] = int(n == 3)
        out.props['_SceneChangeNext'] = int(n == 2)
        return out

    return core.std.ModifyFrame(base, clips=base, selector=fill)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plugin', required=True)
    args = parser.parse_args()
    core = vs.core
    core.num_threads = 2
    core.std.LoadPlugin(path=str(Path(args.plugin).resolve()))
    exercised = set()
    frames = 0
    for format_id in (vs.GRAY8, vs.YUV420P8, vs.YUV420P16, vs.YUV420PH,
                      vs.YUV420PS, vs.RGB24, vs.RGB48, vs.RGBH, vs.RGBS):
        clip = source(core, format_id)
        ref = core.std.FlipHorizontal(clip)
        fmt = clip.format
        cases = [
            ('Median', dict(radius=[1])),
            ('InterQuartileMean', dict(radius=[2])),
            ('SmartMedian', dict(radius=[1], threshold=[50], scalep=True)),
            ('VerticalCleaner', dict(mode=[1])),
            ('RemoveGrain', dict(mode=[2])),
            ('Repair', dict(repairclip=ref, mode=[1])),
            ('Clense', dict(previous=ref, next=ref)),
            ('ForwardClense', {}), ('BackwardClense', {}),
            ('TemporalMedian', dict(radius=2, scenechange=True)),
            ('TemporalSoften', dict(radius=2, scenechange=-1)),
            ('TemporalRepair', dict(repairclip=ref, mode=[0])),
            ('DegrainMedian', dict(limit=[4], scalep=True)),
            ('FluxSmoothT', {}), ('FluxSmoothST', {}),
            ('DCTFilter', dict(factors=[1, .9, .7, .5, .3, .2, .1, 0], planes=[0])),
        ]
        if not (fmt.sample_type == vs.FLOAT and fmt.bits_per_sample == 16):
            cases.append(('TTempSmooth', dict(pfclip=ref, scthresh=-1)))
        if fmt.color_family != vs.GRAY:
            cases.append(('CCD', dict(ref=ref, scale=1.0, temporal_radius=1)))
        if fmt.color_family == vs.YUV and fmt.sample_type == vs.INTEGER:
            cases.append(('Cnr4', dict(ref=ref, scenechange=True)))
        for name, params in cases:
            print(f'{fmt.name}: {name}', flush=True)
            out = getattr(core.neo_smo, name)(clip, **params)
            assert (out.width, out.height, out.num_frames, out.format.id) == (96, 64, 7, fmt.id)
            for n in (0, 3, 6):
                with out.get_frame(n) as frame, clip.get_frame(n) as original:
                    assert frame.props['SmokeFrame'] == n
                    for p in range(fmt.num_planes):
                        data = pixels(frame, p)
                        assert np.all(np.isfinite(data)), (fmt.name, name, n, p)
                        if fmt.sample_type == vs.INTEGER:
                            assert data.max() < (1 << fmt.bits_per_sample)
                        if (name == 'DCTFilter' and p > 0) or (name in ('CCD', 'Cnr4') and fmt.color_family == vs.YUV and p == 0):
                            np.testing.assert_array_equal(data, pixels(original, p))
                frames += 1
            exercised.add(name)
    assert len(exercised) == 19, exercised
    print(f'Passed: {len(exercised)} functions, 9 formats, {frames} frames', flush=True)


if __name__ == '__main__':
    main()
