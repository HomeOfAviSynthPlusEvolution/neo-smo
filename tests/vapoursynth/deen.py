"""Small real-host Deen checks against independent mathematical definitions."""
import argparse
import math
from pathlib import Path
import numpy as np
import vapoursynth as vs


def pixels(f, p):
    a = np.asarray(f[p])
    return a.view(np.float16) if f.format.sample_type == vs.FLOAT and f.format.bits_per_sample == 16 else a


def expected(planes, mode, threshold, minimum=0.5):
    center = planes[0].astype(np.float64)
    num = np.zeros_like(center)
    den = np.zeros_like(center)
    for k, plane in enumerate(planes):
        padded = np.pad(plane.astype(np.float64), 1, mode='edge')
        for dy in range(-1, 2):
            for dx in range(-1, 2):
                s = padded[1+dy:1+dy+center.shape[0], 1+dx:1+dx+center.shape[1]]
                g = 1-(1-minimum)*math.sqrt((dx*dx+dy*dy)/2)
                mask = np.abs(s-center) <= threshold*(g if mode[0] == 'a' else 1)
                weight = g*(2 if k == 0 and len(planes) == 3 else 1) if mode[0] == 'w' else 1
                num += np.where(mask, s, 0 if mode[0] == 'a' else center)*weight
                den += mask if mode[0] == 'a' else weight
    return num/den


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plugin', required=True)
    parser.add_argument('--mode', required=True)
    args = parser.parse_args()
    core = vs.core
    core.num_threads = 2
    core.std.LoadPlugin(path=str(Path(args.plugin).resolve()))
    mode = args.mode
    for fmt in (vs.GRAY8, vs.YUV420P16, vs.YUV420PH, vs.YUV420PS, vs.RGBS):
        base = core.std.BlankClip(width=8, height=6, length=3, format=fmt)
        def fill(n, f):
            out = f.copy()
            for p in range(out.format.num_planes):
                a = pixels(out, p)
                grid = np.arange(a.size).reshape(a.shape)
                v = ((grid*7+n*11+p*19) % 100)/255
                if out.format.sample_type == vs.INTEGER:
                    v = np.floor(v*((1 << out.format.bits_per_sample)-1))
                elif p and out.format.color_family == vs.YUV:
                    v = v-0.5
                a[:] = v
            out.props['DeenFrame'] = n
            return out
        clip = core.std.ModifyFrame(base, clips=base, selector=fill)
        kwargs = dict(mode=mode, thrY=100, thrUV=100, tthY=100, tthUV=100, scenechange=False)
        out = core.neo_smo.Deen(clip, **kwargs)
        # Exercise asynchronous scheduling, reverse order, endpoints, and repeat requests.
        requests = [(n, out.get_frame_async(n)) for n in (2, 0, 1, 1)]
        for n, future in requests:
            f = future.result()
            assert f.props['DeenFrame'] == n
            indices = (n, n-1, n+1) if mode[1] == '3' and n == 1 else (n,)
            frames = [clip.get_frame(i) for i in indices]
            peak = (1 << f.format.bits_per_sample)-1 if f.format.sample_type == vs.INTEGER else 1
            for p in range(f.format.num_planes):
                ref = expected([pixels(src, p) for src in frames], mode, 100*peak/255)
                if f.format.sample_type == vs.INTEGER:
                    ref = np.floor(ref+0.5)
                    tolerance = 1 if mode[0] == 'w' else 0
                else:
                    tolerance = 2e-6
                    if f.format.bits_per_sample == 16:
                        q = np.abs(ref.astype(np.float16).astype(np.float64))
                        tolerance = tolerance+np.exp2(np.maximum(-24, np.floor(np.log2(np.maximum(q, 2**-24)))-10))
                got = pixels(f, p).astype(np.float64)
                assert np.isfinite(got).all()
                assert (np.abs(got-ref) <= tolerance).all(), (mode, fmt, n, p)
        if clip.format.num_planes == 3:
            f = core.neo_smo.Deen(clip, planes=[0], **kwargs).get_frame(1)
            source = clip.get_frame(1)
            for p in (1, 2):
                assert np.array_equal(pixels(f, p), pixels(source, p))
        if mode[1] == '3':
            cut = core.neo_smo.Deen(clip, **dict(kwargs, scenechange=True, scd=0)).get_frame(1)
            spatial = core.neo_smo.Deen(clip, **dict(kwargs, mode=mode[0]+'2d')).get_frame(1)
            for p in range(clip.format.num_planes):
                assert np.array_equal(pixels(cut, p), pixels(spatial, p))
    base = core.std.BlankClip(width=3, height=3, length=1, format=vs.GRAY8, color=[100])
    invalid = [dict(mode='bad'), dict(rad=0), dict(rad=8), dict(thrY=-1),
               dict(thrUV=float('nan')), dict(tthY=float('inf')), dict(min=-0.01),
               dict(min=1.01), dict(scd=-1), dict(scenechange=2), dict(planes=[1]),
               dict(planes=[0, 0]), dict(planes=[])]
    if mode[1] == '3':
        invalid.append(dict(rad=5))
    for override in invalid:
        try:
            core.neo_smo.Deen(base, **dict(dict(mode=mode), **override))
        except vs.Error:
            pass
        else:
            raise AssertionError(('invalid argument accepted', override))
    # A single-frame 3D clip must not request nonexistent temporal dependencies.
    single = core.neo_smo.Deen(base, mode=mode).get_frame(0)
    assert (pixels(single, 0) == 100).all()
    print('VS', mode, 'passed')

if __name__ == '__main__':
    main()
