"""MiniDeen scalar semantics and real VS binding regressions."""
import argparse
from pathlib import Path
import numpy as np
import vapoursynth as vs


def reference(a, radius, threshold, bits):
    threshold = threshold * 2**(bits-8)
    result = np.empty_like(a)
    for y in range(a.shape[0]):
        for x in range(a.shape[1]):
            c = int(a[y, x])
            window = a[max(0,y-radius):y+radius+1, max(0,x-radius):x+radius+1].astype(np.int64)
            accepted = window[np.abs(window-c) < threshold]
            total = 2*c+int(accepted.sum())
            count = 2+accepted.size
            result[y,x] = (2*total+count)//(2*count)
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--plugin', required=True)
    args = ap.parse_args()
    core = vs.core
    core.num_threads = 2
    core.std.LoadPlugin(path=str(Path(args.plugin).resolve()))
    for fmt in (vs.GRAY8, vs.GRAY16, vs.YUV420P8, vs.YUV420P10, vs.YUV420P16, vs.RGB48):
        base = core.std.BlankClip(width=8, height=6, length=3, format=fmt)
        def fill(n, f):
            out = f.copy()
            for p in range(out.format.num_planes):
                a = np.asarray(out[p])
                a[:] = (np.arange(a.size).reshape(a.shape)*3+n*5+p*7)%41 + (1 << (out.format.bits_per_sample-1))
            out.props['MiniDeenFrame'] = n
            return out
        clip = core.std.ModifyFrame(base, clips=base, selector=fill)
        for radius, threshold, planes in (([1], [10], None), ([7,2], [1,11], [0]), ([2,3,1], [0,255,10], None)):
            radius = radius[:clip.format.num_planes]
            threshold = threshold[:clip.format.num_planes]
            out = core.neo_smo.MiniDeen(clip, radius=radius, threshold=threshold, scalep=True, planes=planes)
            for n in (2,0,1,1):
                with out.get_frame(n) as actual, clip.get_frame(n) as source:
                    assert actual.props['MiniDeenFrame'] == n
                    for p in range(source.format.num_planes):
                        a = np.asarray(source[p])
                        want = reference(a, radius[min(p,len(radius)-1)], threshold[min(p,len(threshold)-1)], source.format.bits_per_sample) if planes is None or p in planes else a
                        np.testing.assert_array_equal(np.asarray(actual[p]), want)
        default = core.neo_smo.MiniDeen(clip)
        with default.get_frame(0) as a, clip.get_frame(0) as source:
            for p in range(source.format.num_planes):
                np.testing.assert_array_equal(np.asarray(a[p]), reference(np.asarray(source[p]),1,10,source.format.bits_per_sample))
    # Float formats use the same strict window equation, including signed chroma.
    for fmt in (vs.GRAYH, vs.GRAYS, vs.YUV444PH, vs.YUV420PS, vs.RGBS):
        base = core.std.BlankClip(width=66, height=6, length=1, format=fmt)
        def pix(f, p):
            a = np.asarray(f[p])
            return a.view(np.float16) if f.format.bits_per_sample == 16 else a
        def fill_float(n, f):
            out = f.copy()
            for p in range(out.format.num_planes):
                a = pix(out, p)
                a[:] = ((np.arange(a.size).reshape(a.shape)*7+p*11)%41)/64 - (0.5 if p and out.format.color_family == vs.YUV else 0)
            return out
        clip = core.std.ModifyFrame(base, clips=base, selector=fill_float)
        for r in (1, 7):
            with clip.get_frame(0) as source, core.neo_smo.MiniDeen(clip, radius=r, threshold=0.125).get_frame(0) as actual:
                for p in range(source.format.num_planes):
                    a = pix(source, p).astype(np.float64)
                    want = np.empty_like(a)
                    for y in range(a.shape[0]):
                        for x in range(a.shape[1]):
                            c = a[y,x]
                            window = a[max(0,y-r):y+r+1, max(0,x-r):x+r+1]
                            accepted = window[np.abs(window-c) < 0.125]
                            want[y,x] = (2*c+accepted.sum())/(2+accepted.size)
                    np.testing.assert_allclose(pix(actual,p), want, rtol=0, atol=0.00025 if source.format.bits_per_sample == 16 else 1e-7)
    from deen_api import check_api, check_scene_fallback
    check_api(core)
    check_scene_fallback(core)
    # Tiny planes and half-up center-weight fixture.
    base = core.std.BlankClip(width=2,height=1,format=vs.GRAY8,length=1)
    def fixture(n,f):
        out=f.copy(); np.asarray(out[0])[:]=[[100,110]]; return out
    clip=core.std.ModifyFrame(base,clips=base,selector=fixture)
    with core.neo_smo.MiniDeen(clip,threshold=11).get_frame(0) as f:
        np.testing.assert_array_equal(np.asarray(f[0]), [[103,108]])
    tiny=core.std.BlankClip(width=1,height=1,format=vs.GRAY16,length=1,color=[65535])
    with core.neo_smo.MiniDeen(tiny,radius=7,threshold=255).get_frame(0) as f:
        assert np.asarray(f[0])[0,0] == 65535
    for kwargs in ({'radius':[-1]}, {'radius':[8]}, {'threshold':[-1]}, {'threshold':[65536]},
                   {'radius':[1,2,3,4]}, {'threshold':[1,2,3,4]}, {'planes':[0,0]}, {'planes':[1]}, {'planes':[]}, {'radius':[]}, {'threshold':[]}, {'threshold':[float('nan')]}):
        try: core.neo_smo.MiniDeen(tiny,**kwargs)
        except vs.Error: pass
        else: raise AssertionError(f'accepted {kwargs}')
    print('MiniDeen VS checks passed')

if __name__ == '__main__': main()
