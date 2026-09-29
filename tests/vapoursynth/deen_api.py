"""Public contracts shared by Deen and MiniDeen, exercised in the real host."""
import numpy as np
import vapoursynth as vs


def check_api(core):
    def same(a, b):
        with a.get_frame(1) as x, b.get_frame(1) as y:
            for p in range(x.format.num_planes):
                np.testing.assert_array_equal(np.asarray(x[p]), np.asarray(y[p]))

    def rejects(fn, clip, **params):
        try:
            fn(clip, **params)
        except (vs.Error, TypeError, ValueError):
            return
        raise AssertionError(('accepted invalid parameters', fn.name, params))

    for fmt in (vs.GRAY8, vs.GRAY16, vs.YUV420P16, vs.YUV420PS, vs.RGBS):
        base = core.std.BlankClip(width=18, height=6, length=3, format=fmt)
        fp = base.format.sample_type == vs.FLOAT
        scale = 1/255 if fp else 2**(base.format.bits_per_sample-8)

        def fill(n, f):
            out = f.copy()
            for p in range(out.format.num_planes):
                a = np.asarray(out[p])
                a[:] = ((np.arange(a.size).reshape(a.shape)*3+n*7+p*5)%31)*scale
                if fp and p and out.format.color_family == vs.YUV:
                    a[:] -= 0.5
            return out

        clip = core.std.ModifyFrame(base, base, fill)
        for mini in (False, True):
            fn = core.neo_smo.MiniDeen if mini else core.neo_smo.Deen
            keys = ['radius', 'threshold'] + ([] if mini else ['temporal_threshold', 'minimum'])
            for key in keys:
                rejects(fn, clip, **{key: []})
                rejects(fn, clip, **{key: [1]*(clip.format.num_planes+1)})
            for key in ['threshold'] + ([] if mini else ['temporal_threshold', 'minimum']):
                for bad in (float('nan'), float('inf'), -1):
                    rejects(fn, clip, **{key: bad})
            rejects(fn, clip, radius='invalid')
            rejects(fn, clip, planes=[])
            rejects(fn, clip, planes=[0,0])
            rejects(fn, clip, threshold=256, scalep=True)
            rejects(fn, clip, threshold=2 if fp else 2**clip.format.bits_per_sample)
            same(fn(clip, radius=0), clip)
            same(fn(clip, threshold=7, scalep=True), fn(clip, threshold=7*scale))
            same(fn(clip, radius=1, threshold=7*scale), fn(clip, radius=[1], threshold=[7*scale]))
            defaults = [10 if mini else 9 if p and clip.format.color_family == vs.YUV else 7
                        for p in range(clip.format.num_planes)]
            same(fn(clip), fn(clip, threshold=defaults, scalep=True))
            if clip.format.num_planes == 3:
                same(fn(clip, radius=[2,0], threshold=[7*scale,9*scale]),
                     fn(clip, radius=[2,0,0], threshold=[7*scale,9*scale,9*scale]))
                # Unselected planes are still validated.
                rejects(fn, clip, planes=[0], threshold=[1,-1,1])
            if not mini:
                temporal_defaults = [6 if p and clip.format.color_family == vs.YUV else 4
                                     for p in range(clip.format.num_planes)]
                same(fn(clip), fn(clip, temporal_threshold=temporal_defaults, scalep=True))
                same(fn(clip, scenechange=-1), fn(clip))
                same(fn(clip, 'a3d', 1, 7, 4, 0.5, 0, True, [0]),
                     fn(clip, mode='a3d', radius=1, threshold=7, temporal_threshold=4,
                        minimum=0.5, scenechange=0, scalep=True, planes=[0]))
                if clip.format.color_family == vs.RGB:
                    rejects(fn, clip, scenechange=12)
                elif any(p.identifier == "com.vapoursynth.misc" for p in core.plugins()):
                    same(fn(clip, scenechange=12),
                         fn(clip.misc.SCDetect(threshold=12/255), scenechange=-1))
                else:
                    # With misc absent, exercise the production built-in fallback.
                    for sc in (1, 12, 254):
                        same(fn(clip, scenechange=sc), fn(mark_scene(core, clip, sc), scenechange=-1))
            else:
                same(fn(clip, 1, 10, True, [0]), fn(clip, radius=1, threshold=10, scalep=True, planes=[0]))

    # An integer difference of 10 is accepted by MiniDeen at 10.5, rejected at 10.
    base = core.std.BlankClip(width=2, height=1, length=3, format=vs.GRAY8)
    def pair(n, f):
        out = f.copy()
        np.asarray(out[0])[:] = [[100,110]]
        return out
    clip = core.std.ModifyFrame(base, base, pair)
    same(core.neo_smo.MiniDeen(clip, threshold=10), clip)
    with core.neo_smo.MiniDeen(clip, threshold=10.5).get_frame(1) as f:
        np.testing.assert_array_equal(np.asarray(f[0]), [[103,108]])
    same(core.neo_smo.Deen(clip, mode='c2d', threshold=9.5), clip)
    with core.neo_smo.Deen(clip, mode='c2d', threshold=10).get_frame(1) as f:
        np.testing.assert_array_equal(np.asarray(f[0]), [[103,107]])


def mark_scene(core, clip, threshold):
    def diff(a, b):
        with clip.get_frame(a) as fa, clip.get_frame(b) as fb:
            def luma(f):
                a = np.asarray(f[0])
                if clip.format.sample_type == vs.FLOAT and clip.format.bits_per_sample == 16:
                    a = a.view(np.float16)
                return a.astype(np.float64)
            peak = 1 if clip.format.sample_type == vs.FLOAT else (1 << clip.format.bits_per_sample)-1
            return np.abs(luma(fa)-luma(fb)).mean()/peak > threshold/255
    cuts = [int(diff(n, n+1)) for n in range(clip.num_frames-1)] + [0]
    def props(n, f):
        out = f.copy()
        out.props['_SceneChangePrev'] = cuts[max(n-1, 0)]
        out.props['_SceneChangeNext'] = cuts[n]
        return out
    return core.std.ModifyFrame(clip, clip, props)


def check_scene_fallback(core):
    for fmt in (vs.GRAY8, vs.GRAY16, vs.GRAYH, vs.GRAYS, vs.YUV420P10, vs.YUV420PH, vs.YUV420PS):
        fp = core.get_video_format(fmt).sample_type == vs.FLOAT
        peak = 1 if fp else (1 << core.get_video_format(fmt).bits_per_sample)-1
        frames = []
        for value in (0.1, 0.12, 0.7, 0.72):
            color = [value*peak if fp else int(value*peak)]
            if core.get_video_format(fmt).num_planes == 3:
                color += [0 if fp else int(peak/2)]*2
            frames.append(core.std.BlankClip(width=18, height=6, length=1, format=fmt, color=color))
        clip = core.std.Splice(frames)
        # Positive detection replaces pre-existing flags, even at the endpoints.
        clip = clip.std.SetFrameProps(_SceneChangePrev=1, _SceneChangeNext=1)
        for sc in (1, 12, 254):
            reference = mark_scene(core, clip, sc)
            for mode in ('c3d', 'a3d', 'w3d'):
                options = dict(mode=mode, threshold=255, temporal_threshold=255, scalep=True)
                auto = core.neo_smo.Deen(clip, scenechange=sc, **options)
                explicit = core.neo_smo.Deen(reference, scenechange=-1, **options)
                for n in (3, 0, 2, 1, 1):
                    with auto.get_frame(n) as a, explicit.get_frame(n) as b:
                        for p in range(a.format.num_planes):
                            np.testing.assert_array_equal(np.asarray(a[p]), np.asarray(b[p]))
                        for key in ('_SceneChangePrev', '_SceneChangeNext'):
                            assert a.props[key] == b.props[key], (fmt, sc, n, key)
