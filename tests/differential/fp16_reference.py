"""Match the production FP16 arithmetic policy without widening tolerances."""
import subprocess
from pathlib import Path
import numpy as np
import vapoursynth as vs

class Reference:
    def __init__(self, core, plugin):
        self.core = core
        suffix = ".exe" if Path(plugin).suffix == ".dll" else ""
        probe = Path(plugin).resolve().parent / ("neo_smo_fp16_tests" + suffix)
        mode = subprocess.check_output([str(probe), "--capability"], text=True, timeout=10).strip()
        if mode not in ("native", "fp32"):
            raise RuntimeError("Unexpected FP16 capability: " + mode)
        self.native = mode == "native"
        print("[FP16 reference] " + mode, flush=True)

    def convert(self, clip, bits):
        fmt = clip.format.replace(bits_per_sample=bits)
        base = self.core.std.BlankClip(clip, format=fmt.id)
        def fill(n, f):
            out = f[0].copy()
            for p in range(fmt.num_planes):
                src = np.asarray(f[1][p])
                dst = np.asarray(out[p])
                if clip.format.bits_per_sample == 16:
                    src = src.view(np.float16)
                if bits == 16:
                    dst = dst.view(np.float16)
                dst[:] = src
            for key in f[1].props:
                out.props[key] = f[1].props[key]
            return out
        return self.core.std.ModifyFrame(base, clips=[base, clip], selector=fill)

    def __getattr__(self, name):
        fn = getattr(self.core.zsmooth, name)
        def call(*args, **kwargs):
            clips = [a for a in (*args, *kwargs.values()) if isinstance(a, vs.VideoNode)]
            half = any(c.format.sample_type == vs.FLOAT and c.format.bits_per_sample == 16 for c in clips)
            if self.native or not half:
                return fn(*args, **kwargs)
            def lift(v):
                return self.convert(v, 32) if isinstance(v, vs.VideoNode) and v.format.sample_type == vs.FLOAT and v.format.bits_per_sample == 16 else v
            result = fn(*(lift(a) for a in args), **{k: lift(v) for k, v in kwargs.items()})
            return self.convert(result, 16)
        return call

    def smart_median_threshold_match(self, clip, radius, candidate):
        """Allow only the two decisions at adjacent binary16 thresholds.

        Native half arithmetic and optimized Zig can round the variance to
        opposite sides of a threshold. Do not allow an arbitrary pixel error.
        This check is only for the default thresholds used by the radius sweep.
        """
        if not self.native or clip.format.bits_per_sample != 16 or clip.format.sample_type != vs.FLOAT:
            return False
        radii = list(radius) if isinstance(radius, (list, tuple)) else [radius]
        radii += [radii[-1]] * (clip.format.num_planes - len(radii))
        thresholds = np.asarray([50 / 255 if r == 1 else 128 / 255 for r in radii], dtype=np.float16)
        low = np.nextafter(thresholds, np.float16(0)).astype(float).tolist()
        high = np.nextafter(thresholds, np.float16(1)).astype(float).tolist()
        lo = self.core.zsmooth.SmartMedian(clip, radius=radius, threshold=low, scalep=False).get_frame(0)
        hi = self.core.zsmooth.SmartMedian(clip, radius=radius, threshold=high, scalep=False).get_frame(0)
        nominal = self.core.zsmooth.SmartMedian(clip, radius=radius).get_frame(0)
        changed = 0
        for p in range(clip.format.num_planes):
            a = np.asarray(candidate[p]).view(np.float16)
            l = np.asarray(lo[p]).view(np.float16)
            h = np.asarray(hi[p]).view(np.float16)
            if not np.all((a == l) | (a == h)):
                return False
            changed += int(np.count_nonzero(a != np.asarray(nominal[p]).view(np.float16)))
        print(f"[FP16 threshold rounding] SmartMedian radius={radius}: {changed} pixels, within +/-1 threshold ULP", flush=True)
        return True
