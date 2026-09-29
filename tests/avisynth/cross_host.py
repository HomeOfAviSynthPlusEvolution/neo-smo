"""Compare both host adapters on identical nonconstant frames (NumPy + VS required)."""
import argparse
import ctypes as C
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from isolated import suppress_crash_ui


class Data(C.Union):
    _fields_ = [("pointer", C.c_void_p), ("string", C.c_char_p), ("integer", C.c_int),
                ("floating", C.c_float)]


class Value(C.Structure):
    _fields_ = [("type", C.c_short), ("size", C.c_short), ("data", Data)]


class Runtime:
    def __init__(self, path):
        self.lib = (C.WinDLL if os.name == "nt" else C.CDLL)(path)
        self.bind("create_script_environment", C.c_void_p, C.c_int)
        self.bind("delete_script_environment", None, C.c_void_p)
        self.bind("invoke", Value, C.c_void_p, C.c_char_p, Value, C.c_void_p)
        self.bind("take_clip", C.c_void_p, Value, C.c_void_p)
        self.bind("release_value", None, Value)
        self.bind("release_clip", None, C.c_void_p)
        self.bind("get_frame", C.c_void_p, C.c_void_p, C.c_int)
        self.bind("clip_get_error", C.c_char_p, C.c_void_p)
        self.bind("release_video_frame", None, C.c_void_p)
        for key in ("pitch", "row_size", "height"):
            self.bind("get_" + key + "_p", C.c_int, C.c_void_p, C.c_int)
        self.bind("get_read_ptr_p", C.c_void_p, C.c_void_p, C.c_int)
        self.env = self.create_script_environment(11)
        assert self.env

    def bind(self, name, restype, *argtypes):
        fn = getattr(self.lib, "avs_" + name)
        fn.restype, fn.argtypes = restype, list(argtypes)
        setattr(self, name, fn)

    def evaluate(self, script):
        arg = Value(ord("s"), 0, Data(string=script.encode()))
        value = self.invoke(self.env, b"Eval", arg, None)
        try:
            assert value.type != ord("e"), value.data.string.decode(errors="replace")
            assert value.type == ord("c")
            return self.take_clip(value, self.env)
        finally:
            self.release_value(value)

    def pixels(self, clip, n, ids, dtype):
        import numpy as np
        frame = self.get_frame(clip, n)
        assert frame, self.clip_get_error(clip)
        try:
            result = []
            for plane in ids:
                pitch = self.get_pitch_p(frame, plane)
                height = self.get_height_p(frame, plane)
                row = self.get_row_size_p(frame, plane)
                pointer = self.get_read_ptr_p(frame, plane)
                data = b"".join(C.string_at(pointer + y * pitch, row) for y in range(height))
                result.append(np.frombuffer(data, dtype=dtype).reshape(height, -1))
            return result
        finally:
            self.release_video_frame(frame)


def check(args):
    import numpy as np
    import vapoursynth as vs
    core = vs.core
    core.num_threads = 2
    core.std.LoadPlugin(path=str(Path(args.plugin).resolve()))
    runtime = Runtime(args.runtime)
    pixel = args.case
    rgb = pixel.startswith("RGB")
    floating = pixel.endswith("PS")
    bits = 32 if floating else 16 if pixel.endswith("16") else 10 if pixel.endswith("10") else 8
    dtype = np.float32 if floating else np.uint16 if bits > 8 else np.uint8
    ids = (32, 64, 128) if rgb else (1, 2, 4)
    sub = 1 if "420" in pixel else 0
    fmt = core.query_video_format(vs.RGB if rgb else vs.YUV, vs.FLOAT if floating else vs.INTEGER, bits, sub, sub)
    convert = ".ConvertToPlanarRGB()" if rgb else ".ConvertToYV12()" if sub else ""
    setup = (f'LoadPlugin("{Path(args.plugin).resolve().as_posix()}")\n'
             'a=ColorBars(width=96,height=64,pixel_type="YV24").Trim(0,-1)\n'
             'c=a+a.Invert()+a.FlipHorizontal()+a+a.Invert()+a.FlipVertical()+a\n'
             f'c=c{convert}.ConvertBits({bits}).Crop(2,2,66,50)\n')
    source = runtime.evaluate(setup + "return c")
    inputs = [runtime.pixels(source, n, ids, dtype) for n in range(7)]
    blank = core.std.BlankClip(width=66, height=50, length=7, format=fmt.id)

    def fill(n, f):
        out = f.copy()
        for p in range(3):
            np.copyto(np.asarray(out[p]), inputs[n][p])
        return out

    clip = core.std.ModifyFrame(blank, blank, fill)
    cases = [("Median", {"radius": [1, 2, 3]}), ("InterQuartileMean", {"radius": [2, 1, 3]}),
             ("SmartMedian", {"radius": 2, "threshold": 40, "scalep": True}),
             ("VerticalCleaner", {"mode": [2, 1, 0]}), ("RemoveGrain", {"mode": [20, 13, 24]}),
             ("Repair", {"mode": [13, 16, 24], "repairclip": "REF"}),
             ("Clense", {}), ("ForwardClense", {}), ("BackwardClense", {}),
             ("TemporalMedian", {"radius": 2}), ("TemporalSoften", {"radius": 2, "threshold": 180, "scalep": True}),
             ("TemporalRepair", {"mode": [0, 2, 4], "repairclip": "REF"}),
             ("DegrainMedian", {"mode": [0, 3, 5], "limit": 70, "scalep": True}),
             ("FluxSmoothT", {"temporal_threshold": 120, "scalep": True}),
             ("FluxSmoothST", {"temporal_threshold": 120, "spatial_threshold": 90, "scalep": True}),
             ("TTempSmooth", {"maxr": 2, "thresh": [200, 80, 50], "mdiff": [20, 79, 2], "scthresh": 0}),
             ("TTempSmooth", {"maxr": 2, "pfclip": "REF", "fp": False, "scthresh": 0}),
             ("CCD", {"scale": 1, "temporal_radius": 1, "ref": "REF"}),
             ("DCTFilter", {"factors": [1, 0.9, 0.75, 0.5, 0.4, 0.25, 0.1, 0]})]
    for mode in ("c2d", "c3d", "w2d", "w3d", "a2d", "a3d"):
        cases.append(("Deen", {"mode": mode, "radius": [2, 0, 1], "threshold": [30, 40],
                               "temporal_threshold": [25, 35], "minimum": [0.3, 0.7], "scalep": True, "planes": [0, 2]}))
    cases.append(("Deen", {"mode": "a3d", "scenechange": -1}))
    if not rgb:
        for sc in (1, 254):
            cases.append(("Deen", {"scenechange": sc, "threshold": 255,
                                   "temporal_threshold": 255, "scalep": True}))
    cases.append(("MiniDeen", {"radius": [1, 3, 2], "threshold": [20, 30], "scalep": True, "planes": [0, 2]}))
    cases.append(("MiniDeen", {}))
    if not rgb and not floating:
        for tmode in (0, 1, 4):
            cases.append(("Cnr4", {"scenechange": False, "tmode": tmode, "wmode": tmode % 4, "ref": "REF"}))
    # Flip both axes without changing dimensions; reference pitches differ from cropped source pitches.
    ref = clip.std.FlipHorizontal().std.FlipVertical()

    def avs_value(v):
        if v == "REF":
            return "c.FlipHorizontal().FlipVertical()"
        if isinstance(v, str):
            return '"' + v + '"'
        if isinstance(v, bool):
            return str(v).lower()
        return str(v)

    def resized_guide(expression):
        guide = runtime.evaluate(setup + f"return {expression}.ExtractY().BilinearResize(33,25)")
        planes = [runtime.pixels(guide, n, (1,), dtype)[0] for n in range(7)]
        runtime.release_clip(guide)
        return planes

    def with_resized_luma(node, luma):
        uv = node.std.ShufflePlanes(planes=[1, 1, 2], colorfamily=vs.YUV)
        def replace(n, f):
            out = f.copy()
            np.copyto(np.asarray(out[0]), luma[n])
            return out
        return core.std.ModifyFrame(uv, uv, replace)

    resized_source = with_resized_luma(clip, resized_guide("c")) if sub else None
    resized_ref = with_resized_luma(ref, resized_guide("c.FlipHorizontal().FlipVertical()")) if sub else None
    cases += [("TemporalMedian", {"radius": 2, "scenechange": True}),
              ("TemporalSoften", {"radius": 2, "threshold": 180, "scalep": True, "scenechange": -1}),
              ("TTempSmooth", {"maxr": 2, "thresh": 200, "scthresh": -1, "pfclip": "REF"})]
    if not rgb and not floating:
        cases.append(("Cnr4", {"scenechange": True, "ref": "REF"}))
    total = 0
    for name, params in cases:
        use_scene = bool(params.get("scenechange", False)) or params.get("scthresh", 0) == -1
        scene_setup = 'c=c.propSet("_SceneChangePrev",1).propSet("_SceneChangeNext",1)\n' if use_scene else ""
        actual = runtime.evaluate(setup + scene_setup + "return neo_smo_" + name + "(c" +
                                  "".join(f",{k}={avs_value(v)}" for k, v in params.items()) + ")")
        kwargs = {k: ref if v == "REF" else v for k, v in params.items()}
        # Cnr4's arithmetic uses chroma-sized luma. Feed VS exactly the AVS
        # bilinear guide to isolate the adapter from host resizer rounding.
        aligned_cnr = sub and name == "Cnr4"
        if aligned_cnr:
            kwargs["ref"] = resized_ref
        source_node = resized_source if aligned_cnr else clip
        if use_scene:
            source_node = source_node.std.SetFrameProps(_SceneChangePrev=1, _SceneChangeNext=1)
            for key in ("ref", "pfclip"):
                if key in kwargs:
                    kwargs[key] = kwargs[key].std.SetFrameProps(_SceneChangePrev=1, _SceneChangeNext=1)
        expected = getattr(core.neo_smo, name)(source_node, **kwargs)
        for n in (0, 2, 3, 6):
            avs = runtime.pixels(actual, n, ids, dtype)
            with expected.get_frame(n) as frame:
                for p in range(3):
                    a = avs[p].astype(np.float64)
                    b = (inputs[n][p] if aligned_cnr and p == 0 else np.asarray(frame[p])).astype(np.float64)
                    diff = np.max(np.abs(a - b))
                    tolerance = 2e-6 if floating else 1 if name == "DCTFilter" else 0
                    assert diff <= tolerance, (pixel, name, params, n, p, diff)
            total += 1
        runtime.release_clip(actual)
    runtime.release_clip(source)
    runtime.delete_script_environment(runtime.env)
    print(f"{pixel}: {total} cross-host frames passed", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--plugin", required=True)
    parser.add_argument("--runtime", required=True)
    parser.add_argument("--case")
    args = parser.parse_args()
    suppress_crash_ui()
    if args.case:
        check(args)
        return
    for pixel in ("YUV444P8", "YUV444P10", "YUV444P16", "YUV444PS", "RGBP", "RGBPS", "YUV420P8", "YUV420PS"):
        subprocess.run([sys.executable, __file__, "--plugin", args.plugin, "--runtime", args.runtime,
                        "--case", pixel], check=True, timeout=45)


if __name__ == "__main__":
    main()
