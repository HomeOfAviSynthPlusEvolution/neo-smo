"""Bounded AviSynth host tests; every child inherits crash UI suppression."""
import argparse
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from isolated import suppress_crash_ui


def main():
    parser = argparse.ArgumentParser()
    for name in ("runner", "plugin", "runtime", "work"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--host-api", choices=("c", "cpp"), default="c")
    args = parser.parse_args()
    suppress_crash_ui()
    work = Path(args.work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    header = f'LoadPlugin("{Path(args.plugin).resolve().as_posix()}")\n'
    count = 0

    def run(name, body, *, pixel="YUV420P8", setup=None, frame=4, error=None, extra=(), mode="--video"):
        nonlocal count
        if setup is None:
            setup = f'c=BlankClip(width=64,height=48,length=9,pixel_type="{pixel}")\n'
        script = work / (name + ".avs")
        script.write_text(header + setup + body, encoding="utf-8")
        result = subprocess.run([args.runner, mode, str(script), "--backend", args.host_api,
                                 "--runtime", args.runtime, "--frame", str(frame), *extra],
                                capture_output=True, text=True, timeout=15)
        text = result.stdout + result.stderr
        if error:
            assert result.returncode != 0 and error.lower() in text.lower(), (name, text)
        else:
            assert result.returncode == 0, (name, text)
        count += 1
        return result.stdout

    filters = {
        "Median": "", "InterQuartileMean": "", "SmartMedian": "",
        "VerticalCleaner": ",mode=2", "RemoveGrain": ",mode=20",
        "Repair": ",repairclip=c,mode=13", "Clense": "",
        "ForwardClense": "", "BackwardClense": "", "TemporalMedian": "",
        "TemporalSoften": "", "TemporalRepair": ",repairclip=c",
        "DegrainMedian": "", "FluxSmoothT": "", "FluxSmoothST": "",
        "TTempSmooth": ",scthresh=0", "CCD": ",scale=1",
        "Deen": "", "MiniDeen": "",
        "Cnr4": ",scenechange=false", "DCTFilter": ",factors=[1,1,1,1,1,1,1,1]",
    }
    run("registration", "\n".join(f'Assert(FunctionExists("neo_smo_{f}"))' for f in filters) + "\nreturn c")
    for pixel in ("Y8", "YUV420P8", "YUV422P10", "YUV444P16", "Y32", "YUV420PS", "RGBP", "RGBPS"):
        for f, params in filters.items():
            if f == "CCD" and pixel in ("Y8", "Y32"):
                continue
            if f == "Cnr4" and (not pixel.startswith("YUV") or pixel.endswith("PS")):
                continue
            run(f"{f}-{pixel}", f"return neo_smo_{f}(c{params}).Prefetch(2)", pixel=pixel)

    varying = ('c=ColorBars(width=96,height=64,pixel_type="YV24").Trim(0,3)\n'
               'c=(c+c.Invert()+c.Trim(0,0)).Crop(2,2,64,48)\n')
    for f, params in filters.items():
        for frame in (0, 4, 8):
            run(f"crop-{f}-{frame}", f"return neo_smo_{f}(c{params})", setup=varying, frame=frame)
    run("ttemp-auto", "return neo_smo_TTempSmooth(c)")
    run("soften-auto", "return neo_smo_TemporalSoften(c,scenechange=12)")
    for f, params in (("TemporalMedian", ",scenechange=true"), ("Cnr4", "")):
        run(f"scene-missing-{f}", f"return neo_smo_{f}(c{params})", error="scene change handling requires")
        run(f"scene-present-{f}", f'c=c.propSet("_SceneChangePrev",0).propSet("_SceneChangeNext",0)\nreturn neo_smo_{f}(c{params})')
    for f in ("Median", "TTempSmooth", "CCD", "Cnr4", "DCTFilter", "Deen", "MiniDeen"):
        params = filters[f]
        audio = ('c=BlankClip(width=64,height=48,length=9,pixel_type="YV12")\n'
                 'c=AudioDub(c,Tone(length=1.0,samplerate=48000,channels=2)).AssumeTFF()\n')
        body = f'o=neo_smo_{f}(c{params})\nAssert(AudioRate(o)==AudioRate(c) && AudioChannels(o)==AudioChannels(c))\nAssert(GetParity(o,4)==GetParity(c,4))\n'
        source = run(f + "-audio-source", body + "return c", setup=audio, mode="--audio")
        output = run(f + "-audio-output", body + "return o.Prefetch(2)", setup=audio, mode="--audio")
        assert source == output, (f, source, output)
        run(f + "-properties", f'''c=c.propSet("NeoSMOTest",123)
o=neo_smo_{f}(c{params})
return o.ScriptClip("""Assert(propGetInt(last,"NeoSMOTest")==123)
last""").Prefetch(2)
''')
    invalid = {
        "required": ('Median()', "clip is required"),
        "mode-missing": ('RemoveGrain(c)', "mode is required"),
        "radius": ('Median(c,radius=4)', "radius"),
        "duplicate": ('Median(c,planes=[0,0])', "plane"),
        "plane-index": ('DCTFilter(c,factors=[1,1,1,1,1,1,1,1],planes=3)', "planes"),
        "empty": ('Median(c,planes=[])', "cannot be empty"),
        "type": ('Median(c,radius=[1,true])', "integers"),
        "factors": ('DCTFilter(c,factors=[1,1])', "eight"),
        "ref-size": ('Repair(c,c.Crop(0,0,32,32),1)', "reference format"),
        "short-ref": ('Repair(c,c.Trim(0,1),1)', "requested frame"),
        "cnr-mode": ('Cnr4(c,mode="abc")', "mode only"),
        "rgb-auto": ('TTempSmooth(c.ConvertToPlanarRGB())', "automatic scene detection"),
    }
    for name, (expr, error) in invalid.items():
        run("error-" + name, "return neo_smo_" + expr, error=error)
    for pixel in ("RGB32", "YUY2", "RGBAP", "YUVA444"):
        run("unsupported-" + pixel, "return neo_smo_Median(c)", pixel=pixel, error="only planar")

    temporal = "c=" + " + ".join(
        f'BlankClip(width=16,height=16,length=1,pixel_type="Y8",color_yuv=${value:02x}8080)'
        for value in (10, 90, 30)) + "\n"
    for name, expr, frame, expected in (
        ("median", "TemporalMedian(c)", 1, 30),
        ("clense", "Clense(c)", 1, 30),
        ("soften", "TemporalSoften(c,radius=1,threshold=255)", 1, 43),
        ("soften-edge", "TemporalSoften(c,radius=1,threshold=255)", 0, 50),
        ("median-edge", "TemporalMedian(c)", 0, 10),
        ("soften-scene", "TemporalSoften(c,radius=1,threshold=255,scenechange=12)", 1, 90),
        ("ttemp-scene", "TTempSmooth(c,maxr=1,thresh=256,scthresh=12)", 1, 90),
    ):
        run("value-" + name, "return neo_smo_" + expr, setup=temporal, frame=frame,
            extra=("--expect-y8-sum", str(16 * 16 * expected)))
    for mode in ("c2d", "c3d", "w2d", "w3d", "a2d", "a3d"):
        expected = 90 if mode.endswith("2d") else 55 if mode.startswith("w") else 43
        run("deen-value-" + mode,
            f'return neo_smo_Deen(c,mode="{mode}",threshold=255,temporal_threshold=255,scenechange=0,minimum=1)',
            setup=temporal, frame=1, extra=("--expect-y8-sum", str(256 * expected)))
        run("deen-scene-" + mode,
            f'return neo_smo_Deen(c,mode="{mode}",threshold=255,temporal_threshold=255,scenechange=1)',
            setup=temporal, frame=1, extra=("--expect-y8-sum", str(256 * 90)))
    for name, expr, error in (
        ("deen-mode", 'Deen(c,mode="bad")', "mode"),
        ("deen-radius", 'Deen(c,radius=5)', "radius"),
        ("deen-threshold", 'Deen(c,threshold=-1)', "threshold"),
        ("deen-planes", 'Deen(c,planes=[0,0])', "planes"),
        ("mini-radius", 'MiniDeen(c,radius=[1,8])', "radius"),
        ("mini-threshold", 'MiniDeen(c,threshold=-1)', "threshold"),
        ("mini-planes", 'MiniDeen(c,planes=3)', "planes"),
        ("mini-float-range", 'MiniDeen(c.ConvertBits(32),threshold=2)', "threshold"),
    ):
        run(name, "return neo_smo_" + expr, error=error)
    for name in ("Deen", "MiniDeen"):
        for pixel in ("RGB32", "YUY2", "RGBAP", "YUVA444"):
            run(name + "-reject-" + pixel, "return neo_smo_" + name + "(c)", pixel=pixel, error="only planar")
        run(name + "-empty-planes", "return neo_smo_" + name + "(c,planes=[])", error="planes cannot be empty")
    for name in ("Deen", "MiniDeen"):
        for key in ("radius", "threshold", "planes"):
            run(name + "-empty-" + key, f"return neo_smo_{name}(c,{key}=[])" , error=key)
        for key in ("radius", "threshold"):
            run(name + "-gray-array-" + key, f"return neo_smo_{name}(c,{key}=[1,1])", pixel="Y8", error=key)
        run(name + "-radius-zero", f"return neo_smo_{name}(c,radius=0)", pixel="Y32")
        run(name + "-float-scalep", f"return neo_smo_{name}(c,threshold=10,scalep=true)", pixel="RGBPS")
    for key in ("minimum", "temporal_threshold"):
        run("deen-empty-" + key, f"return neo_smo_Deen(c,{key}=[])" , error=key)
    run("deen-rgb-auto", 'return neo_smo_Deen(c,scenechange=12)', pixel="RGBP", error="Gray or YUV")
    run("deen-scene-missing", 'return neo_smo_Deen(c,scenechange=-1)')
    run("deen-order", 'return neo_smo_Deen(c,"a3d",[1,0],[7,9],[4,6],[0.5],0,true,[0])')
    run("mini-order", 'return neo_smo_MiniDeen(c,[1,0],[10,6],true,[0])')
    print(f"AviSynth acceptance: {count} cases passed")


if __name__ == "__main__":
    main()
