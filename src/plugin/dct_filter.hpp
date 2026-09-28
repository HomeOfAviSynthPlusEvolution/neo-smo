#pragma once
#include <vapoursynth/VapourSynth4.h>

namespace neo_smo::plugin {
void VS_CC dct_filter_create(const VSMap* in, VSMap* out, void*, VSCore* core, const VSAPI* vsapi);
} // namespace neo_smo::plugin
