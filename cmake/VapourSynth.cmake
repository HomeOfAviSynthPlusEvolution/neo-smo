set(NEO_SMO_VS_SDK "" CACHE PATH "Optional VapourSynth SDK")
find_path(NEO_SMO_VS_INCLUDE NAMES vapoursynth/VapourSynth4.h VapourSynth4.h
  HINTS "${NEO_SMO_VS_SDK}" "${NEO_SMO_VS_SDK}/include" "$ENV{VAPOURSYNTH_SDK}/include")
if(NOT NEO_SMO_VS_INCLUDE)
  FetchContent_Declare(neo_smo_vs_sdk
    URL https://github.com/vapoursynth/vapoursynth/archive/refs/tags/R73.zip
    URL_HASH SHA256=7c6b1eb2ec4aeae078675a29cf5e77c6ca1ab8b1dd322677c66e1ca6b76c511d
    SOURCE_SUBDIR neo_smo_headers_only)
  FetchContent_MakeAvailable(neo_smo_vs_sdk)
  set(NEO_SMO_VS_INCLUDE "${neo_smo_vs_sdk_SOURCE_DIR}/include")
endif()
if(NOT EXISTS "${NEO_SMO_VS_INCLUDE}/vapoursynth/VapourSynth4.h")
  configure_file("${NEO_SMO_VS_INCLUDE}/VapourSynth4.h"
    "${CMAKE_CURRENT_BINARY_DIR}/generated/include/vapoursynth/VapourSynth4.h" COPYONLY)
  set(NEO_SMO_VS_INCLUDE "${CMAKE_CURRENT_BINARY_DIR}/generated/include")
endif()
