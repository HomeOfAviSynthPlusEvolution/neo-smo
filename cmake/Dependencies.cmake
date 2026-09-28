include(FetchContent)

# Local source overrides are explicit FETCHCONTENT_SOURCE_DIR_<NAME> options,
# as in neo-mv and neo-fft; do not silently select another project's build tree.

if(NEO_SMO_BUILD_AVISYNTH OR NEO_SMO_BUILD_VAPOURSYNTH)
  set(DS_BUILD_TESTS OFF CACHE BOOL "" FORCE)
  set(DS_BUILD_ACCEPTANCE_PLUGIN OFF CACHE BOOL "" FORCE)
  set(DS_ENABLE_AVISYNTH ${NEO_SMO_BUILD_AVISYNTH} CACHE BOOL "" FORCE)
  set(DS_ENABLE_VAPOURSYNTH ${NEO_SMO_BUILD_VAPOURSYNTH} CACHE BOOL "" FORCE)
  set(_neo_smo_ds_subdir .)
else()
  set(_neo_smo_ds_subdir neo_smo_headers_only)
endif()

FetchContent_Declare(dualsynth2
  GIT_REPOSITORY https://github.com/HomeOfAviSynthPlusEvolution/dualsynth2.git
  GIT_TAG 627127332a7a360f8cbf644be545115bc16585cd
  SOURCE_SUBDIR "${_neo_smo_ds_subdir}")
FetchContent_MakeAvailable(dualsynth2)
unset(_neo_smo_ds_subdir)

foreach(option TESTS EXAMPLES CONTRIB INSTALL)
  set(HWY_ENABLE_${option} OFF CACHE BOOL "" FORCE)
endforeach()
set(HWY_FORCE_STATIC_LIBS ON CACHE BOOL "" FORCE)

FetchContent_Declare(highway
  URL https://github.com/google/highway/archive/refs/tags/1.4.0.tar.gz
  URL_HASH SHA256=e72241ac9524bb653ae52ced768b508045d4438726a303f10181a38f764a453c
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
FetchContent_MakeAvailable(highway)
set_target_properties(hwy PROPERTIES POSITION_INDEPENDENT_CODE ON
  CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN ON)

set(NEO_DCT_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(NEO_DCT_ALLOW_FMA ON CACHE BOOL "" FORCE)
FetchContent_Declare(neo_libdct
  GIT_REPOSITORY https://github.com/HomeOfAviSynthPlusEvolution/neo-libdct.git
  GIT_TAG 52d3712374309c08cc4f724b88f17a0124f17325)
FetchContent_MakeAvailable(neo_libdct)
