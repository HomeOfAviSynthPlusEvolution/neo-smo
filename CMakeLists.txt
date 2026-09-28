cmake_minimum_required(VERSION 3.24)
project(neo_smo VERSION 0.9.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(NEO_SMO_BUILD_AVISYNTH "Build the AviSynth+ plugin" ON)
option(NEO_SMO_BUILD_VAPOURSYNTH "Build the VapourSynth v4 plugin" ON)

include(cmake/Dependencies.cmake)

add_library(neo_smo_core STATIC
  src/algorithms/mini_deen.cpp
  src/kernels/mini_deen.cpp
  src/algorithms/deen.cpp
  src/kernels/deen.cpp
  src/kernels/deen_c.cpp
  src/kernels/deen_w.cpp
  src/kernels/deen_a.cpp
  src/algorithms/plan.cpp
  src/algorithms/dct_filter.cpp
  src/kernels/dct_pixels.cpp
  src/kernels/median.cpp
  src/kernels/median_3x3.cpp
  src/kernels/median_5x5.cpp
  src/kernels/median_7x7.cpp
  src/kernels/vertical_cleaner.cpp
  src/kernels/remove_grain.cpp
  src/kernels/repair.cpp
  src/kernels/repair_int.cpp
  src/kernels/repair_f16.cpp
  src/kernels/repair_f32.cpp
  src/kernels/clense.cpp
  src/kernels/inter_quartile_mean.cpp
  src/kernels/smart_median.cpp
  src/kernels/temporal_median.cpp
  src/kernels/temporal_soften.cpp
  src/kernels/temporal_repair.cpp
  src/kernels/degrain_median.cpp
  src/kernels/degrain_median_u8.cpp
  src/kernels/degrain_median_u16.cpp
  src/kernels/degrain_median_f16.cpp
  src/kernels/degrain_median_f32.cpp
  src/kernels/fluxsmooth.cpp
  src/kernels/ttempsmooth.cpp
  src/kernels/ccd.cpp
  src/kernels/cnr4.cpp
)
target_compile_features(neo_smo_core PUBLIC cxx_std_17)
target_include_directories(neo_smo_core PUBLIC src)
target_include_directories(neo_smo_core SYSTEM PUBLIC "${dualsynth2_SOURCE_DIR}/include")
target_link_libraries(neo_smo_core PUBLIC hwy neo::dct)
set_target_properties(neo_smo_core PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN ON
)

if(MSVC)
  set(NEO_SMO_FP_OPTIONS /fp:fast)
  target_compile_options(neo_smo_core PRIVATE /W4)
else()
  set(NEO_SMO_FP_OPTIONS -ffp-contract=fast)
  target_compile_options(neo_smo_core PRIVATE -Wall -Wextra -Wpedantic)
endif()
target_compile_options(neo_smo_core PRIVATE ${NEO_SMO_FP_OPTIONS})

target_compile_definitions(neo_smo_core PRIVATE HWY_COMPILE_ALL_ATTAINABLE)

# Deen validates finite values and bounds rounding errors before exact boundary
# comparisons. Reassociation would invalidate those bounds and compensated sums.
# SIMD specialization is a separate, deferred change.
if(MSVC)
  set_property(SOURCE src/algorithms/deen.cpp src/kernels/deen_c.cpp src/kernels/deen_w.cpp src/kernels/deen_a.cpp APPEND PROPERTY COMPILE_OPTIONS /fp:precise)
else()
  set_property(SOURCE src/algorithms/deen.cpp src/kernels/deen_c.cpp src/kernels/deen_w.cpp src/kernels/deen_a.cpp APPEND PROPERTY COMPILE_OPTIONS -fno-fast-math -ffp-contract=off)
endif()


# Limit expansion of the EMU128 sorting network in MSVC. Hardware SIMD targets
# retain normal inlining and optimization.
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  set_property(SOURCE src/kernels/median_7x7.cpp
    src/kernels/repair_int.cpp src/kernels/repair_f16.cpp src/kernels/repair_f32.cpp
    src/kernels/inter_quartile_mean.cpp src/kernels/smart_median.cpp
    src/kernels/temporal_median.cpp
    APPEND PROPERTY COMPILE_DEFINITIONS
    "$<$<CONFIG:Release,RelWithDebInfo>:NEO_SMO_OUTLINE_EMU_COMPARE_SWAP>")
  set_property(SOURCE src/kernels/degrain_median_u8.cpp src/kernels/degrain_median_u16.cpp
    src/kernels/degrain_median_f16.cpp src/kernels/degrain_median_f32.cpp APPEND PROPERTY COMPILE_DEFINITIONS
    "$<$<CONFIG:Release,RelWithDebInfo>:NEO_SMO_MSVC_OUTLINE_DGM_EMU>")
  # Reduce integer Repair optimization cost only for the native MSVC compiler.
  set_property(SOURCE src/kernels/repair_int.cpp APPEND PROPERTY COMPILE_DEFINITIONS
    "$<$<CONFIG:Release,RelWithDebInfo>:NEO_SMO_MSVC_REPAIR_ABS_DIFF>")
  # Bound old-SSE integer sorting optimization cost without changing other compilers.
  set_property(SOURCE src/kernels/inter_quartile_mean.cpp src/kernels/smart_median.cpp
    APPEND PROPERTY COMPILE_DEFINITIONS
    "$<$<CONFIG:Release,RelWithDebInfo>:NEO_SMO_MSVC_U16_SORT>")
endif()

# EMU128 half rounding expands into large scalar control flow in MSVC/clang-cl.
# Keep it as a function boundary without restricting hardware SIMD inlining.
if(MSVC)
  set_property(SOURCE src/kernels/remove_grain.cpp src/kernels/repair_f16.cpp
    src/kernels/inter_quartile_mean.cpp src/kernels/smart_median.cpp
    src/kernels/temporal_median.cpp src/kernels/temporal_soften.cpp
    src/kernels/temporal_repair.cpp src/kernels/degrain_median_f16.cpp src/kernels/fluxsmooth.cpp
    src/kernels/ttempsmooth.cpp src/kernels/ccd.cpp src/kernels/cnr4.cpp
    APPEND PROPERTY COMPILE_DEFINITIONS
    "$<$<CONFIG:Release,RelWithDebInfo>:NEO_SMO_OUTLINE_EMU_ROUND_F16>")
endif()

include(CTest)
if(BUILD_TESTING)
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
  foreach(test boundaries fp16 float_accuracy phase3 phase4 phase5 dct deen mini_deen)
    add_executable(neo_smo_${test}_tests tests/unit/${test}.cpp)
    target_include_directories(neo_smo_${test}_tests PRIVATE tests)
    target_link_libraries(neo_smo_${test}_tests PRIVATE neo_smo_core)
    target_compile_definitions(neo_smo_${test}_tests PRIVATE HWY_COMPILE_ALL_ATTAINABLE)
    if(MINGW)
      target_link_options(neo_smo_${test}_tests PRIVATE -static -static-libgcc -static-libstdc++)
    endif()
    if(test STREQUAL "deen")
      foreach(mode c2d c3d w2d w3d a2d a3d)
        add_test(NAME core.deen.${mode} COMMAND "${Python3_EXECUTABLE}"
          "${CMAKE_CURRENT_SOURCE_DIR}/tests/isolated.py" "$<TARGET_FILE:neo_smo_deen_tests>" ${mode})
        set_tests_properties(core.deen.${mode} PROPERTIES LABELS "spec" RUN_SERIAL TRUE TIMEOUT 90)
      endforeach()
    else()
      add_test(NAME core.${test} COMMAND "${Python3_EXECUTABLE}"
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/isolated.py" "$<TARGET_FILE:neo_smo_${test}_tests>")
      set_tests_properties(core.${test} PROPERTIES LABELS "spec" RUN_SERIAL TRUE TIMEOUT 90)
    endif()
  endforeach()
  # Exercise production arithmetic, but generate expectations and check results
  # in a separate translation unit without fast-math assumptions.
  target_sources(neo_smo_fp16_tests PRIVATE tests/unit/fp16_kernel.cpp)
  set_property(SOURCE tests/unit/fp16_kernel.cpp APPEND PROPERTY COMPILE_OPTIONS ${NEO_SMO_FP_OPTIONS})
  if(MSVC)
    set_property(SOURCE tests/unit/fp16.cpp APPEND PROPERTY COMPILE_OPTIONS /fp:precise)
  else()
    set_property(SOURCE tests/unit/fp16.cpp APPEND PROPERTY COMPILE_OPTIONS -fno-fast-math -ffp-contract=off)
  endif()
endif()
option(NEO_SMO_TEST_BLACKBOX "Compare against an external zsmooth binary" OFF)
set(NEO_SMO_REFERENCE_PLUGIN "${CMAKE_CURRENT_SOURCE_DIR}/zsmooth.dll" CACHE FILEPATH "Reference zsmooth plugin for differential tests")
if(BUILD_TESTING AND NEO_SMO_TEST_BLACKBOX AND NOT NEO_SMO_BUILD_VAPOURSYNTH)
  message(FATAL_ERROR "Black-box tests require NEO_SMO_BUILD_VAPOURSYNTH")
endif()

if(NEO_SMO_BUILD_AVISYNTH OR NEO_SMO_BUILD_VAPOURSYNTH)
  add_library(neo_smo MODULE)
  configure_file(src/plugin/version.hpp.in
    "${CMAKE_CURRENT_BINARY_DIR}/generated/neo_smo_version.hpp" @ONLY)
  target_include_directories(neo_smo PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/generated")
  if(WIN32)
    enable_language(RC)
    configure_file(src/plugin/version.rc.in
      "${CMAKE_CURRENT_BINARY_DIR}/generated/neo_smo_version.rc" @ONLY)
    target_sources(neo_smo PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/generated/neo_smo_version.rc")
  endif()
  if(NEO_SMO_BUILD_AVISYNTH)
    include(cmake/AviSynth.cmake)
    target_sources(neo_smo PRIVATE src/plugin/avisynth.cpp src/plugin/avs_filter.cpp src/plugin/avs_weighted.cpp src/plugin/avs_scene.cpp)
    target_include_directories(neo_smo SYSTEM PRIVATE "${NEO_SMO_AVS_INCLUDE}")
  endif()
  if(NEO_SMO_BUILD_VAPOURSYNTH)
    include(cmake/VapourSynth.cmake)
    target_sources(neo_smo PRIVATE src/plugin/vapoursynth.cpp src/plugin/dct_filter.cpp src/plugin/deen.cpp src/plugin/mini_deen.cpp)
    target_include_directories(neo_smo SYSTEM PRIVATE "${NEO_SMO_VS_INCLUDE}")
  endif()
  target_link_libraries(neo_smo PRIVATE neo_smo_core DualSynth::dualsynth)
  if(MSVC)
    target_compile_options(neo_smo PRIVATE /EHsc)
  endif()
  if(BUILD_TESTING AND NEO_SMO_BUILD_VAPOURSYNTH)
    add_executable(neo_smo_frame_errors_tests tests/unit/frame_errors.cpp src/plugin/dct_filter.cpp src/plugin/deen.cpp src/plugin/mini_deen.cpp)
    target_include_directories(neo_smo_frame_errors_tests PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/generated" "${NEO_SMO_VS_INCLUDE}")
    target_link_libraries(neo_smo_frame_errors_tests PRIVATE neo_smo_core)
    if(MSVC)
      target_compile_options(neo_smo_frame_errors_tests PRIVATE /EHsc)
    endif()
    if(MINGW)
      target_link_options(neo_smo_frame_errors_tests PRIVATE -static -static-libgcc -static-libstdc++)
    endif()
    add_test(NAME vs.frame_errors COMMAND "${Python3_EXECUTABLE}"
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/isolated.py" "$<TARGET_FILE:neo_smo_frame_errors_tests>")
    set_tests_properties(vs.frame_errors PROPERTIES LABELS "spec" RUN_SERIAL TRUE TIMEOUT 90)
    add_executable(neo_smo_deen_vs_errors_tests tests/unit/deen_vs_errors.cpp)
    target_include_directories(neo_smo_deen_vs_errors_tests PRIVATE "${NEO_SMO_VS_INCLUDE}")
    target_link_libraries(neo_smo_deen_vs_errors_tests PRIVATE neo_smo_core)
    if(MSVC)
      target_compile_options(neo_smo_deen_vs_errors_tests PRIVATE /EHsc)
    endif()
    if(MINGW)
      target_link_options(neo_smo_deen_vs_errors_tests PRIVATE -static -static-libgcc -static-libstdc++)
    endif()
    add_test(NAME vs.deen.errors COMMAND "${Python3_EXECUTABLE}"
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/isolated.py" "$<TARGET_FILE:neo_smo_deen_vs_errors_tests>")
    set_tests_properties(vs.deen.errors PROPERTIES LABELS "spec" RUN_SERIAL TRUE TIMEOUT 90)
    add_executable(neo_smo_mini_deen_vs_errors_tests tests/unit/mini_deen_vs_errors.cpp)
    target_include_directories(neo_smo_mini_deen_vs_errors_tests PRIVATE "${NEO_SMO_VS_INCLUDE}")
    target_link_libraries(neo_smo_mini_deen_vs_errors_tests PRIVATE neo_smo_core)
    if(MSVC)
      target_compile_options(neo_smo_mini_deen_vs_errors_tests PRIVATE /EHsc)
    endif()
    if(MINGW)
      target_link_options(neo_smo_mini_deen_vs_errors_tests PRIVATE -static -static-libgcc -static-libstdc++)
    endif()
    add_test(NAME vs.mini_deen.errors COMMAND "${Python3_EXECUTABLE}"
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/isolated.py" "$<TARGET_FILE:neo_smo_mini_deen_vs_errors_tests>")
    set_tests_properties(vs.mini_deen.errors PROPERTIES LABELS "spec" RUN_SERIAL TRUE TIMEOUT 90)
  endif()
  set_target_properties(neo_smo PROPERTIES
    PREFIX ""
    OUTPUT_NAME "neo-smo"
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON
  )
  if(MINGW)
    target_link_options(neo_smo PRIVATE -static -static-libgcc -static-libstdc++)
  endif()

  if(BUILD_TESTING AND NEO_SMO_TEST_BLACKBOX)
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    add_test(NAME vs.differential.dct COMMAND "${Python3_EXECUTABLE}"
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/differential/dct.py"
      --plugin "$<TARGET_FILE:neo_smo>"
      --reference "${NEO_SMO_REFERENCE_PLUGIN}")
    set_tests_properties(vs.differential.dct PROPERTIES LABELS "blackbox" RUN_SERIAL TRUE TIMEOUT 60)
    add_test(NAME vs.differential.phase1 COMMAND "${Python3_EXECUTABLE}"
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/differential/run.py"
      --plugin "$<TARGET_FILE:neo_smo>"
      --reference "${NEO_SMO_REFERENCE_PLUGIN}")
    set_tests_properties(vs.differential.phase1 PROPERTIES LABELS "blackbox" RUN_SERIAL TRUE TIMEOUT 300)
    add_test(NAME vs.differential.phase2 COMMAND "${Python3_EXECUTABLE}"
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/differential/phase2.py"
      --plugin "$<TARGET_FILE:neo_smo>"
      --reference "${NEO_SMO_REFERENCE_PLUGIN}")
    set_tests_properties(vs.differential.phase2 PROPERTIES LABELS "blackbox" RUN_SERIAL TRUE TIMEOUT 300)
    add_test(NAME vs.differential.phase3 COMMAND "${Python3_EXECUTABLE}"
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/differential/phase3.py"
      --plugin "$<TARGET_FILE:neo_smo>"
      --reference "${NEO_SMO_REFERENCE_PLUGIN}")
    set_tests_properties(vs.differential.phase3 PROPERTIES LABELS "blackbox" RUN_SERIAL TRUE TIMEOUT 300)
    add_test(NAME vs.differential.phase4 COMMAND "${Python3_EXECUTABLE}"
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/differential/phase4.py"
      --plugin "$<TARGET_FILE:neo_smo>"
      --reference "${NEO_SMO_REFERENCE_PLUGIN}")
    set_tests_properties(vs.differential.phase4 PROPERTIES LABELS "blackbox" RUN_SERIAL TRUE TIMEOUT 300)
    add_test(NAME vs.differential.phase5 COMMAND "${Python3_EXECUTABLE}"
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/differential/phase5.py"
      --neo-smo-dll "$<TARGET_FILE:neo_smo>"
      --zsmooth-dll "${NEO_SMO_REFERENCE_PLUGIN}")
    set_tests_properties(vs.differential.phase5 PROPERTIES LABELS "blackbox" RUN_SERIAL TRUE TIMEOUT 300)
    add_test(NAME vs.differential.regressions COMMAND "${Python3_EXECUTABLE}"
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/differential/regressions.py"
      --plugin "$<TARGET_FILE:neo_smo>" --reference "${NEO_SMO_REFERENCE_PLUGIN}")
    set_tests_properties(vs.differential.regressions PROPERTIES LABELS "blackbox" RUN_SERIAL TRUE TIMEOUT 90)
  endif()
endif()

option(NEO_SMO_TEST_VAPOURSYNTH "Run VapourSynth host smoke tests" OFF)
if(BUILD_TESTING AND NEO_SMO_TEST_VAPOURSYNTH)
  if(NOT NEO_SMO_BUILD_VAPOURSYNTH)
    message(FATAL_ERROR "VapourSynth tests require NEO_SMO_BUILD_VAPOURSYNTH")
  endif()
  add_test(NAME vs.smoke COMMAND "${Python3_EXECUTABLE}"
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/isolated.py" "${Python3_EXECUTABLE}"
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/vapoursynth/smoke.py" --plugin "$<TARGET_FILE:neo_smo>")
  set_tests_properties(vs.smoke PROPERTIES LABELS "spec;vapoursynth" RUN_SERIAL TRUE TIMEOUT 90)
  add_test(NAME vs.mini_deen COMMAND "${Python3_EXECUTABLE}"
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/isolated.py" "${Python3_EXECUTABLE}"
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/vapoursynth/mini_deen.py" --plugin "$<TARGET_FILE:neo_smo>")
  set_tests_properties(vs.mini_deen PROPERTIES LABELS "spec;vapoursynth" RUN_SERIAL TRUE TIMEOUT 90)
  foreach(mode c2d c3d w2d w3d a2d a3d)
    add_test(NAME vs.deen.${mode} COMMAND "${Python3_EXECUTABLE}"
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/isolated.py" "${Python3_EXECUTABLE}"
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/vapoursynth/deen.py" --plugin "$<TARGET_FILE:neo_smo>" --mode ${mode})
    set_tests_properties(vs.deen.${mode} PROPERTIES LABELS "spec;vapoursynth" RUN_SERIAL TRUE TIMEOUT 90)
  endforeach()
endif()

option(NEO_SMO_TEST_AVISYNTH "Run AviSynth host acceptance tests" OFF)
option(NEO_SMO_TEST_CROSS_HOST "Compare AviSynth and VapourSynth adapters on identical input" OFF)
if(BUILD_TESTING AND NEO_SMO_TEST_AVISYNTH)
  if(NOT NEO_SMO_BUILD_AVISYNTH)
    message(FATAL_ERROR "AviSynth tests require NEO_SMO_BUILD_AVISYNTH")
  endif()
  set(NEO_SMO_AVISYNTH_RUNTIME "" CACHE FILEPATH "AviSynth runtime library for host tests")
  if(NOT EXISTS "${NEO_SMO_AVISYNTH_RUNTIME}")
    message(FATAL_ERROR "Set NEO_SMO_AVISYNTH_RUNTIME to the matching AviSynth runtime library")
  endif()
  add_executable(neo_smo_avs_runner "${dualsynth2_SOURCE_DIR}/tests/tools/avisynth_runner.cpp")
  target_include_directories(neo_smo_avs_runner SYSTEM PRIVATE "${NEO_SMO_AVS_INCLUDE}")
  target_link_libraries(neo_smo_avs_runner PRIVATE ${CMAKE_DL_LIBS})
  set(_neo_smo_avs_host_api c)
  if(MSVC)
    set(_neo_smo_avs_host_api cpp)
  endif()
  add_test(NAME avs.acceptance COMMAND "${Python3_EXECUTABLE}"
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/avisynth/acceptance.py"
    --runner "$<TARGET_FILE:neo_smo_avs_runner>" --plugin "$<TARGET_FILE:neo_smo>"
    --runtime "${NEO_SMO_AVISYNTH_RUNTIME}" --host-api "${_neo_smo_avs_host_api}"
    --work "${CMAKE_CURRENT_BINARY_DIR}/avisynth-tests")
  set_tests_properties(avs.acceptance PROPERTIES LABELS "spec;avisynth" RUN_SERIAL TRUE TIMEOUT 180)
  unset(_neo_smo_avs_host_api)
endif()

if(BUILD_TESTING AND NEO_SMO_TEST_CROSS_HOST)
  if(NOT NEO_SMO_BUILD_AVISYNTH OR NOT NEO_SMO_BUILD_VAPOURSYNTH)
    message(FATAL_ERROR "Cross-host tests require both host interfaces")
  endif()
  set(NEO_SMO_AVISYNTH_RUNTIME "" CACHE FILEPATH "AviSynth runtime library for host tests")
  if(NOT EXISTS "${NEO_SMO_AVISYNTH_RUNTIME}")
    message(FATAL_ERROR "Set NEO_SMO_AVISYNTH_RUNTIME to the matching AviSynth runtime library")
  endif()
  add_test(NAME avs.cross_host COMMAND "${Python3_EXECUTABLE}"
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/avisynth/cross_host.py"
    --plugin "$<TARGET_FILE:neo_smo>" --runtime "${NEO_SMO_AVISYNTH_RUNTIME}")
  set_tests_properties(avs.cross_host PROPERTIES LABELS "spec;avisynth" RUN_SERIAL TRUE TIMEOUT 180)
endif()
