cmake_minimum_required(VERSION 3.24)
project(neo_smo VERSION 0.1.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(NEO_SMO_BUILD_VAPOURSYNTH "Build the VapourSynth v4 plugin" ON)

include(cmake/Dependencies.cmake)

add_library(neo_smo_core STATIC
  src/algorithms/plan.cpp
  src/kernels/median.cpp
  src/kernels/median_3x3.cpp
  src/kernels/median_5x5.cpp
  src/kernels/median_7x7.cpp
  src/kernels/vertical_cleaner.cpp
  src/kernels/remove_grain.cpp
)
target_compile_features(neo_smo_core PUBLIC cxx_std_17)
target_include_directories(neo_smo_core PUBLIC src)
target_include_directories(neo_smo_core SYSTEM PUBLIC "${dualsynth2_SOURCE_DIR}/include")
target_link_libraries(neo_smo_core PUBLIC hwy)
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

# Limit expansion of the EMU128 sorting network in MSVC. Hardware SIMD targets
# retain normal inlining and optimization.
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  set_property(SOURCE src/kernels/median_7x7.cpp APPEND PROPERTY COMPILE_DEFINITIONS
    "$<$<CONFIG:Release,RelWithDebInfo>:NEO_SMO_OUTLINE_EMU_COMPARE_SWAP>")
endif()

# EMU128 half rounding expands into large scalar control flow in MSVC/clang-cl.
# Keep it as a function boundary without restricting hardware SIMD inlining.
if(MSVC)
  set_property(SOURCE src/kernels/remove_grain.cpp APPEND PROPERTY COMPILE_DEFINITIONS
    "$<$<CONFIG:Release,RelWithDebInfo>:NEO_SMO_OUTLINE_EMU_ROUND_F16>")
endif()

include(CTest)
if(BUILD_TESTING)
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
  foreach(test boundaries fp16 float_accuracy)
    add_executable(neo_smo_${test}_tests tests/unit/${test}.cpp)
    target_include_directories(neo_smo_${test}_tests PRIVATE tests)
    target_link_libraries(neo_smo_${test}_tests PRIVATE neo_smo_core)
    target_compile_definitions(neo_smo_${test}_tests PRIVATE HWY_COMPILE_ALL_ATTAINABLE)
    if(MINGW)
      target_link_options(neo_smo_${test}_tests PRIVATE -static -static-libgcc -static-libstdc++)
    endif()
    add_test(NAME core.${test} COMMAND "${Python3_EXECUTABLE}"
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/isolated.py" "$<TARGET_FILE:neo_smo_${test}_tests>")
    set_tests_properties(core.${test} PROPERTIES LABELS "spec" RUN_SERIAL TRUE TIMEOUT 90)
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

if(NEO_SMO_BUILD_VAPOURSYNTH)
  add_library(neo_smo MODULE)
  configure_file(src/plugin/version.hpp.in
    "${CMAKE_CURRENT_BINARY_DIR}/generated/neo_smo_version.hpp" @ONLY)
  target_include_directories(neo_smo PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/generated")
  include(cmake/VapourSynth.cmake)
  target_sources(neo_smo PRIVATE src/plugin/vapoursynth.cpp)
  target_include_directories(neo_smo SYSTEM PRIVATE "${NEO_SMO_VS_INCLUDE}")
  target_link_libraries(neo_smo PRIVATE neo_smo_core DualSynth::dualsynth)
  if(BUILD_TESTING)
    add_executable(neo_smo_frame_errors_tests tests/unit/frame_errors.cpp)
    target_include_directories(neo_smo_frame_errors_tests PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/generated" "${NEO_SMO_VS_INCLUDE}")
    target_link_libraries(neo_smo_frame_errors_tests PRIVATE neo_smo_core)
    if(MINGW)
      target_link_options(neo_smo_frame_errors_tests PRIVATE -static -static-libgcc -static-libstdc++)
    endif()
    add_test(NAME vs.frame_errors COMMAND "${Python3_EXECUTABLE}"
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/isolated.py" "$<TARGET_FILE:neo_smo_frame_errors_tests>")
    set_tests_properties(vs.frame_errors PROPERTIES LABELS "spec" RUN_SERIAL TRUE TIMEOUT 90)
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
    add_test(NAME vs.differential.phase1 COMMAND "${Python3_EXECUTABLE}"
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/differential/run.py"
      --plugin "$<TARGET_FILE:neo_smo>"
      --reference "${NEO_SMO_REFERENCE_PLUGIN}")
    set_tests_properties(vs.differential.phase1 PROPERTIES LABELS "blackbox" RUN_SERIAL TRUE TIMEOUT 300)
    add_test(NAME vs.differential.regressions COMMAND "${Python3_EXECUTABLE}"
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/differential/regressions.py"
      --plugin "$<TARGET_FILE:neo_smo>" --reference "${NEO_SMO_REFERENCE_PLUGIN}")
    set_tests_properties(vs.differential.regressions PROPERTIES LABELS "blackbox" RUN_SERIAL TRUE TIMEOUT 90)
  endif()
endif()
