cmake_minimum_required(VERSION 3.16)
project(MinimalResolution CXX)

# ---------------------------------------------------------------------------
# This project was originally built via ten standalone shell scripts at the
# repo root (st_compiling, BPtable_compile, BP_compile, ex_compile,
# mot_compile, mot_comb_compile, mottable_compile, kos_compile,
# tau_mult_compile, tauboc_compile, e2p_compile), each a raw g++ invocation.
# Those scripts are left in place and still work; this CMakeLists.txt builds
# the same ten executables in one reusable, IDE-friendly build.
#
# Toolchain notes (see also README note near the shell scripts):
#   - This build standardizes on AppleClang (the system "g++"/"clang++") for
#     every target, plus Homebrew's libomp for OpenMP support, rather than
#     mixing in Homebrew GCC. AppleClang does not support -fopenmp directly;
#     it needs "-Xpreprocessor -fopenmp" plus libomp's include/lib paths.
#   - Two targets (BPtab, mr_BP) link GMP (gmp/gmpxx) via Qp.h/Qp.cpp.
#     Homebrew's libgmpxx is built against libc++, matching AppleClang's
#     standard library, so AppleClang + Homebrew libomp + Homebrew gmp link
#     together cleanly. Mixing in a libstdc++ toolchain (e.g. Homebrew GCC)
#     would hit ABI-mismatch link errors against Homebrew's gmpxx.
# ---------------------------------------------------------------------------

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# All original scripts use -O2 regardless of CMAKE_BUILD_TYPE; keep that.
add_compile_options(-O2)

# ---------------------------------------------------------------------------
# Optional: widen the packed `exponent` type (polynomial/2.h) from 32 to 64
# bits, to support higher Steenrod-algebra generators/exponents (e.g. h_n
# classes at n>7). Off by default to keep on-disk file formats (mot_gens,
# mot_res, poly_exponents, etc.) byte-identical to all prior output. Files
# from a WIDE_EXPONENT=ON build are NOT binary-compatible with a default
# build; use a separate build directory (e.g. `cmake -B build-wide
# -DWIDE_EXPONENT=ON`) and run a full pipeline end-to-end from binaries of
# the same configuration only.
# ---------------------------------------------------------------------------
option(WIDE_EXPONENT "Use a 64-bit packed exponent type instead of the default 32-bit" OFF)

# Every target whose sources include exponents.cpp and/or ex_exponents.cpp
# needs EXPONENT_WIDTH applied; populated via list(APPEND ...) right after
# each add_executable() below so it can't drift out of sync.
set(EXPONENT_WIDTH_TARGETS "")

# ---------------------------------------------------------------------------
# OpenMP via Homebrew libomp (AppleClang has no built-in OpenMP runtime).
# find_package(OpenMP) does not locate Homebrew's libomp out of the box on
# AppleClang, so we look it up explicitly via `brew --prefix libomp`, with a
# hardcoded fallback to the standard Apple Silicon prefix.
# ---------------------------------------------------------------------------
find_program(BREW_EXECUTABLE brew)
if(BREW_EXECUTABLE)
  execute_process(
    COMMAND ${BREW_EXECUTABLE} --prefix libomp
    OUTPUT_VARIABLE LIBOMP_PREFIX
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
endif()
if(NOT LIBOMP_PREFIX)
  set(LIBOMP_PREFIX /opt/homebrew/opt/libomp)
endif()

find_path(LIBOMP_INCLUDE_DIR omp.h PATHS ${LIBOMP_PREFIX}/include)
find_library(LIBOMP_LIBRARY omp PATHS ${LIBOMP_PREFIX}/lib)

if(NOT LIBOMP_INCLUDE_DIR OR NOT LIBOMP_LIBRARY)
  message(FATAL_ERROR "Could not find Homebrew libomp (looked in ${LIBOMP_PREFIX}). Install it with `brew install libomp`.")
endif()

# Imported target standing in for OpenMP::OpenMP_CXX, since
# find_package(OpenMP) does not reliably resolve Homebrew's libomp with
# AppleClang.
add_library(homebrew_omp INTERFACE)
target_compile_options(homebrew_omp INTERFACE -Xpreprocessor -fopenmp)
target_include_directories(homebrew_omp INTERFACE ${LIBOMP_INCLUDE_DIR})
target_link_libraries(homebrew_omp INTERFACE ${LIBOMP_LIBRARY})

# ---------------------------------------------------------------------------
# GMP via Homebrew (only needed by BPtab and mr_BP, through Qp.h/Qp.cpp).
# CMake ships no FindGMP module, so locate it manually under the Homebrew
# prefix.
# ---------------------------------------------------------------------------
if(BREW_EXECUTABLE)
  execute_process(
    COMMAND ${BREW_EXECUTABLE} --prefix gmp
    OUTPUT_VARIABLE GMP_PREFIX
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
endif()
if(NOT GMP_PREFIX)
  set(GMP_PREFIX /opt/homebrew/opt/gmp)
endif()

find_path(GMP_INCLUDE_DIR gmp.h PATHS ${GMP_PREFIX}/include)
find_library(GMP_LIBRARY gmp PATHS ${GMP_PREFIX}/lib)
find_library(GMPXX_LIBRARY gmpxx PATHS ${GMP_PREFIX}/lib)

if(NOT GMP_INCLUDE_DIR OR NOT GMP_LIBRARY OR NOT GMPXX_LIBRARY)
  message(FATAL_ERROR "Could not find Homebrew gmp/gmpxx (looked in ${GMP_PREFIX}). Install it with `brew install gmp`.")
endif()

add_library(homebrew_gmp INTERFACE)
target_include_directories(homebrew_gmp INTERFACE ${GMP_INCLUDE_DIR})
target_link_libraries(homebrew_gmp INTERFACE ${GMPXX_LIBRARY} ${GMP_LIBRARY})

# ---------------------------------------------------------------------------
# Shared include path: every original script passes -I./ (the repo root).
# ---------------------------------------------------------------------------
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Helper: every script also passes -Wall -Wfatal-errors (kos_compile omits
# -Wfatal-errors but it's harmless to add uniformly).
add_compile_options(-Wall -Wfatal-errors)

# ---------------------------------------------------------------------------
# Targets. Source lists, standards, and output names mirror each shell
# script exactly.
# ---------------------------------------------------------------------------

# st_compiling -> mr_st (-std=c++11)
add_executable(mr_st
  exponents.cpp Fp.cpp mon_index.cpp steenrod.cpp steenrod_init.cpp stmain.cpp
)
set_target_properties(mr_st PROPERTIES CXX_STANDARD 11)
target_link_libraries(mr_st PRIVATE homebrew_omp)
list(APPEND EXPONENT_WIDTH_TARGETS mr_st)

# BPtable_compile -> BPtab (-std=c++11, links gmp/gmpxx)
add_executable(BPtab
  Z2.cpp mon_index.cpp BP.cpp BPQ.cpp BPtable.cpp exponents.cpp Qp.cpp Fp.cpp
)
set_target_properties(BPtab PROPERTIES CXX_STANDARD 11)
target_link_libraries(BPtab PRIVATE homebrew_omp homebrew_gmp)
# Always 32-bit: BP.cpp's etaL/R2L relies on exponent and matrix_index being
# the literal same type, which only holds in the 32-bit build.
target_compile_definitions(BPtab PRIVATE EXPONENT_WIDTH=32)

# BP_compile -> mr_BP (-std=c++11, links gmp/gmpxx via Qp.h)
add_executable(mr_BP
  BPcomplex.cpp streams.cpp algNov.cpp Boc.cpp multiplication.cpp exponents.cpp
  Fp.cpp mon_index.cpp Z2.cpp BP.cpp BP_init.cpp BPmain.cpp
)
set_target_properties(mr_BP PROPERTIES CXX_STANDARD 11)
target_link_libraries(mr_BP PRIVATE homebrew_omp homebrew_gmp)
# Always 32-bit: same BP.cpp constraint as BPtab above.
target_compile_definitions(mr_BP PRIVATE EXPONENT_WIDTH=32)

# ex_compile -> mr_ex (-std=c++14)
add_executable(mr_ex
  ex_steenrod_init.cpp ex_exponents.cpp ex_index.cpp ex_main.cpp ctau_steenrod.cpp
  Fp.cpp streams.cpp
)
set_target_properties(mr_ex PROPERTIES CXX_STANDARD 14)
target_link_libraries(mr_ex PRIVATE homebrew_omp)
list(APPEND EXPONENT_WIDTH_TARGETS mr_ex)

# mot_compile -> mr_mot (no -std specified upstream; use c++14, verified to work)
add_executable(mr_mot
  Fp.cpp streams.cpp mot_main.cpp mot_steenrod.cpp exponents.cpp
)
set_target_properties(mr_mot PROPERTIES CXX_STANDARD 14)
target_link_libraries(mr_mot PRIVATE homebrew_omp)
list(APPEND EXPONENT_WIDTH_TARGETS mr_mot)

# mot_comb_compile -> mot_comb (no -std specified upstream; use c++14)
add_executable(mot_comb
  Fp.cpp streams.cpp mot_combine.cpp mot_steenrod.cpp exponents.cpp
)
set_target_properties(mot_comb PROPERTIES CXX_STANDARD 14)
target_link_libraries(mot_comb PRIVATE homebrew_omp)
list(APPEND EXPONENT_WIDTH_TARGETS mot_comb)

# mottable_compile -> motTab (-std=c++11)
add_executable(motTab
  exponents.cpp mot_steenrod.cpp mot_coactions.cpp mon_index.cpp
)
set_target_properties(motTab PROPERTIES CXX_STANDARD 11)
target_link_libraries(motTab PRIVATE homebrew_omp)
list(APPEND EXPONENT_WIDTH_TARGETS motTab)

# kos_compile -> kos (-std=c++11)
add_executable(kos
  steenrod.cpp steenrod_init.cpp kosul.cpp mon_index.cpp Fp.cpp exponents.cpp
)
set_target_properties(kos PROPERTIES CXX_STANDARD 11)
target_link_libraries(kos PRIVATE homebrew_omp)
list(APPEND EXPONENT_WIDTH_TARGETS kos)

# tau_mult_compile -> mot_mult (-std=c++11)
add_executable(mot_mult
  mot_mult.cpp tao_bockstein.cpp mot_steenrod.cpp exponents.cpp
)
set_target_properties(mot_mult PROPERTIES CXX_STANDARD 11)
target_link_libraries(mot_mult PRIVATE homebrew_omp)
list(APPEND EXPONENT_WIDTH_TARGETS mot_mult)

# tauboc_compile -> tauBoc (-std=c++11)
add_executable(tauBoc
  taubocmain.cpp tao_bockstein.cpp mot_steenrod.cpp exponents.cpp
)
set_target_properties(tauBoc PROPERTIES CXX_STANDARD 11)
target_link_libraries(tauBoc PRIVATE homebrew_omp)
list(APPEND EXPONENT_WIDTH_TARGETS tauBoc)

# e2p_compile -> e2p (-std=c++11)
add_executable(e2p
  ex_exponents.cpp ex_index.cpp ex2poly.cpp
)
set_target_properties(e2p PROPERTIES CXX_STANDARD 11)
target_link_libraries(e2p PRIVATE homebrew_omp)
list(APPEND EXPONENT_WIDTH_TARGETS e2p)

# ---------------------------------------------------------------------------
# Apply the exponent-width compile definition to every target that touches
# exponents.cpp / ex_exponents.cpp.
# ---------------------------------------------------------------------------
if(WIDE_EXPONENT)
  set(EXPONENT_WIDTH_VALUE 64)
else()
  set(EXPONENT_WIDTH_VALUE 32)
endif()
foreach(tgt IN LISTS EXPONENT_WIDTH_TARGETS)
  target_compile_definitions(${tgt} PRIVATE EXPONENT_WIDTH=${EXPONENT_WIDTH_VALUE})
endforeach()
