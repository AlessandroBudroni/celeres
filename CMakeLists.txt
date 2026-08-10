cmake_minimum_required(VERSION 3.16)
project(celeres C ASM)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

set(IMPL "ref" CACHE STRING "Implementation to build: ref or avx2")
set_property(CACHE IMPL PROPERTY STRINGS ref avx2)
set(RING_SIZE "2" CACHE STRING "Static ring size compiled into the library")
set(BUILD_TESTS ON CACHE BOOL "Build protocol tests")
set(BUILD_BENCHMARKS ON CACHE BOOL "Build benchmark programs")
set(BENCH_ITERS "100" CACHE STRING "Iterations used by benchmark programs")
set(CELERES_TAG_MODE "self_orthogonal" CACHE STRING "CELERES tag encoding: self_orthogonal or rref")
set_property(CACHE CELERES_TAG_MODE PROPERTY STRINGS self_orthogonal rref)

string(TOLOWER "${IMPL}" IMPL)
string(TOLOWER "${CELERES_TAG_MODE}" CELERES_TAG_MODE)

if(NOT IMPL STREQUAL "ref" AND NOT IMPL STREQUAL "avx2")
    message(FATAL_ERROR "IMPL must be either 'ref' or 'avx2'")
endif()
if(NOT CELERES_TAG_MODE STREQUAL "self_orthogonal" AND
   NOT CELERES_TAG_MODE STREQUAL "rref")
    message(FATAL_ERROR "CELERES_TAG_MODE must be either 'self_orthogonal' or 'rref'")
endif()
if(RING_SIZE LESS 1)
    message(FATAL_ERROR "RING_SIZE must be at least 1")
endif()
if(IMPL STREQUAL "avx2" AND NOT CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
    message(FATAL_ERROR "IMPL=avx2 requires an x86_64/AMD64 target")
endif()

set(CELERES_COMMON_SOURCES
    src/celeres.c
    src/celeres_known_is.c
    src/celeres_util.c
    src/common/celeres_tag.c
    src/common/monomial.c
    src/common/rng.c
    src/common/seedtree.c)

set(CELERES_REF_SOURCES
    ${CELERES_COMMON_SOURCES}
    src/ref/celeres_internal.c
    src/ref/codes.c
    src/ref/fips202.c
    src/ref/keccakf1600.c
    src/ref/permutation_compress.c
    src/ref/utils.c)

set(CELERES_AVX2_SOURCES
    ${CELERES_COMMON_SOURCES}
    src/avx2/celeres_internal.c
    src/avx2/codes.c
    src/avx2/fips202.c
    src/avx2/keccakf1600.c
    src/avx2/permutation_compress.c
    src/avx2/utils.c
    src/avx2/KeccakP-1600-AVX2.s)

if(APPLE)
    set_source_files_properties(src/avx2/KeccakP-1600-AVX2.s PROPERTIES
        COMPILE_OPTIONS "-Wa,-defsym,no_type=1;-Wa,-defsym,no_size=1;-Wa,-defsym,no_plt=1")
endif()

message(STATUS "Implementation: ${IMPL}")
message(STATUS "Ring size: ${RING_SIZE}")
message(STATUS "CELERES tag mode: ${CELERES_TAG_MODE}")

function(add_celeres_level level category target)
    set(lib celeres_l${level})
    if(IMPL STREQUAL "ref")
        add_library(${lib} STATIC ${CELERES_REF_SOURCES})
        target_include_directories(${lib} PUBLIC src src/common src/ref)
    else()
        add_library(${lib} STATIC ${CELERES_AVX2_SOURCES})
        target_include_directories(${lib} PUBLIC src src/common src/avx2)
        target_compile_definitions(${lib} PUBLIC USE_AVX2)
        target_compile_options(${lib} PRIVATE -mavx2 -mavx -mbmi2 -mbmi)
    endif()

    target_compile_definitions(${lib} PUBLIC
        CATEGORY=${category}
        TARGET=${target}
        CELERES_RING_SIZE=${RING_SIZE})
    if(CELERES_TAG_MODE STREQUAL "self_orthogonal")
        target_compile_definitions(${lib} PUBLIC
            CELERES_TAG_SELF_ORTHOGONAL=1
            CELERES_TAG_MODE_NAME="self_orthogonal"
            CELERES_TAG_COMPRESSION_NAME="diagonal")
    else()
        target_compile_definitions(${lib} PUBLIC
            CELERES_TAG_SELF_ORTHOGONAL=0
            CELERES_TAG_MODE_NAME="rref"
            CELERES_TAG_COMPRESSION_NAME="rref")
    endif()
    target_compile_options(${lib} PRIVATE -Wall -Wextra -Wpedantic -Wno-type-limits)

    if(BUILD_TESTS)
        foreach(test_base IN ITEMS test_celeres test_celeres_link)
            set(test_name ${test_base}_l${level})
            add_executable(${test_name} tests/${test_base}.c)
            target_link_libraries(${test_name} PRIVATE ${lib})
            add_test(NAME ${test_name} COMMAND ${test_name})
        endforeach()
    endif()

    if(BUILD_BENCHMARKS)
        set(bench_name bench_celeres_l${level})
        add_executable(${bench_name} tests/test_celeres_benchmark.c)
        target_link_libraries(${bench_name} PRIVATE ${lib})
        target_compile_definitions(${bench_name} PRIVATE
            IMPL_NAME="${IMPL}"
            BENCH_ITERS=${BENCH_ITERS})
    endif()
endfunction()

enable_testing()

add_celeres_level(1 272 192)
add_celeres_level(3 416 220)
add_celeres_level(5 564 345)
