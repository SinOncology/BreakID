cmake_minimum_required(VERSION 2.8.8)
set(PROJECT_NAME_STR BreakID)
project(${PROJECT_NAME_STR} C CXX)
# external dependencies to be discovered:
set(SAMTOOLS_DIR ${PROJECT_SOURCE_DIR}/thirdparty/samtools/samtools-1.3.1)
set(HTSLIB_DIR ${SAMTOOLS_DIR}/htslib-1.3.1)
set(GTEST_DIR ${PROJECT_SOURCE_DIR}/thirdparty/googletest/googletest-1.8.0)

#set(SAMTOOLS_DIR "/mnt/pipeline-programs/samtools/samtools-1.2")
#set(HTSLIB_DIR "/mnt/pipeline-programs/samtools/samtools-1.2/htslib-1.2.1")
#set(GTEST_DIR "/mnt/pipeline-programs/googletest/googletest-1.8.0")

if (APPLE)
    # For Macs, override the paths here
else ()
    # to compile in aliyun VM using g++-6.3.0
    include_directories(/mnt/sata3_storage/pipe_app/gcc-6.5.0/bin/include)
    set(CMAKE_CXX_COMPILER "/mnt/sata3_storage/pipe_app/gcc-6.5.0/bin/g++")
    set(CMAKE_C_COMPILER "/mnt/sata3_storage/pipe_app/gcc-6.5.0/bin/gcc")

    #    include_directories(/usr/include/c++/4.8.5)
    #
    #    set(CMAKE_CXX_COMPILER "/usr/bin/g++")
    #    set(CMAKE_C_COMPILER "/usr/bin/gcc")

    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath -Wl,/usr/local/lib64")

endif ()

# no version requirement at the moment
find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)
find_package(GTest QUIET)
# Gtest might not be system-available
if (NOT GTest_FOUND)
    message("Didn't find GTest system-wise, checking custom location...")
    find_library(GTEST_MAIN NAMES gtest PATHS ${GTEST_DIR})
    find_library(GTEST_MAIN_MAIN NAMES gtest_main PATHS ${GTEST_DIR} NO_DEFAULT_PATH)
    set(GTEST_INCLUDE_DIR ${GTEST_DIR}/include NO_DEFAULT_PATH)
    set(GTEST_BOTH_LIBRARIES "${GTEST_MAIN};${GTEST_MAIN_MAIN}")
    if (NOT GTEST_MAIN)
        message("No gtest library found in ${GTEST_DIR}")
    endif ()
    if (NOT GTEST_MAIN_MAIN)
        message("No gtest_main library found in ${GTEST_DIR}")
    endif ()
endif ()
include_directories(${GTEST_INCLUDE_DIR})
# Discover samtools and htslib
find_library(SAMTOOLS_BAM_LIB NAMES bam HINTS ${SAMTOOLS_DIR} STATIC IMPORTED)
find_library(HTS_LIB NAMES hts HINTS ${HTSLIB_DIR} STATIC IMPORTED)

# Compiler flags
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    # using Clang
    add_definitions(-Wall -g -O2 -std=c++11)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # using GCC
    add_definitions(-Wall -g -O2  -g -std=c++11)
    #add_definitions(-Wall -g -std=c++11 -fprofile-arcs -ftest-coverage)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    # using Intel C++
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    # using Visual Studio C++
endif ()

# Create the installdir.h include file
set(REFERENCE_DATA_MAKER "./make_ref_data.sh")
#  Set up testing
enable_testing()
# Compile CaGe-A and the unit-tests
add_subdirectory(${PROJECT_SOURCE_DIR}/src)
