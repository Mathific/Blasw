# BSD 3-Clause License
#
# Copyright (c) 2021, Shahriar Rezghi <shahriar25.ss@gmail.com>
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Project
cmake_minimum_required(VERSION 3.5)
project(Blasw VERSION 0.1.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Variables
set(BLASW_MASTER OFF CACHE INTERNAL "")
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(BLASW_MASTER ON CACHE INTERNAL "")
endif()

option(BLASW_INSTALL "Install BLASW Library" ${BLASW_MASTER})
option(BLASW_TESTS "Build BLASW Tests" ${BLASW_MASTER})

# Files
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/src/blasw/config.h.in"
    "${CMAKE_CURRENT_BINARY_DIR}/gen/blasw/config.h")
file(GLOB HEADERS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/blasw/blasw.h"
    "${CMAKE_CURRENT_BINARY_DIR}/gen/blasw/config.h")
list(APPEND CMAKE_MODULE_PATH
    "${CMAKE_CURRENT_LIST_DIR}/cmake/")

# Library
list(APPEND CMAKE_MODULE_PATH
    "${CMAKE_CURRENT_LIST_DIR}/cmake")
find_package(CBLAS REQUIRED)
find_package(LAPACKE)

if(NOT LAPACKE_FOUND)
    add_library(LAPACKE::LAPACKE INTERFACE IMPORTED)
endif()

add_library(BLASW INTERFACE)
target_link_libraries(BLASW INTERFACE
    CBLAS::CBLAS LAPACKE::LAPACKE)
target_include_directories(BLASW INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/>"
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/gen/>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include/>")
add_library(BLASW::BLASW ALIAS BLASW)

if(BLASW_TESTS)
    add_executable(BlaswTests test/tests.cpp ${HEADERS})
    target_link_libraries(BlaswTests BLASW::BLASW)
endif()

# Installation
if(BLASW_INSTALL)
    include(CMakePackageConfigHelpers)
    write_basic_package_version_file(
        "BLASWConfigVersion.cmake"
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion)

    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/BLASWConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/BLASWConfig.cmake" @ONLY)
    file(GLOB MODULES
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindCBLAS.cmake"
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindLAPACKE.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/BLASWConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/BLASWConfigVersion.cmake")
    install(FILES ${HEADERS}
        DESTINATION include/blasw)
    install(FILES ${MODULES}
        DESTINATION lib/cmake/blasw)

    install(TARGETS BLASW
        EXPORT BLASWTargets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include)
    install(EXPORT BLASWTargets
        NAMESPACE BLASW::
        DESTINATION lib/cmake/blasw)
    export(TARGETS BLASW
        NAMESPACE BLASW::
        FILE "${CMAKE_CURRENT_BINARY_DIR}/BLASWTargets.cmake")
endif()
