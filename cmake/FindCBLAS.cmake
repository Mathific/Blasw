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

set(CBLAS_MKL OFF CACHE INTERNAL "")
find_package(BLAS)
unset(CBLAS_LINKER_FLAGS CACHE)

if(BLAS_FOUND)
    include(CheckFunctionExists)
    set(CMAKE_REQUIRED_LIBRARIES ${BLAS_LIBRARIES})
    check_function_exists(cblas_saxpy TEMP_FOUND)

    if(TEMP_FOUND)
        set(CBLAS_LIBRARY ${BLAS_LIBRARIES})
        set(CBLAS_LINKER_FLAGS "")

        foreach(TEMP_NAME ${BLAS_LIBRARIES})
            get_filename_component(TEMP_NAME "${TEMP_NAME}" NAME)
            if(TEMP_NAME MATCHES "libmkl.*.so|libmkl.*.a")
                set(CBLAS_MKL ON CACHE INTERNAL "")
            endif()
        endforeach()
    endif()

    unset(TEMP_FOUND CACHE)
endif()

if(CBLAS_MKL)
    find_path(CBLAS_INCLUDE_DIR
        NAMES mkl_cblas.h
        PATHS /opt/intel/mkl/
        PATH_SUFFIXES include)
endif()

if((NOT CBLAS_INCLUDE_DIR) OR
        (CBLAS_INCLUDE_DIR STREQUAL "CBLAS_INCLUDE_DIR-NOTFOUND"))
    find_path(CBLAS_INCLUDE_DIR
        NAMES cblas.h
        PATHS /usr/include/
        PATH_SUFFIXES)
endif()

if(NOT CBLAS_LIBRARY)
    find_library(CBLAS_LIBRARY
        NAMES cblas
        PATHS ${CBLAS_INCLUDE_DIR}/../
        PATH_SUFFIXES lib lib64)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
    CBLAS
    FOUND_VAR CBLAS_FOUND
    REQUIRED_VARS
    CBLAS_LIBRARY
    CBLAS_INCLUDE_DIR)

if(CBLAS_FOUND)
    set(CBLAS_LIBRARIES ${CBLAS_LIBRARY})
    set(CBLAS_INCLUDE_DIRS ${CBLAS_INCLUDE_DIR})
endif()

if(CBLAS_FOUND AND NOT TARGET CBLAS::CBLAS)
    add_library(CBLAS::CBLAS INTERFACE IMPORTED)
    set_target_properties(CBLAS::CBLAS PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${CBLAS_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES "${CBLAS_LIBRARIES}"
        INTERFACE_LINK_OPTIONS "${CBLAS_LINKER_FLAGS}"
        INTERFACE_COMPILE_DEFINITIONS BLASW_CBLAS_FOUND)

    if(CBLAS_MKL)
        target_compile_definitions(CBLAS::CBLAS INTERFACE BLASW_CBLAS_MKL)
    endif()
endif()

mark_as_advanced(CBLAS_INCLUDE_DIR CBLAS_LIBRARY)
