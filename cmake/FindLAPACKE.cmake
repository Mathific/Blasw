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

set(LAPACKE_MKL OFF CACHE INTERNAL "")
find_package(LAPACK)

if(LAPACK_FOUND)
    include(CheckFunctionExists)
    set(CMAKE_REQUIRED_LIBRARIES ${LAPACK_LIBRARIES})
    set(CMAKE_REQUIRED_FLAGS ${LAPACK_LINKER_FLAGS})
    check_function_exists(LAPACKE_dgels TEMP_FOUND)

    if(TEMP_FOUND)
        set(LAPACKE_LIBRARY ${LAPACK_LIBRARIES})
        set(LAPACKE_LINKER_FLAGS ${LAPACK_LINKER_FLAGS})

        foreach(TEMP_NAME ${LAPACK_LIBRARIES})
            get_filename_component(TEMP_NAME "${TEMP_NAME}" NAME)
            if(TEMP_NAME MATCHES "libmkl.*.so|libmkl.*.a|mkl.*.lib|mkl.*.dll")
                set(LAPACKE_MKL ON CACHE INTERNAL "")
            endif()
        endforeach()
    endif()

    unset(TEMP_FOUND CACHE)
endif()

if(NOT LAPACKE_LIBRARY)
    find_library(LAPACKE_LIBRARY
        NAMES lapacke
        PATHS ${BLASW_PATH} $ENV{BLASW_PATH}
        PATH_SUFFIXES lib lib64)
endif()

if(LAPACKE_LIBRARY)
    foreach(TEMP_DIR ${LAPACKE_LIBRARY})
        get_filename_component(TEMP_DIR "${TEMP_DIR}" DIRECTORY)
        list(APPEND HINT_PATH "${TEMP_DIR}/../")
    endforeach()
endif()

if(LAPACKE_MKL)
    find_path(LAPACKE_INCLUDE_DIR
        NAMES mkl_lapacke.h
        PATHS ${HINT_PATH} "$ENV{MKLROOT}"
        PATH_SUFFIXES include)
endif()

if((NOT LAPACKE_INCLUDE_DIR) OR
        (LAPACKE_INCLUDE_DIR STREQUAL "LAPACKE_INCLUDE_DIR-NOTFOUND"))
    find_path(LAPACKE_INCLUDE_DIR
        NAMES lapacke.h
        PATHS ${HINT_PATH}
        PATH_SUFFIXES include)
endif()

unset(HINT_PATH)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
    LAPACKE
    FOUND_VAR LAPACKE_FOUND
    REQUIRED_VARS
    LAPACKE_LIBRARY
    LAPACKE_INCLUDE_DIR)

if(LAPACKE_FOUND)
    set(LAPACKE_LIBRARIES ${LAPACKE_LIBRARY})
    set(LAPACKE_INCLUDE_DIRS ${LAPACKE_INCLUDE_DIR})
endif()

if(LAPACKE_FOUND AND NOT TARGET LAPACKE::LAPACKE)
    add_library(LAPACKE::LAPACKE INTERFACE IMPORTED)
    set_target_properties(LAPACKE::LAPACKE PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${LAPACKE_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES "${LAPACKE_LIBRARIES}"
        INTERFACE_LINK_OPTIONS "${LAPACKE_LINKER_FLAGS}"
        INTERFACE_COMPILE_DEFINITIONS BLASW_LAPACKE_FOUND)

    if(LAPACKE_MKL)
        target_compile_definitions(LAPACKE::LAPACKE INTERFACE BLASW_LAPACKE_MKL)
    endif()
endif()

mark_as_advanced(LAPACKE_INCLUDE_DIR LAPACKE_LIBRARY)
