############################################################################
#
#   Copyright (c) 2017 ECL Development Team. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in
#    the documentation and/or other materials provided with the
#    distribution.
# 3. Neither the name ECL nor the names of its contributors may be
#    used to endorse or promote products derived from this software
#    without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
# AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
############################################################################

# All cmake commands in this file require at least version 3.2
cmake_minimum_required(VERSION 3.2 FATAL_ERROR)

# The name of the entire project: ECL = Estimation and Control Library.
# CXX is the supported language for this project.
project("ECL" VERSION 0.1 LANGUAGES CXX)

# Must use GNUInstallDirs to install libraries into correct
# locations on all platforms.
include(GNUInstallDirs)

# Check if matrix submodule is available
if(NOT EXISTS ${CMAKE_SOURCE_DIR}/matrixlib/.git)
	message(SEND_ERROR "The git submodules are not available. Please run git submodule update --init --recursive")
endif()

# Initialize built-in CMake variable that differs from default.
# We need to do that only once, the first time cmake is run.
# if (NOT CONFIGURED_ONCE)

	# compiler flag 
	 #set(CMAKE_CXX_FLAGS
	 	#"-Wall -Wextra -Werror"
	  #  CACHE STRING "Flags used by the compiler during all build types." FORCE)
#endif()

# Enable testing for this directory and below using ctest.
enable_testing()

# Add include directory to the build.
#include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory(mathlib)
add_subdirectory(geolib)
#add_subdirectory(attitude_fwlib)
add_subdirectory(matrixlib)
add_subdirectory(tecslib)
#add_subdirectory(tests)

# The first time cmake is run we create a new variable that indicates 
# that there has already been a cmake configuration. 
#set(CONFIGURED_ONCE TRUE CACHE INTERNAL
 #   "A flag showing that CMake has configured at least once.")