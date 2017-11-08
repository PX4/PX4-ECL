# geolib library
project("FixedwingAttitude-Lib" VERSION 0.1 LANGUAGES CXX)

# create library target
add_library(attitude_fw
	src/controller.cpp
	src/pitch_controller.cpp
	src/roll_controller.cpp
	src/wheel_controller.cpp
	src/yaw_controller.cpp
)
add_library(ecl::attitude_fw ALIAS attitude_fw)

# Specify libraries to be used for geolib.
target_link_libraries(attitude_fw
    ecl::mathlib
    ecl::geolib
)
	
# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build
# paths.
target_include_directories(attitude_fw PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE src)
    
# If we have compiler requirements for this library, list them
# here
target_compile_features(attitude_fw
    PUBLIC 
    	cxx_auto_type
    PRIVATE 
    	cxx_variadic_templates)
    
# 'make install' to the correct locations (provided by GNUInstallDirs).
install(TARGETS attitude_fw EXPORT  
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR})  # This is for Windows
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# This makes the project importable from the install directory
# Put config file in per-project dir (name MUST match), can also
# just go into 'cmake'.
#install(EXPORT MathLibConfig DESTINATION share/MathLibrary/cmake)

# This makes the project importable from the build directory
#export(TARGETS mathlib FILE MathLibConfig.cmake)
