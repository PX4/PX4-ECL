# geolib library
project("EKF-Lib" VERSION 0.1 LANGUAGES CXX)

# create library target
add_library(ekf
	src/airspeed_fusion.cpp 
	src/control.cpp
	src/covariance.cpp
	src/drag_fusion.cpp
	src/ekf_helper.cpp
	src/ekf.cpp
	src/estimator_interface.cpp
	src/gps_checks.cpp
	src/mag_fusion.cpp
	src/optflow_fusion.cpp
	src/sideslip_fusion.cpp
	src/terrain_estimator.cpp
	src/vel_pos_fusion.cpp
	)
add_library(ecl::ekflib ALIAS ekf)

# Specify libraries to be used for geolib.
target_link_libraries(ekf
    ecl::mathlib
    ecl::geolib
)
	
# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build
# paths.
target_include_directories(ekf PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE src ../matrixlib/matrix)
    
# If we have compiler requirements for this library, list them
# here
target_compile_features(ekf
    PUBLIC 
    	cxx_auto_type
    PRIVATE 
    	cxx_variadic_templates)
    
# 'make install' to the correct locations (provided by GNUInstallDirs).
install(TARGETS ekf EXPORT  
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
