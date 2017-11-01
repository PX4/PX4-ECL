# Declare library
add_library(libgeo STATIC 
	geo.h geo.cpp
)

# Add include directory to the target.
target_include_directories(libgeo PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

