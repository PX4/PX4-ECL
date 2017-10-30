# Declare library
add_library(libgeo STATIC geo.cpp)

# Specify the included directories exported by thei lib
target_include_directories(libgeo PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

