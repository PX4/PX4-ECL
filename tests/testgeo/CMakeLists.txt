# Add testgeo-exectuable to the project.
add_executable(testgeo main.cpp testgeo.cpp)

# Specify libraries to be used for testgeo.
target_link_libraries(testgeo
    libgeo
    libgtest
    libgmock
)

# Add test such that it can be run with ctest
add_test(Name testgeo
		 Command testgeo)