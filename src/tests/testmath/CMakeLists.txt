# Add testgeo-exectuable to the project.
add_executable(testmath main.cpp testmath.cpp)

# Specify libraries to be used for testgeo.
target_link_libraries(testmath
    ecl::mathlib
    gtestlib
    gmocklib
)

# Add test such that it can be run with ctest
add_test(Name testmath
		 Command testmath)