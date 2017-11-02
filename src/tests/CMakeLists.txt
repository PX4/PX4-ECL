# Test project for ecl
project("ECL-Test" VERSION 0.1 LANGUAGES CXX)

# Googletest needs thread support.
find_package(Threads REQUIRED)

# Enable ExternalProject CMake module for GoogleTest
include(ExternalProject)
ExternalProject_Add(
    gtest
    URL https://github.com/google/googletest/archive/master.zip
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/gtest
    # Disable install step
    INSTALL_COMMAND ""
)

# Get GTest source and binary directories from CMake project
ExternalProject_Get_Property(gtest source_dir binary_dir)

# Create a libgtest target to be used as a dependency by test programs
add_library(gtestlib IMPORTED STATIC GLOBAL)
add_dependencies(gtestlib gtest)


# Set libgtest properties
set_target_properties(gtestlib PROPERTIES
    "IMPORTED_LOCATION" "${binary_dir}/googlemock/gtest/libgtest.a"
    "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
)

# Create a libgmock target to be used as a dependency by test programs
add_library(gmocklib IMPORTED STATIC GLOBAL)
add_dependencies(gmocklib gtest)


# Set libgmock properties
set_target_properties(gmocklib PROPERTIES
    "IMPORTED_LOCATION" "${binary_dir}/googlemock/libgmock.a"
    "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
)

include_directories("${source_dir}/googletest/include"
                    "${source_dir}/googlemock/include")
                    

# Add subdirectories                   
add_subdirectory(testmath)