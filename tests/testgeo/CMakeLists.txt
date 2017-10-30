file(GLOB SRCS *.cpp)

add_executable(testgeo ${SRCS})

target_link_libraries(testgeo
    libgeo
    libgtest
    libgmock
)

install(TARGETS testgeo DESTINATION bin)