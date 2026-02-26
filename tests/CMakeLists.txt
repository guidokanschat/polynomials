find_package(Catch2 3 REQUIRED)

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/testfiles)

add_executable(test_legendre test_legendre.cpp)
target_link_libraries(test_legendre PRIVATE Catch2::Catch2WithMain)

add_executable(test_fraction test_fraction.cpp)
target_link_libraries(test_fraction PRIVATE Catch2::Catch2WithMain)

add_executable(test_helpfunctions test_helpfunctions.cc)
target_link_libraries(test_helpfunctions PRIVATE Catch2::Catch2WithMain)

#find_package(Boost REQUIRED)
#target_link_libraries(test_legendre PRIVATE Boost)

