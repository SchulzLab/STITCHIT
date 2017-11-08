project(LEGTest CXX)

include(LEGTest.cmake)

create_test_config_file()

add_test_case(runGenerateSPANInputTests generateSPANInputTest.cpp)
