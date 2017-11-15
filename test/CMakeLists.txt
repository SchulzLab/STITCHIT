project(LEGTest CXX)

include(LEGTest.cmake)

create_test_config_file()

add_test_case(runGenomeSizeReaderTests genomeSizeReaderTest.cpp)

add_test_case(runGTFReaderTests gtfReaderTest.cpp)

add_test_case(runExpressionReaderTests expressionReaderTest.cpp)

add_test_case(runSPANInputGeneratorTests spanInputGeneratorTest.cpp)
