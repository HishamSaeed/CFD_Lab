cmake_minimum_required(VERSION 3.10)
include_directories (${CFDLAB_SOURCE_DIR}) 

set(OPENSSL_USE_STATIC_LIBS TRUE)
find_package(OpenSSL REQUIRED)

# Test Executables
# add_executable(test_factorial sample_test.cpp ) 
# # add_executable(unit_tests ws1_unit_testing.cpp ../init.cpp ../helper.cpp ../grid.cpp 
# # ../utilities.cpp ../cell.cpp ../boundary_val.cpp ../uvp.cpp ../sor.cpp ./test_utilities.cpp)


# add_executable(unit_tests2 ws2_unit_testing.cpp ../init.cpp ../helper.cpp ../grid.cpp 
# ../utilities.cpp ../cell.cpp ../boundary_val.cpp ../uvp.cpp ../sor.cpp ./test_utilities.cpp)

# # target_link_libraries(unit_tests PRIVATE OpenSSL::Crypto)
# target_link_libraries(unit_tests2 PRIVATE OpenSSL::Crypto)

# # if you want to copy config files  Config Files
# # configure_file(${CFDLAB_SOURCE_DIR}/source.dat ${CMAKE_BINARY_DIR}/tests/target_example.dat COPYONLY)

# # Testing
# add_test (NAME example_test COMMAND test_factorial)
# Example for MPI RUN
# add_test (NAME example_test COMMAND mpirun -np 4 test [2x2])
