find_package(MPI REQUIRED)

function (add_clams_cc_mpi_test test_name)
    set(test_source "${test_name}.cpp")
    set(test_exe "${test_name}")
    add_executable(${test_exe} ${test_source})
    target_include_directories(${test_exe}
                             	PRIVATE ${PROJECT_SOURCE_DIR}/include)
  	target_include_directories(${test_exe} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_link_libraries(${test_exe} PUBLIC ygm::ygm)
    if (${TEST_WITH_SLURM})
        add_test(${test_exe} "srun" "${CMAKE_CURRENT_BINARY_DIR}/${test_exe}")
    else ()
        add_test(${test_exe} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} "4"
                 ${MPIEXEC_PREFLAGS} "${CMAKE_CURRENT_BINARY_DIR}/${test_exe}"
        )
    endif ()
endfunction ()

function (add_clams_cc_test test_name)
    set(test_source "${test_name}.cpp")
    set(test_exe "${test_name}")
    add_executable(${test_exe} ${test_source})
    target_include_directories(${test_exe}
                             	PRIVATE ${PROJECT_SOURCE_DIR}/include)
  	target_include_directories(${test_exe} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_link_libraries(${test_exe} PUBLIC ygm::ygm)
    add_test(${test_exe} "${CMAKE_CURRENT_BINARY_DIR}/${test_exe}")
endfunction ()


add_clams_cc_mpi_test(test_clustering_metrics)