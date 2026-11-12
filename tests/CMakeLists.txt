add_executable(CMFSNTests FVMTests.cpp)
target_link_libraries(CMFSNTests PRIVATE Solver IO)
if(Kokkos_ENABLE_CUDA)
    set_source_files_properties(FVMTests.cpp PROPERTIES LANGUAGE CUDA)
    if(MSVC)
        target_compile_options(CMFSNTests PRIVATE "-Xcompiler=/std:c++20" "-Xcompiler=/Zc:preprocessor")
    endif()
endif()

foreach(test_case IN ITEMS indexing thermodynamics geometry gradients constitutive periodic wall_adiabatic wall_isothermal jst farfield viscous_flux time_step constructors rk_jameson rk_ck export)
    foreach(precision IN ITEMS float double)
        set(test_name "fvm_${test_case}_${precision}")
        add_test(NAME ${test_name} COMMAND CMFSNTests ${test_case} ${precision})
        set_tests_properties(${test_name} PROPERTIES
            LABELS "fvm;${test_case};${precision}"
            TIMEOUT 60
            ENVIRONMENT "OMP_NUM_THREADS=2;OMP_PROC_BIND=false")
        if(Kokkos_ENABLE_CUDA)
            set_tests_properties(${test_name} PROPERTIES RESOURCE_LOCK cmfsn_gpu)
        endif()
    endforeach()
endforeach()
