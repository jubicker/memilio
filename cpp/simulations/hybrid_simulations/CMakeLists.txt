if(MEMILIO_ENABLE_OPENMP)
    add_executable(smm_config sir_metapop/applications/simulation_metapop.cpp)
    target_link_libraries(smm_config PRIVATE memilio ode_sir smm moments)
    target_compile_options(smm_config PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

    add_executable(hybrid1_config sir_metapop/applications/simulation_hybrid_smm_ode.cpp)
    target_link_libraries(hybrid1_config PRIVATE memilio ode_sir smm moments hybrid)
    target_compile_options(hybrid1_config PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

    add_executable(hybrid2_config sir_metapop/applications/simulation_hybrid_smm_set_moments.cpp)
    target_link_libraries(hybrid2_config PRIVATE memilio ode_sir smm moments hybrid)
    target_compile_options(hybrid2_config PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
endif()

add_executable(moment_config sir_metapop/applications/simulation_moments.cpp)
target_link_libraries(moment_config PRIVATE memilio ode_sir smm moments)
target_compile_options(moment_config PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(test_moments sir_metapop/tests/test_moment_equations.cpp)
target_link_libraries(test_moments PRIVATE memilio ode_sir moments)
target_compile_options(test_moments PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(test_moment_closures sir_metapop/tests/test_moment_closures.cpp)
target_link_libraries(test_moment_closures PRIVATE memilio ode_sir moments)
target_compile_options(test_moment_closures PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
