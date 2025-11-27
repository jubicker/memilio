add_library(moments
    sir_metapop/library/moments/model.h
    sir_metapop/library/moments/parameters.h
)
target_link_libraries(moments PUBLIC memilio)
target_include_directories(moments PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(moments PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(smm_config sir_metapop/simulation_metapop.cpp)
target_link_libraries(smm_config PRIVATE memilio ode_sir smm moments)
target_compile_options(smm_config PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(moment_config sir_metapop/simulation_moments.cpp)
target_link_libraries(moment_config PRIVATE memilio ode_sir smm moments)
target_compile_options(moment_config PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(hybrid1_config sir_metapop/simulation_hybrid_smm_ode.cpp)
target_link_libraries(hybrid1_config PRIVATE memilio ode_sir smm moments hybrid)
target_compile_options(hybrid1_config PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(hybrid1_config_2r sir_metapop/simulation_hybrid_smm_ode_two_regions.cpp)
target_link_libraries(hybrid1_config_2r PRIVATE memilio ode_sir smm moments hybrid)
target_compile_options(hybrid1_config_2r PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(test_moments sir_metapop/tests/test_moment_equations.cpp)
target_link_libraries(test_moments PRIVATE memilio ode_sir moments)
target_compile_options(test_moments PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
