add_executable(sir ode_sir.cpp)
target_link_libraries(sir PRIVATE memilio ode_sir)
target_compile_options(sir PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(smm_sir smm_sir.cpp)
target_link_libraries(smm_sir PRIVATE memilio ode_sir smm)
target_compile_options(smm_sir PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
