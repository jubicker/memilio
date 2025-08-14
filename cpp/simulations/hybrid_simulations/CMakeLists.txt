add_library(moments
    moment_equation_array.h
)
target_link_libraries(moments PUBLIC memilio)
target_include_directories(moments PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(moments PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(sir ode_sir.cpp)
target_link_libraries(sir PRIVATE memilio ode_sir)
target_compile_options(sir PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(smm_sir smm_sir.cpp)
target_link_libraries(smm_sir PRIVATE memilio ode_sir smm)
target_compile_options(smm_sir PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(abm_sir abm_sir.cpp)
target_link_libraries(abm_sir PRIVATE memilio abm)
target_compile_options(abm_sir PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(moments_sir moment_equations_sir.cpp)
target_link_libraries(moments_sir PRIVATE memilio ode_sir moments)
target_compile_options(moments_sir PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
