add_library(moments
    model.h
    model.cpp
    parameters.h
    simulation.h
    simulation.cpp
    closure_functions.h
    closure_functions.cpp
)
target_link_libraries(moments PUBLIC memilio)
target_include_directories(moments PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(moments PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
