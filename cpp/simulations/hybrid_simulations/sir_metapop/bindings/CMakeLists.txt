cmake_minimum_required(VERSION 3.13)
project(memilio-python)

set(CMAKE_CXX_STANDARD "20")
set(CMAKE_CXX_STANDARD_REQUIRED "20")

option(MEMILIO_USE_BUNDLED_PYBIND11 "Use pybind11 bundled with this library." ON)
mark_as_advanced(MEMILIO_USE_BUNDLED_PYBIND11)

include(GNUInstallDirs) # set to gnu folders. No cache variable so this is not global

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_PDB_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Windows has no RPATH (CMAKE_SKIP_RPATH defaults to ON there), so set this only on UNIX and APPLE.
# Mach-O (APPLE) uses @loader_path, ELF (other UNIX) uses $ORIGIN as the loader-relative token.
if(APPLE)
    set(MEMILIO_LOADER_PATH_TOKEN "@loader_path")
elseif(UNIX)
    set(MEMILIO_LOADER_PATH_TOKEN "$ORIGIN")
endif()

if(APPLE OR UNIX)
    set(CMAKE_INSTALL_RPATH
        "${MEMILIO_LOADER_PATH_TOKEN}"
        "${MEMILIO_LOADER_PATH_TOKEN}/lib"
    )
endif()

# keep build RPATH separate and avoid putting build paths into installed wheels
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)

# Use a static spdlog so the wheel stays self-contained.
set(SPDLOG_BUILD_SHARED OFF CACHE BOOL "" FORCE)

if(MEMILIO_USE_BUNDLED_PYBIND11)
    # Fetch pybind11
    include(FetchContent)
    FetchContent_Declare(pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11
        GIT_TAG v3.0.0
    )
    FetchContent_GetProperties(pybind11)

    if(NOT pybind11_POPULATED)
        FetchContent_Populate(pybind11)
        add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
    endif()
else()
    find_package(pybind11 REQUIRED)
endif()

# add in C++ library
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../../../ ${CMAKE_CURRENT_BINARY_DIR}/cpp EXCLUDE_FROM_ALL)

# a list of all "LINKED_LIBRARIES" that are given to add_pymio_module. will contain duplicates
# used for wheel installation
set(PYMIO_MEMILIO_LIBS_LIST)

# Function to add a pybind11 module
function(add_pymio_module target_name)
    set(options)
    set(oneValueArgs)
    set(multiValueArgs LINKED_LIBRARIES SOURCES)
    cmake_parse_arguments(PYBIND11_MODULE "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    pybind11_add_module(${target_name} MODULE ${PYBIND11_MODULE_SOURCES})
    target_link_libraries(${target_name} PRIVATE ${PYBIND11_MODULE_LINKED_LIBRARIES})
    target_include_directories(${target_name} PRIVATE memilio/simulation/bindings)

    # Collect linked libraries for later installation
    list(APPEND PYMIO_MEMILIO_LIBS_LIST "${PYBIND11_MODULE_LINKED_LIBRARIES}")
    set(PYMIO_MEMILIO_LIBS_LIST ${PYMIO_MEMILIO_LIBS_LIST} PARENT_SCOPE)

    # Set RPATH so extension can find the C++ libraries in lib subdirectory
    set_target_properties(${target_name} PROPERTIES
        INSTALL_RPATH "${MEMILIO_LOADER_PATH_TOKEN}/lib"
    )

    # Make sure to install shared and static libraries (LIBRARY and ARCHIVE) and the executables (RUNTIME).
    # Destination "." refers to the value of wheel.install-dir in the pyproject.toml.
    install(TARGETS ${target_name}
        LIBRARY DESTINATION .
        ARCHIVE DESTINATION .
        RUNTIME DESTINATION .
    )
endfunction()

# build python extensions
add_pymio_module(_simulation_stochastic
    LINKED_LIBRARIES memilio ode_sir smm moments hybrid
    SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/stochastic_sim.cpp
)
target_include_directories(_simulation_stochastic PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../applications/Influenza_simulations)

# install all shared memilio libraries, which were given as "LINKED_LIBRARIES" to add_pymio_module
list(REMOVE_DUPLICATES PYMIO_MEMILIO_LIBS_LIST)

# Set RPATH for all C++ libraries so they can find each other
foreach(target_name IN LISTS PYMIO_MEMILIO_LIBS_LIST)
    if(TARGET ${target_name})
        set_target_properties(${target_name} PROPERTIES
            INSTALL_RPATH "${MEMILIO_LOADER_PATH_TOKEN};${MEMILIO_LOADER_PATH_TOKEN}/lib"
        )
    endif()
endforeach()

install(TARGETS ${PYMIO_MEMILIO_LIBS_LIST}
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION lib
)

# If jsoncpp comes from the bundled third-party dependencies, copy its shared
# library into the wheel under every file name importers might look for. Wheels
# drop symlinks, so place real copies for each alias.
if(TARGET jsoncpp_lib)
    set_target_properties(jsoncpp_lib PROPERTIES
        INSTALL_RPATH "${MEMILIO_LOADER_PATH_TOKEN}"
    )
    install(FILES
        $<TARGET_FILE:jsoncpp_lib>
        $<TARGET_SONAME_FILE:jsoncpp_lib>
        $<TARGET_LINKER_FILE:jsoncpp_lib>
        DESTINATION lib
    )
endif()

if(TARGET spdlog)
    get_target_property(_SPDLOG_LIB_TYPE spdlog TYPE)

    if(_SPDLOG_LIB_TYPE STREQUAL "SHARED_LIBRARY")
        set_target_properties(spdlog PROPERTIES
            INSTALL_RPATH "${MEMILIO_LOADER_PATH_TOKEN}"
        )
        install(FILES
            $<TARGET_FILE:spdlog>
            $<TARGET_SONAME_FILE:spdlog>
            $<TARGET_LINKER_FILE:spdlog>
            DESTINATION lib
        )
    endif()
endif()
