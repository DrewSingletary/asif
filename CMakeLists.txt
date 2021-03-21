cmake_minimum_required (VERSION 3.10)

# Project name
project(asif++ VERSION 1.0 LANGUAGES CXX)

# Include utilities
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Set various flags
set(default_build_type "Release")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -O3 -DNDEBUG")

# Build examples option
option(BUILD_EXAMPLES "Build example files." ON)

# OSQP version option
option(OSQP_EMBEDDED "Use embedde version of osqp." OFF)
option(USE_ODEINT "Use odeint library" OFF)

# debug option
option(WITH_ASIF_DEBUG "Verbose a lot of stuff" OFF)
if(WITH_ASIF_DEBUG)
	add_definitions(-DASIF_DEBUG)
endif()

# Release type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()


# Chose osqp version
if(OSQP_EMBEDDED)
	message(STATUS "Using embedded version of osqp.")
	set(OSQP_MODE "osqp_embedded")
	set(OSQP_MODE_LIB "osqp_embedded")
else()
	set(OSQP_MODE "osqp")
	set(OSQP_MODE_LIB "osqpstatic")
endif()

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/lib/libaffa")

# Find OSQP library and headers
find_package(${OSQP_MODE} REQUIRED CONFIG)

if (USE_ODEINT)
	find_package(Boost)
	if (NOT DEFINED Boost_FOUND)
		message(STATUS "Boost not found, using default integrator")
	else()
		message(STATUS "Boost found, using odeint for integration")
		add_definitions(-DUSE_ODEINT)
	endif()
endif()

# Add sources
set(SRC "${CMAKE_CURRENT_SOURCE_DIR}/src/asif.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/asif_robust.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/asif_realizable.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/asif_implicit.cpp"
				"${CMAKE_CURRENT_SOURCE_DIR}/src/asif_implicit_robust.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/asif_implicit_tb.cpp"
         "${CMAKE_CURRENT_SOURCE_DIR}/src/qpwrapper_abstract.cpp"
         "${CMAKE_CURRENT_SOURCE_DIR}/src/qpwrapper_osqp.cpp")
file(GLOB "${PROJECT_NAME}_HEADERS" "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h*")

# Make library
add_library(${PROJECT_NAME} STATIC ${SRC} "${${PROJECT_NAME}_HEADERS}")
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

# Set include directory
target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
                                                  "$<INSTALL_INTERFACE:include/${PROJECT_NAME}>")

set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${${PROJECT_NAME}_HEADERS}")


# Link with other libraries
target_link_libraries(${PROJECT_NAME} PUBLIC libaffa
                                             ${OSQP_MODE}::${OSQP_MODE_LIB}
                                             -ldl)

# Activate c++11
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)

# Set installation options
# ----------------------------------------------

install(TARGETS ${PROJECT_NAME} 
        EXPORT  ${PROJECT_NAME}Targets
        ARCHIVE DESTINATION "lib/${PROJECT_NAME}"
        LIBRARY DESTINATION "lib/${PROJECT_NAME}"
        RUNTIME DESTINATION "bin/${PROJECT_NAME}"
        PUBLIC_HEADER DESTINATION "include/${PROJECT_NAME}")

install(EXPORT "${PROJECT_NAME}Targets"
        FILE "${PROJECT_NAME}Targets.cmake"
        NAMESPACE "${PROJECT_NAME}::"
        DESTINATION "share/${PROJECT_NAME}")

write_basic_package_version_file(
	"${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
	VERSION ${PROJECT_VERSION}
	COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
	"${CMAKE_CURRENT_SOURCE_DIR}/cmake/config.cmake.in"
	"${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
	INSTALL_DESTINATION share/${PROJECT_NAME}/cmake
)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
	DESTINATION "share/${PROJECT_NAME}")

# uninstall target
configure_file(
	"${CMAKE_CURRENT_SOURCE_DIR}/cmake/uninstall.cmake.in"
	"${CMAKE_CURRENT_BINARY_DIR}/uninstall.cmake"
	IMMEDIATE @ONLY)

add_custom_target(uninstall COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/uninstall.cmake)

if(BUILD_EXAMPLES)
	add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/examples")
endif()