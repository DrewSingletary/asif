set(LIB_NAME "libaffa")

file(GLOB LIBAFFA_SRC "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
file(GLOB LIBAFFA_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h*")

add_library(${LIB_NAME} ${LIBAFFA_SRC})

target_include_directories(${LIB_NAME}
                           PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>"
                                  "$<INSTALL_INTERFACE:include/${LIB_NAME}>")

set_target_properties(${LIB_NAME} PROPERTIES PUBLIC_HEADER "${LIBAFFA_HEADERS}")

# Specify installation options
install(TARGETS ${LIB_NAME}
        EXPORT  ${LIB_NAME}Targets
        ARCHIVE DESTINATION "lib/${LIB_NAME}"
        LIBRARY DESTINATION "lib/${LIB_NAME}"
        RUNTIME DESTINATION "bin/${LIB_NAME}"
        PUBLIC_HEADER DESTINATION "include/${LIB_NAME}")

install(EXPORT "${LIB_NAME}Targets"
	FILE "${LIB_NAME}Targets.cmake"
	NAMESPACE "${PROJECT_NAME}::"
	DESTINATION "share/${LIB_NAME}")

write_basic_package_version_file(
	"${CMAKE_CURRENT_BINARY_DIR}/${LIB_NAME}ConfigVersion.cmake"
	VERSION ${PROJECT_VERSION}
	COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
	"${PROJECT_SOURCE_DIR}/cmake/config${LIB_NAME}.cmake.in"
	"${CMAKE_CURRENT_BINARY_DIR}/${LIB_NAME}Config.cmake"
	INSTALL_DESTINATION share/${LIB_NAME}/cmake
)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${LIB_NAME}Config.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/${LIB_NAME}ConfigVersion.cmake"
	DESTINATION "share/${LIB_NAME}")
