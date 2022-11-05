include(GNUInstallDirs)

add_library(galley INTERFACE)
# ALIAS target to export a namespaced target even when building in-tree
add_library(Galley::galley ALIAS galley)
# library public header files location
target_include_directories(
    galley INTERFACE
    # different include directory path depending on if using the local or installed version of library
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
# set up version and soversion for the main library object
set_target_properties(
    galley PROPERTIES
    INTERFACE_Galley_MAJOR_VERSION ${PROJECT_VERSION_MAJOR}
    INTERFACE_Galley_MINOR_VERSION ${PROJECT_VERSION_MINOR}
)
set_property(
    TARGET galley
    APPEND PROPERTY COMPATIBLE_INTERFACE_STRING "${Galley_MAJOR_VERSION}.${Galley_MINOR_VERSION}"
)

# install if we're not being built as a sub-project
if (NOT GALLEY_SUBPROJECT)
    # library
    install(
        TARGETS galley
        EXPORT GalleyTargets
        # when a static library
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        # when a shared library on UNIX
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        # when a DLL on Windows
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
    # public headers
    install(
        DIRECTORY "include/" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )
    # CMake export
    install(
        EXPORT GalleyTargets
        FILE GalleyTargets.cmake
        NAMESPACE Galley::
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Galley"
    )
    # CMake package
    include(CMakePackageConfigHelpers)
    configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
        "${CMAKE_CURRENT_BINARY_DIR}/GalleyConfig.cmake"
        INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Galley"
    )
    # generate the version file for the config file
    write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/GalleyConfigVersion.cmake"
        VERSION "${version}"
        COMPATIBILITY SameMinorVersion
    )
    install(
        FILES
        "${CMAKE_CURRENT_BINARY_DIR}/GalleyConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/GalleyConfigVersion.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Galley"
    )
endif()
