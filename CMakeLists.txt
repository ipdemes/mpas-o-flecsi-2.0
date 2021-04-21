#------------------------------------------------------------------------------#
# Copyright (c) 2016, Triad National Security, LLC
# All rights reserved
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.12...3.19)

#------------------------------------------------------------------------------#
# Set the project name.
#------------------------------------------------------------------------------#

project(mpasoflecsi LANGUAGES CXX C)

#------------------------------------------------------------------------------#
# FleCSI currently depends on C++17.
#------------------------------------------------------------------------------#

set(CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)

#------------------------------------------------------------------------------#
# Compiler version requirements.
#------------------------------------------------------------------------------#

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    message(FATAL_ERROR "Version 9.0 of gnu compilers required!")
  endif()
endif()

#------------------------------------------------------------------------------#
# enable testing
#------------------------------------------------------------------------------#
enable_testing()

#------------------------------------------------------------------------------#
# Add local module path.
#------------------------------------------------------------------------------#

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

#------------------------------------------------------------------------------#
# Add flecsi module path.
#------------------------------------------------------------------------------#

find_package(FleCSICMake REQUIRED)
list(APPEND CMAKE_MODULE_PATH ${FleCSICMake_PREFIX})

#----------------------------------------------------------------------------#
# These variables are used to collect library and include dependencies
# for the MPASConfig file below.
#----------------------------------------------------------------------------#

set(MPAS_INCLUDE_DEPENDENCIES)
set(MPAS_LIBRARY_DEPENDENCIES)

#------------------------------------------------------------------------------#
# Find the FleCSI installation.
#
# The easiest way to help CMake find FleCSI is to install it,
# and add the location to the CMAKE_PREFIX_PATH environment variable.
#------------------------------------------------------------------------------#

find_package(FleCSI REQUIRED)

#----------------------------------------------------------------------------#
# Unit tests and logging.
#----------------------------------------------------------------------------#

set(ENABLE_FLOG ${FleCSI_ENABLE_FLOG})
set(FLECSI_RUNTIME_MODEL ${FleCSI_RUNTIME_MODEL})

include(flog)
set(FLECSI_UNIT_MAIN ${FleCSICMake_PREFIX}/unit-main.cc)
include(unit)


message (STATUS "FleCSI runtime model: " ${FleCSI_RUNTIME_MODEL})

list (APPEND MPAS_LIBRARY_DEPENDENCIES FleCSI::FleCSI)

#------------------------------------------------------------------------------#
# HDF5
#------------------------------------------------------------------------------#
#
find_package(HDF5 COMPONENTS C HL REQUIRED)

include_directories( ${HDF5_INCLUDE_DIRS} )
list(APPEND MPAS_LIBRARY_DEPENDENCIES ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
message( STATUS "IO with HDF5 is enabled" )

#------------------------------------------------------------------------------#
# Apps
#------------------------------------------------------------------------------#

add_subdirectory(apps)

#------------------------------------------------------------------------------#
# MPASOFLeCSI library
#------------------------------------------------------------------------------#

include(library)
add_library_target(mpasoflecsi mpasoflecsi EXPORT_TARGET MPASOFleCSITargets)
message(STATUS "MPAS Libs: " "${MPAS_LIBRARY_DEPENDENCIES}")
target_link_libraries(mpasoflecsi LINK_PUBLIC ${MPAS_LIBRARY_DEPENDENCIES})

#----------------------------------------------------------------------------#
# CMake config file: This should be the last thing to happen.
#----------------------------------------------------------------------------#
set(MPAS_CMAKE_DIR ${CMAKE_INSTALL_PREFIX}/${LIBDIR}/cmake/mapsoflecsi)

configure_file(${PROJECT_SOURCE_DIR}/config/mpasoflecsi-config.h.in
  ${CMAKE_BINARY_DIR}/mpasoflecsi-config.h @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/mpasoflecsi-config.h DESTINATION include)

#install(
#    EXPORT MPASOFleCSITargets
#    DESTINATION include
#  )

include_directories(${CMAKE_BINARY_DIR})

