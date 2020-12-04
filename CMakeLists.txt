cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

project(postg
        VERSION 0.1
        LANGUAGES Fortran)
set(postg_AUTHORS      "A. Otero de la Roza and E. R. Johnson")
set(postg_DESCRIPTION  "exchange-hole dipole moment (XDM) dispersion correction for DFT")
set(postg_URL          "https://github.com/aoterodelaroza/postg")
set(postg_LICENSE      "GPLv3")

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

################  Options: Overview and Not Otherwise Mentioned  ###############

#  <<<  CMake build overview  >>>
#
#    >>> ls
#    postg.g90 LICENSE param.f90 ...
#    >>> cmake -H. -Bobjdir -DCMAKE_INSTALL_PREFIX=/path/to/install-postg ...
#    ...
#    -- Generating done
#    -- Build files have been written to: /current/dir/objdir
#    >>> cd objdir && make -j`getconf _NPROCESSORS_ONLN`
#    >>> make install

#  <<<  Required build dependencies that POSTG can't build itself  >>>
#
#    - CMake (e.g., `conda install cmake`)
#    - Fortran compiler

############################  Options: Build How?  #############################

set(CMAKE_BUILD_TYPE Release)

#################################  Main Project  ################################
include(GNUInstallDirs)
include(autocmake_safeguards)

# <<<  Build >>>

set(sources_list atomicdata.f90 meshmod.f90 param.f90 postg.f90 tools_math.f90 wfnmod.f90)

add_executable(postg ${sources_list})
find_package(OpenMP MODULE COMPONENTS Fortran)
target_link_libraries(postg OpenMP::OpenMP_Fortran)

# <<<  Install  >>>

install(TARGETS postg
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# <<< Tests >>>
find_package(NUMDIFF)
option(BUILD_TESTING "Enable the regression tests for the critic2 build." OFF)
if (BUILD_TESTING AND NUMDIFF_FOUND)
  include(CTest)
  add_subdirectory(tests)
endif()
