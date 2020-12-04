## for bash
include(FindUnixCommands)

## test data package info
set(TEST_PACKAGE_NAME "zz_source.tar.xz")
set(TEST_PACKAGE_TARGET "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_PACKAGE_NAME}")
set(TEST_SHA1_HASH 103b3c2ebe12fed6a4e47fac99c593439e1b48d5)

## list of test bundles
set(TESTDIRS 
  001_g16_wfn         ## test Gaussian16 wfn files
  002_g16_wfx         ## test Gaussian16 wfx files
  003_g16_fchk        ## test Gaussian16 fchk files
  004_orca4.0_wfn     ## test orca 4.0 wfn files
  005_orca4.0_molden  ## test orca 4.0 molden files
  006_orca4.2_wfn     ## test orca 4.2 wfn files
  007_orca4.2_wfx     ## test orca 4.2 wfx files
  008_orca4.2_molden  ## test orca 4.2 molden files
  009_psi4_molden     ## test psi4 1.4 molden files
  )

## define the macro that will run the tests
macro (runtests)
  get_filename_component(curdir "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

  foreach(test ${ARGN})
    ## gather info about the test
    include("${CMAKE_CURRENT_SOURCE_DIR}/${test}.cmake")

    ## gen test
    file(TO_NATIVE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/ref/${test}.pgout" _outgen)
    add_test(
      NAME "GEN:${curdir}/${test}" 
      COMMAND ${BASH} -c "$<TARGET_FILE:postg> ${COMMAND_A1} ${COMMAND_A2} ${COMMAND_FILE} ${COMMAND_AHF} > ${_outgen}"
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      CONFIGURATIONS gen
      )
    set_property(TEST "GEN:${curdir}/${test}" PROPERTY LABELS "gen ${label_list}")

    ## run test
    file(TO_NATIVE_PATH "${CMAKE_CURRENT_BINARY_DIR}/${test}.pgout" _outtest)
    add_test(
      NAME "-RUN-:${curdir}/${test}" 
      COMMAND ${BASH} -c "$<TARGET_FILE:postg> ${COMMAND_A1} ${COMMAND_A2} ${COMMAND_FILE} ${COMMAND_AHF} > ${_outtest}"
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
      )
    set_property(TEST "-RUN-:${curdir}/${test}" PROPERTY LABELS "${label_list}")

    ## check tests
    if (NUMDIFF_FOUND)
      add_test(
        NAME "CHECK:${curdir}/${test}" 
        COMMAND ${BASH} -c "${NUMDIFF_EXE} -s ' \t\n=,:;<>[](){}^' ${numdiff_opts} ${_outgen} ${_outtest}"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        )
      set_property(TEST "CHECK:${curdir}/${test}" PROPERTY LABELS "${label_list}")
      set_property(TEST "CHECK:${curdir}/${test}" PROPERTY DEPENDS "-RUN-:${curdir}/${test}")
    endif()
  endforeach()
endmacro()

## if the file exists, calculate the hash and verify
if (EXISTS "${TEST_PACKAGE_TARGET}")
  file(SHA1 "${TEST_PACKAGE_TARGET}" _hash)
  if ("${_hash}" STREQUAL "${TEST_SHA1_HASH}")
    set(TEST_PACKAGE_FOUND "YES")
  endif()
endif()
if(NOT "${TEST_PACKAGE_FOUND}")
  message(FATAL_ERROR "A package file with the correct sha1sum could not be found.")
endif()

## unpack data for the tests (source)
if ("${TEST_PACKAGE_FOUND}")
  if (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/zz_source/")
    message(STATUS "Unpacking the data for tests in ${CMAKE_CURRENT_SOURCE_DIR}/")
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar "xJf" "${TEST_PACKAGE_TARGET}" WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
  endif()

  ## create symlink for the test data (build)
  if (NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/zz_source/")
    message(STATUS "Creating symlink for tests in ${CMAKE_CURRENT_BINARY_DIR}/")
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink "${CMAKE_CURRENT_SOURCE_DIR}/zz_source" "${CMAKE_CURRENT_BINARY_DIR}/zz_source")
  endif()

  ## process all subdirectories
  foreach(dir ${TESTDIRS})
    add_subdirectory(${dir})
  endforeach()
else()
  message(STATUS "Skipping tests because test data package could not be found or downloaded.")
endif()
