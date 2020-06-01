cmake_minimum_required(VERSION 2.8.8)
# CaGe-A sources

file(GLOB SRC_FILES *.cc )
# We don't make a library at the moment
include_directories(${SAMTOOLS_DIR} ${HTSLIB_DIR})
message("Generating installdir.h...")
execute_process(COMMAND ${REFERENCE_DATA_MAKER} WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} RESULT_VARIABLE INSTALLDIR_RES)
if (INSTALLDIR_RES)
    message("Error: ${INSTALLDIR_RES}")
endif ()
message("Done generating installdir.h...")

add_executable(${PROJECT_NAME_STR} ${SRC_FILES}  )
# Link with required pthreads and zlib
target_link_libraries(${PROJECT_NAME_STR} ${ZLIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${SAMTOOLS_BAM_LIB} ${HTS_LIB})
# Static linking to samtools and htslib
#target_link_libraries(${PROJECT_NAME_STR} ${SAMTOOLS_BAM_LIB} ${HTS_LIB})
# Link with gcov for coverage evaluation
if (LINUX)
    target_link_libraries(${PROJECT_NAME_STR} gcov)
    #    target_link_libraries(${PROJECT_NAME_STR} "/usr/local/lib64/libstdc++.so.6.0.22")
endif ()


# We want to install in the project directory for compatibility
set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/" CACHE PATH "default install path" FORCE)
install(TARGETS ${PROJECT_NAME_STR} DESTINATION bin)
