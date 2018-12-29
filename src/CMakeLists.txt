cmake_minimum_required(VERSION 2.8.8)
file(GLOB SRC_FILES *.cc util/*.cc )
include_directories(${SAMTOOLS_DIR} ${HTSLIB_DIR})

add_executable(${PROJECT_NAME_STR} ${SRC_FILES}  )
# Link with required pthreads and zlib
target_link_libraries(${PROJECT_NAME_STR} ${ZLIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
# Static linking to samtools and htslib
target_link_libraries(${PROJECT_NAME_STR} ${SAMTOOLS_BAM_LIB} ${HTS_LIB})

set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/" CACHE PATH "default install path" FORCE)
install(TARGETS ${PROJECT_NAME_STR} DESTINATION bin)
