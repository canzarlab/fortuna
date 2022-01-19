cmake_minimum_required( VERSION 2.8.12 )

file(GLOB sources *.cpp)
file(GLOB headers *.h *.hpp)

list(REMOVE_ITEM sources main.cpp)

include_directories(../ext/htslib)
include_directories(../ext/kallisto/src)

add_library(mt_core ${sources} ${headers})
target_include_directories(mt_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(mt main.cpp)

find_package( Threads REQUIRED )
target_link_libraries(mt mt_core pthread ${CMAKE_CURRENT_SOURCE_DIR}/../ext/htslib/libhts.a)

find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries(mt mt_core ${ZLIB_LIBRARIES})
else()
    message(FATAL_ERROR "zlib not found. Required for to output files" )
endif( ZLIB_FOUND )

install(TARGETS mt DESTINATION "${CMAKE_INSTALL_BINDIR}")
