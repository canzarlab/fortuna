cmake_minimum_required( VERSION 2.8.12 )

project( fortuna )

include(GNUInstallDirs)
include(ExternalProject)

set(CMAKE_CXX_FLAGS_PROFILE "-g")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

if(${CMAKE_VERSION} VERSION_LESS 3.1)
    add_compile_options(-std=c++11)
else()
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
endif()

# add_compile_options(-Wall -Wno-unused-function -Wpedantic)

find_package(ZLIB)
add_library(zlib STATIC IMPORTED)
set_target_properties(zlib PROPERTIES IMPORTED_LOCATION ${ZLIB_LIBRARY})

#set(EXT_PROJECTS_DIR ${PROJECT_SOURCE_DIR}/ext)

ExternalProject_Add(htslib
    PREFIX ${PROJECT_SOURCE_DIR}/ext/htslib
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/ext/htslib
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND autoheader && autoconf && ${PROJECT_SOURCE_DIR}/ext/htslib/configure 
        --prefix=${PREFIX} --disable-bz2 --disable-lzma --disable-libcurl
    BUILD_COMMAND make lib-static
    INSTALL_COMMAND ""
)

if (CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    set(MAKE_COMMAND "$(MAKE)")
else()
    find_program(MAKE_COMMAND NAMES make gmake)
endif()

file(GLOB fortuna_sources src/*.cpp)
file(GLOB fortuna_headers src/*.h src/*.hpp)

file(GLOB ka_sources ext/kallisto/src/*.cpp)
file(GLOB ka_headers ext/kallisto/src/*.h ext/kallisto/src/*.hpp)

list(REMOVE_ITEM fortuna_sources src/main.cpp)
list(REMOVE_ITEM ka_sources ext/kallisto/src/main.cpp)

include_directories(${PROJECT_SOURCE_DIR}/ext/htslib)
include_directories(${PROJECT_SOURCE_DIR}/ext)
include_directories(${PROJECT_SOURCE_DIR}/src)

add_library(fortuna_core ${fortuna_sources} ${fortuna_headers} ${ka_sources} ${ka_headers})
target_include_directories(fortuna_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(fortuna ${PROJECT_SOURCE_DIR}/src/main.cpp)

find_package( Threads REQUIRED )
target_link_libraries(fortuna fortuna_core pthread ${PROJECT_SOURCE_DIR}/ext/htslib/libhts.a)

find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries(fortuna fortuna_core ${ZLIB_LIBRARIES})
else()
    message(FATAL_ERROR "zlib not found. Required for to output files" )
endif( ZLIB_FOUND )

install(TARGETS fortuna DESTINATION "${CMAKE_INSTALL_BINDIR}")

add_dependencies(fortuna htslib)

