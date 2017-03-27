cmake_minimum_required (VERSION 2.8.8)

# If on Ubuntu, get the Ubuntu release codename
find_program   ( LSB_RELEASE lsb_release )
execute_process( COMMAND ${LSB_RELEASE} -cs OUTPUT_VARIABLE LSB_RELEASE_CODE OUTPUT_STRIP_TRAILING_WHITESPACE)



execute_process(
	COMMAND git describe --tags --long
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	OUTPUT_VARIABLE CATH_TOOLS_GIT_VERSION
	OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
	COMMAND git log -1 --date=short --pretty=format:%cd
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	OUTPUT_VARIABLE CATH_TOOLS_GIT_DATE
	OUTPUT_STRIP_TRAILING_WHITESPACE
)

configure_file( "source/cath_tools_git_version.hpp.in"
                "source/cath_tools_git_version.hpp"
                )

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/" )

find_package( GSL REQUIRED )

list                  ( GET GSL_LIBRARIES 0 FIRST_GSL_LIBRARY                 )
get_filename_component( FIRST_GSL_LIBDIR ${FIRST_GSL_LIBRARY} PATH            )
set                   ( GSL_STATIC_LIB      ${FIRST_GSL_LIBDIR}/libgsl.a      )
set                   ( GSLCBLAS_STATIC_LIB ${FIRST_GSL_LIBDIR}/libgslcblas.a )
foreach (loop_var RANGE ${GSL_LIBRARIES})
  list(APPEND GSL_DYN_LINK_FLAGS "-l${loop_var}")
endforeach(loop_var)

# For useful information on getting CMake to work with Boost, type the following command: cmake --help-module FindBoost

if( BUILD_SHARED_LIBS )
	add_definitions( -DBOOST_ALL_DYN_LINK )
	add_definitions( -DBOOST_LOG_DYN_LINK )
	set( GSL_LIB_SUFFIX ${GSL_LIBRARIES} )
else()
	set( Boost_USE_STATIC_LIBS ON )
	SET( GSL_LIB_SUFFIX "${GSL_STATIC_LIB}" "${GSLCBLAS_STATIC_LIB}" )
endif()

if ( ${LSB_RELEASE_CODE} STREQUAL "yakkety" )
	set( GSL_LIB_SUFFIX ${GSL_LIBRARIES} )
endif()

find_library( RT_LIBRARY rt )

find_package(Boost 1.55 REQUIRED filesystem iostreams log_setup log program_options serialization system thread unit_test_framework)

# Default flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y -W -Wall -Werror -Wextra -Wsign-compare -Wsign-conversion -Wcast-qual -Wconversion -Wnon-virtual-dtor -pedantic -ftabstop=2")

# Append debug compiler flags with setting of _GLIBCXX_DEBUG
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_GLIBCXX_DEBUG")

# *** TEMPORARY *** Prevent clang v3.6.0 from failing on Boost Concept checks due to unused local typedefs 
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedef")
endif()

# Note: to use libc++, pass -DCMAKE_CXX_FLAGS="-stdlib=libc++" to cmake

if ( CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_GENERATOR STREQUAL "Ninja" )
	set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -fcolor-diagnostics" )
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fcolor-diagnostics" )
endif()

MESSAGE( STATUS "Using CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS} )

# Add handling of includes and library links for Boost
# 'SYSTEM' tells the compiler it's a system include directory so warnings should be suppressed
include_directories( SYSTEM ${Boost_INCLUDE_DIR} )
include_directories( SYSTEM ${GSL_INCLUDE_DIR}   )
include_directories( SYSTEM rapidjson/include    )
link_directories   (        ${Boost_LIBRARYDIR}  )
link_directories   (        ${Boost_LIBRARYDIR}  )

add_definitions(-DRAPIDJSON_HAS_STDSTRING=1)

# Add handling of includes and library links for bioplib
add_custom_target(
	built_bioplib_libraries
	COMMAND make
	# In CMake v3.2, the BYPRODUCTS field allows specifiecation of the expected resulting files - these can be used by Ninja to avoid calling this step every time
	#BYPRODUCTS ${CMAKE_CURRENT_LIST_DIR}/bioplib/src/libgen.a ${CMAKE_CURRENT_LIST_DIR}/bioplib/src/libbiop.a
	COMMENT Attempting to build bioplib
	WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/bioplib/src
)

link_directories( ${CMAKE_CURRENT_LIST_DIR}/bioplib/src )

set( CMAKE_INCLUDE_CURRENT_DIR      ON )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY .. )

if ( DEFINED ENV{BUILD_REVISION_NUMBER} )
  add_definitions(-DBUILD_REVISION_NUMBER=$ENV{BUILD_REVISION_NUMBER})
endif ( DEFINED ENV{BUILD_REVISION_NUMBER} )
if ( DEFINED ENV{BUILD_BRANCH_NAME} )
  add_definitions(-DBUILD_BRANCH_NAME=$ENV{BUILD_BRANCH_NAME})
endif ( DEFINED ENV{BUILD_BRANCH_NAME} )

add_subdirectory   ( source )
include_directories( source )
