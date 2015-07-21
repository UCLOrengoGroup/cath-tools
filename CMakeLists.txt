cmake_minimum_required (VERSION 2.8.8)

# For useful information on getting CMake to work with Boost, type the following command: cmake --help-module FindBoost

IF(${BUILD_SHARED_LIBS} MATCHES "ON")
	set( Boost_USE_STATIC_LIBS OFF )
	ADD_DEFINITIONS( -DBOOST_ALL_DYN_LINK )
	ADD_DEFINITIONS( -DBOOST_LOG_DYN_LINK )
ELSE()
	set( Boost_USE_STATIC_LIBS ON )
ENDIF()

find_package(Boost 1.55 REQUIRED filesystem iostreams log_setup log program_options serialization system thread unit_test_framework)

# Default flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y -W -Wall -Werror -Wextra -Wsign-compare -Wcast-qual -Wconversion -Wnon-virtual-dtor -pedantic -ftabstop=2")

# Append debug compiler flags with setting of _GLIBCXX_DEBUG
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_GLIBCXX_DEBUG")

# *** TEMPORARY *** Prevent clang v3.6.0 from failing on Boost Concept checks due to unused local typedefs 
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedef")
endif()

if ( CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_GENERATOR STREQUAL "Ninja" )
	set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -fcolor-diagnostics" )
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fcolor-diagnostics" )
endif()

# Add handling of includes and library links for Boost
# 'SYSTEM' tells the compiler it's a system include directory so warnings should be suppressed
include_directories( SYSTEM ${Boost_INCLUDE_DIR} )
link_directories   (        ${Boost_LIBRARYDIR}  )

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
