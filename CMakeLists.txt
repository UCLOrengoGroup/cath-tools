cmake_minimum_required (VERSION 3.2)

# Note: to use libc++, pass -DCMAKE_CXX_FLAGS="-stdlib=libc++" to cmake

MESSAGE( STATUS "Running under CMake v" ${CMAKE_VERSION} )

# Set the default GSL library to be static, If on Linux, get the release codename
SET(USE_STATIC_GSL TRUE CACHE BOOL "Build with the static Gnu Scientific Library")

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
                "source/cath_tools_git_version/cath_tools_git_version.hpp"
                )

include_directories( $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/source/cath_tools_git_version> )


SET( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/" )
IF ( ( ${CMAKE_VERSION} VERSION_GREATER "3.5.0" ) AND ( ${CMAKE_VERSION} VERSION_LESS "3.6.0" ) )
	SET( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/FindBoostModules_in_3.5/" )
ENDIF()
IF ( ${CMAKE_VERSION} VERSION_LESS    "3.3.0" )
	SET( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/FindBoostModules_lt_3.3.0/" )
ENDIF()

# It'd be good to be able to tell FindGSL.cmake that we want static libraries (on all except Ubuntu 16.10 (yakkety))
# and then use the target GSL::gsl but FindGSL.cmake doesn't appear to support that at present.
find_package( GSL REQUIRED )
list                  ( GET GSL_LIBRARIES 0 FIRST_GSL_LIBRARY                 )
get_filename_component( FIRST_GSL_LIBDIR ${FIRST_GSL_LIBRARY} PATH            )
set                   ( GSL_STATIC_LIB      ${FIRST_GSL_LIBDIR}/libgsl.a      )
set                   ( GSLCBLAS_STATIC_LIB ${FIRST_GSL_LIBDIR}/libgslcblas.a )
foreach (loop_var RANGE ${GSL_LIBRARIES})
  list(APPEND GSL_DYN_LINK_FLAGS "-l${loop_var}")
endforeach(loop_var)

if( BUILD_SHARED_LIBS )
	add_definitions( -DBOOST_ALL_DYN_LINK )
	add_definitions( -DBOOST_LOG_DYN_LINK )
	SET( GSL_LIB_SUFFIX ${GSL_LIBRARIES} )
else()
	SET( Boost_USE_STATIC_LIBS ON )
	SET( GSL_LIB_SUFFIX "${GSL_STATIC_LIB}" "${GSLCBLAS_STATIC_LIB}" )
endif()

if ( NOT USE_STATIC_GSL )
	SET( GSL_LIB_SUFFIX ${GSL_LIBRARIES} )
	MESSAGE( STATUS "Directed to use dynamic Gnu Scientific Library, overriding GSL libraries to be : " ${GSL_LIB_SUFFIX} )
endif()

IF( ${CMAKE_SYSTEM_NAME} MATCHES "Linux" )
	find_library( RT_LIBRARY rt )
ENDIF()

# For useful information on getting CMake to work with Boost, type the following command: `cmake --help-module FindBoost`
# or add `SET( Boost_DEBUG "ON" )` to this file to get helpful debug information when running CMake
find_package( Boost 1.60 REQUIRED filesystem iostreams log program_options serialization timer unit_test_framework )

# Compiler options
#
# Not using the following due to excessive false-positives
#  * -Wdouble-promotion
#  * -Wold-style-cast   (fires in common/third_party_code/gnuplot-iostream.h)
#  * -Wuseless-cast
SET( CMAKE_CXX_STANDARD   14    )
SET( CMAKE_CXX_EXTENSIONS "OFF" )
add_compile_options(
	-W
	-Wall
	-Werror
	-Wextra
	-Wcast-qual
	-Wconversion
	-Wformat=2
	-Wnon-virtual-dtor
	-Wshadow
	-Wsign-compare
	-Wsign-conversion
	-pedantic
	-ftabstop=2
	$<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_GREATER:$<CXX_COMPILER_VERSION>,6.0.0>>:-Wduplicated-cond>
	# $<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_GREATER:$<CXX_COMPILER_VERSION>,6.0.0>>:-Wnull-dereference> # Causes problems with superposition_context_from_ptree() - appears to be quite a hassle to resolve
	$<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_GREATER:$<CXX_COMPILER_VERSION>,7.0.0>>:-Wduplicated-branches>
	$<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_GREATER:$<CXX_COMPILER_VERSION>,7.0.0>>:-Wrestrict>
)

# Append _GLIBCXX_DEBUG to GCC debug build to enable the debug mode of libstdc++
#
# ("By default, error messages are formatted to fit on lines of about 78 characters.
#   The environment variable GLIBCXX_DEBUG_MESSAGE_LENGTH can be used to request a different length."
#  Could consider ~150)
add_compile_options( $<$<AND:$<CXX_COMPILER_ID:GNU>,$<CONFIG:DEBUG>>:-D_GLIBCXX_DEBUG> )

# *** TEMPORARY *** Prevent clang v3.6.0 from failing on Boost Concept checks due to unused local typedefs
add_compile_options($<$<AND:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_VERSION:3.6.0>>:-Wno-unused-local-typedef>)

# Add colour diagnostics when using Clang with Ninja
if ( CMAKE_GENERATOR STREQUAL "Ninja" )
	add_compile_options( $<$<CXX_COMPILER_ID:Clang>:-fcolor-diagnostics> )
endif()

# Add handling of includes and library links
# 'SYSTEM' tells the compiler it's a system include directory so warnings should be suppressed
include_directories( SYSTEM ${Boost_INCLUDE_DIR} )
include_directories( SYSTEM ${GSL_INCLUDE_DIR}   )
include_directories( SYSTEM rapidjson/include    )

# Define GSL_RANGE_CHECK_OFF in non-debug builds
add_compile_options( $<$<NOT:$<CONFIG:DEBUG>>:-DGSL_RANGE_CHECK_OFF> )

add_definitions( -DRAPIDJSON_HAS_STDSTRING=1 )

SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY .. )

if ( DEFINED ENV{BUILD_REVISION_NUMBER} )
  add_definitions( -DBUILD_REVISION_NUMBER=$ENV{BUILD_REVISION_NUMBER} )
endif ( DEFINED ENV{BUILD_REVISION_NUMBER} )
if ( DEFINED ENV{BUILD_BRANCH_NAME} )
  add_definitions( -DBUILD_BRANCH_NAME=$ENV{BUILD_BRANCH_NAME} )
endif ( DEFINED ENV{BUILD_BRANCH_NAME} )

add_subdirectory( source )
