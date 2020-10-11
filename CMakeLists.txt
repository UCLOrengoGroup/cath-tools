cmake_minimum_required (VERSION 3.12)

MESSAGE( STATUS "Running under CMake v" ${CMAKE_VERSION} )

project( cath-tools CXX )

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

configure_file( "source/src_external_info/cath/external_info/cath_tools_git_version_impl.hpp.in"
                "source/src_external_info/cath/external_info/cath_tools_git_version_impl.hpp" )

configure_file( "source/src_external_info/cath/external_info/cath_tools_cmake_dirs_impl.hpp.in"
                "source/src_external_info/cath/external_info/cath_tools_cmake_dirs_impl.hpp" )

# Set the default GSL library to be static, If on Linux, get the release codename
SET(USE_STATIC_GSL TRUE CACHE BOOL "Build with the static Gnu Scientific Library")

# It'd be good to be able to tell FindGSL.cmake that we want static libraries (on all except Ubuntu 16.10 (yakkety))
# and then use the target GSL::gsl but FindGSL.cmake doesn't appear to support that at present.
find_package( GSL REQUIRED )
list                  ( GET GSL_LIBRARIES 0 FIRST_GSL_LIBRARY                 )
get_filename_component( FIRST_GSL_LIBDIR ${FIRST_GSL_LIBRARY} PATH            )
set                   ( GSL_STATIC_LIB      ${FIRST_GSL_LIBDIR}/libgsl.a      )
set                   ( GSLCBLAS_STATIC_LIB ${FIRST_GSL_LIBDIR}/libgslcblas.a )
foreach( loop_var ${GSL_LIBRARIES} )
  list(APPEND GSL_DYN_LINK_FLAGS "-l${loop_var}")
endforeach()

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

# For useful information on getting CMake to work with Boost, type the following command: `cmake --help-module FindBoost`
# or add `SET( Boost_DEBUG "ON" )` to this file to get helpful debug information when running CMake
find_package( Boost 1.60 REQUIRED filesystem iostreams log program_options serialization timer unit_test_framework )

# Compiler options
SET( CMAKE_CXX_STANDARD   17    )
SET( CMAKE_CXX_EXTENSIONS "OFF" )

enable_testing()

# Add INTERFACE target for gsl
add_library                ( cath_tools_gsl        INTERFACE                                                 )
target_include_directories ( cath_tools_gsl SYSTEM INTERFACE ${GSL_INCLUDE_DIR}                              )
target_link_libraries      ( cath_tools_gsl        INTERFACE ${GSL_LIB_SUFFIX}                               )
target_compile_options     ( cath_tools_gsl        INTERFACE $<$<NOT:$<CONFIG:DEBUG>>:-DGSL_RANGE_CHECK_OFF> )

# Disallow cycles
set_property( GLOBAL PROPERTY GLOBAL_DEPENDS_NO_CYCLES ON )

add_subdirectory( source )
