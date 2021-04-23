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

configure_file( "source/ct_external_info/cath/external_info/cath_tools_git_version_impl.hpp.in"
                "source/ct_external_info/cath/external_info/cath_tools_git_version_impl.hpp" )

configure_file( "source/ct_external_info/cath/external_info/cath_tools_cmake_dirs_impl.hpp.in"
                "source/ct_external_info/cath/external_info/cath_tools_cmake_dirs_impl.hpp" )

set( CONAN_SYSTEM_INCLUDES ON )

include( ${CMAKE_BINARY_DIR}/conanbuildinfo.cmake )

conan_basic_setup( TARGETS KEEP_RPATHS )

# It'd be good to be able to tell FindGSL.cmake that we want static libraries
# (and then just use GSL::gsl target it gives us) but FindGSL.cmake doesn't
# appear to support that at present.
#
# So here, manaully construct where the .a files are probably located
find_package( GSL REQUIRED )
list                  ( GET GSL_LIBRARIES 0 FIRST_GSL_LIBRARY                 )
get_filename_component( FIRST_GSL_LIBDIR ${FIRST_GSL_LIBRARY} PATH            )
set                   ( GSL_STATIC_LIB      ${FIRST_GSL_LIBDIR}/libgsl.a      )
set                   ( GSLCBLAS_STATIC_LIB ${FIRST_GSL_LIBDIR}/libgslcblas.a )

if( BUILD_SHARED_LIBS AND NOT Boost_USE_STATIC_LIBS )
	# add_definitions( -DBOOST_ALL_DYN_LINK )
	# add_definitions( -DBOOST_LOG_DYN_LINK )
	SET( GSL_LIB_SUFFIX ${GSL_LIBRARIES} )
else()
	# SET( Boost_USE_STATIC_LIBS ON )
	SET( GSL_LIB_SUFFIX "${GSL_STATIC_LIB}" "${GSLCBLAS_STATIC_LIB}" )
endif()

# For useful information on getting CMake to work with Boost, type the following command: `cmake --help-module FindBoost`
# or add `SET( Boost_DEBUG "ON" )` to this file to get helpful debug information when running CMake
find_package( Boost 1.60 REQUIRED filesystem iostreams log program_options serialization timer unit_test_framework )

# Compiler options
SET( CMAKE_CXX_STANDARD   14    )
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
