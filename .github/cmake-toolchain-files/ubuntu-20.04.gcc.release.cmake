set( CMAKE_BUILD_TYPE Release )

set( Boost_USE_STATIC_LIBS ON                                                   )
set( CMAKE_C_COMPILER      "/usr/bin/gcc-10"  CACHE FILEPATH "The C compiler"   )
set( CMAKE_CXX_COMPILER    "/usr/bin/g++-10"  CACHE FILEPATH "The C++ compiler" )
set( CMAKE_CXX_FLAGS_INIT  " ${CMAKE_CXX_FLAGS} -static "                       )

add_compile_options(
	-W
	-Wall
	-Werror
	-Wextra
	-Wcast-qual
	-Wconversion
	-Wformat=2
	-Wno-range-loop-analysis
	-Wnon-virtual-dtor
	-Wold-style-cast
	-Wshadow
	-Wsign-compare
	-Wsign-conversion
	-pedantic
	-ftabstop=2
	$<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_GREATER:$<CXX_COMPILER_VERSION>,4.8.0>>:-Wuseless-cast>
	$<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_GREATER:$<CXX_COMPILER_VERSION>,6.0.0>>:-Wduplicated-cond>
	$<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_GREATER:$<CXX_COMPILER_VERSION>,7.0.0>>:-Wduplicated-branches>
	$<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_GREATER:$<CXX_COMPILER_VERSION>,7.0.0>>:-Wrestrict>
	$<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_GREATER:$<CXX_COMPILER_VERSION>,8.0.0>>:-Wno-maybe-uninitialized>
)
