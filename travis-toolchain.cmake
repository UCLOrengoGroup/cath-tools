
# As of Boost 1.70, the finding becomes something of a mess with tangled
# conflicts between the find module file and the config file, eg see:
#
# https://github.com/boostorg/boost_install/issues/12
#
# This caused the Travis Mac build to start failing with problems like those described here:
#
# https://gitlab.kitware.com/cmake/cmake/issues/19714
#
# So the following line disables the config file mode and reverts behaviour to as before.
set( Boost_NO_BOOST_CMAKE 1 )

add_compile_options(
	-W
	-Wall
	-Werror
	-Wextra
	-Wcast-qual
	-Wconversion
	-Wformat=2
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
