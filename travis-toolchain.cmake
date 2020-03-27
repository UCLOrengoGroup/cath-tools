
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
