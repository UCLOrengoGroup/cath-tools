
Building
========

Compiling from Source
---------------------

You'll need a fairly good (ie recent) C++ compiler (eg on Ubuntu: sudo apt-get install g++). The code's currently being developed against GCC v4.8.2 / v4.9.2 and Clang v3.6.0.

This project currently uses the superposition routine from [bioplib](https://github.com/ACRMGroup/bioplib "Bioplib's GitHub Homepage"), which is included as a Git submodule. There are 2 other dependencies/prerequisites:

 * **CMake** (&ge; v2.8.8  ) : Used to build the software
 * **Boost** (&ge; v1.57.0?) : Used heavily the Boost libraries (both headers and libraries) must be installed on the build machine.

### Boost

    [Ubuntu]
    $ sudo apt-get install libboost-all-dev

### CMake (internal users)

NOTE: You may wish to explicitly specify which cmake to run, particularly if running on a machine that has a local cmake that's more recent than an alternative network cmake, e.g.

    /usr/bin/cmake

(can otherwise produce errors like: ```cmake: error while loading shared libraries: libssl.so.6: cannot open shared object file: No such file or directory```)


Building the Code
-----------------

Once the dependencies are in place, the code can be built with:

    $ cmake -DCMAKE_BUILD_TYPE=RELEASE .
    $ make

NOTE: If you have multiple cores, you can make compiling faster by specifying that it may compile up to N sources simultaneously by appending
      `-j [N]` to the end of the make command.


Compiling Multiple Versions (release/relwithdebinfo/debug)
----------------------------------------------------------

To compile multiple versions, you can put them in individual directories. First make sure that you don't have any
build files in the root directory and then run something like:

    $ mkdir release
    $ cd release
    $ cmake -DCMAKE_BUILD_TYPE=RELEASE ..
    $ cd ..
    $ make -C release

Available types are: release, relwithdebinfo and debug

NOTE: the debug build requires a debug version of Boost built with `_GLIBCXX_DEBUG` enabled otherwise the resulting code will have all sorts of
horrible, non-obvious errors. If you really want to build a debug version without doing this, you can work around this issue by removing all
mentions of `-D_GLIBCXX_DEBUG` from CMakeLists.txt.

Example CMake Commands
----------------------

**clang_debug**

    /usr/bin/cmake -DCMAKE_BUILD_TYPE=DEBUG          -DBUILD_SHARED_LIBS=ON -DBOOST_ROOT=/opt/boost_1_58_0_clang_build -DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DCMAKE_CXX_FLAGS="-stdlib=libc++" ..

**clang_release**

    /usr/bin/cmake -DCMAKE_BUILD_TYPE=RELEASE                               -DBOOST_ROOT=/opt/boost_1_58_0_clang_build -DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DCMAKE_CXX_FLAGS="-stdlib=libc++" ..

**clang_relwithdebinfo**

    /usr/bin/cmake -DCMAKE_BUILD_TYPE=RELWITHDEBINFO                        -DBOOST_ROOT=/opt/boost_1_58_0_clang_build -DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DCMAKE_CXX_FLAGS="-stdlib=libc++" ..

**gcc_debug**

    /usr/bin/cmake -DCMAKE_BUILD_TYPE=DEBUG          -DBUILD_SHARED_LIBS=ON -DBOOST_ROOT=/opt/boost_1_58_0_gcc_build                                                                                                              ..

**gcc_release**

    /usr/bin/cmake -DCMAKE_BUILD_TYPE=RELEASE                               -DBOOST_ROOT=/opt/boost_1_58_0_gcc_build                                                                                                              ..

**gcc_relwithdebinfo**

    /usr/bin/cmake -DCMAKE_BUILD_TYPE=RELWITHDEBINFO                        -DBOOST_ROOT=/opt/boost_1_58_0_gcc_build                                                                                                              ..

Shared Versus Static
----------------

At present, the builds are all completely static. This should be made more configurable in the future so that (at least debug) builds can be run in static mode using the CMake flag `-DBUILD_SHARED_LIBS=ON`.

Building with Clang
-------------------

Basic idea:

    $ cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_CXX_COMPILER=/usr/bin/clang++ .
    $ make

Multiple versions:

    $ mkdir clang_release
    $ cd clang_release
    $ cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_CXX_COMPILER=/usr/bin/clang++ ..
    $ cd ..
    $ make -C clang_release

Building against Clang's C++ library (libc++) rather than GCC's (libstdc++) requires a version of Boost built with Clang and libc++. If you have one, you can use a cmake command like:

    $ cmake -DCMAKE_BUILD_TYPE=RELEASE -DBOOST_ROOT=/opt/boost_1_58_0_clang_build -DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DCMAKE_CXX_FLAGS="-stdlib=libc++" ..

Clang Static Analzyer [For developers]
---------------------------------

Want to search for bugs in the code? These instructions aim to help you run Clang-based static analysis. (TODO: Get this running on a CI server, such as Travis-CI.)

Ensure you have clang installed. Find the analyzer programs `ccc-analyzer` and `c++-analyzer`. For example, on Ubuntu you can do something like:

    dpkg -l | grep clang | awk '{print $2}' | xargs dpkg -L | grep analyzer

Then substitute their locations into the following commands and then run the commands, starting in the root of the cath-tools project (with bioplib already built if you don't want errors from static analysis of bioplib):

    mkdir build-analyze && cd build-analyze
    setenv CCC_CC  clang
    setenv CCC_CXX clang++
    /usr/bin/cmake -DBOOST_ROOT=/opt/boost_1_58_0_clang_build -DCMAKE_C_COMPILER="/usr/share/clang/scan-build-3.6/ccc-analyzer" -DCMAKE_CXX_COMPILER="/usr/share/clang/scan-build-3.6/c++-analyzer" -DCMAKE_CXX_FLAGS="-stdlib=libc++" ..
    scan-build make

To get parallel compilation, you can append ` -j #` to the `scan-build make` (where `#` is the number of threads).

Running the Build Tests
=======================

Once you've built the binaries, run the build tests to sanity check the build. From the root directory of the project, run `build-test` and confirm that all tests pass.

(if you mistakenly run build-test from elsewhere, you'll get lots of `check_required_files_exist` errors).

If your machine has Perl, you can also try running the Perl tests (which includes a run of `build-test` as one of the tests):

 * Set the environment variable CATH_TOOLS_BIN_DIR to the location of the built binaries
 * From the root directory of the project, run `prove -l -v t`
