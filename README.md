README
======

Cloning the cath-tools GitHub Repo
---

This project includes [bioplib](https://github.com/ACRMGroup/bioplib "Bioplib's GitHub Homepage") as a submodule. To ensure the bioplib directory gets populated, clone with:

    git clone --recursive https://github.com/UCLOrengoGroup/cath-tools.git

...or if you've already cloned, then use:

    git submodule update --init --recursive

Building
========

Compiling from Source
---------------------

You'll need a fairly good (ie recent) C++ compiler (eg on Ubuntu: sudo apt-get install g++). The code's currently being developed against GCC v4.8.2 / v4.9.2 and Clang v3.6.0.

cath_binaries is compiled using the CMake build system.

There are 3 dependencies/prerequisites:

 * **CMake** (&ge; v2.8.8  ) : Used to build the software
 * **Boost** (&ge; v1.57.0?) : Used heavily the Boost libraries (both headers and libraries) must be installed on the build machine.
 * **bioplib**             : Currently used for its superposition routine.
   The location is specified with the environment variable BUILDROOT, which must point to a directory that
   contains a directory bioplib/trunk that in turn contains the headers and the two compiled libraries, libbiop.a and libgen.a.

### Boost

    [Ubuntu]
    $ sudo apt-get install libboost-all-dev

### bioplib (external users)

At present, there's an issue that stops us using bioplib &ge; v3.0.
For now, we recommend using v2.1.2.

Download the compressed archive from [v2.1.2 on GitHub](https://github.com/ACRMGroup/bioplib/archive/V2.1.2.tar.gz)
via [the project page](http://www.bioinf.org.uk/software/bioplib/), decompress the archive and cd into the new directory.

Then:

    $ tar -zxvf bioplib-2.1.2.tar.gz
    $ cd bioplib-2.1.2/src
    $ make
    $ mkdir bioplib
    $ ln -s .. bioplib/trunk
    $ setenv BUILDROOT $PWD

If attempting to use &ge; v3.0, be sure to comment out the line in the Makefile that switches on use of libxml2:

    # Include libxml2
    # Note: xml2-config is installed with libxml2
    #       Link to libxml2 with -lxml2
    COPT := $(COPT) -D XML_SUPPORT $(shell xml2-config --cflags)

Otherwise the build against the bioplib libraries will be broken by unresolved dependencies to libxml2.

### bioplib (internal users)

Get a local checkout of the bioplib library and build it

    $ svn co svn://orengosvn/bioplib ~/svn/bioplib
    $ cd ~/svn/bioplib
    $ make

Make this visible to the build by adding the `BUILDROOT` to your environment

    # ~/.bashrc
    $ export BUILDROOT=~/svn

    # ~/.cshrc
    $ setenv BUILDROOT ~/svn


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

Running the Build Tests
=======================

Once you've built the binaries, run the build tests to sanity check the build. From the root directory of the project, run `build-test` and confirm that all tests pass.

(if you mistakenly run build-test from elsewhere, you'll get lots of `check_required_files_exist` errors).

If your machine has Perl, you can also try running the Perl tests (which includes a run of `build-test` as one of the tests):
 * Set the environment variable CATH_BINARIES_BIN_DIR to the location of the built binaries
 * From the root directory of the project, run `prove -l -v t`

Quickstart: Running SSAP
========================

The binary for SSAP is called `cath-ssap` and it'll tell you its usage if you run `cath-ssap --help`.

Most external users will want to use the `--protein_source_files PDB_DSSP_SEC` option (though this may become the default in the future).

You will need to tell it where to find the PDB, wolf and sec files. We recommend you use the path options to manage your list of directories to search. This can be done with command line options, but it may be easier to add suitable environment variables to your profile:

    CATH_BINARIES_PDB_PATH  .:/global/data/directories/pdb
    CATH_BINARIES_DSSP_PATH .:/global/data/directories/dssp
    CATH_BINARIES_SEC_PATH  .:/global/data/directories/sec
