# Build

## Why?

Do you really need to build cath-tools, or can you just use the 64-bit Linux executables from [**DOWNLOADS**](https://github.com/UCLOrengoGroup/cath-tools/releases/latest "The latest CATH Tools release")? Remember to chmod them to be executable (eg `chmod +x cath-ssap`). If you need cath-tools on a different platform (Windows? Mac?), please consider creating a [new GitHub issue](https://github.com/UCLOrengoGroup/cath-tools/issues/new "Open a new GitHub issue for cath-tools"); maybe we can work together to set up an automated build.

## Requirements


You'll need a fairly good (ie recent) C++ compiler (eg on Ubuntu: `sudo apt-get install g++`). The code is currently being developed against GCC v4.9.2 and Clang v3.6.0.

~~~no-highlight
git clone https://github.com/UCLOrengoGroup/cath-tools.git
~~~

There are three further dependencies/prerequisites...

### Boost ( v1.60.0 / v1.61.0 )

This is used heavily throughout the code. Both headers and compiled library files are needed.

~~~no-highlight
[Ubuntu]
$ sudo apt-get install libboost-all-dev
~~~

### CMake ( &ge; v3.2 )

This is used to build the software.

~~~no-highlight
[Ubuntu]
$ sudo apt-get install cmake
~~~

Downloading and running a recent CMake can sometimes be as simple as:

~~~
wget "https://cmake.org/files/v3.9/cmake-3.9.4-Linux-x86_64.tar.gz"
tar -zxvf cmake-3.9.4-Linux-x86_64.tar.gz cmake-3.9.4-Linux-x86_64/bin/cmake
tar -zxvf cmake-3.9.4-Linux-x86_64.tar.gz cmake-3.9.4-Linux-x86_64/share/cmake-3.9/Modules
cmake-3.9.4-Linux-x86_64/bin/cmake --version
~~~

*NOTE for UCL SMB users*: If running cmake gives you errors like:

~~~no-highlight
cmake: error while loading shared libraries: libssl.so.6: cannot open shared object file: No such file or directory
~~~

You probably mixing local/network binaries and libraries so try explicitly running /usr/bin/cmake.

### GSL

The [GNU Scientific Library](https://www.gnu.org/software/gsl/) is used for its Singular Value Decomposition function (`gsl_linalg_SV_decomp()`).

`apt-get install gsl-bin libgsl2:amd64 libgsl-dbg:amd64 libgsl-dev`

Building the Code
-----------------

Once the dependencies are in place, the code can be built with:

~~~no-highlight
$ cmake -DCMAKE_BUILD_TYPE=RELEASE .
$ make
~~~

NOTE: If you have multiple cores, you can make compiling faster by specifying that it may compile up to N sources simultaneously by appending `-j [N]` to the end of the make command.

With default Ubuntu config, this will build with GCC against GCC's standard library (libstdc++). If you instead want to build with Clang against Clang's standard C++ library (libc++), you'll need a version of Boost built with Clang and libc++. If you have one, then you can build with Clang by adding `-DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_CXX_COMPILER=/usr/bin/clang++` to the CMake command and build against libc++ by adding `-DBOOST_ROOT=ROOT_DIR_OF_YOUR_CLANG_BUILD_OF_BOOST` and `-DCMAKE_CXX_FLAGS="-stdlib=libc++"`.

To build against Clang's C++ library (libc++) rather than GCC's (libstdc++), you'll need a version of Boost built with Clang and libc++. If you have one, you can use a cmake command like:

~~~no-highlight
$ cmake -DCMAKE_BUILD_TYPE=RELEASE -DBOOST_ROOT=/opt/boost_1_60_0_clang_build -DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DCMAKE_CXX_FLAGS="-stdlib=libc++" ..
~~~

# Running the Build Tests

Once you've built the binaries, run the build tests to sanity check the build. From the root directory of the project, run `build-test` and confirm that all tests pass.

(if you mistakenly run build-test from elsewhere, you'll get lots of `check_required_files_exist` errors).

If your machine has Perl, you can also try running the Perl tests (which includes a run of `build-test` as one of the tests):

 * Set the environment variable `CATH_TOOLS_BIN_DIR` to the location of the built binaries
 * From the root directory of the project, run `prove -l -v t`

# Building on CentOS 6

Install these packages as root:

~~~no-highlight
yum install bzip2-devel cmake git
yum install centos-release-scl-rh
yum install devtoolset-3-gcc devtoolset-3-gcc-c++
~~~

Then ssh to the build machine as yourself, find some working directory with at least a few Gb of free space, and then substitute it in for `WHATEVER_YOU_HAVE_CHOSEN_AS_YOUR_BUILD_ROOT_DIRECTORY` in these commands:

Setup:

~~~no-highlight
scl enable devtoolset-3 bash
export BUILD_ROOT_DIR=WHATEVER_YOU_HAVE_CHOSEN_AS_YOUR_BUILD_ROOT_DIRECTORY
~~~

Build Boost:

~~~no-highlight
mkdir -p ${BUILD_ROOT_DIR}/boost_1_60_0_build/{include,lib}
cd ${BUILD_ROOT_DIR}/
wget "http://sourceforge.net/projects/boost/files/boost/1.60.0/boost_1_60_0.tar.gz"
tar -zxvf boost_1_60_0.tar.gz
cd ${BUILD_ROOT_DIR}/boost_1_60_0/
./bootstrap.sh --prefix=${BUILD_ROOT_DIR}/boost_1_60_0_build
./b2 -j2         --layout=tagged variant=release
./b2 -j2 install --layout=tagged variant=release
~~~

Build cath-tools:

~~~no-highlight
cd ${BUILD_ROOT_DIR}/
git clone https://github.com/UCLOrengoGroup/cath-tools.git

mkdir -p ${BUILD_ROOT_DIR}/cath-tools/gcc_relwithdebinfo/
cd       ${BUILD_ROOT_DIR}/cath-tools/gcc_relwithdebinfo/
/usr/bin/cmake -DCMAKE_BUILD_TYPE=RELEASE -DBOOST_ROOT=${BUILD_ROOT_DIR}/boost_1_60_0_build ..
cd       ${BUILD_ROOT_DIR}/cath-tools/
make -C gcc_relwithdebinfo -k -j2
ls -l ${BUILD_ROOT_DIR}/cath-tools/gcc_relwithdebinfo/
~~~

# Building for CentOS 5

It's non-trivial got get a modern enough C++ compiler on CentOS 5. The standard downloads of Clang 3.6.2 don't work and it's likely that older versions won't compile the cath-tools code. (Though, for reference, CMake 3.6.3 seems to work cleanly on CentOS 5 from the download.)

Instead, it's possible to build a version on CentOS 6 that's completely statically linked and this should run on CentOS5. You can find such executables in the [v0.16.2](https://github.com/UCLOrengoGroup/cath-tools/releases/tag/v0.16.2) release.

So far, this has only been done in a fairly hacky way, by installing the relevant static libraries (`yum install glibc-static`), modifying the `CMakeLists.txt` file:

~~~diff
diff --git a/CMakeLists.txt b/CMakeLists.txt
index 3badace..06520a3 100644
--- a/CMakeLists.txt
+++ b/CMakeLists.txt
@@ -53,14 +53,18 @@ foreach (loop_var RANGE ${GSL_LIBRARIES})
  list(APPEND GSL_DYN_LINK_FLAGS "-l${loop_var}")
endforeach(loop_var)
 
-if( BUILD_SHARED_LIBS )
-       add_definitions( -DBOOST_ALL_DYN_LINK )
-       add_definitions( -DBOOST_LOG_DYN_LINK )
-       SET( GSL_LIB_SUFFIX ${GSL_LIBRARIES} )
-else()
+SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
+SET(BUILD_SHARED_LIBRARIES OFF)
+SET(CMAKE_EXE_LINKER_FLAGS "-static")
+
+#if( BUILD_SHARED_LIBS )
+#      add_definitions( -DBOOST_ALL_DYN_LINK )
+#      add_definitions( -DBOOST_LOG_DYN_LINK )
+#      SET( GSL_LIB_SUFFIX ${GSL_LIBRARIES} )
+#else()
       SET( Boost_USE_STATIC_LIBS ON )
       SET( GSL_LIB_SUFFIX "${GSL_STATIC_LIB}" "${GSLCBLAS_STATIC_LIB}" )
-endif()
+#endif()
 
if ( ${LSB_RELEASE_CODE} STREQUAL "yakkety" )
       SET( GSL_LIB_SUFFIX ${GSL_LIBRARIES} )
~~~

...and manually repeating the link commands with all `-Wl,-Bdynamic` flags removed.

This process can doubtless be improved with a bit of work.
