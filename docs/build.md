# Build

## Why?

Do you really need to build cath-tools, or can you just use the 64-bit Linux executables from [**DOWNLOADS**](https://github.com/UCLOrengoGroup/cath-tools/releases/latest "The latest CATH Tools release")? Remember to chmod them to be executable (eg `chmod +x cath-ssap`). If you need cath-tools on a different platform (Windows?), please consider creating a [new GitHub issue](https://github.com/UCLOrengoGroup/cath-tools/issues/new "Open a new GitHub issue for cath-tools"); maybe we can work together to set up an automated build.

## Requirements

You'll need a fairly good (ie recent) C++ compiler (eg on Ubuntu: `sudo apt-get install g++`). The code is currently being developed to be buildable with GCC v4.9.2.

~~~no-highlight
git clone https://github.com/UCLOrengoGroup/cath-tools.git
~~~

There are three further dependencies/prerequisites...

### Conan

Conan is used to install some of the C++ dependencies (currently boost and rapidjson)

For Ubuntu:

~~~sh
sudo apt-get install python3-pip
sudo pip3 install --upgrade conan
~~~

### CMake ( &ge; v3.12 )

This is used to build the software.

For Ubuntu:

~~~sh
[Ubuntu]
$ sudo apt-get install cmake
~~~

Downloading and running a recent CMake can sometimes be as simple as:

~~~sh
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

~~~sh
apt-get install gsl-bin libgsl2:amd64 libgsl-dbg:amd64 libgsl-dev
~~~

## Building the dependencies

`cd` into the `cath-tools` directory and then:

~~~sh
mkdir build
conan install --build missing --install-folder build .
cmake -GNinja -B build -S . -D "CMAKE_MODULE_PATH:PATH=${PWD}/build"
~~~

## Building the Code

Once the dependencies are in place, the code can be built with:

~~~sh
cmake -DCMAKE_BUILD_TYPE=Release -B build -S . -D"CMAKE_MODULE_PATH:PATH=${PWD}/build"
make
~~~

NOTE:

* If you have multiple cores, you can make compiling faster by specifying that it may compile up to N sources simultaneously by appending `-j [N]` to the end of the make command.
* If your system does not have the Gnu Scientific Library available as a static library on Linux, you should pass `-D USE_STATIC_GSL:BOOL=FALSE` to cmake.

## Running the Build Tests

Once you've built the binaries, run the build tests to sanity check the build. Run `bin/build-test` in the build directory and confirm that all tests pass. (Alternatively, `cd` into the build directory and run `ctest`, with `-j [N]` for parallelism.)

If your machine has Perl, you can also try running the Perl tests (which includes a run of `build-test` as one of the tests):

* Set the environment variable `CATH_TOOLS_BIN_DIR` to the location of the built binaries
* Make sure your Perl has access to the required dependencies (if you have [cpanm](https://metacpan.org/pod/distribution/Menlo/script/cpanm-menlo) installed then try `cpanm --installdeps ./perl`)
* From the root directory of the project, run `prove -l -v perl/t`

Assuming you have already built the binaries (in the project root):

~~~no-highlight
CATH_TOOLS_BIN_DIR=. prove -l -v ./perl/t
~~~

## Building on CentOS 6

Install these packages as root:

~~~sh
yum install bzip2-devel cmake git
yum install centos-release-scl-rh
yum install devtoolset-3-gcc devtoolset-3-gcc-c++
~~~

Then ssh to the build machine as yourself, find some working directory with at least a few Gb of free space, and then substitute it in for `WHATEVER_YOU_HAVE_CHOSEN_AS_YOUR_BUILD_ROOT_DIRECTORY` in these commands:

Setup:

~~~sh
scl enable devtoolset-3 bash
export BUILD_ROOT_DIR=WHATEVER_YOU_HAVE_CHOSEN_AS_YOUR_BUILD_ROOT_DIRECTORY
~~~

Build Boost:

~~~sh
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

~~~sh
cd ${BUILD_ROOT_DIR}/
git clone https://github.com/UCLOrengoGroup/cath-tools.git

mkdir -p ${BUILD_ROOT_DIR}/cath-tools/gcc_relwithdebinfo/
cd       ${BUILD_ROOT_DIR}/cath-tools/gcc_relwithdebinfo/
/usr/bin/cmake -DCMAKE_BUILD_TYPE=RELEASE -DBOOST_ROOT=${BUILD_ROOT_DIR}/boost_1_60_0_build ..
cd       ${BUILD_ROOT_DIR}/cath-tools/
make -C gcc_relwithdebinfo -k -j2
ls -l ${BUILD_ROOT_DIR}/cath-tools/gcc_relwithdebinfo/
~~~
