
Development
===========

This page contains extra information that may be useful for anyone working on the development of cath-tools. For information on how to set up a standard build of the code, see [Build](build).

Doxygen Code Documentation
--------------------------

Much of the code is documented inline with Doxygen. To view it, install doxygen, run `doxygen` in the cath-tools root directory and then view `doxygen_documentation/html/index.html` in a browser.

<!-- TODO: Set up Doxygen build? -->


Compiling Multiple Versions (release/relwithdebinfo/debug)
----------------------------------------------------------

To compile multiple versions, you can put them in individual directories. First make sure that you don't have any
build files in the root directory and then run something like:


~~~~~no-highlight
$ mkdir release
$ cd release
$ cmake -DCMAKE_BUILD_TYPE=RELEASE ..
$ cd ..
$ make -C release
~~~~~

Available types are: release, relwithdebinfo and debug

NOTE: the debug build requires a debug version of Boost built with `_GLIBCXX_DEBUG` enabled otherwise the resulting code will have all sorts of
horrible, non-obvious errors. If you really want to build a debug version without doing this, you can work around this issue by removing all
mentions of `-D_GLIBCXX_DEBUG` from CMakeLists.txt.

Example CMake Commands
----------------------

**clang_debug**

~~~~~no-highlight
/usr/bin/cmake -DCMAKE_BUILD_TYPE=DEBUG          -DBUILD_SHARED_LIBS=ON -DBOOST_ROOT=/opt/boost_1_58_0_clang_build -DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DCMAKE_CXX_FLAGS="-stdlib=libc++" ..
~~~~~

**clang_release**

~~~~~no-highlight
/usr/bin/cmake -DCMAKE_BUILD_TYPE=RELEASE                               -DBOOST_ROOT=/opt/boost_1_58_0_clang_build -DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DCMAKE_CXX_FLAGS="-stdlib=libc++" ..
~~~~~

**clang_relwithdebinfo**

~~~~~no-highlight
/usr/bin/cmake -DCMAKE_BUILD_TYPE=RELWITHDEBINFO                        -DBOOST_ROOT=/opt/boost_1_58_0_clang_build -DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DCMAKE_CXX_FLAGS="-stdlib=libc++" ..
~~~~~

**gcc_debug**

~~~~~no-highlight
/usr/bin/cmake -DCMAKE_BUILD_TYPE=DEBUG          -DBUILD_SHARED_LIBS=ON -DBOOST_ROOT=/opt/boost_1_58_0_gcc_build                                                                                                              ..
~~~~~

**gcc_release**

~~~~~no-highlight
/usr/bin/cmake -DCMAKE_BUILD_TYPE=RELEASE                               -DBOOST_ROOT=/opt/boost_1_58_0_gcc_build                                                                                                              ..
~~~~~

**gcc_relwithdebinfo**

~~~~~no-highlight
/usr/bin/cmake -DCMAKE_BUILD_TYPE=RELWITHDEBINFO                        -DBOOST_ROOT=/opt/boost_1_58_0_gcc_build                                                                                                              ..
~~~~~

Consider using ninja instead of make
------------------------------------

If you're developing, consider using ninja instead of make. It's a drop-in replacement that's much quicker to get started, which can make a big difference on incremental builds. To use ninja, install a suitable package and add `-GNinja` to your CMake commands.

Shared Versus Static
----------------

At present, the builds are all completely static. This should be made more configurable in the future so that (at least debug) builds can be run in static mode using the CMake flag `-DBUILD_SHARED_LIBS=ON`.

Clang Static Analzyer
--------------------------------------

Want to search for bugs in the code? These instructions aim to help you run Clang-based static analysis. (TODO: Get this running on a CI server, such as Travis-CI.)

Ensure you have clang installed. Find the analyzer programs `ccc-analyzer` and `c++-analyzer`. For example, on Ubuntu you can do something like:

~~~~~no-highlight
dpkg -l | grep clang | awk '{print $2}' | xargs dpkg -L | grep analyzer
~~~~~

Then substitute their locations into the following commands and then run the commands, starting in the root of the cath-tools project (with bioplib already built if you don't want errors from static analysis of bioplib):


~~~~~no-highlight
mkdir build-analyze && cd build-analyze
setenv CCC_CC  clang
setenv CCC_CXX clang++
# For Clang 3.6
/usr/bin/cmake -DBOOST_ROOT=/opt/boost_1_58_0_clang_build -DCMAKE_C_COMPILER="/usr/share/clang/scan-build-3.6/ccc-analyzer"         -DCMAKE_CXX_COMPILER="/usr/share/clang/scan-build-3.6/c++-analyzer"         -DCMAKE_CXX_FLAGS="-stdlib=libc++" ..
# For Clang 3.8
/usr/bin/cmake -DBOOST_ROOT=/opt/boost_1_58_0_clang_build -DCMAKE_C_COMPILER="/usr/share/clang/scan-build-3.8/libexec/ccc-analyzer" -DCMAKE_CXX_COMPILER="/usr/share/clang/scan-build-3.8/libexec/c++-analyzer" -DCMAKE_CXX_FLAGS="-stdlib=libc++" ..
scan-build make
~~~~~

To get parallel compilation, you can append ` -j #` to the `scan-build make` (where `#` is the number of threads).



Checking headers compile independently
--------------------------------------

Clang:

~~~~~no-highlight
find source -iname '*.hpp' | sort | grep third_party_code -v | xargs -P 4 -I VAR clang++ -x c++ -DBOOST_LOG -std=c++1y -stdlib=libc++ -W -Wall -Werror -Wextra -Wno-unused-const-variable -Wno-unused-local-typedef -Wsign-compare -Wcast-qual -Wconversion -Wnon-virtual-dtor -pedantic -ftemplate-backtrace-limit=0 -c -o /tmp/.comp_clang.dummy.header.o -isystem /opt/boost_1_58_0_clang_build/include -I source VAR
~~~~~

GCC:

~~~~~no-highlight
find source -iname '*.hpp' | sort | grep third_party_code -v | xargs -P 4 -I VAR g++     -x c++ -DBOOST_LOG -std=c++1y                -W -Wall -Werror -Wextra -Wno-unused-const-variable -Wno-unused-local-typedef -Wsign-compare -Wcast-qual -Wconversion -Wnon-virtual-dtor -pedantic -ftemplate-backtrace-limit=0 -c -o /tmp/.comp_gcc.dummy.header.o   -isystem /opt/boost_1_58_0_gcc_build/include   -I source VAR
~~~~~


Fixing trailing namespace comments
----------------------------------

~~~~~no-highlight
find source -iname '*.hpp' | sort | grep third_party_code -v | xargs -P 4 -I VAR clang-tidy -fix -checks=llvm-namespace-comment VAR -- -x c++ -std=c++1y -isystem /opt/boost_1_58_0_clang_build/include -I source
~~~~~



Fixing header guards
--------------------

Whilst the header guards start with one underscore (contrary to clang-tidy's preference):

~~~~~no-highlight
find $PWD/source -type f -iname '*.hpp' | sort | grep -vw 'third_party_code' | sed 's/\.hpp$//g' | xargs -I VAR -P 8 ln -s VAR.hpp VAR.h
find $PWD/source -type l -iname '*.h'   | sort | xargs -I VAR -P 4 clang-tidy -fix -checks=llvm-header-guard VAR -- -x c++ -std=c++1y -isystem /opt/boost_1_58_0_clang_build/include -I source
find $PWD/source -type l -iname '*.h'   | sort | xargs rm -f
find $PWD/source -type f -iname '*.h'   | sort | grep -vw 'third_party_code' | xargs sed -i 's/#ifndef __CATH_TOOLS_SOURCE_/#ifndef _CATH_TOOLS_SOURCE_/g'
find $PWD/source -type f -iname '*.h'   | sort | grep -vw 'third_party_code' | xargs sed -i 's/#define __CATH_TOOLS_SOURCE_/#define _CATH_TOOLS_SOURCE_/g'
find $PWD/source -type f -iname '*.h'   | sort | grep -vw 'third_party_code' | sed 's/\.h$//g' | xargs -I VAR diff -C1 VAR.hpp VAR.h
find $PWD/source -type f -iname '*.h'   | sort | grep -vw 'third_party_code' | sed 's/\.h$//g' | xargs -I VAR diff -C1 VAR.hpp VAR.h -q | awk '{print $4}' | xargs -I VAR mv VAR VARpp
find $PWD/source -type f -iname '*.h'   | sort | grep -vw 'third_party_code' | xargs rm -f
~~~~~

...or once all the header guards have been changed to start with two underscores...

~~~~~no-highlight
find $PWD/source -type f -iname '*.hpp' | sort | grep -vw 'third_party_code' | sed 's/\.hpp$//g' | xargs -I VAR -P 8 ln -s VAR.hpp VAR.h
find $PWD/source -type l -iname '*.h'   | sort | xargs -I VAR -P 4 clang-tidy -fix -checks=llvm-header-guard VAR -- -x c++ -std=c++1y -isystem /opt/boost_1_58_0_clang_build/include -I source
find $PWD/source -type l -iname '*.h'   | sort | xargs rm -f
find $PWD/source -type f -iname '*.h'   | sort | grep -vw 'third_party_code' | sed 's/\.h$//g' | xargs -I VAR diff -C1 VAR.hpp VAR.h
find $PWD/source -type f -iname '*.h'   | sort | grep -vw 'third_party_code' | sed 's/\.h$//g' | xargs -I VAR mv       VAR.h   VAR.hpp
~~~~~


Assessing with all clang-tidy checks
--

Prefer using more recent version of clang-tidy; it's changing fast.

Would like to use:
 * `misc-use-override` but (in version 3.6.2), it erroneously fires for method declarations that *do* use override
 * `google-readability-function` but (in version 3.6.2), it fires for declarations, for which I don't want to always name all params

find source -iname '*.?pp' | sort | grep third_party_code -v | xargs -P 4 -I VAR clang-tidy -checks=llvm-include-order VAR -- -std=c++1y -isystem /opt/boost_1_58_0_clang_build/include -I source


~~~~~no-highlight
clang-tidy '-checks=*,-llvm-header-guard,-llvm-namespace-comment,-google-readability-namespace-comments,-google-build-using-namespace,-misc-use-override,-google-readability-function' -list-checks --
clang-tidy '-checks=*,-llvm-header-guard,-llvm-namespace-comment,-google-readability-namespace-comments,-google-build-using-namespace,-misc-use-override,-google-readability-function' -dump-config --
find source -iname '*.?pp' | sort | grep third_party_code -v | xargs -P 4 -I VAR clang-tidy VAR '-checks=*,-llvm-header-guard,-llvm-namespace-comment,-google-readability-namespace-comments,-google-build-using-namespace,-misc-use-override,-google-readability-function' -- -std=c++1y -isystem /opt/boost_1_58_0_clang_build/include -I source

find source -iname '*.?pp' | sort | grep third_party_code -v | xargs -P 4 -I VAR clang-tidy VAR '-checks=*,-llvm-header-guard' - -- -std=c++1y -isystem /opt/boost_1_58_0_clang_build/include -I source
~~~~~


Dumping class/struct memory layouts
-----------------------------------

~~~~~no-highlight
clomp source/file/pdb/pdb_atom.hpp -Xclang -fdump-record-layouts > /tmp/clang_pdb_atom_layout.txt
grep pdb_atom -A60 /tmp/clang_pdb_atom_layout.txt
~~~~~


Simplifying hmmsearch output files for cath-resolve-hits
--------------------------------------------------------

~~~~~no-highlight
grep -A500 '>> O27798_698e555d8cb2c0ea3979641198e527bf' temp_0.hmmsearch.evalcoff0.001 | tr '\n' '@' | sed 's/@>>/\n>>/g' | sed 's/@Internal pipeline statistics summary/\nInternal pipeline statistics summary/g' | grep -v 'Internal pipeline statistics summary' | grep O27798_698e555d8cb2c0ea3979641198e527bf | tr '@' '\n' > bob
echo '\n\nInternal pipeline statistics summary:\n[ok]\n' >> bob
~~~~~

