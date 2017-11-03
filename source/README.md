Layout:

 * cath-superpose-multi-temp-script - Script to make it a bit easier to do multiple superpositions with cath-superpose
 * cath_tools_git_version.hpp.in    - Template from which CMake can make a header in the build directory with the Git version details
 * CMakeLists.txt                   - CMake configuration file for these sources
 * src_*                            - Individual libraries in separate directories to prevent inappropriate #include-ing
 * uni                              - A big mess of interdependent source code that should be separated out into libraries
 * [anything else]                    - Library that isn't quite separated enough as yet to go into its own separate src_* directory

Ideally, code should migrate from uni / [anything else] into src_*.
