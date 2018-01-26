Layout:

 * cath-superpose-multi-temp-script - Script to make it a bit easier to do multiple superpositions with cath-superpose
 * cath_tools_git_version.hpp.in    - Template from which CMake can make a header in the build directory with the Git version details
 * CMakeLists.txt                   - CMake configuration file for these sources
 * src_*                            - Individual libraries in separate directories to prevent inappropriate #include-ing
 * uni                              - A big mess of interdependent source code that should be separated out into individual libraries
 * [anything else]                  - Separate libraries that aren't quite separated enough as yet to go into their own separate src_* directories

Ideally, code should migrate: uni -> [anything else] -> src_*.
