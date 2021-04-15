set( CMAKE_BUILD_TYPE Release )

set( Boost_USE_STATIC_LIBS ON                                                  )
set( CMAKE_C_COMPILER      "/usr/bin/gcc-10" CACHE FILEPATH "The C compiler"   )
set( CMAKE_CXX_COMPILER    "/usr/bin/g++-10" CACHE FILEPATH "The C++ compiler" )
set( CMAKE_CXX_FLAGS_INIT  " ${CMAKE_CXX_FLAGS} -static "                      )
