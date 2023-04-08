# CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8.15)
project (PostProcessing)

# Compiling LIBVORO
add_library (libvoro++ voro++-0.4.6/src/cell.cc voro++-0.4.6/src/common.cc voro++-0.4.6/src/container.cc voro++-0.4.6/src/unitcell.cc voro++-0.4.6/src/v_compute.cc voro++-0.4.6/src/c_loops.cc voro++-0.4.6/src/v_base.cc voro++-0.4.6/src/wall.cc voro++-0.4.6/src/pre_container.cc voro++-0.4.6/src/container_prd.cc)
target_compile_options(libvoro++ PRIVATE -ansi -pedantic -O3) 

# Compiling POSTPROCESSING
add_executable (PostProcessing Calculs.cpp CoarseAdvanced.cpp Compress.cpp FDump.cpp FStep.cpp Ids.cpp LDump.cpp Main.cpp Step.cpp CFStep.cpp Coarse.cpp Dump.cpp Filter.cpp Gunzip.cpp LcfDump.cpp LStep.cpp Multisphere.cpp Statistics.cpp Surface.cpp Writing.cpp)
target_include_directories (PostProcessing PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
#add_compile_definitions(VORONOI)

find_library(Z_LIB z)
if (NOT Z_LIB)
  message(FATAL_ERROR "ZLIB is required to compile")
endif()
target_link_libraries(PostProcessing -lz)

find_library(TIFF_LIB tiff)
if(TIFF_LIB)
  add_compile_definitions(TIFF)
  target_link_libraries(PostProcessing -ltiff)
else() 
  message(WARNING "Tiff library not found, compiled without TIFF support")
endif()


find_package(Boost)
if (Boost_FOUND)
    target_compile_options(PostProcessing PUBLIC -DBOOST)
    include_directories(${Boost_INCLUDE_DIRS})
endif(Boost_FOUND)

find_package(Matlab COMPONENTS MAT_LIBRARY)
if (Matlab_FOUND)
  target_compile_options(PostProcessing PUBLIC -DMATLAB)
  include_directories(${Matlab_INCLUDE_DIRS})
  target_link_libraries(PostProcessing ${Matlab_MX_LIBRARY} ${Matlab_MAT_LIBRARY})
endif()






