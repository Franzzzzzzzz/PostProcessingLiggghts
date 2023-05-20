# CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8.15)
project (PostProcessing)

# ======================== Compiling LIBVORO
add_library (voro++ voro++-0.4.6/src/cell.cc voro++-0.4.6/src/common.cc voro++-0.4.6/src/container.cc voro++-0.4.6/src/unitcell.cc voro++-0.4.6/src/v_compute.cc voro++-0.4.6/src/c_loops.cc voro++-0.4.6/src/v_base.cc voro++-0.4.6/src/wall.cc voro++-0.4.6/src/pre_container.cc voro++-0.4.6/src/container_prd.cc)
target_compile_options(voro++ PRIVATE -ansi -pedantic -O3) 

#========================= Handling internationalisation
find_package (Intl REQUIRED)
if (Intl_FOUND)
    message(STATUS "Internationalization (i18n) found:")
    message(STATUS " INTL_INCLUDE_DIRS: ${Intl_INCLUDE_DIRS}")
    message(STATUS " INTL_LIBRARIES: ${Intl_LIBRARIES}")
    message(STATUS " Version: ${Intl_VERSION}")
    include_directories(${Intl_INCLUDE_DIRS})
    link_directories(${Intl_LIBRARY_DIRS})
else ()
    message(STATUS "Internationalization (i18n) Not found!")
endif ()

find_package(Gettext REQUIRED)
if (Gettext_FOUND)
    message(STATUS "Gettext found:")
    message(STATUS " Version: ${GETTEXT_VERSION_STRING}")
else ()
    message(STATUS "Gettext Not found!")
endif ()

find_program(GETTEXT_XGETTEXT_EXECUTABLE xgettext)

if (GETTEXT_XGETTEXT_EXECUTABLE)

    #message(DEBUG " xgettext: ${GETTEXT_XGETTEXT_EXECUTABLE}")
    file(GLOB_RECURSE CPP_FILES RELATIVE ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/*.cpp)
    add_custom_target(
        pot-update
        COMMENT "pot-update: Done."
        DEPENDS ${CMAKE_SOURCE_DIR}/locale/postproc.pot
    )
    add_custom_command(
        TARGET pot-update
        PRE_BUILD
        COMMAND
            ${GETTEXT_XGETTEXT_EXECUTABLE}
            --from-code=utf-8
            --c++
            --force-po
            --output=${CMAKE_SOURCE_DIR}/locale/postproc.pot
            --keyword=_
            --width=80
            ${CPP_FILES}
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "pot-update: Pot file generated: ${CMAKE_SOURCE_DIR}/locale/postproc.pot"
    )

endif (GETTEXT_XGETTEXT_EXECUTABLE)

if (GETTEXT_MSGMERGE_EXECUTABLE)
    #message(DEBUG " msgmerge: ${GETTEXT_MSGMERGE_EXECUTABLE}")
    add_custom_target(
        pot-merge
        COMMENT "pot-merge: Done."
        DEPENDS ${CMAKE_SOURCE_DIR}/locale/postproc.pot
    )

    file(GLOB PO_FILES ${CMAKE_SOURCE_DIR}/locale/*/postproc.po)
    message(TRACE " PO_FILES: ${PO_FILES}")

    foreach(PO_FILE IN ITEMS ${PO_FILES})
        message(STATUS " Adding msgmerge for: ${PO_FILE}")
        add_custom_command(
            TARGET pot-merge
            PRE_BUILD
            COMMAND
                ${GETTEXT_MSGMERGE_EXECUTABLE} -U ${PO_FILE}
                ${CMAKE_SOURCE_DIR}/locale/postproc.pot
            COMMENT "pot-merge: ${PO_FILE}"
        )
    endforeach()
endif (GETTEXT_MSGMERGE_EXECUTABLE)

if (GETTEXT_MSGFMT_EXECUTABLE)

    #message(DEBUG " msgmerge: ${GETTEXT_MSGFMT_EXECUTABLE}")
    file(GLOB PO_LANGS LIST_DIRECTORIES true ${CMAKE_SOURCE_DIR}/locale/*)
    #message(TRACE " PO_LANGS: ${PO_LANGS}")

    add_custom_target(
        po-compile
        COMMENT "po-compile: Done."
    )

    foreach(PO_LANG IN ITEMS ${PO_LANGS})
        if(IS_DIRECTORY ${PO_LANG})
        message(STATUS " Adding msgfmt for: ${PO_LANG}")
        add_custom_command(
            TARGET po-compile
            PRE_BUILD
            COMMAND
                ${GETTEXT_MSGFMT_EXECUTABLE}
                --output-file=LC_MESSAGES/postproc.mo
                postproc.po
            WORKING_DIRECTORY ${PO_LANG}
            COMMENT "po-compile: ${PO_LANG}"
        )
        endif()
    endforeach()

endif (GETTEXT_MSGFMT_EXECUTABLE)

set(language "" CACHE STRING "Force a locale/language")
set_property(CACHE language PROPERTY STRINGS "fr_FR.UTF-8" "en_GB.UTF-8")


# ======================== Compiling POSTPROCESSING
link_directories(${CMAKE_CURRENT_BINARY_DIR})
add_executable (PostProcessing Calculs.cpp CoarseAdvanced.cpp Compress.cpp FDump.cpp FStep.cpp Ids.cpp LDump.cpp Main.cpp Step.cpp CFStep.cpp Coarse.cpp Dump.cpp Filter.cpp Gunzip.cpp LcfDump.cpp LStep.cpp Multisphere.cpp Statistics.cpp Surface.cpp Writing.cpp)
#target_include_directories (PostProcessing PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
add_compile_definitions(LOCALEPATH="${CMAKE_SOURCE_DIR}/locale/")
add_compile_definitions(VORONOI)
add_compile_definitions(LANGUAGE="${language}")
target_include_directories(PostProcessing PUBLIC ${CMAKE_SOURCE_DIR}/voro++-0.4.6/src/)
target_link_libraries(PostProcessing -lvoro++)

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






