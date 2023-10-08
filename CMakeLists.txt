project(PolyView)

cmake_minimum_required (VERSION 3.15)

# if ("${CONDA_DEPS_DIR}" STREQUAL "")
#   message(FATAL_ERROR "You need to set CONDA_DEPS_DIR")
# endif()

# This block is for Qt
# Tell CMake to run moc when necessary:
set(CMAKE_AUTOMOC ON)
# As moc files are generated in the binary dir, tell CMake
# to always look for includes there:
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Widgets finds its own dependencies (QtGui and QtCore).
find_package(Qt5Widgets REQUIRED)
# The Qt5Widgets_INCLUDES also includes the include directories for
# dependencies QtCore and QtGui
include_directories(${Qt5Widgets_INCLUDES})
# We need add -DQT_WIDGETS_LIB when using QtWidgets in Qt 5.
add_definitions(${Qt5Widgets_DEFINITIONS})
# Executables fail to build with Qt 5 in the default configuration
# without -fPIE. We add that here.

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -fPIC " ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS})

# Compile with OPENMP
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -fopenmp -O2 -DPOLYVIEW_USE_OPENMP" ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS})

find_package(OpenGL REQUIRED)

# The Windows compiler does not understand the line below
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11")

# Find all sources which are not tests
# TODO(oalexan1): Make the code below more compact
file(GLOB SOURCES "geom/*.cpp" "geom/*.h" "gui/*.cpp" "gui/*.h")
file(GLOB TESTS "geom/*test*.cpp" "gui/*test*.cpp")
file(GLOB MAIN_RPOG "gui/*mainProg*")
foreach(test ${TESTS})
  list(REMOVE_ITEM SOURCES ${test})
endforeach()
foreach(test ${MAIN_RPOG})
  list(REMOVE_ITEM SOURCES ${test})
endforeach()

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/geom"
    "${CMAKE_CURRENT_SOURCE_DIR}/gui"
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CONDA_DEPS_DIR}/include")

# Static linking with polyview libs
add_library(polyview_lib STATIC ${SOURCES}) 
target_link_libraries(polyview_lib Qt5::Widgets)

add_executable(polyview "gui/mainProg.cpp")
target_link_libraries(polyview polyview_lib Qt5::Widgets
    ${OPENGL_opengl_LIBRARY}
    ${OPENGL_egl_LIBRARY}
    ${OPENGL_glu_LIBRARY}
    ${OPENGL_glx_LIBRARY}
    #/lib64/libGL.so.1 /home/oalexan1/miniconda3/envs/isis6/lib/libX11.so.6 /home/oalexan1/miniconda3/envs/isis6/lib/libXext.so.6 /lib64/libGLdispatch.so.0 /home/oalexan1/miniconda3/envs/isis6/lib/././libicudata.so.58 /home/oalexan1/miniconda3/envs/isis6/lib/./libxcb.so.1 /home/oalexan1/miniconda3/envs/isis6/lib/././libXau.so.6  /home/oalexan1/miniconda3/envs/isis6/lib/././libXdmcp.so.6 
)

# If conda set the prefix, install there
if (NOT ("$ENV{PREFIX}" STREQUAL ""))
    set(CMAKE_INSTALL_BINDIR $ENV{PREFIX}/bin)
endif()

# TODO(oalexan1): Likely conda will take care of the rpaths, so the block
# below may not be necessary if building with conda.
if (NOT ("${CONDA_DEPS_DIR}" STREQUAL ""))
    if (APPLE)
        set_target_properties(polyview PROPERTIES
            INSTALL_RPATH "@loader_path;${CONDA_DEPS_DIR}/lib")
    elseif(UNIX) # Unix which is not Apple
        set_target_properties(polyview PROPERTIES
            INSTALL_RPATH "$ORIGIN:$ORIGIN/../lib:${CONDA_DEPS_DIR}/lib")
        endif()
endif()

# Install the lib and the tools
install(TARGETS polyview RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
