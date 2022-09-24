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
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS})

find_package(OpenGL REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11")

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
target_link_libraries(polyview polyview_lib Qt5::Widgets ${OPENGL_opengl_LIBRARY})

