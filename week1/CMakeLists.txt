cmake_minimum_required(VERSION 3.15)
project(week1)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES glad.c main.cpp camera.h densityMap.h shader.h shader.cpp densityMap.cpp camera.cpp data.h)
add_compile_options(-O3)
add_executable(week1 ${SOURCE_FILES})
target_link_libraries(week1 GL GLU glfw3 X11 Xxf86vm Xrandr pthread Xi dl Xinerama Xcursor assimp)