cmake_minimum_required(VERSION 3.6)
project(object_tracking_sample)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_subdirectory(object_tracker)
add_subdirectory(object_tracker_binary)
