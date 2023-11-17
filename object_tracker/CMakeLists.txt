add_library(object-tracker
  src/object_tracker.cc
)

target_include_directories(object-tracker
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/include
)
