add_executable(object-tracker-binary
    src/main.cc
)

target_link_libraries(object-tracker-binary
  object-tracker
)