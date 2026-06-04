add_executable(test_expand_config
  test_expand_config.cpp
)

target_link_libraries(test_expand_config
  PRIVATE
    fair_options
)

target_include_directories(test_expand_config
  PRIVATE
    ${CMAKE_SOURCE_DIR}
)

add_test(NAME test_expand_config COMMAND test_expand_config)
