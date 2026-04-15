# reco_alg/CMakeLists.txt
add_subdirectory(reader)
#add_subdirectory(writer)

# Collect all IO readers into a single interface library
add_library(io_readers INTERFACE)
target_link_libraries(io_readers INTERFACE
  RootRawHitReader BinaryRawHitReader SimHitReader
)