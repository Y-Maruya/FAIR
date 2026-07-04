add_subdirectory(PedestalReader)
add_subdirectory(TriggersReader)
add_subdirectory(InterCalibReader)
add_subdirectory(MIPReader)

# Collect all calibration database I/O into a single interface library
add_library(calib_io INTERFACE)
target_link_libraries(calib_io INTERFACE
  PedestalReader
  TriggersReader
  InterCalibReader
  MIPReader
)