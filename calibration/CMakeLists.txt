add_subdirectory(module/pedestal)
add_subdirectory(module/MIP)
add_subdirectory(module/intercalib)
add_subdirectory(module/spe)

# Collect all calibration algorithms into a single interface library
add_library(calib_algs INTERFACE)
target_link_libraries(calib_algs INTERFACE
  PedestalAlg MIPAlg InterCalibAlg SPEAlg
)