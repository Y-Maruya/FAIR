add_library(MIPSimAlg STATIC MIPSimAlg.cpp)
target_include_directories(MIPSimAlg PUBLIC ${FAIR_ROOT_DIR})
target_link_libraries(MIPSimAlg PUBLIC fair_options MIPAlg)
add_library(FAIR::MIPSimAlg ALIAS MIPSimAlg)
