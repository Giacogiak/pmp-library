# Add the executable
add_executable(c_remesh c_remesh.cpp)
set_property(TARGET c_remesh PROPERTY CXX_STANDARD 20)
target_link_libraries(c_remesh pmp)