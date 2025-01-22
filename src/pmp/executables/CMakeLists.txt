# Add the executable
add_executable(remeshApp remeshApp.cpp)
set_property(TARGET remeshApp PROPERTY CXX_STANDARD 20)
target_link_libraries(basics pmp)

