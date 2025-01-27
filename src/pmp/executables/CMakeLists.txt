# Create a custom project for executables
project(executables)

# Specify where executables will be built
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/executables)

# Add executables
add_executable(c_remesh c_remesh.cpp)

# Link required libraries for each executable
# (Example assumes you need to link to pmp-library)
target_link_libraries(c_remesh pmp)

set_target_properties(c_remesh PROPERTIES
    CXX_STANDARD 20
    RUNTIME_OUTPUT_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}
    FOLDER "Executables"
)
