#include <pmp/surface_mesh.h>
#include <pmp/io/io.h>
#include <pmp/algorithms/remeshing.h>
#include "pmp/algorithms/differential_geometry.h"
#include "pmp/algorithms/curvature.h"
#include <iostream>
#include <fstream> // For file output
#include <string>
#include <cmath>     // For sqrt
#include <algorithm> // For std::clamp
#include <memory>
#include <limits>
#include <numbers>
#include <filesystem> // For handling file paths (C++17 and later)

int main()
{
    // Get the path of the executable
    std::filesystem::path exe_path = std::filesystem::current_path();
    std::cout << "Executable directory: " << exe_path << std::endl;

    // Define input and output file names (located in the executable's directory)
    std::string input_file = (exe_path / "input.obj").string();
    std::string output_file = (exe_path / "output.obj").string();
    std::string vals_file = (exe_path / "vals.txt").string();
    // Parameters for remeshing
    double min_length = 5.0;     // Replace with appropriate values
    double max_length = 15.0;     // Replace with appropriate values
    double approx_error = 0.5; // Replace with appropriate values

    // Check if the input file exists
    if (!std::filesystem::exists(input_file))
    {
        std::cerr << "Error: Input file not found: " << input_file << std::endl;
        return 1;
    }


    // Check if the vals file exists
    if (!std::filesystem::exists(vals_file))
    {
        std::cerr << "Error: vals.txt file not found: " << vals_file
                  << std::endl;
        return 1;
    }

    // Read the values from vals.txt into a std::vector<double>
    std::vector<double> vertex_length_values;
    std::ifstream vals_stream(vals_file);
    if (!vals_stream)
    {
        std::cerr << "Error: Could not open vals.txt for reading.\n";
        return 1;
    }

    double value;
    while (vals_stream >> value) // Read doubles line by line
    {
        vertex_length_values.push_back(value);
    }
    vals_stream.close();

     // Print out the values to verify
    std::cout << "Read " << vertex_length_values.size()
              << " values from vals.txt:\n";
    for (const auto& v : vertex_length_values)
    {
        std::cout << v << "\n";
    }
    // Load the mesh
    pmp::SurfaceMesh mesh;
    try
    {
        pmp::read(mesh, input_file);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error loading mesh: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "Loaded mesh with " << mesh.n_vertices() << " vertices and "
              << mesh.n_faces() << " faces.\n";

   // assert(target_edge_lengths.size() == mesh_.n_vertices() &&
          // "Mismatch in size of target_edge_lengths and mesh vertices.");


    std::cout << "Remeshing...\n";

    // Perform remeshing
    try
    {
        pmp::custom_remeshing(mesh, vertex_length_values);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error during remeshing: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "Remeshing complete. Resulting mesh has " << mesh.n_vertices()
              << " vertices and " << mesh.n_faces() << " faces.\n";

    // Save the remeshed output
    try
    {
        pmp::write(mesh, output_file);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error writing output mesh: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "Output written to " << output_file << "\n";

    return 0;
}
