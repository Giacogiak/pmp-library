#include <pmp/surface_mesh.h>
#include <pmp/io/io.h>
#include <pmp/algorithms/remeshing.h>
#include "pmp/algorithms/differential_geometry.h"
#include "pmp/algorithms/curvature.h"
#include <iostream>
#include <fstream> // For file output
#include <string>
#include <cmath> // For sqrt
#include <algorithm> // For std::clamp
#include <cmath>
#include <memory>
#include <limits>
#include <numbers>
#include <random>

int main(int argc, char** argv)
{
    // Check command-line arguments
    if (argc != 7)
    {
        std::cerr << "Usage: " << argv[0] << " <input.obj> <edge_length> <output_curvatures.txt> <output_vsizing.txt> <output.obj> <mode>\n";
        return 1;
    }

    std::string input_file = argv[1];
    double edge_length = std::stod(argv[2]);
    std::string curvatures_file = argv[3];
    std::string vsizing_file = argv[4];
    std::string output_file = argv[5];
    int mode = std::stoi(argv[6]);
    // Load the mesh
    pmp::SurfaceMesh mesh;
    pmp::read(mesh, input_file);

    std::cout << "Loaded mesh with " << mesh.n_vertices() << " vertices and "
              << mesh.n_faces() << " faces.\n";

    // Properties
    auto vfeature_ = mesh.vertex_property<bool>("v:feature", false);
    auto efeature_ = mesh.edge_property<bool>("e:feature", false);
    auto vlocked_ = mesh.add_vertex_property<bool>("v:locked", false);
    auto elocked_ = mesh.add_edge_property<bool>("e:locked", false);
    auto vsizing_ = mesh.add_vertex_property<pmp::Scalar>("v:sizing");

    // Compute curvatures
    pmp::curvature(mesh, pmp::Curvature::max_abs, 0, true, false);
    auto curvatures = mesh.get_vertex_property<pmp::Scalar>("v:curv");

    // Smooth curvatures while considering features
    for (auto v : mesh.vertices())
    {
        if (vfeature_[v])
            continue;
        pmp::Scalar curv = 0, weight = 0, sum_weights = 0;
        for (auto vh : mesh.halfedges(v))
        {
            auto vv = mesh.to_vertex(vh);
            if (vfeature_[vv])
                continue;
            weight = std::max(0.0, pmp::cotan_weight(mesh, mesh.edge(vh)));
            sum_weights += weight;
            curv += weight * curvatures[vv];
        }
        if (sum_weights > 0)
            curvatures[v] = curv / sum_weights;
    }



     // Write curvature_ values to a text file
    std::ofstream offstream(curvatures_file);
    if (!offstream)
    {
        std::cerr << "Error: Could not open file " << curvatures_file << " for writing\n";
        return 1;
    }

    for (auto v : mesh.vertices())
    {
        offstream << curvatures[v] << '\n';
    }
    offstream.close();

    std::cout << "Curvatures values written to " << curvatures_file << "\n";



    // Convert curvatures to edge lengths
    pmp::Scalar min_edge_length_ = edge_length * 0.5;
    pmp::Scalar max_edge_length_ = edge_length * 2.0;
    pmp::Scalar approx_error_ = edge_length * 0.1;
    pmp::Scalar edge_length_ = edge_length;
    for (auto v : mesh.vertices())
    {
        pmp::Scalar c = curvatures[v];

        // Compute edge length from curvature
        const pmp::Scalar r = 1.0 / std::max(c, 1e-6f);
        pmp::Scalar h = sqrt(6.0 * approx_error_ * r - 3.0 * approx_error_ * approx_error_);

        // Clamp to min and max edge lengths
        h = std::clamp(h, min_edge_length_, max_edge_length_);
        vsizing_[v] = h;
    }

    // Write vsizing_ values to a text file
    std::ofstream ofs(vsizing_file);
    if (!ofs)
    {
        std::cerr << "Error: Could not open file " << vsizing_file << " for writing\n";
        return 1;
    }

    for (auto v : mesh.vertices())
    {
        ofs << vsizing_[v] << '\n';
    }
    ofs.close();

    std::cout << "Vertex sizing values written to " << vsizing_file << "\n";

    // Perform remeshing
    if (mode == 1)
    {
        std::cout << "Remeshing Uniform..\n";
        pmp::uniform_remeshing(mesh, edge_length_);
    }
       
    else if (mode == 2)
    {
        std::cout << "Remeshing Adaptive..\n";
        pmp::adaptive_remeshing(mesh, min_edge_length_, max_edge_length_, approx_error_);
    }
        
    else if (mode == 3)
    {
        std::cout << "Remeshing Custom..\n";
        
        // Generate random scalar values for each vertex
        std::vector<double> vertex_scalars(mesh.n_vertices());
        std::random_device rd;  // Seed for random number generation
        std::mt19937 gen(rd()); // Mersenne Twister random number generator
        std::uniform_real_distribution<> dist(min_edge_length_, min_edge_length_);

        for (auto& scalar : vertex_scalars)
        {
            scalar = std::round(dist(gen));
        }
       
        // Perform remeshing using the custom function
        pmp::custom_remeshing(mesh, vertex_scalars);

    }

    else
    {
        std::cerr << "Error: Invalid mode value: " << mode << '\n';
        return 1;
    }
        
    std::cout << "Remeshing complete. Resulting mesh has " << mesh.n_vertices() << " vertices and "
              << mesh.n_faces() << " faces.\n";

    // Save the remeshed output
    pmp::write(mesh, output_file);

    std::cout << "Output written to " << output_file << "\n";

    return 0;
}
