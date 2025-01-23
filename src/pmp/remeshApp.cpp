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
#include <memory>
#include <limits>
#include <numbers>

int main(int argc, char** argv)
{
    // Check command-line arguments
    if (argc != 6)
    {
        std::cerr << "Usage: " << argv[0] << " <input.obj> <edge_length> <output_vsizing.txt> <output.obj> <mode>\n";
        return 1;
    }

    std::string input_file = argv[1];
    double edge_length = std::stod(argv[2]);
    std::string vsizing_file = argv[3];
    std::string output_file = argv[4];
    int mode = std::stoi(argv[5]);
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

     // curvature values for feature vertices and boundary vertices
    // are not meaningful. mark them as negative values.
    for (auto v : mesh.vertices())
    {
        if (mesh.is_boundary(v) || (vfeature_ && vfeature_[v]))
            vsizing_[v] = -1.0;
        else
            vsizing_[v] = curvatures[v];
    }
    // curvature values might be noisy. smooth them.
    // don't consider feature vertices' curvatures.
    // don't consider boundary vertices' curvatures.
    // do this for two iterations, to propagate curvatures
    // from non-feature regions to feature vertices.
     for (int iters = 0; iters < 2; ++iters)
    {
        for (auto v : mesh.vertices())
        {
            pmp::Scalar w, ww = 0.0;
            pmp::Scalar c, cc = 0.0;

            for (auto h : mesh.halfedges(v))
            {
                c = vsizing_[mesh.to_vertex(h)];
                if (c > 0.0)
                {
                    w = std::max(0.0, pmp::cotan_weight(mesh, mesh.edge(h)));
                    ww += w;
                    cc += w * c;
                }
            }

            if (ww)
                cc /= ww;
            vsizing_[v] = cc;
        }
    }

    // Convert curvatures to edge lengths
    pmp::Scalar min_edge_length_ = edge_length * 0.5;
    pmp::Scalar max_edge_length_ = edge_length * 2.0;
    pmp::Scalar approx_error_ = edge_length * 0.1;
    pmp::Scalar edge_length_ = edge_length;
    
      // now convert per-vertex curvature into target edge length
    for (auto v : mesh.vertices())
    {
        pmp::Scalar c = vsizing_[v];

        // get edge length from curvature
        const pmp::Scalar r = 1.0 / c;
        const pmp::Scalar e = approx_error_;
        pmp::Scalar h;
        if (e < r)
        {
            // see mathworld: "circle segment" and "equilateral triangle"
            //h = sqrt(2.0*r*e-e*e) * 3.0 / sqrt(3.0);
            h = sqrt(6.0 * e * r - 3.0 * e * e); // simplified...
        }
        else
        {
            // this does not really make sense
            h = e * 3.0 / sqrt(3.0);
        }

        // clamp to min. and max. edge length
        if (h < min_edge_length_)
            h = min_edge_length_;
        else if (h > max_edge_length_)
            h = max_edge_length_;

        // store target edge length
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

    std::vector<double> vertex_scalars(mesh.n_vertices(), 0.0);
    int cnt = 0;
    for (auto v : mesh.vertices())
    {
        vertex_scalars[cnt] = vsizing_[v];
        cnt++;
    }

    std::cout << "Remeshing Custom..\n";
       
    // Perform remeshing using the custom function
    pmp::adaptive_remeshing(mesh, min_edge_length_, max_edge_length_, approx_error_);

    

 
        
    std::cout << "Remeshing complete. Resulting mesh has " << mesh.n_vertices() << " vertices and "
              << mesh.n_faces() << " faces.\n";

    // Save the remeshed output
    pmp::write(mesh, output_file);

    std::cout << "Output written to " << output_file << "\n";

    return 0;
}
