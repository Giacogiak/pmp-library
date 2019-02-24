//=============================================================================
// Copyright (C) 2011-2019 The pmp-library developers
//
// This file is part of the Polygon Mesh Processing Library.
// Distributed under the terms of the MIT license, see LICENSE.txt for details.
//
// SPDX-License-Identifier: MIT
//=============================================================================

#include "MeshViewer.h"

#include <imgui.h>

#include <cfloat>
#include <iostream>
#include <sstream>

//=============================================================================

namespace pmp {

//=============================================================================

MeshViewer::MeshViewer(const char* title, int width, int height, bool showgui)
    : TrackballViewer(title, width, height, showgui)
{
    // setup draw modes
    clear_draw_modes();
    add_draw_mode("Points");
    add_draw_mode("Hidden Line");
    add_draw_mode("Smooth Shading");
    add_draw_mode("Texture");
    set_draw_mode("Smooth Shading");

    crease_angle_ = 90.0;
}

//-----------------------------------------------------------------------------

MeshViewer::~MeshViewer() = default;

//-----------------------------------------------------------------------------

bool MeshViewer::load_mesh(const char* filename)
{
    // load mesh
    if (mesh_.read(filename))
    {
        // update scene center and bounds
        BoundingBox bb = mesh_.bounds();
        set_scene((vec3)bb.center(), 0.5 * bb.size());

        // compute face & vertex normals, update face indices
        update_mesh();

        // set draw mode
        if (mesh_.n_faces())
        {
            set_draw_mode("Solid Smooth");
        }
        else if (mesh_.n_vertices())
        {
            set_draw_mode("Points");
        }

        // print mesh statistic
        std::cout << "Load " << filename << ": " << mesh_.n_vertices()
                  << " vertices, " << mesh_.n_faces() << " faces\n";

        filename_ = filename;
        crease_angle_ = mesh_.crease_angle();
        return true;
    }

    std::cerr << "Failed to read mesh from " << filename << " !" << std::endl;
    return false;
}

//-----------------------------------------------------------------------------

bool MeshViewer::load_texture(const char* filename, GLint format,
                              GLint min_filter, GLint mag_filter, GLint wrap)
{
    // load texture from file
    if (!mesh_.load_texture(filename, format, min_filter, mag_filter, wrap))
        return false;

    set_draw_mode("Texture");

    // set material
    mesh_.set_ambient(1.0);
    mesh_.set_diffuse(0.9);
    mesh_.set_specular(0.0);
    mesh_.set_shininess(1.0);

    return true;
}

//-----------------------------------------------------------------------------

void MeshViewer::update_mesh()
{
    // re-compute face and vertex normals
    mesh_.update_opengl_buffers();
}

//-----------------------------------------------------------------------------

void MeshViewer::process_imgui()
{
    if (ImGui::CollapsingHeader("Mesh Info", ImGuiTreeNodeFlags_DefaultOpen))
    {
        // output mesh statistics
        ImGui::BulletText("%d vertices", (int)mesh_.n_vertices());
        ImGui::BulletText("%d edges", (int)mesh_.n_edges());
        ImGui::BulletText("%d faces", (int)mesh_.n_faces());

        if (ImGui::Button("Reload Model"))
        {
            load_mesh(filename_.c_str());
        }

        // control crease angle
        ImGui::PushItemWidth(100);
        ImGui::SliderFloat("Crease Angle", &crease_angle_, 0.0f, 180.0f,
                           "%.0f");
        ImGui::PopItemWidth();
        if (crease_angle_ != mesh_.crease_angle())
        {
            mesh_.set_crease_angle(crease_angle_);
        }
    }
}

//-----------------------------------------------------------------------------

void MeshViewer::draw(const std::string& drawMode)
{
    // draw mesh
    mesh_.draw(projection_matrix_, modelview_matrix_, drawMode);
}

//-----------------------------------------------------------------------------

void MeshViewer::keyboard(int key, int scancode, int action, int mods)
{
    if (action != GLFW_PRESS && action != GLFW_REPEAT)
        return;

    switch (key)
    {
        case GLFW_KEY_BACKSPACE: // reload model
        {
            load_mesh(filename_.c_str());
            break;
        }

        case GLFW_KEY_C: // adjust crease angle
        {
            if (mods & GLFW_MOD_SHIFT)
                mesh_.set_crease_angle(mesh_.crease_angle() + 10);
            else
                mesh_.set_crease_angle(mesh_.crease_angle() - 10);
            crease_angle_ = mesh_.crease_angle();

            std::cout << "crease angle: " << mesh_.crease_angle() << std::endl;
            break;
        }

        case GLFW_KEY_O: // write mesh
        {
            mesh_.write("output.off");
            break;
        }

        default:
        {
            TrackballViewer::keyboard(key, scancode, action, mods);
            break;
        }
    }
}

//-----------------------------------------------------------------------------

Vertex MeshViewer::pick_vertex(int x, int y)
{
    Vertex vmin;

    vec3 p;
    Scalar d, dmin(FLT_MAX);

    if (TrackballViewer::pick(x, y, p))
    {
        Point picked_position(p);
        for (auto v : mesh_.vertices())
        {
            d = distance(mesh_.position(v), picked_position);
            if (d < dmin)
            {
                dmin = d;
                vmin = v;
            }
        }
    }
    return vmin;
}

//=============================================================================
} // namespace pmp
//=============================================================================
