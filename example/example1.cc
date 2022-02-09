#include <OpenMesh/Core/IO/MeshIO.hh>

#include "MapsMesh.h"

int main(void) {
    Maps::MapMesh mesh;
    std::vector<Maps::MapMesh::VertexHandle> vertices(8);
    vertices[0] = mesh.add_vertex(Maps::MapMesh::Point(-1, -1, 1));
    vertices[1] = mesh.add_vertex(Maps::MapMesh::Point(1, -1, 1));
    vertices[2] = mesh.add_vertex(Maps::MapMesh::Point(1, 1, 1));
    vertices[3] = mesh.add_vertex(Maps::MapMesh::Point(-1, 1, 1));
    vertices[4] = mesh.add_vertex(Maps::MapMesh::Point(-1, -1, -1));
    vertices[5] = mesh.add_vertex(Maps::MapMesh::Point(1, -1, -1));
    vertices[6] = mesh.add_vertex(Maps::MapMesh::Point(1, 1, -1));
    vertices[7] = mesh.add_vertex(Maps::MapMesh::Point(-1, 1, -1));

    mesh.add_face(
        std::vector<Maps::MapMesh::VertexHandle>({vertices[3], vertices[0], vertices[1]}));
    mesh.add_face(
        std::vector<Maps::MapMesh::VertexHandle>({vertices[3], vertices[1], vertices[2]}));
    mesh.add_face(
        std::vector<Maps::MapMesh::VertexHandle>({vertices[2], vertices[1], vertices[6]}));
    mesh.add_face(
        std::vector<Maps::MapMesh::VertexHandle>({vertices[6], vertices[1], vertices[5]}));
    mesh.add_face(
        std::vector<Maps::MapMesh::VertexHandle>({vertices[6], vertices[3], vertices[2]}));
    mesh.add_face(
        std::vector<Maps::MapMesh::VertexHandle>({vertices[7], vertices[3], vertices[6]}));
    mesh.add_face(
        std::vector<Maps::MapMesh::VertexHandle>({vertices[1], vertices[0], vertices[4]}));
    mesh.add_face(
        std::vector<Maps::MapMesh::VertexHandle>({vertices[1], vertices[4], vertices[5]}));
    mesh.add_face(
        std::vector<Maps::MapMesh::VertexHandle>({vertices[7], vertices[0], vertices[3]}));
    mesh.add_face(
        std::vector<Maps::MapMesh::VertexHandle>({vertices[0], vertices[7], vertices[4]}));
    mesh.add_face(
        std::vector<Maps::MapMesh::VertexHandle>({vertices[4], vertices[7], vertices[5]}));
    mesh.add_face(
        std::vector<Maps::MapMesh::VertexHandle>({vertices[7], vertices[6], vertices[5]}));

    try {
        if (!OpenMesh::IO::write_mesh(mesh, "cube.off")) {
            std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
        }
    } catch (std::exception& x) {
        std::cerr << x.what() << std::endl;
    }

    mesh.Initialize();
    mesh.DownSampling();
    mesh.FaceSubDivision();
    mesh.FaceSubDivision();
    mesh.FaceSubDivision();
    mesh.FaceSubDivision();
    mesh.FaceSubDivision();

    mesh.Remesh();

    mesh.garbage_collection();
    try {
        if (!OpenMesh::IO::write_mesh(mesh, "cube_output.off")) {
            std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
        }
    } catch (std::exception& x) {
        std::cerr << x.what() << std::endl;
    }
}
