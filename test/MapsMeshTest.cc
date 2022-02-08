#include <gtest/gtest.h>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <cstdio>

#include "MapsMesh.h"

class MapsMeshTest : public ::testing::Test {
   protected:
    using Coordinate2Dpair = Maps::MapMesh::Coordinate2DPair;
    using VertexHandle = Maps::MapMesh::VertexHandle;

   protected:
    void SetUp() override {
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
    }

   protected:
    Maps::MapMesh mesh;

   protected:
    bool isEqual(const std::array<VertexHandle, 3>& face1,
                 const std::array<VertexHandle, 3>& face2) {
        if (face1.size() != face2.size()) return false;
        for (const auto& face : face1) {
            if (std::find(face2.begin(), face2.end(), face) == face2.end()) return false;
        }
        return true;
    }
};

TEST_F(MapsMeshTest, MVTTrangle) {  // NOLINT
    Coordinate2Dpair coordinates;
    coordinates.push_back({mesh.vertex_handle(0), {-0.5, -0.5}});
    coordinates.push_back({mesh.vertex_handle(1), {1.5, -0.5}});
    coordinates.push_back({mesh.vertex_handle(2), {1.5, 1.5}});
    coordinates.push_back({mesh.vertex_handle(3), {-0.5, 1.5}});

    std::vector<std::array<VertexHandle, 3>> faces;
    std::array<std::pair<VertexHandle, double>, 3> barycentricCoordinates;
    Maps::MapMesh::MVTTrangle(coordinates, 3, faces, barycentricCoordinates);

    ASSERT_EQ(2, faces.size());
    std::array<VertexHandle, 3> face1{mesh.vertex_handle(0), mesh.vertex_handle(1),
                                      mesh.vertex_handle(3)};
    std::array<VertexHandle, 3> face2{mesh.vertex_handle(3), mesh.vertex_handle(1),
                                      mesh.vertex_handle(2)};
    ASSERT_TRUE(isEqual(faces[0], face1) || isEqual(faces[0], face2));
    ASSERT_TRUE(isEqual(faces[1], face1) || isEqual(faces[1], face2));

    for (const auto& barycentricCoordinate : barycentricCoordinates) {
        if (barycentricCoordinate.first == mesh.vertex_handle(0)) {
            ASSERT_FLOAT_EQ(barycentricCoordinate.second, 0.5);
        } else if (barycentricCoordinate.first == mesh.vertex_handle(1)) {
            ASSERT_FLOAT_EQ(barycentricCoordinate.second, 0.25);
        } else if (barycentricCoordinate.first == mesh.vertex_handle(3)) {
            ASSERT_FLOAT_EQ(barycentricCoordinate.second, 0.25);
        } else {
            ASSERT_TRUE(false);
        }
    }
}

TEST_F(MapsMeshTest, CDTTrangle) {  // NOLINT
    Coordinate2Dpair coordinates;
    coordinates.push_back({mesh.vertex_handle(0), {-0.5, -0.5}});
    coordinates.push_back({mesh.vertex_handle(1), {1.5, -0.5}});
    coordinates.push_back({mesh.vertex_handle(2), {-0.5, 1.5}});

    std::vector<std::array<VertexHandle, 3>> faces;
    std::array<std::pair<VertexHandle, double>, 3> barycentricCoordinates;
    ASSERT_TRUE(Maps::MapMesh::CDTTrangle(coordinates, faces, barycentricCoordinates) != -1);

    for (const auto& barycentricCoordinate : barycentricCoordinates) {
        if (barycentricCoordinate.first == mesh.vertex_handle(0)) {
            ASSERT_FLOAT_EQ(barycentricCoordinate.second, 0.5);
        } else if (barycentricCoordinate.first == mesh.vertex_handle(1)) {
            ASSERT_FLOAT_EQ(barycentricCoordinate.second, 0.25);
        } else if (barycentricCoordinate.first == mesh.vertex_handle(2)) {
            ASSERT_FLOAT_EQ(barycentricCoordinate.second, 0.25);
        } else {
            ASSERT_TRUE(false);
        }
    }
}

TEST_F(MapsMeshTest, FaceSubDivision) {  // NOLINT
    mesh.request_face_status();
    mesh.request_edge_status();
    mesh.request_vertex_status();

    size_t faceCount = mesh.n_faces();
    ASSERT_EQ(12, faceCount);
    mesh.FaceSubDivision();
    mesh.garbage_collection();
    ASSERT_EQ(faceCount * 4, mesh.n_faces());
}

TEST_F(MapsMeshTest, ReCalculate2DCoordinates) {  // NOLINT
    using BaryCoor = Maps::BaryCoor;
    using VertexHandle = Maps::MapMesh::VertexHandle;
    using Point2D = Maps::MapMesh::Point2D;

    auto vertexHandle = mesh.add_vertex(Maps::MapMesh::Point(0, 0, 0));
    auto faceHandle = mesh.face_handle(0);
    auto ringVertex = std::vector<VertexHandle>(mesh.fv_begin(faceHandle), mesh.fv_end(faceHandle));
    ASSERT_EQ(ringVertex.size(), 3);

    BaryCoor baryCoor;
    baryCoor[0] = {ringVertex[0], 0.5};
    baryCoor[1] = {ringVertex[1], 0.25};
    baryCoor[2] = {ringVertex[2], 0.25};
    mesh.data(vertexHandle).baryCoor = baryCoor;
    mesh.data(faceHandle).vertrices.push_back(vertexHandle);

    std::map<Maps::MapMesh::VertexHandle, Maps::MapMesh::Point2D> originPoint2D;
    originPoint2D[ringVertex[0]] = Point2D(0, 0);
    originPoint2D[ringVertex[1]] = Point2D(0, 1);
    originPoint2D[ringVertex[2]] = Point2D(1, 0);

    auto resultOption = mesh.ReCalculate2DCoordinates(originPoint2D, vertexHandle);
    ASSERT_TRUE(resultOption.has_value());
    auto result = resultOption.value();
    ASSERT_FLOAT_EQ(1.0 / 4, result[0]);
    ASSERT_FLOAT_EQ(1.0 / 4, result[1]);
}

TEST_F(MapsMeshTest, TEST1) {  // NOLINT
    mesh.Initialize();
    mesh.DownSampling();

    for (auto vertexIter = mesh.vertices_begin(); vertexIter != mesh.vertices_end(); vertexIter++) {
        auto vertex = *vertexIter;
        const auto& b = mesh.data(vertex).baryCoor;
        if (b.has_value()) {
            const auto& value = b.value();
            printf("%d, %d: %lf, %d: %lf, %d: %lf\n", vertex.idx(), value[0].first.idx(),
                   value[0].second, value[1].first.idx(), value[1].second, value[2].first.idx(),
                   value[2].second);
        } else {
            printf("%d, none\n", vertex.idx());
        }
    }

    std::cout << "new face" << std::endl;
    for (const auto& face : mesh.faces()) {
        std::cout << face.idx() << ':';
        for (const auto& vertexHandle : mesh.fv_range(face)) {
            std::cout << vertexHandle.idx() << ',';
        }
        std::cout << std::endl;
    }

    mesh.Initialize();
    mesh.DownSampling();

    for (auto vertexIter = mesh.vertices_begin(); vertexIter != mesh.vertices_end(); vertexIter++) {
        auto vertex = *vertexIter;
        const auto& b = mesh.data(vertex).baryCoor;
        if (b.has_value()) {
            const auto& value = b.value();
            printf("%d, %d: %lf, %d: %lf, %d: %lf\n", vertex.idx(), value[0].first.idx(),
                   value[0].second, value[1].first.idx(), value[1].second, value[2].first.idx(),
                   value[2].second);
        } else {
            printf("%d, none\n", vertex.idx());
        }
    }
}

TEST_F(MapsMeshTest, TEST2) {  // NOLINT
    mesh.Initialize();
    mesh.DownSampling();
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
