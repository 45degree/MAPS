#include <gtest/gtest.h>

#include "BaseMesh.h"

class BaseMeshTest : public testing::Test {
   protected:
    Maps::BaseMesh<> mesh;

   protected:
    void SetUp() override {
        std::vector<Maps::BaseMesh<>::VertexHandle> vertices(8);
        vertices[0] = mesh.add_vertex(Maps::BaseMesh<>::Point(-1, -1, 1));
        vertices[1] = mesh.add_vertex(Maps::BaseMesh<>::Point(1, -1, 1));
        vertices[2] = mesh.add_vertex(Maps::BaseMesh<>::Point(1, 1, 1));
        vertices[3] = mesh.add_vertex(Maps::BaseMesh<>::Point(-1, 1, 1));
        vertices[4] = mesh.add_vertex(Maps::BaseMesh<>::Point(-1, -1, -1));
        vertices[5] = mesh.add_vertex(Maps::BaseMesh<>::Point(1, -1, -1));
        vertices[6] = mesh.add_vertex(Maps::BaseMesh<>::Point(1, 1, -1));
        vertices[7] = mesh.add_vertex(Maps::BaseMesh<>::Point(-1, 1, -1));

        mesh.add_face(
            std::vector<Maps::BaseMesh<>::VertexHandle>({vertices[3], vertices[0], vertices[1]}));
        mesh.add_face(
            std::vector<Maps::BaseMesh<>::VertexHandle>({vertices[3], vertices[1], vertices[2]}));
        mesh.add_face(
            std::vector<Maps::BaseMesh<>::VertexHandle>({vertices[2], vertices[1], vertices[6]}));
        mesh.add_face(
            std::vector<Maps::BaseMesh<>::VertexHandle>({vertices[6], vertices[1], vertices[5]}));
        mesh.add_face(
            std::vector<Maps::BaseMesh<>::VertexHandle>({vertices[6], vertices[3], vertices[2]}));
        mesh.add_face(
            std::vector<Maps::BaseMesh<>::VertexHandle>({vertices[7], vertices[3], vertices[6]}));
        mesh.add_face(
            std::vector<Maps::BaseMesh<>::VertexHandle>({vertices[1], vertices[0], vertices[4]}));
        mesh.add_face(
            std::vector<Maps::BaseMesh<>::VertexHandle>({vertices[1], vertices[4], vertices[5]}));
        mesh.add_face(
            std::vector<Maps::BaseMesh<>::VertexHandle>({vertices[7], vertices[0], vertices[3]}));
        mesh.add_face(
            std::vector<Maps::BaseMesh<>::VertexHandle>({vertices[0], vertices[7], vertices[4]}));
        mesh.add_face(
            std::vector<Maps::BaseMesh<>::VertexHandle>({vertices[4], vertices[7], vertices[5]}));
        mesh.add_face(
            std::vector<Maps::BaseMesh<>::VertexHandle>({vertices[7], vertices[6], vertices[5]}));
    }
};

TEST_F(BaseMeshTest, CalculateArea) {  // NOLINT
    auto faceHandle = mesh.face_handle(0);
    double area = mesh.CalculateArea(faceHandle);
    ASSERT_FLOAT_EQ(area, 2.0);
}

TEST_F(BaseMeshTest, CalculateAngle) {  // NOLINT
    auto faceHandle = mesh.face_handle(0);
    auto vertexHandle = mesh.vertex_handle(1);
    double angle = mesh.CalculateAngle(vertexHandle, faceHandle);
    ASSERT_FLOAT_EQ(angle, M_PI / 4.0);
}

TEST_F(BaseMeshTest, TryToAddFace) {  // NOLINT
    mesh.request_vertex_status();
    mesh.request_face_status();
    mesh.request_edge_status();

    auto faceHandle = mesh.face_handle(0);
    mesh.delete_face(faceHandle, false);

    bool result = mesh.TryToAddFace(
        {mesh.vertex_handle(3), mesh.vertex_handle(1), mesh.vertex_handle(0)}, faceHandle);

    ASSERT_TRUE(result);
}

TEST_F(BaseMeshTest, TryToAddFaces) {  // NOLINT
    mesh.request_vertex_status();
    mesh.request_face_status();
    mesh.request_edge_status();

    auto faceHandle1 = mesh.face_handle(1);
    auto faceHandle2 = mesh.face_handle(2);
    auto faceHandle3 = mesh.face_handle(3);

    mesh.delete_face(faceHandle1);
    mesh.delete_face(faceHandle2);
    mesh.delete_face(faceHandle3);

    std::vector<std::array<Maps::BaseMesh<>::VertexHandle, 3>> faces;
    faces.push_back({mesh.vertex_handle(3), mesh.vertex_handle(1), mesh.vertex_handle(2)});
    faces.push_back({mesh.vertex_handle(2), mesh.vertex_handle(1), mesh.vertex_handle(6)});
    faces.push_back({mesh.vertex_handle(6), mesh.vertex_handle(1), mesh.vertex_handle(5)});

    std::vector<Maps::BaseMesh<>::FaceHandle> newFaces;
    bool result = mesh.TryToAddFaces(std::move(faces), newFaces);
    ASSERT_TRUE(result);
}

TEST_F(BaseMeshTest, IsInTriangle2D) {  // NOLINT
    using Point2D = Maps::BaseMesh<>::Point2D;
    std::array<Point2D, 3> triangle;
    triangle[0] = Point2D(0, 1);
    triangle[1] = Point2D(-1, -1);
    triangle[2] = Point2D(1, -1);
    ASSERT_TRUE(Maps::BaseMesh<>::IsInTriangle(Point2D(0, 0), triangle));
    ASSERT_FALSE(Maps::BaseMesh<>::IsInTriangle(Point2D(1, 0), triangle));
    ASSERT_FALSE(Maps::BaseMesh<>::IsInTriangle(Point2D(-1, 0), triangle));
}

TEST_F(BaseMeshTest, IsInTriangle3D) {  // NOLINT
    using Point = Maps::BaseMesh<>::Point;
    Point point1(1, 0, 0);
    Point point2(0, 1, 0);
    Point point3(0, 0, 1);
    Point point0(1.0 / 3, 1.0 / 3, 1.0 / 3);
    ASSERT_TRUE(Maps::BaseMesh<>::IsInTriangle(point0, {point1, point2, point3}));

    point1 = Point(-1.0 / 3, 1.0 / 3, 1.0 / 3);
    point2 = Point(1, 1, -1);
    point3 = Point(1, 1, 1);
    point0 = Point(-1, 0, 0);
    ASSERT_FALSE(Maps::BaseMesh<>::IsInTriangle(point0, {point1, point2, point3}));

    point1 = Point(1, 1, -1);
    point2 = Point(-1, 1, 1);
    point3 = Point(1, 1, 1);
    point0 = Point(0, 1, 0);
    ASSERT_TRUE(Maps::BaseMesh<>::IsInTriangle(point0, {point1, point2, point3}));
}

TEST_F(BaseMeshTest, IsCoplanar) {  // NOLINT
    using Point = Maps::BaseMesh<>::Point;
    Point point1(1, 0, 0);
    Point point2(0, 1, 0);
    Point point3(0, 0, 1);
    Point point0(1.0 / 3, 1.0 / 3, 1.0 / 3);

    ASSERT_TRUE(Maps::BaseMesh<>::IsCoplanar(point0, {point1, point2, point3}));

    point1 = Point(-1, 1, 1);
    point2 = Point(-1, -1, 1);
    point3 = Point(1 / 3.0, -1 / 3.0, 1 / 3.0);
    point0 = Point(-1, -1 / 2.0, 1 / 2.0);

    ASSERT_FALSE(Maps::BaseMesh<>::IsCoplanar(point0, {point1, point2, point3}));
}

TEST_F(BaseMeshTest, CalculateBaryCoor2D) {  // NOLINT
    using Point2D = Maps::BaseMesh<>::Point2D;

    Point2D point0(0, 1.0 / 2);
    Point2D point1(-1, 0);
    Point2D point2(1, 0);
    Point2D point3(0, 1);

    auto [alpha, beta, gamma] =
        Maps::BaseMesh<>::CalculateBaryCoor(point0, {point1, point2, point3});
    ASSERT_FLOAT_EQ(alpha, 1.0 / 4);
    ASSERT_FLOAT_EQ(beta, 1.0 / 4);
    ASSERT_FLOAT_EQ(gamma, 1.0 / 2);
}

TEST_F(BaseMeshTest, CalculateBaryCoor3D) {  // NOLINT
    using Point = Maps::BaseMesh<>::Point;

    Point point1(1, 0, 0);
    Point point2(0, 1, 0);
    Point point3(0, 0, 1);
    Point point0(1.0 / 3, 1.0 / 3, 1.0 / 3);
    auto [alpha, beta, gamma] =
        Maps::BaseMesh<>::CalculateBaryCoor(point0, {point1, point2, point3});
    ASSERT_FLOAT_EQ(alpha, 1.0 / 3);
    ASSERT_FLOAT_EQ(beta, 1.0 / 3);
    ASSERT_FLOAT_EQ(gamma, 1.0 / 3);
}

TEST_F(BaseMeshTest, FindFace) {  // NOLINT
    auto vertex1 = mesh.vertex_handle(0);
    auto vertex2 = mesh.vertex_handle(3);
    auto vertex3 = mesh.vertex_handle(1);
    auto face = mesh.FindFace({vertex1, vertex2, vertex3});
    ASSERT_TRUE(face.has_value());
    ASSERT_EQ(face.value(), mesh.face_handle(0));
}

TEST_F(BaseMeshTest, FindFaceForPoint) {  // NOLINT
    auto point0 = mesh.point(mesh.vertex_handle(3));
    auto point1 = mesh.point(mesh.vertex_handle(0));
    auto point2 = mesh.point(mesh.vertex_handle(1));

    auto result = mesh.FindFace(point0 * 1.0 / 3 + point1 * 1.0 / 3 + point2 * 1.0 / 3);
    ASSERT_TRUE(result.has_value());
    ASSERT_EQ(result.value(), mesh.face_handle(0));
}

TEST_F(BaseMeshTest, CalculateBaryCoor) {  // NOLINT
    auto point1 = mesh.point(mesh.vertex_handle(3));
    auto point2 = mesh.point(mesh.vertex_handle(0));
    auto point3 = mesh.point(mesh.vertex_handle(1));
    auto point0 = point1 * 1.0 / 3 + point2 * 1.0 / 3 + point3 * 1.0 / 3;

    auto result = mesh.CalculateBaryCoor(point0, mesh.face_handle(0));
    ASSERT_FLOAT_EQ(result[0].second, 1.0 / 3);
}
