#include "MapsMesh.h"

#include <igl/principal_curvature.h>
#include <poly2tri.h>

#include <Eigen/Dense>
#include <limits>

namespace Maps {

// TODO(45degree): 建立一个函数 判断点是否在三角形内
// TODO(45degree): 建立一个函数 计算重心坐标
// TODO(45degree): 建立一个函数 根据重心坐标还原坐标

template <class Mesh>
static void OpenMesh2IGL(const Mesh *mesh, Eigen::MatrixX3d &V, Eigen::MatrixX3i &F) {

    assert(mesh->is_trimesh());

    auto verticesCount = std::distance(mesh->vertices_sbegin(), mesh->vertices_end());
    auto facesCount = std::distance(mesh->faces_sbegin(), mesh->faces_end());
    V.resize(verticesCount, 3);
    F.resize(facesCount, 3);

    std::map<typename Mesh::VertexHandle, int> newVertexHandleMap;
    int j = 0;
    for (auto i = mesh->vertices_sbegin(); i != mesh->vertices_end(); i++, j++) {
        auto p = mesh->point(*i);
        V(j, 0) = p[0];
        V(j, 1) = p[1];
        V(j, 2) = p[2];
        newVertexHandleMap[*i] = j;
    }

    j = 0;
    for (auto i = mesh->faces_sbegin(); i != mesh->faces_end(); i++, j++) {
        auto fh = *i;
        int vi = 0;
        for (typename Mesh::ConstFaceVertexCCWIter fvi = mesh->cfv_ccwbegin(fh);
             vi < 3 && fvi != mesh->cfv_ccwend(fh); ++fvi) {
            F(j, vi) = newVertexHandleMap[*fvi];
            vi++;
        }
    }
}

bool MapMesh::IsInTriangle(const Point2D &point, const std::array<Point2D, 3> &triangle2D) {

    Point point0(point[0], point[1], 0);
    Point point1(triangle2D[0][0], triangle2D[0][1], 0);
    Point point2(triangle2D[1][0], triangle2D[1][1], 0);
    Point point3(triangle2D[2][0], triangle2D[2][1], 0);

    double z1 = (point0 - point1).cross(point2 - point1)[2];
    double z2 = (point0 - point2).cross(point3 - point2)[2];
    double z3 = (point0 - point3).cross(point1 - point3)[2];

    if ((z1 >= 0 && z2 >= 0 && z3 >= 0) || (z1 <= 0 && z2 <= 0 && z3 <= 0)) return true;

    return false;
}

bool MapMesh::IsInTriangle(const Point &point, const std::array<Point, 3> &triangle) {
    if (!IsCoplanar(point, triangle)) return false;

    double z1 = (point - triangle[0]).cross(triangle[1] - triangle[0])[2];
    double z2 = (point - triangle[1]).cross(triangle[2] - triangle[1])[2];
    double z3 = (point - triangle[2]).cross(triangle[0] - triangle[2])[2];

    if ((z1 >= 0 && z2 >= 0 && z3 >= 0) || (z1 <= 0 && z2 <= 0 && z3 <= 0)) return true;

    return false;
}

bool MapMesh::IsCoplanar(const Point &point, const std::array<Point, 3> &triangle) {
    double esp = 1e-5;
    auto Vec1 = triangle[0] - point;
    auto Vec2 = triangle[1] - point;
    auto Vec3 = triangle[2] - point;

    Eigen::Matrix3d mat;
    mat(0, 0) = Vec1[0], mat(0, 1) = Vec1[1], mat(0, 2) = Vec1[2];
    mat(1, 0) = Vec2[0], mat(1, 1) = Vec2[1], mat(1, 2) = Vec2[2];
    mat(2, 0) = Vec2[0], mat(2, 1) = Vec3[1], mat(2, 2) = Vec3[2];

    double result = mat.determinant();
    return result >= -esp && result <= esp;
}

std::tuple<double, double, double> MapMesh::CalculateBaryCoor(
    const MapMesh::Point2D &point, const std::array<MapMesh::Point2D, 3> &triangle2D) {
    Point point0 = Point(point[0], point[1], 0);
    Point point1 = Point(triangle2D[0][0], triangle2D[0][1], 0);
    Point point2 = Point(triangle2D[1][0], triangle2D[1][1], 0);
    Point point3 = Point(triangle2D[2][0], triangle2D[2][1], 0);

    double area = (point2 - point1).cross(point3 - point1).norm();
    double alpha = (point0 - point2).cross(point0 - point3).norm() / area;
    double beta = (point0 - point1).cross(point0 - point3).norm() / area;
    double gamma = (point0 - point1).cross(point0 - point2).norm() / area;

    return {alpha, beta, gamma};
}

std::tuple<double, double, double> MapMesh::CalculateBaryCoor(
    const Point &point, const std::array<Point, 3> &triangle) {
    Eigen::Matrix3d A;
    // clang-format off
    A << triangle[0][0], triangle[0][1], triangle[0][2],
         triangle[1][0], triangle[1][1], triangle[1][2],
         triangle[2][0], triangle[2][1], triangle[2][2];
    // clang-format on
    Eigen::Vector3d B(point[0], point[1], point[2]);
    Eigen::Vector3d result = A.ldlt().solve(B);

    return {result[0], result[1], result[2]};
}

// TODO(45degree): 移除重心坐标计算
// TODO(45degree): 重新设计借口
int MapMesh::CDTTrangle(const Coordinate2Dpair &coordinates,
                        std::vector<std::array<VertexHandle, 3>> &faces,
                        BarycentricCoordinates &barycentricCoordinates) {

    faces.clear();

    std::vector<p2t::Point *> points(coordinates.size());
    for (int i = 0; i < coordinates.size(); i++) {
        double x = coordinates[i].second[0];
        double y = coordinates[i].second[1];

        points[i] = new p2t::Point(x, y);
    }

    p2t::CDT cdt(points);

    try {
        cdt.Triangulate();
    } catch (const std::runtime_error &exception) {
        std::cout << exception.what() << std::endl;
        throw exception;
    }

    auto triangles = cdt.GetTriangles();
    if (!p2t::IsDelaunay(triangles)) {
        std::cout << "can't delaunay triangles" << std::endl;
        return -1;
    }

    int trangleIdx = -1;
    for (int j = 0; j < triangles.size(); j++) {
        auto &triangle = triangles[j];
        std::array<Point2D, 3> triangle2D;
        triangle2D[0] = Point2D(triangle->GetPoint(0)->x, triangle->GetPoint(0)->y);
        triangle2D[1] = Point2D(triangle->GetPoint(1)->x, triangle->GetPoint(1)->y);
        triangle2D[2] = Point2D(triangle->GetPoint(2)->x, triangle->GetPoint(2)->y);

        std::array<VertexHandle, 3> face;
        for (int i = 0; i < 3; i++) {
            auto position = std::find(points.begin(), points.end(), triangle->GetPoint(i));
            if (position == points.end()) return -1;
            face[i] = coordinates[position - points.begin()].first;
        }
        faces.push_back(face);

        /**
         * 计算原点的重心坐标
         */

        if (!IsInTriangle(Point2D(0, 0), triangle2D)) continue;

        auto [alpha, beta, gamma] = CalculateBaryCoor(Point2D(0, 0), triangle2D);
        barycentricCoordinates[0] = std::make_pair(face[0], alpha);
        barycentricCoordinates[1] = std::make_pair(face[1], beta);
        barycentricCoordinates[2] = std::make_pair(face[2], gamma);
        trangleIdx = j;
    }

    for (auto &point : points) {
        delete point;
    }

    return trangleIdx;
}

int MapMesh::MVTTrangle(const Coordinate2Dpair &coordinates, int startIdx,
                        std::vector<std::array<VertexHandle, 3>> &faces,
                        BarycentricCoordinates &barycentricCoordinates) {
    faces.clear();
    int trangleIdx = -1;
    for (size_t i = 0, ii = 1; i < coordinates.size(); i++, ii++, ii %= coordinates.size()) {
        if (i == startIdx || ii == startIdx) continue;

        std::array<VertexHandle, 3> face(
            {coordinates[i].first, coordinates[ii].first, coordinates[startIdx].first});

        faces.push_back(face);

        std::array<Point2D, 3> triangle;
        triangle[0] = coordinates[startIdx].second;
        triangle[1] = coordinates[i].second;
        triangle[2] = coordinates[ii].second;

        /**
         * 计算原点的重心坐标
         */

        if (!IsInTriangle(Point2D(0, 0), triangle)) continue;

        auto [alpha, beta, gamma] = CalculateBaryCoor(Point2D(0, 0), triangle);
        barycentricCoordinates[0] = std::make_pair(coordinates[startIdx].first, alpha);
        barycentricCoordinates[1] = std::make_pair(coordinates[i].first, beta);
        barycentricCoordinates[2] = std::make_pair(coordinates[ii].first, gamma);

        trangleIdx = static_cast<int>(i);
    }
    return trangleIdx;
}

std::optional<MapMesh::Point2D> MapMesh::ReCalculate2DCoordinates(
    std::map<VertexHandle, Point2D> &originCoor, VertexHandle deleteVertex) {
    if (!IsVertexDeleted(deleteVertex)) return std::nullopt;
    auto baryCoor = data(deleteVertex).barycentricCoordinates.value();

    return originCoor[baryCoor[0].first] * baryCoor[0].second +
           originCoor[baryCoor[1].first] * baryCoor[1].second +
           originCoor[baryCoor[2].first] * baryCoor[2].second;
}

void MapMesh::ReTrangleAndAddFace(const VertexHandle &deleteVertex) {
    // TODO(45degree): 记录1领域面中所有包含的点, 并更新这些点的参数
    // 假设1领域面f上记录了一个点v, v中记录了v在f上的重心坐标(a, b, c), 面f在经过2维映射变成了f'
    // 更新过程: 1. 根据f'和(a, b, c)计算v的二维映射点v'
    //           2. 重新计算包含v'的三角形p‘, 并根据p’重新计算坐标(a', b', c')
    //           3. 删除f, 添加面p, 更新v的重心坐标, 建立p与v的关联

    std::vector<FaceHandle> ringFaces(vf_begin(deleteVertex), vf_end(deleteVertex));
    Coordinate2Dpair coordinates;
    Calculate2D(deleteVertex, coordinates);

    // 重新计算1领域面上所有包含的点的二维坐标
    std::map<VertexHandle, Point2D> currentPoint2DMap(coordinates.begin(), coordinates.end());
    std::map<VertexHandle, Point2D> deletePoint2DMap;
    for (auto &face : ringFaces) {
        for (auto &vertex : data(face).vertrices) {
            auto newCoor = ReCalculate2DCoordinates(currentPoint2DMap, vertex);
            if (!newCoor.has_value()) break;
            deletePoint2DMap[vertex] = newCoor.value();
        }
    }

    delete_vertex(deleteVertex);

    BarycentricCoordinates barycentricCoordinates;
    std::vector<std::array<VertexHandle, 3>> newFaces;

    int trangleIdx = CDTTrangle(coordinates, newFaces, barycentricCoordinates);
    std::vector<FaceHandle> facesHandle;
    if (trangleIdx != -1 && TryToAddFace(newFaces, facesHandle)) {
        data(facesHandle[trangleIdx]).vertrices.push_back(deleteVertex);
        data(deleteVertex).barycentricCoordinates = barycentricCoordinates;

        // TOOD(45degree): 判断被删除的点在那个三角形内部, 并计算重心坐标

        for (auto &[vertex, point2D] : deletePoint2DMap) {
            auto newParma = UpdateParam(newFaces, currentPoint2DMap, point2D);
            if (!newParma.has_value()) continue;
            auto [idx, alpha, beta, gamma] = newParma.value();
            data(facesHandle[idx]).vertrices.push_back(vertex);

            BarycentricCoordinates baryCoor;
            baryCoor[0].first = newFaces[idx][0], baryCoor[0].second = alpha;
            baryCoor[1].first = newFaces[idx][1], baryCoor[1].second = beta;
            baryCoor[2].first = newFaces[idx][2], baryCoor[2].second = gamma;
            data(vertex).barycentricCoordinates = baryCoor;
        }
        return;
    }
    for (int i = 0; i < coordinates.size(); i++) {
        newFaces.clear();
        trangleIdx = MVTTrangle(coordinates, i, newFaces, barycentricCoordinates);
        if (TryToAddFace(newFaces, facesHandle)) {

            // TOOD(45degree): 判断被删除的点在那个三角形内部, 并计算重心坐标

            data(facesHandle[trangleIdx]).vertrices.push_back(deleteVertex);
            data(deleteVertex).barycentricCoordinates = barycentricCoordinates;

            for (auto &[vertex, point2D] : deletePoint2DMap) {
                auto newParma = UpdateParam(newFaces, currentPoint2DMap, point2D);
                if (!newParma.has_value()) continue;
                auto [idx, alpha, beta, gamma] = newParma.value();
                data(facesHandle[idx]).vertrices.push_back(vertex);

                BarycentricCoordinates baryCoor;
                baryCoor[0].first = newFaces[idx][0], baryCoor[0].second = alpha;
                baryCoor[1].first = newFaces[idx][1], baryCoor[1].second = beta;
                baryCoor[2].first = newFaces[idx][2], baryCoor[2].second = gamma;
                data(vertex).barycentricCoordinates = baryCoor;
            }
            return;
        }
    }
    throw std::runtime_error("can't add face");
}

bool MapMesh::TryToAddFace(std::vector<std::array<VertexHandle, 3>> &faces,
                           std::vector<FaceHandle> &addedFaces) {
    addedFaces.clear();
    for (auto &face : faces) {
        for (int i = 0, ii = 1; i < 3; i++, ii++, ii %= 3) {
            auto half_edge = find_halfedge(face[i], face[ii]);
            if (half_edge.is_valid() && !is_boundary(half_edge)) {
                std::swap(face[1], face[2]);
                break;
            }
        }
        auto faceHandle = add_face(std::vector<VertexHandle>(face.begin(), face.end()));
        if (faceHandle.is_valid()) {
            addedFaces.push_back(faceHandle);
        } else {
            for (auto &addedFace : addedFaces) {
                delete_face(addedFace);
            }
            return false;
        }
    }
    return true;
}

void MapMesh::Initialize() {
    request_face_status();
    request_vertex_status();
    request_edge_status();

    for (auto vertexIter = vertices_sbegin(); vertexIter != vertices_end(); vertexIter++) {
        data(*vertexIter).canBeDeleted = true;
    }

    originFaces = std::vector<FaceHandle>(faces_sbegin(), faces_end());
    for (const auto &face : originFaces) {
        originFaceVertices[face] = std::vector<VertexHandle>(fv_begin(face), fv_end(face));
    }

    CalculateCurvature();
    CalculateArea();
    CalculateWeight(0.5);
}

void MapMesh::CalculateCurvature() {
    Eigen::MatrixX3d V, PD1, PD2;
    Eigen::MatrixX3i F;
    Eigen::VectorXd PV1, PV2, PV;

    OpenMesh2IGL(this, V, F);
    igl::principal_curvature(V, F, PD1, PD2, PV1, PV2);
    PV = PV1 + PV2;

    int j = 0;
    maxCurvature = std::numeric_limits<double>::min();
    for (auto i = vertices_sbegin(); i != vertices_end(); i++, j++) {
        data(*i).curvature = PV[j];
        if (PV[j] > maxCurvature) {
            maxCurvature = PV[j];
        }
    }
}

void MapMesh::CalculateWeight(double lambda) {
    while (!curvatureQueue.empty()) {
        curvatureQueue.pop();
    }

    for (auto i = vertices_sbegin(); i != vertices_end(); i++) {
        data(*i).weight = lambda * data(*i).curvature / maxCurvature +
                          (1 - lambda) * data(*i).ringArea / maxRingArea;
        curvatureQueue.push(*i);
    }
}

void MapMesh::CalculateArea() {
    maxRingArea = std::numeric_limits<double>::min();
    for (auto faceIter = faces_sbegin(); faceIter != faces_end(); faceIter++) {
        std::vector<VertexHandle> vertices(fv_begin(*faceIter), fv_end(*faceIter));
        assert(vertices.size() == 3);
        auto Vec1 = point(vertices[1]) - point(vertices[0]);
        auto Vec2 = point(vertices[2]) - point(vertices[0]);
        double area = 0.5 * Vec1.cross(Vec2).norm();
        data(*faceIter).area = area;
    }

    for (auto vertexIter = vertices_sbegin(); vertexIter != vertices_end(); vertexIter++) {
        double ringArea = 0;
        for (auto faceIter = vf_begin(*vertexIter); faceIter != vf_end(*vertexIter); faceIter++) {
            ringArea += data(*faceIter).area;
        }
        if (ringArea > maxRingArea) {
            maxRingArea = ringArea;
        }
        data(*vertexIter).ringArea = ringArea;
    }
}

double MapMesh::CalculateAngle(VertexHandle vertexHandle, FaceHandle face) {
    auto point1 = point(vertexHandle);
    Eigen::Vector3d vertex;
    vertex[0] = point1[0];
    vertex[1] = point1[1];
    vertex[2] = point1[2];

    std::vector<Eigen::Vector3d> points;
    for (auto _pointIter = fv_begin(face); _pointIter != fv_end(face); _pointIter++) {
        if (*_pointIter != vertexHandle) {
            auto _point = point(*_pointIter);
            Eigen::Vector3d p;
            p[0] = _point[0];
            p[1] = _point[1];
            p[2] = _point[2];
            points.emplace_back(p);
        }
    }
    assert(points.size() == 2);

    auto Vec1 = points[0] - vertex;
    auto Vec2 = points[1] - vertex;

    return std::acos(Vec1.normalized().dot(Vec2.normalized()));
}

void MapMesh::MapFaceFromOriginMesh(const FaceHandle &face, std::array<Point, 3> &mapFace) {
    const std::vector<VertexHandle> &vertices = originFaceVertices[face];
    assert(vertices.size() == 3);
    for (int i = 0; i < 3; i++) {
        auto vertexHandle = vertices[i];
        if (!IsVertexDeleted(vertexHandle)) {
            mapFace[i] = point(vertexHandle);
            continue;
        }

        auto baryCoor = data(vertexHandle).barycentricCoordinates.value();
        mapFace[i] = point(baryCoor[0].first) * baryCoor[0].second +
                     point(baryCoor[1].first) * baryCoor[1].second +
                     point(baryCoor[2].first) * baryCoor[2].second;
    }
}

void MapMesh::FaceSubDivision() {
    std::vector<FaceHandle> originFaceHandle(faces_sbegin(), faces_end());

    for (const auto &faceHandle : originFaceHandle) {
        std::vector<VertexHandle> vertices(fv_begin(faceHandle), fv_end(faceHandle));
        assert(vertices.size() == 3);

        for (int i = 0, ii = 1; i < 3; i++, ii++, ii %= 3) {
            auto newPoint = (point(vertices[i]) + point(vertices[ii])) / 2.0;
            auto newPointHandle = add_vertex(newPoint);
            data(newPointHandle).isNew = true;
            vertices.push_back(newPointHandle);
        }
        delete_face(faceHandle);
        add_face({vertices[0], vertices[3], vertices[5]});
        add_face({vertices[3], vertices[1], vertices[4]});
        add_face({vertices[5], vertices[4], vertices[2]});
        add_face({vertices[3], vertices[4], vertices[5]});
    }
}

void MapMesh::DownSampling() {
    while (!curvatureQueue.empty()) {
        auto vertexHandle = curvatureQueue.top();
        curvatureQueue.pop();
        if (!data(vertexHandle).canBeDeleted) continue;

        for (auto ringVertex = vv_begin(vertexHandle); ringVertex != vv_end(vertexHandle);
             ringVertex++) {
            data(*ringVertex).canBeDeleted = false;
        }

        // 计算2维坐标
        Coordinate2Dpair coordinates;
        std::vector<std::vector<VertexHandle>> faces;
        Calculate2D(vertexHandle, coordinates);

        ReTrangleAndAddFace(vertexHandle);
    }
}

void MapMesh::Remesh() {

    for (const auto &vertex : vertices()) {
        if (!data(vertex).isNew) continue;

        for (const auto &face : originFaces) {
            std::array<Point, 3> mapFace;
            const std::vector<VertexHandle> &fv = originFaceVertices[face];
            MapFaceFromOriginMesh(face, mapFace);
            if (IsInTriangle(point(vertex), mapFace)) {
                auto [alpha, beta, gamma] = CalculateBaryCoor(point(vertex), mapFace);
                Point newPoint = alpha * point(fv[0]) + beta * point(fv[1]) + gamma * point(fv[2]);
                point(vertex) = newPoint;
            }
        }
    }
}

std::optional<std::tuple<int, double, double, double>> MapMesh::UpdateParam(
    const std::vector<std::array<VertexHandle, 3>> &faces,
    const std::map<VertexHandle, Point2D> &point2DMap, const Point2D &point) {

    for (int i = 0; i < faces.size(); i++) {
        const auto &face = faces[i];
        std::array<Point2D, 3> triangle;
        triangle[0] = point2DMap.at(face[0]);
        triangle[1] = point2DMap.at(face[1]);
        triangle[2] = point2DMap.at(face[2]);
        if (!IsInTriangle(point, triangle)) continue;

        auto [alpha, beta, gamma] = CalculateBaryCoor(point, triangle);
        return std::make_tuple(i, alpha, beta, gamma);
    }

    return std::nullopt;
}

void MapMesh::Calculate2D(VertexHandle vertex, Coordinate2Dpair &coordinates) {
    coordinates.clear();
    auto vertex3D = point(vertex);

    double K_i = 0;
    for (auto faceIter = vf_begin(vertex); faceIter != vf_end(vertex); faceIter++) {
        K_i += CalculateAngle(vertex, *faceIter);
    }
    double a = 2 * M_PI / K_i;

    K_i = 0;

    VertexHandle lastVertex;
    for (auto ringPointIter = vv_begin(vertex); ringPointIter != vv_end(vertex); ringPointIter++) {
        if (ringPointIter != vv_begin(vertex)) {
            auto Vec1 = point(lastVertex) - point(vertex);
            auto Vec2 = point(*ringPointIter) - point(vertex);

            K_i += std::acos(Vec1.normalized().dot(Vec2.normalized()));
        }

        auto point3D = point(*ringPointIter);
        double r_k = (point3D - vertex3D).norm();

        Point2D point2D{std::pow(r_k, a) * std::cos(K_i * a), std::pow(r_k, a) * std::sin(K_i * a)};

        coordinates.emplace_back(*ringPointIter, point2D);
        lastVertex = *ringPointIter;
    }
}

}  // namespace Maps
