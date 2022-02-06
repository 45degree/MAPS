#include "MapsMesh.h"

#include <igl/gaussian_curvature.h>
#include <poly2tri.h>
#include <tbb/tbb.h>

#include <Eigen/Dense>
#include <algorithm>
#include <limits>

namespace Maps {

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

void MapMesh::ReadMeshFromLibigl(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F) {
    std::vector<VertexHandle> vertexHandle(V.rows());
    for (int i = 0; i < V.rows(); i++) {
        vertexHandle[i] = add_vertex(Point(V(i, 0), V(i, 1), V(i, 2)));
    }

    for (int i = 0; i < F.rows(); i++) {
        add_face({vertexHandle[F(i, 0)], vertexHandle[F(i, 1)], vertexHandle[F(i, 2)]});
    }
}

int MapMesh::CDTTrangle(const Coordinate2DPair &coordinates,
                        std::vector<std::array<VertexHandle, 3>> &faces,
                        BaryCoor &barycentricCoordinates) {

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

int MapMesh::MVTTrangle(const Coordinate2DPair &coordinates, int startIdx,
                        std::vector<std::array<VertexHandle, 3>> &faces,
                        BaryCoor &barycentricCoordinates) {
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

        trangleIdx = static_cast<int>(faces.size() - 1);
    }
    return trangleIdx;
}

std::optional<MapMesh::Point2D> MapMesh::ReCalculate2DCoordinates(
    std::map<VertexHandle, Point2D> &originCoor, VertexHandle deleteVertex) {
    if (!IsVertexDeleted(deleteVertex)) return std::nullopt;
    auto baryCoor = data(deleteVertex).baryCoor.value();

    return originCoor[baryCoor[0].first] * baryCoor[0].second +
           originCoor[baryCoor[1].first] * baryCoor[1].second +
           originCoor[baryCoor[2].first] * baryCoor[2].second;
}

void MapMesh::ReTrangleAndAddFace(const VertexHandle &deleteVertex) {
    // 记录1领域面中所有包含的点, 并更新这些点的参数
    // 假设1领域面f上记录了一个点v, v中记录了v在f上的重心坐标(a, b, c), 面f在经过2维映射变成了f'
    // 更新过程: 1. 根据f'和(a, b, c)计算v的二维映射点v'
    //           2. 重新计算包含v'的三角形p‘, 并根据p’重新计算坐标(a', b', c')
    //           3. 删除f, 添加面p, 更新v的重心坐标, 建立p与v的关联

    std::vector<FaceHandle> ringFaces(vf_begin(deleteVertex), vf_end(deleteVertex));
    Coordinate2DPair coordinates;
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

    delete_vertex(deleteVertex, false);

    BaryCoor barycentricCoordinates;
    std::vector<std::array<VertexHandle, 3>> newFaces;

    int trangleIdx = CDTTrangle(coordinates, newFaces, barycentricCoordinates);
    std::vector<FaceHandle> facesHandle;
    if (trangleIdx != -1 && TryToAddFaces(std::move(newFaces), facesHandle)) {
        data(facesHandle[trangleIdx]).vertrices.push_back(deleteVertex);
        data(deleteVertex).baryCoor = barycentricCoordinates;

        // 判断被删除的点在那个三角形内部, 并计算重心坐标

        for (auto &[vertex, point2D] : deletePoint2DMap) {
            auto newParma = UpdateParam(newFaces, currentPoint2DMap, point2D);
            if (!newParma.has_value()) continue;
            auto [idx, alpha, beta, gamma] = newParma.value();
            data(facesHandle[idx]).vertrices.push_back(vertex);

            BaryCoor baryCoor;
            baryCoor[0].first = newFaces[idx][0], baryCoor[0].second = alpha;
            baryCoor[1].first = newFaces[idx][1], baryCoor[1].second = beta;
            baryCoor[2].first = newFaces[idx][2], baryCoor[2].second = gamma;
            data(vertex).baryCoor = baryCoor;
        }

        return;
    }
    for (int i = 0; i < coordinates.size(); i++) {
        newFaces.clear();
        trangleIdx = MVTTrangle(coordinates, i, newFaces, barycentricCoordinates);
        if (TryToAddFaces(std::move(newFaces), facesHandle)) {

            data(facesHandle[trangleIdx]).vertrices.push_back(deleteVertex);
            data(deleteVertex).baryCoor = barycentricCoordinates;

            for (auto &[vertex, point2D] : deletePoint2DMap) {
                auto newParma = UpdateParam(newFaces, currentPoint2DMap, point2D);
                if (!newParma.has_value()) continue;
                auto [idx, alpha, beta, gamma] = newParma.value();
                data(facesHandle[idx]).vertrices.push_back(vertex);

                BaryCoor baryCoor;
                baryCoor[0].first = newFaces[idx][0], baryCoor[0].second = alpha;
                baryCoor[1].first = newFaces[idx][1], baryCoor[1].second = beta;
                baryCoor[2].first = newFaces[idx][2], baryCoor[2].second = gamma;
                data(vertex).baryCoor = baryCoor;
            }
            return;
        }
    }
    throw std::runtime_error("can't add face");
}

void MapMesh::Initialize() {
    originFaces = std::vector<FaceHandle>(faces_sbegin(), faces_end());
    for (const auto &face : originFaces) {
        originFaceVertices[face] = std::vector<VertexHandle>(fv_begin(face), fv_end(face));
    }
}

void MapMesh::CalculateCurvature() {
    Eigen::MatrixX3d V;
    Eigen::MatrixX3i F;
    Eigen::VectorXd K;

    OpenMesh2IGL(this, V, F);
    igl::gaussian_curvature(V, F, K);

    int j = 0;
    maxCurvature = (std::numeric_limits<double>::min)();
    for (auto i = vertices_sbegin(); i != vertices_end(); i++, j++) {
        data(*i).curvature = K[j];
        if (K[j] > maxCurvature) {
            maxCurvature = K[j];
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

void MapMesh::CalculateAreas() {
    maxRingArea = (std::numeric_limits<double>::min)();
    for (auto faceIter = faces_sbegin(); faceIter != faces_end(); faceIter++) {
        double area = CalculateArea(*faceIter);
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

void MapMesh::MapFaceFromOriginMesh(const FaceHandle &face, std::array<Point, 3> &mapFace) {
    const std::vector<VertexHandle> &vertices = originFaceVertices[face];
    assert(vertices.size() == 3);
    for (int i = 0; i < 3; i++) {
        auto vertexHandle = vertices[i];
        if (!IsVertexDeleted(vertexHandle)) {
            mapFace[i] = point(vertexHandle);
            continue;
        }

        auto baryCoor = data(vertexHandle).baryCoor.value();
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
    request_face_status();
    request_vertex_status();
    request_edge_status();

    for (auto vertexIter = vertices_sbegin(); vertexIter != vertices_end(); vertexIter++) {
        data(*vertexIter).canBeDeleted = true;
    }

    CalculateCurvature();
    CalculateAreas();
    CalculateWeight(0.5);
    while (!curvatureQueue.empty()) {
        auto vertexHandle = curvatureQueue.top();
        curvatureQueue.pop();
        if (!data(vertexHandle).canBeDeleted) continue;

        for (auto ringVertex = vv_begin(vertexHandle); ringVertex != vv_end(vertexHandle);
             ringVertex++) {
            data(*ringVertex).canBeDeleted = false;
        }

        ReTrangleAndAddFace(vertexHandle);
    }
}

void MapMesh::Remesh() {

    int N = static_cast<int>(std::distance(vertices_begin(), vertices_end()));

    int i = 0;
    std::mutex _mutex;

    auto IsMapedInASigleFace = [&](const FaceHandle &faceHandle) -> bool {
        std::set<VertexHandle> verticesHandle;
        for (auto vertexIter = fv_begin(faceHandle); vertexIter != fv_end(faceHandle);
             vertexIter++) {
            if (data(*vertexIter).baryCoor.has_value()) {
                auto baryCoor = data(*vertexIter).baryCoor.value();
                verticesHandle.insert(baryCoor[0].first);
                verticesHandle.insert(baryCoor[1].first);
                verticesHandle.insert(baryCoor[2].first);
            }
        }

        return verticesHandle.size() == 3;
    };

    auto findMaxCommonVertex = [&](const std::vector<FaceHandle> &facesHandle) {
        std::multiset<VertexHandle> verticesHandle;
        for (const auto &faceHandle : facesHandle) {
            for (const auto &vertexHandle : fv_range(faceHandle)) {
                verticesHandle.insert(vertexHandle);
            }
        }

        VertexHandle maxCommonVertex;
        int maxCount = -1;
        for (const auto &vertex : verticesHandle) {
            if (verticesHandle.count(vertex) > maxCount) {
                maxCount = static_cast<int>(verticesHandle.count(vertex));
                maxCommonVertex = vertex;
            }
        }

        return maxCommonVertex;
    };

    oneapi::tbb::parallel_for_each(vertices_sbegin(), vertices_end(),
                                   [&](const VertexHandle &vertex) {
        if (!data(vertex).isNew) {
            _mutex.lock();
            i++;
            _mutex.unlock();
            return;
        }
        oneapi::tbb::parallel_for_each(originFaces.begin(), originFaces.end(),
                                       [&](const FaceHandle &face) {
            std::array<Point, 3> mapFace;
            const std::vector<VertexHandle> &fv = originFaceVertices[face];
            if (IsMapedInASigleFace(face)) {  // 原面能够完全映射到基面上
                MapFaceFromOriginMesh(face, mapFace);
                if (IsInTriangle(point(vertex), mapFace)) {
                    auto [alpha, beta, gamma] = CalculateBaryCoor(point(vertex), mapFace);
                    Point newPoint =
                        alpha * point(fv[0]) + beta * point(fv[1]) + gamma * point(fv[2]);
                    point(vertex) = newPoint;
                    return;
                }
            } else {
                // 在基面上找到公共顶点，并展开到2维平面计算重心坐标
                std::vector<FaceHandle> baseFaces;
                auto vertexFace = FindFace(point(vertex)).value();
                for (const auto &vertexHandle : fv) {
                    if (!IsVertexDeleted(vertexHandle)) continue;

                    auto baryCoor = data(vertexHandle).baryCoor.value();
                    auto face = FindFace({baryCoor[0].first, baryCoor[1].first, baryCoor[2].first});
                    if (face.has_value()) {
                        baseFaces.push_back(face.value());
                    }
                }
                auto commonVertex = findMaxCommonVertex(baseFaces);
                Coordinate2DPair coordinates;
                Calculate2D(commonVertex, coordinates);
                Point2D vertex2D;

                std::array<std::pair<VertexHandle, Point2D>, 3> face2D;
                for (int i = 0; i < 3; i++) {
                    auto vertexHandle = fv[i];
                    auto v = std::find_if(
                        coordinates.begin(), coordinates.end(),
                        [&vertexHandle](const std::pair<VertexHandle, Point2D> &cooridinate) {
                        return cooridinate.first == vertexHandle;
                        });
                    face2D[i] = *v;
                }
            }
        });
        _mutex.lock();
        i++;
        _mutex.unlock();
    });
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

void MapMesh::Calculate2D(VertexHandle vertex, Coordinate2DPair &coordinates) {
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
