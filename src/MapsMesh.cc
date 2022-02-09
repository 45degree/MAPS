#include "MapsMesh.h"

#include <igl/gaussian_curvature.h>
#include <poly2tri.h>
#include <tbb/tbb.h>

#include <Eigen/Dense>
#include <algorithm>
#include <limits>
#include <progressbar/progressbar.hpp>

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
    std::map<VertexHandle, Point2D> &originCoor, VertexHandle deleteVertex) const {
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
    if (trangleIdx != -1 && TryToAddFaces(newFaces, facesHandle)) {
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
        if (TryToAddFaces(newFaces, facesHandle)) {

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

void MapMesh::MapFaceFromOriginMesh(const FaceHandle &face, std::array<Point, 3> &mapFace) const {
    if (originFaceVertices.find(face) == originFaceVertices.end()) return;

    const std::vector<VertexHandle> &vertices = originFaceVertices.at(face);
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
    SaveBaseLevelMesh();
}

void MapMesh::Remesh() {

    int N = static_cast<int>(std::distance(vertices_sbegin(), vertices_end()));
    progressbar bar(N);

    std::mutex _mutex;

    auto IsMapedInASingleFace = [&](const FaceHandle &faceHandle) -> bool {
        std::set<VertexHandle> verticesHandle;
        for (const auto &vertexHandle : originFaceVertices[faceHandle]) {
            if (data(vertexHandle).baryCoor.has_value()) {
                auto baryCoor = data(vertexHandle).baryCoor.value();
                verticesHandle.insert(baryCoor[0].first);
                verticesHandle.insert(baryCoor[1].first);
                verticesHandle.insert(baryCoor[2].first);
            } else {
                verticesHandle.insert(vertexHandle);
            }
        }

        return verticesHandle.size() == 3;
    };

    auto findMaxCommonVertex = [&](const std::vector<FaceHandle> &facesHandle) {
        std::multiset<VertexHandle> verticesHandle;
        for (const auto &faceHandle : facesHandle) {
            for (const auto &vertexHandle : baseLevelMesh->fv_range(faceHandle)) {
                verticesHandle.insert(vertexHandle);
            }
        }

        VertexHandle maxCommonVertex;
        int maxCount = -1;
        for (const auto &vertex : verticesHandle) {
            int count = static_cast<int>(verticesHandle.count(vertex));
            if (count > maxCount) {
                maxCount = static_cast<int>(verticesHandle.count(vertex));
                maxCommonVertex = vertex;
            }
        }

        return maxCommonVertex;
    };

    auto isRingFace = [&](const VertexHandle &vertexHandle, const FaceHandle &faceHandle) {
        for (const auto &face : baseLevelMesh->vf_range(vertexHandle)) {
            if (face == faceHandle) return true;
        }
        return false;
    };

    auto isRingPoint =
        [&](const VertexHandle &vertexHandle, const std::vector<VertexHandle> &ringVertices) {
        for (const auto &ringVertex : ringVertices) {
            bool found = false;
            for (const auto &realRingVertex : baseLevelMesh->vv_range(vertexHandle)) {
                if (realRingVertex == ringVertex) {
                    found = true;
                }
            }
            if (!found) return false;
        }
        return true;
    };

    oneapi::tbb::parallel_for_each(vertices_sbegin(), vertices_end(),
                                   [&](const VertexHandle &vertex) {
        if (!data(vertex).isNew) {
            _mutex.lock();
            bar.update();
            _mutex.unlock();
            /* continue; */
            return;
        }

        // 在基面上找到公共顶点，并展开到2维平面计算重心坐标
        auto vertexFaceOption = baseLevelMesh->FindFace(point(vertex));
        if (!vertexFaceOption.has_value()) {
            _mutex.lock();
            bar.update();
            _mutex.unlock();
            /* continue; */
            return;
        }
        auto vertexFace = vertexFaceOption.value();

        for (const auto &face : originFaces) {
            std::array<Point, 3> mapFace;
            const std::vector<VertexHandle> &fv = originFaceVertices[face];

            if (IsMapedInASingleFace(face)) {  // 原面能够完全映射到基面上
                MapFaceFromOriginMesh(face, mapFace);
                if (IsInTriangle(point(vertex), mapFace)) {
                    auto [alpha, beta, gamma] = CalculateBaryCoor(point(vertex), mapFace);
                    Point newPoint =
                        alpha * point(fv[0]) + beta * point(fv[1]) + gamma * point(fv[2]);
                    point(vertex) = newPoint;
                    break;
                }
            } else {

                std::vector<FaceHandle> baseFaces;
                std::vector<VertexHandle> fixVertex;
                for (const auto &vertexHandle : fv) {
                    if (!IsVertexDeleted(vertexHandle)) {
                        fixVertex.push_back(vertexHandle);
                        continue;
                    }

                    auto baryCoor = data(vertexHandle).baryCoor.value();
                    auto face = baseLevelMesh->FindFace(
                        {baryCoor[0].first, baryCoor[1].first, baryCoor[2].first});

                    if (face.has_value()) {
                        baseFaces.push_back(face.value());
                    }
                }

                VertexHandle commonVertex;
                if (baseFaces.size() == 3) {
                    if (std::find(baseFaces.begin(), baseFaces.end(), vertexFace) ==
                        baseFaces.end()) {
                        continue;
                    }
                    commonVertex = findMaxCommonVertex(baseFaces);
                } else {
                    bool foundCommonVertex = false;
                    for (const auto &baseFace : baseFaces) {
                        for (const auto &baseVert : baseLevelMesh->fv_range(baseFace)) {
                            if (isRingFace(baseVert, vertexFace) &&
                                isRingPoint(baseVert, fixVertex)) {
                                commonVertex = baseVert;
                                foundCommonVertex = true;
                                break;
                            }
                        }
                        if (foundCommonVertex) break;
                    }
                    if (!foundCommonVertex) continue;
                }

                Coordinate2DPair coorPair;
                baseLevelMesh->Calculate2D(commonVertex, coorPair);

                Coordinate2DMap coorMap(coorPair.begin(), coorPair.end());
                Point2D vertex2D(0, 0);
                auto vertex2DBaryCoor = baseLevelMesh->CalculateBaryCoor(point(vertex), vertexFace);
                for (const auto &vertex2DBaryCoorItem : vertex2DBaryCoor) {
                    if (coorMap.find(vertex2DBaryCoorItem.first) == coorMap.end()) {
                        continue;
                    }
                    vertex2D +=
                        coorMap.at(vertex2DBaryCoorItem.first) * vertex2DBaryCoorItem.second;
                }

                std::array<std::pair<VertexHandle, Point2D>, 3> face2D;
                bool canMapTo2D = true;
                for (int i = 0; i < 3; i++) {
                    face2D[i].first = fv[i];
                    face2D[i].second = Point2D(0, 0);
                    if (!IsVertexDeleted(fv[i])) {
                        face2D[i].second = coorMap[fv[i]];
                        continue;
                    }

                    auto baryCoor = data(fv[i]).baryCoor.value();
                    for (const auto &baryCoorItem : baryCoor) {
                        if (baryCoorItem.first == commonVertex) continue;
                        if (coorMap.find(baryCoorItem.first) == coorMap.end()) {
                            canMapTo2D = false;
                            break;
                            /* throw std::runtime_error("can't found item"); */
                        }
                        face2D[i].second += coorMap.at(baryCoorItem.first) * baryCoorItem.second;
                    }
                }
                if (!canMapTo2D) continue;

                std::array<Point2D, 3> triangle{face2D[0].second, face2D[1].second,
                                                face2D[2].second};

                if (!IsInTriangle(vertex2D, triangle)) continue;

                auto [alpha, beta, gamma] = CalculateBaryCoor(vertex2D, triangle);

                point(vertex) = alpha * point(face2D[0].first) + beta * point(face2D[1].first) +
                                gamma * point(face2D[2].first);
                break;
            }
        }
        _mutex.lock();
        bar.update();
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

void MapMesh::Calculate2D(VertexHandle vertex, Coordinate2DPair &coordinates) const {
    coordinates.clear();
    auto vertex3D = point(vertex);

    double K_i = 0;
    for (const auto &face : vf_range(vertex)) {
        K_i += CalculateAngle(vertex, face);
    }
    double a = 2 * M_PI / K_i;

    K_i = 0;

    VertexHandle lastVertex;
    int i = 0;
    for (const auto &ringPoint : vv_range(vertex)) {
        if (i != 0) {
            auto Vec1 = point(lastVertex) - point(vertex);
            auto Vec2 = point(ringPoint) - point(vertex);

            K_i += std::acos(Vec1.normalized().dot(Vec2.normalized()));
        }
        auto point3D = point(ringPoint);
        double r_k = (point3D - vertex3D).norm();

        Point2D point2D{std::pow(r_k, a) * std::cos(K_i * a), std::pow(r_k, a) * std::sin(K_i * a)};

        coordinates.emplace_back(ringPoint, point2D);
        lastVertex = ringPoint;

        i++;
    }
}

}  // namespace Maps
