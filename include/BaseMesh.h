#ifndef MAPS_BASE_MESH_H
#define MAPS_BASE_MESH_H

#include <Eigen/Dense>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <optional>

namespace Maps {

template <typename Trait = OpenMesh::DefaultTraits>
class BaseMesh : public OpenMesh::TriMesh_ArrayKernelT<Trait> {
   public:
    using VertexHandle = typename OpenMesh::TriMesh_ArrayKernelT<Trait>::VertexHandle;
    using FaceHandle = typename OpenMesh::TriMesh_ArrayKernelT<Trait>::FaceHandle;
    using Point = typename OpenMesh::TriMesh_ArrayKernelT<Trait>::Point;
    using HalfedgeHandle = typename OpenMesh::TriMesh_ArrayKernelT<Trait>::HalfedgeHandle;

   public:
    using Point2D = OpenMesh::VectorT<double, 2>;
    using Coordinate2DPair = std::vector<std::pair<VertexHandle, Point2D>>;
    using Coordinate2DMap = std::map<VertexHandle, Point2D>;

    using CoordinatePair = std::vector<std::pair<VertexHandle, Point>>;
    using Triangles = std::vector<std::array<std::pair<VertexHandle, Point2D>, 3>>;

   public:
    BaseMesh() = default;
    ~BaseMesh() = default;
    BaseMesh(const BaseMesh& mesh) = default;

   public:
    /**
     * @brief 计算一个面的面积
     */
    double CalculateArea(const FaceHandle& faceHandle);

    /**
     * @brief 计算一个面中一个顶点的角度
     */
    double CalculateAngle(const VertexHandle& vertex, const FaceHandle& face);

    /**
     * @brief 尝试添加一个面
     *
     * @param[in]  faces
     * @param[out] facesHandle
     *
     * @return true: 添加成功
     *         false: 添加失败
     */
    bool TryToAddFace(std::array<VertexHandle, 3>&& face, FaceHandle& faceHandle);

    /**
     * @brief 尝试添加多个面
     *
     * @param[in]  faces
     * @param[out] facesHandle
     *
     * @return true: 添加成功
     *         false: 添加失败
     */
    bool TryToAddFaces(std::vector<std::array<VertexHandle, 3>>&& faces,
                       std::vector<FaceHandle>& facesHandle);

    std::optional<FaceHandle> FindFace(const std::array<VertexHandle, 3>& faceVertices) const {
        auto half_edge = this->find_halfedge(faceVertices[0], faceVertices[1]);
        if (half_edge.next().to() == faceVertices[2]) {
            return half_edge.face();
        } else if (half_edge.opp().next().to() == faceVertices[2]) {
            return half_edge.opp().face();
        }
        return std::nullopt;
    }

    std::optional<FaceHandle> FindFace(const Point& point) {
        for (const auto& face : this->faces()) {
            std::vector<VertexHandle> triangle(this->fv_begin(face), this->fv_end(face));
            if (IsInTriangle(point, {this->point(triangle[0]), this->point(triangle[1]),
                                     this->point(triangle[2])})) {
                return face;
            }
        }
        return std::nullopt;
    }

    std::vector<std::pair<VertexHandle, double>> CalculateBaryCoor(const Point& point,
                                                                   const FaceHandle& face);

    /**
     * @brief 判断2维点是否在三角形内部
     */
    static bool IsInTriangle(const Point2D& point, const std::array<Point2D, 3>& triangle2D);

    /**
     * @brief 判断3维点是否在三角形内部(如果4点不共面，也不算在内部)
     */
    static bool IsInTriangle(const Point& point, const std::array<Point, 3>& triangle);

    /**
     * @brief 判断点是否和三角形在同一平面上
     */
    static bool IsCoplanar(const Point& point, const std::array<Point, 3>& triangle);

    /**
     * @brief 在二维平面上计算重心坐标
     *
     * @param point
     * @param triangle2D
     *
     * @return
     */
    static std::tuple<double, double, double> CalculateBaryCoor(
        const Point2D& point, const std::array<Point2D, 3>& triangle2D);

    /**
     * @brief 在三维平面上计算重心坐标
     *
     * @param point
     * @param triangle
     *
     * @return
     */
    static std::tuple<double, double, double> CalculateBaryCoor(
        const Point& point, const std::array<Point, 3>& triangle);
};

template <typename T>
double BaseMesh<T>::CalculateArea(const FaceHandle& faceHandle) {

    std::vector<VertexHandle> vertices(this->fv_begin(faceHandle), this->fv_end(faceHandle));
    assert(vertices.size() == 3);
    auto Vec1 = this->point(vertices[1]) - this->point(vertices[0]);
    auto Vec2 = this->point(vertices[2]) - this->point(vertices[0]);
    double area = 0.5 * Vec1.cross(Vec2).norm();

    return area;
}

template <typename T>
double BaseMesh<T>::CalculateAngle(const VertexHandle& vertexHandle, const FaceHandle& faceHandle) {

    auto point1 = this->point(vertexHandle);
    Eigen::Vector3d vertex;
    vertex[0] = point1[0];
    vertex[1] = point1[1];
    vertex[2] = point1[2];

    std::vector<Eigen::Vector3d> points;
    for (auto _pointIter = this->fv_begin(faceHandle); _pointIter != this->fv_end(faceHandle);
         _pointIter++) {

        if (*_pointIter != vertexHandle) {
            auto _point = this->point(*_pointIter);
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

template <typename T>
bool BaseMesh<T>::TryToAddFace(std::array<VertexHandle, 3>&& face, FaceHandle& facesHandle) {

    for (int i = 0, ii = 1; i < 3; i++, ii++, ii %= 3) {
        auto half_edge = this->find_halfedge(face[i], face[ii]);
        if (half_edge.is_valid() && !this->is_boundary(half_edge)) {
            std::swap(face[1], face[2]);
            break;
        }
    }
    facesHandle = this->add_face(std::vector<VertexHandle>(face.begin(), face.end()));
    return facesHandle.is_valid();
}

template <typename T>
bool BaseMesh<T>::TryToAddFaces(std::vector<std::array<VertexHandle, 3>>&& faces,
                                std::vector<FaceHandle>& facesHandle) {
    facesHandle.clear();
    for (auto& face : faces) {
        for (int i = 0, ii = 1; i < 3; i++, ii++, ii %= 3) {
            auto half_edge = this->find_halfedge(face[i], face[ii]);
            if (half_edge.is_valid() && !this->is_boundary(half_edge)) {
                std::swap(face[1], face[2]);
                break;
            }
        }
        auto faceHandle = this->add_face(std::vector<VertexHandle>(face.begin(), face.end()));
        if (faceHandle.is_valid()) {
            facesHandle.push_back(faceHandle);
        } else {
            for (auto& addedFace : facesHandle) {
                this->delete_face(addedFace);
            }
            return false;
        }
    }
    return true;
}

template <typename T>
bool BaseMesh<T>::IsInTriangle(const Point2D& point, const std::array<Point2D, 3>& triangle2D) {
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

template <typename T>
bool BaseMesh<T>::IsInTriangle(const Point& point, const std::array<Point, 3>& triangle) {

    if (!IsCoplanar(point, triangle)) return false;

    auto Vec1 = (point - triangle[0]).cross(triangle[1] - triangle[0]);
    auto Vec2 = (point - triangle[1]).cross(triangle[2] - triangle[1]);
    auto Vec3 = (point - triangle[2]).cross(triangle[0] - triangle[2]);

    double z1 = Vec1.dot(Vec2);
    double z2 = Vec2.dot(Vec3);
    double z3 = Vec3.dot(Vec1);

    return z1 >= 0 && z2 >= 0 && z3 >= 0;
}

template <typename T>
bool BaseMesh<T>::IsCoplanar(const Point& point, const std::array<Point, 3>& triangle) {
    double esp = 1e-5;
    auto Vec1 = triangle[2] - triangle[0];
    auto Vec2 = triangle[1] - triangle[0];
    auto Vec3 = point - triangle[0];

    double result = Vec1.cross(Vec2).dot(Vec3);
    return result >= -esp && result <= esp;
}

template <typename T>
std::tuple<double, double, double> BaseMesh<T>::CalculateBaryCoor(
    const Point2D& point, const std::array<Point2D, 3>& triangle2D) {

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

template <typename T>
std::tuple<double, double, double> BaseMesh<T>::CalculateBaryCoor(
    const Point& point, const std::array<Point, 3>& triangle) {

    double area = (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]).norm();
    double area1 = (triangle[1] - point).cross(triangle[2] - point).norm();
    double area2 = (triangle[0] - point).cross(triangle[2] - point).norm();
    double area3 = (triangle[1] - point).cross(triangle[0] - point).norm();

    return {area1 / area, area2 / area, area3 / area};
}

template <typename T>
std::vector<std::pair<typename BaseMesh<T>::VertexHandle, double>> BaseMesh<T>::CalculateBaryCoor(
    const Point& point, const FaceHandle& face) {

    std::array<Point, 3> triangle;
    std::array<VertexHandle, 3> vertexHandle;
    auto vertexIter = this->fv_begin(face);
    for (int i = 0; i < 3; i++, vertexIter++) {
        triangle[i] = this->point(*vertexIter);
        vertexHandle[i] = *vertexIter;
    }

    auto [alpha, beta, gamma] = CalculateBaryCoor(point, triangle);

    return {{vertexHandle[0], alpha}, {vertexHandle[1], beta}, {vertexHandle[2], gamma}};
}

}  // namespace Maps

#endif
