#ifndef MAPS_MAPS_MESH_H
#define MAPS_MAPS_MESH_H

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <array>
#include <functional>
#include <optional>
#include <queue>

namespace Maps {

using BarycentricCoordinates = std::array<std::pair<OpenMesh::VertexHandle, double>, 3>;

struct MapTrait : public OpenMesh::DefaultTraits {
   public:
    FaceTraits {
        double area;                                    // 面积
        std::vector<OpenMesh::VertexHandle> vertrices;  // 所有与该面相关连的顶点
    };

    VertexTraits {
        double curvature;           // 主曲率
        bool canBeDeleted = false;  // 是否能够被删除
        double weight;              // 权重
        double ringArea;            // 1领域三角面片的面积之和
        double isDeleted = false;   // 是否已经被删除
        bool isTranversed = false;  // 在更新重心坐标时是否被遍历过
        std::optional<BarycentricCoordinates> barycentricCoordinates;  // 重心坐标
    };
};

using BaseMesh = OpenMesh::TriMesh_ArrayKernelT<MapTrait>;
class MapMesh : public BaseMesh {
   private:
    using Point2D = OpenMesh::VectorT<double, 2>;
    using Coordinate2Dpair = std::vector<std::pair<VertexHandle, Point2D>>;
    using Triangles = std::vector<std::array<std::pair<VertexHandle, Point2D>, 3>>;

   private:
    struct CompareFun {
       private:
        MapMesh* _mesh;

       public:
        explicit CompareFun(MapMesh* mesh) : _mesh(mesh){};
        bool operator()(const OpenMesh::VertexHandle& L, const OpenMesh::VertexHandle& R) {
            return _mesh->data(L).weight < _mesh->data(R).weight;
        }
    };

   public:
    MapMesh() : curvatureQueue(CompareFun(this)) {}
    ~MapMesh() override = default;

    explicit MapMesh(const char* filename) : curvatureQueue(CompareFun(this)) {
        OpenMesh::IO::read_mesh(*this, filename);
    }
    MapMesh(const MapMesh& mapMesh) = default;

   public:
    void Initialize();

    void DownSampling();

    // TODO(45degree): need to impl
    void Remesh();

   private:
    /**
     * @brief 计算所有未删除点的曲率
     */
    void CalculateCurvature();

    /**
     * @brief 计算每个面的面积以及每个点的1领域面的面积之和
     */
    void CalculateArea();

    /**
     * @brief 计算每个点的权重
     */
    void CalculateWeight(double lambda = 0.5);

    /**
     * @brief 计算一个面中一个顶点的角度
     */
    double CalculateAngle(VertexHandle vertex, FaceHandle face);

    /**
     * @brief 重新三角化并将三角化后的结果添加到面中
     *
     * @param coordinates
     * @param deleteVertex
     */
    void ReTrangleAndAddFace(const VertexHandle& deleteVertex);

    /**
     * @brief 尝试添加面
     *
     * @param[in]  faces
     * @param[out] facesHandle
     *
     * @return true: 添加成功
     *         false: 添加失败
     */
    bool TryToAddFace(std::vector<std::array<VertexHandle, 3>>& faces,
                      std::vector<FaceHandle>& facesHandle);

    /**
     * @brief 计算1领域点的2维坐标
     *
     * @param[in]  vertex
     * @param[out] coordinates
     */
    void Calculate2D(VertexHandle vertex, Coordinate2Dpair& coordinates);

    /**
     * @brief 将每个面细分为四个面
     */
    void FaceSubDivision();

    /**
     * @brief 更新重心坐标
     */
    void UpdateBarycentricCoordinates();

    /**
     * @brief 判断顶点是否被删除
     */
    [[nodiscard]] bool IsVertexDeleted(const MapMesh::VertexHandle& vertexHandle) const {
        return data(vertexHandle).barycentricCoordinates.has_value();
    }

    /**
     * @brief 重新计算被删除点的2维坐标
     *
     * @param originCoordinates
     * @param deleteVertex
     *
     * @return
     */
    std::optional<Point2D> ReCalculate2DCoordinates(
        std::map<VertexHandle, Point2D>& originCoordinates, VertexHandle deleteVertex);

    /** TODO(45degree): 需要重新实现
     * @brief CDT三角化
     *
     * @param[in]  coordinates
     * @param[out] faces
     * @param[out] barycentricCoordinates
     *
     * @return i != -1: 原点存在于三角化后的第i个面中
     *         i == -1: 三角化失败
     */
    static int CDTTrangle(const Coordinate2Dpair& coordinates,
                          std::vector<std::array<VertexHandle, 3>>& faces,
                          BarycentricCoordinates& barycentricCoordinates);

    /** TODO(45degree): 需要重新实现
     * @brief MVT三角化
     *
     * @param[in]  coordinates
     * @param[in]  startIdx
     * @param[out] faces
     * @param[out] barycentricCoordinates
     *
     * @return i != -1: 原点存在于三角化后的第i个面中
     *         i == -1: 三角化失败
     */
    static int MVTTrangle(const Coordinate2Dpair& coordinates, int startIdx,
                          std::vector<std::array<VertexHandle, 3>>& faces,
                          BarycentricCoordinates& barycentricCoordinates);

    static bool IsInTriangle(const Point2D& point, const std::array<Point2D, 3>& triangle2D);

    static std::tuple<double, double, double> CalculateBaryCoor(
        const Point2D& point, const std::array<Point2D, 3>& triangle2D);

    static std::optional<std::tuple<int, double, double, double>> UpdateParam(
        const std::vector<std::array<VertexHandle, 3>>& faces,
        const std::map<VertexHandle, Point2D>& point2DMap, const Point2D& point);

   private:
    double maxRingArea;
    double maxCurvature;
    std::priority_queue<VertexHandle, std::vector<VertexHandle>, CompareFun> curvatureQueue;
};
}  // namespace Maps

#endif
