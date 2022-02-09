#ifndef MAPS_MAPS_MESH_H
#define MAPS_MAPS_MESH_H

#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <array>
#include <functional>
#include <optional>
#include <queue>

#include "BaseMesh.h"

namespace Maps {

using BaryCoor = std::array<std::pair<OpenMesh::VertexHandle, double>, 3>;

struct MapTrait : public OpenMesh::DefaultTraits {
   public:
    FaceTraits {
        double area;                                    // 面积
        std::vector<OpenMesh::VertexHandle> vertrices;  // 所有与该面相关连的顶点
    };

    VertexTraits {
        double curvature;                  // 主曲率
        bool canBeDeleted = false;         // 是否能够被删除
        double weight;                     // 权重
        double ringArea;                   // 1领域三角面片的面积之和
        bool isTranversed = false;         // 在更新重心坐标时是否被遍历过
        bool isNew = false;                // 是否是细分出来的点
        std::optional<BaryCoor> baryCoor;  // 重心坐标
    };
};

class MapMesh : public BaseMesh<MapTrait> {

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
    ~MapMesh() override {
        if (baseLevelMesh != nullptr) {
            delete baseLevelMesh;
        }
    }

    MapMesh(const MapMesh& mapMesh) = default;

   public:
    void Initialize();

    void DownSampling();

    /**
     * @brief 将每个面细分为四个面
     */
    void FaceSubDivision();

    void Remesh();

   private:
    /**
     * @brief 计算所有未删除点的曲率
     */
    void CalculateCurvature();

    /**
     * @brief 计算每个面的面积以及每个点的1领域面的面积之和
     */
    void CalculateAreas();

    /**
     * @brief 计算每个点的权重
     */
    void CalculateWeight(double lambda = 0.5);

    /**
     * @brief 重新三角化并将三角化后的结果添加到面中
     *
     * @param coordinates
     * @param deleteVertex
     */
    void ReTrangleAndAddFace(const VertexHandle& deleteVertex);

    inline void SaveBaseLevelMesh() {
        auto* newMesh = new MapMesh(*this);
        newMesh->ClearBaseLevelMesh();
        if (baseLevelMesh != nullptr) {
            delete baseLevelMesh;
        }
        baseLevelMesh = newMesh;
    }

   public:
    inline void ClearBaseLevelMesh() { baseLevelMesh = nullptr; }

    /**
     * @brief 计算1领域点的2维坐标
     *
     * @param[in]  vertex
     * @param[out] coordinates
     */
    void Calculate2D(VertexHandle vertex, Coordinate2DPair& coordinates) const;

    /**
     * @brief 判断顶点是否被删除
     */
    [[nodiscard]] bool IsVertexDeleted(const MapMesh::VertexHandle& vertexHandle) const {
        return data(vertexHandle).baryCoor.has_value();
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
        std::map<VertexHandle, Point2D>& originCoordinates, VertexHandle deleteVertex) const;

    /**
     * @brief 将原始网格中的一个面映射到现网格中
     *
     * @param[in] face
     * @param[out] mapFace
     */
    void MapFaceFromOriginMesh(const FaceHandle& face, std::array<Point, 3>& mapFace) const;

    /**
     * @brief CDT三角化
     *
     * @param[in]  coordinates
     * @param[out] faces
     * @param[out] barycentricCoordinates
     *
     * @return i != -1: 原点存在于三角化后的第i个面中
     *         i == -1: 三角化失败
     */
    static int CDTTrangle(const Coordinate2DPair& coordinates,
                          std::vector<std::array<VertexHandle, 3>>& faces,
                          BaryCoor& barycentricCoordinates);

    /**
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
    static int MVTTrangle(const Coordinate2DPair& coordinates, int startIdx,
                          std::vector<std::array<VertexHandle, 3>>& faces,
                          BaryCoor& barycentricCoordinates);

    /**
     * @brief 重新计算某一点的参数
     *
     * @param faces
     * @param point2DMap
     * @param point
     *
     * @return
     */
    static std::optional<std::tuple<int, double, double, double>> UpdateParam(
        const std::vector<std::array<VertexHandle, 3>>& faces,
        const std::map<VertexHandle, Point2D>& point2DMap, const Point2D& point);

   private:
    double maxRingArea;
    double maxCurvature;
    MapMesh* baseLevelMesh = nullptr;
    std::vector<FaceHandle> originFaces;  // 原始网格中所有的面
    std::map<FaceHandle, std::vector<VertexHandle>> originFaceVertices;
    std::priority_queue<VertexHandle, std::vector<VertexHandle>, CompareFun> curvatureQueue;
};
}  // namespace Maps

#endif
