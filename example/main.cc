#include <OpenMesh/Core/IO/MeshIO.hh>

#include "MapsMesh.h"
#include "igl/gaussian_curvature.h"
#include "igl/principal_curvature.h"
#include "igl/readOBJ.h"
#include "igl/writeOBJ.h"

int main(void) {
    Maps::MapMesh mesh;
    if (!OpenMesh::IO::read_mesh(mesh, "model/newmodel.obj")) {
        std::cout << "error load mesh from file";
        return -1;
    }

    mesh.Initialize();
    for (int i = 0; i < 3; i++) {
        mesh.DownSampling();
    }
    mesh.FaceSubDivision();
    mesh.FaceSubDivision();
    std::cout << "subdivided finished!" << std::endl;
    mesh.Remesh();

    mesh.garbage_collection();
    try {
        if (!OpenMesh::IO::write_mesh(mesh, "model_output.off")) {
            std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
            return 1;
        }
    } catch (std::exception& x) {
        std::cerr << x.what() << std::endl;
        return 1;
    }

    return 0;
}
