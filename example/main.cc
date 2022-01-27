#include "MapsMesh.h"
#include "igl/principal_curvature.h"
#include "igl/readPLY.h"
#include "igl/writeOBJ.h"

int main(void) {
    Maps::MapMesh mesh("model/Armadillo.obj");

    for (int i = 0; i < 1; i++) {
        mesh.Initialize();
        mesh.DownSampling();
        mesh.garbage_collection();
    }

    try {
        if (!OpenMesh::IO::write_mesh(mesh, "Armadillo_output.off")) {
            std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
            return 1;
        }
    } catch (std::exception& x) {
        std::cerr << x.what() << std::endl;
        return 1;
    }

    // Eigen::MatrixX3d V, PD1, PD2;
    // Eigen::MatrixX3i F;
    // Eigen::VectorXd PV1, PV2, PV;
    // igl::readPLY("model/Armadillo.ply", V, F);
    // igl::writeOBJ("model/Armdillo.obj", V, F);

    // igl::principal_curvature(V, F, PD1, PD2, PV1, PV2);

    return 0;
}
