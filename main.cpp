#include <igl/opengl/glfw/Viewer.h>
#include <igl/triangle/triangulate.h>
#include "MeshConnectivity.h"
#include "AmpSolver.h"
#include <igl/jet.h>

int exampleid = 0;

Eigen::Vector3d vectorField(Eigen::Vector3d pos)
{
    switch (exampleid)
    {
    case 0:
    {
        // whirlpool

        Eigen::Vector3d ans(-pos[1], pos[0], 0);
        ans /= pos.squaredNorm();
        return ans;
    }
    case 1:
    {
        // plane wave
        return Eigen::Vector3d(1, 1, 0);
    }
    default:
        // constant
        return Eigen::Vector3d(0, 0, 0);
    }
}

void createSquare(Eigen::MatrixXd& V, Eigen::MatrixXi& F, double triarea)
{
    int totpts = 4;
    Eigen::MatrixXd planeV(totpts, 2);
    Eigen::MatrixXi planeE(totpts, 2);
    planeV(0, 0) = -1.0;
    planeV(0, 1) = -1.0;
    planeV(1, 0) = 1.0;
    planeV(1, 1) = -1.0;
    planeV(2, 0) = 1.0;
    planeV(2, 1) = 1.0;
    planeV(3, 0) = -1.0;
    planeV(3, 1) = 1.0;
    planeE(0, 0) = 0;
    planeE(0, 1) = 1;
    planeE(1, 0) = 1;
    planeE(1, 1) = 2;
    planeE(2, 0) = 2;
    planeE(2, 1) = 3;
    planeE(3, 0) = 3;
    planeE(3, 1) = 0;

    Eigen::MatrixXi H(0, 2);

    Eigen::MatrixXd V3d;
    Eigen::MatrixXi F3d;
    std::stringstream ss;
    ss << "a" << triarea << "q";
    igl::triangle::triangulate(planeV, planeE, H, ss.str(), V3d, F3d);

    F = F3d;

    int nverts = V3d.rows();
    V.resize(nverts, 3);
    for (int i = 0; i < nverts; i++)
    {
        V(i, 0) = V3d(i, 0);
        V(i, 1) = V3d(i, 1);
        V(i, 2) = 0;
    }

}

int main(int argc, char *argv[])
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    createSquare(V, F, 0.01);

    MeshConnectivity mesh(F);
    int nedges = mesh.nEdges();
    Eigen::MatrixXd omegas(nedges, 2);
    for (int i = 0; i < nedges; i++)
    {
        int v0 = mesh.edgeVertex(i, 0);
        int v1 = mesh.edgeVertex(i, 1);
        Eigen::Vector3d end0 = V.row(v0).transpose();
        Eigen::Vector3d end1 = V.row(v1).transpose();
        Eigen::Vector3d pos = 0.5 * (end0 + end1);
        Eigen::Vector3d vec = vectorField(pos);
        double omega = vec.dot(end1 - end0);
        omegas(i, 0) = omega;
        omegas(i, 1) = omega;
    }
    Eigen::VectorXd amp;
    ampSolver(V, mesh, omegas, amp);

    //std::cout << amp << std::endl;
    int nverts = V.rows();
    Eigen::MatrixXd C(nverts, 3);
    igl::jet(amp, false, C);

    std::cout << "minimum amplitude: " << amp.minCoeff() << ", maximum: " << amp.maxCoeff() << std::endl;

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.data().set_colors(C);
  viewer.data().set_face_based(false);
  viewer.launch();
}
