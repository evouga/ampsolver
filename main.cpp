#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include <igl/triangle/triangulate.h>
#include "MeshConnectivity.h"
#include "AmpSolver.h"
#include <igl/jet.h>
#include <igl/readOBJ.h>
#include <fstream>
#include <igl/cotmatrix_entries.h>
#include <igl/massmatrix.h>

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
    polyscope::init();

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
        Eigen::Vector3d vec1 = vectorField(end0);
        Eigen::Vector3d vec2 = vectorField(end1);
        double omega1 = vec1.dot(end1 - end0);
        double omega2 = vec2.dot(end0 - end1);        
        omegas(i, 0) = std::isnan(omega1) ? 1e6 : omega1;
        omegas(i, 1) = std::isnan(omega2) ? 1e6 : omega2;
    }

    Eigen::VectorXd amp;
    ampSolver(V, mesh, omegas, amp);

    std::map<std::pair<int, int>, int > polyordering;

    int idx = 0;
    for (int i = 0; i < F.rows(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            int v0 = F(i, j);
            int v1 = F(i, (j + 1) % 3);
            auto it = polyordering.find({ v0,v1 });
            if (it != polyordering.end())
                continue;
            polyordering[{v0, v1}] = idx;
            polyordering[{v1, v0}] = idx;
            idx++;
        }
    }

    std::vector<int> perm(nedges);
    for (int i = 0; i < nedges; i++)
    {
        //perm[i] = polyordering[{mesh.edgeVertex(i, 0), mesh.edgeVertex(i, 1)}];
        perm[polyordering[{mesh.edgeVertex(i, 0), mesh.edgeVertex(i, 1)}]] = i;
    }




    auto* pmesh = polyscope::registerSurfaceMesh("Mesh", V, F);
    pmesh->setEdgePermutation(perm);
    pmesh->addEdgeScalarQuantity("Omega 1", omegas.col(0));
    pmesh->addEdgeScalarQuantity("Omega 2", omegas.col(1));
    pmesh->addVertexScalarQuantity("Amp", amp);

    // visualize!
    polyscope::show();
}
