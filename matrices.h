#ifndef MATRICES_H
#define MATRICES_H
#include <meshVFE.h>
#include <Eigen/Sparse>

Eigen::Vector3d position(Vertex & myvertex);
double faceArea( Mesh & mymesh, face_index f_i);
face_index findNextAdjFace(const Mesh & mymesh, vertex_index myvertex, face_index ref_face);

Eigen::SparseMatrix<double> meshToMX(Mesh & mymesh);
Eigen::SparseMatrix<double> meshToMX_inv(Mesh & mymesh);
Eigen::SparseMatrix<double> meshToMS(Mesh & mymesh);
Eigen::SparseMatrix<double> meshToMS_inv(Mesh & mymesh);
Eigen::SparseMatrix<double> meshToMSstar(Mesh & mymesh);
Eigen::SparseMatrix<double> meshToMSstar_inv(Mesh & mymesh);
Eigen::SparseMatrix<double> meshToGS(Mesh & mymesh);
Eigen::SparseMatrix<double> meshToGSstar(Mesh & mymesh);
Eigen::SparseMatrix<double> meshToJ(Mesh & mymesh);

Eigen::SparseMatrix<double> meshToDiv_h(Mesh & mymesh);
Eigen::SparseMatrix<double> meshToCurl_h(Mesh & mymesh);
Eigen::SparseMatrix<double> meshToDiv_hs(Mesh & mymesh);
Eigen::SparseMatrix<double> meshToCurl_hs(Mesh & mymesh);

Eigen::SparseMatrix<double> meshToHodgeLaplace(Mesh & mymesh);
Eigen::SparseMatrix<double> meshToDirichlet(Mesh & mymesh);
Eigen::SparseMatrix<double> meshToBiharmonic(Mesh & mymesh);

Eigen::SparseMatrix<double> meshToDirichlet_S(Mesh & mymesh);
Eigen::SparseMatrix<double> meshToLaplaceBeltrami_S(Mesh & mymesh);
Eigen::SparseMatrix<double> meshToDirichlet_Ss(Mesh & mymesh);
Eigen::SparseMatrix<double> meshToLaplaceBeltrami_Ss(Mesh & mymesh);
#endif // MATRICES_H
