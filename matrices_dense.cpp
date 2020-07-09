#include "meshVFE.h"
#include <tsvet.h>

#include <QFile>
#include <QString>
#include <QDebug>
#include <QTextStream>
#include <cmath>

#include <queue>
#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Geometry> // for cross product


face_index findNextAdjFace(const Mesh & mymesh, vertex_index myvertex, face_index ref_face)
{
    face_index target;
    for(int k=0; k<3; k++)
    {
                if(mymesh.fs()[ref_face].vis()[k] == myvertex)
                {
                    target = mymesh.fs()[ref_face].fis()[(k+1)%3];
                }
    }
    return target;
}

Eigen::Vector3d position(Vertex & myvertex)
{   // return the position of the vertex
    return Eigen::Vector3d(myvertex.x(), myvertex.y(), myvertex.z());
}

double faceArea( Mesh & mymesh, face_index f_i)
{
    Eigen::Vector3d a = position(mymesh.vs()[ mymesh.fs()[f_i].vis()[0]]);
    Eigen::Vector3d b = position(mymesh.vs()[ mymesh.fs()[f_i].vis()[1]]);
    Eigen::Vector3d c = position(mymesh.vs()[ mymesh.fs()[f_i].vis()[2]]);
    return (b-a).cross(c-a).norm();
}


//----------------------------
// Basics Matrices
//----------------------------

typedef Eigen::Triplet<double> T;

Eigen::SparseMatrix<double> expandMatrix(Eigen::MatrixXd & mymatrix)
{   // expand every coefficient of mymatrix to a 2x2 block
    long n_rows = mymatrix.rows();
    long n_cols = mymatrix.cols();
    Eigen::MatrixXd new_M = Eigen::MatrixXd::Zero(n_rows*2,n_cols*2);
    for(long i=0; i<n_rows; i++)
    {
        for(long j=0; j<n_cols; j++)
        {
            new_M(2*i,2*j) = mymatrix(i,j);
            new_M(2*i+1,2*j+1) = mymatrix(i,j);
        }
    }
    return new_M;
}

// L2 scalar product in chi
Eigen::SparseMatrix<double> meshToMX(Mesh & mymesh)
{
    int n_face = mymesh.fs().size();
    QVector<T> tripletList;
    tripletList.reserve(n_face);
    for(int i=0; i<n_face; i++)
    {
      tripletList.push_back(T(i,i,faceArea(mymesh, i)));
    }
    Eigen::SparseMatrix<double> new_M(n_face,n_face);
    new_M.setFromTriplets(tripletList.begin(), tripletList.end());

/*
    int n_face = mymesh.fs().size();
    Eigen::MatrixXd new_M = Eigen::MatrixXd::Zero(n_face,n_face);
    for(int i=0; i<n_face; i++)
        {
        new_M(i,i) = faceArea(mymesh, i);
    }
*/
    return new_M;
}
// double area(Face & myface);


// L2 scalar product in S (vertex)
Eigen::SparseMatrix<double> meshToMS(Mesh & mymesh)
{
    int n_vertex = mymesh.vs().size();
    Eigen::MatrixXd new_M = Eigen::MatrixXd::Zero(n_vertex,n_vertex);
    for(int i=0; i<n_vertex; i++)
    {
        face_index start_face = mymesh.vs()[i].fi();
        new_M(start_face, start_face) = faceArea(mymesh,start_face)/3;
        face_index next_face = findNextAdjFace(mymesh, i, start_face);
        while(next_face!=start_face)
            {
            new_M(next_face,next_face) = faceArea(mymesh,next_face)/3;
            next_face = findNextAdjFace(mymesh, i,next_face);
            }
    }
    return new_M;
}

// L2 scalar product in S* (edges)
Eigen::SparseMatrix<double> meshToMSstar(Mesh & mymesh)
{
    int n_edge = mymesh.es().size();
    Eigen::MatrixXd new_M = Eigen::MatrixXd::Zero(n_edge,n_edge);
    for(int i=0; i<n_edge; i++){
        face_index j1 = mymesh.es()[i].fi1();
        face_index j2 = mymesh.es()[i].fi2();
        new_M(i,i) = faceArea(mymesh, j1)/3 + faceArea(mymesh, j2)/3;
    }
    return new_M;
}


// gradient in S (vertex)
Eigen::SparseMatrix<double> meshToGS(Mesh & mymesh)
{
    unsigned int n_face = mymesh.fs().size();
    unsigned int n_vertex = mymesh.vs().size();
    Eigen::MatrixXd new_M = Eigen::MatrixXd::Zero(n_vertex,n_face);
    for(face_index i=0; i<n_face; i++)
    {
        Eigen::Vector3d x0 = position(mymesh.vs()[mymesh.fs()[i].vis()[0]]);
        Eigen::Vector3d x1 = position(mymesh.vs()[mymesh.fs()[i].vis()[1]]);
        Eigen::Vector3d x2 = position(mymesh.vs()[mymesh.fs()[i].vis()[2]]);
        double area = faceArea(mymesh, i);
        Eigen::Vector3d normal = ((x1-x0).cross(x2-x0)).normalized(); //don't know qhy normalize() doesn't work here.
        Eigen::Vector3d xixkp = Eigen::AngleAxisd(M_PI/2., normal).toRotationMatrix()*(x0-x2);
        Eigen::Vector3d xjxi = x1-x0;
        double nxjxi = xjxi.norm();
        Eigen::Vector3d xjxip = Eigen::AngleAxisd(M_PI/2., normal).toRotationMatrix()*(xjxi);
        double a1 = xixkp.dot(xjxi)/2/area/nxjxi;
        double a2 = xixkp.dot(xjxip)/2/area/nxjxi;
        double a3 = nxjxi/2/area;
        new_M(mymesh.fs()[i].vis()[0],i) = -a1-a2-a3;
        new_M(mymesh.fs()[i].vis()[1],i) = a1+a2;
        new_M(mymesh.fs()[i].vis()[2],i) = a3;
    }
    return new_M;
}

// gradient in S* (vertex)
Eigen::SparseMatrix<double> meshToGSstar(Mesh & mymesh)
{
    int n_face = mymesh.fs().size();
    int n_edge = mymesh.es().size();
    Eigen::MatrixXd new_M = Eigen::MatrixXd::Zero(n_edge,n_face);
    for(face_index i=0; i<n_face; i++)
    {
        Eigen::Vector3d xv0 = position(mymesh.vs()[mymesh.fs()[i].vis()[0]]);
        Eigen::Vector3d xv1 = position(mymesh.vs()[mymesh.fs()[i].vis()[1]]);
        Eigen::Vector3d xv2 = position(mymesh.vs()[mymesh.fs()[i].vis()[2]]);
        Eigen::Vector3d x0 = (xv1+xv2)/2;
        Eigen::Vector3d x1 = (xv2+xv0)/2;
        Eigen::Vector3d x2 = (xv0+xv1)/2;
        double area = faceArea(mymesh, i)/4.;
        Eigen::Vector3d normal = ((x1-x0).cross(x2-x0)).normalized(); //don't know qhy normalize() doesn't work here.
        Eigen::Vector3d xixkp = Eigen::AngleAxisd(M_PI/2., normal).toRotationMatrix()*(x0-x2);
        Eigen::Vector3d xjxi = x1-x0;
        double nxjxi = xjxi.norm();
        Eigen::Vector3d xjxip = Eigen::AngleAxisd(M_PI/2., normal).toRotationMatrix()*(xjxi);
        double a1 = xixkp.dot(xjxi)/2/area/nxjxi;
        double a2 = xixkp.dot(xjxip)/2/area/nxjxi;
        double a3 = nxjxi/2/area;
        new_M(mymesh.fs()[i].vis()[0],i) = -a1-a2-a3;
        new_M(mymesh.fs()[i].vis()[1],i) = a1+a2;
        new_M(mymesh.fs()[i].vis()[2],i) = a3;
    }
    return new_M;
}



Eigen::SparseMatrix<double> meshToJ(Mesh & mymesh)
{   // rotation matrix J
    int n_face = mymesh.fs().size();
    Eigen::MatrixXd new_M = Eigen::MatrixXd::Zero(2*n_face,2*n_face);
    for(int i=0; i<n_face; i++)
    {
         new_M(2*i+1,2*i) = -1;
         new_M(2*i,2*i+1) =  1;
    }
    return new_M;
}

Eigen::SparseMatrix<double> meshToDiv_h(Mesh & mymesh)
{
    Eigen::MatrixXd MS_inv = meshToMS(mymesh).inverse();
    Eigen::MatrixXd GS_T = meshToGS(mymesh).transpose();
    Eigen::MatrixXd MX = meshToMX(mymesh);
    Eigen::MatrixXd div = MS_inv*GS_T*MX;
    return expandMatrix(div);
}

Eigen::SparseMatrix<double> meshToCurl_h(Mesh & mymesh)
{
    Eigen::MatrixXd MS_inv0 = meshToMS(mymesh).inverse();
    Eigen::MatrixXd MS_inv = expandMatrix(MS_inv0);
    Eigen::MatrixXd GS_T0 = meshToGS(mymesh).transpose();
    Eigen::MatrixXd GS_T = expandMatrix(GS_T0);
    Eigen::MatrixXd J = meshToJ(mymesh);
    Eigen::MatrixXd MX0 = meshToMX(mymesh);
    Eigen::MatrixXd MX = expandMatrix(MX0);
    return -MS_inv*GS_T*J*MX;
}

Eigen::SparseMatrix<double> meshToDiv_hs(Mesh & mymesh)
{
    Eigen::MatrixXd MS_inv = meshToMSstar(mymesh).inverse();
    Eigen::MatrixXd GS_T = meshToGSstar(mymesh).transpose();
    Eigen::MatrixXd MX = meshToMX(mymesh);
    Eigen::MatrixXd div = MS_inv*GS_T*MX;
    return expandMatrix(div);
}

Eigen::SparseMatrix<double> meshToCurl_hs(Mesh & mymesh)
{
    Eigen::MatrixXd MS_inv0 = meshToMSstar(mymesh).inverse();
    Eigen::MatrixXd MS_inv = expandMatrix(MS_inv0);
    Eigen::MatrixXd GS_T0 = meshToGSstar(mymesh).transpose();
    Eigen::MatrixXd GS_T = expandMatrix(GS_T0);
    Eigen::MatrixXd J = meshToJ(mymesh);
    Eigen::MatrixXd MX0 = meshToMX(mymesh);
    Eigen::MatrixXd MX = expandMatrix(MX0);
    return -MS_inv*GS_T*J*MX;
}

Eigen::SparseMatrix<double> meshToHodgeLaplace(Mesh & mymesh)
{   // HodgeLaplace operator
    Eigen::MatrixXd MS_inv0 = meshToMS(mymesh).inverse();
    Eigen::MatrixXd MS_inv = expandMatrix(MS_inv0);

    Eigen::MatrixXd MSs_inv0 = meshToMSstar(mymesh).inverse();
    Eigen::MatrixXd MSs_inv = expandMatrix(MSs_inv0);

    Eigen::MatrixXd GS0 = meshToGS(mymesh);
    Eigen::MatrixXd GS = expandMatrix(GS0);
    Eigen::MatrixXd GS_T0 = GS.transpose();
    Eigen::MatrixXd GS_T = expandMatrix(GS_T0);

    Eigen::MatrixXd GS0s = meshToGSstar(mymesh);
    Eigen::MatrixXd GSs = expandMatrix(GS0);
    Eigen::MatrixXd GSs_T0 = GS.transpose();
    Eigen::MatrixXd GSs_T = expandMatrix(GS_T0);

    Eigen::MatrixXd J = meshToJ(mymesh);

    Eigen::MatrixXd MX0 = meshToMX(mymesh);
    Eigen::MatrixXd MX = expandMatrix(MX0);

    return (GS*MS_inv*GS_T - J*GSs*MSs_inv*GSs_T*J)*MX;
}

Eigen::SparseMatrix<double> meshToDirichlet(Mesh & mymesh)
{   // Dirichlet Energy of vector fields
    Eigen::MatrixXd HL = meshToHodgeLaplace(mymesh);
    Eigen::MatrixXd MX0 = meshToMX(mymesh);
    Eigen::MatrixXd MX = expandMatrix(MX0);
    return MX*HL;
}


Eigen::SparseMatrix<double> meshToBiharmonic(Mesh & mymesh)
{   // Biharmonic Energy
    Eigen::MatrixXd HL = meshToHodgeLaplace(mymesh);
    Eigen::MatrixXd HL_T = HL.transpose();
    Eigen::MatrixXd MX0 = meshToMX(mymesh);
    Eigen::MatrixXd MX = expandMatrix(MX0);
    return HL_T*MX*HL;
}


Eigen::SparseMatrix<double> meshToDirichlet_S(Mesh & mymesh)
{   // Dirichlet Energy of functions on S
    Eigen::MatrixXd GS0 = meshToGS(mymesh);
    Eigen::MatrixXd GS = expandMatrix(GS0);
    Eigen::MatrixXd GS_T0 = GS.transpose();
    Eigen::MatrixXd GS_T = expandMatrix(GS_T0);

    Eigen::MatrixXd MX0 = meshToMX(mymesh);
    Eigen::MatrixXd MX = expandMatrix(MX0);
    return GS_T*MX*GS;
}

Eigen::MatrixXd meshToLaplaceBeltrami_S(Mesh & mymesh)
{   // Laplace Beltrami operator on S
    Eigen::MatrixXd S = meshToDirichlet_S(mymesh);
    Eigen::MatrixXd MS_inv0 = meshToMS(mymesh).inverse();
    Eigen::MatrixXd MS_inv = expandMatrix(MS_inv0);
    return MS_inv*S;
}

Eigen::MatrixXd meshToDirichlet_Ss(Mesh & mymesh)
{   // Dirichlet Energy of functions on S
    Eigen::MatrixXd GS0 = meshToGSstar(mymesh);
    Eigen::MatrixXd GS = expandMatrix(GS0);
    Eigen::MatrixXd GS_T0 = GS.transpose();
    Eigen::MatrixXd GS_T = expandMatrix(GS_T0);

    Eigen::MatrixXd MX0 = meshToMX(mymesh);
    Eigen::MatrixXd MX = expandMatrix(MX0);
    return GS_T*MX*GS;
}

Eigen::MatrixXd meshToLaplaceBeltrami_Ss(Mesh & mymesh)
{   // Laplace Beltrami operator on S*
    Eigen::MatrixXd S = meshToDirichlet_Ss(mymesh);
    Eigen::MatrixXd MS_inv0 = meshToMSstar(mymesh).inverse();
    Eigen::MatrixXd MS_inv = expandMatrix(MS_inv0);
    return MS_inv*S;
}

//----------------------------
// Computation of Eigenfields
//----------------------------

void eigenVectorField(Mesh mymesh)
{
    Eigen::MatrixXd LB_S = meshToLaplaceBeltrami_S(mymesh);
    Eigen::MatrixXd LB_Ss = meshToLaplaceBeltrami_Ss(mymesh);

    Eigen::EigenSolver<Eigen::MatrixXd> es_S(LB_S);
    std::cout << "The eigenvalues of m are:" << std::endl << es_S.eigenvalues() << std::endl;
    std::cout << "The matrix of eigenvectors, V, is:" << std::endl << es_S.eigenvectors() << std::endl << std::endl;

    Eigen::EigenSolver<Eigen::MatrixXd> es_Ss(LB_Ss);
    std::cout << "The eigenvalues of m are:" << std::endl << es_Ss.eigenvalues() << std::endl;
    std::cout << "The matrix of eigenvectors, V, is:" << std::endl << es_Ss.eigenvectors() << std::endl << std::endl;

}



