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
//#include <Eigen/Dense>
#include <Eigen/Geometry> // for cross product

face_index findNextAdjFace(const Mesh & mymesh, vertex_index myvertex, face_index ref_face)
{
    face_index target = mymesh.fs()[ref_face].fis()[1];
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


Eigen::SparseMatrix<double> meshToMX(Mesh & mymesh)
{   // L2 scalar product in chi
    int n_face = mymesh.fs().size();
    QVector<T> tripletList;
    tripletList.reserve(2*n_face);
    for(int i=0; i<n_face; i++)
    {
      tripletList.push_back(T(2*i,2*i,faceArea(mymesh, i)));
      tripletList.push_back(T(2*i+1,2*i+1,faceArea(mymesh, i)));
    }
    Eigen::SparseMatrix<double> new_M(2*n_face,2*n_face);
    new_M.setFromTriplets(tripletList.begin(), tripletList.end());
    return new_M;
}
Eigen::SparseMatrix<double> meshToMX_inv(Mesh & mymesh)
{   // Inverse of L2 scalar product in chi
    int n_face = mymesh.fs().size();
    QVector<T> tripletList;
    tripletList.reserve(2*n_face);
    for(int i=0; i<n_face; i++)
    {
      tripletList.push_back(T(2*i,2*i,1./faceArea(mymesh, i)));
      tripletList.push_back(T(2*i+1,2*i+1,1./faceArea(mymesh, i)));
    }
    Eigen::SparseMatrix<double> new_M(2*n_face,2*n_face);
    new_M.setFromTriplets(tripletList.begin(), tripletList.end());
    return new_M;
}

Eigen::SparseMatrix<double> meshToMS(Mesh & mymesh)
{   // L2 scalar product in S (vertex)
    int n_vertex = mymesh.vs().size();
    QVector<T> tripletList;
    tripletList.reserve(n_vertex);
    for(int i=0; i<n_vertex; i++)
    {
        double area = 0;
        face_index start_face = mymesh.vs()[i].fi();
        area += faceArea(mymesh, start_face)/double(3);
        face_index next_face = findNextAdjFace(mymesh, i, start_face);
        while(next_face!=start_face)
            {
            area += faceArea(mymesh, next_face)/double(3);
            next_face = findNextAdjFace(mymesh, i,next_face);
            }
        tripletList.push_back(T(i,i,area));
    }
    Eigen::SparseMatrix<double> new_M(n_vertex,n_vertex);
    new_M.setFromTriplets(tripletList.begin(), tripletList.end());
    return new_M;
}
Eigen::SparseMatrix<double> meshToMS_inv(Mesh & mymesh)
{   // Inverse of L2 scalar product in S (vertex)
    int n_vertex = mymesh.vs().size();
    QVector<T> tripletList;
    tripletList.reserve(n_vertex);
    for(int i=0; i<n_vertex; i++)
    {
        double area = 0;
        face_index start_face = mymesh.vs()[i].fi();
        area+= 3/faceArea(mymesh, start_face);
        face_index next_face = findNextAdjFace(mymesh, i, start_face);
        while(next_face!=start_face)
            {
            area+= 3./faceArea(mymesh, next_face);
            next_face = findNextAdjFace(mymesh, i,next_face);
            }
        tripletList.push_back(T(i,i,area));
    }
    Eigen::SparseMatrix<double> new_M(n_vertex,n_vertex);
    new_M.setFromTriplets(tripletList.begin(), tripletList.end());
    return new_M;
}

Eigen::SparseMatrix<double> meshToMSstar(Mesh & mymesh)
{   // L2 scalar product in S* (edges)
    int n_edge = mymesh.es().size();
    QVector<T> tripletList;
    tripletList.reserve(n_edge);
    for(int i=0; i<n_edge; i++){
        face_index j1 = mymesh.es()[i].fi1();
        face_index j2 = mymesh.es()[i].fi2();
        tripletList.push_back(T(i,i,faceArea(mymesh, j1)/3. + faceArea(mymesh, j2)/3.));
    }
    Eigen::SparseMatrix<double> new_M(n_edge,n_edge);
    new_M.setFromTriplets(tripletList.begin(), tripletList.end());
    return new_M;
}
Eigen::SparseMatrix<double> meshToMSstar_inv(Mesh & mymesh)
{   // Inverse L2 scalar product in S* (edges)
    int n_edge = mymesh.es().size();
    QVector<T> tripletList;
    tripletList.reserve(n_edge);
    for(int i=0; i<n_edge; i++){
        face_index j1 = mymesh.es()[i].fi1();
        face_index j2 = mymesh.es()[i].fi2();
        tripletList.push_back(T(i,i,3./(faceArea(mymesh, j1) + faceArea(mymesh, j2))));
    }
    Eigen::SparseMatrix<double> new_M(n_edge,n_edge);
    new_M.setFromTriplets(tripletList.begin(), tripletList.end());
    return new_M;
}

Eigen::SparseMatrix<double> meshToGS(Mesh & mymesh)
{   // gradient in S (vertex)
    unsigned int n_face = mymesh.fs().size();
    unsigned int n_vertex = mymesh.vs().size();
    QVector<T> tripletList;
    tripletList.reserve(n_face*6);
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
        double a3 = nxjxi/2./area;
        tripletList.push_back(T(2*i,mymesh.fs()[i].vis()[0],-a1));
        tripletList.push_back(T(2*i+1,mymesh.fs()[i].vis()[0],-a2-a3));
        tripletList.push_back(T(2*i,mymesh.fs()[i].vis()[1], a1));
        tripletList.push_back(T(2*i+1,mymesh.fs()[i].vis()[1], a2));
        tripletList.push_back(T(2*i,mymesh.fs()[i].vis()[2], 0));
        tripletList.push_back(T(2*i+1,mymesh.fs()[i].vis()[2], a3));
    }
    Eigen::SparseMatrix<double> new_M(2*n_face,n_vertex);
    new_M.setFromTriplets(tripletList.begin(), tripletList.end());
    return new_M;
}
QVector<edge_index> findEdgeIndex(Mesh & mymesh, face_index fi)
{   // find three edges that surround the face fi
    edge_index n_edge = mymesh.es().size();
    QMap<face_index,edge_index> list_of_edges;
    for(edge_index ei=0; ei<n_edge; ei++)
    {
        Edge myedge = mymesh.es()[ei];
        if(myedge.fi1()==fi) list_of_edges[myedge.fi2()] = ei;
        if(myedge.fi2()==fi) list_of_edges[myedge.fi1()] = ei;
    }
    QVector<face_index> fiTab = mymesh.fs()[fi].fis();
    QVector<edge_index> eiTab;
    for(int i=0; i<3; i++)
    {
        eiTab.push_back(list_of_edges[fiTab[i]]);
    }
    return eiTab;
}
Eigen::SparseMatrix<double> meshToGSstar(Mesh & mymesh)
{   // gradient in S* (vertex)
    int n_face = mymesh.fs().size();
    int n_edge = mymesh.es().size();
    QVector<T> tripletList;
    tripletList.reserve(n_face*6);
    for(face_index i=0; i<n_face; i++)
    {
        Eigen::Vector3d x0 = position(mymesh.vs()[mymesh.fs()[i].vis()[0]]);
        Eigen::Vector3d x1 = position(mymesh.vs()[mymesh.fs()[i].vis()[1]]);
        Eigen::Vector3d x2 = position(mymesh.vs()[mymesh.fs()[i].vis()[2]]);

        double area = faceArea(mymesh, i);
        Eigen::Vector3d normal = ((x1-x0).cross(x2-x0)).normalized();
        Eigen::Vector3d xixkp = Eigen::AngleAxisd(M_PI/2., normal).toRotationMatrix()*(x0-x2);
        Eigen::Vector3d xjxi = x1-x0;
        double nxjxi = xjxi.norm();
        Eigen::Vector3d xjxip = Eigen::AngleAxisd(M_PI/2., normal).toRotationMatrix()*(xjxi);
        double a1 = -xixkp.dot(xjxi)/area/nxjxi;
        double a2 = -xixkp.dot(xjxip)/area/nxjxi;
        double a3 = -nxjxi/area;
        QVector<edge_index> eis = findEdgeIndex(mymesh, i);
        //qDebug() <<"face" << i << " " << eis[0]<<" " << eis[1]<<" " << eis[2];
        tripletList.push_back(T(2*i, eis[0],-a1));
        tripletList.push_back(T(2*i+1, eis[0],-a2-a3));

        tripletList.push_back(T(2*i, eis[1], a1));
        tripletList.push_back(T(2*i+1, eis[1], a2));
        tripletList.push_back(T(2*i, eis[2], 0));
        tripletList.push_back(T(2*i+1, eis[2], a3));
    }
    Eigen::SparseMatrix<double> new_M(2*n_face,n_edge);
    new_M.setFromTriplets(tripletList.begin(), tripletList.end());
    return new_M;
}

Eigen::SparseMatrix<double> meshToJ(Mesh & mymesh)
{   // rotation matrix J
    int n_face = mymesh.fs().size();
    QVector<T> tripletList;
    tripletList.reserve(n_face*2);
    for(int i=0; i<n_face; i++)
    {
         tripletList.push_back(T(2*i+1,2*i,-1));
         tripletList.push_back(T(2*i,2*i+1, 1));
    }
    Eigen::SparseMatrix<double> new_M(n_face*2,n_face*2);
    new_M.setFromTriplets(tripletList.begin(), tripletList.end());
    return new_M;
}
Eigen::SparseMatrix<double> meshToDiv_h(Mesh & mymesh)
{
    Eigen::SparseMatrix<double> MS_inv = meshToMS_inv(mymesh);
    Eigen::SparseMatrix<double> GS_T = meshToGS(mymesh).transpose();
    Eigen::SparseMatrix<double> MX = meshToMX(mymesh);
    Eigen::SparseMatrix<double> div = MS_inv*GS_T*MX;
    return div;
}
Eigen::SparseMatrix<double> meshToCurl_h(Mesh & mymesh)
{
    Eigen::SparseMatrix<double> MS_inv = meshToMS_inv(mymesh);
    Eigen::SparseMatrix<double> GS_T = meshToGS(mymesh).transpose();
    Eigen::SparseMatrix<double> J = meshToJ(mymesh);
    Eigen::SparseMatrix<double> MX = meshToMX(mymesh);
    return -MS_inv*GS_T*J*MX;
}
Eigen::SparseMatrix<double> meshToDiv_hs(Mesh & mymesh)
{
    Eigen::SparseMatrix<double> MS_inv = meshToMSstar_inv(mymesh);
    Eigen::SparseMatrix<double> GS_T = meshToGSstar(mymesh).transpose();
    Eigen::SparseMatrix<double> MX = meshToMX(mymesh);
    Eigen::SparseMatrix<double> div = MS_inv*GS_T*MX;
    return div;
}
Eigen::SparseMatrix<double> meshToCurl_hs(Mesh & mymesh)
{
    Eigen::SparseMatrix<double> MS_inv = meshToMSstar_inv(mymesh);
    Eigen::SparseMatrix<double> GS_T = meshToGSstar(mymesh).transpose();
    Eigen::SparseMatrix<double> J = meshToJ(mymesh);
    Eigen::SparseMatrix<double> MX = meshToMX(mymesh);
    return -MS_inv*GS_T*J*MX;
}

Eigen::SparseMatrix<double> meshToHodgeLaplace(Mesh & mymesh)
{   // HodgeLaplace operator
    Eigen::SparseMatrix<double> MX = meshToMX(mymesh);
    Eigen::SparseMatrix<double> MS_inv = meshToMS_inv(mymesh);
    Eigen::SparseMatrix<double> MSs_inv = meshToMSstar_inv(mymesh);

    Eigen::SparseMatrix<double> GS = meshToGS(mymesh);
    Eigen::SparseMatrix<double> GS_T = GS.transpose();
    Eigen::SparseMatrix<double> GSs = meshToGSstar(mymesh);
    Eigen::SparseMatrix<double> GSs_T = GSs.transpose();
    Eigen::SparseMatrix<double> J = meshToJ(mymesh);

    return (GS*MS_inv*GS_T - J*GSs*MSs_inv*GSs_T*J)*MX;
}
Eigen::SparseMatrix<double> meshToDirichlet(Mesh & mymesh)
{   // Dirichlet Energy of vector fields
    Eigen::SparseMatrix<double> HL = meshToHodgeLaplace(mymesh);
    Eigen::SparseMatrix<double> MX = meshToMX(mymesh);
    return MX*HL;
}
Eigen::SparseMatrix<double> meshToBiharmonic(Mesh & mymesh)
{   // Biharmonic Energy
    Eigen::SparseMatrix<double> Dir = meshToDirichlet(mymesh);
    Eigen::SparseMatrix<double> MX_inv = meshToMX_inv(mymesh);
    return Dir*MX_inv*Dir;
}



Eigen::SparseMatrix<double> meshToDirichlet_S(Mesh & mymesh)
{   // Dirichlet Energy of functions on S
    Eigen::SparseMatrix<double> GS = meshToGS(mymesh);
    Eigen::SparseMatrix<double> GS_T = GS.transpose();
    Eigen::SparseMatrix<double> MX = meshToMX(mymesh);
    return GS_T*MX*GS;
}
Eigen::SparseMatrix<double> meshToLaplaceBeltrami_S(Mesh & mymesh)
{   // Laplace Beltrami operator on S
    Eigen::SparseMatrix<double> S = meshToDirichlet_S(mymesh);
    Eigen::SparseMatrix<double> MS_inv = meshToMS_inv(mymesh);
    return MS_inv*S;
}
Eigen::SparseMatrix<double> meshToDirichlet_Ss(Mesh & mymesh)
{   // Dirichlet Energy of functions on S
    Eigen::SparseMatrix<double> GS = meshToGSstar(mymesh);
    Eigen::SparseMatrix<double> GS_T = GS.transpose();
    Eigen::SparseMatrix<double> MX = meshToMX(mymesh);
    return GS_T*MX*GS;
}
Eigen::SparseMatrix<double> meshToLaplaceBeltrami_Ss(Mesh & mymesh)
{   // Laplace Beltrami operator on S*
    Eigen::SparseMatrix<double> S = meshToDirichlet_Ss(mymesh);
    Eigen::SparseMatrix<double> MS_inv = meshToMSstar_inv(mymesh);
    return MS_inv*S;
}



