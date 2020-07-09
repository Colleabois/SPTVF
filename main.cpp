#include "mainwindow.h"
#include <QApplication>
#include <iostream>
//#include <Eigen/Sparse>
//#include <Eigen/Eigenvalues>
//#include "meshVFE.h"
//#include "matrices.h"
#include "vectorfield.h"
#include "weakmesh.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    //construct_model();
    //vertexUpdateFile(":/output/model_2139-13032020.txt", 1);
    //test_eigens();
    //constructEigenSpace();
    return a.exec();
}
