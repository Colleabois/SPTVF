
#include "gldisplaywidget.h"
#ifdef __APPLE__
    #include <glu.h>
#else
    #include <GL/glu.h>
#endif

//#include <glm-0.9.9.6/glm/glm.hpp> //complete glu

#include "QDebug"
#include "meshVFE.h"
#include "matrices.h"
#include "vectorfield.h"
#include "weakmesh.h"


GLDisplayWidget::GLDisplayWidget(QWidget *parent) : QGLWidget(parent),
    _choice (1),_cameraX(0),_cameraY(0),_cameraZ(10), _theta(0), _distance(45)
{
    // Update the scene
    connect( &_timer, SIGNAL(timeout()), this, SLOT(updateGL()));
    _timer.start(16);
}

void GLDisplayWidget::initializeGL()
{
    // background color
    glClearColor(0.2f, 0.2f, 0.2f, 1.0f);

    // Shader
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHT0);
    glEnable(GL_COLOR_MATERIAL);
}

void GLDisplayWidget::paintGL(){

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Center the camera
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    //gluLookAt(0,0,200,  0,0,0,   0,1,0); //eye, center, up
    _cameraX = _distance*cos(_theta/10*M_PI);
    _cameraY = _distance*sin(_theta/10*M_PI);
    gluLookAt(_cameraX,_cameraY,_cameraZ,  0., 0.,_cameraZ, 0., 0., 1.); //eye, center, up // also deprecated

    // Translation
    glTranslated(_X, _Y, _Z);

    // Rotation
    glRotatef(_angle, 1.0f, 1.0f, 0.0f);

    //QString filename = ":/output/model_2139-13032020.txt";
    //QString updated =":/output/Vertex_Updated_2154-13032020.txt";
    //const QMap<QString, Eigen::MatrixXd> map = read_model(filename);
    //Mesh modmesh = read_mesh_from_file(updated);
    //VectorField myvf = solveSytem_read(filename);

    //QString filename = ":/output/Eigen_Space_1411-26052020.txt";
    QString filename = ":/output/Eigen_Space_2146-26052020.txt";
    Mesh mymesh = read_mesh_from_file(filename);
    const QMap<QString, Eigen::MatrixXd> map = read_eigen(filename);
    VectorField myvf1 = 25*showEigenField(filename, 1, 61);
    VectorFieldOnMesh myvfom1 = VectorFieldOnMesh(mymesh, myvf1);
    VectorField myvf2 = 25*showEigenField(filename, 2, 61);
    VectorFieldOnMesh myvfom2 = VectorFieldOnMesh(mymesh, myvf2);

    // for tests' use
    if(_choice% 2==0)
    {
        showVfOnMesh2(map, myvfom1, myvfom2);
        //showVfOnMesh(map, myvfom);
    }
    // possible display
    if(_choice% 3==0)
    { //1
        //showVfWithTriangles(map, myvfom);

    }
    if(_choice% 5==0)
    { //2
        //plot_curves_as_conditions(map, mymesh);
        //plotMeshNaive(modmesh, grassColor()[3]);
        //Eigen::MatrixXd U = map.value("U");
        //VectorField my_eigen_vector = U.col(5);
        //VectorFieldOnMesh myvfom2 = VectorFieldOnMesh(mymesh, my_eigen_vector);
        //showVfOnMesh(map, myvfom2);
    }

}

void GLDisplayWidget::resizeGL(int width, int height){

    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0f, GLdouble(width)/GLdouble(height), 0.1f,300.0f);//fovy,aspect,zN,zF
// fovy: Specifies the field of view angle, in degrees, in the y direction.
// aspect: Specifies the aspect ratio that determines the field of view in the x direction.
// zNear: Specifies the distance from the viewer to the near clipping plane (always positive).
// zFar: Specifies the distance from the viewer to the far clipping plane (always positive).
/* gluPerspective is deprecated */
    updateGL();
}


// - - - - - - - - - - - - Mouse Management  - - - - - - - - - - - - - - - -
// When you click, the position of your mouse is saved
void GLDisplayWidget::mousePressEvent(QMouseEvent *event)
{
    if( event != nullptr )
        _lastPosMouse = event->pos();
}

// Mouse movement management
void GLDisplayWidget::mouseMoveEvent(QMouseEvent *event)
{
    int dx = event->x() - _lastPosMouse.x();
    // int dy = event->y() - lastPosMouse.y();

    if( event != nullptr )
    {
        _angle += dx;
        _lastPosMouse = event->pos();

        updateGL();
    }
}

// Mouse Management for the zoom
void GLDisplayWidget::wheelEvent(QWheelEvent *event) {
    QPoint numDegrees = event->angleDelta();
    double stepZoom = 0.25;
    if (!numDegrees.isNull())
    {
      _Z = (numDegrees.x() > 0 || numDegrees.y() > 0) ? _Z + stepZoom : _Z - stepZoom;
    }
}

// all the followings are added at TP3

void GLDisplayWidget::selectPositionAndColorOn3DModel() {
    int xMouse = _lastPosMouse.x();
    int yMouse = height() - 1 - _lastPosMouse.y();

    float OnObjectZ;

    unsigned char color[4];
    color[3] = 0;
    glReadPixels(xMouse, yMouse, 1, 1, GL_RGB,GL_UNSIGNED_BYTE, &color[0]);

    if(color[0] != 0 || color[1] != 0 || color[2] != 0)
    {
        glReadPixels(xMouse,yMouse,1,1,GL_DEPTH_COMPONENT,GL_FLOAT,&OnObjectZ);
        QVector3D posClickSelected = unproject(xMouse,yMouse,OnObjectZ); // Pos of click
        qDebug() << posClickSelected;
        qDebug() << color[0]<< ", " << color[1] << ", "<< color[2]; // Color of face
    }
}

// UnProjection: **Model** => GL_MODELVIEW_MATRIX => **View** => GL_PROJECTION_MATRIX ==> **Projection**
// Return the position on the face

QVector3D GLDisplayWidget::unproject(int X, int Y, float Depth)
{
  std::vector<GLdouble> ModelView(16),Projection(16);
  glGetDoublev(GL_MODELVIEW_MATRIX,&ModelView[0]);
  glGetDoublev(GL_PROJECTION_MATRIX,&Projection[0]);
  std::vector<GLint> Viewport(4);
  glGetIntegerv(GL_VIEWPORT,&Viewport[0]);
  GLdouble x,y,z;
  gluUnProject(X,Y,Depth,&ModelView[0],&Projection[0],&Viewport[0],&x,&y,&z);
// gluUnProject is deprecated use "glm::unProject"
  return QVector3D(x,y,z);
}

