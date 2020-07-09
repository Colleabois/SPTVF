#ifndef GLDISPLAYWIDGET_H
#define GLDISPLAYWIDGET_H

#include <QGLWidget>
#include <QtWidgets>
#include <QTimer>

class GLDisplayWidget : public QGLWidget
{
public:
    //for display changing
    int _choice;
    int _subchoice;

    // Translation (was in private)
    GLdouble _X, _Y, _Z;
    GLdouble _cameraX,_cameraY,_cameraZ,_theta,_distance;
    explicit GLDisplayWidget(QWidget *parent= nullptr); //before was O as null pointer

    void initializeGL();
    void paintGL(); // Display Gl
    void resizeGL(int width, int height);

    void selectPositionAndColorOn3DModel();         //TP3
    QVector3D unproject(int X, int Y, float Depth); //TP3

protected:
    // Mouse Management
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void wheelEvent(QWheelEvent *event);

private:

    QTimer _timer; // To update the scene

    float _angle; // Rotation

    QPoint _lastPosMouse; // To keep the last position of the mouse

};

#endif // GLDISPLAYWIDGET_H
