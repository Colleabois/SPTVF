#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

private:
    Ui::MainWindow *ui;

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_Exit_released();
    void on_checkBox0_stateChanged(int arg1);
    void on_checkBox1_stateChanged(int arg1);
    void on_checkBox2_stateChanged(int arg1);

    void on_resizer_sliderMoved(int position);
    void on_resizerY_sliderMoved(int position);
    void on_resizerX_sliderMoved(int position);
    void on_Reset_released();
};

#endif // MAINWINDOW_H
