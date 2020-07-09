#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    //connect(ui-> Reset, SIGNAL(on_reset_released());ui-> resizerX ,SLOT(setHidden(bool));
    ui->resizer ->setMinimum(-10);   ui->resizer ->setMaximum(10);
    ui->resizerX->setMinimum(0);   ui->resizerX->setMaximum(100);
    ui->resizerY->setMinimum(-100);   ui->resizerY->setMaximum(100);
}

MainWindow::~MainWindow()
{

    delete ui;
}

void MainWindow::on_Exit_released()
{
    close();
}

//-------------------------------------------------------------
// Swapping displays
//-------------------------------------------------------------
void MainWindow::on_checkBox0_stateChanged(int arg1)
{
    int & mychoice = ui-> mydisplay ->_choice;
    if(arg1) mychoice = mychoice*2;
    else mychoice = mychoice/2;
}
void MainWindow::on_checkBox1_stateChanged(int arg1)
{
    int & mychoice = ui-> mydisplay ->_choice;
    if(arg1) mychoice = mychoice*3;
    else mychoice = mychoice/3;
}
void MainWindow::on_checkBox2_stateChanged(int arg1)
{
    int & mychoice = ui-> mydisplay ->_choice;
    if(arg1) mychoice = mychoice*5;
    else mychoice = mychoice/5;
}


//-------------------------------------------------------------
// Scrollbars
//-------------------------------------------------------------
void MainWindow::on_resizer_sliderMoved(int position){ui->mydisplay-> _theta = position;}
void MainWindow::on_resizerY_sliderMoved(int position){ui->mydisplay-> _cameraZ = position;}
void MainWindow::on_resizerX_sliderMoved(int position){ui->mydisplay-> _distance = position;}


//-------------------------------------------------------------
// RESET
//-------------------------------------------------------------

void MainWindow::on_Reset_released()
{
    ui->mydisplay-> _cameraX = 0;
    ui->mydisplay-> _cameraY = 0;
    ui->mydisplay-> _cameraZ = 0;
    //ui-> mydisplay ->_choice = 1;
}



