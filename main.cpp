#include "mainwindow.h"
#include <QApplication>

// Polynomial approximation of discrete function
// Solving overdetermined system of linear equations with SVD matrix decomposition
// author: Wojciech Mormul
// date 14.01.2017

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    return a.exec();
}
