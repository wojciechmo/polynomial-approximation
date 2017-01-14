#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QPainter>
#include <QMouseEvent>
#include <cmath>
#include "matrix.h"
#include "svd.h"

using namespace std;

struct Point{float x;float y;};

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    void paintEvent(QPaintEvent *event);
    void mousePressEvent(QMouseEvent * ev);
    void draw_origin();
    void draw_crosses();
    void draw_polynomial();

private slots:
    void on_exit_clicked();
    void on_clear_clicked();
    void on_approximate_clicked();

private:
    Ui::MainWindow *ui;
    QPen pen;
    QFont font;

    Matrix *A;
    Matrix *b;
    Matrix *c;
    Matrix *A_pinv;

    SVD svd_deomposition;

    Point origin_centre;
    float scale;
    int x_axis_length;
    int y_axis_length;
    int cross_length;
    QPointF *x_axis_arrow;
    QPointF *y_axis_arrow;
    vector <Point> drawing_points;
    vector <Point> real_points;
    int polynomial_draw=0;
};

#endif // MAINWINDOW_H
