#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    origin_centre.x=650;origin_centre.y=350;
    x_axis_length=620;y_axis_length=320;
    scale=60;
    cross_length=12;

    x_axis_arrow=new QPointF[3];
    x_axis_arrow[0]= QPointF(origin_centre.x+x_axis_length, origin_centre.y+10);
    x_axis_arrow[1]=QPointF(origin_centre.x+x_axis_length, origin_centre.y-10);
    x_axis_arrow[2]=QPointF(origin_centre.x+x_axis_length+20, origin_centre.y);

    y_axis_arrow=new QPointF[3];
    y_axis_arrow[0]= QPointF(origin_centre.x+10, origin_centre.y-y_axis_length);
    y_axis_arrow[1]=QPointF(origin_centre.x-10, origin_centre.y-y_axis_length);
    y_axis_arrow[2]=QPointF(origin_centre.x, origin_centre.y-y_axis_length-20);

    ui->degree_spin_box->setValue(1);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_approximate_clicked()
{
    int M=real_points.size();
    int N=ui->degree_spin_box->value()+1;

    if(M>=N)
    {
        A=new Matrix(M,N,"A",Matrix::as_matrix);
        b=new Matrix(M,1,"b",Matrix::as_vector);
        c=new Matrix(N,1,"c",Matrix::as_vector);

        for(int i=0;i<M;i++)
        {
            b->mat[i]=real_points[i].y;

            for(int j=0;j<N;j++)
                A->mat[i*A->cols+j]=pow(real_points[i].x,j);
        }

        A_pinv = svd_deomposition.pinv_compute(A,A->rows,A->cols);
        svd_deomposition.multiply(A_pinv,b,c);

        A->draw();
        b->draw();

        svd_deomposition.U->draw();
        svd_deomposition.S->draw();
        svd_deomposition.VT->draw();

        A_pinv->draw();
        c->draw();

        polynomial_draw=1;
        update();
    }
}

void MainWindow::paintEvent(QPaintEvent *event)
{
    draw_origin();

    if(polynomial_draw==1)
        draw_polynomial();

    draw_crosses();
}

void MainWindow::draw_origin()
{
    QPainter painter(this);

    pen.setColor(QColor(100,100,100,255));
    pen.setWidth(1);
    painter.setPen(pen);

    for(int i=-10;i<=10;i=i+1)
    {
        painter.drawLine(origin_centre.x-x_axis_length,origin_centre.y+i*scale,origin_centre.x+x_axis_length,origin_centre.y+i*scale);
        painter.drawLine(origin_centre.x+i*scale,origin_centre.y-y_axis_length,origin_centre.x+i*scale,origin_centre.y+y_axis_length);
    }

    pen.setColor(QColor(0,0,0,255));
    pen.setWidth(3);
    painter.setPen(pen);

    painter.drawLine(origin_centre.x-x_axis_length,origin_centre.y,origin_centre.x+x_axis_length,origin_centre.y);
    painter.drawLine(origin_centre.x,origin_centre.y-y_axis_length,origin_centre.x,origin_centre.y+y_axis_length);

    for(int i=-10;i<=10;i=i+1)
    {
        painter.drawLine(origin_centre.x-10,origin_centre.y+i*scale,origin_centre.x+10,origin_centre.y+i*scale);
        painter.drawLine(origin_centre.x+i*scale,origin_centre.y-10,origin_centre.x+i*scale,origin_centre.y+10);
    }

    painter.setBrush(QBrush(QColor(0,0,0,255)));
    painter.drawConvexPolygon(x_axis_arrow, 3);

    font.setPointSize (20);
    painter.setFont(font);
    painter.drawText(QPointF(origin_centre.x+x_axis_length+8,origin_centre.y-10),"X");
    painter.drawConvexPolygon(y_axis_arrow,3);
    painter.drawText(QPointF(origin_centre.x+10,origin_centre.y-y_axis_length),"Y");
}

void MainWindow::draw_crosses()
{
    QPainter painter(this);

    pen.setColor(QColor(0,100,255,255));
    pen.setWidth(4);
    painter.setPen(pen);

    for(int i=0;i<drawing_points.size();i++)
    {
      painter.drawLine(drawing_points[i].x-cross_length,drawing_points[i].y,drawing_points[i].x+cross_length,drawing_points[i].y);
      painter.drawLine(drawing_points[i].x,drawing_points[i].y-cross_length,drawing_points[i].x,drawing_points[i].y+cross_length);
    }

    pen.setColor(QColor(0,0,0,255));
    pen.setWidth(4);
    painter.setPen(pen);
    painter.setBrush(QBrush(QColor(239,235,231,255)));
    painter.drawRect(QRect(1140,480,143,200));
}

void MainWindow::draw_polynomial()
{
    QPainter painter(this);

    pen.setColor(QColor(255,0,0,255));
    pen.setWidth(4);
    painter.setPen(pen);

    float value;
    float prev_value=0;
    int degree=ui->degree_spin_box->value()+1;

    for(int n=0;n<degree;n++)
        prev_value=prev_value-c->mat[n]*pow(-10.21,n);

    for(float i=-10.2;i<10.2;i=i+0.01)
    {
        value=0;
        for(int n=0;n<degree;n++)
            value=value-c->mat[n]*pow(i,n);

        painter.drawLine((i-0.01)*scale+origin_centre.x,prev_value*scale+origin_centre.y,i*scale+origin_centre.x,value*scale+origin_centre.y);
        prev_value=value;
    }
}

void MainWindow::mousePressEvent(QMouseEvent * ev)
{
    Point mouse_pos;
    mouse_pos.x=ev->x();
    mouse_pos.y=ev->y();

    drawing_points.push_back(mouse_pos);
    mouse_pos.x=(mouse_pos.x-origin_centre.x)/scale;
    mouse_pos.y=-(mouse_pos.y-origin_centre.y)/scale;
    real_points.push_back(mouse_pos);

    update();
}

void MainWindow::on_exit_clicked()
{
    exit(0);
}

void MainWindow::on_clear_clicked()
{
    drawing_points.resize(0);
    real_points.resize(0);

    polynomial_draw=0;
    update();
}
