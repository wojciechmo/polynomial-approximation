QT       += core gui
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = aproximation
TEMPLATE = app
CONFIG += console

SOURCES += main.cpp\
        mainwindow.cpp \
    matrix.cpp \
    svd.cpp

HEADERS  += mainwindow.h \
    matrix.h \
    svd.h

FORMS    += mainwindow.ui
