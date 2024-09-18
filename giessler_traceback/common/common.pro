TEMPLATE = lib
CONFIG += staticlib
TARGET = traceback
#CONFIG += console qt
#QT -= gui

SOURCES += \
    starparam.cpp \
    VecNorm.cpp

HEADERS += \
    Vector.h \
    Matrix.h \
    MatNorm.h \
    VecNorm.h \
    starparam.h \
    rng.h

#LIBS += -lgsl -lgslcblas -lrt
