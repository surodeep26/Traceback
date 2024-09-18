TEMPLATE = app
CONFIG += console qt
QT -= gui

INCLUDEPATH += ../../common

DEPENDPATH += ../../common

PRE_TARGETDEPS = ../../common/libtraceback.a

LIBS += -L../../common -ltraceback -lgsl -lgslcblas -lrt

#QMAKE_RPATHDIR = $$replace(PWD,example,common)
