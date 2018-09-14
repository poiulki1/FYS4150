TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib

SOURCES += \
        main.cpp \
    eigenvalue_solver.cpp \
    jacobi_test.cpp

LIBS += -larmadillo -llapack -lblas

HEADERS += \
    eigenvalue_solver.h
