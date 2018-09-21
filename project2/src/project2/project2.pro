TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib

SOURCES += \
        main.cpp \
    jacobi.cpp \
    test.cpp \
    quantum.cpp

LIBS += -larmadillo -llapack -lblas

HEADERS += \
    test.h \
    jacobi.h \
    quantum.h
