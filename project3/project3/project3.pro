TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib

SOURCES += \
        main.cpp \
    body.cpp \
    vec3.cpp \
    solver.cpp

LIBS += -larmadillo -llapack -lblas

HEADERS += \
    body.h \
    vec3.h \
    solver.h
