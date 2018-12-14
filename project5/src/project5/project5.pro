TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib
SOURCES += \
        main.cpp \
        solver.cpp \
    savetofile.cpp \
    read_parameters.cpp
        LIBS += -larmadillo -llapack -lblas

HEADERS += \
    solver.h \
    savetofile.h \
    read_parameters.h
