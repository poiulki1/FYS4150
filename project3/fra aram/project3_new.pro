TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    body.cpp \
    solver.cpp \
    quantities.cpp \
    store_values.cpp \
    vec3.cpp \
    solver.cpp

HEADERS += \
    body.h \
    solver.h \
    quantities.h \
    store_values.h \
    vec3.h
