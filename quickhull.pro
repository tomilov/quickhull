CONFIG -= qt
TEMPLATE = app

CONFIG *= exceptions rtti warn_on
QMAKE_CXXFLAGS_CXX11 = -std=gnu++1y
CONFIG *= c++11

DESTDIR = "bin/"

TARGET = "quickhull"

INCLUDEPATH += "include/"
DEPENDPATH += "include/"
LIBS += -L"lib/"
CONFIG(release, debug|release) {
}
CONFIG(debug, debug|release) {
    DEFINES += BOOST_UBLAS_NDEBUG=1
}

QMAKE_CXXFLAGS_DEBUG += -fno-default-inline -fno-inline
QMAKE_LFLAGS_DEBUG   += -fno-default-inline -fno-inline

win32 {
    CONFIG += console thread
}

SOURCES += "src/quickhull.cpp"
HEADERS += "include/quickhull.hpp"
