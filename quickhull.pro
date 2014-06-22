CONFIG -= qt
TEMPLATE = app

CONFIG *= exceptions rtti warn_on
QMAKE_CXXFLAGS_CXX11 = -std=gnu++1y
CONFIG *= c++11

DESTDIR = "bin/"

TARGET = "quickhull"

INCLUDEPATH += "include/"
DEPENDPATH += "include/"

CONFIG(release, debug|release) {
    DEFINES += NDEBUG=1
}
CONFIG(debug, debug|release) {
    DEFINES += _DEBUG=1 DEBUG=1
    DEFINES += _GLIBCXX_DEBUG=1
}

QMAKE_CXXFLAGS_DEBUG += -fno-default-inline -fno-inline
QMAKE_LFLAGS_DEBUG   += -fno-default-inline -fno-inline
#QMAKE_CXXFLAGS_RELEASE += -g -pg
#QMAKE_LFLAGS_RELEASE   += -g -pg

win32 {
    CONFIG += console
#thread
}

SOURCES += "src/quickhull.cpp"
HEADERS += "include/quickhull.hpp"
