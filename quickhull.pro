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
    #LIBS += -linsituc
}
CONFIG(debug, debug|release) {
    #LIBS += -linsitucd
}

win32 {
    CONFIG += console thread

    BOOST_INCLUDE="D:/libs/boost"
    INCLUDEPATH += $$BOOST_INCLUDE
    QMAKE_CXXFLAGS += -isystem $$BOOST_INCLUDE
    DEPENDPATH += $$BOOST_INCLUDE
    #QMAKE_INCDIR_QT =
}

SOURCES += "src/quickhull.cpp"
HEADERS += "include/quickhull.hpp"
