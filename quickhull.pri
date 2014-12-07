CONFIG -= qt
TEMPLATE = app

CONFIG *= exceptions rtti warn_on optimize_full
QMAKE_CXXFLAGS_CXX11 = -std=gnu++1z
CONFIG *= c++11

DESTDIR = "bin/"

INCLUDEPATH += "include/"
DEPENDPATH += "include/"

CONFIG(release, debug|release) {
    DEFINES += NDEBUG=1
}
CONFIG(debug, debug|release) {
    DEFINES += _DEBUG=1 DEBUG=1
    DEFINES += _GLIBCXX_DEBUG=1
}

QMAKE_CXXFLAGS_DEBUG += -fno-inline -fno-omit-frame-pointer -fno-optimize-sibling-calls
#QMAKE_LFLAGS_DEBUG   +=

win32 {
    CONFIG += console
}
