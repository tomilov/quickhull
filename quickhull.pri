CONFIG -= qt
TEMPLATE = app

CONFIG *= exceptions rtti warn_on optimize_full
QMAKE_CXXFLAGS_CXX11 = -std=gnu++1z
CONFIG *= c++11

linux:DESTDIR = /tmp/
linux:OBJECTS_DIR = /tmp/

INCLUDEPATH += "include/"
DEPENDPATH += "include/"

CONFIG(release, debug|release) {
    DEFINES += NDEBUG=1
}
CONFIG(debug, debug|release) {
    DEFINES += _DEBUG=1 DEBUG=1
    DEFINES += _GLIBCXX_DEBUG=1
}

#QMAKE_CXXFLAGS_RELEASE += -fno-inline -fno-omit-frame-pointer -fno-optimize-sibling-calls
#QMAKE_CXXFLAGS_RELEASE += -gline-tables-only
#QMAKE_LFLAGS_RELEASE   += -gline-tables-only
#QMAKE_LFLAGS_RELEASE   += -Wl,--no-as-needed -lprofiler -Wl,--as-needed

QMAKE_CXXFLAGS_WARN_ON = \
    -W -Weverything -pedantic -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded \
    -ftemplate-backtrace-limit=0 -fdiagnostics-color=always
