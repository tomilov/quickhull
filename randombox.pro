include(quickhull.pri)

CONFIG(release, debug|release) {
    LIBS += -lboost_program_options-mt-s
}
CONFIG(debug, debug|release) {
    LIBS += -lboost_program_options-mt-sd
}

TARGET = "randombox"

SOURCES += "src/randombox.cpp"
