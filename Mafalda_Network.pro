HEADERS += \
	bcr.h \
	cell.h \
	chemokines3d.h \
	dynarray.h \
	GC3D.h \
	lattice.h \
	mafalda.h \
	main.h \
	output.h \
	parameters.h \
	random.h \
	setparam.h \
	trackball.h \
	vector3d.h \
        SHM.h\
	events.h\
	network.h\
	odesolver.h

SOURCES += \
	bcr.cpp \
	cell.cpp \
	chemokines3D.cpp \
	GC3D.cpp \
	lattice.cpp \
	mafalda.cpp \
	main.cpp \
	output.cpp \
	parameters.cpp \
	random.cpp \
	setparam.cpp \
	trackball.cpp \
	vector3d.cpp \
        SHM.cpp\
	events.cpp\
	network.cpp\
	odesolver.cpp


QT += opengl

mac: LIBS += -lil -framework GLUT

unix:!mac {LIBS += -lglut -lGLU -lGL -lSOIL -lILU -lIL -lILUT -lstdc++fs}

win32: LIBS += -L$$PWD/freeglut/lib/ -lfreeglut -lopengl32
INCLUDEPATH += $$PWD/freeglut/include
#DEPENDPATH += $$PWD/freeglut/lib

win32: LIBS += -L$$PWD/soil/ -lSOIL
INCLUDEPATH += $$PWD/soil
DEPENDPATH += $$PWD/soil

win32: LIBS += -L$$PWD/glStatic/ -lglu32
INCLUDEPATH += $$PWD/glStatic
DEPENDPATH += $$PWD/glStatic

win32: LIBS += -L$$PWD/glStatic/ -lopengl32
INCLUDEPATH += $$PWD/glStatic
DEPENDPATH += $$PWD/glStatic

#to use vector<doule> = {1,2};
QMAKE_CXXFLAGS += -std=c++17

mac: QMAKE_CXXFLAGS += -Wno-deprecated


mac: LIBS += -L$$PWD/../../../../../../../usr/local/Cellar/devil/1.8.0_1/lib/ -lIL -lILU -lILUT
mac:INCLUDEPATH += $$PWD/../../../../../../../usr/local/Cellar/devil/1.8.0_1/include/
mac:DEPENDPATH += $$PWD/../../../../../../../usr/local/Cellar/devil/1.8.0_1

CONFIG += c++17
