ROOTCFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS   = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS  = $(shell $(ROOTSYS)/bin/root-config --glibs)

CXXFLAGS   = ${ROOTCFLAGS} -Wall #-g
LIBS       = ${ROOTLIBS} ${ROOTGLIBS}

SRCS	= $(wildcard *.cc)
HEADERS	= $(wildcard *.h)
EXESP   = main
ALLOBJS = $(addprefix bin/obj/, $(SRCS:.cc=.o))
OBJSP	= $(filter-out $(addsuffix .o, $(EXESP)), $(SRCS:.cc=.o))
OBJS    = $(addprefix bin/obj/, $(OBJSP))
EXES    = $(addprefix bin/, $(addsuffix .exe, $(EXESP)))
EXELN   = EMuBackgroundStudy.exe

all: ${EXES}
	@ln -sf bin/main.exe ${EXELN}

bin/obj/%.o: %.cc ${HEADERS}
	@mkdir -p bin/obj
	g++ ${CXXFLAGS} -c $< -o $@

bin/%.exe: ${OBJS} bin/obj/%.o
	g++ ${CXXFLAGS} ${LIBS} $^ -o $@

clean:
	rm ${EXELN}
	rm -f bin/*.exe bin/obj/*.o gmon.out
	rmdir bin/obj bin

.SECONDARY: ${ALLOBJS}
