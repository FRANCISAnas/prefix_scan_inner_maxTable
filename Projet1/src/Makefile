#=======================================================================
# Basic C: makefile example to use as a reminder or as a template
#-----------------------------------------------------------------------
# Anas Francis
#=======================================================================

#only ok for project with no templated classes


#compiler
COMPILER = gcc
#linker
LINKER = gcc

#options for linker and compiler
FLAGS =-Wall -lm -fopenmp

EXE_NAME= out.mp

#can have several ones separated by spaces, only cpp files
SOURCES = main.c

#PATH to extra header used in the project (when using libs not installed in the OS)
INC_PATHS= -I

#extra libs used in the project (for instance -lpthread)
LIBS	= -L 

LINKOBJ	= $(SOURCES:.c=.o)
SOURCEHEADERS = $(SOURCES:.c=.h)

$(EXE_NAME):
	$(LINKER) -g -o $(EXE_NAME) $(SOURCES) $(FLAGS)
	valgrind time ./$(EXE_NAME) 512

%.o: %.c %.h
	$(COMPILER) $(FLAGS) -o $@ $(INC_PATHS) 

# cleanup
clean:
	rm -f *.o
	rm $(EXE_NAME)
