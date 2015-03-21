############################################################################
#
# Sample Makefile for GNU g++ under Linux
#
# From: 
# 1) Steve Oualline, "Practical C++ Programming", O'Reilly, 2003,
# Chapters 7 and 23;
#
# 2) Daoqi Yang, "C++ and Object Oriented Numeric Computing", Springer, 2001,
# section 4.4.4
#
# 3) Sample make program, Computational Physics 780, Ohio State University, 2009
# www.physics.ohio-state.edu/~ntg/780/programs/make_program
#
############################################################################
#
# This file contains a set of rules used by the "make"
#   command.  This makefile $(MAKEFILE) tells "make" how the
#   executable program $(COMMAND) should be created
#   from the source files $(SRCS) via the object
#   files $(OBJS) and the header files $(HDRS) (or equivalently, the required 
#   header files in $(INCLS)); type the command:
#        "make -f make_program"
#   where "make_program" should be replaced by the name of
#   the makefile, if not the default "Makefile" (in which case the -f option 
#   is not required)
# 
# This file also tells make what to do if arguments
#   follow the "make" command.
#
# To remove the OBJS files; type the command:
#        "make -f make_program clean"
#
# To create a gzipped tar file with name $(COMMAND).tar.gz
#  containing this makefile and the SRCS and HDRS files,
#  type the command:
#        "make -f make_program tarz"

#
############################################################################

# where to look for files if not in directory containing Makefile
INCLS  	= -I./include
SOURCE 	= ./source
MAIN  	= .

# macro for name of compiler
CC = g++

# macro for compilation options, if not default
CFLAGS = -g -O3
# -c       compile only, generate corresponding .o files
# -g       generate code suitable for the debugger
# -p       generate code suitable for the profiler
# -O       attemp code optimizations

# optional macro for compiler warnings, if not default
WARNFLAGS = -Wall

# name of executable file
COMMAND = hw2

# optional directory to find additional or user defined libraries
#LIBS	= -L/home/feehan/gcc/usr/local/lib # for user defined libraries

# optional math and gsl libraries
#LDFLAGS= -lgsl -lgslcblas -lm

# list of main and required source files
SRCS = 	$(MAIN)/hw2main.cpp \
        $(SOURCE)/Simulator.cpp \
				$(SOURCE)/Random2.cpp \
				$(SOURCE)/Arrays.cpp \
				$(SOURCE)/Random3.cpp \
				$(SOURCE)/Normals.cpp \
				$(SOURCE)/Statistics.cpp \
				$(SOURCE)/Option.cpp

# list of object files to be created
OBJS = 	hw2main.o \
        Simulator.o \
				Arrays.o \
				Random3.o \
				Random2.o \
				Normals.o \
				Statistics.o \
				Option.o


# optional list of header files for tarz command
#HDRS =

# GSL and math libraries (uncomment if required)
#LDFLAGS= -lgsl -lgslcblas -lm

# make the executable from the object files; add $(LIBS) or $(LDFLAGS) if required
$(COMMAND): $(OBJS) 
	@echo "linking ..."
	$(CC) $(CFLAGS) -o $(COMMAND) $(OBJS)

# make the object files; $(SRCS) can be replaced by $(HDRS) in first line
$(OBJS): $(SRCS)
	$(CC) $(CFLAGS) $(WARNFLAGS) $(INCLS) -c $(SRCS)

# dependencies list
hw2main.o : $(SRCS)
#Simulator.o : $(MAIN)/Simulator.cpp
#Arrays.o : $(SOURCE)/Arrays.cpp
#PayOffBarrier.o : $(MAIN)/PayOffBarrier.cpp

# remove the object files
clean:
	rm -f $(OBJS) $(COMMAND)

# create a gzipped tar file containing this makefile and the SRCS and HDRS files
tarz:
	tar zcf - $(MAKEFILE) $(SRCS) $(HDRS) > $(COMMAND).tar.gz
