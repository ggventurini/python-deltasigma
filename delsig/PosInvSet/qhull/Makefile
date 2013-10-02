# Makefile for qhull and rbox
#Commands added for compilation on HP computers and for the MEX interface. 
# This makefile should be called from ../.. in order to get definitions
# for INCLUDE and SFX
# The MEX file is automatically moved to the parent directory.
#See the qhull.m file for a MATLAB synopsis of the function.
#(R. Schreier, July 1995)
#
#       see README
#
#	make           to produce qhull, rbox, qhull.1, and rbox.1
#	make qhullx    to produce qhull w/o using libqhull.a
#	make doc       to print documentation
#	make install   to copy qhull, rbox, qhull.1, rbox.1 to BINDIR, MANDIR
#	make new       to rebuild qhull, rbox, qhull.1, and rbox.1 from source
#
#	make printall  to print all files
#	make user_eg   to produce user_eg
#	make clean     to remove object files and core
#	make cleanall  to remove all generated files
#
#	PRINTMAN --  command for printing manual pages
#	PRINTC --  command for printing C files
#       BINDIR -- directory where to copy executables
#       MANDIR -- directory where to copy manual pages
#       CC --     ANSI C compiler
#       CCOPTS1 - options used to compile .c files
#       CCOPTS2 -- options used to link .o files
#
#       CFILES -- .c files for printing
#       HFILES -- .h files for printing
#       DFILES -- documentation files
#	MFILES -- man pages
#	FILES -- all other files
#       OBJS -- specifies the object files of libqhull.a
#
BINDIR  = ./bin
MANDIR  = ./man

# if you do not have enscript, try a2ps or just use lpr.  The files are text.
PRINTMAN = enscript -2rl
PRINTC = enscript -2r
# PRINTMAN = lpr
# PRINTC = lpr

#for gcc compiler -O2 for optimization, -g for debugging, -Wall for check
#
#CC     = gcc
#CCOPTS1 = -O2 -ansi 

#for HP's cc compiler -Aa for ANSI mode, -O for optimization(NOT YET)
# +z for PIC (needed for cmex to work)
#
CC     = cc
CCOPTS1 = -Aa +z

# for Sun's cc compiler, -fast or O2 for optimization, -g for debugging, -Xc for ANSI
#CC = cc
#CCOPTS1 = -Xc -v -fast

# for Silicon Graphics cc compiler, -O2 for optimization, -g for debugging
#CC = cc
#CCOPTS1 = -ansi -O2

# for Next cc compiler with fat executable
#CC = cc
#CCOPTS1 = -ansi -O2 -arch m68k -arch i386 -arch hppa

# for loader, ld
CCOPTS2 = $(CCOPTS1)

# OBJS in execution frequency order.  CFILES after qhull.c are alphabetical
OBJS = global.o stat.o geom2.o poly2.o io.o merge.o qhull.o \
       mem.o set.o poly.o geom.o

CFILES= unix.c qhull.c geom.c geom2.c global.c io.c mem.c merge.c poly.c \
        poly2.c set.c stat.c
HFILES= user.h qhull.h qhull_a.h geom.h io.h mem.h merge.h poly.h set.h stat.h
DFILES=  README Announce Changes REGISTER COPYING
FILES=  Makefile qhull.man rbox.man rbox.c user_eg.c q_test q_egtest q_eg
MFILES= qhull.1 rbox.1

# These are overridden when called from ../../Makefile
SFX = mexhp7
INCLUDE = /usr/local/apps/matlab/extern/include

qhull.${SFX}: mex.o libqhull.a
	cmex -I${INCLUDE} mex.o libqhull.a
	mv mex.${SFX} ../$@
 
mex.o: mex.c
	cc -c -I${INCLUDE} $(CCOPTS1) mex.c

all: rbox qhull 

unix.o:   qhull.h user.h
qhull.o:  qhull.h stat.h qhull_a.h user.h
geom.o:   qhull.h stat.h geom.h qhull_a.h user.h
geom2.o:   qhull.h stat.h geom.h qhull_a.h user.h
global.o: qhull.h stat.h qhull_a.h user.h
io.o:     qhull.h stat.h io.h qhull_a.h user.h
mem.o:    mem.h 
merge.o:  qhull.h stat.h merge.h qhull_a.h user.h
poly.o:   qhull.h stat.h poly.h qhull_a.h user.h
poly2.o:   qhull.h stat.h poly.h qhull_a.h user.h
set.o:    set.h mem.h 
stat.o:   qhull.h stat.h qhull_a.h user.h

.c.o:
	$(CC) -c $(CCOPTS1) $<

clean:
	rm -f *.o core 

cleanall: clean
	rm -f rbox qhull user_eg libqhull.a >/dev/null

doc:
	$(PRINTMAN) $(DFILES) $(MFILES)

install: all
	cp qhull $(BINDIR)/qhull
	cp rbox $(BINDIR)/rbox
	cp qhull.man $(MANDIR)/qhull.man
	cp rbox.man $(MANDIR)/rbox.man

new:    cleanall all

printall: doc printh printc printf

printc:
	$(PRINTC) $(CFILES)

printf:
	$(PRINTC) $(FILES) 

printh:
	$(PRINTC) $(HFILES)

libqhull.a: $(OBJS)
	#echo if 'ar' or 'ranlib' fails, try 'make qhullx' >/dev/null
	ar r libqhull.a $(OBJS)
	#echo the next line may need to be removed. >/dev/null
	#-test -x /bin/ranlib -o -x /usr/bin/ranlib && ranlib libqhull.a

qhull: unix.o libqhull.a
	$(CC) -o qhull $(CCOPTS2) unix.o -L. -lqhull -lm 
	-chmod +x q_test q_eg q_egtest
	-./rbox D4 | ./qhull

# compile qhull without using libqhull.a
qhullx: unix.o $(OBJS)
	$(CC) -o qhull $(CCOPTS2) unix.o $(OBJS) -lm 
	-chmod +x q_test q_eg q_egtest
	-./rbox D4 | ./qhull

qhull.1: qhull.man
	nroff -man qhull.man | col -b >qhull.1

rbox: rbox.o
	$(CC) -o rbox rbox.o $(CCOPTS2) -lm

rbox.1: rbox.man
	nroff -man rbox.man | col -b >rbox.1

user_eg: user_eg.o libqhull.a 
	$(CC)  -o user_eg $(CCOPTS2) user_eg.o  -L. -lqhull -lm 

# end of Makefile








