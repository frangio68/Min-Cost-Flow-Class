##############################################################################
################################# makefile ###################################
##############################################################################
#									     #
#   makefile of MCFSolve						     #
#									     #
#   'make clean' cleans up						     #
#   'make' or 'make MCFSolve' builds the module				     #
#									     #
#                                VERSION 3.10				     #
#                	        16 - 05 - 2012				     #
#									     #
#			      Antonio Frangioni				     #
#   			   Operations Research Group			     #
#			  Dipartimento di Informatica			     #
#   			     Universita' di Pisa			     #
#									     #
##############################################################################

# basic directory
DIR = ./

# basic libreries
LIB = -lm

# debug switches
SW = -Wall -g -std=c++14

# production switches
#SW = -O3

# compiler
CC = g++

# default target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

default: MCFSolve TestMCF

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -
# if a module is not used in the current configuration, just comment out the
# corresponding include line
# each module outputs some macros to be used here:
# *OBJ is the final object/library
# *LIB external libreries + -L<libdirs>
# *H   is the list of all include files
# *INC is the -I<include directories>

# OPT****.h
OPTH	= $(DIR)OPTUtils.h
OPTINC	= -I$(DIR)

# the whole libMCF, with all required external libreries
libMCFDIR = ../
include $(libMCFDIR)lib/makefile-c

# main module (linking phase) - - - - - - - - - - - - - - - - - - - - - - - -

MCFSolve: Main.o $(libMCFOBJ)
	$(CC) -o MCFSolve Main.o $(libMCFOBJ) $(libMCFLIB) $(LIB) $(SW)

TestMCF: TestMain.o $(libMCFOBJ)
	$(CC) -o TestMCF TestMain.o $(libMCFOBJ) $(libMCFLIB) $(LIB) $(SW)

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(DIR)Main.o: $(DIR)Main.C $(libMCFH)
	$(CC) -c $*.C -o $@ $(libMCFINC) $(OPTINC) $(SW)

$(DIR)TestMain.o: $(DIR)TestMain.C $(libMCFH)
	$(CC) -c $*.C -o $@ $(libMCFINC) $(OPTINC) $(SW)

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -Rf MCFSolve TestMCF *.o *~ *.dSYM

############################ End of makefile #################################
