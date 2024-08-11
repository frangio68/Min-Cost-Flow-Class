##############################################################################
################################# makefile ###################################
##############################################################################
#                                                                            #
#   makefile of MCFSolve                                                     #
#                                                                            #
#                              Antonio Frangioni                             #
#                          Dipartimento di Informatica                       #
#                              Universita' di Pisa                           #
#                                                                            #
##############################################################################

# basic directory
DIR = .

# basic libraries
LIB = -lm

# debug switches
SW = -g -std=c++14 -Wall
# production switches
#SW = -O3 -DNDEBUG -std=c++14

# compiler
CC = clang++

# default target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

default: $(DIR)/MCFSolve $(DIR)/TestMCF

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -
# if a module is not used in the current configuration, just comment out the
# corresponding include line
# each module outputs some macros to be used here:
# *OBJ is the final object/library
# *LIB external libraries + -L< libdirs >
# *H   is the list of all include files
# *INC is the -I< include directories >

# the whole libMCF, with all required external libraries
libMCFClDIR = ..
include $(libMCFClDIR)/lib/makefile-inc

# main module (linking phase) - - - - - - - - - - - - - - - - - - - - - - - -

$(DIR)/MCFSolve: $(DIR)/Main.o $(libMCFClOBJ)
	$(CC) -o $(DIR)/MCFSolve $(DIR)/Main.o $(libMCFClOBJ) \
	$(libMCFClLIB) $(LIB) $(SW)

$(DIR)/TestMCF: $(DIR)/TestMain.o $(libMCFOBJ)
	$(CC) -o $(DIR)/TestMCF $(DIR)/TestMain.o $(libMCFClOBJ) \
	$(libMCFClLIB) $(LIB) $(SW)

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(DIR)/Main.o: $(DIR)/Main.C $(libMCFClH)
	$(CC) -c $*.C -o $@ $(libMCFClINC) $(SW)

$(DIR)/TestMain.o: $(DIR)/TestMain.C $(libMCFClH)
	$(CC) -c $*.C -o $@ $(libMCFClINC) $(SW)

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -Rf $(DIR)/MCFSolve $(DIR)/TestMCF $(DIR)/*.o $(DIR)/*~ \
	$(DIR)/*.dSYM

############################ End of makefile #################################
