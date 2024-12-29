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

# common flags
COMMON_SW = -std=c++14

# debug switches
SW_DEBUG = -g -Wall $(COMMON_SW)
# production switches
SW_RELEASE = -O3 -DNDEBUG $(COMMON_SW)

# compiler
CC = clang++

# default target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

default: release

# debug target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

debug: SW = $(SW_DEBUG)
debug: $(DIR)/MCFSolve $(DIR)/TestMCF

# release target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

release: SW = $(SW_RELEASE)
release: $(DIR)/MCFSolve $(DIR)/TestMCF

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
include $(libMCFClDIR)/lib/makefile-lib

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

# clean target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -Rf $(DIR)/MCFSolve $(DIR)/TestMCF $(DIR)/*.o $(DIR)/*~ \
	$(DIR)/*.dSYM

# distclean target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

distclean: clean

# phony targets - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.PHONY: debug release clean distclean

############################ End of makefile #################################
