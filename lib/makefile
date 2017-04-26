##############################################################################
################################# makefile ###################################
##############################################################################
#									     #
#   makefile of libMCF							     #
#   all modules are given in input                                           #
#									     #
#   Input:  $(CC)        = compiler command				     #
#           $(SW)        = compiler options                                  #
#           $(OPTH)      = the OPT****.h files to include		     #
#           $(OPTINC)    = the -I$(OPT****.h directory)			     #
#           $(libCPXINC) = the -I$(include directories)	for libCPX	     #
#           $(libMCFDIR) = the basic directory for the library		     #
#									     #
#   Output: $(MCFH)      = the include files for MCFClass		     #
#           $(MCFINC)    = the -I$(MCFClass.h directory)		     #
#           $(libMCFOBJ) = the library (that must be built)		     #
#           $(libMCFLIB) = the library + external libreries + -L<libdirs>    #
#           $(libMCFH)   = all the .h files of the MCF library		     #
#           $(libMCFINC) = all the -I$(include dirs) of the MCF library	     #
#									     #
#                                VERSION 2.10				     #
#                	        24 - 04 - 2017				     #
#									     #
#			      Antonio Frangioni				     #
#   			   Operations Research Group			     #
#			  Dipartimento di Informatica			     #
#   			     Universita' di Pisa			     #
#									     #
##############################################################################

# module name
libMCFOBJ = $(libMCFDIR)lib/libMCF.a

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -
# if a module is not used in the current configuration, just comment out the
# corresponding include line
# each module outputs some macros to be used here:
# *OBJ is the final object/library
# *LIB external libreries + -L<libdirs>
# *H   is the list of all include files
# *INC is the -I<include directories>

# MCFClass
MCFDIR = $(libMCFDIR)MCFClass/
include $(MCFDIR)makefile

# CS2
#CS2DIR = $(libMCFDIR)CS2/
#include $(CS2DIR)makefile

# MCFClone
MCLDIR = $(libMCFDIR)MCFClone/
include $(MCLDIR)makefile

# MCFCplex
MCXDIR = $(libMCFDIR)MCFCplex/
include $(MCXDIR)makefile

# MCFSimplex
MFSDIR = $(libMCFDIR)MCFSimplex/
include $(MFSDIR)makefile

# MCFZIB
#ZIBDIR = $(libMCFDIR)MCFZIB/
#include $(ZIBDIR)makefile

# RelaxIV
#RIVDIR = $(libMCFDIR)RelaxIV/
#include $(RIVDIR)makefile

# SPTree
SPTDIR = $(libMCFDIR)SPTree/
include $(SPTDIR)makefile

# main module (archiving phase) - - - - - - - - - - - - - - - - - - - - - - -

# object files (this is local)
lcllibMCFOBJ =	$(CS2OBJ)\
		$(MCXOBJ)\
		$(MFSOBJ)\
		$(ZIBOBJ)\
		$(SPTOBJ)\
		$(RIVOBJ)
# libreries
libMCFLIB =	$(CS2LIB)\
		$(MCXLIB)\
		$(MFSLIB)\
		$(ZIBLIB)\
		$(SPTLIB)\
		$(RIVLIB)
# includes
libMCFH   =	$(OPTH)\
		$(MCFH)\
		$(CS2H)\
		$(MCLH)\
		$(MCXH)\
		$(MFSH)\
		$(ZIBH)\
		$(RIVH)\
		$(SPTH)
# include dirs
libMCFINC =	$(MCFINC)\
		$(CS2INC)\
		$(MCLINC)\
		$(MCXINC)\
		$(MFSINC)\
		$(ZIBINC)\
		$(RIVINC)\
		$(SPTINC)

$(libMCFOBJ): $(lcllibMCFOBJ)
	ar -ruv $(libMCFOBJ) $(lcllibMCFOBJ)

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(libMCFOBJ)

############################ End of makefile #################################
