##############################################################################
################################## makefile ##################################
##############################################################################
#									     #
#   makefile of MCFSimplex 						     #
#									     #
#   Input:  $(CC)     = compiler command				     #
#           $(SW)     = compiler options                                     #
#           $(OPTH)   = the OPT****.h files to include			     #
#           $(OPTINC) = the -I$(OPT****.h directory)			     #
#           $(MCFH)   = the include files for MCFClass			     #
#           $(MCFINC) = the -I$(MCFClass.h directory)			     #
#           $(MFSDIR) = the directory where the source is		     #
#									     #
#   Output: $(MFSOBJ) = the final object(s) / library			     #
#           $(MFSLIB) = external libreries + -L<libdirs>		     #
#           $(MFSH)   = the .h files to include				     #
#           $(MFSINC) = the -I$(include directories), excluding OPTDIR	     #
#									     #
#                                VERSION 1.00				     #
#                	        09 - 05 - 2012				     #
#									     #
#			      Antonio Frangioni				     #
#   			   Operations Research Group			     #
#			  Dipartimento di Informatica			     #
#   			     Universita' di Pisa			     #
#									     #
##############################################################################

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(MFSDIR)*.o $(MFSDIR)*~

# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

MFSOBJ = $(MFSDIR)MCFSimplex.o
MFSLIB =
MFSH   = $(MFSDIR)MCFSimplex.h
MFSINC = -I$(MFSDIR)

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(MFSOBJ): $(MFSDIR)MCFSimplex.C $(MFSH) $(MCFH) $(OPTH)
	$(CC) -c $*.C -o $@ $(MCFINC) $(OPTINC) $(SW)

############################# End of makefile ################################
