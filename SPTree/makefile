##############################################################################
################################## makefile ##################################
##############################################################################
#									     #
#   makefile of SPTree 							     #
#									     #
#   Input:  $(CC)     = compiler command				     #
#           $(SW)     = compiler options                                     #
#           $(OPTH)   = the OPT****.h files to include			     #
#           $(OPTINC) = the -I$(OPT****.h directory)			     #
#           $(MCFH)   = the include files for MCFClass			     #
#           $(MCFINC) = the -I$(MCFClass.h directory)			     #
#           $(SPTDIR) = the directory where the source is		     #
#									     #
#   Output: $(SPTOBJ) = the final object(s) / library			     #
#           $(SPTLIB) = external libreries + -L<libdirs>		     #
#           $(SPTH)   = the .h files to include				     #
#           $(SPTINC) = the -I$(include directories), excluding OPTDIR	     #
#									     #
#                                VERSION 3.00				     #
#                	        07 - 05 - 2012				     #
#									     #
#			      Antonio Frangioni				     #
#   			   Operations Research Group			     #
#			  Dipartimento di Informatica			     #
#   			     Universita' di Pisa			     #
#									     #
##############################################################################

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(SPTDIR)*.o $(SPTDIR)*~

# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

SPTOBJ = $(SPTDIR)SPTree.o
SPTLIB =
SPTH   = $(SPTDIR)SPTree.h
SPTINC = -I$(SPTDIR)

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(SPTOBJ): $(SPTDIR)SPTree.C $(MCFH) $(OPTH)
	$(CC) -c $*.C -o $@ $(MCFINC) $(OPTINC) $(SW)

############################# End of makefile ################################
