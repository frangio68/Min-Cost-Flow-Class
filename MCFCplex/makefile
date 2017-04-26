##############################################################################
################################## makefile ##################################
##############################################################################
#									     #
#   makefile of MCFCplex 						     #
#   all modules are given in input                                           #
#									     #
#   Input:  $(CC)        = compiler command				     #
#           $(SW)        = compiler options                                  #
#           $(OPTH)      = the OPT****.h files to include		     #
#           $(OPTINC)    = the -I$(OPT****.h directory)			     #
#           $(MCFH)      = the include files for MCFClass		     #
#           $(MCFINC)    = the -I$(MCFClass.h directory)		     #
#           $(libCPXINC) = the -I$(include directories)	for libCPX	     #
#           $(MCXDIR)    = the directory where the source is		     #
#									     #
#   Output: $(MCXOBJ)    = the final object(s) / library		     #
#           $(MCXLIB)    = external libreries + -L<libdirs>		     #
#           $(MCXH)      = the .h files to include			     #
#           $(MCXINC)    = the -I$(include directories)			     #
#									     #
#                                VERSION 3.00				     #
#                	        12 - 05 - 2012				     #
#									     #
#			      Antonio Frangioni				     #
#   			   Operations Research Group			     #
#			  Dipartimento di Informatica			     #
#   			     Universita' di Pisa			     #
#									     #
##############################################################################

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(MCXDIR)*.o $(MCXDIR)*~

# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

MCXOBJ = $(MCXDIR)MCFCplex.o
MCXLIB = 
MCXH   = $(MCXDIR)MCFCplex.h
MCXINC = -I$(MCXDIR)

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(MCXOBJ): $(MCXDIR)MCFCplex.C $(MCXH) $(OPTH)
	$(CC) -c $*.C -o $@ $(MCFINC) $(libCPXINC) $(OPTINC) $(SW)

############################# End of makefile ################################
