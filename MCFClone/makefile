##############################################################################
################################## makefile ##################################
##############################################################################
#									     #
#   makefile of MCFClone 						     #
#									     #
#   Input:  $(MCLDIR) = the directory where the source is		     #
#									     #
#   Output: $(MCLH)   = the .h files to include				     #
#           $(MCLINC) = the -I$(MCFClone directory)			     #
#									     #
#                                VERSION 1.60				     #
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
	rm -f $(MCLDIR)*~

# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

MCLH   = $(MCLDIR)MCFClone.h
MCLINC = -I$(MCLDIR)

# dependencies: every .o from its .C + every recursively included .h- - - - -

############################# End of makefile ################################
