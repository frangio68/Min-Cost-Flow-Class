##############################################################################
################################## makefile ##################################
##############################################################################
#                                                                            #
#   makefile of MCFCplex                                                     #
#   all modules are given in input                                           #
#                                                                            #
#   Input:  $(CC)          = compiler command                                #
#           $(SW)          = compiler options                                #
#           $(MCFClH)      = the include files for MCFClass                  #
#           $(MCFClINC)    = the -I$(MCFClass.h directory)                   #
#           $(libCPLEXINC) = the -I$(include directories) for libCPLEX       #
#           $(MCFCxDIR)    = the directory where the source is               #
#                                                                            #
#   Output: $(MCFCxOBJ)    = the final object(s) / library                   #
#           $(MCFCxLIB)    = external libraries + -L< libdirs >              #
#           $(MCFCxH)      = the .h files to include                         #
#           $(MCFCxINC)    = the -I$(include directories)                    #
#                                                                            #
#                              Antonio Frangioni                             #
#                          Dipartimento di Informatica                       #
#                              Universita' di Pisa                           #
#                                                                            #
##############################################################################

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(MCFCxDIR)/*.o $(MCFCxDIR)/*~

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

MCFCxOBJ = $(MCFCxDIR)/MCFCplex.o
MCFCxLIB = 
MCFCxH   = $(MCFCxDIR)/MCFCplex.h
MCFCxINC = -I$(MCFCxDIR)

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(MCFCxOBJ): $(MCFCxDIR)/MCFCplex.C $(MCFCxH) $(MCFClH)
	$(CC) -c $*.C -o $@ $(MCFClINC) $(libCPLEXINC) $(SW)

############################# End of makefile ################################
