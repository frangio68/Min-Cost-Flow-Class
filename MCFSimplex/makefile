##############################################################################
################################## makefile ##################################
##############################################################################
#                                                                            #
#   makefile of MCFSimplex                                                   #
#                                                                            #
#   Input:  $(CC)       = compiler command                                   #
#           $(SW)       = compiler options                                   #
#           $(MCFClH)   = the include files for MCFClass                     #
#           $(MCFClINC) = the -I$( MCFClass.h directory )                    #
#           $(MCFSxDIR) = the directory where the source is                  #
#                                                                            #
#   Output: $(MCFSxOBJ) = the final object(s) / library                      #
#           $(MCFSxLIB) = external libraries + -L< libdirs >                 #
#           $(MCFSxH)   = the .h files to include                            #
#           $(MCFSxINC) = the -I$( include directories )                     #
#                                                                            #
#                              Antonio Frangioni                             #
#                          Dipartimento di Informatica                       #
#                              Universita' di Pisa                           #
#                                                                            #
##############################################################################

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(MCFSxDIR)/*.o $(MCFSxDIR)/*~

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

MCFSxOBJ = $(MCFSxDIR)/MCFSimplex.o
MCFSxLIB =
MCFSxH   = $(MCFSxDIR)/MCFSimplex.h
MCFSxINC = -I$(MCFSxDIR)

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(MCFSxOBJ): $(MCFSxDIR)/MCFSimplex.C $(MCFSxH) $(MCFClH)
	$(CC) -c $*.C -o $@ $(MCFClINC) $(SW)

############################# End of makefile ################################
