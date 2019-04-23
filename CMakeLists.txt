include_directories(MCFClass OPTUtils)
add_library(MCFClass STATIC
        MCFSimplex/MCFSimplex.C
        SPTree/SPTree.C
        )
