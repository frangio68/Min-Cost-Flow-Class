/*--------------------------------------------------------------------------*/
/*---------------------------- File Manual ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * User Manual for the "pure LGPL part" of the MCFClass project, a set of C++
 * solvers for (Linear or Convex Quadratic Separable) Min Cost Flow Problem
 * solvers under the same interface deriving from the base class MCFClass.
 *
 *  \author Antonio Frangioni \n
 *          Dipartimento di Informatica \n
 *          Universita' di Pisa \n
 *
 *  \author Alessandro Bertolini \n
 *          Dipartimento di Informatica \n
 *          Universita' di Pisa \n
 */

/** \mainpage The MCFClass Project Documentation

\section The MCFClass Project

\subsection intro Introduction

This is the Doxygen documentation for the "pure LGPL part" of MCFClass
project, a set of C++ solvers for (Linear or Convex Quadratic Separable)
Min Cost Flow Problem solvers under the same interface seriving from the
base class MCFClass.

The aim of MCFClass is to provide an abstraction layer between practitioners
who need to solve MCF problems within complex applications and developers of
MCF software. The idea is to provide an interface which caters for all the
needs that a practitioner can have, thereby allowing him/her to use whichever
algorithm - among those that have an implementation conforming to this
interface - without bothering with the details of the implementation, and to
easily switch between different algorithms.

MCFClass defines and "exports" the types of "flows" (MCFClass::FNumber),
"costs" (MCFClass::CNumber) and so on, together with a set of comparison
operators (ETZ, GTZ, ...) which automatically detect whether or not the
underlying types are integers or floats, inserting appropriate "epsilons" in
the latter case and avoiding them (for speed) in the former; these things are
sorted out at compile time without user intervention.

\subsection license License

This code is provided free of charge under the "GNU Lesser General Public 
License". There are actually two more complete solvers available under the
MCFClass interface, namely CS2 and MCFZIB. These are, however, distributed
under a more restrictive academic license, which has to be explicitly
accepted before getting hold of the code. Request forms are available in
the req folder


\subsection contents Contents

This release comprises:

-  docs/: this documentation, also available at

    https://frangio68.github.io/Min-Cost-Flow-Class/

-  doxygen/: doxygen files to produce the documentation

-  License.md: the text of the "GNU Lesser General Public
   License", Version 3.0, under which most of this code is distributed
   (but not all of it, see RelaxIV below)

-  MCFClass/: definition of the base class

-  MCFClone/: implements a "fake" MCF solver that takes two "real" ones
   and does everything on both; useful for testing the solvers (either for
   correctness or for efficiency) when used within "complex" approaches

-  MCFCplex/: implements a MCF solver conforming to the MCFClass interface
   based on calls to the commercial (but free for academic purposes)
   IBM/ILOG Cplex solver 

-  MCFSimplex/: implements a MCF solver conforming to the MCFClass interface
   based on the primal and dual revised network simplex algorithm

-  OPTUtils/: contains the `OPTUtils.h` file with a few minor utility
   functions

-  ReadMe.md: this file

-  RelaxIV/: implements a MCF solver conforming to the MCFClass interface
   based on the RELAXIV code by D. Bertsekas and P. Tseng, as described in
       Bertsekas, Dimitri P., and Paul Tseng.
       "RELAX-IV: A faster version of the RELAX code for solving minimum
       cost flow problems." (1994), Report LIDS-P-2276, MIT.
   be aware that RelaxIV is distributed under a less permissive academic
   license than the rest of the code (see RelaxIV/academicl.txt for details),
   which only applies to researchers of noncommercial and academic
   institutions, e.g., universities; if the license does not apply to you,
   you must not download the source or delete it immediately

-  req/: folder containing two request forms for two other solvers derived
   from the MCFClass interface that cannot be directly distributed in the
   repo due to licensing restrictions, see below

-  SPTree/: implements a MCF solver partly conforming to the MCFClass
   interface, in the sense that is only able to solve MCF instances that
   are in fact Shortest Path Tree ones (that is, only one source node and
   no arc capacities), but then does so using SPT algorithms (both
   label-setting and label-correcting variants can be used) that are much
   faster than complete MCF ones

-  pyMCFSimplex-0.9/: a Python-Wrapper for the MCFSimplex solver by Johannes
   from the G#.Blog, check the README.txt for details

-  test/: contains two example Main files to use the library. One solves
   a given MCF instance with any one MCF solver, which can be chosen by
   just changing two lines of code. The other compares the results of two
   solvers in order to verify that they agree. See the comments in both
   files for more details

There are two more complete solvers available under the MCFClass interface,
namely CS2 and MCFZIB. These are, however, distributed under a more
restrictive academic license, which has to be explicitly accepted before
getting hold of the code. Request forms are available in the req folder.

A further solver, MCFIntPnt (based on Interior-Point algorithms) has been 
developed, but it has not yet reached a sufficient maturuty to be distributed;
its principles are discussed at

  http://pages.di.unipi.it/frangio/abstracts.html#SIOPT04

  http://pages.di.unipi.it/frangio/abstracts.html#COAP06

  http://pages.di.unipi.it/frangio/abstracts.html#OMS06

For installation instruction, check the ReadMe. */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File Manual -------------------------------*/
/*--------------------------------------------------------------------------*/
