# Min-Cost-Flow-Class

ReadMe for the MCFClass project, a set of C++ solvers for (Linear or Convex
Quadratic Separable) Min Cost Flow Problem solvers under the same interface
deriving from the base class `MCFClass`.

The aim of `MCFClass` is to provide an abstraction layer between practitioners
who need to solve MCF problems within complex applications and developers of
MCF software. The idea is to provide an interface which caters for all the
needs that a practitioner can have, thereby allowing him/her to use whichever
algorithm - among those that have an implementation conforming to this
interface - without bothering with the details of the implementation, and to
easily switch between different algorithms.

`MCFClass` defines and "exports" the types of "flows" (`MCFClass::FNumber`),
"costs" (`MCFClass::CNumber`) and so on, together with a set of comparison
operators (ETZ, GTZ, ...) which automatically detect whether or not the
underlying types are integers or floats, inserting appropriate "epsilons" in
the latter case and avoiding them (for speed) in the former; these things are
sorted out at compile time without user intervention.

This release comprises:

-  [`docs/`](docs): HTML doxygen documentation, also available at

    https://frangio68.github.io/Min-Cost-Flow-Class/

-  [`doxygen/`](doxygen): files to produce the documentation

-  [`License.md`](License.md): the text of the "GNU Lesser General Public License",
   Version 3.0, under which most of this code is distributed
   (but not all of it, see RelaxIV below)

-  [`MCFClass/`](MCFClass): definition of the base class

-  [`MCFClone/`](MCFClone): implements a "fake" MCF solver that takes two "real" 
   ones  and does everything on both; useful for testing the solvers (either for
   correctness or for efficiency) when used within "complex" approaches

-  [`MCFCplex/`](MCFCplex): implements a MCF solver conforming to the `MCFClass`
   interface based on calls to the commercial (but free for academic purposes)
   IBM/ILOG Cplex solver 

-  [`MCFSimplex/`](MCFSimplex): implements a MCF solver conforming to the `MCFClass`
   interface based on the primal and dual revised network simplex algorithm

-  [`OPTUtils/`](OPTUtils): contains the `OPTUtils.h` file with a few minor utility
   functions

-  [`ReadMe.md`](ReadMe.md): this file

-  [`RelaxIV/`](RelaxIV): implements a MCF solver conforming to the `MCFClass`
   interface based on the RELAXIV code by D. Bertsekas and P. Tseng, as described
   in

       Bertsekas, Dimitri P., and Paul Tseng.
       "RELAX-IV: A faster version of the RELAX code for solving minimum
       cost flow problems." (1994), Report LIDS-P-2276, MIT

   be aware that `RelaxIV` is distributed under a [less permissive academic
   license](RelaxIV/academicl.txt), which only applies to researchers of
   noncommercial and academic institutions, e.g., universities; if the license
   does not apply to you, you must not download the source or delete it immediately

-  [`req/`](req): folder containing two request forms for two other solvers derived
   from the `MCFClass` interface that cannot be directly distributed in the
   repo due to licensing restrictions, see below
 
-  [`SPTree/`](SPTree): implements a MCF solver partly conforming to the `MCFClass`
   interface, in the sense that is only able to solve MCF instances that are in
   fact Shortest Path Tree ones (that is, only one source node and no arc
   capacities), but then does so using SPT algorithms (both label-setting and
   label-correcting variants can be used) that are much faster than complete MCF
   ones

-  [`pyMCFSimplex-0.9/`](pyMCFSimplex-0.9): a Python-Wrapper for the MCFSimplex
   solver by Johannes from the G#.Blog, check `README.txt` for details

-  [`test/`](test): contains two example Main files to use the library. One solves
   a given MCF instance with any one MCF solver, which can be chosen by just
   changing two lines of code. The other compares the results of two solvers in
   order to verify that they agree. See the comments in both files for more details

There are two more complete solvers available under the `MCFClass` interface,
namely CS2 and MCFZIB. These are, however, distributed under a more
restrictive academic license, which has to be explicitly accepted before
getting hold of the code. Request forms are available in the [`req/`](req) folder.


## Build and install

You can either use [CMake](https://cmake.org) or plain makefiles to build the
library, your choice. CMake compiles off-source and it is therefore perhaps
better suited to one-off, compile-and-forget installations, whereby the
provided makefiles compile on-source and we find that they are better suited
while developing and testing the code (if that's your cup of tea; it is ours).

In both cases, all external dependencies should be automatically dealt with if
they are installed in their default paths, as specified in the `*_ROOT` values
of [`extlib/makefile-default-paths`](extlib/makefile-default-paths). If not,
the suggested way to change them is to copy the file into
[`extlib/makefile-paths`](extlib/makefile-paths) and edit it. The file (if
present) is automatically read and the values found there replace the
corresponding non-default definitions. The rationale for not changing
makefile-default-paths is that makefile-paths file is .gitignore-d. Hence, it
should not be necessary to re-change the makefiles (or stash/restore the
changes) each time the project is pulled, or manually ignore the changes when
it is pushed, which is very convenient for anyone who actually develops
`MCFClass` components (anyone there?). However, note that the reading of both
[`extlib/makefile-default-paths`](extlib/makefile-default-paths) and
[`extlib/makefile-paths`](extlib/makefile-paths) can be disabled; see the
`MCFClass_READ_PATHS` option in the CMake section and the `MCFC_NO_PATHS` macro
in the makefile section below.


### Using CMake

Configure and build the library with:

```sh
mkdir build
cd build
cmake ..
make
```

If CPLEX is not available or you are not interested in `MCFCplex`, run

```sh
cmake -DMCFClass_USE_CPLEX=OFF ..
```

- Optionally, you can install the library with:

```sh
sudo make install
```

- After the library is built, you can use it in your CMake project with:

```cmake
find_package(MCFClass)
target_link_libraries(<my_target> MCFClass::MCFClass)
```

- If you use the `MCFClass` project inside some other project that already
  properly defines the `*_ROOT` values, you can avoid them being read by
  setting the option `MCFClass_READ_PATHS` to `OFF`, e.g., by adding

```cmake
set(MCFClass_READ_PATHS OFF CACHE BOOL
		                "Whether MCFClass will read locations for
		                dependencies or not." FORCE)
```

  in your CMake project.


### Using makefiles

- To create the library, go into `lib/` and type

```sh
make -f makefile-lib
```

- To test the library, go into `test` and type `make`.

- If you want to use the MCFCplex class, which comes excluded by default,
  uncomment the two lines in `lib/makefile`:

```makefile
MCFCxDIR = $(libMCFClDIR)/MCFCplex
include $(MCFCxDIR)/makefile
```

You can similarly enable (or disable) any solver, both the LGPL ones and
those under the academic license, if you have obtained them, by commenting
out (or commenting) the corresponding two lines in `lib/makefile`.

- If you want to use `MCFClass` as a part of some larger project you can
  just include [`lib/makefile-c`](lib/makefile-c) or
  [`lib/makefile-inc`](lib/makefile-inc) in the "main" makefile, provided that
  you have properly defined all the necessary input macros; see, e.g.,
  [`test/makefile`](test/makefile) for an example.

- If you use the `MCFClass` project inside some other project that already
  properly defines the `*_ROOT` values, you can avoid them being read by
  defining the macro `MCFC_NO_PATHS` in your "main" makefiles prior to
  including [`lib/makefile-c`](lib/makefile-c) or
  [`lib/makefile-inc`](lib/makefile-inc)


## Other stuff

More information about (some of) the implemented algorithms can be found at

  http://pages.di.unipi.it/frangio/abstracts.html#JOC06

A further solver, MCFIntPnt (based on Interior-Point algorithms) has been
developed, but it has not yet reached a sufficient maturuty to be
distributed; its principles are discussed at

  http://pages.di.unipi.it/frangio/abstracts.html#SIOPT04
  http://pages.di.unipi.it/frangio/abstracts.html#COAP06
  http://pages.di.unipi.it/frangio/abstracts.html#OMS06
