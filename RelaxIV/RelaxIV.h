/*--------------------------------------------------------------------------*/
/*------------------------- File RelaxIV.h ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Linear Min Cost Flow problems solver, based on the RELAXIV code by
 * D. Bertsekas and P. Tseng, as described in
 *
 *    Bertsekas, Dimitri P., and Paul Tseng.
 *    "RELAX-IV: A faster version of the RELAX code for solving minimum
 *     cost flow problems." (1994), Report LIDS-P-2276, MIT.
 *
 * Conforms to the standard (MCF) interface defined in MCFClass.h.
 *
 * RelaxIV is based on a primal-dual algorithm which essentially operates
 * as follows: a pseudoflow (a flow vector which satisfies bound and
 * non-negativity constraints but not necessarily flow conservation
 * constraints) is kept which satisfies complementarity slackness conditions
 * with the current vector of potentials; that is, only the flow on arcs
 * whose reduced cost
 * \f[
 *  RC[ i , j ] =  C[ i , j ] - Pi[ j ] + Pi[ i ]
 * \f]
 * is zero can be chosen to any value between 0 and the capacity, while arcs
 * with reduced cost < 0 are saturated (fixed to their capacity) and arcs
 * with reduced cost > 0 are empty (fixed to 0).
 *
 * The algorithm attempts to convert the pseudoflow into a flow (i.e.,
 * to satisfy the flow conservation constraints) by essentially running
 * a max-flow algorithm of the augmenting path type. If the flow is found then
 * this is an optimal solution of the problem and the algorithm is stopped.
 * Otherwise, a saturated cut is identified which separates the origins (nodes
 * not yet producing enough flow) to the destinations (nodes not yet consuming
 * enough flow); this cut is used to modify the potentials, thereby creating
 * new arcs with zero reduced cost, which can be used to push further flow
 * from the origins to the destinations. If no such arcs can be created the
 * problem is declared unfeasible. Much care is devoted to stop the max-flow
 * computation as soon as a proof that the set of potentials is not optimal,
 * in order to reach as soon as possible a dual optimal solution, and to
 * re-use all available information to "warm start" the max-flow computation
 * after a change in the potentials.
 *
 * \warning The original code has been written for integer data only.
 *          By properly setting the flow and cost tolerances [see
 *          SetEps****() in MCFClass.h] we have always been able to solve
 *          any MCF that we could throw at the solver, but in principle this
 *          kind of algorithm may fail to converge with nonintegral data, so
 *          consider yourselves warned.
 *
 * \warning A private type SIndex is defined which is intended to hold arc
 *          and node indices "with a sign", used to represent orientation.
 *          This has to be "in sync" with Index, in the sense that for every
 *          unsigned index value in Index, the two signed values should be
 *          feasible in SIndex. In other words, either Index is not using at
 *          least half of its feasible values, or SIndex has to be a "bigger"
 *          data type than Index. The default value for SIndex is int.
 *
 * \author <b>(original FORTRAN code)</b> \n
 *         Dimitri P. Bertsekas \n
 *         Lab. for Information and Decision Systems \n
 *         Massachusetts Institute of Technology \n
 *
 * \author <b>(original FORTRAN code)</b> \n
 *         Paul Tseng \n
 *         Department of Mathematics \n
 *         University of Washington \m
 *
 * \author <b>(C++ porting and polishing)</b> \n
 *         Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author <b>(C++ porting and polishing)</b> \n
 *         Claudio Gentile \n
 *         Istituto di Analisi di Sistemi e Informatica \n
 *         Consiglio Nazionale delle Ricerche \n
 *
 * Copyright &copy by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __RelaxIV
 #define __RelaxIV  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MCFClass.h"

/*--------------------------------------------------------------------------*/
/*--------------------------------- MACROS ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup RELAXIV_MACROS Compile-time switches in RelaxIV.h
    These macros control some important details of the implementation.
    Although using macros for activating features of the implementation is
    not very C++, switching off some unused features may make the code
    more efficient in running time or memory.
    @{ */

/*----------------------------- DYNMC_MCF_RIV ------------------------------*/

#define DYNMC_MCF_RIV 3

/**< Decides if the graph topology (arcs, nodes) can be changed.
   If DYNMC_MCF_RIV > 0, the methods of the public interface of class that
   allow to change the topology of the underlying network are actually
   implemented. Possible values of this macro are:

   - 0 => the topology of the graph cannot be changed;

   - 1 => the methods that "close" arcs and delete nodes are implemented;

   - 2 => the methods that "open" previously closed arcs and add nodes are
          implemented;

   - 3 => the methods that change the start and end node of a (possibly
          "closed") arc, delete and create new arcs are implemented. */

/*-------------------------------- AUCTION ---------------------------------*/

#define AUCTION 0

/**< Decides if the auction/shortest paths inizialization procedure is used.
   If AUCTION == 1, then an auction/shortest paths inizialization procedure
   is provided [see SetPar() below] that has been reported to make the
   RelaxIV algorithm run faster on some classes of instances.
   The auction initialization essentially "spreads" the imbalances around
   the graph by performing some steps of the "pure" epsilon-relaxation
   method: this should produce "short" augmenting steps, that seem to be
   the best situation for RelaxIV.

   By setting AUCTION == 0, some memory is saved. */

/*-------------------------- RELAXIV_STATISTICS ----------------------------*/

#define RELAXIV_STATISTICS 0

/**< If RELAXIV_STATISTICS > 0, then statistic information about the behaviour
   of the Relaxation algorithm is computed. */

/** @} ---------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace MCFClass_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS RelaxIV --------------------------------*/
/*--------------------------------------------------------------------------*/
/** The RelaxIV class derives from the abstract base class MCFClass, thus
    sharing its (standard) interface, and implements a Relaxation algorithm
    for solving (Linear) Min Cost Flow problems. */

class RelaxIV : public MCFClass {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** Public enum describing the possible parameters of the MCF solver,
    "extended" from MCFClass::MCFParam, to be used with the methods
    SetPar() and GetPar(). */

  enum MCFRParam { kAuction = kLastParam     ///< crash initialization
                   };

/*--------------------------------------------------------------------------*/
/** Public enum describing the more file formats in RelaxIV::WriteMCF(). */

  enum RIVFlFrmt { kCLP = kMPS + 1 ,  ///< the "LP" format
		   kRIV               ///< RelaxIV-specific format
                   };

/*--------------------------------------------------------------------------*/

  typedef bool           *Bool_Vec;        ///< vector of booleans

/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

   RelaxIV( Index nmx = 0 , Index mmx = 0 );

/**< Constructor of the class, as in MCFClass::MCFClass(). */

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

   void LoadNet( Index nmx = 0 , Index mmx = 0 , Index pn = 0 , Index pm = 0 ,
		 cFRow pU = 0 , cCRow pC = 0 , cFRow pDfct = 0 ,
		 cIndex_Set pSn = 0 , cIndex_Set pEn = 0 ) override;

/**< Inputs a new network, as in MCFClass::LoadNet().

   Arcs with pC[ i ] == Inf< CNumber >() do not "exist". If DYNMC_MCF_RIV > 0,
   these arcs are "closed".

   If DYNMC_MCF_RIV == 0, these arcs are just removed from the formulation.
   However, they have some sort of a "special status" (after all, if the user
   wants to remove them completely he/she can just change the data), in that
   they are still counted into the number of arcs of the graph and they will
   always have 0 flow and Inf< CNumber >() reduced cost as "closed" or
   "deleted" arcs. */

/*--------------------------------------------------------------------------*/
/// set integer parameters of the algorithm
/** Set integer parameters of the algorithm.

   @param par   is the parameter to be set;

   @param value is the value to assign to the parameter.  

   Apart from the parameters of the base class, this method handles:

   - kAuction: if set to kYes, the auction/shortest paths initialization is
               used in SolveMCF() to generate the starting solution; if
	       set to kNo (default), then the default initialization based on
	       special single-node relaxation iterations is used instead.
	       Note that this parameter is *ignored* if AUCTION == 0. */

   void SetPar( int par , int val ) override
   {
    if( par == kAuction ) {
     #if( AUCTION )
      crash = ( val == kYes ) ? TRUE : FALSE;
     #else
      if( val == kYes )
       throw( MCFException( "Auction initialization not available" ) );
     #endif
     }
    else
     MCFClass::SetPar( par , val );
  }

/*--------------------------------------------------------------------------*/
// set double parameters of the algorithm
/* Set double parameters of the algorithm.
 *
 * This should in princible not be necessary, as RelaxIV has no double
 * parameters to set. However, without this being well-defined, template
 * classes having RelaxIV as template type may fail to be able to use the
 * base class method in its stead and resort to wrongly calling the
 * SetPar( , int ) version instead (no idea why), so this useless method has
 * to be kept here. */

   void SetPar( int par , double val ) override {
    MCFClass::SetPar( par , val );
    }

/*--------------------------------------------------------------------------*/
/** Returns one of the integer parameters of the algorithm.

   @param par  is the parameter to return [see SetPar( int ) for comments];

   @param val  upon return, it will contain the value of the parameter.

   Apart from the parameters of the base class, this method handles kAuction.
   */

   void GetPar( int par , int &val ) const override {
    if( par == kAuction )
     #if( AUCTION )
      val = crash ? kYes : kNo;
     #else
      val = kNo;
     #endif
    else
     MCFClass::GetPar( par , val );
    }

/*--------------------------------------------------------------------------*/
/** Returns one of the double parameters of the algorithm
 *
 * This should in princible not be necessary, as RelaxIV has no double
 * parameters to report. However, without this being well-defined, template
 * classes having RelaxIV as template type may fail to be able to use the
 * base class method in its stead (no idea why), so this useless method has
 * to be kept here. */

 void GetPar( int par , double &val ) const override {
  MCFClass::GetPar( par , val );
  }

/*--------------------------------------------------------------------------*/
/** If this method is called, a preprocessing phase is performed trying to
   reduce the arc capacities. This may sometimes help in speeding up the
   solution of the problem, but may also change the capacities returned by
   MCFUCap[s]() [see below].

   For this method to work properly, arc capacities, node deficits and the
   topology of the graph must have already been provided with LoadNet()
   [see above].

   This method can be called more than once, for instance whenever the
   capacities of some arcs or the deficits of some nodes are changed;
   however, it destroys the provious optimal solution (if any), forcing the
   algorithm to restart from scratch. */

   void PreProcess( void ) override;

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

   void SolveMCF( void ) override;

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

   void MCFGetX( FRow F , Index_Set nms = 0 ,
		 Index strt = 0 , Index stp = Inf< Index >() ) const override;

   cFRow MCFGetX( void ) const override { return( X + 1 ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   void MCFGetRC( CRow CR , cIndex_Set nms = 0 ,
		  Index strt = 0 , Index stp = Inf< Index >() )
    const override;

   cCRow MCFGetRC( void ) const override { return( RC + 1 ); }

   CNumber MCFGetRC( Index i ) const override { return( RC[ i + 1 ] ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   void MCFGetPi( CRow P , cIndex_Set nms = NULL ,
		  Index strt = 0 , Index stp = Inf< Index >() )
    const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   cCRow MCFGetPi( void ) const override { return( Pi + 1 ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   FONumber MCFGetFO( void ) const override { return( FO ); }

/*--------------------------------------------------------------------------*/

   MCFStatePtr MCFGetState( void ) const override;

/**< Same meaning as MCFClass::MCFGetState().

   The state of the algorithm is the pair S = ( X[] , RC[] ) of the arc
   flows and reduced costs. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   void MCFPutState( MCFClass::MCFStatePtr S ) override;

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

   void MCFArcs( Index_Set Startv , Index_Set Endv , cIndex_Set nms = 0 ,
		 Index strt = 0 , Index stp = Inf< Index >() ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   Index MCFSNde( Index i ) const override {
    return( Startn[ i + 1 ] - USENAME0 );
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   Index MCFENde( Index i ) const override {
    return( Endn[ i + 1 ] - USENAME0 );
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/** Same meaning as MCFClass::MCFSNdes().

   \note MCFSNdes() returns a pointers to a (read-only) vector containing
         the arc start nodes *only if USENAME0 == 0*; otherwise, it returns
	 NULL. */

   cIndex_Set MCFSNdes( void ) const override {
    #if( USENAME0 )
     return( 0 );
    #else
     return( Startn + 1 );
    #endif
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/** Same meaning as MCFClass::MCFENdes().

   \note MCFENdes() returns a pointers to a (read-only) vector containing
         the arc end nodes *only if USENAME0 == 0*; otherwise, it returns
	 NULL. */

   cIndex_Set MCFENdes( void ) const override {
    #if( USENAME0 )
     return( 0 );
    #else
     return( Endn + 1 );
    #endif
    }

/*--------------------------------------------------------------------------*/

   void MCFCosts( CRow Costv , cIndex_Set nms = 0  ,
		  Index strt = 0 , Index stp = Inf< Index >() )
    const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   CNumber MCFCost( Index i ) const override { return( C[ i + 1 ] ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   cCRow MCFCosts( void ) const override { return( C + 1 ); }

/*--------------------------------------------------------------------------*/

   void MCFUCaps( FRow UCapv , cIndex_Set nms = 0 ,
		  Index strt = 0 , Index stp = Inf< Index >() )
    const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   FNumber MCFUCap( Index i ) const override { return( Cap[ i + 1 ] ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   cFRow MCFUCaps( void ) const override { return( Cap + 1 ); }

/*--------------------------------------------------------------------------*/

   void MCFDfcts( FRow Dfctv , cIndex_Set nms = 0  ,
		  Index strt = 0 , Index stp = Inf< Index >() )
    const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   FNumber MCFDfct( Index i ) const override { return( B[ i + 1 ] ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   cFRow MCFDfcts( void ) const override { return( B + 1 ); }

/*--------------------------------------------------------------------------*/

   void WriteMCF( ostream &oStrm , int frmt = 0 ) const override;

/**< Extends MCFClass::WriteMCF() to support two new formats:

   - kCLP is the "LP" format read by several LP solvers;

   - kRIV is the following RelaxIV-specific format:

          - < number of nodes > < number of arcs >

	  - for( < each arc > )
	    < start node > < end node > < reduced_capacity > < reduced_cost >

	  - for( < each node > )
	    < reduced flow deficit at node >

	  \note the data of the problem in this format is not that of the
	        original problem, but rather that of the "reduced" problem
	        corresponding to the current pair (flow, potential) of the
		relaxation algorithm. */

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/*----- Changing the costs, deficits and upper capacities of the (MCF) -----*/
/*--------------------------------------------------------------------------*/

   void ChgCosts( cCRow NCost , cIndex_Set nms = 0 ,
		  Index strt = 0 , Index stp = Inf< Index >() ) override;

   void ChgCost( Index arc , CNumber NCost ) override;

/*--------------------------------------------------------------------------*/

   void ChgDfcts( cFRow NDfct , cIndex_Set nms = 0 ,
		  Index strt = 0 , Index stp = Inf< Index >() ) override;

   void ChgDfct( Index nod , FNumber NDfct ) override;

/*--------------------------------------------------------------------------*/

   void ChgUCaps( cFRow NCap , cIndex_Set nms = 0 ,
		  Index strt = 0 , Index stp = Inf< Index >() ) override;

   void ChgUCap( Index arc , FNumber NCap  ) override;

/*--------------------------------------------------------------------------*/
/*--------------- Modifying the structure of the graph ---------------------*/
/*--------------------------------------------------------------------------*/

   void CloseArc( Index name ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   bool IsClosedArc( Index name ) const override {
    #if( DYNMC_MCF_RIV > 2 )
     return( ( RC[ name + 1 ] == Inf< CNumber >() ) &&
	     ( Startn[ name + 1 ] < Inf< Index >() ) );
    #elif( DYNMC_MCF_RIV )
     return( RC[ name + 1 ] == Inf< CNumber >() );
    #else
     return( false );
    #endif
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   void DelNode( Index name ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   void OpenArc( Index name ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   Index AddNode( FNumber aDfct ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   void ChangeArc( Index name , Index nSS = Inf< Index >() ,
		   Index nEN = Inf< Index >() ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   void DelArc( Index name ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   bool IsDeletedArc( Index name ) const override {
    #if( DYNMC_MCF_RIV > 2 )
     return( Startn[ name + 1 ] == Inf< Index >() );
    #else
     return( false );
    #endif
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   Index AddArc( Index Start , Index End , FNumber aU , CNumber aC ) override;

/*--------------------------------------------------------------------------*/
/*------------------------ SPECIALIZED INTERFACE ---------------------------*/
/*--------------------------------------------------------------------------*/
   /// total number of (single-node or multinode) iterations

   int MCFiter( void ) const { return( iter ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
   /// number of flow augmentations

   int MCFaug( void ) const { return( num_augm ); }

/*--------------------------------------------------------------------------*/

#if( RELAXIV_STATISTICS )
   /// number of multinode iterations

   int RelaxIV::MCFmulti( void ) const { return( nmultinode ); }
 
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
   /// number of dual ascent steps

   int RelaxIV::MCFascnt( void ) const { return( num_ascnt ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 #if( AUCTION )
   /// number of iterations in the Auction() initialization

   int RelaxIV::MCFauct( void ) const { return( nsp ); }
 #endif
#endif

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

  virtual ~RelaxIV();

/*--------------------------------------------------------------------------*/
/*------------------------ PUBLIC DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*-- Nobody should ever look at this part: everything that is under this  --*/
/*-- advice may be changed without notice in any new release of the code. --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE TYPES -------------------------------*/
/*--------------------------------------------------------------------------*/

  typedef int		  SIndex;            ///< an index with a sign
  typedef SIndex          *SIndex_Set;       ///< set (array) of SIndex
  typedef const SIndex    cSIndex;           ///< a read-only SIndex
  typedef cSIndex        *cSIndex_Set;       ///< read-only SIndex array

/*--------------------------------------------------------------------------*/

   class RIVState : public MCFClass::MCFState {
    public:

     RIVState( Index m );
     ~RIVState();

     FRow Flow;
     CRow RedCost;
     };

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------- called in SolveMCF() ---------------------------*/
/*--------------------------------------------------------------------------*/

   void init_tree( void );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   void init_standard( void );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   FNumber svblncdarcs( Index node , cIndex_Set tfst1 , cIndex_Set tnxt1 ,
			             cIndex_Set tfst2 , cIndex_Set tnxt2 );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   FNumber dascnt( Index node , CNumber &delprc , cIndex_Set F1 ,
		   cIndex_Set Nxt1 , cIndex_Set F2 , cIndex_Set Nxt2 );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   void relist( Index node );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   void AugFlow( Index augnod , Index root , Index node_p , Index node_n ,
		 cIndex_Set Term1 , cIndex_Set Term2 );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   bool Ascnt( FNumber sdm , FNumber delx , Index &nlabel ,
	       bool &Switch , Index &nscan , Index &curnode ,
	       cIndex_Set Term1 , cIndex_Set Term2 ,  cIndex_Set F1 ,
	       cIndex_Set Nxt1 , cIndex_Set F2 , cIndex_Set Nxt2 );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 #if( AUCTION )
   void auction( void );
 #endif

/*--------------------------------------------------------------------------*/
/*----------------------- called in init_standard --------------------------*/
/*--------------------------------------------------------------------------*/

   CNumber nxtbrkpt( cIndex_Set t_St1 , cIndex_Set NSt1 ,
		     cIndex_Set t_St2 , cIndex_Set NSt2 );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   CNumber mvflw1( Index arc , FRow tDfct , FRow tDDNeg ,
		   cIndex_Set Term , FRow Flow1 , FRow Flow2 );

   CNumber mvflw2( Index arc , FRow tDfct , FRow tDDPos ,
		   cIndex_Set Term , FRow Flow1 , FRow Flow2 );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   void decrsRC( Index arc , CNumber trc , CNumber delprc , CNumber &nxtbrk ,
		 FRow tDD1 , FRow DD2 , cIndex_Set Term );

   void incrsRC( Index arc , CNumber trc , CNumber delprc , CNumber &nxtbrk ,
		 FRow tDD1 , FRow DD2 , cIndex_Set Term );

/*--------------------------------------------------------------------------*/
/*--------------------------- called in Chg**** ----------------------------*/
/*--------------------------------------------------------------------------*/

   void chgcsti( Index i , CNumber NCost );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   void chgcapi( Index i , FNumber NCap );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

#if( DYNMC_MCF_RIV )

   void delarci( Index arc );

 #if( DYNMC_MCF_RIV > 1 )

   void addarci( Index arc );

 #endif
#endif

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   void cmptprices( void );

/*--------------------------------------------------------------------------*/

   void MemAlloc( void );

   void MemDeAlloc( void );

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

 FRow X;                  ///< arc Flows
 FRow U;                  ///< arc residual capacities
 FRow Cap;                ///< arc Capacities

 CRow C;                  ///< arc Costs
 CRow RC;                 ///< arc Reduced Costs

 FRow B;                  ///< node deficits vector
 FRow Dfct;               ///< node residual deficits

 FONumber FO;             ///< Objective Function value

 Index_Set tfstou;        ///< first forward balanced arc
 Index_Set tnxtou;        ///< next forward balanced arc
 Index_Set tfstin;        ///< first backward balanced arc
 Index_Set tnxtin;        ///< next backward balanced arc

 Index nb_pos;            ///< number of "directed" balanced arcs
 Index nb_neg;            ///< number of "inverse" balanced arcs

 #if( DYNMC_MCF_RIV > 2 )
 Index ffp;  ///< first free arc name, InINF if none
	     /**< ffp, if not-InINF, is the head of a queue of available
	      * arc names implemented in Endn[]. That is, Endn[ ffp ] is
	      * the next available name, Endn[ Endn[ ffp ] ] is the one
	      * after, and so on. The queue is kept ordered by arc name,
	      * which requires O( number of deleted arcs ) in DelArc() but
	      * O( 1 ) in AddArc(), and it is InINF-terminated. */
 #endif

 #if( AUCTION )
  bool crash;             /**< true => initialization is perfomed by the
			   * auction routine, false => it is performed by
			   * single node relaxation iterations */
 #endif

 int iter;                ///< number of iterations (of both types)
 int num_augm;            ///< number of flow augmentation steps
 #if( RELAXIV_STATISTICS )
  int nmultinode;         ///< number of multinode iterations
  int num_ascnt;          ///< number of multinode ascent steps
  #if( AUCTION )
   int nsp;               ///< n. of auction/shortest path iterations
  #endif
 #endif

 Index error_node;  ///< node where unfeasibility/unboundednedd is detected

 Index error_info;  /**< 1 unfeasibility detected in PreProcessing(): out
                     *     capacity of node error_node < - deficit
                     *   2 unfeasibility detected in PreProcessing(): in
                     *     capacity of node error_node < deficit
                     *   3 exit during initialization by single node
                     *     iterations: dual ascent feasible ray was found
                     *     while increasing price of node error_node;
                     *   4 exit during initialization by single node
                     *     iterations: dual ascent feasible ray was found
                     *     while decreasing price of node error_node;
                     *   5 dual ascent feasible ray was found during a
                     *     relaxation iterazion at node error_node with
                     *     positive deficit;
                     *   6 dual ascent feasible ray was found during a
                     *     relaxation iterazion at node error_node with
                     *     negative deficit;
                     *   7 dual ascent feasible ray was found during a
                     *     multinode relaxation iteration, error_node is
                     *     the starting node of the iteration;
                     *   8 problem has been detected unfeasible in
                     *     Auction() initialization. */

 CRow Pi;          ///< node Potentials

 Bool_Vec mark;      ///< temporary for multinode iterations
 Index_Set save;     ///< temporary for multinode iterations
 Index_Set label;    ///< temporary for multinode iterations
 SIndex_Set Prdcsr;  ///< temporary for multinode iterations

 Bool_Vec scan;    ///< which node belongs to S in multinode iteration
 Index_Set queue;  ///< queue of non zero deficit nodes
 Index lastq;      ///< index of the last element in the queue
 Index prvnde;     ///< index of the element preceding lastqueue

 FRow DDNeg;       ///< positive directional derivative at nodes
 FRow DDPos;       ///< negative directional derivative at nodes

 #if( AUCTION )
  CRow SB_level;          ///< temporary used in Auction()
  SIndex_Set extend_arc;  ///< temporary used in Auction()
  SIndex_Set SB_arc;      ///< temporary used in Auction()
  Index_Set FpushF;       ///< temporary used in Auction()
  Index_Set NxtpushF;     ///< temporary used in Auction()
  Index_Set FpushB;       ///< temporary used in Auction()
  Index_Set NxtpushB;     ///< temporary used in Auction()
 #endif

 Index_Set Startn;  ///< Start node of each arc
 Index_Set Endn;    ///< End node of each arc

 Index_Set FOu;     ///< first arc exiting from node
 Index_Set NxtOu;   ///< next arc exiting from Startn[ a ]
 Index_Set FIn;     ///< first arc entering into node
 Index_Set NxtIn;   ///< next arc entering into Endn[ a ]

/*--------------------------------------------------------------------------*/

 };  // end( class RelaxIV )

/*--------------------------------------------------------------------------*/

}  // end( namespace MCFClass_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* RelaxIV.h included */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File RelaxIV.h ----------------------------*/
/*--------------------------------------------------------------------------*/
