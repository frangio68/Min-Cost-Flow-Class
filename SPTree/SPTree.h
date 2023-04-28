/*--------------------------------------------------------------------------*/
/*----------------------------- File SPTree.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Definition of SPTree, a class deriving from MCFClass, and therefore
 * conforming to the standard MCF interface defined therein, and implementing
 * several "classic" Shortest Path Tree algorithms to solve uncapacitated
 * single-source Min Cost Flow problems. The actual algorithm can be chosen
 * at compile time by setting a proper switch.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni.
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef _SPTree
 #define _SPTree  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MCFClass.h"

/*--------------------------------------------------------------------------*/
/*---------------------------- MACROS --------------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup SPTREE_MACROS Compile-time switches in SPTree.h
 *  These macros control some important details of the implementation.
 *  Although using macros for activating features of the implementation is
 *  not very C++, switching off some unused features may make the code
 *  more efficient in running time or memory.
 *  @{ */

/*------------------------------ SPT_ALGRTM --------------------------------*/

#define SPT_ALGRTM 4

/** This macro decides which SPT algorithm has to be used.
   Possible values are:

   - 0  =>  LQueue
   - 1  =>  LDeque
   - 2  =>  (currently unused)
   - 3  =>  Dijkstra
   - 4  =>  Heap

   For algorithms based on priority lists, the macro LABEL_SETTING [see
   below] can be set to 1 to say that the algorithm is of the
   "label-setting" (nodes only exit from Q once) rather than of the
   "label-correcting" (nodes may exit from Q more than once) type. */

#if( SPT_ALGRTM <= 2 )
 #define LABEL_SETTING 0
 ///< this is a label-correcting SPT algorithm
#else
 #define LABEL_SETTING 1

 /**< This macro decides if the "label-setting" style is used.

    With a priority lists, the SPT algorithm applied to SPT problems
    with *all nonnegative arc costs* has the "label-setting" property:
    nodes only exit from Q once, hence when a node exits from Q its
    label is permanently set.

    If LABEL_SETTING > 0 the code will assume that this property holds
    and implement some things accordingly; in particular, the algorithm
    is terminated when the last destination is extracted from Q even
    though Q is still nonempty.

    \warning Solving a SPT algorithm with negative arc costs with
             LABEL_SETTING > 0 may produce a suboptimal solution. */

 #if( SPT_ALGRTM == 4 )
  #define HeapCard 2

  /**< Number of sons of each node in the heap.
     SPT_ALGRTM == 4 means using a C-ary heap to hold the node set Q: each
     HeapCard is the ariety of the heap, i.e. the max number of sons of a
     node in the heap. Special treatment is deserved to the case
     HeapCard == 2. */
 #endif
#endif

/*------------------------------ ORDRD_NMS ---------------------------------*/

#define ORDRD_NMS 1

/**< Decides if arc names in MCFGetX() are ordered.
   If ORDRD_NMS > 0, and MCFGetX() [see below] is asked for a "sparse" flow
   solution (i.e., nms != 0), then the set of indices returned at the end
   of the method is ordered in increasing sense. If ORDRD_NMS == 0 instead,
   the set of indices may not be ordered.

   ORDRD_NMS > 0 may be useful for some applications, but it is more costly
   (basically, it requires either to compute the "dense" flow solution or to
   sort a vector). Also, "sparse" flow solutions in this class are guaranteed
   to contain no more than n - 1 nonzeroes, hence if ORDRD_NMS == 0 then the
   parameter `F' in MCFGetX( F , nms ) can actually point to a (n - 1)-vector,
   while if ORDRD_NMS > 0 it must point to a m-vector anyway. */

/*----------------------------- DYNMC_MCF_SPT ------------------------------*/

#define DYNMC_MCF_SPT 0

/**< Decides if the graph topology (arcs, nodes) can be changed.
   If DYNMC_MCF_SPT > 0, some of the methods of the public interface of
   class that allow to change the topology of the underlying network are
   actually implemented. Possible values of this macro are:

   - 0 => the topology of the graph cannot be changed;

   - 1 => all the methods that change the topology of the graph are
          implemented. */

/** @} ---------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace MCFClass_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASSES -------------------------------------*/
/*--------------------------------------------------------------------------*/
/** The SPTree class derives from the abstract base class MCFClass, thus
 *  sharing its (standard) interface, and implements Shortest Path Tree
 *  algorithms for solving "uncapacitated" (Linear) Min Cost Flow
 *  problems with one source node.
 *
 *  \warning The SPT algorithm will enter in an infinite loop if a directed
 *           cycle of negative cost exists in the graph: there is no check
 *	     about this in the code. */

class SPTree : public MCFClass
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  The following methods and data are the actual interface of the      --*/
/*--  class: the standard user should use these methods and data only.    --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

   SPTree( Index nmx = 0 , Index mmx = 0 , bool Drctd = true );

/**< Constructor of the class.

   For the meaning of nmx and mmx see MCFClass::MCFClass().

   The parameter `Drctd' tells if the given graph has really to be
   understood as directed (default), i.e., if the i-th arc is
   Sn[ i ] --> En[ i ], or undirected, i.e., the i-th arc is
   Sn[ i ] <--> En[ i ]. Undirected graphs are internally implemented by
   doubling each arc, but this is completely hidden by the interface. */

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

   void LoadNet( Index nmx = 0 , Index mmx = 0 , Index pn = 0 , Index pm = 0 ,
		 cFRow pU = 0 , cCRow pC = 0 , cFRow pDfct = 0 ,
		 cIndex_Set pSn = 0 , cIndex_Set pEn = 0 ) override;

/**< Inputs a new network, as in MCFClass::LoadNet().

   Arcs with pC[ i ] == Inf< CNumber >() do not "exist". If
   DYNMC_MCF_SPT > 0, these arcs are "closed".

   If DYNMC_MCF_SPT == 0, these arcs are just removed from the formulation.
   However, they have some sort of a "special status" (after all, if the
   user wants to remove them completely he/she can just change the data), in
   that they are still counted into the number of arcs of the graph and they
   will always have 0 flow and Inf< CNumber >() reduced cost as "closed" or
   "deleted" arcs. */

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

   void SolveMCF( void ) override;

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

   using MCFClass::MCFGetX;  // the ( void ) method, which is otherwise hidden

   void MCFGetX( FRow F , Index_Set nms = 0  ,
		 Index strt = 0 , Index stp = Inf< Index >() ) const override;

/*--------------------------------------------------------------------------*/

   using MCFClass::MCFGetRC;  // the ( void ) method, which is otherwise hidden

   void MCFGetRC( CRow CR , cIndex_Set nms = 0  ,
		  Index strt = 0 , Index stp = Inf< Index >() ) const override;

   CNumber MCFGetRC( Index i ) const override;

/*--------------------------------------------------------------------------*/

   void MCFGetPi( CRow P , cIndex_Set nms = 0  ,
		  Index strt = 0 , Index stp = Inf< Index >() ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** Same meaning as MCFClass::MCFGetPi().

   \note Some of the potentials may be + Inf< CNumber >(): this means that

   - the node is *not* a destination and it cannot be reached from the Origin
     (however, this does *not* mean that the problem is unfeasible);

   - if LABEL_SETTING == 1, the node is *not* a destination and it has not
     been reached during the algorithm. */

   cCRow MCFGetPi( void ) const override { return( Pi + 1 ); }

/*--------------------------------------------------------------------------*/
/** Same meaning as MCFClass::MCFGetFO().

   \note if not all the specified destinations can be reached from the
         Origin, returns Inf< FONumber >(). */

   FONumber MCFGetFO( void ) const override { return( FO ); }

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

   void MCFArcs( Index_Set Startv , Index_Set Endv , cIndex_Set nms = 0  ,
		 Index strt = 0 , Index stp = Inf< Index >() ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   Index MCFSNde( Index i ) const override {
    if( DirSPT ) 
     return( Startn[ i ] );
    else
     return( FS[ DictM1[ 2 * i + 1 ] ].Nde );
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   Index MCFENde( Index i ) const override {
    if( DirSPT ) 
     return( FS[ DictM1[ i ] ].Nde - USENAME0 );
    else
     return( FS[ DictM1[ 2 * i ] ].Nde - USENAME0 );
    }

/*--------------------------------------------------------------------------*/

   void MCFCosts( CRow Costv , cIndex_Set nms = 0  ,
		  Index strt = 0 , Index stp = Inf< Index >() )
    const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   CNumber MCFCost( Index i ) const override {
    if( DirSPT ) 
     return( FS[ DictM1[ i ] ].Cst );
    else
     return( FS[ DictM1[ 2 * i ] ].Cst );
    }

/*--------------------------------------------------------------------------*/

   void MCFUCaps( FRow UCapv , cIndex_Set nms = 0  ,
		  Index strt = 0 , Index stp = Inf< Index >() )
    const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   FNumber MCFUCap( Index i ) const override {
    #if( ! DYNMC_MCF_SPT )
     cIndex pos = DictM1[ i ];
     if( pos >= StrtFS[ n + 1 ] )
       return( 0 );
      else
    #endif
       return( - B[ Origin ] );
    }

/*--------------------------------------------------------------------------*/

   void MCFDfcts( FRow Dfctv , cIndex_Set nms = 0  ,
		  Index strt = 0 , Index stp = Inf< Index >() )
    const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   FNumber MCFDfct( Index i ) const override { return( B[ i + 1 ] ); }

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/*----- Changing the costs, deficits and upper capacities of the (MCF) -----*/
/*--------------------------------------------------------------------------*/

   void ChgCosts( cCRow NCost , cIndex_Set nms = 0  ,
		  Index strt = 0 , Index stp = Inf< Index >() ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   void ChgCost( Index arc , CNumber NCost ) override;

/*--------------------------------------------------------------------------*/

   void ChgDfcts( cFRow NDfct , cIndex_Set nms = 0  ,
		  Index strt = 0 , Index stp = Inf< Index >() ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   void ChgDfct( Index nod , FNumber NDfct ) override;

/*--------------------------------------------------------------------------*/

   void ChgUCaps( cFRow NCap , cIndex_Set nms = 0  ,
		  Index strt = 0 , Index stp = Inf< Index >() ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   void ChgUCap( Index arc , FNumber NCap ) override;

/*--------------------------------------------------------------------------*/
/*--------------- Modifying the structure of the graph ---------------------*/
/*--------------------------------------------------------------------------*/

  void CloseArc( Index name ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  bool IsClosedArc( Index name ) const override {
   #if( DYNMC_MCF_SPT )
    cIndex pos = DictM1[ name ];  // current position of arc name
    Index nde = Startn[ name ];   // start node of arc name
    return( pos < StrtFS[ nde ] + LenFS[ nde ] );
   #else
    return( false );
   #endif
   }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  void DelNode( Index name ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  void OpenArc( Index name ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  Index AddNode( FNumber aDfct ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  void ChangeArc( Index name , Index nSS = Inf< Index >() ,
		  Index nEN = Inf< Index >() ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  void DelArc( Index name ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  bool IsDeletedArc( Index name ) const override {
   return( SPTree::IsClosedArc( name ) );
   }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  Index AddArc(cIndex Start , Index End , FNumber aU , CNumber aC ) override; 

/*--------------------------------------------------------------------------*/
/*------------------------ SPECIALIZED INTERFACE ---------------------------*/
/*--------------------------------------------------------------------------*/

   void ShortestPathTree( void );

/**< Solver of the Shortest Path Tree Problem from the current Origin.
   (specified in the constructor or by SetOrigin(), see below)

   If LABEL_SETTING == 0, or if no Destination is speficied (Dst ==
   Inf< Index >() in SetDest() [see below]), the whole Shortest Path Tree (at
   least, the SPT of the component of the graph connected with Origin) is
   computed, otherwise the code stops as soon as the shortest path between
   Origin and Dest is computed.

   Note that methods such as MCFGetX(), MCFGetRC() and MCFGetFO() may need
   some complicate calculations in order to put the solution of the Shortest
   Path in the correct format; since these calculations change some of the
   internal data structures, it is not permitted to call again
   ShortestPathTree() after that any of these methods have been called. */

/*--------------------------------------------------------------------------*/
/// Changes the Origin from which Shortest Paths are computed.

   void SetOrigin( Index NewOrg ) {
    if( Origin != NewOrg + USENAME0 ) {
     Origin = NewOrg + USENAME0;
     status = MCFClass::kUnSolved;
     }
    }

/*--------------------------------------------------------------------------*/
/** Changes the Destination node of Shotest Paths. If LABEL_SETTING == 0, it
   has no influence since label correcting methods cannot stop before the
   whole SPT has been computed. Conversely, label setting algorithms can solve
   Origin-Dest Shortest Path Problems; therefore, it is possible to obtain
   shortest paths between Origin and a subset of the nodes, by calling
   ShortestPathTree() with one of the destinations, and controlling upon
   completion that all the desidered nodes have been visited (see Reached()
   below). If this is not the case, ShortestPathTree() can be invoked again
   with one of the unreached nodes, until they are all visited.

   If no Dest is given, or if Dest is set to Inf< Index >(), the whole Shortest
   Path Tree (at least, the SPT of the component of the graph connected with
   Origin) is computed. */

   void SetDest( Index NewDst ) {
    if( Dest != NewDst + USENAME0 ) {
     #if( LABEL_SETTING )
      Dest = NewDst + USENAME0;
     #endif
     status = MCFClass::kUnSolved;
     }
    }

/*--------------------------------------------------------------------------*/

   void MCFGetX( Index ND , cIndex_Set DB , FRow F , Index_Set nms = 0 ,
		 Index strt = 0 , Index stp = Inf< Index >() ) const;

/**< Like SPTree::MCFGetX( FRow , Index_Set , cIndex , Index ), except that
   the primal solution that is returned is relative only to the subset of
   destinations whose names are contained in the first ND entries of the
   vector DB.

   Note: node names in ND must be in 1 ... n irrespective of USENAME0. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   FONumber MCFGetFO( Index ND , cIndex_Set DB ) const;

/**< Like SPTree::MCFGetFO( void ), except that the cost that is returned is
   that of the primal solution relative only to the subset of destinations
   whose names are contained in the first ND entries of the vector DB.

   Note: node names in ND must be in 1 ... n irrespective of USENAME0. */

/*--------------------------------------------------------------------------*/
/** Returns true if a shortest path from Origin to i have already been
   computed; this can be used when LABEL_SETTING == 1 to determine if a
   shortest from Origin to i have been obtained as a by-product of the
   calculation of the shortest path between Origin and some other Dest. */

   bool Reached( Index i ) const {
    return( ( Pi[ i ] < Inf< CNumber >() ) && ( Q[ i ] == Inf< Index >() ) );
    }

/*--------------------------------------------------------------------------*/
/** Return a cIndex* vector p[] such that p[ i ] is the predecessor of node
   i in the shortest path tree. If a node i has no predecessor, i.e.,
   i == Origin, i does not belong to the connected component of the origin or
   the computation have been stopped before reaching i, then p[ i ] == 0.

   \note if the name "0" is used for nodes, (USENAME0 == 1) then node names
         are internally "translated" of +1 to avoid it being used - the
	 the names reported in this vector will follow the same rule.

   For this reason, the first entry of p (*p) is not significative. */

   cIndex_Set Predecessors( void ) const { return( NdePrd ); }

/*--------------------------------------------------------------------------*/
/** Return a cIndex* vector a[] such that a[ i ] is the index of the arc
   ( p[ i ] , i ), being p[] the vector returned by the above method, and
   with the same structure. If p[ i ] == 0, then a[ i ] is not significative:
   for the Origin (that has p[ Origin ] == 0), however, it is guaranteed that
   a[ Origin ] == Inf< Index >(). */

   cIndex_Set ArcPredecessors( void );

/*--------------------------------------------------------------------------*/
/// returns the root of the SPT problem

   Index Orig( void ) const { return( Origin ); }

/*--------------------------------------------------------------------------*/
/// returns the number of destination nodes in the SPT problem

   Index DestN( void ) const { return( NDsts ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/** Returns the DestN()-vector containig the names of destination nodes in
   the SPT problem; the names are in increasing order and INF-terminated. */

   cIndex_Set Dests( void ) const { return( DstBse ); }

/*--------------------------------------------------------------------------*/
/// returns the size of the Forward Star of node i

   Index LenFS( Index i ) const {
    #if( DYNMC_MCF_SPT )
     return( LenFS[ i ] );
    #else
     return( StrtFS[ i + 1 ] - StrtFS[ i ] );
    #endif
    }

/*--------------------------------------------------------------------------*/
/// returns the h-th arc in FS( i ) for h = 0, ... , LenFS( i ) - 1

   Index ReadFS( Index i , Index h ) const {
    return( Dict[ StrtFS[ i ] + h ] );
    }

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

   virtual ~SPTree();

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  The standard user should not care about the following part: users   --*/
/*--  who need to extend the code by deriving a new class may use these   --*/
/*--  methods and data structures. It is *dangerous* to *modify* the      --*/
/*--  data structures, while it safe to read them.                        --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/

 struct FSElmnt {               // one entry of the Forward Star
                  CNumber Cst;  // cost of the arc
                  Index   Nde;  // end node of the arc
	          };

 typedef FSElmnt *FrwdStr;

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED DATA STRUCTURES --------------------------*/
/*--------------------------------------------------------------------------*/

 Index Origin;       // the source
 Index Dest;         // the sink
 Index NDsts;        // total number of destinations;
 Index_Set DstBse;   // array of indices of the destinations

 Index_Set NdePrd;   // NdePrd[ i ] = predecessor of i in the shortest path
                     // NdePrd[ Origin ] = 0

 Index_Set ArcPrd;   // ArcPrd[ i ] = index of arc ( NdePrd[ i ] , i )
                     // ArcPrd[ Origin ] = 0

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
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

   void Initialize( void );

/* Initialize the data structures for a "cold start". */

/*--------------------------------------------------------------------------*/

   void ScanFS( cIndex mi );

/* Scans the Forward Star of mi and puts in Q those nodes whose distance
   label can be decreased by using an arc emanating from mi. */

/*--------------------------------------------------------------------------*/

   Index ExtractQ( void );

/* Extracts an element (depending on the particular algoritm) from the set Q:
   if Q is empty, returns 0. */

/*--------------------------------------------------------------------------*/

#if( ( SPT_ALGRTM == 0 ) || ( SPT_ALGRTM == 3 ) )

   void InsertQ( cIndex j );

#else

   void InsertQ( cIndex j , cCNumber label );

#endif

/* Inserts the node with name j and label label somewhere in Q: the label is
   not needed for LQueue and Djkstra algorithms. */

/*--------------------------------------------------------------------------*/

   void CalcArcP( void );

/* Calculates the ArcPrd[] vector. */

/*--------------------------------------------------------------------------*/

   void MemAlloc( void );

   void MemDeAlloc( void );

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

 CRow Pi;            // node Potentials
 FRow B;             // node deficits vector

 FONumber FO;        // Objective Function value

 bool ReadyArcP;     // if the "arc predecessor" data structure has already
                     // been updated after a "final" ShortestPathTree() call 
 FrwdStr FS;         // the Forward Star (itself)

 Index_Set Q;        // the set of scanned nodes: Q[ i ] = INF ==> i \notin Q
 #if( SPT_ALGRTM <= 3 )
                     // here, Q is an array pointer implementation of a list,
                     // and *Q is the head of the list (node names are >= 1)
 #else
  Index_Set H;       // here, Q[ i ] tells the position of node i in the
		     // vector implementing the heap, and H is that vector
 #endif

 Index cFS;          // cardinality of the FS (m if DirSPT, 2m otherwise)
 Index tail;         // the tail element of the list, or the first free
                     // position in the heap

 Index_Set Startn;
 Index_Set StrtFS;   // position in FS[] where FS[ i ] begins
 #if( DYNMC_MCF_SPT )
  Index_Set LenFS;   // how many arcs there are in FS[ i ]
 #endif

 Index_Set Dict;     // arc dictionary: for each position in FS[], tells
                     // which arc is that one
 Index_Set DictM1;   // inverse of Dict: for each arc, tells where it stands
                     // in FS[] - if the graph is undirected, the two
		     // consecutive entries 2 * i and 2 * i + 1 tells the
		     // two positions of arc i in FS[]
 bool DirSPT;        // true if the graph is directed

/*--------------------------------------------------------------------------*/

 };  // end( class SPTree )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

};  // end( namespace MCFClass_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* SPTree.h included */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File SPTree.h -----------------------------*/
/*--------------------------------------------------------------------------*/
