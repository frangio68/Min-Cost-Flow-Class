/*--------------------------------------------------------------------------*/
/*------------------------- File MCFCplex.h --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Quadratic Min Cost Flow problems solver, based on calls to the Cplex
 * Callable Libraries. Conforms to the standard MCF interface
 * defined in MCFClass.h
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Antonio Manca \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Matteo Sammartino \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Antonio Manca, Matteo Sammartino
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __MCFCplex
 #define __MCFCplex  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*----------------------------- INCLUDES -----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MCFClass.h"

#include <ilcplex/cplex.h>

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS -----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup MCFCPLEX_MACROS Compile-time switches in MCFCplex.h
    These macros control some important details of the implementation.
    Although using macros for activating features of the implementation is
    not very C++, switching off some unused features may make the code
    more efficient in running time or memory.
    @{ */

#define CNUMBER_IS_DOUBLE 1

/**< Tells if CNumber is in fact double.
   Although the MCFClass interface is designed to work seamlessly for every
   possibile choice of the basic types FNumber and CNumber, Cplex only works
   with doubles. Thus, when CNumber != double some conversions have to be done
   between the internal Cplex format and the CNumbers that are received in
   input/expected in output; when CNumber == double this can be saved, and
   setting CNUMBER_IS_DOUBLE == 1 heqps the code to perform some operations
   faster and using less memory. */

#define FNUMBER_IS_DOUBLE 1

/**< Tells if FNumber is in fact double.
   Although the MCFClass interface is designed to work seamlessly for every
   possibile choice of the basic types FNumber and CNumber, Cplex only works
   with doubles. Thus, when FNumber != double some conversions have to be done
   between the internal Cplex format and the FNumbers that are received in
   input/expected in output; when FNumber == double this can be saved, and
   setting FNUMBER_IS_DOUBLE == 1 heqps the code to perform some operations
   faster and using less memory. */

#define INDEX_IS_UINT 1

/**< Tells if Index is in fact [unsigned] int.
   Indices in Cplex are ints; for Index == [unsigned] int, some conversions
   can be avoided [note: this assumes that unsigned ints and ints are
   phisically the same type], and setting INDEX_IS_UINT == 1 does that. */

/*----------------------------- DYNMC_MCF_CPX ------------------------------*/

#define DYNMC_MCF_CPX 1

/**< Decides if the graph topology (arcs, nodes) can be changed.
   If DYNMC_MCF_CPX > 0, some of the methods of the public interface of
   class that allow to change the topology of the underlying network are
   actually implemented. Possible values of this macro are:

   - 0 => arcs cannot be added or deleted, closed arcs cannot be reopened;
     all the other operations are possible;

   - 1 => all the methods that change the topology of the graph are
          implemented. */

/** @} ---------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace MCFClass_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS MCFCplex --------------------------------*/
/*--------------------------------------------------------------------------*/
/** The MCFCplex class derives from the abstract base class MCFClass, thus
 *  sharing its (standard) interface, and solves (Linear) Min Cost Flow
 *  problems via calls to Cplex Callable Library functions. */

class MCFCplex: public MCFClass {

/*--------------------------------------------------------------------------*/
/*-------------------PUBLIC PART OF THE CLASS-------------------------------*/
/*--                                                                     ---*/
/*--  The following methods and data are the actual interface of the     ---*/
/*--  class: the standard user should use these methods and data only    ---*/
/*--                                                                     ---*/
/*--------------------------------------------------------------------------*/
   
 public:

/*--------------------------------------------------------------------------*/
/*------------------------- PUBLIC DATA STRUCTURES -------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/** Public enum describing the possible parameters of the MCF solver,
 *  "extended" from MCFClass::MCFParam, to be used with the methods
 *  SetPar() and GetPar(). */

  enum MCFCParam { kQPMethod = kLastParam     ///< solution method
                   };

/*--------------------------------------------------------------------------*/
/** enum describing possible ways for solving QP problem: see the Cplex
 *  manual for details */

   enum QPMethod { qpAutomatic = 0 ,
		   qpPSimplex  = 1 ,
		   qpDSimplex  = 2 ,
		   qpNSimplex  = 3 ,
		   qpBarrier   = 4
                   }; 

/*--------------------------------------------------------------------------*/
/*------------------------PUBLIC METHODS------------------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------CONSTRUCTOR--------------------------------------*/
/*--------------------------------------------------------------------------*/
/** Constructor of the class.
 *
 * For the meaning of nmx and mmx see MCFClass::MCFClass(). */

   MCFCplex( cIndex nmx = 0 , cIndex mmx = 0 );

/*--------------------------------------------------------------------------*/
/*---------------------- OTHER INITIALIZATIONS -----------------------------*/
/*--------------------------------------------------------------------------*/
/** Inputs a new network, as in MCFClass::LoadNet().
 *
 * Passing pC[ i ] == C_INF means that the arc `i' does not exist in the
 * problem. These arcs are just "closed" and their cost is set to 0: this is
 * done for being (if DYNMC_MCF_CPX > 0) subsequently capable of "opening"
 * them back with OpenArc(). If the corresponding pU[ i ] is == F_INF then
 * the arc is just "deleted". */

   void LoadNet( cIndex nmx = 0 , cIndex mmx = 0 , cIndex pn = 0 ,
		 cIndex pm = 0 , cFRow pU = NULL , cCRow pC = NULL ,
		 cFRow pDfct = NULL , cIndex_Set pSn = NULL ,
		 cIndex_Set pEn = NULL ) override;

/*--------------------------------------------------------------------------*/
/** Set integer parameters of the algorithm.

   @param par   is the parameter to be set;

   @param value is the value to assign to the parameter.  

   Apart from the parameters of the base class, this method handles:

   - kQPMethod:   the alorithm used to solve the QP, possible values are
                  defined in the enum QPMethod.

   - <any other>: any unrecognized value is taken to be one of the the many
                  "int" algorithmic parameters of Cplex and passed right
                  away via CPXsetintparam() [see the documentation in the
                  Cplex manual for details. */

   void SetPar( int par , int val ) override {
    try {
     MCFClass::SetPar( par , val );  // is it handled by the base class?

     if( par == kMaxIter )           // let MaxIter to be handled by Cplex
      CPXsetintparam( env , CPX_PARAM_NETITLIM , val > 0 ? val : 2100000000 );
     }
    catch( MCFException &e ) {      // it is *not* handled by the base class
     if( par == kQPMethod )
      QPMthd = (QPMethod) val;
     else
      CPXsetintparam( env , par , val );
     }
    }

/*--------------------------------------------------------------------------*/
/** Set float parameters of the algorithm.

   @param par   is the parameter to be set;

   @param value is the value to assign to the parameter.  

   Apart from the parameters of the base class, this method handles:

   - <any other>: any unrecognized value is taken to be one of the the many
                  "int" algorithmic parameters of Cplex and passed right
                  away via CPXsetintparam(); see the documentation in the
                  Cplex manual for details. */

   void SetPar( int par , double val ) override {
    try {
     MCFClass::SetPar( par , val );  // is it handled by the base class?

     if( par == kMaxTime )           // let MaxTime to be handled by Cplex
      CPXsetdblparam( env , CPX_PARAM_TILIM , val > 0 ? val : 1e+75 );
     }
    catch( MCFException &e ) {      // it is *not* handled by the base class
     CPXsetdblparam( env , par , val );
     }
    }

/*--------------------------------------------------------------------------*/
/** Returns one of the integer parameters of the algorithm
 *
 * @param par  is the parameter to return [see SetPar( int ) for comments];
 *
 * @param val  upon return, it will contain the value of the parameter.
 *
 * Apart from the parameters of the base class, this method handles
 * kQPMethod. */

   void GetPar( int par , int &val ) const override {
    if( par == kQPMethod )
     val = (int) QPMthd;
    else
     MCFClass::GetPar( par , val );
    }

/*--------------------------------------------------------------------------*/
/** Returns one of the double parameters of the algorithm
 *
 * This should in princible not be necessary, as MCFCplex has no double
 * parameters to report. However, without this being well-defined, template
 * classes having MCFCplex as template type may fail to be able to use the
 * base class method in its stead (no idea why), so this useless method has
 * to be kept here. */
   
   void GetPar( int par , double & val ) const override {
    MCFClass::GetPar( par, val );
    }

/*--------------------------------------------------------------------------*/
/// returns a pointer to the internal Cplex environment

   CPXENVptr GetCplexEnv( void ) const { return( env ); }

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

   void SolveMCF( void ) override;

/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR READING RESULTS  -------------------------*/
/*--------------------------------------------------------------------------*/

  using MCFClass::MCFGetX;  // the ( void ) method, which is otherwise hidden

  void MCFGetX( FRow F , Index_Set nms = 0 , Index strt = 0 ,
		Index stp = Inf< Index >() ) const override;

/*--------------------------------------------------------------------------*/

  using MCFClass::MCFGetRC;  // the ( void ) method, which is otherwise hidden

  void MCFGetRC( CRow CR , cIndex_Set nms = 0 , Index strt = 0 ,
		 Index stp = Inf< Index >() ) const override;

  CNumber MCFGetRC( Index i ) const override;

/*--------------------------------------------------------------------------*/

  using MCFClass::MCFGetPi;  // the ( void ) method, which is otherwise hidden

  void MCFGetPi( CRow P , cIndex_Set nms = 0 , Index strt = 0 ,
		 Index stp = Inf< Index >() ) const override;

/*--------------------------------------------------------------------------*/

  FONumber MCFGetFO( void ) const override;

/*--------------------------------------------------------------------------*/
/*---------- METHODS FOR READING THE DATA OF THE PROBLEM -------------------*/
/*--------------------------------------------------------------------------*/

  void MCFArcs( Index_Set Startv , Index_Set Endv , cIndex_Set nms = 0 ,
		Index strt = 0 , Index stp = Inf< Index >() ) const override;

/*--------------------------------------------------------------------------*/

  Index MCFSNde( Index i ) const override {
   int strtn;
   if( net )
    CPXNETgetarcnodes( env , net , &strtn , NULL , int( i ) , int( i ) );
   else
    strtn = Startn[ i ];

   #if( ! USENAME0 )
    strtn++;
   #endif

   return( strtn );
   }

/*--------------------------------------------------------------------------*/

  Index MCFENde( Index i ) const override {
   int endn;
   if( net )
    CPXNETgetarcnodes( env , net , NULL , &endn , int( i ) , int( i ) );
   else
    endn = Endn[ i ];

   #if( ! USENAME0 )
    endn++;
   #endif

   return( endn );
   }

/*--------------------------------------------------------------------------*/
 
  void MCFCosts( CRow Costv , cIndex_Set nms = 0 , Index strt = 0 ,
		 Index stp = Inf< Index >() ) const override;

/*--------------------------------------------------------------------------*/

  CNumber MCFCost( Index i ) const override {
   double cst;
   if( net )  
    CPXNETgetobj( env , net , &cst , int( i ) , int( i ) );
   else
    CPXgetobj( env , qp , &cst , int( i ) , int( i ) );

   return( cst );
   }

/*--------------------------------------------------------------------------*/

  void MCFQCoef( CRow Qv , cIndex_Set nms = 0 , Index strt = 0 ,
		 Index stp = Inf< Index >() ) const override;

/*--------------------------------------------------------------------------*/

  CNumber MCFQCoef( Index i ) const override {
   double qcoef = 0;
   if( qp )
    CPXgetqpcoef( env , qp , int( i ) , int( i ) , &qcoef );

   return( qcoef );
   }

/*--------------------------------------------------------------------------*/

  void MCFUCaps( FRow UCapv , cIndex_Set nms = 0 , Index strt = 0 ,
		 Index stp = Inf< Index >() ) const override;

/*--------------------------------------------------------------------------*/

  FNumber MCFUCap( Index i ) const override {
   #if( DYNMC_MCF_CPX )
    if( ( ArcPos[ i ] >= 0 ) && ( ArcPos[ i ] < Inf< FNumber >() ) )
     return( ArcPos[ i ] );
   #endif

   double ucap;
   if( net )
    CPXNETgetub( env , net , &ucap , int( i ) , int( i ) );
   else
    CPXgetub( env , qp , &ucap , int( i ) , int ( i ) );

   return( ucap );
   }

/*--------------------------------------------------------------------------*/
     
  void MCFDfcts( FRow Dfctv , cIndex_Set nms = 0 , Index strt = 0 ,
		 Index stp = Inf< Index >() ) const override;

/*--------------------------------------------------------------------------*/

  FNumber MCFDfct( Index i ) const override {
   double dfct;
   if( net )
    CPXNETgetsupply( env , net , &dfct , int( i ) , int( i ) );
   else
    CPXgetrhs( env , qp , &dfct , int( i ) , int( i ) );

   return( - dfct ); 
   }

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

  void ChgCosts( cCRow NCost , cIndex_Set nms = 0 ,
		 Index strt = 0 , Index stp = Inf< Index >() ) override;

  void ChgCost( Index arc , CNumber NCost ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  void ChgQCoef( cCRow NQCoef, cIndex_Set nms = 0 ,
		Index strt = 0 , Index stp = Inf< Index >() ) override;

  void ChgQCoef( Index arc , CNumber NQCoef ) override;

/*--------------------------------------------------------------------------*/
  
  void ChgDfcts( cFRow NDfct , cIndex_Set nms = 0 ,
		 Index strt = 0 , Index stp = Inf< Index >() ) override;

  void ChgDfct( Index node , cFNumber NDfct ) override;

/*--------------------------------------------------------------------------*/

  void ChgUCaps( cFRow NCap , cIndex_Set nms = 0 ,
		 Index strt = 0 , Index stp = Inf< Index >() ) override;

  void ChgUCap( Index arc , FNumber NCap  ) override;

/*--------------------------------------------------------------------------*/
/*----------------- MODIFYING THE STRUCTURE OF THE GRAPH -------------------*/
/*--------------------------------------------------------------------------*/ 

  void CloseArc( Index name ) override;
     
  bool IsClosedArc( Index name ) const override {
   #if( DYNMC_MCF_CPX )
    return( ( ArcPos[ name ] >= 0 ) &&
	    ( ArcPos[ name ] < Inf< FNumber >() ) );
   #else
    return( MCFCplex::MCFUCap( name ) == 0 );
   #endif
   }

  void DelNode( Index name ) override;

  void OpenArc( Index name ) override;

  Index AddNode( FNumber aDfct ) override;

  void ChangeArc( Index name ,
		  Index nSN = Inf< Index >() , Index nEN = Inf< Index >() )
   override;

  void DelArc( Index name ) override;

  bool IsDeletedArc( Index name ) const override {
   #if( DYNMC_MCF_CPX )
    return( ArcPos[ name ] == Inf< FNumber >() );
   #else
    return( false );
   #endif
   }

  Index AddArc( Index Start , Index End , FNumber aU , CNumber aC ) override; 

/*--------------------------------------------------------------------------*/
/*---------------------------- DESTRUCTOR ----------------------------------*/
/*--------------------------------------------------------------------------*/

  ~MCFCplex();  

/*--------------------------------------------------------------------------*/ 
/*----------------- PRIVATE PART OF THE CLASS ------------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

 void MemAlloc( void );

 void MemDeAlloc( void );

 void TurnToQP( void );

 void QPchgarcnode( int name , int sn , int en );

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

 CPXENVptr env;      // Cplex environment pointer 
 CPXNETptr net;      // network pointer
 CPXLPptr qp;        // QP pointer 

 int* Startn;        // arcs' start nodes
 int* Endn;          // arcs' end nodes
 
 QPMethod QPMthd;    // QP solving method

 #if( DYNMC_MCF_CPX )
  FRow ArcPos;       // ArcPos[ i ] < 0 means that arc i exists, == F_INF
                     // means that the position is available for creating
                     // a new arc, anything in between means that the arc
                     // is closed and that is its original capacity
  Index FreePos;     // first position available for creating a new arc, i.e.
                     // smallest index i s.t. ArcPos[ i ] == F_INF
 #endif

/*-------------------------------------------------------------------------*/

 };  // end( class MCFCplex )

/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/

}  // end( namespace MCFClass_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* MCFCplex.h included */

/*--------------------------------------------------------------------------*/
/*-------------------- End File MCFCplex.h ---------------------------------*/
/*--------------------------------------------------------------------------*/
