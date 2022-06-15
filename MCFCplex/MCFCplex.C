/*--------------------------------------------------------------------------*/
/*------------------------- File MCFCplex.C --------------------------------*/
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
 * Copyright &copy by Antonio Frangioni.
 */
/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MCFCplex.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS -----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( CPX_VERSION >= 800 )
 #define CPX_OPTIMAL CPX_STAT_OPTIMAL
 #define CPX_UNBOUNDED CPX_STAT_UNBOUNDED
 #define CPX_INFEASIBLE CPX_STAT_INFEASIBLE

 #define CPX_INForUNBD CPX_STAT_INForUNBD
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace MCFClass_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------- TYPES & CONSTEXPRS -----------------------------*/
/*--------------------------------------------------------------------------*/

using Index = MCFCplex::Index;
using cIndex = MCFCplex::cIndex;
using Index_Set = MCFCplex::Index_Set;
using cIndex_Set = MCFCplex::cIndex_Set;

static constexpr Index IInf = Inf< Index >();

using FNumber = MCFCplex::FNumber;
static constexpr FNumber FInf = Inf< FNumber >();

using CNumber = MCFCplex::CNumber;
static constexpr CNumber CInf = Inf< CNumber >();

/*--------------------------------------------------------------------------*/
/*----------------------------- FUNCTIONS ----------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  These functions are not implemented as methods of the class, since  --*/
/*--  they don't use directly its data structures.                        --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

#if( ( ! FNUMBER_IS_DOUBLE ) || ( ! CNUMBER_IS_DOUBLE ) )

template<class T>
static Index_Set SparseAssign( T* g1 , const double *g2 , Index_Set nms ,
			       Index n , Index Bs = 0 )
{
 // takes the dense n-vector g2 and assigns its nonzeroes to the sparse vector
 // g1, meanwhile constructiong its set of nonzeroes in nms, with names from
 // Bs onwards; returns a pointer to the first element in nms after the last
 // index vritten: this can be used for computing the number of nonzeroes in
 // the "sparsified" vector and/or for InINF-terminating nms[]

 for( ; n-- ; Bs++ ) {
  const double tg2 = *(g2++);
  if( tg2 ) {
   *(nms++) = Bs;
   *(g1++) = T( tg2 );
   }
  }

 return( nms );
 }

/*--------------------------------------------------------------------------*/

template< class T >
static Index_Set SparseAssign( T* g1 , const double *g2 , Index_Set nms ,
			       Index n , double eps , Index Bs = 0 )
{
 // as SparseAssign(), but elements are considered nonzero only if they are
 // >= eps (the idea is that all elements are >= 0)

 for( ; n-- ; Bs++ ) {
  const double tg2 = *(g2++);
  if( tg2 >= eps ) {
   *(nms++) = Bs;
   *(g1++) = T( tg2 );
   }
  }

 return( nms );
 }

#endif

/*--------------------------------------------------------------------------*/

static Index VectLength( cIndex_Set nms , cIndex stp = IInf )
{
 // count the number of elements in the vector nms that are < stp; stops
 // as soon as the first element >= stp is found

 MCFCplex::cIndex_Set tnms = nms;
 while( *tnms < stp )
  tnms++;

 return( tnms - nms );
 }

/*--------------------------------------------------------------------------*/

static void VectFill( int *nms , int strt ,  Index n )
{
 // fills nms with the indices strt, strt + 1, ..., strt + n - 1

 for( ; n-- ; )
  *(nms++) = strt++;
 }

/*--------------------------------------------------------------------------*/

template< class T >
static void VectAssign( T *const g , T x , Index n )
{
 // g[ i ] = x for each i = 0 .. n - 1

 for( T *tg = g + n ; tg > g ; )
  *(--tg) = x;
 }

/*--------------------------------------------------------------------------*/

template< class T >
static void VectAssign( T *const g1 , T *const g2 , Index n )
{
 // g1[ i ] = g2[ i ] for each i = 0 .. n - 1

 const T *tg2 = g2 + n;
 for( T *tg1 = g1 + n ; tg1 > g1 ; )
  *(--tg1) = *(--tg2);
 }

/*--------------------------------------------------------------------------*/

template< class T1 , class T2 >
static void VectMAssign( T1 *g1 , T2 *g2 , Index n )
{
 // g1 := - g2

 for( ; n-- ; )
  *(g1++) = - *(g2++);
 }

/*--------------------------------------------------------------------------*/

template< class T >
static Index_Set Sparsify( T* g , Index_Set B , Index n , Index Bs = 0 )
{
 // turns g from a "dense" n-vector to a "sparse" one, eliminating all items
 // that are exactly == 0; writes the set of nonzero items in B, with names
 // from Bs onwards, ordered in increasing sense; returns a pointer to the
 // first element in B after the last index vritten: this can be used for
 // computing the number of nonzeroes in the "sparsified" vector and/or for
 // InINF-terminating the set

 for( ; n ; n-- , g++ )
  if( *g )
   *(B++) = Bs++;
  else
   break;

 if( n ) {
  T* tg = g++;
  for( Bs++ ; --n ; g++ , Bs++ )
   if( *g ) {
    *(tg++) = *g;
    *(B++) = Bs;
    }
  }

 return( B );

 }  // end( Sparsify )

/*--------------------------------------------------------------------------*/

template< class T >
static Index_Set SparsifyT( T* g , Index_Set B , Index n , T eps ,
			    Index Bs = 0 )
{
 // as Sparsify(), but elements are considered nonzero only if they are >=
 // eps (the idea is that all elements are >= 0)

 for( ; n ; n-- , g++ )
  if( *g >= eps )
   *(B++) = Bs++;
  else
   break;

 if( n ) {
  T* tg = g++;
  for( Bs++ ; --n ; g++ , Bs++ )
   if( *g >= eps ) {
    *(tg++) = *g;
    *(B++) = Bs;
    }
  }

 return( B );

 }  // end( SparsifyT )

/*--------------------------------------------------------------------------*/

template<class T>
static void VectSum( T *const g , const T x , cIndex n )
{
 // g[ i ] += x for each i = 0 .. n - 1

 for( T *tg = g + n ; tg > g ; )
  *(--tg) += x;
 }

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS -----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( USENAME0 )
 #define MINUSONE
 #define PLUSONE
#else
 #define MINUSONE -1
 #define PLUSONE +1
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS MCFCplex --------------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- PUBLIC METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

MCFCplex::MCFCplex( Index nmx , Index mmx ) : MCFClass( nmx , mmx )
{
 // open the environment- - - - - - - - - - - - - - - - - - - - - - - - - - -

 int ts;
 env = CPXopenCPLEX( &ts );
 if( ( ! env ) || ts )
  throw( MCFException( "Problem opening Cplex environment" ) );

 // allocate memory - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( nmax && mmax )
  MemAlloc();
 else
  nmax = mmax = 0;

 }  // end( MCFCplex )

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void MCFCplex::LoadNet( Index nmx , Index mmx , Index pn , Index pm ,
			cFRow pU , cCRow pC , cFRow pDfct , cIndex_Set pSn ,
			cIndex_Set pEn )
{
 // allocating and deallocating memory- - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( mmax && nmax ) && ( ( nmx != nmax ) || ( mmx != mmax ) ) ) {
  MemDeAlloc();
  nmax = mmax = 0;
  }

 if( ( mmx && nmx ) && ( ( nmx != nmax ) || ( mmx != mmax ) ) ) {
  nmax = nmx;
  mmax = mmx;
  MemAlloc();
  }

 if( ( ! nmax ) || ( ! mmax ) )  // just sit down in the corner and wait
  return;

 // now setting up data - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 n = pn;
 m = pm;

 // setup data structures for arc creation/deletion- - - - - - - - - - - - -

 #if( DYNMC_MCF_CPX )
  VectAssign( ArcPos , FNumber( -1 ) , FreePos = m );
  VectAssign( ArcPos + m , FInf , mmax - m );
 #endif

 // create and set up temporary data structures - - - - - - - - - - - - - - -

 #if( ( ! INDEX_IS_UINT ) || ( ! USENAME0 ) )
  int* stn = new int[ m ];
  int* enn = new int[ m ];

  for( Index i = m ; i-- ; ) {
   stn[ i ] = pSn[ i ] MINUSONE;
   enn[ i ] = pEn[ i ] MINUSONE;
   }
 #else
  #define stn (int *) pSn
  #define enn (int *) pEn
 #endif

 double* sup = new double[ n ];

 if( pDfct )
  VectMAssign( sup , pDfct , n );  // invert the sign of deficits
 else
  VectAssign( sup , double( 0 ) , n );

 double* upc = new double[ m ];

 if( pU )
  for( Index i = m ; i-- ; )
   if( pU[ i ] == FInf )
    upc[ i ] = CPX_INFBOUND;
   else
    upc[ i ] = double( pU[ i ] );
 else
  VectAssign( upc , CPX_INFBOUND , m );

 double* obj = new double[ m ];

 if( pC )
  for( Index i = m ; i-- ; )
   if( pC[ i ] == CInf ) {
    #if( DYNMC_MCF_CPX )
     ArcPos[ i ] = pU[ i ];
    #endif
    obj[ i ] = upc[ i ] = 0;
    }
   else
    obj[ i ] = double( pC[ i ] );
 else
  VectAssign( obj , double( 0 ) , m );

 #if( DYNMC_MCF_CPX )
  while( FreePos && ( ArcPos[ FreePos - 1 ] == FInf ) )
   FreePos--;
 #endif

 // load internal structure of Cplex- - - - - - - - - - - - - - - - - - - - -

 status = CPXNETcopynet( env , net , CPX_MIN , n , sup , NULL , m , stn ,
			 enn , NULL , upc , obj , NULL );
 if( status )
  throw( MCFException( "Problem loading data" ) );

 // setup QP data - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 qp = NULL;             // problem is initially Network
 QPMthd = qpNSimplex;   // set default QP solving method to Network Simplex
 

 // delete temporaries- - - - - - - - - - - - - - - - - - - - - - - - - - - -

 delete[] obj;
 delete[] upc;
 delete[] sup;

 #if( ( ! INDEX_IS_UINT ) || ( ! USENAME0 ) )
  delete[] enn;
  delete[] stn;
 #endif

 status = MCFClass::kUnSolved;

 }  // end( MCFCplex::LoadNet )

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

void MCFCplex::SolveMCF( void )
{
 if( MCFt )
  MCFt->Start();
 
 if( net ) {
  CPXNETprimopt( env , net );  // call the network simplex- - - - - - - - - - 
  status = CPXNETgetstat( env , net );
  }
 else {
  CPXqpopt( env , qp );       // call the QP solver - - - - - - - - - - - - -
  status = CPXgetstat( env , qp );
  }

 switch( status ) {
   case( CPX_OPTIMAL ):       status = MCFClass::kOK;
                              break;
   case( CPX_INForUNBD ):
   case( CPX_INFEASIBLE ):    status = MCFClass::kUnfeasible;
                              break;
   case( CPX_UNBOUNDED ):     status = MCFClass::kUnbounded;
                              break;
   default:                   status = MCFClass::kStopped;
   }

 if( MCFt )
  MCFt->Stop();

 }  // end( MCFCplex::SolveMCF )

/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR READING RESULTS  -------------------------*/
/*--------------------------------------------------------------------------*/

void MCFCplex::MCFGetX( FRow F , Index_Set nms , Index strt , Index stp )
 const {
 if( stp > m )
  stp = m;

 #if( FNUMBER_IS_DOUBLE )
  if( net ) 
   CPXNETgetx( env , net , F , strt , stp - 1 );
  else
   CPXgetx( env , qp , F , strt , stp - 1 );
   
  if( nms )
   #if( EPS_FLOW )
    *SparsifyT( F , nms , stp - strt , EpsFlw , strt ) = IInf;
   #else
    *Sparsify( F , nms , stp - strt , strt ) = IInf;
   #endif
 #else
  auto val = new double[ stp - strt ];
  if( net ) 
   CPXNETgetx( env , net , val , strt , stp - 1 );
  else
   CPXgetx( env , qp , val , strt , stp - 1 );

  if( nms )
   #if( EPS_FLOW )
    *SparseAssign( F , val , nms , stp - strt , EpsFlw , strt ) = IInf;
   #else
    *SparseAssign( F , val , nms , stp - strt , strt ) = IInf;
   #endif
  else
   VectAssign( F , val , stp - strt );
  delete[] val;
 #endif

 }  // end( MCFCplex::MCFGetX( F ) )

/*--------------------------------------------------------------------------*/

void MCFCplex::MCFGetRC( CRow CR , cIndex_Set nms ,
			 Index strt , Index stp ) const
{
 if( stp > m )
  stp = m;

 if( strt >= stp )
  return;

 if( nms ) {
  auto val = new double[ stp - strt ];

  if( net )
   CPXNETgetdj( env , net , val , strt , stp - 1 );
  else
   CPXgetdj( env , qp , val , strt , stp - 1 );

  while( *nms < strt )
   nms++;

  for( Index h ; ( h = *(nms++) ) < stp ; )
   *(CR++) = val[ h - strt ];

  delete[] val;
  }
 else {
  #if( CNUMBER_IS_DOUBLE )
   if( net )
    CPXNETgetdj( env , net , CR , strt , stp - 1 );
   else
    CPXgetdj( env , qp , CR , strt , stp - 1);
  #else
   auto val = new double[ stp - strt ];
   if( net )
    CPXNETgetdj( env , net , val , strt , stp - 1 );
   else
    CPXgetdj( env , qp , val , strt , stp - 1);

   VectAssign( CR , val , stp - strt );
   delete[] val;
  #endif
  }
 }  // end( MCFCplex::MCFGetRC( CR ) )

/*--------------------------------------------------------------------------*/

MCFCplex::CNumber MCFCplex::MCFGetRC( Index i ) const
{
 double temp;
 if( net )
  CPXNETgetdj( env , net , &temp , int( i ) , int( i ) );
 else
  CPXgetdj( env , qp , &temp , int( i ) , int( i ) );

 return( CNumber( temp ) );

 }  // end( MCFCplex::MCFGetRC( i ) )

/*--------------------------------------------------------------------------*/

void MCFCplex::MCFGetPi( CRow P , cIndex_Set nms ,
			 Index strt , Index stp ) const
{
 if( stp > n )
  stp = n;

 if( nms ) {
  auto val = new double[ stp - strt ];

  if( net )
   CPXNETgetpi( env , net , val , strt , stp - 1 ); 
  else
   CPXgetpi( env , qp , val , strt , stp - 1 ); 

  while( *nms < strt )
   nms++;

  for( Index h ; ( h = *(nms++) ) < stp ; )
   *(P++) = - val[ h - strt ];

  delete[] val;
  }
 else {
  #if( CNUMBER_IS_DOUBLE ) 
   if( net )
    CPXNETgetpi( env , net , P , strt , stp - 1 );          
   else
    CPXgetpi( env , qp , P , strt , stp - 1 );          

   for( Index i = 0 ; i < stp - strt ; i++ )
    P[ i ] = - P[ i ];
  #else
   auto val = new double[ stp - strt ];
   if( net )	   
    CPXNETgetpi( env , net , val , strt , stp - 1 );
   else
    CPXgetpi( env , qp , val , strt , stp - 1 );

   VectMAssign( P , val , stp - strt );
   delete[] val;
  #endif
  }    
 }  // end( MCFGetPi( P ) )

/*--------------------------------------------------------------------------*/

MCFCplex::FONumber MCFCplex::MCFGetFO( void ) const
{
 if( status == MCFClass::kOK ) {
  double objval;
  if( net )
   CPXNETgetobjval( env , net , &objval );
  else
   CPXgetobjval(env, qp, &objval);

  return( FONumber( objval ) );  
  }
 else
  if( status == MCFClass::kUnbounded )
   return( - Inf< FONumber >() );
  else
   return( Inf< FONumber >() );

 }  // end( MCFCplex::MCFGetFO )

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

void MCFCplex::MCFArcs( Index_Set Startv , Index_Set Endv ,
			cIndex_Set nms , Index strt , Index stp ) const
{
 if( stp > m )
  stp = m;

 if( strt >= stp )
  return;
 
 if( nms ) {
  while( *nms < strt )
   nms++;

  for( Index i ; ( i = *(nms++) ) < stp ; ) {
   int st , en;
   if( net ) 
    CPXNETgetarcnodes( env , net , &st , &en , int( i ) , int( i ) );
   else {
    st = Startn[ i ];
    en = Endn[ i ];
    }

   if( Startv )
    *(Startv++) = Index( st PLUSONE );
   if( Endv )
    *(Endv++) = Index( en PLUSONE );
   }
  }
 else {
  if( net ) {
   #if( INDEX_IS_UINT )
    CPXNETgetarcnodes( env , net , (int *) Startv , (int *) Endv , strt ,
		       stp - 1 );
   #else
    auto ind = new int[ stp - strt ];
    if( Startv ) {
     CPXNETgetarcnodes( env , net , ind , NULL , strt , stp - 1 );
     VectAssign( Startv , ind , stp - strt );
     }

    if( Endv ) {
     CPXNETgetarcnodes( env , net , NULL , ind , strt , stp - 1 );
     VectAssign( Endv , ind , stp - strt );
     }

    delete[] ind;
   #endif
   }
  else {
   if( Startv )
    VectAssign( Startv , Index_Set( Startn + strt ) , stp - strt );
  
   if( Endv )
    VectAssign( Endv , Index_Set( Endn + strt ) , stp - strt );
   }

  #if( ! USENAME0 )
   if( Startv )
    VectSum( Startv , Index( 1 ) , stp - strt );

   if( Endv )
    VectSum( Endv , Index( 1 ) , stp - strt );
  #endif
  }
 }  // end( MCFArcs )

/*-------------------------------------------------------------------------*/

void MCFCplex::MCFCosts( CRow Costv , cIndex_Set nms ,
			 Index strt , Index stp ) const
{
 if( stp > m )
  stp = m;

 if( strt >= stp )
  return;

 if( nms ) {
  while( *nms < strt )
   nms++;

  for( Index i ; ( i = *(nms++) ) < stp ; ) {
   double cst;
   if( net )
    CPXNETgetobj( env , net , &cst , int( i ) , int( i ) );
   else
    CPXgetobj( env , qp , &cst , int( i ) , int( i ) );

   *(Costv++) = cst;
   }
  }
 else {
  #if( CNUMBER_IS_DOUBLE )
   if( net )
    CPXNETgetobj( env , net , Costv , strt , stp - 1 );
   else
    CPXgetobj( env , qp , Costv , strt , stp - 1 );
  #else
   auto val = new double[ stp - strt ];
   if( net )
    CPXNETgetobj( env , net , val , strt , stp - 1 );
   else
    CPXgetobj( env , qp , val , strt , stp - 1 );

   VectAssign( Costv , val , stp - strt );
   delete[] val;
  #endif
  }
 }  // end( MCFCplex::MCFCosts )

/*-------------------------------------------------------------------------*/

void MCFCplex::MCFQCoef( CRow Qv , cIndex_Set nms ,
			 Index strt , Index stp ) const
{
 if( stp > m )
  stp = m;

 if( strt >= stp )
  return;

 double qcoef = 0;
 if( nms ) {
  while( *nms < strt )
   nms++;

  for( Index i ; ( i = *(nms++) ) < stp ; ) {
   if( qp )
    CPXgetqpcoef( env , qp , int( i ) , int( i ) , &qcoef );

   *(Qv++) = qcoef;
   }    
  }
 else {
  for( Index i = strt ; i < stp ; i++ ) {
   if( qp )
    CPXgetqpcoef( env , qp , int( i ) , int( i ) , &qcoef );

   *(Qv++) = qcoef;
   }
  }
 } // end( MCFCplex::MCFQCoef )

/*-------------------------------------------------------------------------*/

void MCFCplex::MCFUCaps( FRow UCapv , cIndex_Set nms ,
			 Index strt , Index stp ) const
{
 if( stp > m )
  stp = m;

 if( strt >= stp )
  return;

 if( nms ) {
  while( *nms < strt )
   nms++;

  for( Index i ; ( i = *(nms++) ) < stp ; )
   #if( DYNMC_MCF_CPX )
    if( ( ArcPos[ i ] >= 0 ) && ( ArcPos[ i ] < FInf ) )
     *(UCapv++) = ArcPos[ i ];
    else
   #endif
    {
     double cap;
     if( net )
      CPXNETgetub( env , net , &cap , int( i ) , int( i ) );
     else
      CPXgetub( env , qp , &cap , int( i ) , int( i ) );

     *(UCapv++) = cap;
     }
  }
 else {
  #if( FNUMBER_IS_DOUBLE )
   if( net )
    CPXNETgetub( env , net , UCapv , strt , stp - 1 );  
   else
    CPXgetub( env , qp , UCapv , strt , stp - 1 );
  #else
   auto val = new double[ stp - strt ];
   if( net )
    CPXNETgetub( env , net , val , strt , stp - 1 );
   else
    CPXgetub( env , qp , val , strt , stp - 1 );

   VectAssign( UCapv , val , stp - strt );
   delete[] val;
  #endif

  #if( DYNMC_MCF_CPX )
   for( Index i = strt ; i < stp ; i++ )
    if( ( ArcPos[ i ] >= 0 ) && ( ArcPos[ i ] < FInf ) )
     UCapv[ i ] = ArcPos[ i ];
  #endif
  }
 }  // end( MCFCplex::MCFUCaps )

/*-------------------------------------------------------------------------*/

void MCFCplex::MCFDfcts( FRow Dfctv , cIndex_Set nms ,
			 Index strt , Index stp ) const
{
 if( stp > n )
  stp = n;

 if( strt >= stp )
  return;

 if( nms ) {
  while( *nms < strt )
   nms++;

  for( Index i ; ( i = *(nms++) ) < stp ; ) {
   double dfct;
   if( net )
    CPXNETgetsupply( env , net , &dfct , int( i ) , int( i ) );
   else
    CPXgetrhs( env , qp , &dfct , int( i ) , int( i ) );
   
   *(Dfctv++) = - dfct;
   }
  }
 else {
  #if( FNUMBER_IS_DOUBLE )
   if( net )
    CPXNETgetsupply( env , net , Dfctv , strt , stp - 1 ); 
   else
    CPXgetrhs( env , qp  , Dfctv , strt , stp - 1 ); 

   for( FRow tDf = Dfctv + stp - strt ; tDf-- > Dfctv ; )
    *tDf = - *tDf;
  #else
   auto val = new double[ stp - strt ];
   if( net )
    CPXNETgetsupply( env , net , val , strt , stp - 1 );
   else
    CPXgetrhs( env , qp , val , strt , stp - 1 );

   VectMAssign( Dfctv , val , stp - strt );
   delete[] val;
  #endif
  }
 }  // end( MCFCplex::MCFDfcts )

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

void MCFCplex::ChgCosts( cCRow NCost , cIndex_Set nms ,
			 Index strt , Index stp )
{
 if( stp > m )
  stp = m;

 if( strt >= stp )
  return;

 Index cnt;
 if( nms ) {
  while( *nms < strt ) {
   nms++;
   NCost++;
   }

  cnt = VectLength( nms , stp );
  if( ! cnt )
   return;
  }
 else
  cnt = stp - strt;

 #undef MCFCplex_VALUEPASSED
 #undef MCFCplex_INDEXPASSED

 #if( CNUMBER_IS_DOUBLE )
  #define MCFCplex_VALUEPASSED CRow( NCost )
 #else
  auto val = new double[ cnt ];
  VectAssign( val , NCost , cnt );
  #define MCFCplex_VALUEPASSED val
 #endif

 if( nms ) {
  #if( INDEX_IS_UINT )
   #define MCFCplex_INDEXPASSED (int *) nms
  #else
   auto ind = new int[ cnt ];
   VectAssign( ind , nms , cnt );
   #define MCFCplex_INDEXPASSED ind
  #endif

  if( net )
   CPXNETchgobj( env , net , cnt , MCFCplex_INDEXPASSED ,
		 MCFCplex_VALUEPASSED );
  else
   CPXchgobj( env , qp , cnt , MCFCplex_INDEXPASSED , MCFCplex_VALUEPASSED );

  #if( ! INDEX_IS_UINT )
   delete[] ind;
  #endif
  }
 else {
  auto ind = new int[ cnt ];
  VectFill( ind , strt , cnt );
  
  if( net )
   CPXNETchgobj( env , net , cnt , ind , MCFCplex_VALUEPASSED );
  else
   CPXchgobj( env , qp , cnt , ind , MCFCplex_VALUEPASSED );

  delete[] ind;
  }

 #if( ! CNUMBER_IS_DOUBLE )
  delete[] val;
 #endif

 }  // end( MCFCplex::ChgCosts )

/*--------------------------------------------------------------------------*/

void MCFCplex::ChgCost( Index arc , CNumber NCost )
{
 if( arc >= MCFm() )
  throw( MCFException( "MCFCplex::ChgCost: invalid arc name" ) );

 double cost = NCost;
 int which = arc;

 if( net )
  CPXNETchgobj( env , net , 1 , &which , &cost );
 else
  CPXchgobj( env , qp  , 1 , &which , &cost);

 }  // end( MCFCplex::ChgCost )

/*--------------------------------------------------------------------------*/

void MCFCplex::ChgQCoef( cCRow NQCoef , cIndex_Set nms ,
			 Index strt , Index stp )
{
 if( stp > m )
  stp = m;

 if( strt >= stp )
  return;

 if( nms ) {
  while( *nms < strt ) {
   nms++;
   if( NQCoef )
    NQCoef++;
   }

  if( ! VectLength( nms , stp ) )
   return;
  }

 if( ( ! qp ) && NQCoef )  // the problem in not QP
  TurnToQP();              // turn it into QP

 double qcoef = 0;
 if( nms ) {
  for( Index arc ; ( arc = *(nms++) ) < stp ; ) {
   if( NQCoef )
    qcoef = *(NQCoef++);

   CPXchgqpcoef( env , qp , int( arc ) , int( arc ) , qcoef );
   }
  }
 else {
  for( Index arc = strt ; arc < stp ; ++arc ) {
   if( NQCoef )
    qcoef = *(NQCoef++);

   CPXchgqpcoef( env , qp , int( arc ) , int( arc ) , qcoef );
   }
  }
 }  // end( MCFCplex::ChgQCoef )

/*--------------------------------------------------------------------------*/

void MCFCplex::ChgQCoef( Index arc , CNumber NQCoef ) 
{
 if( arc >= MCFm() )
  throw( MCFException( "MCFCplex::ChgQCoef: invalid arc name" ) );

 if( ( ! qp ) && NQCoef )
  TurnToQP();

 CPXchgqpcoef( env , qp , int( arc ) , int( arc ) , NQCoef );

 }  // end( MCFCplex::ChgQCoef )

/*--------------------------------------------------------------------------*/
   
void MCFCplex::ChgDfcts( cFRow NDfct , cIndex_Set nms ,
			 Index strt , Index stp ) 
{
 if( stp > n )
  stp = n;

 if( strt >= stp )
  return;

 Index cnt;
 if( nms ) {
  while( *nms < strt ) {
   nms++;
   NDfct++;
   }

  cnt = VectLength( nms , stp );
  if( ! cnt )
   return;
  }
 else
  cnt = stp - strt;

 #undef MCFCplex_INDEXPASSED

 auto val = new double[ cnt ];
 VectMAssign( val , NDfct , cnt );

 if( nms ) {
  #if( INDEX_IS_UINT )
   #define MCFCplex_INDEXPASSED (int *) nms
  #else
   auto ind = new int[ cnt ];
   VectAssign( ind , nms , cnt );
   #define MCFCplex_INDEXPASSED ind
  #endif

  if( net )
   CPXNETchgsupply( env , net , cnt , MCFCplex_INDEXPASSED , val );
  else
   CPXchgrhs( env , qp , cnt , MCFCplex_INDEXPASSED , val );

  #if( ! INDEX_IS_UINT )
   delete[] ind;
  #endif
  }
 else {
  auto ind = new int[ cnt ];
  VectFill( ind , strt , cnt );

  if( net )
   CPXNETchgsupply( env , net , cnt , ind , val );
  else
   CPXchgrhs( env , qp , cnt , ind , val );

  delete[] ind;
  }

 delete[] val;
 
 }  // end( MCFCplex::ChgDfcts )

/*--------------------------------------------------------------------------*/

void MCFCplex::ChgDfct( Index node , FNumber NDfct )
{
 if( node >= MCFn() )
  throw( MCFException( "MCFCplex::ChgDfct: invalid node name" ) );

 double dfct = - NDfct;
 int which = node;

 if( net ) 
  CPXNETchgsupply( env , net , 1 , &which , &dfct );
 else
  CPXchgrhs( env , qp , 1 , &which , &dfct );

 }  // end( MCFCplex::ChgDfct )

/*--------------------------------------------------------------------------*/
 
void MCFCplex::ChgUCaps( cFRow NCap , cIndex_Set nms ,
			 Index strt , Index stp )
{
 if( stp > m )
  stp = m;

 if( strt >= stp )
  return;

 Index cnt;
 if( nms ) {
  while( *nms < strt ) {
   nms++;
   NCap++;
   }

  cnt = VectLength( nms , stp );
  if( ! cnt )
   return;
  }
 else
  cnt = stp - strt;

 auto ChangeUB = new char[ cnt ];
 VectAssign( ChangeUB , 'U' , cnt );

 #undef MCFCplex_VALUEPASSED
 #undef MCFCplex_INDEXPASSED

 #if( FNUMBER_IS_DOUBLE && ( ! DYNMC_MCF_CPX ) )
  #define MCFCplex_VALUEPASSED FRow( NCap )
 #else
  auto val = new double[ cnt ];
  #if( DYNMC_MCF_CPX )
   // for all arcs that are closed, store the new value in ArcPos[] and
   // leave the "new" capacity just being set to 0
   auto tval = val;
   if( nms ) {
    auto tnms = nms;
    for( Index i ; ( i = *(tnms++) ) < stp ; ) {
     if( ( ArcPos[ i ] >= 0 ) && ( ArcPos[ i ] < FInf ) ) {
      ArcPos[ i ] = *(NCap++);
      *(tval++) = 0;
      }
     else
      *(tval++) = *(NCap++);
     }
    }
   else
    for( Index i = strt ; i < stp ; ++i ) {
     if( ( ArcPos[ i ] >= 0 ) && ( ArcPos[ i ] < FInf ) ) {
      ArcPos[ i ] = *(NCap++);
      *(tval++) = 0;
      }
     else
      *(tval++) = *(NCap++);
     }
  #else
   VectAssign( val , NCap , cnt );
  #endif
  #define MCFCplex_VALUEPASSED val
 #endif

 if( nms ) {
  #if( INDEX_IS_UINT )
   #define MCFCplex_INDEXPASSED (int *) nms
  #else 
   auto ind = new int[ cnt ];
   VectAssign( ind , nms , cnt );
   #define MCFCplex_INDEXPASSED ind
  #endif

  if( net )
   CPXNETchgbds( env , net , cnt , MCFCplex_INDEXPASSED , ChangeUB ,
		 MCFCplex_VALUEPASSED );
  else
   CPXchgbds( env , qp , cnt , MCFCplex_INDEXPASSED , ChangeUB ,
	      MCFCplex_VALUEPASSED );

  #if( ! INDEX_IS_UINT )
   delete[] ind;
  #endif
  }
 else {
  auto ind = new int[ cnt ];
  VectFill( ind , strt , cnt );

  if( net )
   CPXNETchgbds( env , net , cnt , ind , ChangeUB , MCFCplex_VALUEPASSED );
  else
   CPXchgbds( env , qp , cnt , ind , ChangeUB , MCFCplex_VALUEPASSED );

  delete[] ind;
  }

 #if( ( ! FNUMBER_IS_DOUBLE ) || DYNMC_MCF_CPX )
  delete[] val;
 #endif

 delete[] ChangeUB;
 
 }  // end( MCFCplex::ChgUCaps )

/*--------------------------------------------------------------------------*/

void MCFCplex::ChgUCap( Index arc , FNumber NCap )
{
 if( arc >= MCFm() )
  throw( MCFException( "MCFCplex::ChgUCap: invalid arc name" ) );

 #if( DYNMC_MCF_CPX )
  if( ArcPos[ arc ] == FInf )  // the arc is deleted
   return;                   // silently return

  // if the arc is closed, change the capacity stored in ArcPos but leave
  // it in its closed state
  if( ArcPos[ arc ] >= 0 ) {
   ArcPos[ arc ] = NCap;
   return;
   }
 #endif

 double cap = NCap;
 int which = arc;

 if( net )
  CPXNETchgbds( env , net , 1 , &which , "U", &cap );
 else
  CPXchgbds( env , qp , 1 , &which , "U" , &cap );

 }  // end( MCFCplex::ChgUCap )

/*--------------------------------------------------------------------------*/
/*----------------- MODIFYING THE STRUCTURE OF THE GRAPH -------------------*/
/*--------------------------------------------------------------------------*/ 

void MCFCplex::CloseArc( Index name )
{
 if( name >= MCFm() )
  throw( MCFException( "MCFCplex::CloseArc: invalid arc name" ) );

 #if( DYNMC_MCF_CPX )
  // if the arc is closed already, or there is no arc in that position
  if( ArcPos[ name ] >= 0 )
   return;  // nothing to do
 #endif

 double ub;
 int temp = name;

 #if( DYNMC_MCF_CPX )
  if( net )
   CPXNETgetub( env , net , &ub , temp , temp );
  else
   CPXgetub( env , qp , &ub , temp , temp );

  ArcPos[ name ] = ub;  // save current upper bound
 #endif

 // set upper bound to 0
 ub = 0;
 if( net )
  CPXNETchgbds( env , net , 1 , &temp , "U" , &ub );
 else
  CPXchgbds( env , qp , 1 , &temp , "U" , &ub );

 }  // end( MCFCplex::CloseArc )

/*--------------------------------------------------------------------------*/

void MCFCplex::DelNode( Index name )
{
 int index = int( name MINUSONE );
 double newdfct = 0;

 if( net ) { 
  int arcbeg = 0 , surplus = 0 , arccount_p = 0;
  CPXNETchgsupply( env , net , 1 , &index , &newdfct );  // set supply to 0

  auto ind = new int[ m ];
  CPXNETgetnodearcs( env , net , &arccount_p , &arcbeg , ind , m , &surplus ,
		     index , index );  

  for( int i = 0 ; i < arccount_p ; )
   CloseArc( Index( ind[ i++ ] ) );  // incident arcs are "closed"

  delete[] ind;
  }
 else {
  for( Index arc = 0 ; arc < m ; arc++ )
   if( ( Startn[ arc ] == index ) || ( Endn[ arc ] == index ) )
    CloseArc( arc );  // incident arcs are "closed"
    
  CPXchgrhs( env , qp , 1 , &index , &newdfct );  // set rhs to 0
  }
 }  // end( MCFCplex::DelNode )

/*--------------------------------------------------------------------------*/

void MCFCplex::OpenArc( Index name )
{
 if( name >= MCFm() )
  throw( MCFException( "MCFCplex::OpenArc: invalid arc name" ) );

 #if( DYNMC_MCF_CPX )
  // if the arc exists and it is opened already
  if( ArcPos[ name ] < 0 )
   return;  // nothing to do

  if( ArcPos[ name ] == FInf )  // the arc is deleted
   throw( MCFException( "MCFCplex::OpenArc: cannot open a deleted arc" ) );

  int temp = name;
  double ub = ArcPos[ name ];
  ArcPos[ name ] = -1;

  // restore upper bound
  if( net )
   CPXNETchgbds( env , net , 1 , &temp , "U" , &ub ); 
  else
   CPXchgbds( env , qp , 1 , &temp , "U" , &ub ); 

 #else
  throw( MCFException(
	     "MCFCplex::OpenArc() not implemented if DYNMC_MCF_CPX == 0" ) );
 #endif

 }  // end( MCFCplex::OpenArc )
 
/*--------------------------------------------------------------------------*/

MCFCplex::Index MCFCplex::AddNode( FNumber aDfct )
{
 if( n < nmax ) {
  double dfct = aDfct;
 
  if( net ) 
   CPXNETaddnodes( env , net , 1 , &dfct , NULL );
  else
   CPXnewrows( env, qp , 1 , &dfct , "E" , NULL , NULL ); // add new row in
                                                          // constraint matrix 
  n++;
  }

 return( n );

 }  // end( MCFCplex::AddNode )

/*--------------------------------------------------------------------------*/

void MCFCplex::DelArc( Index name ) 
{
 if( name >= MCFm() )
  throw( MCFException( "MCFCplex::DelArc: invalid arc name" ) );

 #if( DYNMC_MCF_CPX )
  if( ArcPos[ name ] == FInf )  // the arc is deleted already
   return;                      // nothing to do

  bool clsd = ( ArcPos[ name ] >= 0 );  // was it closed already?
  ArcPos[ name ] = FInf;        // position now is available for a new arc
  if( name < FreePos )
   FreePos = name;

  int which = name;

  if( name == m - 1 ) {  // deleting the last arc(s) - - - - - - - - - - - -
   do {
    m--;                 // decrement number of arcs
 
    if( net )
     CPXNETdelarcs( env , net , which , which );
    else
     CPXdelcols( env , qp , which , which );  // delete matching column in
                                              // constraint matrix

    } while( ArcPos[ m - 1 ] == FInf );

   if( FreePos > m ) 
    FreePos = m; 
   }
  else                 // just set the bound to 0- - - - - - - - - - - - - -
   if( ! clsd ) {      // unless the arc wass closed, the bound is 0 already
    double cpct = 0; 
    if( net )
     CPXNETchgbds( env , net , 1 , &which , "U" , &cpct );
    else
     CPXchgbds( env , qp , 1 , &which , "U" , &cpct );
    }
 #else
  throw( MCFException(
	     "MCFCplex::DelArc() not implemented if DYNMC_MCF_CPX == 0" ) );
 #endif

 }  // end( MCFCplex::DelArc )

/*-------------------------------------------------------------------------*/

MCFCplex::Index MCFCplex::AddArc( Index Start , Index End , FNumber aU ,
				  CNumber aC )
{
 #if( DYNMC_MCF_CPX )
  int strn = int( Start MINUSONE );
  int endn = int( End MINUSONE );
 
  double cst = aC;
  double cpct = aU;

  if( FreePos >= m ) {  // the first free position is at the end
   if( m >= mmax )      // but there is no space left
    return( IInf );     // operation failed
   FreePos = ++m;       // increase number of arcs

   int NewArc = m - 1;  // now physically create the new arc
   ArcPos[ NewArc ] = -1;  // mark the arc as existent

   if( net ) 
    CPXNETaddarcs( env , net , 1 , &strn , &endn , NULL , &cpct , &cst ,
		   NULL );
   else {
    // add new column in objective function and constraint matrix
    CPXnewcols( env , qp , 1 , &cst , NULL, &cpct , NULL , NULL );
    
    // modify arc-node incidence matrix
    CPXchgcoef( env , qp , strn , NewArc , 1 );   
    CPXchgcoef( env , qp , endn , NewArc , -1 );

    Startn[ NewArc ] = strn;
    Endn[ NewArc ]   = endn;
    }

   return( Index( NewArc ) );
   }
  else {                // the first free position is in the middle
   int temp = FreePos;  // just fill the hole

   if( net ) {
    CPXNETchgarcnodes( env , net , 1 , &temp , &strn , &endn );
    CPXNETchgobj( env , net , 1 , &temp , &cst );
    CPXNETchgbds( env , net , 1 , &temp , "U" , &cpct );
    }
   else {
    QPchgarcnode( temp, strn, endn );
    CPXchgobj( env , qp , 1 , &temp , &cst );
    CPXchgbds( env , qp , 1 , &temp , "U" , &cpct );

    Startn[ temp ] = strn;
    Endn[ temp ] = endn;
    }

   ArcPos[ FreePos ] = -1;  // mark the arc as existent
   // look for the first non-existent one
   while( ( FreePos < mmax ) && ( ArcPos[ FreePos ] < FInf ) )
    FreePos++;

   return( temp );
   }
 #else
  throw( MCFException(
	      "MCFCplex::AddArc() not implemented if DYNMC_MCF_CPX == 0" ) );

  return( IInf );
 #endif

 }  // end( MCFCplex::AddArc )

/*-------------------------------------------------------------------------*/

void MCFCplex::ChangeArc( Index name , Index nSN , Index nEN ) 
{
 if( name >= MCFm() )
  throw( MCFException( "MCFCplex::ChangeArc: invalid arc name" ) );

 int temp = int( name );
 int sn , en;
 if( ( nSN == IInf ) || ( nEN == IInf ) ) {
  if( net ) {
   CPXNETgetarcnodes( env, net, &sn, &en, temp, temp );
   }
  else {
   sn = Startn[ name ];
   en = Endn[ name ];
   }
  }

 if( nSN < IInf )
  sn = int( nSN MINUSONE );

 if( nEN < IInf )
  en = int( nEN MINUSONE );

 if( net ) 
  CPXNETchgarcnodes( env , net , 1 , &temp , &sn , &en ) ;
 else
  QPchgarcnode( temp , sn , en );
       
 }  // end( MCFCplex::ChangeArc )

/*--------------------------------------------------------------------------*/
/*---------------------------- DESTRUCTOR ----------------------------------*/
/*-------------------------------------------------------------------------*/

MCFCplex::~MCFCplex()
{
 if( mmax && nmax )
  MemDeAlloc();

 CPXcloseCPLEX( &env );

 }  // end( MCFCplex::~MCFCplex )

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

void MCFCplex::MemAlloc( void )
{
 // create a new Cplex problem- - - - - - - - - - - - - - - - - - - - - - - -

 int ts;
 net = CPXNETcreateprob( env , &ts , "NET" );
 if( ( ! net ) || ts )
  throw( MCFException( "Problem creating Cplex problem" ) );

 // create data structures for arc creation/deletion- - - - - - - - - - - - -

 #if( DYNMC_MCF_CPX )
  ArcPos = new FNumber[ mmax ];
 #endif

 }  // end( MemAlloc )

/*--------------------------------------------------------------------------*/

void MCFCplex::MemDeAlloc( void )
{
 #if( DYNMC_MCF_CPX )
  delete[] ArcPos;
 #endif

 if( net ) 
  CPXNETfreeprob( env , &net ); 
 else {
  CPXfreeprob( env , &qp );
  delete[] Startn;
  delete[] Endn;
  }
 }  // end( MemDeAlloc )

/*--------------------------------------------------------------------------*/

void MCFCplex::TurnToQP( void )
{
 int status_p;
 qp = CPXcreateprob( env , &status_p , "QP" );
 if( ( ! qp ) || status_p )
  throw( MCFException( "Problem creating Cplex LP" ) );

 Startn = new int[ mmax ];
 Endn   = new int[ mmax ];

 CPXNETgetarcnodes( env , net , Startn , Endn , 0 , m - 1 );

 status_p = CPXcopynettolp( env , qp , net );
 if( ( ! qp ) || status_p )
  throw( MCFException( "Problem creating copying NET to LP" ) );

 CPXNETfreeprob( env , &net );
 net = NULL;

 status_p = CPXchgprobtype( env , qp , CPXPROB_QP );
 if( status_p )
  throw( MCFException( "Problem changing problem type to QP" ) );
 
 status_p = CPXsetintparam( env , CPX_PARAM_QPMETHOD , QPMthd );
 if( status_p )
  throw( MCFException( "Problem setting QP solving method" ) );
 }

/*--------------------------------------------------------------------------*/

void MCFCplex::QPchgarcnode( int name , int sn , int en ) 
{
 // Modify node-arc incidence matrix - - - - - - - - - - - - - - - - - - - - -

 CPXchgcoef( env , qp , Startn[ name ] , name , 0 );
 CPXchgcoef( env , qp , Endn[ name ] , name , 0 );
 // set to 0 components of the node-arc incidence matrix corresponding to
 // arc "name"

 CPXchgcoef( env , qp , sn , name , 1 );  // setup new
 CPXchgcoef( env , qp , en , name , -1 ); // components

 Startn[ name ] = sn; 
 Endn[ name ]   = en;
 }

/*--------------------------------------------------------------------------*/
/*-------------------- End File MCFCplex.C ---------------------------------*/
/*--------------------------------------------------------------------------*/
