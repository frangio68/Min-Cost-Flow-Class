/*--------------------------------------------------------------------------*/
/*-------------------------- File MCFClone.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Testing template class for Min Cost Flow Problem solvers deriving from
 * MCFClass. It is a template deriving from Master and holding an object of
 * class Slave (both deriving from MCFClass), that does whatever Master does
 * but also does whatever Slave does. The output is that of Master.
 *
 * \version 2.00
 *
 * \date 07 - 05 - 2012
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright(C) 1992 - 2012 Antonio Frangioni
 */

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __MCFClone
 #define __MCFClone  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MCFClass.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE and USING ----------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
namespace MCFClass_di_unipi_it
{
 using namespace OPTtypes_di_unipi_it;
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS MCFClone --------------------------------*/
/*--------------------------------------------------------------------------*/
/** Testing template class for Min Cost Flow Problem solvers deriving from
 * MCFClass. It is a template deriving from Master and holding an object of
 * class Slave (both deriving from MCFClass), that does whatever Master does
 * but also does whatever Slave does. The output is that of Master. */

template<class Master, class Slave>
class MCFClone : public Master {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

   MCFClone( cIndex nmx = 0 , cIndex mmx = 0 )
             :
             Master( nmx , mmx )
   {
    SlvMCF = new Slave( nmx , mmx );

    }  // end( MCFClone )

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

   void LoadNet( cIndex nmx = 0 , cIndex mmx = 0 , cIndex pn = 0 ,
		 cIndex pm = 0 , cFRow pU = 0 , cCRow pC = 0 ,
		 cFRow pDfct = 0 , cIndex_Set pSn = 0 ,
		 cIndex_Set pEn = 0 )
   {
    Master::LoadNet( nmx , mmx , pn , pm , pU , pC , pDfct , pSn , pEn );
    SlvMCF->LoadNet( nmx , mmx , pn , pm , pU , pC , pDfct , pSn , pEn );
    }

/*--------------------------------------------------------------------------*/

   void PreProcess( void )
   {
    Master::PreProcess();
    SlvMCF->PreProcess();
    }

/*--------------------------------------------------------------------------*/

   void SetPar( int par , int val )
   {
    Master::SetPar( par , val );
    SlvMCF->SetPar( par , val );
    }

   void SetPar( int par , double val )
   {
    Master::SetPar( par , val );
    SlvMCF->SetPar( par , val );
    }

/*--------------------------------------------------------------------------*/

   void SetMCFTime( bool TimeIt = TRUE )
   {
    Master::SetMCFTime( TimeIt );
    SlvMCF->SetMCFTime( TimeIt );
    }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR SOLVING THE PROBLEM -------------------*/
/*--------------------------------------------------------------------------*/

   void SolveMCF( void )
   {
    Master::SolveMCF();
    SlvMCF->SolveMCF();
    }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

   void TimeMCF( double &t_us , double &t_ss )
   {
    Master::TimeMCF( t_us , t_ss );

    double st__us , st_ss;
    SlvMCF->TimeMCF( st__us , st_ss );

    t_us += st__us;
    t_ss += st__ss;
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   double TimeMCF( void )
   {
    return( Master::TimeMCF() + SlvMCF->TimeMCF() );
    }

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

   void ChgCosts( cCRow NCost , cIndex_Set nms = 0 ,
		  cIndex strt = 0 , Index stp = Inf<Index>() )
   {
    Master::ChgCosts( NCost , nms , strt , stp );
    SlvMCF->ChgCosts( NCost , nms , strt , stp );
    }

   void ChgCost( Index arc , cCNumber NCost )
   {
    Master::ChgCost( arc , NCost );
    SlvMCF->ChgCost( arc , NCost );
    }

/*--------------------------------------------------------------------------*/

   void ChgQCoef( cCRow NQCoef = 0 , cIndex_Set nms = 0 ,
		  cIndex strt = 0 , Index stp = Inf<Index>() )
   {
    Master::ChgQCoef( NQCoef , nms , strt , stp );
    SlvMCF->ChgQCoef( NQCoef , nms , strt , stp );
    }

   void ChgQCoef( Index arc , cCNumber NQCoef )
   {
    Master::ChgQCoef( arc , NQCoef );
    SlvMCF->ChgQCoef( arc , NQCoef );
    }

/*--------------------------------------------------------------------------*/

   void ChgUCaps( cFRow NCap  , cIndex_Set nms = 0 ,
		  cIndex strt = 0 , Index stp = Inf<Index>() )
   {
    Master::ChgUCaps( NCap  , nms , strt , stp );
    SlvMCF->ChgUCaps( NCap  , nms , strt , stp );
    }

   void ChgUCap( Index arc  , cFNumber NCap  )
   {
    Master::ChgUCap( arc  , NCap  );
    SlvMCF->ChgUCap( arc  , NCap  );
    }

/*--------------------------------------------------------------------------*/

   void ChgDfcts( cFRow NDfct , cIndex_Set nms = 0 ,
		  cIndex strt = 0 , Index stp = Inf<Index>() )
   {
    Master::ChgDfcts( NDfct , nms , strt , stp );
    SlvMCF->ChgDfcts( NDfct , nms , strt , stp );
    }

   void ChgDfct( Index node , cFNumber NDfct )
   {
    Master::ChgDfct( node , NDfct );
    SlvMCF->ChgDfct( node , NDfct );
    }

/*--------------------------------------------------------------------------*/

   void CloseArc( cIndex name )
   {
    Master::CloseArc( name );
    SlvMCF->CloseArc( name );
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   void DelNode( cIndex name )
   {
    Master::DelNode( name );
    SlvMCF->DelNode( name );
    }

/*--------------------------------------------------------------------------*/

   void OpenArc( cIndex name )
   {
    Master::OpenArc( name );
    SlvMCF->OpenArc( name );
    }

   Index AddNode( cFNumber aDfct )
   {
    SlvMCF->AddNode( aDfct );
    return( Master::AddNode( aDfct ) );
    }

/*--------------------------------------------------------------------------*/

   void ChangeArc( cIndex name , cIndex nSN = Inf<Index>() ,
		   cIndex nEN = Inf<Index>() )
   {
    Master::ChangeArc( name , nSN , nEN );
    SlvMCF->ChangeArc( name , nSN , nEN );
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   void DelArc( cIndex name )
   {
    Master::DelArc( name );
    SlvMCF->DelArc( name );
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   Index AddArc( cIndex Start , cIndex End , cFNumber aU , cCNumber aC ) 
   {
    SlvMCF->AddArc( Start , End , aU , aC );
    return( Master::AddArc( Start , End , aU , aC ) );
    }

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

   ~MCFClone()
   {
    delete SlvMCF;
    }

/*--------------------------------------------------------------------------*/
/*------------------------ PUBLIC DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

   Slave *SlvMCF;     /**< pointer to the "slave" object; it is public in
			 order to allow calling the methods of the
			 specialized interface of the slave class */

 };   // end( class MCFClone )

/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 };  // end( namespace MCFClass_di_unipi_it )
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* MCFClone.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File MCFClone.h ----------------------------*/
/*--------------------------------------------------------------------------*/
