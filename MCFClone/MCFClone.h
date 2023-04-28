/*--------------------------------------------------------------------------*/
/*-------------------------- File MCFClone.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Testing template class for Min Cost Flow Problem solvers deriving from
 * MCFClass. It is a template deriving from Master and holding an object of
 * class Slave (both deriving from MCFClass), that does whatever Master does
 * but also does whatever Slave does. The output is that of Master.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
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

namespace MCFClass_di_unipi_it
{
 using namespace OPTtypes_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS MCFClone --------------------------------*/
/*--------------------------------------------------------------------------*/
/** MCFClone derives from MCFClass and is a template class for testing Min
 * Cost Flow Problem solvers deriving from MCFClass. It is a template
 * deriving from Master and holding an object of class Slave, both deriving
 * from MCFClass, that does whatever Master does but also does whatever
 * Slave does. The output is that of Master. Note that since it derives from
 * Master, all the methods for readins stuff need not be defined since these
 * of Master are automatically available and do rhe right thing. Yet, one
 * may want to re-define some of them to check that the results agree. */

template< class Master , class Slave >
class MCFClone : public Master {
 
 static_assert( std::is_base_of< MCFClass , T >::value );
 static_assert( std::is_base_of< MCFClass , T >::value );

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

   MCFClone( cIndex nmx = 0 , cIndex mmx = 0 ) : Master( nmx , mmx ) {
    SlvMCF = new Slave( nmx , mmx );
    }

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

   void LoadNet( Index nmx = 0 , Index mmx = 0 , Index pn = 0 ,
		 Index pm = 0 , cFRow pU = 0 , cCRow pC = 0 ,
		 cFRow pDfct = 0 , cIndex_Set pSn = 0 ,
		 cIndex_Set pEn = 0 ) override {
    Master::LoadNet( nmx , mmx , pn , pm , pU , pC , pDfct , pSn , pEn );
    SlvMCF->LoadNet( nmx , mmx , pn , pm , pU , pC , pDfct , pSn , pEn );
    }

/*--------------------------------------------------------------------------*/

   void PreProcess( void ) override {
    Master::PreProcess();
    SlvMCF->PreProcess();
    }

/*--------------------------------------------------------------------------*/

   void SetPar( int par , int val ) override {
    Master::SetPar( par , val );
    SlvMCF->SetPar( par , val );
    }

   void SetPar( int par , double val ) override {
    Master::SetPar( par , val );
    SlvMCF->SetPar( par , val );
    }

/*--------------------------------------------------------------------------*/

   void SetMCFTime( bool TimeIt = true ) override {
    Master::SetMCFTime( TimeIt );
    SlvMCF->SetMCFTime( TimeIt );
    }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR SOLVING THE PROBLEM -------------------*/
/*--------------------------------------------------------------------------*/

   void SolveMCF( void ) override {
    Master::SolveMCF();
    SlvMCF->SolveMCF();
    }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

   void TimeMCF( double &t_us , double &t_ss ) override {
    Master::TimeMCF( t_us , t_ss );

    double st__us , st_ss;
    SlvMCF->TimeMCF( st__us , st_ss );

    t_us += st__us;
    t_ss += st__ss;
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   double TimeMCF( void ) override {
    return( Master::TimeMCF() + SlvMCF->TimeMCF() );
    }

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

   void ChgCosts( cCRow NCost , cIndex_Set nms = 0 ,
		  Index strt = 0 , Index stp = Inf< Index >() ) override {
    Master::ChgCosts( NCost , nms , strt , stp );
    SlvMCF->ChgCosts( NCost , nms , strt , stp );
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   void ChgCost( Index arc , CNumber NCost ) {
    Master::ChgCost( arc , NCost );
    SlvMCF->ChgCost( arc , NCost );
    }

/*--------------------------------------------------------------------------*/

   void ChgQCoef( cCRow NQCoef = 0 , cIndex_Set nms = 0 ,
		  Index strt = 0 , Index stp = Inf< Index >() ) override {
    Master::ChgQCoef( NQCoef , nms , strt , stp );
    SlvMCF->ChgQCoef( NQCoef , nms , strt , stp );
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   void ChgQCoef( Index arc , CNumber NQCoef ) override {
    Master::ChgQCoef( arc , NQCoef );
    SlvMCF->ChgQCoef( arc , NQCoef );
    }

/*--------------------------------------------------------------------------*/

   void ChgUCaps( cFRow NCap  , cIndex_Set nms = 0 ,
		  Index strt = 0 , Index stp = Inf< Index >() ) override {
    Master::ChgUCaps( NCap  , nms , strt , stp );
    SlvMCF->ChgUCaps( NCap  , nms , strt , stp );
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   void ChgUCap( Index arc  , FNumber NCap  ) override {
    Master::ChgUCap( arc  , NCap  );
    SlvMCF->ChgUCap( arc  , NCap  );
    }

/*--------------------------------------------------------------------------*/

   void ChgDfcts( cFRow NDfct , cIndex_Set nms = 0 ,
		  Index strt = 0 , Index stp = Inf< Index >() ) override {
    Master::ChgDfcts( NDfct , nms , strt , stp );
    SlvMCF->ChgDfcts( NDfct , nms , strt , stp );
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   void ChgDfct( Index node , FNumber NDfct ) override {
    Master::ChgDfct( node , NDfct );
    SlvMCF->ChgDfct( node , NDfct );
    }

/*--------------------------------------------------------------------------*/

   void CloseArc( Index name ) override {
    Master::CloseArc( name );
    SlvMCF->CloseArc( name );
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   void DelNode( Index name ) override {
    Master::DelNode( name );
    SlvMCF->DelNode( name );
    }

/*--------------------------------------------------------------------------*/

   void OpenArc( Index name ) override {
    Master::OpenArc( name );
    SlvMCF->OpenArc( name );
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   Index AddNode( FNumber aDfct ) override {
    SlvMCF->AddNode( aDfct );
    return( Master::AddNode( aDfct ) );
    }

/*--------------------------------------------------------------------------*/

   void ChangeArc( Index name , Index nSN = Inf< Index >() ,
		   Index nEN = Inf< Index >() ) override {
    Master::ChangeArc( name , nSN , nEN );
    SlvMCF->ChangeArc( name , nSN , nEN );
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   void DelArc( Index name ) override {
    Master::DelArc( name );
    SlvMCF->DelArc( name );
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   Index AddArc( Index Start , Index End , FNumber aU , CNumber aC )
    override {
    SlvMCF->AddArc( Start , End , aU , aC );
    return( Master::AddArc( Start , End , aU , aC ) );
    }

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

   ~MCFClone() { delete SlvMCF; }

/*--------------------------------------------------------------------------*/
/*------------------------ PUBLIC DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

   Slave * SlvMCF;     /**< pointer to the "slave" object; it is public in
			* order to allow calling the methods of the
			* specialized interface of the slave class */

/*--------------------------------------------------------------------------*/

 };   // end( class MCFClone )

/*--------------------------------------------------------------------------*/

};  // end( namespace MCFClass_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* MCFClone.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File MCFClone.h ----------------------------*/
/*--------------------------------------------------------------------------*/
