/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2016-2017
                  
                Copyright (C) 2016 Stefano Brambilla
======================================================================*/
/*! 
  @file   transport3d1d.cpp
  @author Stefano Brambilla <s.brambilla93@gmail.com>
  @date   September 2016.
  @brief  Definition of the aux class for algorithm description strings for transport problem.
 */
 

 
/** @defgroup input User-defined parameters  */

#ifndef M3D1D_DESCR3D1D_TRANSP_HPP_
#define M3D1D_DESCR3D1D_TRANSP_HPP_

#include <string>

namespace getfem {

//! Class to import the descriptors of the coupled 3D/1D solver
/*!
	\ingroup input
 */
struct descr3d1d_transp {

	//! Absolute path to the vessel mesh file
	std::string MESH_FILEV;
	//! Identifier of tissue concentration's FEM type
	std::string FEM_TYPET_C;
	//! Identifier of vessel concentration's FEM type
	std::string FEM_TYPEV_C;
	//! Identifier of vessel integration method type
	std::string IM_TYPEV;
	//! Output directory for transport problem
	std::string OUTPUT;
	//! Flag to make the problem stationary
	size_type STATIONARY;	
	// Solver information
	//! Identifief of the monolithic solver for transport problem
	std::string SOLVE_METHOD;
	//! Maximum number of iterations (iterative solvers)
	size_type   MAXITER;
	//! Mamimum residual (iterative solvers)
	scalar_type RES; 
	//! Number of target points for the tissue-to-vessel average
	size_type   NInt;
	// Utils
	//! File .param
	ftool::md_param FILE_;
	//! Import algorithm specifications from file .param
	void import(ftool::md_param & fname) 
	{
		FILE_ = fname;
		
		MESH_FILEV  = FILE_.string_value("MESH_FILEV_TRANSP","1D points file for transport problem");
		FEM_TYPET_C   = FILE_.string_value("FEM_TYPET_C","FEM 3D tissue - concentration");
		
		FEM_TYPEV_C   = FILE_.string_value("FEM_TYPEV_C","FEM 1D vessel - concentration");
		IM_TYPEV 	= FILE_.string_value("IM_TYPEV_TRANSP","Name of integration method");
		SOLVE_METHOD = FILE_.string_value("SOLVE_METHOD", "Monolithic Solver"); 
		if (SOLVE_METHOD != "SuperLU") { // iterative solver
			MAXITER  = FILE_.int_value("MAXITER", "Max number of sub-iterations");
			RES = FILE_.real_value("RES"); if (RES == 0.) RES = 2.0e-10;
		}
		NInt = size_type(FILE_.int_value("NInt", "Node numbers on the circle for the nonlocal term"));  
		OUTPUT = FILE_.string_value("OUTPUT","Output Directory");
		STATIONARY = size_type(FILE_.int_value("STATIONARY", "Flag to make the problem stationary")); 
	}

	//! Overloading of the output operator
	friend std::ostream & operator << (
		std::ostream & out, const descr3d1d_transp & descr
		)
	{ 
		cout << "---- TRANSPORT PROBLEM DESCRIPTORS--------------------------" << endl;
		
		cout << " FEM TYPE  3D concentration     : " << descr.FEM_TYPET_C   << endl;
		cout << " FEM TYPE  1D concentration     : " << descr.FEM_TYPEV_C   << endl;
		cout << "--------------------------------------------------" << endl;

		return out;            
	}

}; /* end of class */

} /* end of namespace */

#endif
