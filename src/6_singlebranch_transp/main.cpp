/* -*- c++ -*- (enables emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2014-2015
                  
                Copyright (C) 2015 Domenico Notaro
======================================================================*/
/*! 
  @file   main.cpp
  @author Domenico Notaro <domenico.not@gmail.com>
  @date   April 2015.
  @brief  Main program for test simulations.
  @details
    We solve the coupled 3D/1D problem of fluid exchange between a 1D 
    network \Lambda and the 3D interstitial tissue \Omega
    
    *****************************************
      Benchmark : single-vessel network 
      Mixed finite elements approximation
      Monolithic resolution by SuperLU 3.0
    *****************************************
    
	See Section "Code verification: test-cases"
 */
#include <iostream>
#include <getfem/bgeot_config.h> // for FE_ENABLE_EXCEPT
#include <problem3d1d.hpp>
#include <transport3d1d.hpp> 
using namespace getfem;

/* default 4 Byte integer types */
#ifndef APPL_INT
#define APPL_INT int
#endif
/* end of integer.h */

//! main program
int main(int argc, char *argv[]) 
{

	GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
	FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

	try {   
		// Declare a new problem 
		getfem::transport3d1d p; 

		/// fluid problem: velocity field and pressure
		// Initialize the problem

		p.problem3d1d::init(argc, argv);
//FEDE
/*
		if(!p.LINEAR_LYMPH())
			{
			// Build the monolithic system
			p.problem3d1d::assembly_fixpoint();
			// Solve the problem
			if (!p.problem3d1d::solve_fixpoint()) GMM_ASSERT1(false, "solve procedure has failed");
			}
		else
			{
			// Build the monolithic system
			p.problem3d1d::assembly();
			// Solve the problem
			if (!p.problem3d1d::solve()) GMM_ASSERT1(false, "solve procedure has failed");
			}
		// Save results in .vtk format
		p.problem3d1d::export_vtk();
		// Display some global results: mean pressures, total flow rate
		std::cout << "--- FINAL RESULTS -------------------------" << std::endl; 
		std::cout << "  Pt average            = " << p.mean_pt()   << std::endl;
		std::cout << "  Pv average            = " << p.mean_pv()   << std::endl;
		std::cout << "  Network-to-Tissue TFR = " << p.flow_rate() << std::endl;
		std::cout << "  Lymphatic FR          = " << p.lymph_flow_rate() << std::endl;
		std::cout << "  FR from the cube      = " << p.cube_flow_rate() << std::endl;
		std::cout << "-------------------------------------------" << std::endl; 	
*/
		//transport problem: concentration  
		
		//initialize 
		p.init(argc, argv);
		//assemble        
		p.assembly();    
		//solve     
		// if (!p.solve()) GMM_ASSERT1(false, "solve procedure has failed");  // the export is in the solve at each time step 
	 	 if (!p.solve_samg()) GMM_ASSERT1(false, "solve procedure has failed");  // the export is in the solve at each time step 
				      
		   
	}

	GMM_STANDARD_CATCH_ERROR;

	return 0; 
	
} /* end of main program */

