 /* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2016-2017
                  
                Copyright (C) 2016 Annagiulia Tiozzo
======================================================================*/
/*! 
  @file   assembling1d_transp_nano.hpp
  @author Annagiulia Tiozzo <annagiulia.tiozzo@mail.polimi.com>
  @date   September 2016.
  @brief  Miscelleanous assembly routines for the 1D transport problem.
 */
/** @defgroup asm Assembly routines */

#ifndef M3D1D_ASSEMBLING_1D_TRANSP_NANO_HPP_
#define M3D1D_ASSEMBLING_1D_TRANSP_NANO_HPP_
#include <defines.hpp>
#include <utilities.hpp>

namespace getfem { //NANO

/*! Build the adhesion term
    @f$ A=\int_{\Lambda}  M c_v b_v~d\sigma@f$
 */
/*!
	@param A         Computed mass matrix
	@param mim       The integration method to be used
	@param mf_c      The finite element method for the concentration
	@param mf_u  	 The finite element method for the adhesion parameter M
	@param adh	 The value for the adhesion parameter M
	@param rg        The region where to integrate

	@ingroup asm
 */

template<typename MAT, typename VEC>
void 
asm_network_nano_transp
	(MAT & A,  
	 const mesh_im & mim,
	 const mesh_fem & mf_c,
	 const mesh_fem & mf_u,
	 const VEC & adh,
	 const mesh_region & rg = mesh_region::all_convexes()
	 ) 		
{
	GMM_ASSERT1(mf_c.get_qdim() == 1 ,
		"invalid data mesh fem (Qdim=1 required)");
	
	generic_assembly 
	assem("m=data$1(#2);"
		  "t=comp(Base(#2).Base(#1).Base(#1));"
		  "M$1(#1,#1)+=t(k,:,:).m(k);"); 
	assem.push_mi(mim);
	assem.push_mf(mf_c);
	assem.push_mf(mf_u);
	assem.push_data(adh);
	assem.push_mat(A);
	assem.assembly(rg);

}


template<typename MAT, typename VEC>
void 
asm_network_nano_transp
	(MAT & A,  
	 const mesh_im & mim,
	 const mesh_fem & mf_c,
	 const VEC & adh,
	 const mesh_region & rg = mesh_region::all_convexes()
	 ) 		
{
	GMM_ASSERT1(mf_c.get_qdim() == 1 ,
		"invalid data mesh fem (Qdim=1 required)");
	
	generic_assembly 
	assem("m=data$1(#1);"
		  "t=comp(Base(#1).Base(#1).Base(#1));"
		  "M$1(#1,#1)+=t(k,:,:).m(k);"); 
	assem.push_mi(mim);
	assem.push_mf(mf_c);
	assem.push_data(adh);
	assem.push_mat(A);
	assem.assembly(rg);

}


} /* end of namespace */

#endif
