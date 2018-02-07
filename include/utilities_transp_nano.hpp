/* -*- c++ -*- (enables emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2014-2015
                  
                Copyright (C) 2017 Annagiulia Tiozzo
======================================================================*/
/*! 
  @file   utilities_transp_nano.hpp
  @author Annagiulia Tiozzo <annagiulia.tiozzo@mail.polimi.com>
  @date   March 2017.
  @brief  Declaration of some miscelleanous auxiliary routines.
 */
#ifndef M3D1D_UTILITIES_TRANSP_NANO_HPP_
#define M3D1D_UTILITIES_TRANSP_NANO_HPP_

#include <getfem/getfem_assembling.h> 
#include <getfem/getfem_mesh.h>
#include <gmm/gmm.h>
#include <math.h>

namespace getfem {
//! Compute the Wall Shear Stress: 	
//! @f$ W = 4 vel / Rad @f$
/*
 @param WSR  wall shear stress
 @param vel velocity 
 @param Rad radius of the branch
*/

template <typename VEC>
void wall_shear_stress	
	(VEC & W,
	const scalar_type & Rad,
	VEC & vel
	)
{	
	int lengthW = gmm::mat_nrows(gmm::col_vector(W));
	int lengthvel= gmm::mat_nrows(gmm::col_vector(vel));

	gmm::add(vel, W);
	gmm::scale(W, 4.0/Rad); 
	//for(int i=0; i<lengthW; i++){
	//std::cout<<"wWSR al nodo: "<<i<<" = "<<W[i]<<std::endl;
	//}
	//std::cout<<"raggio per iol WSR: "<<Rad<<std::endl;
};


//Compute the probability of adhesion P_a (s): 	
//! @f$ P=c1*exp(c2*WSR) @f$
//deve restituire un vettore
/*
 @param P 	probability of adhesion
 @param WSR  	wall shear stress
 @param c1 	la costante fuori dall'esponenziale (m_l K_a^0 alpha2 pi r_0^2)
 @param c2 	la costante dentrol'esponenziale (-beta * mu/alpha2)
*/

template <typename VEC>
void probability_adhesion0
	(VEC & P,
	VEC & W,
	const scalar_type & ml,
	const scalar_type & Ka,
	const scalar_type & alpha,
	const scalar_type & r0,
	const scalar_type & beta,
	const scalar_type & mu,
	const scalar_type & mr
	)

{ 	int lengthP= gmm::mat_nrows(gmm::col_vector(P));
	int lengthW= gmm::mat_nrows(gmm::col_vector(W));
	GMM_ASSERT1(lengthP==lengthW, "impossible to operate on the vectors: different lengths");

	vector_type temporary(lengthP); gmm::clear(temporary); 
	vector_type temporary2(lengthP); gmm::clear(temporary2); 
	gmm::add(W,temporary);

	gmm::scale(temporary, -(beta*mu/alpha));

	for(int i=0; i<lengthP; i++){
	temporary2[i]=exp(temporary[i]);
	}

	//gmm::scale(temporary2, (ml*Ka*alpha*pi*pow(r0,2.0))); //In Pa serve dp/2 al posto di r0
	gmm::scale(temporary2, (ml*Ka*mr*pi*pow(r0,2.0)));
	gmm::add(temporary2,P);
	gmm::clear(temporary);
	gmm::clear(temporary2);	


};

//Compute the probability of adhesion P_a (s): 	
//! come funzione interpolata linearmente dei dati di Alessandro
//deve restituire un vettore
/*
 @param P 	probability of adhesion
 @param Re  	Reynolds number
 @param Pa_min 	P_a(Re_min)
 @param Pa_med 	P_a(Re_med)
 @param Pa_max 	P_a(Re_max)
*/

template <typename VEC>
void probability_adhesion
	(VEC & P,
	VEC & Re,
	scalar_type Pa_min,
	scalar_type Pa_med,
	scalar_type Pa_max
	)

{ 	int lengthP= gmm::mat_nrows(gmm::col_vector(P));
	int lengthRe= gmm::mat_nrows(gmm::col_vector(Re));
	GMM_ASSERT1(lengthP==lengthRe, "impossible to operate on the vectors: different lengths");

	scalar_type Re_min=0.01;
	scalar_type Re_med=0.1;
	scalar_type Re_max=1.0;


	for(int i=0; i<lengthP; i++){
		if(Re[i]<=Re_min)
			P[i]=1.0; //Questo non è vero nel caso di Rhol=0.30, sistemare
		if(Re[i]<=Re_med && Re[i]>Re_min)
			P[i]=(Pa_med-Pa_min)/(Re_med-Re_min)*(Re[i]-Re_min)+Pa_min;
		if(Re[i]<=Re_max && Re[i]>Re_med)
			P[i]=(Pa_max-Pa_med)/(Re_max-Re_med)*(Re[i]-Re_med)+Pa_med;
		if(Re[i]>Re_max)
			P[i]=0.0; //Vedere se sistemare questo a valori diversi da zero
	}




};

//Compute the Reynolds number: 	
//! @f$ Re = vel \rho Rad / \mu @f$

/*
 @param Re 	Reynolds number
 @param vel  	velocity
 @param mu 	viscosity
 @param rho 	density
*/
template <typename VEC>
void reynolds
	(VEC & Re,
	 VEC & vel,
	const scalar_type & Rad,
	const scalar_type & mu,
	const scalar_type & rho
	)
{
	//std::cout<<"mu: "<<mu<<std::endl;
	//std::cout<<"rho: "<<rho<<std::endl;
	//std::cout<<"Raggio: "<<Rad<<std::endl;
	
	int lengthRe= gmm::mat_nrows(gmm::col_vector(Re));
	int lengthVel= gmm::mat_nrows(gmm::col_vector(vel));

	//for(int i=0; i<lengthRe; i++){
	//std::cout<<"Velocita al nodo: "<<i<<" = "<<vel[i]<<std::endl;
	//}
	
	gmm::add(vel, Re);
	gmm::scale(Re, rho*Rad/mu);

	//for(int i=0; i<lengthRe; i++){
	//std::cout<<"Reynolds al nodo: "<<i<<" = "<<Re[i]<<std::endl;
	//}

};
/*
template <typename VEC>
void compute_Psi(void)
 	{
    	size_type nb_dof_v = mf_uv.nb_dof();
    	std::vector<scalar_type> temp(nb_dof_v, 0.0);
    	scalar_type cost;
    	cost = DT*constTempo/2;   //dimensional time
    	gmm::copy(Psi, Psiold);
    	gmm::copy(gmm::scaled(Cv, cost), temp);
    	gmm::add(gmm::scaled(Cvold, cost), temp);

    	for(size_type i=0; i<nb_dof_v; i++) {
      	Psi[i] = constPsi*vectPiGreco[i]*temp[i];
      	}
 	   gmm::add(Psiold, Psi);
 	}
*/

} /* end of namespace */

	//MAGARI CI SERVE PER SCORRERE GLI ELEMENTI DI UNA SINGOLA REGIONE
	//for (size_type b=0; b<nb_branches; ++b)
	//	for (getfem::mr_visitor mrv(mf_data.linked_mesh().region(b)); !mrv.finished(); ++mrv)
	//		for (auto i : mf_data.ind_basic_dof_of_element(mrv.cv()))
	//			Radius[i] = Rdata[b];


#endif
