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
  @brief  Definition of the main class for the 3D/1D coupled transport problem.
  */

#include <transport3d1d.hpp>
#include <assembling1d_transp_nano.hpp>

//SAMG
#define CSC_INTERFACE 
//#define SPARSE_INTERFACE
//#define CSR_INTERFACE 
// ------------------------------------
#define DIRECT_SOLVER 
//#define AMG_STAND_ALONE
//#define AMG_ACCELERATED

#ifdef USE_SAMG
#include "samg.h"
#endif

/* default 4 Byte integer types */
#ifndef APPL_INT
#define APPL_INT int
#endif
/* end of integer.h */


namespace getfem {

	void transport3d1d::init (int argc, char *argv[]) 
	{
		std::cout << "initialize transport problem..."<<std::endl<<std::endl;

		import_data();
		build_mesh();
		set_im_and_fem();
		build_param();
		build_vessel_boundary();
		build_tissue_boundary();

	}; // end of init


	// Aux methods for init

	//! Import algorithm specifications
	void transport3d1d::import_data(void)
	{
		std::cout<<"init part 1: import data!......" <<std::endl;
#ifdef M3D1D_VERBOSE_
		cout << "Importing descriptors for tissue and vessel problems ..." << endl;
#endif
		descr_transp.import(PARAM);
#ifdef M3D1D_VERBOSE_
		cout << descr_transp;
#endif


	};



	//! Import mesh for tissue (3D) and vessel (1D)  
	void transport3d1d::build_mesh(void){

		//but, in order to have the boundary conditions for the nodes
		//we need to build again the 1D mesh from another pts file
		mesht.clear();
		bool test = 0;
		test = PARAM.int_value("TEST_GEOMETRY");
		if(test==0){
#ifdef M3D1D_VERBOSE_
			cout << "Importing the 3D mesh for the tissue ...  "   << endl;
#endif
			import_msh_file(descr.MESH_FILET, mesht);
		}else{
#ifdef M3D1D_VERBOSE_
			cout << "Building the regular 3D mesh for the tissue ...  "   << endl;
#endif
			string st("GT='" + PARAM.string_value("GT_T") + "'; " +
					"NSUBDIV=" + PARAM.string_value("NSUBDIV_T") + "; " +  
					"ORG=" + PARAM.string_value("ORG_T") + "; " +  
					"SIZES=" + PARAM.string_value("SIZES_T") + "; " +  
					"NOISED=" + PARAM.string_value("NOISED_T")); 
			cout << "mesht description: " << st << endl;
			regular_mesh(mesht, st);
		}

#ifdef M3D1D_VERBOSE_
		cout << "Importing the 1D mesh for the vessel (transport problem)... "   << endl;
#endif
		std::ifstream ifs(descr_transp.MESH_FILEV);
		meshv.clear();
		//	mesh meshv_transp;
		//	vector_size_type nb_vertices_transp;
		GMM_ASSERT1(ifs.good(), "impossible to read from file " << descr_transp.MESH_FILEV);
		import_pts_file(ifs, meshv, BCv_transp, nb_vertices, descr.MESH_TYPEV);
		nb_branches = nb_vertices.size();
		ifs.close();

		/*
		   cout<<"BC per il problema di trasporto"<<endl;
		   for (size_type bc=0; bc < BCv_transp.size(); bc++) {
		   cout<<"bc:"<<bc<<endl;
		   cout <<BCv_transp[bc]<<endl; 
		   }
		   cout<<"BC per il problema di stokes"<<endl;
		   for (size_type bc=0; bc < BCv.size(); bc++) {
		   cout<<"bc:"<<bc<<endl;
		   cout <<BCv[bc]<<endl;
		   }
		   */

	};
	//! Set finite elements methods and integration methods 
	void transport3d1d::set_im_and_fem(void)
	{
		std::cout<<"init part 2: set fem methods!......" <<std::endl;


#ifdef M3D1D_VERBOSE_
		cout << "Setting FEMs for tissue and vessel problems ..." << endl;
#endif


		pfem pf_Ct = fem_descriptor(descr_transp.FEM_TYPET_C);
		pfem pf_Cv = fem_descriptor(descr_transp.FEM_TYPEV_C);

#ifdef M3D1D_VERBOSE_
		cout << "Setting IMs and FEMs for tissue ..." << endl;
#endif


		mf_Ct.set_finite_element(mesht.convex_index(), pf_Ct);

#ifdef M3D1D_VERBOSE_
		cout << "Setting IMs and FEMs for vessel branches ..." << endl;
#endif

		mf_Cv.set_finite_element(meshv.convex_index(), pf_Cv);


#ifdef M3D1D_VERBOSE_
		cout << "Setting FEM dimensions for tissue and vessel problems ..." << endl;
#endif

		dof_transp.set(mf_Ct, mf_Cv);
#ifdef M3D1D_VERBOSE_
		cout << std::scientific << dof_transp;
#endif


		//HO DOVUTO CANCELLARE MESHV!!!!!!!!!!!!!!
		// quindi metodi fem e di integrazione definiti sulla rete non valgono più! 
		// li devo caricare di nuovo!
		mimv.clear();
		mf_Uvi.clear();
		mf_Pv.clear();
		mf_coefv.clear();
		mf_coefvi.clear();

		mimt.clear();
		mf_Ut.clear();
		mf_Pt.clear();
		mf_coeft.clear();


		problem3d1d::set_im_and_fem();
		/*
#ifdef M3D1D_VERBOSE_
cout << "Setting IMs for tissue and vessel problems ..." << endl;
#endif
pintegration_method pim_v = int_method_descriptor(descr.IM_TYPEV);
mimv.set_integration_method(meshv.convex_index(), pim_v);

#ifdef M3D1D_VERBOSE_
cout << "Setting FEMs for tissue and vessel problems ..." << endl;
#endif
bgeot::pgeometric_trans pgt_v = bgeot::geometric_trans_descriptor(descr.MESH_TYPEV);
pfem pf_Uv = fem_descriptor(descr.FEM_TYPEV);
pfem pf_Pv = fem_descriptor(descr.FEM_TYPEV_P);
pfem pf_coefv = fem_descriptor(descr.FEM_TYPEV_DATA);

#ifdef M3D1D_VERBOSE_
cout << "Setting IMs and FEMs for vessel branches ..." << endl;
#endif
mf_Uvi.reserve(nb_branches);
mf_coefvi.reserve(nb_branches);
for(size_type i=0; i<nb_branches; ++i){

mesh_fem mf_tmp(meshv);
mf_tmp.set_finite_element(meshv.region(i).index(), filepf_coefv);
mf_coefvi.emplace_back(mf_tmp);
mf_tmp.clear();

mf_tmp.set_finite_element(meshv.region(i).index(), pf_Uv);
mf_Uvi.emplace_back(mf_tmp);
mf_tmp.clear();
}
mf_Pv.set_finite_element(meshv.convex_index(), pf_Pv);
mf_coefv.set_finite_element(meshv.convex_index(), pf_coefv);

*/


		};


//! Build problem parameters
void transport3d1d::build_param(void)
{
	std::cout<<"init part 3: build dimensionless parameters!" <<std::endl;

#ifdef M3D1D_VERBOSE_
	cout << "Building parameters for tissue and vessel problems ..." << endl;
#endif
	//FEDE PER FIXED SOURCE: in input devo dare anche mf_Cv perchè ho modificato la funzione in param3d1d_transp.hpp
	param_transp.build(PARAM, mf_coeft, mf_coefv, mf_Cv);
#ifdef M3D1D_VERBOSE_
	cout << param_transp ;
#endif
	cout<<param_transp;
};



	void
transport3d1d::build_tissue_boundary (void) 
{
#ifdef M3D1D_VERBOSE_
	cout << "Building tissue boundary ..." << endl;
#endif
	BCt_transp.clear();
	BCt_transp.reserve(2*DIMT);
	// Parse BC data
	string label_in = PARAM.string_value("BClabel_transp", "Array of tissue boundary labels");
	string value_in = PARAM.string_value("BCvalue_transp", "Array of tissue boundary values");
	vector<string> labels = split(label_in, ' ');
	vector<string> values = split(value_in, ' ');
	GMM_ASSERT1(labels.size()==2*DIMT, "wrong number of BC labels");
	GMM_ASSERT1(values.size()==2*DIMT, "wrong number of BC values");
	for (unsigned f=0; f<2*DIMT; ++f) {
		BCt_transp.emplace_back(labels[f], std::stof(values[f]), 0, f);
#ifdef M3D1D_VERBOSE_
		cout << "  face " << f << " : " << BCt_transp.back() << endl;
#endif
	} 

	for (size_type bc=0; bc < BCt_transp.size(); bc++)
		cout<<BCt_transp[bc]<<endl;

	// Build mesht regions
	mesh_region border_faces;
	outer_faces_of_mesh(mesht, border_faces);

	for (mr_visitor i(border_faces); !i.finished(); ++i) {

		assert(i.is_face());

		// Unit outward normal : used to identify faces
		//! \todo Use getfem 5.0's function select_faces_of_normal?
		base_node un = mesht.normal_of_face_of_convex(i.cv(), i.f());
		un /= gmm::vect_norm2(un);

		if (gmm::abs(un[0] + 1.0) < 1.0E-7)      // back
			mesht.region(0).add(i.cv(), i.f());
		else if (gmm::abs(un[0] - 1.0) < 1.0E-7) // front
			mesht.region(1).add(i.cv(), i.f());
		else if (gmm::abs(un[1] + 1.0) < 1.0E-7) // left
			mesht.region(2).add(i.cv(), i.f());
		else if (gmm::abs(un[1] - 1.0) < 1.0E-7) // right
			mesht.region(3).add(i.cv(), i.f());
		else if (gmm::abs(un[2] + 1.0) < 1.0E-7) // bottom
			mesht.region(4).add(i.cv(), i.f());
		else if (gmm::abs(un[2] - 1.0) < 1.0E-7) // top
			mesht.region(5).add(i.cv(), i.f());

	}

}
	void 
transport3d1d::build_vessel_boundary(void)
{
#ifdef M3D1D_VERBOSE_
	cout << "Building vessel boundary ..." << endl;
#endif
	try {

		dal::bit_vector junctions; // global idx of junctions vertices in meshv
		dal::bit_vector extrema;   // global idx of extreme vertices in meshv

		Jv.clear();
		nb_extrema=0; 
		nb_junctions=0;

		size_type fer = nb_branches; // first empty region
		GMM_ASSERT1(meshv.has_region(fer)==0, 
				"Overload in meshv region assembling!");

		// List all the convexes
		dal::bit_vector nn = meshv.convex_index();
		bgeot::size_type cv;
		for (cv << nn; cv != bgeot::size_type(-1); cv << nn) {

			bgeot::pconvex_structure cvs = meshv.structure_of_convex(cv);
			if (cvs->nb_points()>2) 
				cerr << "Error: convex #" << cv << "has more than 2 vertices!" << endl;
			if (cvs->nb_faces()>2)  
				cerr << "Error: convex #" << cv << "has more than 2 faces!" << endl;

			// Build regions for BCs and junctions
			// Global idx of mesh vertices
			size_type i0 = meshv.ind_points_of_convex(cv)[cvs->ind_points_of_face(1)[0]];
			size_type i1 = meshv.ind_points_of_convex(cv)[cvs->ind_points_of_face(0)[0]];
			// Identify vertex type
			if (meshv.convex_to_point(i0).size()==1){ /* inflow extremum */
				// Update information
				extrema.add(i0);
				nb_extrema++;
				// Build a new region made by a single face
				GMM_ASSERT1(meshv.has_region(fer)==0, 
						"Overload in meshv region assembling!");
				meshv.region(fer).add(cv, 1);
				// Store the current index and then update it
				size_type bc = 0; 
				bool found = false;
				while (!found && (bc<BCv_transp.size())) {
					found = (i0 == BCv_transp[bc].idx);
					if (!found) bc++;
				}
				GMM_ASSERT1(found=true, "Miss a boundary node in BCv list!");
				BCv_transp[bc].rg = fer; 
				fer++;
				// Store the containing branch index
				size_type branch = 0; 
				bool contained = false;
				while (!contained && branch<nb_branches ) {
					contained = meshv.region(branch).is_in(cv);
					if (!contained) branch++;
				}
				GMM_ASSERT1(contained=true, "No branch region contains node i0!");
				BCv_transp[bc].branches.emplace_back(branch); 
			}
			else if (meshv.convex_to_point(i0).size()==2){ /* trivial inflow junction */
				// DO NOTHING
			}
			else if (meshv.convex_to_point(i0).size()>=2){ /* non-trivial inflow junction */
				// Check if jucntion has been already stored, 
				// if not add to the junction list (J) and build a new region
				dal::bit_vector tmp; tmp.add(i0);
				if(!junctions.contains(tmp)){
					// Store the junction vertex
					junctions.add(i0);
					nb_junctions++;
					GMM_ASSERT1(meshv.has_region(fer)==0, 
							"Overload in meshv region assembling!");
					// Build a new region with idx "first empty region"
					meshv.region(fer).add(cv, 1); // single-face region
					// Create a new junction node
					Jv.emplace_back("JUN", 0, i0, fer);
					fer++;
				}
				// Search for index of containing branch (\mathcal{P}^{in}_j)
				size_type branch = 0; 
				bool contained = false;
				while (!contained && branch<nb_branches ) {
					contained = meshv.region(branch).is_in(cv);
					if (!contained) branch++;
				}
				GMM_ASSERT1(contained=true, "No branch region contains node i0!");
				// Add the inflow branch (to the right junction node)
				size_type jj = 0;
				bool found = false;
				while (!found && jj < nb_junctions){
					found = (i0 == Jv[jj].idx);
					if (!found) jj++;
				}
				//cout << "Branch -" << branch << " added to junction " << jj << endl;
				Jv[jj].value += param.R(mimv, branch);
				Jv[jj].branches.emplace_back(-branch);
				GMM_ASSERT1(branch>0, 
						"Error in network labeling: -0 makes no sense");
			}

			if (meshv.convex_to_point(i1).size()==1){ 
				size_type bc = 0; 
				bool found = false;
				while (!found && (bc<BCv_transp.size())) {
					found = (i1 == BCv_transp[bc].idx);
					if (!found) bc++;
				}
				if (found){ /* outlow extremum */
					extrema.add(i1); 
					nb_extrema++; 
					// Build a new region made by a single face
					GMM_ASSERT1(meshv.has_region(fer)==0, 
							"Overload in meshv region assembling!");
					meshv.region(fer).add(cv, 0);
					// Store the current index and then update it
					BCv_transp[bc].value *= +1.0;
					BCv_transp[bc].rg = fer; 
					fer++;
					// Store the containing branch index
					size_type branch = 0; 
					bool contained = false;
					while (!contained && branch<nb_branches ) {
						contained = meshv.region(branch).is_in(cv);
						if (!contained) branch++;
					}
					GMM_ASSERT1(contained=true, "No branch region contains node i1!");
					BCv_transp[bc].branches.emplace_back(branch); 
				}
				/*else { // interior -> Mixed point
				// "MIX" label via post-processing
				// Build a new region made by a single face
				GMM_ASSERT1(meshv.has_region(fer)==0, 
				"Overload in meshv region assembling!");
				meshv.region(fer).add(cv, 0);
				BCv_transp.emplace_back("MIX", 0.0, i1, fer);
				fer++;
				// Store the containing branch index
				size_type branch = 0; 
				bool contained = false;
				while (!contained && branch<nb_branches ) {
				contained = meshv.region(branch).is_in(cv);
				if (!contained) branch++;
				}
				GMM_ASSERT1(contained=true, "No branch region contains node i1!");
				BCv_transp.back().branches.emplace_back(branch); 
				}*/
			}
			else if (meshv.convex_to_point(i1).size()==2){ /* trivial outflow junction */

				// Search for index of first containing branch (\mathcal{P}^{out}_j)
				size_type firstbranch = 0; 
				bool contained = false;
				while (!contained && firstbranch<nb_branches ) {
					contained = meshv.region(firstbranch).is_in(cv);
					if (!contained) firstbranch++;
				}
				GMM_ASSERT1(contained=true, "No branch region contains node i1!");

				// Check if i1 is a trivial junction (or a INT point)
				size_type cv1 = meshv.convex_to_point(i1)[0];
				size_type cv2 = meshv.convex_to_point(i1)[1];
				bool is_junc = (meshv.region(firstbranch).is_in(cv1) < 1 ||
						meshv.region(firstbranch).is_in(cv2) < 1 );

				if (is_junc){
					cout << "Found a trivial junction at i1 = " << i1 << endl;
					// Check if jucntion has been already stored, 
					// if not add to the junction list (J) and build a new region
					dal::bit_vector tmp; tmp.add(i1);
					if(!junctions.contains(tmp)){
						// Store the junction vertex
						junctions.add(i1);
						nb_junctions++;
						GMM_ASSERT1(meshv.has_region(fer)==0, 
								"Overload in meshv region assembling!");
						// Build a new region with idx "first empty region"
						meshv.region(fer).add(cv, 0);
						// Create a new junction node
						Jv.emplace_back("JUN", 0, i1, fer);
						fer++;
					}
					// Search for index of second containing branch (\mathcal{P}^{out}_j)
					size_type secondbranch = firstbranch+1; 
					size_type secondcv = (( cv1 == cv) ? cv2 : cv1);
					contained = false;
					while (!contained && secondbranch<nb_branches ) {
						contained = meshv.region(secondbranch).is_in(secondcv);
						if (!contained) secondbranch++;
					}
					GMM_ASSERT1(contained=true, "No branch region contains node i1!");
					// Add the two branches
					Jv.back().branches.emplace_back(+firstbranch);
					Jv.back().branches.emplace_back(-secondbranch);
					Jv.back().value += param.R(mimv, firstbranch);
					Jv.back().value += param.R(mimv, secondbranch);
				}
			}
			else if (meshv.convex_to_point(i1).size()>=2){ /* non-trivial outflow junction */

				// Search for index of containing branch (\mathcal{P}^{out}_j)
				size_type branch = 0; 
				bool contained = false;
				while (!contained && branch<nb_branches ) {
					contained = meshv.region(branch).is_in(cv);
					if (!contained) branch++;
				}
				GMM_ASSERT1(contained=true, "No branch region contains node i0!");

				// Check if jucntion has been already stored, 
				// if not add to the junction list (J) and build a new region
				dal::bit_vector tmp; tmp.add(i1);
				if(!junctions.contains(tmp)){
					// Store the junction vertex
					junctions.add(i1);
					nb_junctions++;
					GMM_ASSERT1(meshv.has_region(fer)==0, 
							"Overload in meshv region assembling!");
					// Build a new region with idx "first empty region"
					meshv.region(fer).add(cv, 0);
					// Create a new junction node
					Jv.emplace_back("JUN", 0, i1, fer);
					// Add the outflow branch
					Jv.back().branches.emplace_back(+branch);
					Jv.back().value += param.R(mimv, branch);
					//cout << "Branch " << branch << " added to junction " << i1 << endl;
					fer++;
				}
				else {
					// Add the outflow branch (to the right junction node)
					size_type jj = 0;
					bool found = false;
					while (!found && jj < nb_junctions){
						found = (i1 == Jv[jj].idx);
						if (!found) jj++;
					}
					Jv[jj].branches.emplace_back(+branch);
					Jv[jj].value += param.R(mimv, branch);
					//cout << "Branch " << branch << " added to junction " << jj << endl;
				}
			}

		} /* end of convexes loop */

		// Ckeck network assembly
#ifdef M3D1D_VERBOSE_
// 		cout << "--- NETWORK ASSEMBLY ------------------ "   << endl;
// 		cout << "  Branches:   " << nb_branches << endl
// 			<< "  Vertices:   " << nn.size()+1 << endl;
// 		cout << "  Extrema:    " << extrema << endl;	  
// 		for (size_type i=0; i<BCv_transp.size(); ++i)
// 			cout << "    -  label=" << BCv_transp[i].label 
// 				<< ", value=" << BCv_transp[i].value << ", ind=" << BCv_transp[i].idx 
// 				<< ", rg=" << BCv_transp[i].rg << ", branches=" << BCv_transp[i].branches << endl; 
// 		cout << "  Junctions: " << junctions << endl;
// 		for (size_type i=0; i<Jv.size(); ++i)
// 			cout << "    -  label=" << Jv[i].label 
// 				<< ", value=" << Jv[i].value << ", ind=" << Jv[i].idx 
// 				<< ", rg=" << Jv[i].rg << ", branches=" << Jv[i].branches << endl; 
// 		cout << "---------------------------------------- "   << endl;
#endif

	} 
	GMM_STANDARD_CATCH_ERROR; // catches standard errors


} /* end of build_vessel_boundary */



void transport3d1d::assembly (void)
{
	std::cout<<"assemble transport problem"<<std::endl;
	//1 Build the monolithic matrix AM
	assembly_mat(); 
	//2 Build the monolithic rhs FM
	//assembly_rhs();
} // end of assembly

	void 
transport3d1d::assembly_mat(void)
{

#ifdef M3D1D_VERBOSE_
	cout << "Allocating AM, UM, FM ..." << endl;
#endif
	gmm::resize(AM_transp, dof_transp.tot(), dof_transp.tot()); gmm::clear(AM_transp);
	gmm::resize(UM_transp, dof_transp.tot()); gmm::clear(UM_transp);
	gmm::resize(FM_transp, dof_transp.tot()); gmm::clear(FM_transp);
#ifdef M3D1D_VERBOSE_
	cout << "Assembling the monolithic matrix AM ..." << endl;
#endif

	sparse_matrix_type Dt(dof_transp.Ct(), dof_transp.Ct());gmm::clear(Dt);
	sparse_matrix_type Dv(dof_transp.Cv(), dof_transp.Cv());gmm::clear(Dv);	


	// Junction compatibility matrix for the network problem
	//sparse_matrix_type Jvv(dof.Pv(), dof.Uv());
	// Tissue-to-tissue exchange matrix
	sparse_matrix_type Btt(dof_transp.Ct(), dof_transp.Ct());gmm::clear(Btt);
	// Vessel-to-tissue exchange matrix
	sparse_matrix_type Btv(dof_transp.Ct(), dof_transp.Cv());gmm::clear(Btv);
	// Tissue-to-vessel exchange matrix
	sparse_matrix_type Bvt(dof_transp.Cv(), dof_transp.Ct());gmm::clear(Bvt);
	// Vessel-to-vessel exchange matrix
	sparse_matrix_type Bvv(dof_transp.Cv(), dof_transp.Cv());gmm::clear(Bvv);

	// Aux tissue-to-vessel averaging matrix
	sparse_matrix_type Mbar(dof_transp.Cv(), dof_transp.Ct());gmm::clear(Mbar);
	// Aux tissue-to-vessel interpolation matrix
	sparse_matrix_type Mlin(dof_transp.Cv(), dof_transp.Ct());gmm::clear(Mlin);


#ifdef M3D1D_VERBOSE_
	cout << "  Assembling  Dt ..." << endl;
#endif

	bool DIFFUSION_T = PARAM.int_value("DIFFUSION_T", "flag for diffusion term in tissue");
	if(DIFFUSION_T ==1){
		getfem::asm_stiffness_matrix_for_laplacian(Dt,mimt,mf_Ct, mf_coeft, param_transp.At()); //Divergence matrix	
		gmm::add(Dt,
				gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_transp.Ct()), 
					gmm::sub_interval(0, dof_transp.Ct()))); 

	} 
#ifdef M3D1D_VERBOSE_
	cout << "  Assembling  Dv ..." << endl;
#endif
	bool FIXED_SOURCE = PARAM.int_value("FIXED_SOURCE", "flag for fixed source ");
	if(FIXED_SOURCE==0){

		bool DIFFUSION_V = PARAM.int_value("DIFFUSION_V", "flag for diffusion term in vessel");
		if(DIFFUSION_V ==1){
			getfem::asm_stiffness_matrix_for_laplacian(Dv,mimv,mf_Cv, mf_coefv, param_transp.Av()); //Divergence matrix	
			gmm::add(Dv,
					gmm::sub_matrix(AM_transp, 
						gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()), 
						gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()))); 

		} 
	}



#ifdef M3D1D_VERBOSE_
	cout << "  Assembling aux exchange matrices Mbar and Mlin ..." << endl;
#endif
	asm_exchange_aux_mat(Mbar, Mlin, 
			mimv, mf_Ct, mf_Cv, param.R(), descr.NInt);
#ifdef M3D1D_VERBOSE_
	cout << "  Assembling exchange matrices ..." << endl;
#endif
	bool NEWFORM = PARAM.int_value("NEW_FORMULATION", "flag for the new formulation");
	vector_type coeff(dof.Pv());

	//ricalcolo Rv_coeff=(1-sigma)*Q*((p_v-\bar(p_t) -sigma*(pi_v-pi_t)))*0.5


	sparse_matrix_type Mbar_(dof.Pv(), dof.Pt());
	sparse_matrix_type Mlin_(dof.Pv(), dof.Pt());

	asm_exchange_aux_mat(Mbar_, Mlin_, mimv, mf_Pt, mf_Pv, param.R(), descr.NInt);

	asm_exchange_mat_transp(Btt, Btv, Bvt, Bvv,
			mimv, mf_Cv, mf_coefv, Mbar, Mlin, param_transp.Y(), NEWFORM);
	

	

	bool COUPLING = PARAM.int_value("COUPLING", "flag for coupling-exchange term ");
	if(COUPLING==1){
		
		gmm::add(gmm::scaled(Btt, 2.0*pi*param.R(0)),			 

			gmm::sub_matrix(AM_transp, 
				gmm::sub_interval(0, dof_transp.Ct()), 
				gmm::sub_interval(0, dof_transp.Ct()))); 

		gmm::MatrixMarket_IO::write("Btt_fede.mm",gmm::scaled(Btt, 2.0*pi*param.R(0)));
		
		gmm::add(gmm::scaled(Btv, -2.0*pi*param.R(0)),									

				gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_transp.Ct()),
					gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()))); 

		gmm::add(gmm::scaled(Bvt, -2.0/param.R(0)),  	
				gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()),
					gmm::sub_interval(0, dof_transp.Ct())));


		gmm::add(gmm::scaled(Bvv, 2.0/param.R(0)),								
				gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()), 
					gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()))); 

	}

	if(COUPLING==0 && FIXED_SOURCE==1){

		gmm::add(gmm::scaled(Btt, 2.0*pi*param.R(0)),			 

			gmm::sub_matrix(AM_transp, 
				gmm::sub_interval(0, dof_transp.Ct()), 
				gmm::sub_interval(0, dof_transp.Ct()))); 

		gmm::MatrixMarket_IO::write("Btt_fede.mm",gmm::scaled(Btt, 2.0*pi*param.R(0)));

		gmm::scale(Btv, 2.0*pi*param.R(0));


		//SE VOGLIO SOURCE COEF COSTANTE DA INPUT	

		//gmm::mult(Btv,param_transp.Source_coef(),TTmp);

		//gmm::add(TTmp,	gmm::sub_vector(FM_transp, 
		//			gmm::sub_interval(0,dof_transp.Ct()))); 
		//gmm::clear(TTmp);
		//gmm::scale(Bvt1, -2.0/param.R(0));
		//gmm::mult(Btv1,param_transp.Source_coef(),TTmp);

		//gmm::add(TTmp,	gmm::sub_vector(FM_transp, 
		//			gmm::sub_interval(0,dof_transp.Ct()))); 
		//gmm::clear(TTmp);


		//costruisco il vettore dei coefficienti della sorgente	perchè non lo voglio costante da inpur
		vector_type x(dof_transp.Cv());
		vector_type Source_coef(dof_transp.Cv());
		for(scalar_type k=0;k<dof_transp.Cv();k++){
			x[k]=(k+1)*(1.0/(dof_transp.Cv()));
			Source_coef[k]=1.0-x[k];

		}
		gmm::clear(x);	

		vector_type TTmp(dof_transp.Ct());
		gmm::mult(Btv,Source_coef,TTmp);

		gmm::add(TTmp,	gmm::sub_vector(FM_transp, 
					gmm::sub_interval(0,dof_transp.Ct()))); 
		gmm::clear(TTmp);

	}

	// De-allocate memory

	gmm::clear(Dt); gmm::clear(Dv); 
	gmm::clear(Mbar);  gmm::clear(Mlin);
	gmm::clear(Btt);  gmm::clear(Btv);
	gmm::clear(Bvt);  gmm::clear(Bvv);
}

	void 
transport3d1d::assembly_rhs(void)
{


#ifdef M3D1D_VERBOSE_
	cout << "Assembling the monolithic rhs FM ... " << endl;
#endif
#ifdef M3D1D_VERBOSE_
	cout << "  Initializing RHS for FM ..." << endl;
#endif

bool NONNULL_RHS = PARAM.int_value("NONNULL_RHS", "flag for non null rhs ");
	double source=PARAM.real_value("Source_coef","coefficient for fixed source");
	if(NONNULL_RHS==1){
		/*//non null RHS for tissue problem
		vector_type RHS_coef(dof_transp.Ct());
		for (int i=0;i<dof_transp.Ct();i++){
			RHS_coef[i]=1.0;	
		}
		vector_type RHS(dof_transp.Ct());
		getfem::asm_source_term(RHS,mimt,mf_Ct,mf_Ct,RHS_coef);
		std::cout<<"ho aggiumnto rhs=1"<<std::endl;
		gmm::add(RHS, 
				gmm::sub_vector(FM_temp,
					gmm::sub_interval(0, dof_transp.Ct())));	*/

		//non null RHS for vessel problem
		vector_type RHS_coef(dof_transp.tot());
		for (int i=0;i<dof_transp.Ct();i++){
			RHS_coef[i]=0.0;	
		}
		for (int i=dof_transp.Ct();i<dof_transp.tot();i++){
			RHS_coef[i]=source;	
		}
		vector_type RHS(dof_transp.tot());
		getfem::asm_source_term(gmm::sub_vector(RHS,gmm::sub_interval(0, dof_transp.Ct())),mimt,mf_Ct,mf_Ct,gmm::sub_vector(RHS_coef,gmm::sub_interval(0, dof_transp.Ct())));
		getfem::asm_source_term(gmm::sub_vector(RHS,gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())),mimv,mf_Cv,mf_Cv,gmm::sub_vector(RHS_coef,gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())));
		std::cout<<"ho aggiumnto rhs=1"<<std::endl;
		gmm::add(RHS, FM_temp);	
	}



#ifdef M3D1D_VERBOSE_
	cout << "  Building tissue boundary term ..." << endl;
#endif

	cout<<"assemble BC for tissue"<<endl;

	sparse_matrix_type Att(dof_transp.Ct(),   dof_transp.Ct());
	vector_type Ft(dof_transp.Ct());
	gmm::add(	gmm::sub_matrix(AM_temp,
				gmm::sub_interval(0,dof_transp.Ct()),
				gmm::sub_interval(0,dof_transp.Ct()))
			, Att);
	gmm::scale(	gmm::sub_matrix(AM_temp,
				gmm::sub_interval(0,dof_transp.Ct()),
				gmm::sub_interval(0,dof_transp.Ct()))
			, 0.0);	

	gmm::add(	 gmm::sub_vector(FM_temp, gmm::sub_interval(0,dof_transp.Ct()))
			,Ft);	 
	gmm::scale(	 gmm::sub_vector(FM_temp, gmm::sub_interval(0,dof_transp.Ct()))
			,0.0);


	scalar_type beta_t  = PARAM.real_value("BETAtissue_transp", "Coefficient for mixed BC for transport problem in tissue");
	asm_tissue_bc_transp(Att, Ft, mimt, mf_Ct, mf_coeft, BCt_transp,beta_t);
	gmm::add(Att, 
			gmm::sub_matrix(AM_temp,
				gmm::sub_interval(0,dof_transp.Ct()),
				gmm::sub_interval(0,dof_transp.Ct())));
	gmm::add(Ft, 
			gmm::sub_vector(FM_temp,
				gmm::sub_interval(0,dof_transp.Ct())));
	// De-allocate memory
	gmm::clear(Att);
	gmm::clear(Ft);

	bool FIXED_SOURCE = PARAM.int_value("FIXED_SOURCE", "flag for fixed source ");
	if(FIXED_SOURCE==0){
		//costruisco le condizioni al contorno: influenzeranno AM_temp e Fv
		sparse_matrix_type Avv(dof_transp.Cv(), dof_transp.Cv());
		vector_type Fv(dof_transp.Cv());
		gmm::add(	gmm::sub_matrix(AM_temp,
					gmm::sub_interval(dof_transp.Ct(),dof_transp.Cv()),
					gmm::sub_interval(dof_transp.Ct(),dof_transp.Cv()))
				, Avv);
		gmm::scale(	gmm::sub_matrix(AM_temp,
					gmm::sub_interval(dof_transp.Ct(),dof_transp.Cv()),
					gmm::sub_interval(dof_transp.Ct(),dof_transp.Cv()))
				, 0.0);	

		gmm::add(	 gmm::sub_vector(FM_temp, gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()))
				,Fv);	
		gmm::scale(	 gmm::sub_vector(FM_temp, gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()))
				,0.0);

		scalar_type beta_v  = PARAM.real_value("BETAvessel_transp", "Coefficient for mixed BC for transport problem in vessels");
		asm_network_bc_transp(Avv, Fv, mimv, mf_Cv, mf_coefv, BCv_transp, beta_v );
		gmm::add(Avv, 
				gmm::sub_matrix(AM_temp,
					gmm::sub_interval(dof_transp.Ct(),dof_transp.Cv()),
					gmm::sub_interval(dof_transp.Ct(),dof_transp.Cv())));
		gmm::add(Fv, 
				gmm::sub_vector(FM_temp,
					gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())));
		// De-allocate memory
		gmm::clear(Avv);
		gmm::clear(Fv);
	}
	
	

	
}

void transport3d1d::update (vector_type Pigreco){

	// ho la matrice AM assemblata con tutti i termini (da non modificare mai più!)
	// ho  il vettore F ancora vuoto (da aggiungere: il termine temporale; le condizioni al contorno 1d)
	gmm::copy(AM_transp, AM_temp);
	gmm::copy(FM_transp, FM_temp);

	bool SATURATION  = PARAM.int_value("SATURATION","Flag to use the saturation per the adhesive term");

	if(SATURATION==1){
		sparse_matrix_type Adhv(dof_transp.Cv(), dof_transp.Cv()); gmm::clear(Adhv);
		asm_network_nano_transp(Adhv, mimv, mf_Cv, Pigreco); 	
		masslumping(Adhv);
		gmm::add(Adhv, gmm::sub_matrix(AM_temp, 
					gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()),
					gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())));
	}



	// faccio l'update temporale: sommo a F le concentrazioni al passo precedente

	// update rhs (there is the time step mass term)

	vector_type TFt(dof_transp.Ct());
	vector_type TFv(dof_transp.Cv());
	//vector_type Ct(dof_transp.Ct());
	//vector_type Cv(dof_transp.Cv());
	//gmm::add(gmm::sub_vector(UM_transp, gmm::sub_interval(0, dof_transp.Ct())), Ct);
	//gmm::add(gmm::sub_vector(UM_transp, gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())), Cv);
	asm_source_term(TFt,mimt, mf_Ct, mf_Ct,gmm::sub_vector(UM_transp, gmm::sub_interval(0, dof_transp.Ct()))); // tentativo con generic_assembly: lo trovi in 										assembling3d_transp, la funzione è asm_time_rhs_transp
	asm_source_term(TFv,mimv, mf_Cv, mf_Cv,gmm::sub_vector(UM_transp, gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())));
	gmm::scale(TFt, (1.0/param_transp.dt())); // dt time step
	gmm::scale(TFv, (1.0/param_transp.dt())); // dt time step
	gmm::add(TFt, gmm::sub_vector(FM_temp, gmm::sub_interval(0, dof_transp.Ct())));
	gmm::add(TFv, gmm::sub_vector(FM_temp, gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())));
	gmm::clear(UM_transp);
	gmm::clear(TFt); gmm::clear(TFv);

	assembly_rhs();


}

bool transport3d1d::solve (void)
{
	std::cout<<"solve transport problem"<<std::endl<<std::endl;
#ifdef M3D1D_VERBOSE_
	cout << "Solving the monolithic system ... " << endl;
#endif
	gmm::resize(AM_temp, dof_transp.tot(), dof_transp.tot()); gmm::clear(AM_temp);
	gmm::resize(FM_temp, dof_transp.tot()); gmm::clear(FM_temp);

	double time = gmm::uclock_sec();
	double time_count = 0;	
	int iteraz = 0;

	std::ofstream outMeanCt("MeanCt.txt"); //contiene la media di Ct sul dominio ad ogni istante

	for(double t=0;t<=param_transp.T()*(!descr_transp.STATIONARY) ; t = t + param_transp.dt() + (param_transp.dt()==0) ){ 
		time_count++; 
		iteraz++; 
		std::cout<<"iteration number:"<<time_count<<std::endl;
		std::cout<<"time = "<<t<<" s"<<std::endl;	

		//FEDE aggungo qui quello che mi serve da update
		gmm::copy(AM_transp, AM_temp);
		gmm::copy(FM_transp, FM_temp);
		assembly_rhs();
               

		gmm::csc_matrix<scalar_type> A_transp;
		gmm::clean(AM_transp, 1E-12);
		gmm::copy(AM_temp, A_transp);
		//FEDE esporto la matrice in formato MatrixMarket
		gmm::MatrixMarket_IO::write("A_fede.mm",A_transp);

		vector_type F_transp(gmm::vect_size(FM_transp));
		gmm::clean(FM_transp, 1E-12);
		gmm::copy(FM_temp, F_transp);

		std::ofstream outF("F.txt");
		outF << gmm::col_vector(FM_temp);
		outF.close();


		//gmm::clear(AM_transp); // to be postponed for preconditioner





		if ( descr_transp.SOLVE_METHOD == "SuperLU" ) { // direct solver //
			bool FIXED_SOURCE = PARAM.int_value("FIXED_SOURCE", "flag for fixed source ");

			//FEDE: calcolo il tempo effettivo di risoluzione del sistema con Super LU
			double time_LU = gmm::uclock_sec();

			if(FIXED_SOURCE==0){
#ifdef M3D1D_VERBOSE_
				cout << "  Applying the SuperLU method ... " << endl;
#endif	
				scalar_type cond;
				gmm::SuperLU_solve(A_transp, UM_transp, F_transp, cond);
				cout << "  Condition number (transport problem): " << cond << endl;
			}
			else{
#ifdef M3D1D_VERBOSE_
				cout << "  Applying the SuperLU method ... " << endl;
#endif
				scalar_type cond;
				gmm::SuperLU_solve(gmm::sub_matrix(A_transp, 
							gmm::sub_interval(0, dof_transp.Ct()),
							gmm::sub_interval(0, dof_transp.Ct())),
						gmm::sub_vector(UM_transp, gmm::sub_interval(0, dof_transp.Ct())),
						gmm::sub_vector(F_transp, gmm::sub_interval(0, dof_transp.Ct())),
						cond);

				//costruisco il vettore dei coefficienti della sorgente	perchè non lo voglio costante da input
				vector_type x(dof_transp.Cv());
				vector_type Source_coef(dof_transp.Cv());
				for(scalar_type k=0;k<dof_transp.Cv();k++){
					x[k]=(k+1)*(1.0/(dof_transp.Cv()));
					Source_coef[k]=1.0-x[k];
				}
				gmm::clear(x);	

				for(scalar_type k=0;k<dof_transp.Cv();k++)
					UM_transp[dof_transp.Ct()+k]=Source_coef[k];

				cout << "  Condition number (transport problem): " << cond << endl;
			}	
			//FEDE
			cout << endl<<"... time to solve the system using Super LU method: " << gmm::uclock_sec() - time_LU << " seconds\n";	
		}
		else { // Iterative solver //

			// Iterations
			gmm::iteration iter(descr_transp.RES);  // iteration object with the max residu
			iter.set_noisy(1);               // output of iterations (2: sub-iteration)
			iter.set_maxiter(descr_transp.MAXITER); // maximum number of iterations

			// Preconditioners
			//! \todo Add preconditioner choice to param file
			// See \link http://download.gna.org/getfem/html/homepage/gmm/iter.html
			gmm::identity_matrix PM; // no precond
			//gmm::diagonal_precond<sparse_matrix_type> PM(AM); // diagonal preocond
			//gmm::ilu_precond<sparse_matrix_type> PM(AM);
			// ...
			//gmm::clear(AM);
			// See <http://download.gna.org/getfem/doc/gmmuser.pdf>, pag 15

			double time_iter = gmm::uclock_sec();
			if ( descr_transp.SOLVE_METHOD == "CG" ) {
#ifdef M3D1D_VERBOSE_
				cout << "  Applying the Conjugate Gradient method ... " << endl;
#endif
				gmm::identity_matrix PS;  // optional scalar product
				gmm::cg(AM_transp, UM_transp, F_transp, PS, PM, iter);
			}
			else if ( descr_transp.SOLVE_METHOD == "BiCGstab" ) {
#ifdef M3D1D_VERBOSE_
				cout << "  Applying the BiConjugate Gradient Stabilized method ... " << endl;
#endif
				gmm::bicgstab(AM, UM, FM, PM, iter);
			}
			else if ( descr_transp.SOLVE_METHOD == "GMRES" ) {
#ifdef M3D1D_VERBOSE_
				cout << "  Applying the Generalized Minimum Residual method ... " << endl;
#endif
				size_type restart = 50;
				gmm::gmres(A_transp, UM, FM, PM, restart, iter);
			}
			else if ( descr_transp.SOLVE_METHOD == "QMR" ) {
#ifdef M3D1D_VERBOSE_
				cout << "  Applying the Quasi-Minimal Residual method ... " << endl;
#endif
				gmm::qmr(AM, UM, FM, PM, iter);
			}
			else if ( descr_transp.SOLVE_METHOD == "LSCG" ) {
#ifdef M3D1D_VERBOSE_
				cout << "  Applying the unpreconditionned Least Square CG method ... " << endl;
#endif
				gmm::least_squares_cg(AM, UM, FM, iter);
			}
			// Check
			if (iter.converged())
				cout << "  ... converged in " << iter.get_iteration() << " iterations." << endl;
			else if (iter.get_iteration() == descr_transp.MAXITER)
				cerr << "  ... reached the maximum number of iterations!" << endl;

			std::cout << endl<<"... time to solve the system iterative  method: " << gmm::uclock_sec() - time_iter << " seconds\n";	
		}


		//export solution
		std::cout<<"solved! going to export..."<<std::endl;
		string time_suff = "";
		std::ostringstream convert;
		convert << time_count;
		time_suff = convert.str();
		//Compute the total c_t in the tissue



export_vtk("1"); 
export_vtk("2"); 

		outMeanCt<<mean_ct()<<std::endl;


		if(t==0){
			export_vtk(time_suff); 
			std::cout<<"exported! now new iteration..."<<std::endl;
		}
		//if(iteraz%30==0){ 
		if(iteraz%10==0){ 
			export_vtk(time_suff); 
		}

	} //end of cycle over time 

	outMeanCt.close();

	cout << endl<<"... total time to solve : " << gmm::uclock_sec() - time << " seconds\n";

	return true;
	}; // end of solve








	bool transport3d1d::solve_samg (void)
	{
#ifdef USE_SAMG
		std::cout<<"solve transport problem"<<std::endl<<std::endl;
#ifdef M3D1D_VERBOSE_
		cout << "Solving the monolithic system ... " << endl;
#endif
		gmm::resize(AM_temp, dof_transp.tot(), dof_transp.tot()); gmm::clear(AM_temp);
		gmm::resize(FM_temp, dof_transp.tot()); gmm::clear(FM_temp);

		double time = gmm::uclock_sec();
		double time_count = 0;	
		int iteraz = 0;


		for(double t=0;t<=param_transp.T()*(!descr_transp.STATIONARY) ; t = t + param_transp.dt() + (param_transp.dt()==0) ){ 
			time_count++; 
			iteraz++; 
			std::cout<<"iteration number:"<<time_count<<std::endl;
			std::cout<<"time = "<<t<<" s"<<std::endl;	

			//aggungo qui quello che mi serve da update
			gmm::copy(AM_transp, AM_temp);
			gmm::copy(FM_transp, FM_temp);
			assembly_rhs();


			gmm::csc_matrix<scalar_type> A_transp;
			gmm::clean(AM_transp, 1E-12);
			gmm::copy(AM_temp, A_transp);

			vector_type F_transp(gmm::vect_size(FM_transp));
			gmm::clean(FM_transp, 1E-12);
			gmm::copy(FM_temp, F_transp);

			//gmm::clear(AM_transp); // to be postponed for preconditioner

			float told_ass;
			SAMG_CTIME(&told_ass);
			double time_ass=gmm::uclock_sec();	

			//Interface with csc matrix
			//
#ifdef CSC_INTERFACE

			bool COUPLING = PARAM.int_value("COUPLING", "flag for coupling-exchange term ");
			bool FIXED_SOURCE = PARAM.int_value("FIXED_SOURCE", "flag for fixed source ");

			gmm::csc_matrix<scalar_type> A_csc;gmm::clean(AM_temp, 1E-12);// the clean remove the very low values
			std::vector<scalar_type> U_2, B_2; /* main unknown, and right hand side  */

			if(COUPLING==0 && FIXED_SOURCE==1 ){			

				// gmm::clean(A_2, 1E-12);
				// gmm::copy(AM_temp, A_2);
				gmm::copy(gmm::transposed(gmm::sub_matrix(AM_temp, 
							gmm::sub_interval(0, dof_transp.Ct()),
							gmm::sub_interval(0, dof_transp.Ct()))), A_csc);
				//
				//

				gmm::resize(U_2,dof_transp.Ct());gmm::clean(U_2, 1E-12);
				gmm::copy(gmm::sub_vector(UM_transp, gmm::sub_interval(0, dof_transp.Ct())),U_2);

				gmm::resize(B_2,dof_transp.Ct());gmm::clean(B_2, 1E-12);
				gmm::copy(gmm::sub_vector(F_transp, gmm::sub_interval(0, dof_transp.Ct())),B_2);

			}
			
		bool NONNULL_RHS = PARAM.int_value("NONNULL_RHS", "flag for non null rhs ");
			if(COUPLING==0 && NONNULL_RHS==1){			

				// gmm::clean(A_2, 1E-12);
				// gmm::copy(AM_temp, A_2);
				gmm::copy(gmm::transposed(gmm::sub_matrix(AM_temp, 
							gmm::sub_interval(0, dof_transp.Ct()),
							gmm::sub_interval(0, dof_transp.Ct()))), A_csc);
				//
				//

				gmm::resize(U_2,dof_transp.Ct());gmm::clean(U_2, 1E-12);
				gmm::copy(gmm::sub_vector(UM_transp, gmm::sub_interval(0, dof_transp.Ct())),U_2);

				gmm::resize(B_2,dof_transp.Ct());gmm::clean(B_2, 1E-12);
				gmm::copy(gmm::sub_vector(F_transp, gmm::sub_interval(0, dof_transp.Ct())),B_2);

			}

			if(COUPLING==1){			

				// gmm::clean(A_2, 1E-12);
				// gmm::copy(AM_temp, A_2);

				//gmm::clean(AM_temp, 1E-12)
				gmm::copy(gmm::transposed(AM_temp), A_csc);

				gmm::resize(U_2,dof_transp.tot()); gmm::clean(U_2, 1E-12);
				gmm::copy(UM_transp,U_2);

				gmm::resize(B_2,dof_transp.tot());gmm::clean(B_2, 1E-12);
				gmm::copy(F_transp,B_2);


			}

   			std::cout<<"*** parameters SAMG matrix   "<<A_csc.nrows()<<" parameters gmm matrix "<<A_transp.nrows()<<std::endl;	
			std::cout<<"*** parameters SAMG matrix   "<<A_csc.ncols()<<" parameters gmm matrix "<<A_transp.ncols()<<std::endl;	

			//Print values of csc matrix    
			  std::ofstream output_file("./ir.txt");
			  std::ofstream output_filev("./v.txt");
			  std::ofstream output_filec("./jc.txt");
			  std::ostream_iterator<int> output_iterator(output_file, "\n");
			  std::copy(A_csc.ir.begin(), A_csc.ir.end(), output_iterator);    
			  std::ostream_iterator<int> output_iteratorc(output_filec, "\n");
			  std::copy(A_csc.jc.begin(), A_csc.jc.end(), output_iteratorc);    
			  std::ostream_iterator<double> output_iteratorv(output_filev, "\n");
			  std::copy(A_csc.pr.begin(), A_csc.pr.end(), output_iteratorv);    

			//creating samg a matrix

			APPL_INT nnu,nna;
			nnu=A_csc.nc; nna=A_csc.pr.size();
			double *a_samg;a_samg=new double[nna];
			double *u_samg;u_samg=new double[nnu];
			double *b_samg;b_samg=new double[nnu];
			APPL_INT *ja_samg;ja_samg=new APPL_INT[nna];
			APPL_INT *ia_samg;ia_samg=new APPL_INT[nnu+1];
			ia_samg[0]=1;
			unsigned int offset=0;
			for(int ia=0;ia<A_csc.jc.size()-1;ia++){
				unsigned int nonzero=A_csc.jc.at(ia+1)-A_csc.jc.at(ia);
				ia_samg[ia+1] = nonzero + ia_samg[ia] ;
				u_samg[ia] = U_2[ia];
				b_samg[ia] = B_2[ia];
				// std::cout<<"Going from "<<offset<<" to "<<offset + nonzero<<std::endl;
				for (int innz=offset; innz<offset+nonzero; innz++ )
					if( ia == (int) A_csc.ir.at(innz))
					{
						//	std::cout<<"diagonal term"<<std::endl; 
						a_samg[offset]=A_csc.pr.at(innz);
						ja_samg[offset]=A_csc.ir.at(innz)+1;

					}
				int shift=1;
				for (int innz=offset; innz<offset+nonzero; innz++ )
					if( ia != (int) A_csc.ir.at(innz))
					{
						//	std::cout<<"non diagonal term "<< ia << " " << A_csc.ir.at(innz) <<std::endl; 
						a_samg[offset+shift]=A_csc.pr.at(innz);
						ja_samg[offset+shift]=A_csc.ir.at(innz)+1;
						shift++;
					}
				//	std::cout<<std::endl;
				// std::cout<<"non zero in i "<<ia<<" are "<< nonzero<<std::endl; 
				offset+=nonzero;
			}
			
#endif //cscinterface

#ifdef CSR_INTERFACE

			bool COUPLING = PARAM.int_value("COUPLING", "flag for coupling-exchange term ");
			bool FIXED_SOURCE = PARAM.int_value("FIXED_SOURCE", "flag for fixed source ");

			gmm::csr_matrix<scalar_type> A_csr; gmm::clean(AM_temp, 1E-12);// the clean remove the very low values
			std::vector<scalar_type> U_2, B_2; /* main unknown, and right hand side  */

			if(COUPLING==0 && FIXED_SOURCE==1 ){			

				// gmm::clean(A_2, 1E-12);
				// gmm::copy(AM_temp, A_2);
				gmm::copy(gmm::sub_matrix(AM_temp, 
							gmm::sub_interval(0, dof_transp.Ct()),
							gmm::sub_interval(0, dof_transp.Ct())), A_csr);
				//
				//

				gmm::resize(U_2,dof_transp.Ct());gmm::clean(U_2, 1E-12);
				gmm::copy(gmm::sub_vector(UM_transp, gmm::sub_interval(0, dof_transp.Ct())),U_2);

				gmm::resize(B_2,dof_transp.Ct());gmm::clean(B_2, 1E-12);
				gmm::copy(gmm::sub_vector(F_transp, gmm::sub_interval(0, dof_transp.Ct())),B_2);

			}
			
			bool NONNULL_RHS = PARAM.int_value("NONNULL_RHS", "flag for non null rhs ");

			if(COUPLING==0 && NONNULL_RHS==1){			

				// gmm::clean(A_2, 1E-12);
				// gmm::copy(AM_temp, A_2);
				gmm::copy(gmm::sub_matrix(AM_temp, 
							gmm::sub_interval(0, dof_transp.Ct()),
							gmm::sub_interval(0, dof_transp.Ct())), A_csr);
				//
				//

				gmm::resize(U_2,dof_transp.Ct());gmm::clean(U_2, 1E-12);
				gmm::copy(gmm::sub_vector(UM_transp, gmm::sub_interval(0, dof_transp.Ct())),U_2);

				gmm::resize(B_2,dof_transp.Ct());gmm::clean(B_2, 1E-12);
				gmm::copy(gmm::sub_vector(F_transp, gmm::sub_interval(0, dof_transp.Ct())),B_2);

			}

			if(COUPLING==1){			

				// gmm::clean(A_2, 1E-12);
				// gmm::copy(AM_temp, A_2);

				
				// gmm::copy(gmm::transposed(gmm::transposed(AM_temp)), A_csr);
gmm::copy(AM_temp, A_csr);
				gmm::resize(U_2,dof_transp.tot()); gmm::clean(U_2, 1E-12);
				gmm::copy(UM_transp,U_2);

				gmm::resize(B_2,dof_transp.tot());gmm::clean(B_2, 1E-12);
				gmm::copy(F_transp,B_2);


			}

   			std::cout<<"*** parameters SAMG matrix   "<<A_csr.nrows()<<" parameters gmm matrix "<<A_transp.nrows()<<std::endl;	
			std::cout<<"*** parameters SAMG matrix   "<<A_csr.ncols()<<" parameters gmm matrix "<<A_transp.ncols()<<std::endl;	

			//Print values of csr matrix    
			   std::ofstream output_file("./ir_csr.txt");
			     std::ofstream output_filev("./v_csr.txt");
			     std::ofstream output_filec("./jc_csr.txt");
			     std::ostream_iterator<int> output_iterator(output_file, "\n");
			     std::copy(A_csr.ir.begin(), A_csr.ir.end(), output_iterator);    
			     std::ostream_iterator<int> output_iteratorc(output_filec, "\n");
			     std::copy(A_csr.jc.begin(), A_csr.jc.end(), output_iteratorc);    
			     std::ostream_iterator<double> output_iteratorv(output_filev, "\n");
			     std::copy(A_csr.pr.begin(), A_csr.pr.end(), output_iteratorv);    

			//creating samg a matrix

			APPL_INT nnu,nna;
			nnu=A_csr.nc; nna=A_csr.pr.size();
			double *a_samg;a_samg=new double[nna];
			double *u_samg;u_samg=new double[nnu];
			double *b_samg;b_samg=new double[nnu];
			APPL_INT *ja_samg;ja_samg=new APPL_INT[nna];
			APPL_INT *ia_samg;ia_samg=new APPL_INT[nnu+1];
			ia_samg[0]=1;
			unsigned int offset=0;
			for(int ia=0;ia<nnu+1-1;ia++){
				unsigned int nonzero=A_csr.jc.at(ia+1)-A_csr.jc.at(ia);
				ia_samg[ia+1] = nonzero + ia_samg[ia] ;
				u_samg[ia] = U_2[ia];
				b_samg[ia] = B_2[ia];
				// std::cout<<"Going from "<<offset<<" to "<<offset + nonzero<<std::endl;
				for (int innz=offset; innz<offset+nonzero; innz++ )
					if( ia == (int) A_csr.ir.at(innz))
					{	//WARNING: not working for a matrix with zero on the diagonal
						//	std::cout<<"diagonal term"<<std::endl; 
						if(fabs(A_csr.pr.at(innz))<1E-30)
							std::cout<<"****************************diagonal zero"<<std::endl;
						a_samg[offset]=A_csr.pr.at(innz);
						ja_samg[offset]=A_csr.ir.at(innz)+1;

					}
				int shift=1;
				for (int innz=offset; innz<offset+nonzero; innz++ )
					if( ia != (int) A_csr.ir.at(innz))
					{
						//	std::cout<<"non diagonal term "<< ia << " " << A_csr.ir.at(innz) <<std::endl; 
												if(fabs(A_csr.pr.at(innz))<1E-30)
							std::cout<<"****************************non diagonal zero "<<A_csr.pr.at(innz)<<std::endl;
						a_samg[offset+shift]=A_csr.pr.at(innz);
						ja_samg[offset+shift]=A_csr.ir.at(innz)+1;
						shift++;
					}
				//	std::cout<<std::endl;
				// std::cout<<"non zero in i "<<ia<<" are "<< nonzero<<std::endl; 
				offset+=nonzero;
			}
			
#endif //csrinterface


			//std::copy(A_csc.pr.begin(), A_csc.pr.end(), a_samg);	
#ifdef SPARSE_INTERFACE	
			bool COUPLING = PARAM.int_value("COUPLING", "flag for coupling-exchange term ");
			bool FIXED_SOURCE = PARAM.int_value("FIXED_SOURCE", "flag for fixed source ");
			sparse_matrix_type AM_1;
			std::vector<scalar_type> U_1, B_1;      /* main unknown, and right hand side  */
			if(COUPLING==0 && FIXED_SOURCE==1){
				//==========================interfce against sparse matrix

				gmm::resize(AM_1,dof_transp.Ct(), dof_transp.Ct());gmm::clean(AM_1, 1E-12);
				std::cout<<"--------------- number of dof ---------" << dof_transp.Ct()<<std::endl;
				gmm::copy(gmm::sub_matrix(A_transp, 
							gmm::sub_interval(0, dof_transp.Ct()),
							gmm::sub_interval(0, dof_transp.Ct())),AM_1);


				gmm::resize(U_1,dof_transp.Ct());gmm::clean(U_1, 1E-12);
				gmm::copy(gmm::sub_vector(UM_transp, gmm::sub_interval(0, dof_transp.Ct())),U_1);

				gmm::resize(B_1,dof_transp.Ct());gmm::clean(B_1, 1E-12);
				gmm::copy(gmm::sub_vector(F_transp, gmm::sub_interval(0, dof_transp.Ct())),B_1);
			}			
			
			bool NONNULL_RHS = PARAM.int_value("NONNULL_RHS", "flag for non null rhs ");
			if(COUPLING==0 && NONNULL_RHS==1){			

				gmm::resize(AM_1,dof_transp.Ct(), dof_transp.Ct());gmm::clean(AM_1, 1E-12);
				std::cout<<"--------------- number of dof ---------" << dof_transp.Ct()<<std::endl;
				gmm::copy(gmm::sub_matrix(A_transp, 
							gmm::sub_interval(0, dof_transp.Ct()),
							gmm::sub_interval(0, dof_transp.Ct())),AM_1);


				gmm::resize(U_1,dof_transp.Ct());gmm::clean(U_1, 1E-12);
				gmm::copy(gmm::sub_vector(UM_transp, gmm::sub_interval(0, dof_transp.Ct())),U_1);

				gmm::resize(B_1,dof_transp.Ct());gmm::clean(B_1, 1E-12);
				gmm::copy(gmm::sub_vector(F_transp, gmm::sub_interval(0, dof_transp.Ct())),B_1);
			}

			if(COUPLING==1){			

				//==========================interfce against sparse matrix

				gmm::resize(AM_1,dof_transp.tot(), dof_transp.tot());gmm::clean(AM_1, 1E-12);
				std::cout<<"--------------- number of dof ---------" << dof_transp.tot()<<std::endl;
				gmm::copy(A_transp,AM_1);


				gmm::resize(U_1,dof_transp.tot());gmm::clean(U_1, 1E-12);
				gmm::copy(UM_transp,U_1);

				gmm::resize(B_1,dof_transp.tot());gmm::clean(B_1, 1E-12);
				gmm::copy(F_transp,B_1);
			}



			////SAMG 
			//
			//
			std::vector<int> ja_v ;
			std::vector<double> a_v ; 
			int non_zero=0;
			int total_non_zero=0;
			//
			double * a, * u, * f;
			int nrows=AM_1.nrows();
			f= new double[nrows];	u= new double[nrows];
			APPL_INT * ia, * ja;
			ia = new APPL_INT[nrows+1];
			ia[0]=(APPL_INT) 1;	
			for(int i = 0 ; i < nrows ; i++ ){
				non_zero=0;
				for(int j = 0 ; j < nrows ; j++ ){
					if (AM_1(i,j)*AM_1(i,j) > 1.e-40 || i==j ) {
						a_v.push_back(AM_1(i,j));
						//		if (i==j) {
						//		 a_v.push_back(10);
						ja_v.push_back(j);	
						non_zero++;
					}
					}
					f[i]=B_1[i];
					// f[i]=U[i];
					// u[i]=U[i];
					u[i]=0;	ia[i+1]=ia[i]+non_zero;
					total_non_zero+=non_zero;
				}
				std::cout<< "We have "<< total_non_zero << " non zero values "<<std::endl;
				ja=new APPL_INT[total_non_zero];
				a=new double [total_non_zero];
				APPL_INT nnu,nna;
				nnu=nrows; nna=total_non_zero;

				for(int i = 0 ; i < nnu ; i++ ){ //looping over the rows
					non_zero = ia[i+1] - ia[i]; // row number of non zero
					int diag_idx=ia[i]-1;//diagonal idex
					a[diag_idx] = AM_1(i,i);   // diagonal value in proper position
					ja[diag_idx]=i+1;	// diagonal index 
					//	std::cout<< "We have "<< non_zero << " in line "<< i <<std::endl;
					//	 std::cout<< "going from "<< ia[i] << " to  "<< ia[i+1]<<"/"<<ja_v.size() <<std::endl;
					int shift=1;//shift for of diagonal node of the row	
					for(int j = 0  ; j < non_zero  ; j++ ){ //looping on non diagonal node
						//		std::cout<< "i   "<< i << " ia "<< ia[i] <<" j "<<  j << " ja_v " << ja_v[j] <<std::endl;
						if (i != ja_v[diag_idx + j]) // skip diagonal term
						{
							a[diag_idx + shift ]  = AM_1(i,ja_v[diag_idx + j]);// store off diagonal term
							ja[diag_idx+shift] =ja_v[diag_idx + j]+1;// store off diagonal index
							shift++;// increment shift
						}
					}
				}
				std::cout<< "Non zero = "<< a_v.size()<<std::endl;
				//	for (int i =0; i<nnu; i++)
				//		std::cout<< ia[i]<< std::endl;
				//	for (int i =0; i<nna; i++)
				//		std::cout<< ja[i]<< std::endl;
				// matrix=11:  A is symmetric and rowsum zero6 
				//	matrix=12:  A is symmetric and not rowsum zero 
				//	matrix=21:  A is non‐symmetric and rowsum zero 
				//	matrix=22:  A is non‐symmetric and not rowsum zero 
#endif //sparse interface

				APPL_INT npnt,nsys,matrix;
				matrix=22; nsys=1;npnt=0;
				APPL_INT ndiu      = 1;        // dimension of (dummy) vector iu
				APPL_INT ncyc;

				float told,tnew,tamg;

				//SAMG configuration
				//i

				 APPL_INT * iu;
			if (nsys==1) iu= new APPL_INT[1];
			else {
 			iu  = new APPL_INT[nnu];
 			ndiu   = nnu;
			 for(int iiu=0;iiu<nnu;iiu++){
				//only working for nsys = 2 check it out for larger systems				 
				if(iiu<dof_transp.Ct())iu[iiu]=1;
				else iu[iiu]=2;
				}//end for iiu
			}
				
			//=====================================================
			/// printing format
			char ch[] = "f";int len=1;
			SAMG_SET_IOFORM(&ch[0],&len);

			
			//==================================================			
			// int nrc	;
			// int nptmax; // max number of points for which we switch to nrc_emergency
			// int clsolver_finest;	// 
			// SAMG_GET_NRC(&nrc); // retreive nrc=solver for coarser levels page 67 user guide
			// std::cout<<"..value of nrc======"<<nrc<<std::endl;
			// ================================
			// 1  Iterative application of currently selected smoother 
			// 2  ILU(0) preconditioned CG 
			// 3  ILUT preconditioned BiCGstab 
			// 4  Diagonally preconditioned CG 
			// 5  Diagonally preconditioned BiCGstab 
			// 6  Dense direct solver (standard)  
			// 7  Sparse direct solver (standard)  
			// 8  Least squares solver (LSQ; robust but very expensive!) 
			// 9  Dummy interface routine (→ Section 8.3.1)  
			// 10  Dense direct solver (highly tuned, → cldense_ctrl below) 
			// 11  Sparse direct solver (Intel’s pardiso) if MKL has been linked 
			// 99  Another instance of SAMG 
			// ================================



			//==================================================	





#ifdef DIRECT_SOLVER	
	
			int levelx;
			SAMG_GET_LEVELX(&levelx); // retreive levelx=number of maximum coarsening levels 
			 std::cout<<"..value of levelx======"<<levelx<<std::endl;
			levelx=1;
			SAMG_SET_LEVELX(&levelx);// change levelx=number of maximum coarsening levels 
			SAMG_GET_LEVELX(&levelx);// retreive levelx=number of maximum coarsening levels 
			std::cout<<"..check if change value of levelx======"<<levelx<<std::endl;
			
			int  nrc=11;int  nrc_emergency=11;int nptmax=5000;int clsolver_finest=1;
			  SAMG_SET_NRC(&nrc);// change  nrc=solver for coarser levels page 67 user guide
			SAMG_SET_NRC_EMERGENCY(&nrc_emergency);
			  SAMG_SET_NPTMAX(&nptmax);// change  nptmax
			  SAMG_SET_CLSOLVER_FINEST(&clsolver_finest);// change  nptmax
			// SAMG_GET_NRC(&nrc);// retreive  nrc=solver for coarser levels page 67 user guide
			// std::cout<<"..check if change value of nrc======"<<nrc<<std::endl;
			
			ncyc      = 11050;    // V-cycle as pre-conditioner for CG; at most 50 iterations

#endif //direct_solver



#ifdef AMG_STAND_ALONE
			//Both ncgrad (the "nd number in ncyc) and ncgrad_default must be equal to 0 to use SAMG solver as a stand-alone solver (not as a preconditioner)
 				
			int ncgrad_default=0;
			SAMG_SET_NCGRAD_DEFAULT(&ncgrad_default);
			ncyc      = 10050;    // V-cycle as pre-conditioner for CG; at most 50 iterations
#endif //amg_stand_alone
			


#ifdef AMG_ACCELERATED

			ncyc      = 11050;    // V-cycle as pre-conditioner for CG; at most 50 iterations	
			
#endif //amg_accelerated


				APPL_INT ndip      = 1;        // dimension of (dummy) vector ip
				// APPL_INT * iscale = new APPL_INT[1];
				// this vector (iscale) indicates which uknowns require scaling if 0 no scaling
				APPL_INT * iscale = new APPL_INT[nsys]; for(int i_sys=0; i_sys<nsys; i_sys++) iscale[i_sys]=0; 
				APPL_INT * ip     = new APPL_INT[1];
				APPL_INT nsolve    = 2;        // results in scalar approach (current system is scalar)
				APPL_INT ifirst    = 1;        // first approximation = zero
				double  eps       = 1.0e-6;   // required (relative) residual reduction
				
				//ncyc      = 50050;
				APPL_INT n_default = 20;       // select default settings for secondary parameters

				// CURRENTLY AVAILABLE: 10-13, 15-18, 20-23, 25-28
				// NOTE: the higher the respective SECOND digit, the
				// more aggressive the coarsening (--> lower memory at
				// the expense of slower convergence)
				APPL_INT iswtch    = 5100+n_default; // complete SAMG run ....
				// ... memory de-allocation upon return ....
				// ... memory extension feature activated ....
				// ... residuals measured in the L2-norm .....
				// ... secondary parameter default setting # n_default

				// ===> Secondary parameters which have to be set if n_default=0
				//      (at the same time a demonstration of how to access secondary or hidden parameters)

				double a_cmplx   = 2.2;      // estimated dimensioning
				double g_cmplx   = 1.7;      // estimated dimensioning
				double w_avrge   = 2.4;      // estimated dimensioning
				double p_cmplx   = 0.0;      // estimated dimensioning (irrelevant for scalar case)
				double  chktol    = -1.0;    // input checking de-activated (we know it's ok!)
				
				//============================================
				//     idump controls the matrix dumping of SAMG				
				// 1  Standard print output, no matrix dump. 
				// 2‐6  Write matrices to disk: level 2 up to level idmp. 
				// 7  Write matrices to disk: level 2 up to the coarsest level. 
				// 8  Write finest‐level matrix to disk (incl. right hand side etc.). 
				// 9  Write all matrices to disk. 
				APPL_INT idump     = 1;       // minimum output during setup
				//============================================
				// iout page 44 Userguide. it controls display outpu. default 2 very verbose 43
				APPL_INT iout      = 2;        // display residuals per iteration and work statistics



				// output:
				APPL_INT ierr,ierrl,ncyc_done;
				double res_out,res_in;
				//end

				// for(int xx=0;xx<200;xx++)std::cout<<a_samg[xx]<<" "<<a[xx]<< " -- "<<ja_samg[xx]<<" "<<ja[xx] <<std::endl;
				SAMG_CTIME(&told);
				double time_SAMG=gmm::uclock_sec();
				std::cout<<"start solving with samg " << std::endl;
#ifdef SPARSE_INTERFACE  
				SAMG(&nnu,&nna,&nsys,
						&ia[0],&ja[0],&a[0],&f[0],&u[0],&iu[0],&ndiu,&ip[0],&ndip,&matrix,&iscale[0],
						&res_in,&res_out,&ncyc_done,&ierr,
						&nsolve,&ifirst,&eps,&ncyc,&iswtch,
						&a_cmplx,&g_cmplx,&p_cmplx,&w_avrge,
						&chktol,&idump,&iout);
#endif
#ifdef CSC_INTERFACE

				SAMG(&nnu,&nna,&nsys,
						&ia_samg[0],&ja_samg[0],&a_samg[0],&b_samg[0],&u_samg[0],&iu[0],&ndiu,&ip[0],&ndip,&matrix,&iscale[0],
						&res_in,&res_out,&ncyc_done,&ierr,
						&nsolve,&ifirst,&eps,&ncyc,&iswtch,
						&a_cmplx,&g_cmplx,&p_cmplx,&w_avrge,
						&chktol,&idump,&iout);
 
 // std::cout << std::string(90, '=') << std::endl;
 // std::cout << "New samg iteration" << std::endl;
 // std::cout << std::string(90, '=') << std::endl;
// clsolver_finest=111;levelx=1;
//			SAMG_SET_LEVELX(&levelx);// change levelx=number of maximum coarsening levels 

#endif

#ifdef CSR_INTERFACE

				SAMG(&nnu,&nna,&nsys,
						&ia_samg[0],&ja_samg[0],&a_samg[0],&b_samg[0],&u_samg[0],&iu[0],&ndiu,&ip[0],&ndip,&matrix,&iscale[0],
						&res_in,&res_out,&ncyc_done,&ierr,
						&nsolve,&ifirst,&eps,&ncyc,&iswtch,
						&a_cmplx,&g_cmplx,&p_cmplx,&w_avrge,
						&chktol,&idump,&iout);
 
 // std::cout << std::string(90, '=') << std::endl;
 // std::cout << "New samg iteration" << std::endl;
 // std::cout << std::string(90, '=') << std::endl;
// clsolver_finest=111;levelx=1;
//			SAMG_SET_LEVELX(&levelx);// change levelx=number of maximum coarsening levels 

#endif


cout << endl<<"*** time to solve the system using SAMG with gmm time: " << gmm::uclock_sec() - time_SAMG << " seconds\n";	
	cout << endl<<"... time to solve the system using SAMG with assembling********= " << gmm::uclock_sec() - time_ass<< " seconds\n";	

				if (ierr > 0) {
					cout << endl << " SAMG terminated with error code " 
						<< ierr << " **** " << endl;
				}
				else if (ierr < 0) {
					cout << endl << " SAMG terminated with warning code " 
						<< ierr << " **** " << endl;
				}

				SAMG_CTIME(&tnew);
				tamg=tnew-told;
				cout << endl << " ***** total run time: " << tamg << " ***** " << endl; 
				cout << endl << " ***** total run time with assembling: " << tnew-told_ass << " ***** " << endl; 






#ifdef SPARSE_INTERFACE
				for(int i = 0 ; i < nrows ; i++ ){U_1[i]=u[i];
					UM_transp[i]=u[i];	
				}
				gmm::copy(U_1, UM_transp);
				if(FIXED_SOURCE==1){
					vector_type x(dof_transp.Cv());
					vector_type Source_coef(dof_transp.Cv());
					for(scalar_type k=0;k<dof_transp.Cv();k++){
						x[k]=(k+1)*(1.0/(dof_transp.Cv()));
						Source_coef[k]=1.0-x[k];
					}
					gmm::clear(x);	
					for(scalar_type k=0;k<dof_transp.Cv();k++)
					UM_transp[dof_transp.Ct()+k]=Source_coef[k];
				}
#endif
#ifdef CSC_INTERFACE
				for(int i = 0 ; i < nnu ; i++ ){
					U_2[i]=u_samg[i];UM_transp[i]=u_samg[i];}
				gmm::copy(U_2,UM_transp);
				if(FIXED_SOURCE==1){
					vector_type x(dof_transp.Cv());
					vector_type Source_coef(dof_transp.Cv());
					for(scalar_type k=0;k<dof_transp.Cv();k++){
						x[k]=(k+1)*(1.0/(dof_transp.Cv()));
						Source_coef[k]=1.0-x[k];
					}
					gmm::clear(x);	
					for(scalar_type k=0;k<dof_transp.Cv();k++)
					UM_transp[dof_transp.Ct()+k]=Source_coef[k];
				}
#endif
				
#ifdef CSR_INTERFACE
				for(int i = 0 ; i < nnu ; i++ ){
					U_2[i]=u_samg[i];UM_transp[i]=u_samg[i];}
				gmm::copy(U_2,UM_transp);
				if(FIXED_SOURCE==1){
					vector_type x(dof_transp.Cv());
					vector_type Source_coef(dof_transp.Cv());
					for(scalar_type k=0;k<dof_transp.Cv();k++){
						x[k]=(k+1)*(1.0/(dof_transp.Cv()));
						Source_coef[k]=1.0-x[k];
					}
					gmm::clear(x);	
					for(scalar_type k=0;k<dof_transp.Cv();k++)
					UM_transp[dof_transp.Ct()+k]=Source_coef[k];
				}
#endif
	

				//export solution
				std::cout<<"solved! going to export..."<<std::endl;
				string time_suff = "";
				std::ostringstream convert;
				convert << time_count;
				time_suff = convert.str();
				//Compute the total c_t in the tissue

export_vtk("1"); 
export_vtk("2"); 

			// 	if(t==0){
			// 		export_vtk(time_suff); 
			// 		std::cout<<"exported! now new iteration..."<<std::endl;
			// 	}
			// 	//if(iteraz%30==0){ 
			// 	if(iteraz%10==0){ 
			// 		export_vtk(time_suff); 
			// 	}

			} //end of cycle over time 


			cout << endl<<"... time to solve : " << gmm::uclock_sec() - time << " seconds\n";
#endif //USE_SAMG
			return true;
			}; // end of solve


















			void transport3d1d::export_vtk (const string & time_suff,const string & suff)
			{
				if (PARAM.int_value("VTK_EXPORT"))
				{
#ifdef M3D1D_VERBOSE_
					cout << "Exporting the solution (vtk format) to " << descr.OUTPUT << " ..." << endl;
#endif
#ifdef M3D1D_VERBOSE_
					cout << "  Saving the results from the monolithic unknown vector ... " << endl;
#endif

					// Array of unknown dof of the interstitial velocity
					vector_type Ct(dof_transp.Ct()); 

					// Array of unknown dof of the network velocity
					vector_type Cv(dof_transp.Cv()); 

					//Copy solution
					gmm::copy(gmm::sub_vector(UM_transp, 
								gmm::sub_interval(0, dof_transp.Ct())), Ct);
					gmm::copy(gmm::sub_vector(UM_transp, 
								gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())), Cv);


#ifdef M3D1D_VERBOSE_
					cout << "  Exporting Ct ..." << endl;
#endif
					vtk_export exp_Ct(descr_transp.OUTPUT+"Ct"+suff+"_t"+time_suff+".vtk");
					exp_Ct.exporting(mf_Ct);
					exp_Ct.write_mesh();
					exp_Ct.write_point_data(mf_Ct, Ct, "Ct");



#ifdef M3D1D_VERBOSE_
					cout << "  Exporting Cv ..." << endl;
#endif
					vtk_export exp_Cv(descr_transp.OUTPUT+"Cv"+suff+"_t"+time_suff+".vtk");
					exp_Cv.exporting(mf_Cv);
					exp_Cv.write_mesh();
					exp_Cv.write_point_data(mf_Cv, Cv, "Cv");



					//esporto le concentrazioni su un file .txt
					//std::ofstream outCv("Cv"+suff+"_t"+time_suff+".txt");
					//outCv << gmm::col_vector(Cv);
					//outCv.close();



#ifdef M3D1D_VERBOSE_
					cout << "... export done, visualize the data file with (for example) Paraview " << endl; 
#endif
				}
			}; // end of export


			void transport3d1d::test(void){

				// risolvi il problema: (At + Btt)Ct = Btv Cv
				// stazionario
				//con Cv considerato costante e uguale a 1 su tutta la mesh
				//At è la matrice di stiffness per il laplaciano
				//Btt e Btv sono le matrici di accoppiamento tra le mesh 3d - 1d


				//definizioni di matrici e vettori
				vector_type Ct(dof_transp.Ct()); gmm::clear(Ct);
				sparse_matrix_type At(dof_transp.Ct(), dof_transp.Ct());gmm::clear(At);
				sparse_matrix_type Btt(dof_transp.Ct(), dof_transp.Ct());gmm::clear(Btt);
				sparse_matrix_type Btv(dof_transp.Ct(), dof_transp.Cv());gmm::clear(Btv);
				sparse_matrix_type Bvt(dof_transp.Cv(), dof_transp.Ct());gmm::clear(Bvt); // non servono, ma devo comunque passarle alla funzione asm_exchange_mat
				sparse_matrix_type Bvv(dof_transp.Cv(), dof_transp.Cv());gmm::clear(Bvv); //	
				vector_type Cv(dof_transp.Cv(), 1.0);
				vector_type F(dof_transp.Ct()); gmm::clear(F);  //F= Btv Cv

				//stiffness matrix
				getfem::generic_assembly
					assem("M$1(#1,#1) += sym(comp(vGrad(#1).vGrad(#1)) (:, i,k, : ,i,k) )");
				assem.push_mi(mimt);
				assem.push_mf(mf_Ct);
				assem.push_mat(At);
				assem.assembly();

				//matrici di accoppiamento
				sparse_matrix_type Mbar(dof_transp.Cv(), dof_transp.Ct());gmm::clear(Mbar);
				sparse_matrix_type Mlin(dof_transp.Cv(), dof_transp.Ct());gmm::clear(Mlin);
				asm_exchange_aux_mat(Mbar, Mlin, mimv, mf_Ct, mf_Cv, param.R(), descr.NInt);

				bool NEWFORM = true;	
				asm_exchange_mat_transp(Btt, Btv, Bvt, Bvv,
						mimv, mf_Cv, mf_coefv, Mbar, Mlin, param_transp.Y(), NEWFORM);

				// Copying Btt and add it to At
				gmm::add(Btt, At); 

				// F è il termine noto, F= Btv*Cv
				gmm::mult(Btv, Cv, F);

				//costruisci una matrice del tipo adatto per SuperLU_solver
				gmm::csc_matrix<scalar_type> A;
				gmm::clean(At, 1E-12); 
				gmm::copy(At, A);

				gmm::MatrixMarket_IO::write("At.mm",At);
				gmm::MatrixMarket_IO::write("Btt.mm",Btt);
				gmm::MatrixMarket_IO::write("Btv.mm",Btv);
				gmm::MatrixMarket_IO::write("Bvt.mm",Bvt);
				gmm::MatrixMarket_IO::write("Bvv.mm",Bvv);	


				//std::ofstream outF("F.txt");
				//outF << gmm::col_vector(F);
				//outF.close();	

				int a;
				cout<<"type 0 for superLU; 1 for GMRES"<<endl;
				cin>>a;

				if(a==0){
					cout<<"	superLU"<<endl;
					scalar_type cond;
					gmm::SuperLU_solve(At, Ct, F, cond);
					cout << "  Condition number (test diffusion problem): " << cond << endl;}
				else if(a==1){
					cout<<"	GMRES"<<endl;
					gmm::iteration iter (0.0000001, 1,40000);
					//gmm::ilutp_precond<sparse_matrix_type> P(At, 20, 1E-6);
					gmm::identity_matrix P;
					gmm::gmres (At, Ct, F, P, 50, iter);

				}



				//std::ofstream outUU("Ct.txt");
				//outUU << gmm::col_vector(Ct);
				//outUU.close();		

				vtk_export exp_Ct("test_Ct.vtk");
				exp_Ct.exporting(mf_Ct);
				exp_Ct.write_mesh();
				exp_Ct.write_point_data(mf_Ct, Ct, "Ct");


	};//end of test

} // end of namespace
