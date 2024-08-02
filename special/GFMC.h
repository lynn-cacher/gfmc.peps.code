#ifndef GFMC_H
#define GFMC_H

#include "include.h"
#include "spcfcs.h"
#include "Hamiltonian.h"

class GFMC
{

   public:
   
		//constructor with input trialwavefunction
		GFMC(const MyMpi &mm_i, const Pars &p_i, Model &mdl_i);
      
		//Destructor
		~GFMC();
      
		//Let the walkers propagate for n_steps steps
		void walk();

		void reconfig();

		void population_ctrl();

		void calc_eloc();

		void calc_Gn();

		void calc_energy_avg();

		bool if_include_current_thread(Int index);

		void set_new_cfg(OccCfg cfg_in, Int pos_in);

		void print_stat_label(std::ostream &os = std::cout);
		void print_stat(Int itr, std::ostream &os = std::cout);

		void backup_cfg_pars(const Str &file_name_append = empty_str) const;


		//for test
		void calc_lnw_avg_estimation();

		void calc_bin_avg();

		Real calc_Energy_target();

		void check_num_of_same_cfg();

   private:

		const MyMpi &mm;
		const Pars  &p;
		Model &mdl;

		//!The total desired number of walkers
		Int Nw;
		Int Nw_per_prog;
		UniformUnitRandomSFMT rng;
		Bool ifreconfig;

		//!Distribution of possible final states given a walker state
		WalkersControl m_walkersctrl;

		//Nw weights for current step
		VecReal curr_lnw_lt;
		VecReal curr_w_lt;

		//nP average weights
		vector< Real > w_lt;
		vector< Real > lnw_avg_lt;

		vector< Real > Gn_lt;
		vector< Real > eloc_lt;

		vector< Real > eloc_bin_lt;
		vector< Real > eloc_sqr_bin_lt;

		vector<VecReal> ss_lt;

		vector<VecReal> ss_bin_lt;
		vector<VecReal> ss_sqr_bin_lt;

		//In oder to reduce G, we have to divide w each step by r_w_avg_estimation
		Real  lnw_avg_estimation;


		//for test
		vector<Real> lnw_var_lt;
		vector<Real> beta_one_flip_lt;

		//20230414
		vector<VecReal> op_lt;
		vector<VecReal> op_bin_lt;
		vector<VecReal> op_sqr_bin_lt;

		//20230414
		vector<VecReal> dimer_lt;
		vector<VecReal> dimer_bin_lt;
		vector<VecReal> dimer_sqr_bin_lt;

		//20230613
		vector< OccCfg > cfg_samples_lt;
		vector< Real > w_samples_lt;
		vector< Real > eloc_samples_lt;

		//20231031
		vector<VecReal> all_ss_lt;
		vector<VecReal> all_ss_bin_lt;
		vector<VecReal> all_ss_sqr_bin_lt;

		//20231211
		vector<VecReal> SziSzj_for_different_k_lt;
		vector<VecReal> SziSzj_for_different_k_bin_lt;
		vector<VecReal> SziSzj_for_different_k_sqr_bin_lt;

		vector<VecReal> dimer_for_different_k_lt;
		vector<VecReal> dimer_for_different_k_bin_lt;
		vector<VecReal> dimer_for_different_k_sqr_bin_lt;

		vector<VecReal> get_Si_lt;
		vector<VecReal> get_Si_bin_lt;
		vector<VecReal> get_Si_sqr_bin_lt;


};

#endif
