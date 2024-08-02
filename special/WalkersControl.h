#ifndef WALKERS_CONTROL_H
#define WALKERS_CONTROL_H

#include <vector>
#include "Walker.h"

class WalkersControl
{
   public:

	   WalkersControl(const MyMpi &mm_i, const Pars &p_i, Model &mdl_i, Int n_w_i);
	   ~WalkersControl();
	   void init_walkers_lt();
	   void reset_lnw();
	   void propagate();
	   VecReal gather_lnw();
	   VecReal gather_eloc();
	   Vec<OccCfg> gather_cfg();
	   void set_new_cfg(OccCfg cfg_in, Int index);
	   void backup_cfg_pars(const Str &file_name_append = empty_str) const;
	   void reset_n_index();
   private:

		const MyMpi &mm;
		const Pars  &p;
		Model &mdl;
		Int n_w;
		Bool ifreconfig;

		std::vector< Walker* > walkers_lt;
		Real w_sum;
		Real eloc_sum;
		VecReal ss_sum;
		Int all_ss_bonds_num;
		Int ns;
};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
