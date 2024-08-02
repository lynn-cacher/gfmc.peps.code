#ifndef WALKER_H
#define WALKER_H

#include <Python.h>
#include "vector"
#include "spcfcs.h"
#include "pars.h"
#include "Hamiltonian.h"


//save the Intermediate process information of spin correlation
class Walker 
{
   public:

      Walker(const MyMpi &mm_i, const Pars &p_i, Model &mdl_i);
      
	  ~Walker();

	  //operate the weight
      Real get_lnw() const;

      void set_lnw(Real);

	  //operate the config
	  OccCfg get_cfg() const;

	  void set_cfg( OccCfg cfg );

	  void init_cfg(Int nId = 0);

	  bool read_cfg(OccCfg &cfg_read) const;

	  void write_cfg() const;

	  void adjust_cfg();

	  Int rng_Int();

	  //evolution of the walker
	  void calc_Hxx_Psi();

	  void calc_Hfnxx();

	  void calc_tau();
	 
	  void get_txx_lt();

	  Int	draw();

	  void get_new_cfg(OccCfg &cfg_i, VEC< VEC<Int> >& flips);

	  Real get_factor();

	  void propagate();
	  
	  Real get_eloc() const;

	  StdVecReal get_psip_o_psi(const OccCfg& cfg, const VEC<OccCfg> &cfg_lt);

	  //20230614
	  Real save_cfg_w_eloc();

	  void reset_n_index();

  private:
	  const MyMpi &mm;
	  const Pars  &p;
	  Model &mdl;
	  PRng prng;
	  Int nv;
	  Int ne;
	  bool ifnec;
      Real lnw;
	  OccCfg cfg;

	  vector<OccCfg> cfg_lt;
	  StdVecReal hxx_lt;
	  StdVecReal psixp_o_psix_lt;
	  StdVecReal hfnxx_lt;
	  StdVecReal pxx_lt;
	  StdVecReal txx_lt;
	  Real beta;
	  Real beta_left;
	  Real tau;
	  Real gamma;
	  Real e_loc;
	  Int n_cfg_id;
	  Int n_index;
	  PyObject* pModule;
	  PyObject* pDict;
	  PyObject* pClass;
	  PyObject* pInstance;
	  PyObject* result;
};

#endif
