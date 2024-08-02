#include "GFMC.h"
#include "file.h"

GFMC::GFMC(const MyMpi &mm_i, const Pars &p_i, Model &mdl_i)
	: mm(mm_i)
	, p(p_i)
	, mdl(mdl_i)
	, ifreconfig(p.ifreconfig)
	, Nw(p.Nw)
	, Nw_per_prog(p.Nw_per_prog)
	, m_walkersctrl(mm, p, mdl_i, Nw_per_prog)
	, rng((unsigned)time(NULL)+mm.id())
	, lnw_avg_estimation(0.0)
	, curr_lnw_lt(Nw)
	, curr_w_lt(Nw)
{
	if (!ifreconfig)
	{
		if (1 != Nw_per_prog)
		{
			ERR("GFMC::GFMC for the non-reconfigration case, there is only one walker for each process!");
		}

		curr_lnw_lt.reset(1);
		curr_w_lt.reset(1);
	}
}

GFMC::~GFMC(){}

void GFMC::backup_cfg_pars(const Str &file_name_append) const
{
	m_walkersctrl.backup_cfg_pars();
}

void GFMC::walk()
{ 
	print_stat_label();
	//preparation and thermalization
	Int n_thermal_step = p.nStep / 10;
	for_Int(i, 0, 10/*n_thermal_step*/)
	{
		if(mm) std::cout << "thermal_step: " << i << std::endl;
		reconfig();
	}
	//Estimate the average value of global weight
	calc_lnw_avg_estimation();
	lnw_avg_lt.clear();
	w_lt.clear();
	//sampling process
	Int n_total_step = p.nP + p.nStep;
	Int n_stat_id = 0;
	Int n_ss_total = 0;

	m_walkersctrl.reset_n_index();

	//add 20200520
	for_Int(i, 0, n_total_step)
	{
		if (mm) std::cout << "total_step: " << i << std::endl;
		reconfig();
		if (i < p.nP) continue;
		calc_Gn();
		calc_eloc();
		print_stat(n_stat_id++);
	}

	calc_bin_avg();
	calc_energy_avg();
}


void GFMC::reconfig()
{
	m_walkersctrl.reset_lnw();
	m_walkersctrl.propagate();

	curr_lnw_lt = m_walkersctrl.gather_lnw() - VecReal(curr_lnw_lt.size(), lnw_avg_estimation);
	curr_w_lt = EXP(curr_lnw_lt);

	Real lnw_avg = AVG(curr_lnw_lt);
	lnw_avg_lt.push_back(lnw_avg);
	w_lt.push_back(AVG(curr_w_lt));

	population_ctrl();
}

void GFMC::calc_lnw_avg_estimation()
{
	lnw_avg_estimation = 0.;
	for_Int(i, 0, lnw_avg_lt.size())
	{
		lnw_avg_estimation += lnw_avg_lt[i];
	}

	lnw_avg_estimation *= INV(lnw_avg_lt.size());


	VecReal vec_lnw_avg(mm.np());
	vec_lnw_avg = mm.Allgather(lnw_avg_estimation);
	lnw_avg_estimation = AVG(vec_lnw_avg);
}

void GFMC::calc_eloc()
{
	VecReal vec_eloc(m_walkersctrl.gather_eloc());

	if (curr_w_lt.size() != vec_eloc.size()) ERR("void GFMC::calc_eloc(): size error!");

	Real top = 0., bottom = 0.;
	for_Int(i, 0, curr_w_lt.size())
	{
		top += curr_w_lt[i] * vec_eloc[i];
		bottom += curr_w_lt[i];
	}

	eloc_lt.push_back(top / bottom);
	
}


void GFMC::calc_Gn()
{
	if (w_lt.size() < p.nP) ERR("void GFMC::calc_Gn()---w_lt.size() error!!!");


	if (w_lt.size() >= p.nP && Gn_lt.empty())
	{
		Real Gn = 1.;
		for (Int i = 0; i < p.nP; i++)
		{
			Gn *= w_lt[w_lt.size() - (i + 1)];
		}
		Gn_lt.push_back(Gn);
	}
	else
	{
		Real Gn = Gn_lt.back();
		Gn *= w_lt.back();
		Gn /= w_lt[w_lt.size() - (p.nP+1)];
		Gn_lt.push_back(Gn);
	}
}

void GFMC::population_ctrl()
{
	VecReal vec_w_prob( INV(SUM(curr_w_lt)) * curr_w_lt );

	//get the cfg vector
	Vec<OccCfg> vec_total_cfg = m_walkersctrl.gather_cfg();

	//generates a random number between [0,1)
	Real random = 0.;
	
	for_Int(i, 0, Nw)
	{
		random = rng();
		mm.Bcast(random);

		Int ch_id = 0;
		Real sum = 0.;
		Real Z = (i + random) / Nw;
		for_Int(j, 0, vec_w_prob.size())
		{
			sum += vec_w_prob[j];
			if (sum > Z)
			{
				ch_id = j;
				break;
			}
		}

		if (i == ch_id)
		{
			continue;
		}
		else if (if_include_current_thread(i))
		{
			//walker include in this progress
			set_new_cfg(vec_total_cfg[ch_id], i % Nw_per_prog);
		}
	}

	//debug
	//check_num_of_same_cfg();

}

bool GFMC::if_include_current_thread(Int index)
{
	Int pro_id  = mm.id();
	if (index >= pro_id * Nw_per_prog && index < (pro_id+1) * Nw_per_prog)
	{
		return true;
	}
	return false;
}

void GFMC::set_new_cfg( OccCfg cfg_in, Int pos_in)
{
	m_walkersctrl.set_new_cfg(cfg_in, pos_in);
}


void GFMC::check_num_of_same_cfg()
{
	vector<OccCfg> vRes;
	vector<Real> vReal;
	//get the cfg vector
	Vec<OccCfg> vec_total_cfg = m_walkersctrl.gather_cfg();
	for_Int(i, 0, vec_total_cfg.size())
	{
		OccCfg& cfg = vec_total_cfg[i];
		bool flag = true;
		for_Int(k, 0, vRes.size())
		{
			if (cfg == vRes[k])
			{
				vReal[k] += 1;
				flag = false;
				break;
			}
		}

		if (flag) {
			vReal.push_back(1.);
			vRes.push_back(cfg);
		}
	}

	Real total = 0.;
	for_Int(i, 0, vReal.size())
	{
		if(mm) std::cout << vReal[i] << std::endl;
		//total += vReal[i];
	}
	//if (mm) std::cout << "total: " << total << std::endl;
}

void GFMC::print_stat_label(std::ostream &os)
{
	using namespace std;
	Int w = TOSTRLEN(p.nStep);

	if (mm) os << setw(w) << "i" << " >  "
		<< std::uppercase << std::scientific << std::setprecision(6)
		<< setw(w_Real_def) << "Gp" << "  "
		<< setw(w_Real_def) << "Eloc" << "  "
		<< present() << endl << iofmt("def");
}

void GFMC::print_stat(Int itr, std::ostream &os)
{
	using namespace std;
	const Int w = TOSTRLEN(p.nStep);
	Real eloc = eloc_lt[itr];
	Real Gn = Gn_lt[itr];

	if (mm) os << setw(w) << itr << " >  "
		<< std::uppercase << std::scientific << std::setprecision(6)
		<< setw(w_Real_def) << Gn << "  "
		<< setw(w_Real_def) << eloc * INV(p.ns) << "  "
		<< present() << endl << iofmt("def");
}

void GFMC::calc_bin_avg()
{
	if (eloc_lt.size() != Gn_lt.size())
	{
		ERR("GFMC::compute_energy(void) : size error !");
	}

	Int nN = p.nStep / p.nbin;
	if ( nN <= 1)
	{
		ERR("GFMC::compute_energy(void) : bin size error !");
	}

	for_Int(i, 0, p.nbin)
	{
		Real topEnergy = 0., buttom = 0.;

		for_Int(j, 0, nN)
		{
			Int index = i * nN + j;
			topEnergy += eloc_lt[index] * Gn_lt[index];
			buttom += Gn_lt[index];
		}
		topEnergy *= INV(buttom);

		eloc_bin_lt.push_back(topEnergy);
		eloc_sqr_bin_lt.push_back(SQR(topEnergy));
	}
}

Real GFMC::calc_Energy_target()
{
	if (eloc_lt.size() != Gn_lt.size())
	{
		ERR("GFMC::compute_energy(void) : size error !");
	}

	Real topEnergy = 0., buttom = 0.;

	for_Int(i, 0, eloc_lt.size())
	{
		topEnergy += eloc_lt[i] * Gn_lt[i];
		buttom += Gn_lt[i];
	}

	if (!ifreconfig)
	{
		topEnergy = SUM(mm.Allgather(topEnergy));
		buttom = SUM(mm.Allgather(buttom));
	}

	return topEnergy * INV(buttom);;
}

void GFMC::calc_energy_avg()
{
	if (eloc_bin_lt.size() != p.nbin)
	{
		ERR("GFMC::compute_energy(void) : size error !");
	}
	Real bin = p.nbin;
	Real E_avg = AVG(VecReal(eloc_bin_lt));
	Real E_sqr_avg = AVG(VecReal(eloc_sqr_bin_lt));

	if (!ifreconfig)
	{
		E_avg = AVG(mm.Allgather(E_avg));
		E_sqr_avg = AVG(mm.Allgather(E_sqr_avg));
		bin = SUM(mm.Allgather(bin));
	}

	//D(x) = E(x^2) - (E(x))^2
	Real erro_bar = SQRT(ABS(E_sqr_avg - nsqr(E_avg)) * INV(bin));

	Real Energy = calc_Energy_target();

	//print and output to file
	if (mm) std::cout << std::uppercase<< std::scientific << std::setprecision(6)
		<< setw(20) << "Enegry:"	<< setw(w_Real_def) << Energy * INV(p.ns)
		<< setw(20) << "error bar:" << setw(w_Real_def) << erro_bar * INV(p.ns) << "\n"
		<< std::endl;
}
