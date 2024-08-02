#include "Walker.h"

Walker::Walker(const MyMpi &mm_i, const Pars &p_i, Model &mdl_i)
	: mm(mm_i)
	, p(p_i)
	, mdl(mdl_i)
	, prng(mm_i)
	, nv(p_i.ns)
	, cfg(nv)
	, ifnec(p_i.ifnec)
	, beta(p_i.beta)
	, beta_left(p_i.beta)
	, gamma(p_i.gamma)
	, ne(p_i.ne)
	, lnw(0.0)
	, tau(0.0)
	, n_cfg_id(0)
	, n_index(0)
{
	Py_Initialize();

	pModule = NULL;
	pDict = NULL;
	pClass = NULL;
	pInstance = NULL;
	result = NULL;

	PyRun_SimpleString("import sys");
	PyRun_SimpleString("sys.path.append('./peps/')");
	pModule = PyImport_ImportModule("amplitute");
	if( pModule == NULL )
	{
		cout <<"module not found" << endl;
		exit(0);
	}
	pDict = PyModule_GetDict(pModule);
	if( pDict == NULL )
	{
		cout <<"dictionary not found" << endl;
		exit(0);
	}
	pClass = PyDict_GetItemString(pDict, "Amplitute");
	if( pClass == NULL )
	{
		cout <<"class not found" << endl;
		exit(0);
	}

	pInstance = PyObject_CallObject(pClass, NULL);
	if( pInstance == NULL )
	{
		cout <<"instance not found" << endl;
		exit(0);
	}
}

Walker::~Walker()
{
	Py_Finalize();
}

void Walker::set_lnw(Real lnw_i)
{
	lnw = lnw_i;
}

void Walker::set_cfg(OccCfg cfg_i)
{
	cfg = cfg_i;
}

Int Walker::rng_Int()
{ 
	return Int(prng() * nv) % nv; 
}

void Walker::init_cfg(Int nId_i)
{
	n_cfg_id = nId_i;
	for_Int(i, 0, nv)
	{
		cfg[i] = prng() < 0.5 ? 0 : 1;
	}
	
	if (ifnec) adjust_cfg();

	OccCfg cfg_read(cfg);
	if (read_cfg(cfg_read))
	{
		cfg = cfg_read;
	}
}

bool Walker::read_cfg(OccCfg &cfg_read) const
{
	Str file = bix + "cfg" + STR(n_cfg_id);
	IFS ifs(file.c_str(), std::ios::binary);
	if (ifs) {
		ifs.read(cfg_read.p(), cfg.szof());
	}
	return mm.Allreduce(Int(!ifs.fail())) == mm.np();
}

void Walker::write_cfg() const
{
	Str file = bix + "cfg_" + STR(n_cfg_id) + "_" + STR(n_index);
	OFS ofs(file.c_str(), std::ios::binary);
	if (!ofs) ERR(NAV1(file));
	ofs.write(cfg.p(), cfg.szof());
	ofs.close();
}

void Walker::adjust_cfg()
{
	Int ne_cfg = 0;
	for_Int(i, 0, nv) ne_cfg += cfg[i];
	while (ne_cfg > ne) {
		Int rs = rng_Int();
		while (cfg[rs] == 0) {
			rs = rng_Int();
		}
		cfg[rs] = 0;
		ne_cfg -= 1;
	}
	while (ne_cfg < ne) {
		Int rs = rng_Int();
		while (cfg[rs] == 1) {
			rs = rng_Int();
		}
		cfg[rs] = 1;
		ne_cfg += 1;
	}
#ifdef _ASSERTION_
	ASSERT_EQ(ne_cfg, ne);
#endif
}

void Walker::get_new_cfg(OccCfg &cfg_i, VEC< VEC<Int> >& flips)
{
	cfg_lt.clear();

	for_Int(i, 0, flips.size())
	{
		OccCfg cfg = cfg_i;
		VEC<Int> &fp = flips[i];
		for_Int(j, 0, fp.size()) cfg[fp[j]] ^= 1;
		cfg_lt.push_back(cfg);
	}
}


void Walker::calc_Hxx_Psi()
{
	hxx_lt.clear();
	psixp_o_psix_lt.clear();

	vector<StdVecInt > conn;

 	mdl.find_conn(cfg, hxx_lt, conn);
	get_new_cfg(cfg, conn);

	psixp_o_psix_lt = get_psip_o_psi(cfg, cfg_lt);
}

StdVecReal Walker::get_psip_o_psi(const OccCfg& cfg, const VEC<OccCfg> &cfg_lt)
{
	Real psi_x;
	result = PyObject_CallMethod(pInstance, "get_amplitute", "s", OccCfg2String(cfg).c_str());
	PyArg_Parse(result, "d", &psi_x);

	const Int nconn = cfg_lt.size();
	StdVecReal psi_xp(nconn);
	for_Int(k, 0, nconn)
	{
		Real res;
		result = PyObject_CallMethod(pInstance, "get_amplitute", "s", OccCfg2String(cfg_lt[k]).c_str());
		PyArg_Parse(result, "d", &res);
		psi_xp[k] = res / psi_x;
	}
	return psi_xp;
}


void Walker::calc_Hfnxx()
{
	if (hxx_lt.size() != psixp_o_psix_lt.size())
	{
		ERR("void Walker::calc_Hfnxx()");
	}
	hfnxx_lt.clear();
	hfnxx_lt.resize(hxx_lt.size());
	Real sxx = 0, Vsf = 0.;
	for_Int(i, 1, hfnxx_lt.size())
	{
		sxx = hxx_lt[i] * psixp_o_psix_lt[i];

		if (sxx > 0)
		{
			hfnxx_lt[i] = -1.* gamma * sxx;
			Vsf += sxx;
		}
		else
		{
			hfnxx_lt[i] = sxx;
		}
	}
	hfnxx_lt[0] = hxx_lt[0] + (1 + gamma) * Vsf;
}

void Walker::calc_tau(void)
{
	Real xi = prng();

	e_loc = 0.;
	for_Int(i, 0, hfnxx_lt.size())
	{
		e_loc += hfnxx_lt[i];
	}

	Real temp = std::log(1.- xi) / (e_loc - hfnxx_lt[0]);

	if (beta_left < temp)
	{
		tau = beta_left;
	}
	else
	{
		tau = temp;
	}
}

void Walker::get_txx_lt()
{
	txx_lt.clear();
	txx_lt.resize(cfg_lt.size());
	txx_lt[0] = 0.;
	for_Int(i, 1, txx_lt.size())
	{
		txx_lt[i] = hfnxx_lt[i] / (e_loc - hfnxx_lt[0]);
	}
}

Int Walker::draw()
{
	get_txx_lt();

	//random number
	Real rd = prng();
	Real sum = 0.;
	for_Int(i, 0, txx_lt.size())
	{
		sum += txx_lt[i];
		if (sum > rd)
		{
			return i;
		}
	}
	return 0;
}


Real Walker::get_factor()
{
	return -tau * e_loc;
}

void Walker::propagate()
{
	beta_left = beta;

	Int icout = 0;
	//return until beta_left=0
	//std::cout << "propagate" << std::endl;
	while (true)
	{
		icout++;
		calc_Hxx_Psi();
		calc_Hfnxx();

		//get the time of evolution
		calc_tau();
		beta_left = beta_left - tau;
		//update the weight
		lnw += get_factor();
		if(0 >= beta_left) break;

		//get config of next step
		Int idx = draw();
		set_cfg(cfg_lt[idx]);
	}

	save_cfg_w_eloc();
	n_index += 1;
}

Real Walker::get_lnw() const
{
	return lnw;
}

Real Walker::get_eloc() const
{
	return e_loc;
}

OccCfg Walker::get_cfg() const
{
	return cfg;
}


//20230614
Real Walker::save_cfg_w_eloc()
{
	if (!p.if_savesamplings)
	{
		return 0;
	}
	write_cfg();
	
	using namespace std;
	Str file = iox + "w_eloc_" + STR(n_cfg_id) + "_" + STR(n_index);
	OFS ofs(file.c_str());
	if (!ofs) ERR(NAV1(file));

	ofs << iofmt("sci");
	ofs << setw(w_Real) << lnw << "\t" << setw(w_Real) << e_loc << endl;
	ofs << iofmt("def");
	return 1;
}

void Walker::reset_n_index()
{
	n_index = 0;
}