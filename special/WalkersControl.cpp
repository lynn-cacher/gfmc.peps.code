#include "WalkersControl.h"

WalkersControl::WalkersControl(const MyMpi &mm_i, const Pars &p_i, Model &mdl_i, Int n_w_i)
	: mm(mm_i)
	, p(p_i)
	, mdl(mdl_i)
	, n_w(n_w_i)
	, ns(p_i.ns)
	, ifreconfig(p.ifreconfig)
{
	init_walkers_lt();
}

WalkersControl::~WalkersControl()
{
	for_Int(i, 0, walkers_lt.size())
	{
		delete walkers_lt[i];
	}
}

void WalkersControl::init_walkers_lt()
{
	for_Int(i, 0, n_w)
	{
		Walker* wal = new Walker(mm, p, mdl);
		wal->init_cfg(mm.id() * n_w + i);
		walkers_lt.push_back(wal);
	}
}

void WalkersControl::reset_lnw()
{
	for_Int(i, 0, walkers_lt.size())
	{
		walkers_lt[i]->set_lnw(0.0);
	}
}

void WalkersControl::propagate()
{
	for_Int(i, 0, walkers_lt.size())
	{
		walkers_lt[i]->propagate();
	}
}

void WalkersControl::set_new_cfg(OccCfg cfg_in, Int index)
{
	if (index >= walkers_lt.size())
	{
		ERR("void WalkersControl::set_new_cfg");
	}
	walkers_lt[index]->set_cfg(cfg_in);
}

VecReal WalkersControl::gather_lnw()
{
	VecReal vec(walkers_lt.size());
	for_Int(i, 0, walkers_lt.size())
	{
		vec[i] = walkers_lt[i]->get_lnw();
	}

	if (!ifreconfig) return vec;

	return  mm.Allgather(vec);
}

VecReal WalkersControl::gather_eloc()
{
	VecReal vec(walkers_lt.size());
	for_Int(i, 0, walkers_lt.size())
	{
		vec[i] = walkers_lt[i]->get_eloc();
	}

	if (!ifreconfig) return vec;

	return mm.Allgather(vec);
}

Vec<OccCfg> WalkersControl::gather_cfg()
{
	Vec<OccCfg> vec_curr_cfg(walkers_lt.size());
	for_Int(i, 0, walkers_lt.size())
	{
		vec_curr_cfg[i] = walkers_lt[i]->get_cfg();
	}

	if (!ifreconfig) return vec_curr_cfg;

	Vec<OccCfg> vec_total_cfg(mm.np() * vec_curr_cfg.size());

	for_Int(i, 0, vec_curr_cfg.size())
	{
		Vec<Char> vectemp = mm.Allgather(vec_curr_cfg[i]);
		Mat<Char> mat(mm.np(), p.ns, vectemp);

		for_Int(j, 0, mat.nrows())
		{
			Int id = j * vec_curr_cfg.size() + i;
			vec_total_cfg[id] = mat[j];
		}
	}

	return vec_total_cfg;
}

void WalkersControl::backup_cfg_pars(const Str &file_name_append) const
{
	for_Int(i, 0, walkers_lt.size())
	{
		walkers_lt[i]->write_cfg();
	}
}

void WalkersControl::reset_n_index()
{
	for_Int(i, 0, walkers_lt.size())
	{
		walkers_lt[i]->reset_n_index();
	}
}