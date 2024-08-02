#include "Hamiltonian.h"

/////////////////////////////////////////////////////////////////////////////
/////////////////////////     Heisenberg2d        ///////////////////////////
////////////////////////////////////////////////////////////////////////////
void Heisenberg2d::get_all_bonds()
{
	if (ltt.n0 < 2 || ltt.n1 < 2)
	{
		// this is not a two-dimensional lattice
		ERR(NAV2(ltt.n0, ltt.n1));
	}

	//Defining the nearest-neighbors,0-left 1-right 2-top 3-buttom
	//Defining the next-nearest-neighbors,0-upper left 1-upper right 2-buttom left 3-buttom right
	Int i = 0, j = 0;
	for_Int(y, 0, ltt.n1)
	{
		for_Int(x, 0, ltt.n0)
		{
			Int i = ltt.index(LttVec(x, y));

			//nearest-neighbors
			j = ltt.left(i);
			if (i < j)
			{
				if (0 != x || isPBC) bonds.push_back(VEC<Int>{i, j});
			}

			j = ltt.right(i);
			if (i < j)
			{
				bonds.push_back(VEC<Int>{i, j});
			}
			j = ltt.top(i);
			if (i < j)
			{
				if (0 != y || isPBC) bonds.push_back(VEC<Int>{i, j});
			}
			j = ltt.buttom(i);
			if (i < j)
			{
				bonds.push_back(VEC<Int>{i, j});
			}

			//next-nearest-neighbors
			j = ltt.left(ltt.top(i));
			if (i < j)
			{
				if (0 != y || isPBC) next_bonds.push_back(VEC<Int>{i, j});		
			}
			j = ltt.right(ltt.top(i));
			if (i < j)
			{
				if (0 != y || isPBC) next_bonds.push_back(VEC<Int>{i, j});
			}
			j = ltt.left(ltt.buttom(i));
			if (i < j)
			{
				if (0 != x || isPBC) next_bonds.push_back(VEC<Int>{i, j});
			}
			j = ltt.right(ltt.buttom(i));
			if (i < j)
			{
				if ((ltt.n0 - 1) != x || isPBC) next_bonds.push_back(VEC<Int>{i, j});
			}
		}
	}
}

void Heisenberg2d::find_conn(const OccCfg &cfg, VEC<Real> &elem, VEC<VEC<Int> > &conn) const
{
	conn.clear();
	conn.resize(1);
	elem.resize(1);

	//computing interaction part Sz*Sz
	elem[0] = 0.;
	conn[0].resize(0);

	//Looks for possible spin flips
	double dX = 0., Z = 0., Zp = 0.;
	for_Int(i, 0, bonds.size())
	{
		dX = sgmz(cfg[bonds[i][0]]) * sgmz(cfg[bonds[i][1]]);
		Z += dX;
		elem.push_back(J - J * dX);
		conn.push_back(bonds[i]);
	}

	for_Int(j, 0, next_bonds.size())
	{
		dX = sgmz(cfg[next_bonds[j][0]]) * sgmz(cfg[next_bonds[j][1]]);
		Zp += dX;
		elem.push_back(Jp - Jp * dX);
		conn.push_back(next_bonds[j]);
	}

	elem[0] = J*Z + Jp*Zp;
}

