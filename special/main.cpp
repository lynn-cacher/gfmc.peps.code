#include "include.h"
#include "spcfcs.h"
#include "GFMC.h"

int main(int argc,char *argv[])
{
	using namespace std;
	MPI_Init(&argc, &argv);
	MyMpi mm;
	if (mm) cout << NAV(pwd()) << endl;
	use_mkl(mm);
	if (mm) io_init();
	clock_t t_program_bgn;
	if (mm) TIME_BGN("program", t_program_bgn);

	Pars p(mm);
	Lattice ltt(p.lx, p.ly, p.nsu);

	Heisenberg2d mdl(ltt, p);
	
	GFMC gfmc(mm, p, mdl);
	gfmc.walk();

	if (mm) TIME_END("program", t_program_bgn);
	MPI_Finalize();
	return 0;
}
