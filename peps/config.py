import attr
import torch
import os
from os.path import dirname, abspath
# os.environ["CUDA_VISIBLE_DEVICES"] = "0"

def _validate_directory(obj, attribute, value):
    """Validates value is a directory."""
    del obj
    if value and not os.path.isdir(value):
        raise ValueError(f'{attribute.name} is not a directory')


class Config:
    # device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    device = torch.device('cpu')
    path = dirname(abspath(__file__))
    cpu: str = 'cpu'
    dtype = torch.float64

    Hamiltonian: str = "nH"            # ["H", "nH"]   no sH
    lx: int = 4
    ly: int = 4
    size = lx * ly
    n_occ: int = size // 2
    J: torch.float64 = 1.0
    Jp: torch.float64 = 0.5

    beta = 0.02                  #the evolution of time for the Continuous-Time limit process
    nP = 16                      #the number of reconfig processes. beta * nP should equal to 20~40.
    nStep = 100                  #the number of sample step

    ifreconfig = False          #if false, then Nw must equal to the number of processes

    Nw = 1                      #the number of walkers
    ntotalSample = Nw * nStep   #the number of total samples
    nbin = 20
    #gamma = 0.0
    gamma = 1
    max_relnpsi = 0.0

    # Nw_per_prog = (Nw - 1) / np + 1
    # Nw = Nw_per_prog * np

    # guiding wave function model, 0: ffnn,   1: PEPS
    psi_g_mode = 1

    # PEPS
    bond_dim = 4
    chi = 60
    D = 20
    Dnorm = 20
    # 0: hotrg, 1: tebd, 2: not cutting
    cut_mode = 1

    # 1SITE, STRIPE
    tiling = "4SITE"

    # ffnn
    ifcanonical: bool = False

