import torch
import os
from peps.ipeps.ipeps_c4v import *
from peps.models import j1j2
import argparse
import peps.config as cfg
import time
from peps.optim_j1j2_c4v import optim_j1j2_c4v
from peps.configuration_input import configuration_input
from peps.contract_net_cut import contract_net_cut
from peps.contract_net_tebd import contract_net_tebd
from peps.TorNcon import ncon

if __name__=='__main__':
   # parse command line args and build necessary configuration objects
   parser= cfg.get_args_parser()
   # additional model-dependent arguments
   parser.add_argument("--folder", default='../data/PEPS/',help="where to store results")
   parser.add_argument("--D", type=int, default=20, help="bond dimension of network contraction")
   parser.add_argument("--Dnorm", type=int, default=10, help="bond dimension of norm calculation")
   parser.add_argument("--j1", type=float, default=1., help="nearest-neighbour coupling")
   parser.add_argument("--j2", type=float, default=0., help="next nearest-neighbour coupling")
   parser.add_argument("--j3", type=float, default=0., help="next-to-next nearest-neighbour coupling")
   parser.add_argument("--hz_stag", type=float, default=0., help="staggered mag. field")
   parser.add_argument("--delta_zz", type=float, default=1., help="easy-axis (nearest-neighbour) anisotropy")
   parser.add_argument("--top_freq", type=int, default=-1, help="freuqency of transfer operator spectrum evaluation")
   parser.add_argument("--top_n", type=int, default=2, help="number of leading eigenvalues of transfer operator to compute")
   parser.add_argument("--force_cpu", action='store_true', help="evaluate energy on cpu")                  
   args, unknown_args = parser.parse_known_args()
   outputstatefile= args.out_prefix+"_state.json"

   if len(unknown_args)>0:
      print("args not recognized: "+str(unknown_args))
      raise Exception("Unknown command line arguments")
   
   curPath=os.getcwd()
   targetPath=curPath+'/data/PEPS'

   dataPath = './data/PEPS/J1_' + str(args.j1) + 'J2_' + str(args.j2) +\
              '_D' + str(args.bond_dim) + '_chi' + str(args.chi)

   if not os.path.exists(dataPath+'_state.json'):
      print('the target local tensor does not exist, need to calculate it!')
   else:
      print('the target local tensor already exists, just load it, its name is:', dataPath+'xxx')
   state= read_ipeps_c4v(dataPath+'_state.json', aux_seq=[1,0,3,2], peps_args=cfg.peps_args, global_args=cfg.global_args)

##############################################
########  peps generated from adpeps  ########
########     u s                      ########
########     |/                       ########
########  l--A--r  <=> a[s,l,u,r,d]   ########
########     |                        ########
########     d                        ########
##############################################
   A = state.sites[(0,0)].clone()
   DA = A.size()
   
   configuration = configuration_input()
   sys_size = configuration.size()
   Nx, Ny = sys_size[0], sys_size[1] ## Nx: row, Ny: column
   print('Row, Column:', Nx, Ny)

   if os.path.exists(dataPath+'Nx'+str(Nx)+'Ny'+str(Ny)+'_normalized.tensor'):
      print('the near normalized local tensor already exists, just load it, its name is:', dataPath+'xxx')
      A = torch.load(dataPath+'Nx'+str(Nx)+'Ny'+str(Ny)+'_normalized.tensor')

   ## 2. component
   t1 = time.time()
   T_amp = torch.zeros(Nx, Ny, DA[1], DA[2], DA[3], DA[4], dtype=A.dtype, device=A.device)
   for j in range(Nx):
       for k in range(Ny):
           if (j+k)%2 == 0:
              T_amp[j, k] = A[configuration[j, k]].clone()
           else:
              tempA_jk = ((-1)**(configuration[j, k]+1)) * A[1-configuration[j, k]]
              T_amp[j, k] = tempA_jk.clone()

   amplitute = contract_net_cut(T_amp, args.D, norm_flag = False)
   #amplitute2 = contract_net_tebd(T_amp, args.D)
   t2 = time.time()
   print('configuration:', configuration)
   print('amplitute_cut:', amplitute.item())
   #print('amplitute_tebd:', amplitute2.item())
   print('time:', t2-t1)
