import torch
import os
import time
from peps.ipeps.ipeps import *
# from peps.models import j1j2
import peps.config as cfg
import config as cfg1
from peps.contract_net import contract_net
from peps.contract_net_ed import contract_net_ed
from peps.contract_net_cut import contract_net_cut
from peps.contract_net_tebd import contract_net_tebd
from peps.contract_net_transfer import contract_net_transfer
from peps.TorNcon import ncon


class Amplitute:
    def __init__(self):
        self.config = cfg1.Config()
        self.init_A()
         
    def init_A(self):
        dataPath = './peps/peps/data/PEPS/J1_' + str(self.config.J) + 'J2_' + str(self.config.Jp) +\
                '_D' + str(round(self.config.bond_dim)) + '_chi' + str(self.config.chi)
        print(dataPath, flush=True)

        if not os.path.exists(dataPath+'_state.json'):
            print('the target local tensor does not exist, need to calculate it!', flush = True)
        else:
            print('the target local tensor already exists, just load it, its name is:', dataPath+'xxx', flush = True)
        state= read_ipeps(dataPath+'_state.json', aux_seq=[1,0,3,2], peps_args=cfg.peps_args, global_args=cfg.global_args)

        
        if self.config.tiling == "1SITE" or self.config.tiling == "1SITE-STRIPE":
            self.CNx, self.CNy = 1, 1
        elif self.config.tiling == "4SITE":
            self.CNx, self.CNy = 2, 2
        elif self.config.tiling == "STRIPE" or self.config.tiling == "2SITE":
            self.CNx, self.CNy = 1, 2
        self.A = state.sites
        self.DA = self.A[(0, 0)].size()
        

        # if os.path.exists(dataPath+'Nx'+str(self.config.lx)+'Ny'+str(self.config.ly)+'_normalized.tensor'):
        #     print('the near normalized local tensor already exists, just load it, its name is:', dataPath+'xxx', flush = True)
        #     print("111", flush=True)
        #     self.A = torch.load(dataPath+'Nx'+str(self.config.lx)+'Ny'+str(self.config.ly)+'_normalized.tensor').to(self.config.device)
        #     print("222", flush=True)
        # else:
        T = torch.zeros(self.config.lx, self.config.ly, self.DA[1]**2, self.DA[2]**2, self.DA[3]**2, self.DA[4]**2, dtype=self.A[(0,0)].dtype, device=self.A[(0,0)].device)
        for j in range(self.config.lx):
            for k in range(self.config.ly):
                ix, iy = k%self.CNy, j%self.CNx
                Aik = self.A[(ix, iy)]
                AAik = torch.einsum('aijkl,amnpq->imjnkplq', Aik, Aik)
                AAik = AAik.contiguous().reshape(self.DA[1]**2, self.DA[2]**2, self.DA[3]**2, self.DA[4]**2)
                T[j, k] = AAik.clone()
        if self.config.bond_dim**(2*self.config.lx) < 5000:
            Tnorm = contract_net_ed(T)
            Coef = torch.exp(torch.log(Tnorm)/(self.config.lx*self.config.ly))
            T = T/Coef
        else:
            Tnorm, Coef, T = contract_net(T, self.config.Dnorm, norm_flag = True)

        for ix in range(self.CNy):
            for iy in range(self.CNx):
                self.A[(ix,iy)] = self.A[(ix,iy)]/Coef.sqrt()
        if self.config.bond_dim**(2*self.config.lx) < 5000:
            Tnorm = contract_net_ed(T)
        else:      
            Tnorm = contract_net(T, self.config.Dnorm, norm_flag = False)
        print('===================================================', flush = True)
        print('test the normalization.', flush = True)
        print('attention: the normalization is only approximated when the truncation error is non-zero!!!', flush = True)
        print('===================================================', flush = True)
        print('Tnorm, norm_diff, Coef', Tnorm.item(), (1-Tnorm).item(), Coef.item(), flush = True)
        torch.save(self.A, dataPath+'Nx'+str(self.config.lx)+'Ny'+str(self.config.ly)+'_normalized.tensor')
        
        # self.index = torch.arange(self.config.ly - 1, -1, -1).to(torch.long)
        # print("self.index: ", self.index, flush = True)

    def get_amplitute_peps(self, cfg):
        lx = self.config.lx
        ly = self.config.ly
        cut_mode = self.config.cut_mode
        if cfg.shape[0] != lx * ly:
            print("amplitute.py: cfg size erro!!!", flush=True)
        cfg = torch.logical_not(cfg).reshape((lx, ly)).to(torch.int)
        ## 2. component
        # print(N, flush=True)
        t1 = time.time()
        amplitute = 0.
        T_amp = torch.zeros(lx, ly, self.DA[1], self.DA[2], self.DA[3], self.DA[4], dtype=self.A[(0,0)].dtype, device=self.A[(0,0)].device)
        for j in range(lx):
            for k in range(ly):
                ix, iy = k%self.CNy, j%self.CNx
                Aik = self.A[(ix, iy)]
                if self.config.tiling == "1SITE":
                    if (j+k)%2 == 0:
                        T_amp[j, k] = Aik[cfg[j, k]].clone()
                    else:
                        T_amp[j, k] = (((-1)**(cfg[j, k]+1)) * Aik[1-cfg[j, k]]).clone()
                elif self.config.tiling == "1SITE-STRIPE":
                    if k%2 == 0:
                        T_amp[j, k] = Aik[cfg[j, k]].clone()
                    else:
                        T_amp[j, k] = (((-1)**(cfg[j, k]+1)) * Aik[1-cfg[j, k]]).clone()
                elif self.config.tiling == "4SITE" or self.config.tiling == "STRIPE" or self.config.tiling == "2SITE":
                    T_amp[j, k] = Aik[cfg[j, k]].clone()
                else:
                    print("tiling form is not supported now", flush=True)
        # if self.config.bond_dim**(2*lx) < 5000:
        #     amplitute = contract_net_ed(T_amp)
        # else:
        #     amplitute = contract_net(T_amp, self.config.D, norm_flag = False)
        if 0 == cut_mode:
            amplitute = contract_net_cut(T_amp, self.config.D, norm_flag = False)
        elif 1 == cut_mode:
            amplitute = contract_net_tebd(T_amp, self.config.D)
        elif 2 == cut_mode:
            amplitute = contract_net_ed(T_amp)
        else:
            print("ampitute.py cut_mode is wrong!", flush=True)
            exit(0)
        t2 = time.time()
        #print("amplitute time: ", t2-t1, flush=True)   
        return amplitute
        

    def get_amplitute(self, cfg):
        cfg_in = torch.tensor(list(map(int, list(cfg))))
        res = self.get_amplitute_peps(cfg_in)
	    #print(res)
        return res
