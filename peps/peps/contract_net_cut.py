import torch
import math
import numpy as np
import scipy.io as sio
from peps.CalGaugeU import InitGaugeU
from peps.TorNcon import ncon
import time

def contract_net_cut(T, Dcut, norm_flag = False):
    ###    Ny: column
    ###  -|--|-  Nx: row
    ###  -|--|-     
    logTnorm = torch.zeros(1, dtype=T.dtype, device=T.device)
    tempT = T.clone()
    DT = tempT.size()
    Ntot = DT[0]*DT[1]

    for iter_step in range(4*Ntot):
        if DT[0] != 1:
           t1 = time.time()
           tempT, logCoef = OneContract(tempT, Dcut)
           t2 = time.time()
           #print('time:', Dcut, logCoef, t2-t1)
           #input()
           logTnorm = logTnorm + logCoef

        ## clock-wise rotation
        tempT = tempT.permute([1,0,5,2,3,4])
        tempT2 = tempT.clone()
        ny = DT[0]//2 + DT[0]%2
        for g in range(ny):
            tempT[:, g] = tempT2[:, ny-g-1].clone()
        DT = tempT.size()
        if DT[0] == 1 and DT[1] == 1:
           break
    Tfnorm = torch.einsum('ijij', tempT[0,0])

    if norm_flag:
       DT = T.size()
       Coef = torch.exp((logTnorm + Tfnorm.log())/Ntot)
       for j in range(DT[0]):
           for k in range(DT[1]):
               T[j, k] = T[j, k] / Coef
       return Tfnorm*torch.exp(logTnorm), Coef, T
    else:
       return Tfnorm*torch.exp(logTnorm)

def OneContract(Told, Dcut):
    DT = Told.size()
    if DT[0]%2 == 0:
       Dnew = [DT[0]//2, DT[1], min(DT[2]**2, Dcut), DT[3], min(DT[4]**2, Dcut), DT[5]]
    else:
       Dnew = [(DT[0]+1)//2, DT[1], min(DT[2]**2, Dcut), DT[3], min(DT[4]**2, Dcut), DT[5]]
    Tnew = torch.zeros(Dnew, dtype=Told.dtype, device=Told.device)
    logCoef = torch.zeros(1, dtype=Told.dtype, device=Told.device)

    for j in range(DT[0]//2):
        UL, UR = [], []
        ## 1. calculate gauge U
        for k in range(DT[1]):
            Ta, Tb = Told[2*j, k%DT[1]].clone(), Told[2*j+1, k%DT[1]].clone()
            Tc, Td = Told[2*j, (k+1)%DT[1]].clone(), Told[2*j+1, (k+1)%DT[1]].clone()
            ULk, URk, Spectra, TrunError = InitGaugeU( Ta, Tb, Tc, Td, Dcut )
            UL.append(ULk)
            UR.append(URk)
            # print('size, TrunError:', DT[0], DT[1], j, TrunError)

        ## 2. renormalization
        for k in range(DT[1]):
            Tak, Tbk = Told[2*j, k%DT[1]].clone(), Told[2*j+1, k%DT[1]].clone()
            ULk, URk = UL[(k-1)%DT[1]], UR[k%DT[1]]
            tempT = ncon([ULk,Tak,Tbk,URk],[[1,2,-3],[1,-4,5,6],[2,6,7,-9],[5,7,-8]],[1,7,2,6,5])
            logCoef = logCoef + torch.log(tempT.norm())
            Tnew[j, k] = tempT/tempT.norm()
    
    if DT[0]%2 == 1:
       for k in range(DT[1]):
           logCoef = logCoef + torch.log(Told[DT[0]-1, k].norm())
           Tnew[DT[0]//2, k, :DT[2], :DT[3], :DT[4], :DT[5]] = Told[DT[0]-1, k] / Told[DT[0]-1, k].norm()
           #Tnew[DT[0]//2, k, :DT[2], :DT[3], :DT[4], :DT[5]] = Told[DT[0]-1, k]
    return Tnew, logCoef
    
