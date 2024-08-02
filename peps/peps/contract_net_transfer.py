import torch
import math
import numpy as np
import scipy.io as sio
from peps.CalGaugeU import InitGaugeU
from peps.TorNcon import ncon
import math

def contract_net_transfer(T, Dcut, norm_flag = False):
    ###    Ny: column
    ###  -|--|-  Nx: row
    ###  -|--|-     
    DT = T.size()
    logTnorm, TMT = transferM(T, Dcut)
    TM = TMT[0, 0]

    for ny in range(1, DT[1]):
        TMny = TMT[0, ny]
        TM = TM@TMny
        logTnorm = logTnorm + torch.log(TM.norm())
        TM = TM / TM.norm()

    Tfnorm = TM.trace()
    if norm_flag:
       Ntot = DT[0]*DT[1]
       Coef = torch.exp((logTnorm + Tfnorm.log())/Ntot)
       for j in range(DT[0]):
           for k in range(DT[1]):
               T[j, k] = T[j, k] / Coef
       return Tfnorm*torch.exp(logTnorm), Coef, T
    else:
       return Tfnorm*torch.exp(logTnorm)


def transferM(T, Dcut):
    ###    Ny: column
    ###  -|--|-  Nx: row
    ###  -|--|-     
    logTnorm = torch.zeros(1, dtype=T.dtype, device=T.device)
    DT = T.size()
    Dmin = min(DT[2]**DT[0], Dcut)
    TMT = torch.zeros(1, DT[1], Dmin, Dmin, dtype=T.dtype, device=T.device)

    for iter_step in range(DT[0]):
        if DT[0] == 1:
           for ny in range(DT[1]):
               TMT[0,ny] = torch.einsum('ijkj->ik', T[0,ny])
           break
        else:
           T, logCoef = OneContract(T, Dcut)
           logTnorm = logTnorm + logCoef
           DT = T.size()

    return logTnorm, TMT

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
            print('size, TrunError:', DT[0], DT[1], j, TrunError)

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
    
