import torch
import math
import numpy as np
import scipy.io as sio
from peps.CalGaugeU import InitGaugeU
from peps.TorNcon import ncon
import time

def contract_net_tebd(T, Dcut):
    ###    Ny: column
    ###  -|--|-  Nx: row
    ###  -|--|-     
    logTnorm = torch.zeros(1, dtype=T.dtype, device=T.device)
    tempT = T[0].clone()
    DT0 = T.size()

    for iter_step in range(DT0[0]-4):
        t1 = time.time()
        tempT, logCoef = OneContract(tempT, T[iter_step+1].clone(), Dcut)
        t2 = time.time()
        #print('time:', Dcut, logCoef, t2-t1)
        #input()
        logTnorm = logTnorm + logCoef

    TTM0 = ncon([tempT[0],T[DT0[0]-3,0],T[DT0[0]-2,0],T[DT0[0]-1,0]],[[-1,1,-5,2],[-2,2,-6,3],[-3,3,-7,4],[-4,4,-8,1]],[3,4,1,2])
    DTTM0 = TTM0.size()
    TTM0 = TTM0.contiguous().reshape(DTTM0[0]*DTTM0[1]*DTTM0[2]*DTTM0[3], -1)
    for k in range(1, DT0[1]):
        TTM = ncon([tempT[k],T[DT0[0]-3,k],T[DT0[0]-2,k],T[DT0[0]-1,k]],[[-1,1,-5,2],[-2,2,-6,3],[-3,3,-7,4],[-4,4,-8,1]],[3,4,1,2])
        TTM = TTM.contiguous().reshape(DTTM0[0]*DTTM0[1]*DTTM0[2]*DTTM0[3], -1)
        TTM0 = TTM0@TTM
    Tfnorm = torch.trace(TTM0)

    return Tfnorm*torch.exp(logTnorm)

def OneContract(T1, T2, Dcut):
    DT1 = T1.size()
    DT2 = T2.size()
    Dnew = [DT1[0], min(DT1[1]*DT2[1], Dcut), DT1[2], min(DT1[3]*DT2[3], Dcut), DT2[4]]
    Tnew = torch.zeros(Dnew, dtype=T1.dtype, device=T1.device)
    logCoef = torch.zeros(1, dtype=T1.dtype, device=T1.device)

    UL, UR = [], []
    ## 1. calculate gauge U
    for k in range(DT1[0]):
       Ta, Tb = T1[k%DT1[0]].clone(), T2[k%DT2[0]].clone()
       Tc, Td = T1[(k+1)%DT1[0]].clone(), T2[(k+1)%DT2[0]].clone()
       ULk, URk, Spectra, TrunError = InitGaugeU( Ta, Tb, Tc, Td, Dcut )
       UL.append(ULk)
       UR.append(URk)
       # print('step, TrunError:', k, TrunError)

    ## 2. renormalization
    for k in range(DT1[0]):
       Tak, Tbk = T1[k%DT1[0]].clone(), T2[k%DT2[0]].clone()
       ULk, URk = UL[(k-1)%DT1[0]], UR[k%DT2[0]]
       tempT = ncon([ULk,Tak,Tbk,URk],[[1,2,-3],[1,-4,5,6],[2,6,7,-9],[5,7,-8]],[1,7,2,6,5])
       logCoef = logCoef + torch.log(tempT.norm())
       Tnew[k] = tempT/tempT.norm()
    
    return Tnew, logCoef
    
