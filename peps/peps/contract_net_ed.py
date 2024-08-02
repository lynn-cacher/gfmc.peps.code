import torch
import math
import numpy as np
import scipy.io as sio
from peps.TorNcon import ncon
import math

def contract_net_ed(T):
    ###    Ny: column
    ###  -|--|-  Nx: row
    ###  -|--|-     
    DT = T.size()
    TM = transferM(T, 0)

    for ny in range(1, DT[1]):
        TMny = transferM(T, ny)
        TM = TM@TMny

    Tnorm = TM.trace()
    return Tnorm

def transferM(T, ny):
    ## transfer matrix of ny-th column, ny = 0,1,...,Ny-1
    DT = T.size()
    if DT[0] == 1:
       TMny = torch.einsum('ijkj->ik', T[0,ny])
    else:
       if DT[0] > 2:
          TMny = ncon([T[0,ny],T[1,ny]],[[-1,-3,-4,7],[-2,7,-5,-6]],[7])
          DTMny = TMny.size()
          TMny = TMny.contiguous().reshape(DTMny[0]*DTMny[1], DTMny[2], DTMny[3]*DTMny[4], DTMny[5])

          for nx in range(DT[0]-3):
              TMny = ncon([TMny,T[nx+2,ny]],[[-1,-3,-4,7],[-2,7,-5,-6]],[7])
              DTMny = TMny.size()
              TMny = TMny.contiguous().reshape(DTMny[0]*DTMny[1], DTMny[2], DTMny[3]*DTMny[4], DTMny[5])
       if DT[0] == 2:
          TMny = T[0,ny]
       TMny = ncon([TMny, T[DT[0]-1,ny]],[[-1,5,-3,6],[-2,6,-4,5]],[5,6])
       DTMny = TMny.size()
       TMny = TMny.contiguous().reshape(DTMny[0]*DTMny[1], DTMny[2]*DTMny[3])
    
    return TMny
