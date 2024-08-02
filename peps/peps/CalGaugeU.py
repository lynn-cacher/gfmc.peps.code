import torch
from peps.TorNcon import ncon
###  Ta&Tb:     2               U:
###             |                    1 ---|
###         1 ----- 3                     |---3
###             |                    2 ---|
###             4

def InitGaugeU( Ta, Tb, Tc, Td, Dcut ):
    DTa, DTb, DTc, DTd = Ta.size(), Tb.size(), Tc.size(), Td.size()
    ## 1. QR
    TTem = ncon([Ta, Tb, Ta, Tb], ([1,2,-3,4], [5,4,-6,7], [1,2,-8,9], [5,9,-10,7]), [1,2,5,7,4,9])
    DTem = TTem.size()
    TTem = TTem.contiguous().reshape(DTem[0]*DTem[1], DTem[2]*DTem[3])
    TTem = (TTem + TTem.t())/2
    ER, UR = torch.linalg.eigh(TTem)
    R = torch.diag(ER.abs().sqrt()) @ UR.t()

    ## 2. LQ 
    TTem = ncon([Tc, Td, Tc, Td], ([-1,2,3,4], [-5,4,6,7], [-8,2,3,9], [-10,9,6,7]), [2,3,6,7,4,9])
    DTem = TTem.size()
    TTem = TTem.contiguous().reshape(DTem[0]*DTem[1], DTem[2]*DTem[3])
    TTem = (TTem + TTem.t())/2
    EL, UL = torch.linalg.eigh(TTem)
    L = UL @ torch.diag(EL.abs().sqrt())

    ## 3. R@L = U@torch.diag(S)@Vh (linalg.svd)
    U, S, Vh = torch.linalg.svd(R@L)
    sorted, indices = torch.sort(S, descending = True)
    Spectra = S[indices]/S[indices[0]]
    U = U[:, indices]
    Vh = Vh[indices, :]
    S = S[indices]
    U = U[:, :Dcut]
    Vh = Vh[:Dcut, :]
    S = S[:Dcut]
    Sinv = 1/(S + 1e-12)
    Sinv_sqrt = Sinv.sqrt()

    ## 4. Calculate P & Q
    UR = L @ Vh.t() @ torch.diag(Sinv_sqrt)
    UL = R.t() @ U @ torch.diag(Sinv_sqrt)

    TrunError = torch.sum(Spectra[Dcut:])/torch.sum(Spectra)
    UL = UL.contiguous().reshape(DTc[0], DTd[0], -1)
    UR = UR.contiguous().reshape(DTa[2], DTb[2], -1)
    
    return UL, UR, Spectra, TrunError
