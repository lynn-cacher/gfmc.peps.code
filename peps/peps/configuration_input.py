import torch

def configuration_input():
    ## two-dimensional configuration array, where 0 is spin up and 1 is spin down
    ## Nx (number of rows) \times Ny (number of columns) array
    ## one should modify the configuration tensor,
    ## if he/she wants to calculate coefficient of other configuration
    '''
    configuration = torch.tensor([[0,0,0,0,0,0],
                                  [0,0,0,0,0,0],
                                  [0,0,0,0,0,0],
                                  [0,0,0,0,0,0],
                                  [0,0,0,0,0,0],
                                  [0,0,0,0,1,1]])
    '''
    '''
    configuration = torch.tensor([[0,1,0,1,0,1],
                                  [1,0,1,0,1,0],
                                  [0,1,0,1,0,1],
                                  [1,0,1,0,1,0],
                                  [0,1,0,1,0,1],
                                  [1,0,1,0,1,0]])
    '''
    configuration = torch.tensor([[0,1,0,1,0,1,0,1,0,1],
                                  [1,0,1,0,1,0,1,0,1,0],
                                  [0,1,0,1,0,1,0,1,0,1],
                                  [1,0,1,0,1,0,1,0,1,0],
                                  [0,1,0,1,0,1,0,1,0,1],
                                  [1,0,1,0,1,0,1,0,1,0],
                                  [0,1,0,1,0,1,0,1,0,1],
                                  [1,0,1,0,1,0,1,0,1,1],
                                  [0,1,0,1,0,1,0,1,0,1],
                                  [1,0,1,0,1,0,1,0,1,0]])
    return configuration
