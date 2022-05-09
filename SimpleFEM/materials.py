import numpy as np

class Elastic:
    def __init__(self, Y, nu):
        self.Y = Y
        self.nu = nu
        self.E = Y / (1-nu**2) * np.array([[1, nu, 0],
                                           [nu, 1, 0],
                                           [0, 0, 0.5*(1-nu)]])