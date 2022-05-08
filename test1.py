import numpy as np
import PlateFEM as pfem

#test from pdf 2D triangular Elemnt pdf

if __name__ == "__main__":
    nodes = np.array([[3, 3, 0, 0], [0, 2, 2, 0]])
    connectivity = np.array([[0, 1, 3], [2, 3, 1]])
    boundaryConditions = {
        "fixed" : {0 : 1, 2 : "all", 3 : "all" },
        "applied force" : {1: [0, -1000]}
    }
    alumininium = pfem.Elastic(30E6, 0.25)
    t = 0.5
    plate = pfem.Body(pfem.createTriangles(connectivity, nodes, alumininium, t), np.shape(nodes)[1])
    plate.assembleGlobalStiffness()
    plate.applyBcs(boundaryConditions)
    plate.solve()
    assert abs(plate.sol[0] - 1.9e-5) < 1e-5
    assert abs(plate.sol[2] - 8.73e-6) < 1e-5
    assert abs(plate.sol[3] + 7.415e-5) < 1e-5
