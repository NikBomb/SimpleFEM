import numpy as np
import SimpleFEM as simFem
import matplotlib.pyplot as plt
plt.interactive(True)

#test from pdf 2D triangular Element from ->  Development of Membrane, Plate and Flat Shell Elements in Java by Kaushalkumar Kansara


def test_2():
    nodes = np.array([[0,12,24,36,48,0,12,24,36,48],
                      [0,0,0,0,0,12,12,12,12,12]])
    connectivity = np.array([[0,6,5],
                            [0,1,6],
                            [1,7,6],
                            [1,2,7],
                            [2,8,7],
                            [2,3,8],
                            [3,9,8],
                            [3,4,9]])

    alumininium = simFem.Elastic(30E6, 0.25)
    t = 1.0
    boundaryConditions = {
        "fixed": {0: "all", 5: "all"},
        "applied force": {4: [0, 20000], 9: [0, 20000] }
    }
    plate = simFem.Body(simFem.createTriangles(connectivity, nodes, alumininium, t), np.shape(nodes)[1])
    plate.assembleGlobalStiffness()
    plate.applyBcs(boundaryConditions)
    plate.solve()
    assert abs(plate.sol[18] + 0.014159) < 1e-5
    assert abs(plate.sol[19] - 0.090347) < 1e-5
    assert abs(plate.sol[14] + 0.010825) < 1e-5
    assert abs(plate.sol[15] - 0.030403) < 1e-5
    plt.triplot(nodes[0,:], nodes[1,:], connectivity)
    plt.show()



