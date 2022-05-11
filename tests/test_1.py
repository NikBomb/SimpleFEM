import numpy as np
import SimpleFEM as simFem

#test from pdf 2D triangular Element from ->  Development of Membrane, Plate and Flat Shell Elements in Java by Kaushalkumar Kansara

def test_1():
    nodes = np.array([[3, 3, 0, 0], [0, 2, 2, 0]])
    connectivity = np.array([[0, 1, 3], [2, 3, 1]])
    boundaryConditions = {
        "fixed" : {0 : 1, 2 : "all", 3 : "all" },
        "applied force" : {1: [0, -1000]}
    }
    alumininium = simFem.Elastic(30E6, 0.25)
    t = 0.5
    plate = simFem.Body(simFem.createTriangles(connectivity, nodes, alumininium, t), np.shape(nodes)[1])
    plate.assembleGlobalStiffness()
    plate.applyBcs(boundaryConditions)
    plate.solve()
    assert abs(plate.sol[0] - 1.9e-5) < 1e-5
    assert abs(plate.sol[2] - 8.73e-6) < 1e-5
    assert abs(plate.sol[3] + 7.415e-5) < 1e-5
    assert np.linalg.norm(plate.sigma - np.array([[-93, -1136, -62], [93, 23, -297]])) < 1
