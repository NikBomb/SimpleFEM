import numpy as np
import SimpleFEM as simFem
import matplotlib.pyplot as plt
plt.interactive(True)

#test from pdf 2D triangular Element from ->  Development of Membrane, Plate and Flat Shell Elements in Java by Kaushalkumar Kansara


def test_3():
    nodes = np.array([[0,12,24,36,48,
                       0,12,24,36,48,
                       6,6,0,6,12,
                       18,18,18,24,30,
                       30,30,36,42,42,
                       42,48],
                      [0,0,0,0,0,
                       12,12,12,12,12,
                       0,6,6,12,6,
                       0,6,12,6,0,
                       6,12,6,0,6,
                       12,6]])

    connectivity = np.array([[0,12,10],
                             [11,10,12],
                             [12,11,13],
                             [12,5,13],
                             [10,1,11],
                             [1,14,11],
                             [11,13,6],
                             [11,14,6],
                             [14,15,1],
                             [14,16,15],
                             [14,16,17],
                             [14,6,17],
                             [16,15,2],
                             [16,18,2],
                             [16,18,7],
                             [16,17,7],
                             [18,2,19],
                             [18,20,19],
                             [18,20,21],
                             [18,7,21],
                             [20,19,3],
                             [20,22,3],
                             [20,22,8],
                             [20,21,8],
                             [22,3,23],
                             [22,24,23],
                             [22,24,25],
                             [22,25,8],
                             [24,23,4],
                             [24,26,4],
                             [24,26,9],
                             [24,25,9],
                             ])

    alumininium = simFem.Elastic(30E6, 0.25)
    t = 1.0
    boundaryConditions = {
        "fixed": {0: "all", 5: "all", 12: "all"},
        "applied force": {4: [0, 6.67 * 1000], 9: [0, 6.67 * 1000], 26: [0, 26.67 * 1000]}
    }
    plate = simFem.Body(simFem.createTriangles(connectivity, nodes, alumininium, t), np.shape(nodes)[1])
    plate.assembleGlobalStiffness()
    plate.applyBcs(boundaryConditions)
    plate.solve()

    fig, ax = plt.subplots()
    ax.scatter(nodes[0, :] + plate.sol[0::2], nodes[1, :] + plate.sol[1::2])
    for i in range(np.shape(nodes)[1]):
        ax.annotate(i, (nodes[0, i], nodes[1, i]))
    plt.triplot(nodes[0, :], nodes[1, :], connectivity)
    plt.show()



