import numpy as np
from .elements import Triangle

def createTriangles(connectivity, nodes, material, t):
    triangles = []
    for triangle in connectivity:
        x = np.empty(3)
        y = np.empty(3)
        labels = []
        for local, label in zip(range(3), triangle):
            x[local] = nodes[0][label]
            y[local] = nodes[1][label]
            labels.append(label)
        triangles.append(Triangle(x, y, labels, material, t))
    return triangles