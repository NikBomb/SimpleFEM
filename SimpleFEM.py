import numpy as np
# vectors are horz
# kg	mm	ms	kN	GPa	kN-mm	7.83e-06	2.07e+02	15.65	9.806e-03
#       in  s   lbf psi



def nodeToMatrix(node1: int, dof1: int, node2:int, dof2:int ) -> list[int] :
    return [node1 * 2 + dof1, node2 * 2 + dof2]


class Triangle:
    def __init__(self, x, y, labels, elastic, t):
        self.x = x
        self.y = y
        self.elastic = elastic
        self.labels = labels
        self.t = t
        x1 = x[0]
        x2 = x[1]
        x3 = x[2]
        y1 = y[0]
        y2 = y[1]
        y3 = y[2]
        self.detJ = (x1 - x3)*(y2-y3) - (y1-y3)*(x2 - x3)
        self.A = 0.5*abs(self.detJ)
        self.B = 1 / self.detJ * np.array([[y2 - y3, 0, y3 - y1, 0, y1 - y2, 0],
                                              [0, x3 - x2, 0, x1 - x3, 0, x2 - x1],
                                              [x3 - x2, y2 - y3, x1 - x3, y3 - y1, x2- x1, y1 - y2] ])
        self.K = self.t * self.A * np.matmul(np.matmul(self.B.transpose(),  elastic.E), self.B)


class Body:
     def __init__(self, triangles : list[Triangle], nnodes:int):
         self.triangles = triangles
         ndof = 2
         self.K = np.zeros([nnodes * ndof , nnodes * ndof])
         self.f = np.zeros([nnodes * ndof])
         self.sol = np.zeros([nnodes * ndof])
         self.solInd = range(0, nnodes*ndof)

     def assembleGlobalStiffness(self):
          for triangle in self.triangles:
              for row in range(0,6):
                  for column in range(0,6):
                     globalRow = triangle.labels[int(row / 2)] * 2 + int(row % 2)
                     globalCol = triangle.labels[int(column / 2)] * 2 + int(column % 2)
                     self.K[globalRow, globalCol] = self.K[globalRow, globalCol] + triangle.K[row, column]


     def applyBcs(self, bcs : dict):
         toDeleteDof = []
         for type, dictApplied in bcs.items():
            if type == "fixed":
                for node, dofs in dictApplied.items():
                    if (dofs == "all"):
                        for dof in range(0, 2):
                            row = nodeToMatrix(node, dof, node, dof)[0]
                            toDeleteDof.append(row)
                    else:
                          row = nodeToMatrix(node, dofs, node, dofs)[0]
                          toDeleteDof.append(row)
            elif type == "applied force":
                for node, force in dictApplied.items():
                        rows = [ nodeToMatrix(node, 0, node, 0)[0], nodeToMatrix(node, 1, node, 1)[1]]
                        self.f[rows[0]] = force[0]
                        self.f[rows[1]] = force[1]
            else:
                print ("Boundary condition of type {} not implemented".format(type))
         print("To Delete: {}".format(toDeleteDof) )
         self.solInd = [x for x in self.solInd if x not in toDeleteDof]
         self.KR = self.K[np.ix_(self.solInd, self.solInd)]
         self.fR = self.f[self.solInd]

     def solve(self):
         sol = np.linalg.solve(self.KR, self.fR)
         for dof, idx in zip(sol, self.solInd):
             self.sol[idx] = dof





class Elastic:
    def __init__(self, Y, nu):
        self.Y = Y
        self.nu = nu
        self.E = Y / (1-nu**2) * np.array([[1, nu, 0],
                                           [nu, 1, 0],
                                           [0, 0, 0.5*(1-nu)]])


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