import numpy as np
from .elements import Triangle

# vectors are horz
# kg	mm	ms	kN	GPa	kN-mm	7.83e-06	2.07e+02	15.65	9.806e-03
#       in  s   lbf psi



def nodeToMatrix(node1: int, dof1: int, node2:int, dof2:int ) -> list[int] :
    return [node1 * 2 + dof1, node2 * 2 + dof2]


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
             self.sigma = np.zeros([len(self.triangles), 3])
         for [triangle, idx] in zip(self.triangles, range(0, len(self.triangles))):
             u1 = self.sol[triangle.labels[0] * 2]
             v1 = self.sol[triangle.labels[0] * 2 + 1]
             u2 = self.sol[triangle.labels[1] * 2]
             v2 = self.sol[triangle.labels[1] * 2 + 1]
             u3 = self.sol[triangle.labels[2] * 2]
             v3 = self.sol[triangle.labels[2] * 2 + 1]
             epsilon = np.matmul(triangle.B, np.array([u1, v1, u2, v2, u3, v3]))
             self.sigma[idx, :] = np.matmul(triangle.elastic.E, epsilon).transpose()





