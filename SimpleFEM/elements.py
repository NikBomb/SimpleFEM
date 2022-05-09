import numpy as np

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