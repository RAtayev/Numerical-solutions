import numpy as np
import qr

A = np.array([[5, 6, 3],[-1, 0, 1],[1, 2, -1]])
res = qr.qr(A, 0.00001, 11)
print(np.diag(res))