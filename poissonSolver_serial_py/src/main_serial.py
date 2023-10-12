import numpy as np


class gridType:
    def __init__(self):
        self.p = None


nP = [51, 51, 48]
pStart = [0.0, 0.0, 0.0]
pEnd = [1.0, 1.0, 1.0]
tol = 1e-4
w = 1.75

x = [gridType() for _ in range(3)]
filename = ""

i, j, k, cnt = 0, 0, 0, 0
u = np.zeros((nP[0], nP[1], nP[2]))
f = np.zeros((nP[0], nP[1], nP[2]))
t1, t2, t3, res, totRes = 0.0, 0.0, 0.0, 0.0, [0.0, 0.0]
lhs, rhs, it1, it2, dx = 0.0, 0.0, 0.0, 0.0, [0.0, 0.0, 0.0]

x[0].p = np.zeros(nP[0])
x[1].p = np.zeros(nP[1])
x[2].p = np.zeros(nP[2])

x[0].p[0] = pStart[0]
dx[0] = (pEnd[0] - pStart[0]) / (nP[0] - 1)
for i in range(nP[0] - 1):
    x[0].p[i + 1] = x[0].p[i] + dx[0]

x[1].p[0] = pStart[1]
dx[1] = (pEnd[1] - pStart[1]) / (nP[1] - 1)
for i in range(nP[1] - 1):
    x[1].p[i + 1] = x[1].p[i] + dx[1]

x[2].p[0] = pStart[2]
dx[2] = (pEnd[2] - pStart[2]) / (nP[2] - 1)
for i in range(nP[2] - 1):
    x[2].p[i + 1] = x[2].p[i] + dx[2]

for k in range(nP[2]):
    for j in range(nP[1]):
        for i in range(nP[0]):
            f[i][j][k] = -5.0
            u[i][j][k] = 0.0

t1 = 1.0 / (dx[0] * dx[0])
t2 = 1.0 / (dx[1] * dx[1])
t3 = 1.0 / (dx[2] * dx[2])
cnt = 0

while True:
    cnt += 1

    for k in range(1, nP[2] - 1):
        for j in range(1, nP[1] - 1):
            for i in range(1, nP[0] - 1):
                u[i][j][k] = (1 - w) * u[i][j][k] + w * ((u[i - 1][j][k] + u[i + 1][j][k]) * t1 + (u[i][j - 1][k] +
                                                                                                   u[i][j + 1][k]) * t2 + (u[i][j][k - 1] + u[i][j][k + 1]) * t3 - f[i][j][k]) / (2 * (t1 + t2 + t3))

    res = 0.0
    totRes[1] = totRes[0]

    for k in range(1, nP[2] - 1):
        for j in range(1, nP[1] - 1):
            for i in range(1, nP[0] - 1):
                lhs = (u[i - 1][j][k] + u[i + 1][j][k]) * t1 + (u[i][j - 1][k] + u[i][j + 1][k]) * \
                    t2 + (u[i][j][k - 1] + u[i][j][k + 1]) * \
                    t3 - 2 * u[i][j][k] * (t1 + t2 + t3)
                rhs = -5.0
                res += (lhs - rhs)

    totRes[0] = abs(res / (nP[0] * nP[1] * nP[2]))
    if totRes[0] < tol and totRes[0] < totRes[1]:
        print("Poisson converged at", cnt, "iteration : Residual", totRes[0])
        break

    if cnt > 50 and totRes[0] > totRes[1]:
        print("Poisson diverged at", cnt, "iteration : Residual", totRes[0])
        break

filename = "../output/3D.dat"
file = open(filename, "w")
file.write("Title = \"Poisson Solution\"\n")
file.write("Variables = \"x\",\"y\",\"z\",\"u\"\n")
file.write("Zone k=" + str(nP[2]) + ",j=" + str(nP[1]) +
           ",i=" + str(nP[0]) + ", DATAPACKING=\"POINT\"\n")

for k in range(nP[2]):
    for j in range(nP[1]):
        for i in range(nP[0]):
            file.write(str(x[0].p[i]) + " " + str(x[1].p[j]) +
                       " " + str(x[2].p[k]) + " " + str(u[i][j][k]) + "\n")

file.close()
print("Program executed successfully")
del x[0].p
del x[1].p
del x[2].p
