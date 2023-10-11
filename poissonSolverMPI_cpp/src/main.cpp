#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <bits/stdc++.h>

using namespace std;

int main(int argc, char **argv)
{

  // bool period[3] = {false, false, false};
  int period[3] = {0, 0, 0};
  int domMPI[3] = {2, 2, 2};
  int nP[3] = {51, 51, 48};
  double pStart[3] = {0.0, 0.0, 0.0};
  double pEnd[3] = {1.0, 1.0, 1.0};
  double tol = 1e-4;
  double w = 1.75;

  struct gridType
  {
    double *p;
  };

  gridType x[3];

  string filename;
  int myperiod[3];
  bool bottom, top, east, west, north, south;
  int liney, planeyz, planexz, planexy;
  int i, j, k, cnt, nproc, ierr, id, flag, comm3d;
  MPI_Aint sizDP;
  int bs[3];
  int es[3];
  int bn[3];
  int en[3];
  int siz[3];
  int myDim[3];
  int myCoord[3];
  int src[3];
  int dest[3];
  int **tArr;
  int **bb;
  int **ee;
  MPI_Status STATUS;

  double ***u;
  double ***f;
  double t1, t2, t3, res, totRes[2], lhs, rhs, it1, it2, dx[3];

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Cart_create(MPI_COMM_WORLD, 3, domMPI, period, false, &comm3d);
  MPI_Cart_get(comm3d, 3, myDim, myperiod, myCoord);
  MPI_Cart_shift(comm3d, 0, 1, &src[0], &dest[0]);
  MPI_Cart_shift(comm3d, 1, 1, &src[1], &dest[1]);
  MPI_Cart_shift(comm3d, 2, 1, &src[2], &dest[2]);

  tArr = new int *[*max_element(myDim, myDim + 2)];
  for (i = 0; i < *max_element(myDim, myDim + 2); i++)
  {
    tArr[i] = new int[3];
  }

  bb = new int *[3];
  ee = new int *[3];
  for (i = 0; i < 3; i++)
  {
    bb[i] = new int[nproc];
    ee[i] = new int[nproc];
  }

  east = (myCoord[0] == myDim[0] - 1);
  west = (myCoord[0] == 0);
  north = (myCoord[1] == myDim[1] - 1);
  south = (myCoord[1] == 0);
  top = (myCoord[2] == myDim[2] - 1);
  bottom = (myCoord[2] == 0);

  for (k = 0; k < 3; k++)
  {
    it1 = nP[k] / myDim[k];
    it2 = myDim[k] - (nP[k] % myDim[k]);
    tArr[0][0] = 1;
    for (i = 0; i < myDim[k]; i++)
    {
      if (i == it2)
      {
        it1 = it1 + 1;
      }
      tArr[i][1] = tArr[i][0] + it1 - 1;
      tArr[i][2] = tArr[i][1] - tArr[i][0] + 1;
      if (i == myDim[k] - 1)
      {
        break;
      }
      tArr[i + 1][0] = tArr[i][1] + 1;
    }

    for (i = 0; i < myDim[k]; i++)
    {
      if (i == myCoord[k])
      {
        bs[k] = tArr[i][0];
        es[k] = tArr[i][1];
        siz[k] = tArr[i][2];
      }
    }
    bn[k] = bs[k];
    en[k] = es[k];
    if ((west && k == 0) || (south && k == 1) || (bottom && k == 2))
    {
      bn[k] = bs[k] + 1;
    }
    if ((east && k == 0) || (north && k == 1) || (top && k == 2))
    {
      en[k] = es[k] - 1;
    }
  }

  MPI_Type_vector(siz[1], 1, nP[0], MPI_DOUBLE_PRECISION, &liney);
  MPI_Type_commit(&liney);
  MPI_Type_extent(MPI_DOUBLE_PRECISION, &sizDP);
  MPI_Type_hvector(siz[2], 1, nP[0] * nP[1] * sizDP, liney, &planeyz);
  MPI_Type_vector(siz[2], siz[0], nP[0] * nP[1], MPI_DOUBLE_PRECISION, &planexz);
  MPI_Type_vector(siz[1], siz[0], nP[0], MPI_DOUBLE_PRECISION, &planexy);
  MPI_Type_commit(&planeyz);
  MPI_Type_commit(&planexz);
  MPI_Type_commit(&planexy);
  MPI_Barrier(MPI_COMM_WORLD);

  x[0].p = new double[nP[0]];
  x[1].p = new double[nP[1]];
  x[2].p = new double[nP[2]];

  x[0].p[0] = pStart[0];
  dx[0] = (pEnd[0] - pStart[0]) / (nP[0] - 1);
  for (i = 1; i < nP[0]; i++)
  {
    x[0].p[i] = x[0].p[i - 1] + dx[0];
  }

  x[1].p[0] = pStart[1];
  dx[1] = (pEnd[1] - pStart[1]) / (nP[1] - 1);
  for (i = 1; i < nP[1]; i++)
  {
    x[1].p[i] = x[1].p[i - 1] + dx[1];
  }

  x[2].p[0] = pStart[2];
  dx[2] = (pEnd[2] - pStart[2]) / (nP[2] - 1);
  for (i = 1; i < nP[2]; i++)
  {
    x[2].p[i] = x[2].p[i - 1] + dx[2];
  }

  u = new double **[nP[0]];
  f = new double **[nP[0]];
  for (i = 0; i < nP[0]; i++)
  {
    u[i] = new double *[nP[1]];
    f[i] = new double *[nP[1]];
    for (j = 0; j < nP[1]; j++)
    {
      u[i][j] = new double[nP[2]];
      f[i][j] = new double[nP[2]];
    }
  }

  for (i = bs[0]; i <= es[0]; i++)
  {
    for (j = bs[1]; j <= es[1]; j++)
    {
      for (k = bs[2]; k <= es[2]; k++)
      {
        f[i][j][k] = -5.0;
        u[i][j][k] = 0.0;
      }
    }
  }

  t1 = 1.0 / pow(dx[0], 2.0);
  t2 = 1.0 / pow(dx[1], 2.0);
  t3 = 1.0 / pow(dx[2], 2.0);
  flag = 0;
  cnt = 0;

  while (true)
  {
    cnt = cnt + 1;
    for (k = bn[2]; k <= en[2]; k++)
    {
      for (j = bn[1]; j <= en[1]; j++)
      {
        for (i = bn[0]; i <= en[0]; i++)
        {
          u[i][j][k] = (1 - w) * u[i][j][k] + w * ((u[i - 1][j][k] + u[i + 1][j][k]) * t1 + (u[i][j - 1][k] + u[i][j + 1][k]) * t2 + (u[i][j][k - 1] + u[i][j][k + 1]) * t3 - f[i][j][k]) / (2 * (t1 + t2 + t3));
        }
      }
    }

    double sendBuf[2];
    double recvBuf[2];

    MPI_Sendrecv(&u[en[0]][bn[1]][bn[2]], 1, planeyz, dest[0], 50, &u[bn[0] - 1][bn[1]][bn[2]], 1, planeyz, src[0], 50, MPI_COMM_WORLD, &STATUS);
    MPI_Sendrecv(&u[bn[0]][bn[1]][bn[2]], 1, planeyz, src[0], 50, &u[en[0] + 1][bn[1]][bn[2]], 1, planeyz, dest[0], 50, MPI_COMM_WORLD, &STATUS);
    MPI_Sendrecv(&u[bn[0]][en[1]][bn[2]], 1, planexz, dest[1], 50, &u[bn[0]][bn[1] - 1][bn[2]], 1, planexz, src[1], 50, MPI_COMM_WORLD, &STATUS);
    MPI_Sendrecv(&u[bn[0]][bn[1]][bn[2]], 1, planexz, src[1], 50, &u[bn[0]][en[1] + 1][bn[2]], 1, planexz, dest[1], 50, MPI_COMM_WORLD, &STATUS);
    MPI_Sendrecv(&u[bn[0]][bn[1]][en[2]], 1, planexy, dest[2], 50, &u[bn[0]][bn[1]][bn[2] - 1], 1, planexy, src[2], 50, MPI_COMM_WORLD, &STATUS);
    MPI_Sendrecv(&u[bn[0]][bn[1]][bn[2]], 1, planexy, src[2], 50, &u[bn[0]][bn[1]][en[2] + 1], 1, planexy, dest[2], 50, MPI_COMM_WORLD, &STATUS);

    res = 0.0;
    totRes[1] = totRes[0];
    for (k = bn[2]; k <= en[2]; k++)
    {
      for (j = bn[1]; j <= en[1]; j++)
      {
        for (i = bn[0]; i <= en[0]; i++)
        {
          lhs = (u[i - 1][j][k] + u[i + 1][j][k]) * t1 +
                (u[i][j - 1][k] + u[i][j + 1][k]) * t2 +
                (u[i][j][k - 1] + u[i][j][k + 1]) * t3 -
                2 * u[i][j][k] * (t1 + t2 + t3);
          rhs = -5.0;
          res = res + (lhs - rhs);
        }
      }
    }

    MPI_Reduce(&res, &totRes[0], 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD);
    totRes[0] = abs(totRes[0] / (nP[0] * nP[1] * nP[2]));

    if (id == 0)
    {
      if (totRes[0] < tol && totRes[0] < totRes[1])
      {
        cout << "Poisson converged at " << cnt << " iteration : Residual " << totRes[0] << endl;
        flag = 1;
      }
      if (cnt > 50 && totRes[0] > totRes[1])
      {
        cout << "Poisson diverged at " << cnt << " iteration : Residual " << totRes[0] << endl;
        break;
      }
    }

    MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (flag == 1)
    {
      break;
    }
  }

  filename = "../output/3D" + to_string(id + 1) + ".dat";
  ofstream file(filename);

  file << "Title = \"Poisson Solution\"" << endl;
  file << "Variables = \"x\",\"y\",\"z\",\"u\"" << endl;
  file << "Zone k=" << es[2] - bs[2] + 1 << ",j=" << es[1] - bs[1] + 1 << ",i=" << es[0] - bs[0] + 1 << ", DATAPACKING=\"POINT\"" << endl;

  for (k = bs[2]; k <= es[2]; k++)
  {
    for (j = bs[1]; j <= es[1]; j++)
    {
      for (i = bs[0]; i <= es[0]; i++)
      {
        file << x[0].p[i] << " " << x[1].p[j] << " " << x[2].p[k] << " " << u[i][j][k] << endl;
      }
    }
  }

  file.close();

  MPI_Gather(&bs[1], 3, MPI_INT, bb[0], 3, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&es[1], 3, MPI_INT, ee[0], 3, MPI_INT, 0, MPI_COMM_WORLD);

  if (id == 0)
  {
    ofstream domInfo("domInfo.txt");
    domInfo << nP[0] << " " << nP[1] << " " << nP[2] << endl;
    domInfo << nproc << endl;
    for (i = 0; i < nproc; i++)
    {
      domInfo << bb[0][i] << " " << ee[0][i] << " " << bb[1][i] << " " << ee[1][i] << " " << bb[2][i] << " " << ee[2][i] << endl;
    }
    domInfo.close();
  }

  MPI_Type_free(&planeyz);
  MPI_Type_free(&planexz);
  MPI_Type_free(&planexy);
  MPI_Type_free(&liney);
  MPI_Finalize();

  if (id == 0)
  {
    cout << "Program executed successfully" << endl;
  }

  return 0;
}
