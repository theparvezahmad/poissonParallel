#include <iostream>
#include <fstream>
#include <cmath>
#include <omp.h>

using namespace std;

struct gridType
{
  double *p{nullptr};
};

int main()
{
  const int nP[3] = {101, 101, 101};
  const double pStart[3] = {0.0, 0.0, 0.0};
  const double pEnd[3] = {1.0, 1.0, 1.0};
  const double tol = 1e-8;
  const double w = 1.75;

  gridType x[3];
  string filename;
  int i, j, k, cnt;
  bool isConverged = false;
  bool isDiverged = false;

  // double u[nP[0]][nP[1]][nP[2]];
  double f;
  double t1, t2, t3, res, totRes[2], lhs, rhs, it1, it2, dx[3];

  // Dynamically allocate the 3D array
  double ***u = new double **[nP[0]];
  for (int i = 0; i < nP[0]; ++i)
  {
    u[i] = new double *[nP[1]];
    for (int j = 0; j < nP[1]; ++j)
    {
      u[i][j] = new double[nP[2]];
    }
  }

  x[0].p = new double[nP[0]];
  x[1].p = new double[nP[1]];
  x[2].p = new double[nP[2]];

  x[0].p[0] = pStart[0];
  dx[0] = (pEnd[0] - pStart[0]) / (nP[0] - 1);
  for (i = 0; i < nP[0] - 1; i++)
  {
    x[0].p[i + 1] = x[0].p[i] + dx[0];
  }

  x[1].p[0] = pStart[1];
  dx[1] = (pEnd[1] - pStart[1]) / (nP[1] - 1);
  for (i = 0; i < nP[1] - 1; i++)
  {
    x[1].p[i + 1] = x[1].p[i] + dx[1];
  }

  x[2].p[0] = pStart[2];
  dx[2] = (pEnd[2] - pStart[2]) / (nP[2] - 1);
  for (i = 0; i < nP[2] - 1; i++)
  {
    x[2].p[i + 1] = x[2].p[i] + dx[2];
  }

  for (k = 0; k < nP[2]; k++)
  {
    for (j = 0; j < nP[1]; j++)
    {
      for (i = 0; i < nP[0]; i++)
      {
        // f[i][j][k] = -5.0;
        f = -5.0;
        u[i][j][k] = 0.0;
      }
    }
  }

  t1 = 1.0 / (dx[0] * dx[0]);
  t2 = 1.0 / (dx[1] * dx[1]);
  t3 = 1.0 / (dx[2] * dx[2]);
  cnt = 0;

  // #pragma omp parallel
  // {
  while (!isConverged && !isDiverged)
  {
    cnt++;
#pragma omp for
    for (k = 1; k < nP[2] - 1; k++)
    {
      for (j = 1; j < nP[1] - 1; j++)
      {
        for (i = 1; i < nP[0] - 1; i++)
        {
          u[i][j][k] = (1 - w) * u[i][j][k] + w * ((u[i - 1][j][k] + u[i + 1][j][k]) * t1 + (u[i][j - 1][k] + u[i][j + 1][k]) * t2 + (u[i][j][k - 1] + u[i][j][k + 1]) * t3 - f) / (2 * (t1 + t2 + t3));
        }
      }
    }

    // #pragma omp barrier
    res = 0.0;
    totRes[1] = totRes[0];
    // #pragma omp for reduction(+ : res)
    for (k = 1; k < nP[2] - 1; k++)
    {
      for (j = 1; j < nP[1] - 1; j++)
      {
        for (i = 1; i < nP[0] - 1; i++)
        {
          lhs = (u[i - 1][j][k] + u[i + 1][j][k]) * t1 +
                (u[i][j - 1][k] + u[i][j + 1][k]) * t2 +
                (u[i][j][k - 1] + u[i][j][k + 1]) * t3 -
                2 * u[i][j][k] * (t1 + t2 + t3);
          rhs = -5.0;
          res += (lhs - rhs);
        }
      }
    }

    totRes[0] = abs(res) / (nP[0] * nP[1] * nP[2]);

    // #pragma omp single
    // {
    if (totRes[0] < tol && totRes[0] < totRes[1])
    {
      cout << "Poisson converged at " << cnt << " iteration : Residual " << totRes[0] << endl;
      isConverged = true;
    }

    if (cnt > 50 && totRes[0] > totRes[1])
    {
      cout << "Poisson diverged at " << cnt << " iteration : Residual " << totRes[0] << endl;
      isDiverged = true;
    }
    // }
  }

  filename = "../output/3D.dat";
  ofstream file(filename);
  file << "Title = \"Poisson Solution\"" << endl;
  file << "Variables = \"x\",\"y\",\"z\",\"u\"" << endl;
  file << "Zone k=" << nP[2] << ",j=" << nP[1] << ",i=" << nP[0] << ", DATAPACKING=\"POINT\"" << endl;

  for (k = 0; k < nP[2]; k++)
  {
    for (j = 0; j < nP[1]; j++)
    {
      for (i = 0; i < nP[0]; i++)
      {
        file << x[0].p[i] << " " << x[1].p[j] << " " << x[2].p[k] << " " << u[i][j][k] << endl;
      }
    }
  }

  file.close();

  cout << "Program executed successfully" << endl;

  delete[] x[0].p;
  delete[] x[1].p;
  delete[] x[2].p;

  return 0;
}
