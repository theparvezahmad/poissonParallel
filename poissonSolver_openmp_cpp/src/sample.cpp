#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <omp.h>

const int N = 32;               // Grid size in each dimension
const double L = 1.0;           // Domain size
const double dx = L / N;        // Grid spacing
const double D = 1.0;           // Diffusion coefficient
const double tol = 1e-5;        // Tolerance for convergence
const int maxIterations = 1000; // Maximum iterations

void saveSolution(const std::vector<std::vector<std::vector<double>>> &grid)
{
  std::ofstream outputFile("solution.txt");
  if (!outputFile.is_open())
  {
    std::cerr << "Error: Failed to open the output file." << std::endl;
    return;
  }

  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < N; j++)
    {
      for (int k = 0; k < N; k++)
      {
        outputFile << i << " " << j << " " << k << " " << grid[i][j][k] << std::endl;
      }
    }
  }

  outputFile.close();
}

int main()
{
  std::vector<std::vector<std::vector<double>>> grid(N, std::vector<std::vector<double>>(N, std::vector<double>(N, 0.0)));

  int iteration = 0;
  double error = 0.0;
  double maxError = 0.0;
  bool converged = false;

#pragma omp parallel
  {
    while (!converged && iteration < maxIterations)
    {
      maxError = 0.0;

#pragma omp for reduction(max : maxError)
      for (int i = 1; i < N - 1; i++)
      {
        for (int j = 1; j < N - 1; j++)
        {
          for (int k = 1; k < N - 1; k++)
          {
            double newValue = 0.5 * (grid[i - 1][j][k] + grid[i + 1][j][k] + grid[i][j - 1][k] + grid[i][j + 1][k] + grid[i][j][k - 1] + grid[i][j][k + 1]) / 6.0;

            double localError = std::abs(newValue - grid[i][j][k]);
            maxError = std::max(maxError, localError);

            grid[i][j][k] = newValue;
          }
        }
      }

#pragma omp barrier

#pragma omp single
      {
        iteration++;
        if (maxError < tol)
        {
          converged = true;
        }
      }
    }
  }

#pragma omp barrier

#pragma omp single
  {
    if (converged)
    {
      std::cout << "Converged in " << iteration << " iterations with error " << maxError << std::endl;
    }
    else
    {
      std::cout << "Maximum iterations reached without convergence." << std::endl;
    }

    saveSolution(grid); // Save the solution to a file
  }

  return 0;
}
