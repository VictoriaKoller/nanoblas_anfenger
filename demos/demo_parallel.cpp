#include <iostream>
#include <matrix.hpp>
#include <inverse.hpp>
#include <lapack_interface.hpp>
#include <algorithm>
#include "../HP_Anfenger/src/taskmanager.hpp"
#include "../HP_Anfenger/src/timer.hpp"

using namespace std;
using namespace nanoblas;

int main() {
  const size_t n = 1000;
    Matrix<double, RowMajor> A(n,n);
    Matrix<double, RowMajor> B(n,n);
    Matrix<double, RowMajor> C(n,n);


  StartWorkers(3);                    // 3 Worker + Hauptthread = 4 Threads

  // Vite-Timer starten
  // TimerRegion reg("matmat_parallel");  // oder wie es im ASC-HPC-Paket heißt

  addMatMat_parallel(A, B, C);

  // TimerRegion endet beim Destruktor -> Trace für Vite wird geschrieben

  StopWorkers();

  std::cout << C(0,0) << std::endl;   // Plausibilitätscheck (sollte 2*n sein)
}
