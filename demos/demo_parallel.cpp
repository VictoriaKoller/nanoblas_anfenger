#include <iostream>
#include <matrix.hpp>
#include <inverse.hpp>
#include <lapack_interface.hpp>
#include <algorithm>
#include <taskmanager.hpp>
#include <timer.hpp>

using namespace std;
using namespace nanoblas;
using namespace ASC_HPC;

int main() {
  const size_t n = 1000;
    Matrix<double, nanoblas::ColMajor> A(n,n);
    Matrix<double, nanoblas::ColMajor> B(n,n);
    Matrix<double, nanoblas::ColMajor> C(n,n);


    A = 2.0;
    B = 3.0;
    C = 0.0;


  StartWorkers(3);                    // 3 Worker + Hauptthread = 4 Threads

  // Vite-Timer starten
  // TimerRegion reg("matmat_parallel");  // oder wie es im ASC-HPC-Paket heißt

  addMatMat_parallel(A, B, C);

  // TimerRegion endet beim Destruktor -> Trace für Vite wird geschrieben

  StopWorkers();

  std::cout << C(0,0) << std::endl;   // Plausibilitätscheck 
  return 0;
}
