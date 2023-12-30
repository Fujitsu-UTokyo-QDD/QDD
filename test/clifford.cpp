#include "algorithms/grover.hpp"
#include "task.h"
#include <boost/fiber/buffered_channel.hpp>
#include <boost/fiber/future/future.hpp>
#include <boost/fiber/future/packaged_task.hpp>
#include <boost/fiber/operations.hpp>
#include <chrono>
#include <condition_variable>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <mutex>
#include <thread>
#if OPENMP
#include <omp.h>
#endif
using namespace std::chrono_literals;
using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;

static mEdge random_clifford(QubitCount qnum, int rand, int target,
                             int control) {
  switch (rand) {
  case 0:
    return H(qnum, target);
  case 1:
    return S(qnum, target);
  case 2:
    return T(qnum, target);
  case 3:
    return CX(qnum, target, control);
  default:
    return H(qnum, target);
  }
}

int main(int argc, char *argv[]) {

  srand(0);
#if OPENMP
  omp_set_num_threads(std::atoi(argv[1]));
  std::cout << "OMP threads = " << omp_get_max_threads() << std::endl;
#endif

  const int nworkers = std::stoi(argv[1]);
  const int nqubits = std::stoi(argv[2]);
  const int gcfreq = std::stoi(argv[3]);
  const int ngates = std::stoi(argv[4]);
  std::cout << "run with " << nworkers << " workers, " << ngates << " gates, "
            << nqubits << " qubits" << std::endl;

  Scheduler s(nworkers, gcfreq);
  for (auto i = 0; i < ngates; i++) {
    int target = rand() % nqubits;
    int control = target;
    while (control == target) {
      control = rand() % nqubits;
    }
    int gate = rand() % 4;
    s.addGate(random_clifford(nqubits, gate, target, control));
  }
  auto t1 = std::chrono::high_resolution_clock::now();
  vEdge result = s.buildCircuit(makeZeroState(nqubits));
  auto t2 = std::chrono::high_resolution_clock::now();
  duration<double, std::micro> ms = t2 - t1;
  std::cout << ms.count() / 1000000 << " seconds" << std::endl;
  return 0;
}
