#include "dd.h"
#include <iostream>
int main(int argc, char *argv[]) {
    mEdge e1 = makeHybridGate(4, Hmat, 2, 1, Controls{Control{0}, Control{1}});
    mEdge e2 = makeGate(4, Hmat, 2, Controls{Control{0}, Control{1}});
    e1.printMatrix();
    e2.printMatrix();

    vEdge v = makeHybridOneState(4, 2);
    v.printVector();
}
