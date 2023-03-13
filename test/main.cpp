#include "dd.h"
#include <iostream>
int main(int argc, char *argv[]) {
    mEdge e1 = makeHybridGate(4, Hmat, 2, 1);
    mEdge e2 = makeGate(4, Hmat, 2);
    e1.printMatrix();
    e2.printMatrix();
}
