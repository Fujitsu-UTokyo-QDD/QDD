#include "dd.h"
#include <iostream>
int main(int argc, char *argv[]) {
    mEdge e1 = makeHybridGate(3, Hmat, 1, 2);
    mEdge e2 = makeGate(3, Hmat, 1);
    e1.printMatrix();
    e2.printMatrix();
}
