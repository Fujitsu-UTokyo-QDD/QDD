#include <cstdio>
#include <vector>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include "gtest/gtest.h"
#include "dd.h"


TEST(MPITest, MPIAllTest){
    int argc=0;
    char **argv;
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;

    {
        std::cout << world.rank() << " / " << world.size() << std::endl;
    }
    world.barrier();
    {
        vEdge e0 = makeZeroStateMPI(4, world);
        e0.printVectorMPI(world);
        if(world.rank()==0)
            std::cout << "----" << std::endl;
        vEdge e1 = makeOneStateMPI(4, world);
        e1.printVectorMPI(world);
    }
    if(world.rank()==0)
            std::cout << "----" << std::endl;
    world.barrier();
    {
            assert(world.size() > 1);
            mEdge cnot = CX(2, 0, 1);
            mEdge rz = RZ(2, 0, 1.0);
            mEdge h = makeGate(2, Hmat, 1);
            mEdge tmp0 = mm_multiply(rz, cnot);
            mEdge tmp = mm_multiply(h, tmp0);
            if (world.rank() == 0)
            {
            tmp.printMatrix();
            std::cout << "----" << std::endl;
        }
        mEdge mat00 = getMPIGate(tmp, 0, 0, world.size());
        mEdge mat01 = getMPIGate(tmp, 0, 1, world.size());
        mEdge mat10 = getMPIGate(tmp, 1, 0, world.size());
        mEdge mat11 = getMPIGate(tmp, 1, 1, world.size());
        if(world.rank()==0){
            mat00.printMatrix();
            mat01.printMatrix();
            mat10.printMatrix();
            mat11.printMatrix();
        }
    }
}
