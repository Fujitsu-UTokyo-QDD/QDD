#include <cstdio>
#include <vector>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include "gtest/gtest.h"
#include "dd.h"


TEST(MPITest, MPIAllTest){
    // With GoogleTest, all tests with MPI must be included in the single test,
    // unless MPI library makes an error.

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
    world.barrier();
    {
        if(world.rank()==0)
            std::cout << "mv_multiply_MPI" << std::endl;
        vEdge v = makeZeroStateMPI(3, world);
        mEdge h = makeGate(3, Hmat,2);
        v = mv_multiply_MPI(h, v, world, 3, 2);
        v.printVectorMPI(world);
        world.barrier();
        if (world.rank() == 0){
            std::cout << "mv_multipl" << std::endl;
            vEdge v_single = makeZeroState(3);
            mv_multiply(h,v_single).printVector();
        }
    }
    world.barrier();
    {
        std::random_device rd;
        std::mt19937_64 mt(rd());
        for (int i = 0; i<10; i++){
            vEdge v = makeZeroStateMPI(3, world);
            mEdge h0 = makeGate(3, Hmat, 0);
            v = mv_multiply_MPI(h0, v, world, 3, 0);
            mEdge h1 = makeGate(3, Hmat, 1);
            v = mv_multiply_MPI(h1, v, world, 3, 1);
            auto result1 = measureOneCollapsingMPI(world, v, 0, 3, mt);
            mEdge cx = CX(3, 2, 0);
            v = mv_multiply_MPI(cx, v, world, 3, 2);
            auto result2 = measureOneCollapsingMPI(world, v, 2, 3, mt);
            assert(result1 == result2);
        }
    }
    world.barrier();
    {
        std::random_device rd;
        std::mt19937_64 mt(rd());
        for (int i = 0; i<10; i++){
            vEdge v = makeZeroStateMPI(3, world);
            mEdge h0 = makeGate(3, Hmat, 2);
            v = mv_multiply_MPI(h0, v, world, 3, 2);
            mEdge h1 = makeGate(3, Hmat, 1);
            v = mv_multiply_MPI(h1, v, world, 3, 1);
            auto result1 = measureOneCollapsingMPI(world, v, 2, 3, mt);
            mEdge cx = CX(3, 0, 2);
            v = mv_multiply_MPI(cx, v, world, 3, 2);
            auto result2 = measureOneCollapsingMPI(world, v, 0, 3, mt);
            assert(result1 == result2);
        }
    }
}
