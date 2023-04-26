#include <cstdio>
#include <vector>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include "gtest/gtest.h"
#include "dd.h"

/*
TEST(MPITest, BasicTest){
    int argc=0;
    char **argv;
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;
    std::cout << world.rank() << " / " << world.size() << std::endl;
}
*/

TEST(MPITest, InitTest){
    int argc=0;
    char **argv;
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;

    vEdge e0 = makeZeroStateMPI(4, world);
    e0.printVectorMPI(world);
    if(world.rank()==0)
        std::cout << "----" << std::endl;
    vEdge e1 = makeOneStateMPI(4, world);
    e1.printVectorMPI(world);
}
