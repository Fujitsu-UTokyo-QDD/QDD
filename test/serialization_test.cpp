
#include "dd.h"
#include "gtest/gtest.h"
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <cstdio>
#include <fstream>
#include <vector>

TEST(SerializationTest, InitTest2) {

    std::ofstream ofs("filename");

    // create class instance
    vEdge v = makeZeroState(3);

    // save data to archive
    {
        boost::archive::text_oarchive oa(ofs);
        // write class instance to archive
        oa << v;
        // archive and stream closed when destructors are called
    }

    // ... some time later restore the class instance to its orginal state
    vEdge newv;
    {
        // create and open an archive for input
        std::ifstream ifs("filename");
        boost::archive::text_iarchive ia(ifs);
        // read class state from archive
        ia >> newv;
        // archive and stream closed when destructors are called
    }
    std::cout << "input:" << std::endl;
    v.printVector(); std::cout << "output:" << std::endl;
    newv.printVector();
}
