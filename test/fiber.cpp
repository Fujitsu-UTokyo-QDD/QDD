#include <boost/fiber/all.hpp>
#include <iostream>

using boost::fibers::fiber;

struct Bar{
    
    Bar(int ii): i(ii){}
    void operator()() const{
        std::cout<<"Bar::i = "<< i<<std::endl;
    
    }
    int i;

};

void fn(const Bar& b) {
    b();
}

int main(int argc, char* argv[]){
    Bar b(26);
    fiber f(std::ref(b));

    return 0;
}
