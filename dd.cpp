#include "dd.h"
#include "common.h"
#include "engine.h"
#include <algorithm>


static mEdge normalize(Worker* w, mNode& node){

    
    // check for redundant node
    const Complex& w0 = node.children[0].w;
    const Index n0 = node.children[0].n;
    if(std::all_of(node.children.begin(), node.children.end(), [&](const mEdge& e){ return e.w.isApproximatelyEqual(w0) && e.n == n0;} ) ){
        return {w->getComplexFromCache(w0.getValuePair()), n0}; 
    }
    
    // check for all zero weights
    if(std::all_of(node.children.begin(), node.children.end(), [](const mEdge& e){ return e.w.isApproximatelyZero();})){
        return {w->getComplexFromCache({0.0,0.0}), TERMINAL};
    }

    auto result = std::max_element(node.children.begin(), node.children.end(), [](const mEdge& lhs, const mEdge& rhs){
        return lhs.w < rhs.w;   
    });
    const Complex& max_weight = result->w;
    assert(!max_weight.isApproximatelyZero() && max_weight.mag2()!=0.0);
    const std::size_t idx = std::distance(node.children.begin(), result);

    for(int i = 0; i < 4; i++){
        double_pair r = Complex::div(node.children[i].w, max_weight);
        // here we insert into the worker-local ComplexTable
        node.children[i].w = w->getComplexFromTable(r);
    }

    assert(std::all_of(node.children.begin(), node.children.end(), [](const mEdge& e){
        return e.w.src == Complex::Table && e.w.mag() <= 1.0;
    }));
    
    //write the node into the UniqueTable
    
    Index n = w->uniquefy(node);
    return {w->getComplexFromCache(max_weight.getValuePair()), n};
    

}

mEdge makeEdge(Worker* w, Qubit q, const std::array<mEdge, 4> c){
    
    //safety check
    std::for_each(c.begin(), c.end(), [](const mEdge& e){
        assert(e.w.src == Complex::Cache);
    });
    
    mNode node{ .v = q, .children = c };

    
    mEdge e =  normalize(w, node); 

    assert(e.w.src == Complex::Cache);
    return e;
}
mNode* mEdge::getNode() const {
    return uniqueTable.get_data(n);
}

Qubit mEdge::getVar() const {
    return this->getNode()->v;
}

bool mEdge::isTerminal() const {
    return n == 0;
}

static double_pair mul(const double_pair& lhs, const double_pair& rhs) {
    
    const auto lr = lhs.first;
    const auto li = lhs.second;
    const auto rr = rhs.first;
    const auto ri = rhs.second;

    return { lr * rr - li * ri, lr * ri + li * rr};
}

static void printMatrix2(const mEdge& edge, size_t col, size_t row, double_pair w, uint64_t left,  double_pair** m){
    
    double_pair wp = mul(edge.w.getValuePair(), w);  
    
    if(edge.isTerminal() && left == 0){
        m[row][col] = wp;
        return;
    }else if (edge.isTerminal()){
        col = col << left;
        row = row << left;

        for(std::size_t i = 0; i < (1<<left); i++){
            for(std::size_t j = 0; j < (1<<left); j++){
                m[row|i][col|j] = wp;  
            }
        }
        return;
    }

    mNode* node = uniqueTable.get_data(edge.n);
    printMatrix2(node->getEdge(0), (col<<1)|0, (row<<1)|0, wp, left-1, m);
    printMatrix2(node->getEdge(1), (col<<1)|1, (row<<1)|0, wp, left-1, m);
    printMatrix2(node->getEdge(2), (col<<1)|0, (row<<1)|1, wp, left-1, m);
    printMatrix2(node->getEdge(3), (col<<1)|1, (row<<1)|1, wp, left-1, m);
    

}
void mEdge::printMatrix() const {
    Qubit q = this->getVar();   
    std::size_t dim = 1 << (q+1);

    double_pair** matrix = new double_pair*[dim];
    for(std::size_t i = 0; i < dim; i++) matrix[i] = new double_pair[dim];


    printMatrix2(*this, 0, 0, {1.0,0.0}, q+1, matrix);

    for(size_t i = 0 ; i < dim; i++){
        for(size_t j = 0; j < dim; j++){
            std::cout<<matrix[i][j].first<<"+"<<matrix[i][j].second<<"i ";
        }
        std::cout<<"\n";

    }

    for(size_t i = 0; i < dim; i++){
        delete[] matrix[i];
    }
    delete[] matrix;
}



mEdge makeIdent(Worker* w, Qubit q){
    Index n = TERMINAL; 
    Complex one = w->getComplexFromCache({1.0, 0.0});

    Complex zero = w->getComplexFromCache({0.0, 0.0});
    mEdge e; e.n = n;
    for(Qubit i = 0; i <= q; i++){
       e = makeEdge(w, i, {{{one, e.n},{zero, TERMINAL},{zero, TERMINAL},{one, e.n}}}); 
       if(i != q) w->returnComplexToCache(std::move(e.w));
    }


    w->returnComplexToCache(std::move(one));
    w->returnComplexToCache(std::move(zero));

    return e;


}

