#include "dd.h"
#include "common.h"
#include "engine.h"
#include <algorithm>

std::vector<mEdge> identityTable;

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
    Complex max_weight = result->w;
    assert(!max_weight.isApproximatelyZero() && max_weight.mag2()!=0.0);
    const std::size_t idx = std::distance(node.children.begin(), result);

    for(int i = 0; i < 4; i++){
        double_pair r = Complex::div(node.children[i].w, max_weight);
        // here we insert into the worker-local ComplexTable
        node.children[i].w = w->getComplexFromTable(r);
    }
    /*
    for(auto it = node.children.begin(); it != node.children.end(); ++it){
        if(!(it->w.src == Complex::Table && it->w.mag() <= 1.0)){
            std::cout<<it->w.mag()<<std::endl;
            std::cout<<it->w.src<<std::endl;
            assert(false);
        }
    
    }
    */

    
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
    if(this->isTerminal()) {
        std::cout<<this->w<<std::endl;
        return;
    }
    Qubit q = this->getVar();   
    std::size_t dim = 1 << (q+1);

    double_pair** matrix = new double_pair*[dim];
    for(std::size_t i = 0; i < dim; i++) matrix[i] = new double_pair[dim];


    printMatrix2(*this, 0, 0, {1.0,0.0}, q+1, matrix);

    for(size_t i = 0 ; i < dim; i++){
        for(size_t j = 0; j < dim; j++){
            //std::cout<<matrix[i][j].first<<"+"<<matrix[i][j].second<<"i ";
            std::cout<<matrix[i][j].first<<" ";
        }
        std::cout<<"\n";

    }

    for(size_t i = 0; i < dim; i++){
        delete[] matrix[i];
    }
    delete[] matrix;
}



mEdge makeIdent(Worker* w, Qubit q){
    if(!identityTable[q].isTerminal()) {
        return identityTable[q];
    }
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
    identityTable[q] = e;
    return e;


}

mEdge makeGate(Worker* w, GateMatrix g, QubitCount q, Qubit target, const Controls& c){
    std::array<mEdge, 4> edges;
    std::array<Complex, 4> matrix_entries;
    for(auto i = 0; i < 4; i++) {
        matrix_entries[i] = w->getComplexFromCache({g[i].r, g[i].i});
    }
    for(auto i = 0; i < 4; i++) edges[i] = mEdge{matrix_entries[i], TERMINAL}; 
    
    auto it = c.begin();
    
    mEdge zero = {w->getComplexFromCache({0.0, 0.0}), TERMINAL};
    mEdge one = {w->getComplexFromCache({1.0, 0.0}), TERMINAL};
    Qubit z = 0;
    for(; z < target; z++){
        for(int b1 = 0; b1 < 2; b1++){
            for(int b0 = 0; b0 < 2; b0++){
                std::size_t i = (b1<<1) | b0; 
                if(it != c.end() && *it == z){
                    //positive control 
                    if(z == 0)
                        edges[i] = makeEdge(w, z, { (b1 == b0)? one : zero , zero, zero, edges[i]});
                    else
                        edges[i] = makeEdge(w, z, { (b1 == b0)? makeIdent(w, z-1) : zero , zero, zero, edges[i]});

                }else{
                    edges[i] = makeEdge(w, z, {edges[i], zero, zero, edges[i]});
                }
            }
        }
        if(it != c.end() && *it == z) ++it;
    }

    auto e = makeEdge(w, z, edges);

    for( z = z+1 ; z < q; z++){
       if(it != c.end() && *it == z){
           e = makeEdge(w, z, {makeIdent(w, z-1), zero, zero, e});
           ++it;
       }else{
            e = makeEdge(w, z, {e, zero, zero, e}); 
       }
    }

    w->returnComplexToCache(std::move(zero.w));
    for(auto i = 0; i < 4; i++) w->returnComplexToCache(std::move(matrix_entries[i]));
    return e;

}


static Qubit rootVar(const mEdge& lhs, const mEdge& rhs) {
    assert(!(lhs.isTerminal() && rhs.isTerminal()));
    
    return (lhs.isTerminal() || (!rhs.isTerminal() && ( rhs.getVar()> lhs.getVar() )))?  rhs.getVar() : lhs.getVar();
}


mEdge add2(Worker* w, const mEdge& lhs, const mEdge& rhs, int32_t current_var){
    
    if(current_var == -1) {
        assert(lhs.isTerminal() && rhs.isTerminal());
        double_pair r = Complex::add(lhs.w, rhs.w);
        return {w->getComplexFromCache(r), TERMINAL};
    }


    mEdge x, y;

    Qubit lv = lhs.getVar();
    Qubit rv = rhs.getVar();
    mNode* lnode = lhs.getNode();
    mNode* rnode = rhs.getNode();

    std::array<mEdge, 4> edges;

    for(auto i = 0; i < 4; i++){
        if(lv == current_var && !lhs.isTerminal()){
            x = lnode->getEdge(i);
            double_pair xv = Complex::mul(lhs.w, x.w);
            x.w = w->getComplexFromCache(xv);
        }else{
            x = lhs;
        }
        if(rv == current_var && !rhs.isTerminal()){
            y = rnode->getEdge(i);
            double_pair yv = Complex::mul(rhs.w, y.w);
            y.w = w->getComplexFromCache(yv);
        }else{
            y = rhs;
        }

        edges[i] = add2(w, x, y, current_var - 1);
    }

    return makeEdge(w, current_var, edges);


}

mEdge add(Worker* w, const mEdge& lhs, const mEdge& rhs){
    if(lhs.isTerminal() && rhs.isTerminal()){
        double_pair d = Complex::add(lhs.w, rhs.w);
        return {w->getComplexFromCache(d), TERMINAL};
    }


    Qubit root = rootVar(lhs, rhs);
    return add2(w, lhs, rhs, root);
}

mEdge addSerial(Worker* w, const std::vector<Job*> jobs, std::size_t start, std::size_t end){
    //std::cout<<"enter add serial" <<std::endl; 
    end = std::min(jobs.size(), end);

    mEdge result = jobs[start]->getResult();

    for(auto i = start+1; i < end; i++){
       result = add(w, result, jobs[i]->getResult()); 
    }

    //std::cout<<"leave serial "<< start<<std::endl;
    return result;

}
