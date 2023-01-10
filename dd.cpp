#include "dd.h"
#include "common.h"
#include "engine.h"
#include <algorithm>

std::vector<mEdge> identityTable;
AddTable addTable(16* AddTable::region_size, 256*AddTable::region_size);
MulTable mulTable(16* MulTable::region_size, 128*MulTable::region_size);


static mEdge normalize(Worker* w, mNode& node){

    auto d0 = node.children[0].w.getValuePair();
    auto d1 = node.children[1].w.getValuePair();
    auto d2 = node.children[2].w.getValuePair();
    auto d3 = node.children[3].w.getValuePair();
    
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
    //Complex max_weight = w->getComplexFromCache(result->w.getValuePair());
    double_pair d = result->w.getValuePair();
    //assert(!max_weight.isApproximatelyZero() && max_weight.mag2()!=0.0);
    const std::size_t idx = std::distance(node.children.begin(), result);
    for(int i = 0; i < 4; i++){
        double_pair r = Complex::div(node.children[i].w, d);
        //w->returnComplexToCache(node.children[i].w);
        node.children[i].w = w->getComplexFromTable(r);
    }


    
    //write the node into the UniqueTable
    
    Index n = w->uniquefy(node);
    return {w->getComplexFromCache(d), n};
    

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

    mEdge e; e.n = TERMINAL;
    for(Qubit i = 0; i <= q; i++){
       e = makeEdge(w, i, {{{Complex::one, e.n},{Complex::zero, TERMINAL},{Complex::zero, TERMINAL},{Complex::one, e.n}}}); 
       ComplexReturner r(w, e.w);
    }
    
    e.w = Complex::one;


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
    
    mEdge zero = {Complex::zero, TERMINAL};
    mEdge one = {Complex::one, TERMINAL};

    Qubit z = 0;
    for(; z < target; z++){
        for(int b1 = 0; b1 < 2; b1++){
            for(int b0 = 0; b0 < 2; b0++){
                std::size_t i = (b1<<1) | b0; 
                if(it != c.end() && *it == z){
                    ComplexReturner r(w, edges[i].w);
                    //positive control 
                    if(z == 0)
                        edges[i] = makeEdge(w, z, { (b1 == b0)? one : zero , zero, zero, edges[i]});
                    else
                        edges[i] = makeEdge(w, z, { (b1 == b0)? makeIdent(w, z-1) : zero , zero, zero, edges[i]});

                }else{
                    ComplexReturner r(w, edges[i].w);
                    edges[i] = makeEdge(w, z, {edges[i], zero, zero, edges[i]});
                }
            }
        }
        if(it != c.end() && *it == z) ++it;
    }

    auto e = makeEdge(w, z, edges);
    {
        ComplexReturner r0(w, edges[0].w);
        ComplexReturner r1(w, edges[1].w);
        ComplexReturner r2(w, edges[2].w);
        ComplexReturner r3(w, edges[3].w);
    }

    for( z = z+1 ; z < q; z++){
       if(it != c.end() && *it == z){
           ComplexReturner r(w, e.w);
           e = makeEdge(w, z, {makeIdent(w, z-1), zero, zero, e});
           ++it;
       }else{
           ComplexReturner r(w, e.w);
            e = makeEdge(w, z, {e, zero, zero, e}); 
       }
    }

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
    AddQuery* q = addTable.get_data(addTable.find_or_insert({.lhs = lhs, .rhs = rhs, .current_var = current_var}));
    mEdge result;
    if(q->load_result(w,result)) return result;

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

    result =  makeEdge(w, current_var, edges);

    q->set_result(w,result);

    return result;


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
    end = std::min(jobs.size(), end);

    mEdge result = jobs[start]->getResult();

    for(auto i = start+1; i < end; i++){
       result = add(w, result, jobs[i]->getResult()); 
    }

    return result;

}

mEdge multiply2(Worker* w, const mEdge& lhs, const mEdge& rhs, int32_t current_var){
    if(current_var == -1) {
        assert(lhs.isTerminal() && rhs.isTerminal());
        double_pair r = Complex::mul(lhs.w, rhs.w);
        return {w->getComplexFromCache(r), TERMINAL};
    }
    MulQuery* q = mulTable.get_data(mulTable.find_or_insert({.lhs = lhs.n, .rhs = rhs.n, .current_var = current_var}));
    Index n;
    if(q->load_result(w,n)){
        double_pair p = Complex::mul(lhs.w, rhs.w);
        return {w->getComplexFromCache(p), n};
    }

    mEdge x, y;

    Qubit lv = lhs.getVar();
    Qubit rv = rhs.getVar();
    mNode* lnode = lhs.getNode();
    mNode* rnode = rhs.getNode();

    std::array<mEdge, 4 > edges;
    for(auto i = 0; i < 4; i++){
        std::size_t row = i >> 1;
        std::size_t col = i & 0x1;
        
        std::array<mEdge, 2> product;
        for(auto k = 0; k < 2; k++){
            if(lv == current_var && !lhs.isTerminal()){
                x = lnode->getEdge((row<<1) | k);
                double_pair xv = Complex::mul(lhs.w, x.w);
                x.w = w->getComplexFromCache(xv);
            }else{
                x = lhs; 
            }


            if(rv == current_var && !rhs.isTerminal()){
                y = rnode->getEdge((k<<1) | col);
                double_pair yv = Complex::mul(rhs.w, y.w);
                y.w = w->getComplexFromCache(yv);
            
            }else{
                y = rhs;     
            }
            
            product[k] = multiply2(w, x, y, current_var - 1); 
            
        }

        edges[i] = add2(w, product[0], product[1], current_var - 1);
    }

    mEdge result = makeEdge(w, current_var, edges);
    q->set_result(w, result.n);
    return result;
    
}


mEdge multiply(Worker* w, const mEdge& lhs, const mEdge& rhs){
    if(lhs.isTerminal() && rhs.isTerminal()){
        double_pair d = Complex::mul(lhs.w, rhs.w);
        return {w->getComplexFromCache(d), TERMINAL};
    }


    Qubit root = rootVar(lhs, rhs);
    return multiply2(w, lhs, rhs, root);

}

mEdge mulSerial(Worker* w, const std::vector<Job*> jobs, std::size_t start, std::size_t end){
    end = std::min(jobs.size(), end);

    mEdge result = jobs[start]->getResult();

    for(auto i = start+1; i < end; i++){
       result = multiply(w, result, jobs[i]->getResult()); 
    }

    return result;

}
