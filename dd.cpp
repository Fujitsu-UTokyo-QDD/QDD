#include "dd.h"
#include "common.h"
#include "engine.h"
#include <algorithm>

std::vector<mEdge> identityTable;
AddTable addTable;
MulTable mulTable;

static mEdge normalize(Worker* w, mNode& node){


    
    // check for redundant node
    const std_complex& w0 = node.children[0].w;
    const Index n0 = node.children[0].n;
    if(std::all_of(node.children.begin(), node.children.end(), [&](const mEdge& e){ return e.w == w0 && e.n == n0;} ) ){
        
        return {w0, n0}; 
    }
    
    // check for all zero weights
    if(std::all_of(node.children.begin(), node.children.end(), [](const mEdge& e){ return std::norm(e.w) == 0.0;})){
        return {{0.0,0.0}, TERMINAL};
    }

    auto result = std::max_element(node.children.begin(), node.children.end(), [](const mEdge& lhs, const mEdge& rhs){
            return std::norm(lhs.w) < std::norm(rhs.w);
    });
    std_complex d = result->w;
    const std::size_t idx = std::distance(node.children.begin(), result);
    for(int i = 0; i < 4; i++){
        std_complex r = node.children[i].w/d;
        node.children[i].w = r;
    }


    
    //write the node into the UniqueTable
    
    Index n = w->uniquefy(node);
    return {d, n};
    

}

mEdge makeEdge(Worker* w, Qubit q, const std::array<mEdge, 4> c){
    

    
    mNode node{ .v = q, .children = c };

    
    mEdge e =  normalize(w, node); 

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



static void printMatrix2(const mEdge& edge, size_t col, size_t row, const std_complex& w, uint64_t left,  std_complex** m){
    
    std_complex wp = edge.w * w;
    
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

    std_complex** matrix = new std_complex*[dim];
    for(std::size_t i = 0; i < dim; i++) matrix[i] = new std_complex[dim];


    printMatrix2(*this, 0, 0, {1.0,0.0}, q+1, matrix);

    for(size_t i = 0 ; i < dim; i++){
        for(size_t j = 0; j < dim; j++){
            std::cout<<matrix[i][j]<<" ";
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
       e = makeEdge(w, i, {{{{1.0,0.0}, e.n},{{0.0,0.0}, TERMINAL},{{0.0,0.0}, TERMINAL},{{1.0,0.0}, e.n}}}); 
    }
    
    identityTable[q] = e;
    return e;


}

mEdge makeGate(Worker* w, GateMatrix g, QubitCount q, Qubit target, const Controls& c){
    std::array<mEdge, 4> edges;

    for(auto i = 0; i < 4; i++) edges[i] = mEdge{{g[i].r, g[i].i}, TERMINAL}; 
    
    auto it = c.begin();
    
    mEdge zero = {{0.0,0.0}, TERMINAL};
    mEdge one = {{1.0, 0.0}, TERMINAL};

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

    return e;

}


static Qubit rootVar(const mEdge& lhs, const mEdge& rhs) {
    assert(!(lhs.isTerminal() && rhs.isTerminal()));
    
    return (lhs.isTerminal() || (!rhs.isTerminal() && ( rhs.getVar()> lhs.getVar() )))?  rhs.getVar() : lhs.getVar();
}


mEdge add2(Worker* w, const mEdge& lhs, const mEdge& rhs, int32_t current_var){
    
    if(current_var == -1) {
        assert(lhs.isTerminal() && rhs.isTerminal());
        return {lhs.w + rhs.w, TERMINAL};
    }

    AddQuery query(lhs, rhs, current_var);
    AddQuery* q = addTable.find_or_overwrite(std::hash<AddQuery>()(query), query);
    mEdge result;
    if(q->load_result(w,result)) {
        return result;
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
            x.w = lhs.w * x.w;
        }else{
            x = lhs;
        }
        if(rv == current_var && !rhs.isTerminal()){
            y = rnode->getEdge(i);
            y.w = rhs.w * y.w;
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
        return {lhs.w + rhs.w, TERMINAL};
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
        return {lhs.w * rhs.w, TERMINAL};
    }

    MulQuery query(lhs.n,  rhs.n, current_var);
    MulQuery* q = mulTable.find_or_overwrite(std::hash<MulQuery>()(query), query);
    Index n;
    if(q->load_result(w,n)){
        return {lhs.w * rhs.w, n};
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
                x.w = lhs.w * x.w;
            }else{
                x = lhs; 
            }


            if(rv == current_var && !rhs.isTerminal()){
                y = rnode->getEdge((k<<1) | col);
                y.w = rhs.w * y.w;
            
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
        return {lhs.w * rhs.w, TERMINAL};
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
