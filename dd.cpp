#include "dd.h"
#include "common.h"
#include "engine.h"
#include <algorithm>
#include "table.hpp"


mNodeTable mUnique(20);

std::vector<mEdge> identityTable;
AddTable addTable;
MulTable mulTable;


mEdge mEdge::one{{1.0, 0.0}, mNode::terminal};
mEdge mEdge::zero{{0.0,0.0}, mNode::terminal};

mNode mNode::terminalNode{.v = -1, .children = {}, .next = nullptr};

static mEdge normalize(Worker* w,  const mEdge& e){
    /*

    const mEdge& e0 = e.n->children[0]; 

    if(std::all_of(e.n->children.begin(), e.n->children.end(), [&](const mEdge& e){ return e == e0;} ) ){
        mNode* n = e0.n;
        mUnique.returnNode(e.n);
        return {e.w * e0.w ,  n}; 
    }
    */
    // check for all zero weights
    if(std::all_of(e.n->children.begin(), e.n->children.end(), [](const mEdge& e){ return norm(e.w) == 0.0;})){
        mUnique.returnNode(e.n);
        return mEdge::zero;
    }

    auto result = std::max_element(e.n->children.begin(), e.n->children.end(), [](const mEdge& lhs, const mEdge& rhs){
            return norm(lhs.w) < norm(rhs.w);
    });

    
    std_complex max_weight = result->w;
    const std::size_t idx = std::distance(e.n->children.begin(), result);
    

    for(int i = 0; i < 4; i++){
        std_complex r = e.n->children[i].w/max_weight;
        e.n->children[i].w = r;
    }


    mNode* n = mUnique.lookup(e.n);
    return {max_weight * e.w, n};
    

}

mEdge makeEdge(Worker* w, Qubit q, const std::array<mEdge, 4>& c){
    

    

    mNode* node = mUnique.getNode();
    node->v = q;
    node->children = c;

    for(int i = 0; i < 4; i++){
        assert(&node->children[i] != &mEdge::one && (&node->children[i] != &mEdge::zero));
    }

    
    mEdge e =  normalize(w, {{1.0,0.0}, node}); 

    assert(e.getVar() == q || e.isTerminal());

    return e;
}

mNode* mEdge::getNode() const {
    return n;
}

vNode* vEdge::getNode() const {
    return nullptr;
}

Qubit mEdge::getVar() const {
    return n->v;
}

Qubit vEdge::getVar() const {
    return 0;
}

bool mEdge::isTerminal() const {
    return n == mNode::terminal;
}

bool vEdge::isTerminal() const {
    return false;
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
    
    mNode* node = edge.getNode();
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
    std::cout<<std::endl;

    for(size_t i = 0; i < dim; i++){
        delete[] matrix[i];
    }
    delete[] matrix;
}



mEdge makeIdent(Worker* w, Qubit q){
    if(!identityTable[q].isTerminal()) {
        return identityTable[q];
    }

    mEdge e = makeEdge(w, 0, {mEdge::one, mEdge::zero, mEdge::zero, mEdge::one});
    for(Qubit i = 1; i <= q; i++){
       e = makeEdge(w, i, {{e,mEdge::zero,mEdge::zero,e}}); 
    }
    
    identityTable[q] = e;
    return e;


}

mEdge makeGate(Worker* w, GateMatrix g, QubitCount q, Qubit target, const Controls& c){
    std::array<mEdge, 4> edges;

    for(auto i = 0; i < 4; i++) edges[i] = mEdge{{g[i].r, g[i].i}, mNode::terminal}; 
    
    auto it = c.begin();
    


    Qubit z = 0;
    for(; z < target; z++){
        for(int b1 = 0; b1 < 2; b1++){
            for(int b0 = 0; b0 < 2; b0++){
                std::size_t i = (b1<<1) | b0; 
                if(it != c.end() && *it == z){
                    //positive control 
                    if(z == 0)
                        edges[i] = makeEdge(w, z, { (b1 == b0)? mEdge::one : mEdge::zero , mEdge::zero, mEdge::zero, edges[i]});
                    else
                        edges[i] = makeEdge(w, z, { (b1 == b0)? makeIdent(w, z-1) : mEdge::zero , mEdge::zero, mEdge::zero, edges[i]});

                }else{
                    edges[i] = makeEdge(w, z, {edges[i], mEdge::zero, mEdge::zero, edges[i]});
                }
            }
        }
        if(it != c.end() && *it == z) ++it;
    }

    auto e = makeEdge(w, z, edges);


    for( z = z+1 ; z < q; z++){
       if(it != c.end() && *it == z){
           e = makeEdge(w, z, {makeIdent(w, z-1), mEdge::zero, mEdge::zero, e});
           ++it;
       }else{
            e = makeEdge(w, z, {e, mEdge::zero, mEdge::zero, e}); 
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
        return {lhs.w + rhs.w, mNode::terminal};
    }

    Query query(lhs, rhs, current_var);
    Query* q = addTable.find_or_overwrite(std::hash<Query>()(query), query);
    mEdge result;
    if(q->load_result(w,result)) {
        if(w!=nullptr)
            w->addCacheHit++;
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
        return {lhs.w + rhs.w, mNode::terminal};
    }


    Qubit root = rootVar(lhs, rhs);
    return add2(w, lhs, rhs, root);
}

mEdge addSerial(Worker* w,  const std::vector<Job*>& jobs, std::size_t start, std::size_t end){
    //end = std::min(jobs.size(), end);

    mEdge result = jobs[start]->getResult();

    for(auto i = start+1; i < end; i++){
       result = add(w, result, jobs[i]->getResult()); 
    }

    return result;

}

mEdge multiply2(Worker* w, const mEdge& lhs, const mEdge& rhs, int32_t current_var){
    if(current_var == -1) {
        assert(lhs.isTerminal() && rhs.isTerminal());
        return {lhs.w * rhs.w, mNode::terminal};
    }

    mEdge lcopy = lhs;
    mEdge rcopy = rhs;
    lcopy.w = {1.0,0.0};
    rcopy.w = {1.0, 0.0};

    Query query(lcopy,  rcopy, current_var);
    Query* q = mulTable.find_or_overwrite(std::hash<Query>()(query), query);
    mEdge result;
    if(q->load_result(w,result)){
        if(w!= nullptr)
            w->mulCacheHit++;
        return {result.w * lhs.w * rhs.w, result.n};
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
                //x.w = lhs.w * x.w;
            }else{
                x = lhs; 
            }


            if(rv == current_var && !rhs.isTerminal()){
                y = rnode->getEdge((k<<1) | col);
                //y.w = rhs.w * y.w;
            
            }else{
                y = rhs;     
            }
            
            product[k] = multiply2(w, x, y, current_var - 1); 
            
        }

        edges[i] = add2(w, product[0], product[1], current_var - 1);
    }

    result = makeEdge(w, current_var, edges);

    q->set_result(w, result);
    return {result.w * lhs.w * rhs.w, result.n};
    
}


mEdge multiply(Worker* w, const mEdge& lhs, const mEdge& rhs){
    if(lhs.isTerminal() && rhs.isTerminal()){
        return {lhs.w * rhs.w, mNode::terminal};
    }


    Qubit root = rootVar(lhs, rhs);
    return multiply2(w, lhs, rhs, root);

}

mEdge mulSerial(Worker* w,  const std::vector<Job*>& jobs, std::size_t start, std::size_t end){
    //end = std::min(jobs.size(), end);

    mEdge result = jobs[start]->getResult();

    for(auto i = start+1; i < end; i++){
       result = multiply(w, result, jobs[i]->getResult()); 
    }

    return result;

}



static void printVector2(const vEdge& e, std::size_t i, const std_complex& w, uint64_t left, std_complex* m){
    
    const std_complex ww = w * e.w;
    if(e.isTerminal() && left == 0){
        m[i] = ww;
        return;
    }else if(e.isTerminal()){
        i = i << left;

        for(auto b = 0; b < (1<<left); b++){
            m[i|b] = ww;
        }

        return;
    }
    /* 
    vNode* node = v_uniqueTable.get_data(e.n);
    printVector2(node->getEdge(0), (i<<1)|0, ww, left-1, m );
    printVector2(node->getEdge(1), (i<<1)|1, ww, left-1, m );
*/


}
mEdge kronecker2(Worker* w, const mEdge& lhs, const mEdge& rhs){
    if (lhs.isTerminal()){
        return {lhs.w * rhs.w, rhs.n};
    }

    std::array<mEdge, 4> edges;
    mEdge x;
    mNode* lnode = lhs.getNode();


    Qubit lv = lhs.getVar();
    Qubit rv = rhs.getVar();
    for(auto i = 0; i < 4; i++){
        x = lnode->getEdge(i);
        x.w = lhs.w * x.w;
        edges[i] = kronecker2(w, x, rhs);
    }

    mEdge ret = makeEdge(w,  lv + rv + 1, edges);
    return ret;
}

mEdge kronecker(Worker* w, const mEdge& lhs, const mEdge& rhs){
   if(lhs.isTerminal() && rhs.isTerminal()){
        return {lhs.w * rhs.w, mNode::terminal};
   }

   return kronecker2(w, lhs, rhs);

}



void vEdge::printVector() const {
    if(this->isTerminal()){
        std::cout<<this->w<<std::endl;
        return;
    }
    Qubit q = this->getVar();   
    std::size_t dim = 1 << (q+1);

    std_complex* matrix = new std_complex[dim];
    
    printVector2(*this, 0, {1.0,0.0}, q+1, matrix);
    
    for(auto i = 0; i < dim; i++){
        std::cout<<matrix[i]<<std::endl;
    }

    delete[] matrix;
}

vEdge mEdge::get_column(std::size_t col)  const {
    return vEdge();

}
