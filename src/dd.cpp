#include "dd.h"
#include "common.h"
#include <algorithm>
#include "table.hpp"
#include "graph.hpp"
#include <bitset>

#define SUBTASK_THRESHOLD 5

mNodeTable mUnique(40);
vNodeTable vUnique(40);

std::vector<mEdge> identityTable(40);


mEdge mEdge::one{{1.0, 0.0}, mNode::terminal};
mEdge mEdge::zero{{0.0,0.0}, mNode::terminal};
vEdge vEdge::one{{1.0, 0.0}, vNode::terminal};
vEdge vEdge::zero{{0.0,0.0}, vNode::terminal};

mNode mNode::terminalNode{.v = -1, .children = {}, .next = nullptr};
vNode vNode::terminalNode{.v = -1, .children = {}, .next = nullptr};

static mEdge normalizeM(const mEdge& e){

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

static vEdge normalizeV(const vEdge& e){

    // check for all zero weights
    if(std::all_of(e.n->children.begin(), e.n->children.end(), [](const vEdge& e){ return norm(e.w) == 0.0;})){
        vUnique.returnNode(e.n);
        return vEdge::zero;
    }

    auto result = std::max_element(e.n->children.begin(), e.n->children.end(), [](const vEdge& lhs, const vEdge& rhs){
            return norm(lhs.w) < norm(rhs.w);
    });

    
    std_complex max_weight = result->w;
    const std::size_t idx = std::distance(e.n->children.begin(), result);
    

    for(int i = 0; i < 2; i++){
        std_complex r = e.n->children[i].w/max_weight;
        e.n->children[i].w = r;
    }


    vNode* n = vUnique.lookup(e.n);
    return {max_weight * e.w, n};
    

}

mEdge makeMEdge(Qubit q, const std::array<mEdge, 4>& c){
    

    mNode* node = mUnique.getNode();
    node->v = q;
    node->children = c;

    for(int i = 0; i < 4; i++){
        assert(&node->children[i] != &mEdge::one && (&node->children[i] != &mEdge::zero));
    }

    
    mEdge e =  normalizeM({{1.0,0.0}, node}); 

    assert(e.getVar() == q || e.isTerminal());

    return e;
}

vEdge makeVEdge(Qubit q, const std::array<vEdge, 2>& c){
    

    

    vNode* node = vUnique.getNode();
    node->v = q;
    node->children = c;

    for(int i = 0; i < 2; i++){
        assert(&node->children[i] != &vEdge::one && (&node->children[i] != &vEdge::zero));
    }

    
    vEdge e =  normalizeV({{1.0,0.0}, node}); 

    assert(e.getVar() == q || e.isTerminal());

    return e;
}

mNode* mEdge::getNode() const {
    return n;
}

vNode* vEdge::getNode() const {
    return n;
}

Qubit mEdge::getVar() const {
    return n->v;
}

Qubit vEdge::getVar() const {
    return n->v;
}

bool mEdge::isTerminal() const {
    return n == mNode::terminal;
}

bool vEdge::isTerminal() const {
    return n == vNode::terminal;
}



static void fillMatrix(const mEdge& edge, size_t col, size_t row, const std_complex& w, uint64_t left,  std_complex** m){
    
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
    fillMatrix(node->getEdge(0), (col<<1)|0, (row<<1)|0, wp, left-1, m);
    fillMatrix(node->getEdge(1), (col<<1)|1, (row<<1)|0, wp, left-1, m);
    fillMatrix(node->getEdge(2), (col<<1)|0, (row<<1)|1, wp, left-1, m);
    fillMatrix(node->getEdge(3), (col<<1)|1, (row<<1)|1, wp, left-1, m);
    
    

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


    fillMatrix(*this, 0, 0, {1.0,0.0}, q+1, matrix);

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

std_complex** mEdge::getMatrix(std::size_t* dim) const {
    assert(!this->isTerminal());

    Qubit q = this->getVar();   
    std::size_t d = 1 << (q+1);

    std_complex** matrix = new std_complex*[d];
    for(std::size_t i = 0; i < d; i++) matrix[i] = new std_complex[d];


    fillMatrix(*this, 0, 0, {1.0,0.0}, q+1, matrix);
    if(dim != nullptr) *dim = d;
    return matrix;

}

struct MatrixGuard{
    MatrixGuard(std_complex** m, std::size_t dim): _m(m), _dim(dim){} 
    ~MatrixGuard(){
        for(size_t i = 0; i < _dim; i++){
            delete[] _m[i];
        }
        delete[] _m;
    
    }

    std_complex** _m;
    const std::size_t _dim;
};

bool mEdge::compareNumerically(const mEdge& other) const noexcept {
    if(this->getVar() != other.getVar()) return false;
    
    auto m1 = this->getMatrix(nullptr);
    auto m2 = other.getMatrix(nullptr);
    Qubit q = this->getVar();
    std::size_t dim = 1<<(q+1);
    MatrixGuard mg1(m1 , dim);
    MatrixGuard mg2(m2 , dim);


    for(auto i = 0; i < dim; i++){
        for(auto j = 0; j < dim; j++){
            if(m1[i][j] != m2[i][j])
                return false;
        }
    }
    return true;

}



mEdge makeIdent( Qubit q){

    if(q < 0) return mEdge::one;

    if(identityTable[q].n != nullptr) {
        return identityTable[q];
    }

    mEdge e = makeMEdge(0, {mEdge::one, mEdge::zero, mEdge::zero, mEdge::one});
    for(Qubit i = 1; i <= q; i++){
       e = makeMEdge(i, {{e,mEdge::zero,mEdge::zero,e}}); 
    }
    
    identityTable[q] = e;
    return e;


}

vEdge makeZeroState( QubitCount q){
    vEdge e = makeVEdge(0, {vEdge::one, vEdge::zero});
    for(Qubit i = 1; i < q; i++){
       e = makeVEdge(i, {{e,vEdge::zero}}); 
    }
    return e;
}

vEdge makeOneState(Worker *w, QubitCount q){
    vEdge e = makeVEdge( 0, {vEdge::zero, vEdge::one});
    for(Qubit i = 1; i < q; i++){
       e = makeVEdge( i, {{vEdge::zero, e}}); 
    }
    return e;
}
mEdge makeGate(QubitCount q, GateMatrix g,Qubit target){
    return makeGate(q, g, target, {});
}

mEdge makeGate(QubitCount q,  GateMatrix g, Qubit target, const Controls& c){
    std::array<mEdge, 4> edges;

    for(auto i = 0; i < 4; i++) edges[i] = mEdge{{g[i].real(), g[i].imag()}, mNode::terminal}; 
    
    auto it = c.begin();
    


    Qubit z = 0;
    for(; z < target; z++){
        for(int b1 = 0; b1 < 2; b1++){
            for(int b0 = 0; b0 < 2; b0++){
                std::size_t i = (b1<<1) | b0; 
                if(it != c.end() && it->qubit == z){
                    if(it->type == Control::Type::neg)
                        edges[i] = makeMEdge(z, { edges[i], mEdge::zero, mEdge::zero, (b1 == b0)? makeIdent(z-1): mEdge::zero});
                    else
                        edges[i] = makeMEdge(z, { (b1 == b0)? makeIdent(z-1) : mEdge::zero , mEdge::zero, mEdge::zero, edges[i]});

                }else{
                    edges[i] = makeMEdge(z, {edges[i], mEdge::zero, mEdge::zero, edges[i]});
                }
            }
        }
        if(it != c.end() && it->qubit == z) ++it;
    }

    auto e = makeMEdge(z, edges);


    for( z = z+1 ; z < q; z++){
       if(it != c.end() && it ->qubit == z){
            if(it->type == Control::Type::neg)
                e = makeMEdge(z, { e, mEdge::zero, mEdge::zero, makeIdent(z-1)});
            else
                e = makeMEdge(z, { makeIdent(z-1) , mEdge::zero, mEdge::zero, e});
            ++it;
       }else{
            e = makeMEdge(z, {e, mEdge::zero, mEdge::zero, e}); 
       }
    }

    return e;

}


static Qubit rootVar(const mEdge& lhs, const mEdge& rhs) {
    assert(!(lhs.isTerminal() && rhs.isTerminal()));
    
    return (lhs.isTerminal() || (!rhs.isTerminal() && ( rhs.getVar()> lhs.getVar() )))?  rhs.getVar() : lhs.getVar();
}


mEdge mm_add2(Worker* w, const mEdge& lhs, const mEdge& rhs, int32_t current_var){
    if(lhs.w.isZero()){
        return rhs;
    }else if(rhs.w.isZero()){
        return lhs;
    }
    
    if(current_var == -1) {
        assert(lhs.isTerminal() && rhs.isTerminal());
        return {lhs.w + rhs.w, mNode::terminal};
    }

    mEdge result = w->_addCache.find(lhs,rhs);
    if(result.n != nullptr){
        if(result.w.isApproximatelyZero()){
            return mEdge::zero;
        }else{
            return result;
        }
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

        edges[i] = mm_add2(w, x, y, current_var - 1);
    }


    result =  makeMEdge(current_var, edges);
    w->_addCache.set(lhs, rhs, result);
    

    return result;


}

mEdge mm_add2_no_worker( const mEdge& lhs, const mEdge& rhs, int32_t current_var){
    if(lhs.w.isZero()){
        return rhs;
    }else if(rhs.w.isZero()){
        return lhs;
    }
    
    if(current_var == -1) {
        assert(lhs.isTerminal() && rhs.isTerminal());
        return {lhs.w + rhs.w, mNode::terminal};
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

        edges[i] = mm_add2_no_worker( x, y, current_var - 1);
    }


    mEdge result =  makeMEdge(current_var, edges);
    

    return result;


}
mEdge mm_add(Worker* w, const mEdge& lhs, const mEdge& rhs){
    if(lhs.isTerminal() && rhs.isTerminal()){
        return {lhs.w + rhs.w, mNode::terminal};
    }


    Qubit root = rootVar(lhs, rhs);
    return mm_add2(w, lhs, rhs, root);
}


mEdge mm_add_no_worker( const mEdge& lhs, const mEdge& rhs){
    if(lhs.isTerminal() && rhs.isTerminal()){
        return {lhs.w + rhs.w, mNode::terminal};
    }


    Qubit root = rootVar(lhs, rhs);
    return mm_add2_no_worker( lhs, rhs, root);
}

static void brp(){
}

mEdge mm_multiply2(Worker* w, const mEdge& lhs, const mEdge& rhs, int32_t current_var){

    if(lhs.w.isZero() || rhs.w.isZero()){
        return mEdge::zero;
    }

    if(current_var == -1) {
        if(!lhs.isTerminal() || !rhs.isTerminal()){
            brp();
        }
        assert(lhs.isTerminal() && rhs.isTerminal());
        return {lhs.w * rhs.w, mNode::terminal};
    }
    
    mEdge result = w->_mulCache.find(lhs.n, rhs.n);
    if(result.n != nullptr){
        if(result.w.isApproximatelyZero()){
            return mEdge::zero;
        }else{
            result.w = result.w * lhs.w * rhs.w;
            if(result.w.isApproximatelyZero()) return mEdge::zero;
            else return result;
        }
    }
    

    Qubit lv = lhs.getVar();
    Qubit rv = rhs.getVar();
    assert(lv <= current_var && rv <= current_var);
    mNode* lnode = lhs.getNode();
    mNode* rnode = rhs.getNode();
    mEdge x,y;
    mEdge lcopy = lhs;
    mEdge rcopy = rhs;
    lcopy.w = {1.0,0.0};
    rcopy.w = {1.0, 0.0};


    std::array<mEdge, 4 > edges;

    for(auto i = 0; i < 4; i++){

        std::size_t row = i >> 1;
        std::size_t col = i & 0x1;
        
        std::array<mEdge, 2> product;
        for(auto k = 0; k < 2; k++){
            if(lv == current_var && !lhs.isTerminal()){
                x = lnode->getEdge((row<<1) | k);
            }else{
                x = lhs; 
            }


            if(rv == current_var && !rhs.isTerminal()){
                y = rnode->getEdge((k<<1) | col);
            }else{
                y = rhs;     
            }

            

            product[k] = mm_multiply2(w, x, y, current_var - 1); 
            
        }
        edges[i] = mm_add2(w, product[0], product[1], current_var - 1);
        
    }


    result = makeMEdge(current_var, edges);
    w->_mulCache.set(lhs.n, rhs.n, result);

    result.w = result.w * lhs.w * rhs.w;
    if(result.w.isApproximatelyZero()) return mEdge::zero;
    else return result;

    
}

mEdge mm_multiply2_no_worker(const mEdge& lhs, const mEdge& rhs, int32_t current_var){

    if(lhs.w.isZero() || rhs.w.isZero()){
        return mEdge::zero;
    }

    if(current_var == -1) {
        if(!lhs.isTerminal() || !rhs.isTerminal()){
            brp();
        }
        assert(lhs.isTerminal() && rhs.isTerminal());
        return {lhs.w * rhs.w, mNode::terminal};
    }
    

    

    Qubit lv = lhs.getVar();
    Qubit rv = rhs.getVar();
    assert(lv <= current_var && rv <= current_var);
    mNode* lnode = lhs.getNode();
    mNode* rnode = rhs.getNode();
    mEdge x,y;
    mEdge lcopy = lhs;
    mEdge rcopy = rhs;
    lcopy.w = {1.0,0.0};
    rcopy.w = {1.0, 0.0};


    std::array<mEdge, 4 > edges;

    for(auto i = 0; i < 4; i++){

        std::size_t row = i >> 1;
        std::size_t col = i & 0x1;
        
        std::array<mEdge, 2> product;
        for(auto k = 0; k < 2; k++){
            if(lv == current_var && !lhs.isTerminal()){
                x = lnode->getEdge((row<<1) | k);
            }else{
                x = lhs; 
            }


            if(rv == current_var && !rhs.isTerminal()){
                y = rnode->getEdge((k<<1) | col);
            }else{
                y = rhs;     
            }

            

            product[k] = mm_multiply2_no_worker(x, y, current_var - 1); 
            
        }
        edges[i] = mm_add2_no_worker(product[0], product[1], current_var - 1);
        
    }


    mEdge result = makeMEdge(current_var, edges);

    result.w = result.w * lhs.w * rhs.w;
    if(result.w.isApproximatelyZero()) return mEdge::zero;
    else return result;

    
}

mEdge mm_multiply(Worker* w, const mEdge& lhs, const mEdge& rhs){

    if(lhs.isTerminal() && rhs.isTerminal()){
        return {lhs.w * rhs.w, mNode::terminal};
    }


    Qubit root = rootVar(lhs, rhs);
    mEdge result = mm_multiply2(w, lhs, rhs, root);
    return result;
}


mEdge mm_multiply_no_worker(const mEdge& lhs, const mEdge& rhs){

    if(lhs.isTerminal() && rhs.isTerminal()){
        return {lhs.w * rhs.w, mNode::terminal};
    }


    Qubit root = rootVar(lhs, rhs);
    mEdge result = mm_multiply2_no_worker(lhs, rhs, root);
    return result;
}


static void printVector2(const vEdge& edge, std::size_t row, const std_complex& w, uint64_t left, std_complex* m){
        
    std_complex wp = edge.w * w;
    
    if(edge.isTerminal() && left == 0){
        m[row] = wp;
        return;
    }else if (edge.isTerminal()){
        row = row << left;

        for(std::size_t i = 0; i < (1<<left); i++){
            m[row|i] = wp;  
        }
        return;
    }
    
    vNode* node = edge.getNode();
    printVector2(node->getEdge(0), (row<<1)|0, wp, left-1, m);
    printVector2(node->getEdge(1), (row<<1)|1, wp, left-1, m);

}


mEdge mm_kronecker2(Worker* w, const mEdge& lhs, const mEdge& rhs){
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
        edges[i] = mm_kronecker2(w, x, rhs);
    }

    mEdge ret = makeMEdge(lv + rv + 1, edges);
    return ret;
}


mEdge mm_kronecker(Worker* w, const mEdge& lhs, const mEdge& rhs){
   if(lhs.isTerminal() && rhs.isTerminal()){
        return {lhs.w * rhs.w, mNode::terminal};
   }

   return mm_kronecker2(w, lhs, rhs);

}



vEdge vv_add2(Worker* w, const vEdge& lhs, const vEdge& rhs, int32_t current_var){
    if(lhs.w.isZero()){
        return rhs;
    }else if(rhs.w.isZero()){
        return lhs;
    }
    
    if(current_var == -1) {
        assert(lhs.isTerminal() && rhs.isTerminal());
        return {lhs.w + rhs.w, vNode::terminal};
    }

    vEdge result = w->_addCache.find(lhs,rhs);
    if(result.n != nullptr){
        return result;
    }


    vEdge x, y;

    Qubit lv = lhs.getVar();
    Qubit rv = rhs.getVar();
    vNode* lnode = lhs.getNode();
    vNode* rnode = rhs.getNode();
    std::array<vEdge, 2> edges;

    for(auto i = 0; i < 2; i++){
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

        edges[i] = vv_add2(w, x, y, current_var - 1);
    }

    result =  makeVEdge( current_var, edges);
    w->_addCache.set(lhs, rhs, result);
    

    return result;


}

vEdge vv_add(Worker* w, const vEdge& lhs, const vEdge& rhs){
    if(lhs.isTerminal() && rhs.isTerminal()){
        return {lhs.w + rhs.w, vNode::terminal};
    }

    // assume lhs and rhs are the same size vector.
    assert(lhs.getVar()==rhs.getVar());
    return vv_add2(w, lhs, rhs, lhs.getVar());
}

vEdge vv_kronecker2(Worker* w, const vEdge& lhs, const vEdge& rhs){
    if (lhs.isTerminal()){
        return {lhs.w * rhs.w, rhs.n};
    }
    if (rhs.isTerminal()){
        return {lhs.w * rhs.w, lhs.n};
    }

    std::array<vEdge, 2> edges;
    vEdge x;
    vNode* lnode = lhs.getNode();


    Qubit lv = lhs.getVar();
    Qubit rv = rhs.getVar();
    for(auto i = 0; i < 2; i++){
        x = lnode->getEdge(i);
        x.w = lhs.w * x.w;
        edges[i] = vv_kronecker2(w, x, rhs);
    }

    vEdge ret = makeVEdge(lv + rv + 1, edges);
    return ret;
}

vEdge vv_kronecker(Worker* w, const vEdge& lhs, const vEdge& rhs){
    if(lhs.isTerminal() && rhs.isTerminal()){
        return {lhs.w * rhs.w, vNode::terminal};
   }

   return vv_kronecker2(w, lhs, rhs);
}


void vEdge::printVector() const {
    if(this->isTerminal()) {
        std::cout<<this->w<<std::endl;
        return;
    }
    Qubit q = this->getVar();   
    std::size_t dim = 1 << (q+1);

    std_complex* vector = new std_complex[dim];

    printVector2(*this, 0, {1.0,0.0}, q+1, vector);

    for(size_t i = 0 ; i < dim; i++){
        std::cout<<vector[i]<<" ";
        std::cout<<"\n";

    }
    std::cout<<std::endl;

    delete[] vector;
}

static void fillVector(const vEdge& edge, std::size_t row, const std_complex& w, uint64_t left, std_complex* m){
        
    std_complex wp = edge.w * w;
    
    if(edge.isTerminal() && left == 0){
        m[row] = wp;
        return;
    }else if (edge.isTerminal()){
        row = row << left;

        for(std::size_t i = 0; i < (1<<left); i++){
            m[row|i] = wp;  
        }
        return;
    }
    
    vNode* node = edge.getNode();
    fillVector(node->getEdge(0), (row<<1)|0, wp, left-1, m);
    fillVector(node->getEdge(1), (row<<1)|1, wp, left-1, m);

}

std_complex* vEdge::getVector(std::size_t* dim) const{
    assert(!this->isTerminal());


    Qubit q = this->getVar();   
    std::size_t d = 1 << (q+1);

    std_complex* vector = new std_complex[d];
    fillVector(*this, 0, {1.0,0.0}, q+1, vector);
    if(dim != nullptr) *dim = d;
    return vector;

}



vEdge mv_multiply2(Worker* w, const mEdge& lhs, const vEdge& rhs, int32_t current_var){
    /*
    std::cout<<"lhs\n";
    lhs.printMatrix();
    std::cout<<"rhs\n";
    rhs.printVector();
    std::cout<<std::endl;
*/
    if(lhs.w.isZero() || rhs.w.isZero()){
        return vEdge::zero;
    }

    if(current_var == -1) {
        assert(lhs.isTerminal() && rhs.isTerminal());
        return {lhs.w * rhs.w, vNode::terminal};
    }
    
    vEdge result = w->_mulCache.find(lhs.n, rhs.n);
    if(result.n != nullptr){
        return {result.w * lhs.w * rhs.w, result.n};
    }
    

    Qubit lv = lhs.getVar();
    Qubit rv = rhs.getVar();
    mNode* lnode = lhs.getNode();
    vNode* rnode = rhs.getNode();
    mEdge x;
    vEdge y;
    mEdge lcopy = lhs;
    vEdge rcopy = rhs;
    lcopy.w = {1.0, 0.0};
    rcopy.w = {1.0, 0.0};


    std::array<vEdge, 2> edges;

    for(auto i = 0; i < 2; i++){
        std::array<vEdge, 2> product;
        for(auto k = 0; k < 2; k++){
            if(lv == current_var && !lhs.isTerminal()){
                x = lnode->getEdge((i<<1) | k);
            }else{
                x = lhs; 
            }


            if(rv == current_var && !rhs.isTerminal()){
                y = rnode->getEdge(k);
            }else{
                y = rhs;     
            }

            
            product[k] = mv_multiply2(w, x, y, current_var - 1); 
            
        }

        edges[i] = vv_add2(w, product[0], product[1], current_var - 1);
    }

    result = makeVEdge(current_var, edges);
    w->_mulCache.set(lhs.n, rhs.n, result);

    return {result.w * lhs.w * rhs.w, result.n};
    
}

vEdge mv_multiply(Worker *w, const mEdge& lhs, const vEdge& rhs){
    if(lhs.isTerminal() && rhs.isTerminal()){
        return {lhs.w * rhs.w, vNode::terminal};
    }

    // assume lhs and rhs are the same length.
    assert(lhs.getVar() == rhs.getVar());
    vEdge v = mv_multiply2(w, lhs, rhs, lhs.getVar());
    return v;
}




vEdge vv_multiply(Worker* w, const vEdge& lhs, const vEdge& rhs){
    return vEdge{}; 
}




std::string measureAll(vEdge& rootEdge, const bool collapse, std::mt19937_64& mt, float epsilon) {
    std::size_t dim; 
    std_complex* vs = rootEdge.getVector(&dim);
    std::size_t max_idx = 0;
    float max_val = vs[max_idx].mag2();
    for(auto i = 1; i < dim; i++){
        if(vs[i].mag2() > max_val){
            max_val = vs[i].mag2();
            max_idx = i;
        }
    }
    std::bitset<10> b(max_idx);
    return b.to_string();

}

mEdge makeSwap(QubitCount q, Qubit target0, Qubit target1){
    Controls c1{Control{target0, Control::Type::pos}};
    mEdge e1 = makeGate(q, Xmat, target1, c1);

    Controls c2{Control{target1, Control::Type::pos}};
    mEdge e2 = makeGate(q, Xmat, target0, c2);

    mEdge e3 = mm_multiply_no_worker(e2, e1);
    e3 = mm_multiply_no_worker(e1, e3);

    return e3;
}
