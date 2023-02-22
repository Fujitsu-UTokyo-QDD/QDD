#include "dd.h"
#include "common.h"
#include <algorithm>
#include "table.hpp"
#include "graph.hpp"
#include <bitset>
#include <oneapi/tbb/enumerable_thread_specific.h>
#include <boost/fiber/future/future.hpp>
#include <boost/fiber/future/packaged_task.hpp>

#define SUBTASK_THRESHOLD 5

mNodeTable mUnique(40);
vNodeTable vUnique(40);

std::vector<mEdge> identityTable(40);


mEdge mEdge::one{{1.0, 0.0}, mNode::terminal};
mEdge mEdge::zero{{0.0,0.0}, mNode::terminal};
vEdge vEdge::one{{1.0, 0.0}, vNode::terminal};
vEdge vEdge::zero{{0.0,0.0}, vNode::terminal};

//mNode mNode::terminalNode{.v = -1, .children = {}, .next = nullptr, .ref = MAX_REF };
//vNode vNode::terminalNode{.v = -1, .children = {}, .next = nullptr, .ref = MAX_REF };

mNode mNode::terminalNode = mNode(-1, {}, nullptr, MAX_REF);
vNode vNode::terminalNode = vNode(-1, {}, nullptr, MAX_REF);

oneapi::tbb::enumerable_thread_specific<AddCache> _aCache(40);
oneapi::tbb::enumerable_thread_specific<MulCache> _mCache(40);

static int LIMIT =8;
const int MINUS = 3;

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
    assert(n->v >= -1);
    

    
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

MatrixXcf mEdge::getEigenMatrix(){
    std::size_t dim;
    auto m = getMatrix(&dim);
    MatrixGuard g(m, dim);
    MatrixXcf M(dim,dim);

    for(auto i = 0; i < dim; i++){
        for(auto j = 0; j < dim; j++ ){
           M(i,j) = std::complex<double>{m[i][j].r, m[i][j].i}; 
        }
    }
    return M;
}


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
        assert(identityTable[q].n->v > -1);
        return identityTable[q];
    }

    mEdge e = makeMEdge(0, {mEdge::one, mEdge::zero, mEdge::zero, mEdge::one});
    for(Qubit i = 1; i <= q; i++){
       e = makeMEdge(i, {{e,mEdge::zero,mEdge::zero,e}}); 
    }
    
    identityTable[q] = e;
    //e.incRef();
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
    if(lhs.w.isApproximatelyZero()){
        return rhs;
    }else if(rhs.w.isApproximatelyZero()){
        return lhs;
    }
    
    if(current_var == -1) {
        assert(lhs.isTerminal() && rhs.isTerminal());
        return {lhs.w + rhs.w, mNode::terminal};
    }

    AddCache& local_aCache = _aCache.local();
    mEdge result = local_aCache.find(lhs,rhs);
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
    local_aCache.set(lhs, rhs, result);
    

    return result;


}

mEdge mm_add_fiber2(const mEdge& lhs, const mEdge& rhs, int32_t current_var){
    if(lhs.w.isApproximatelyZero()){
        return rhs;
    }else if(rhs.w.isApproximatelyZero()){
        return lhs;
    }
    
    if(current_var == -1) {
        assert(lhs.isTerminal() && rhs.isTerminal());
        return {lhs.w + rhs.w, mNode::terminal};
    }

    AddCache& local_aCache = _aCache.local();
    mEdge result = local_aCache.find(lhs,rhs);
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

        edges[i] = mm_add_fiber2( x, y, current_var - 1);
    }


    result =  makeMEdge(current_var, edges);
    local_aCache.set(lhs, rhs, result);
    

    return result;


}
mEdge mm_add2_no_worker( const mEdge& lhs, const mEdge& rhs, int32_t current_var){
    if(lhs.w.isApproximatelyZero()){
        return rhs;
    }else if(rhs.w.isApproximatelyZero()){
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

    if(lhs.w.isApproximatelyZero() || rhs.w.isApproximatelyZero()){
        return mEdge::zero;
    }

    if(current_var == -1) {
        if(!lhs.isTerminal() || !rhs.isTerminal()){
            brp();
        }
        assert(lhs.isTerminal() && rhs.isTerminal());
        return {lhs.w * rhs.w, mNode::terminal};
    }
    
    MulCache& local_mCache = _mCache.local();
    mEdge result = local_mCache.find(lhs.n, rhs.n);
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
                x = lcopy; 
            }


            if(rv == current_var && !rhs.isTerminal()){
                y = rnode->getEdge((k<<1) | col);
            }else{
                y = rcopy;     
            }

            

            product[k] = mm_multiply2(w, x, y, current_var - 1); 
            
        }
        edges[i] = mm_add2(w, product[0], product[1], current_var - 1);
        
    }


    result = makeMEdge(current_var, edges);
    local_mCache.set(lhs.n, rhs.n, result);

    result.w = result.w * lhs.w * rhs.w;
    if(result.w.isApproximatelyZero()) return mEdge::zero;
    else return result;

    
}

mEdge mm_test(mEdge lhs, mEdge rhs, int32_t current_var){
    return mEdge();
}
mEdge mm_multiply2_no_worker(const mEdge& lhs, const mEdge& rhs, int32_t current_var){

    if(lhs.w.isApproximatelyZero() || rhs.w.isApproximatelyZero()){
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
mEdge mm_multiply_fiber2(const mEdge& lhs, const mEdge& rhs, int32_t current_var){

    if(lhs.w.isApproximatelyZero() || rhs.w.isApproximatelyZero()){
        return mEdge::zero;
    }

    if(current_var == -1) {

        assert(lhs.isTerminal() && rhs.isTerminal());
        return {lhs.w * rhs.w, mNode::terminal};
    }
    
    MulCache& local_mCache = _mCache.local();
    mEdge result = local_mCache.find(lhs.n, rhs.n);
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

    std::vector<boost::fibers::future<mEdge>> products;
    std::vector<mEdge> products_nofuture;


    for(auto i = 0; i < 4; i++){

        std::size_t row = i >> 1;
        std::size_t col = i & 0x1;
        
        std::array<mEdge, 2> product;
        for(auto k = 0; k < 2; k++){
            if(lv == current_var && !lhs.isTerminal()){
                x = lnode->getEdge((row<<1) | k);
            }else{
                x = lcopy; 
            }


            if(rv == current_var && !rhs.isTerminal()){
                y = rnode->getEdge((k<<1) | col);
            }else{
                y = rcopy;     
            }

            if(current_var > LIMIT){ 
                boost::fibers::packaged_task<mEdge()> pt(std::bind(mm_multiply_fiber2, x, y, current_var - 1)); 
                products.emplace_back(pt.get_future());
                boost::fibers::fiber f(std::move(pt));
                f.detach();
            }else{
                products_nofuture.emplace_back(mm_multiply_fiber2(x,y, current_var - 1));
            }

            
        }
        
    }
    if(current_var > LIMIT){
        assert(products.size() == 8);

        for(int i = 0; i < 8; i += 2){
            edges[i/2] = mm_add_fiber2(products[i].get(), products[i+1].get(), current_var - 1);
        
        }
    }else{
        assert(products_nofuture.size() == 8);
        for(int i = 0; i < 8; i += 2){
            edges[i/2] = mm_add_fiber2(products_nofuture[i], products_nofuture[i+1], current_var - 1);
        
        }

    
    }


    result = makeMEdge(current_var, edges);
    local_mCache.set(lhs.n, rhs.n, result);

    result.w = result.w * lhs.w * rhs.w;
    if(result.w.isApproximatelyZero()) return mEdge::zero;
    if(result.w.isApproximatelyOne()) result.w = {1.0,0.0};
    return result;

    
}

mEdge mm_multiply(Worker* w, const mEdge& lhs, const mEdge& rhs){

    if(lhs.isTerminal() && rhs.isTerminal()){
        return {lhs.w * rhs.w, mNode::terminal};
    }


    Qubit root = rootVar(lhs, rhs);
    mEdge result = mm_multiply2(w, lhs, rhs, root);
    return result;
}

mEdge mm_multiply_fiber(const mEdge& lhs, const mEdge& rhs){

    if(lhs.isTerminal() && rhs.isTerminal()){
        return {lhs.w * rhs.w, mNode::terminal};
    }


    Qubit root = rootVar(lhs, rhs);
    LIMIT = root-MINUS;
    mEdge result = mm_multiply_fiber2( lhs, rhs, root);
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
    if(lhs.w.isApproximatelyZero()){
        return rhs;
    }else if(rhs.w.isApproximatelyZero()){
        return lhs;
    }
    
    if(current_var == -1) {
        assert(lhs.isTerminal() && rhs.isTerminal());
        return {lhs.w + rhs.w, vNode::terminal};
    }
    
    AddCache& local_aCache = _aCache.local();
    vEdge result = local_aCache.find(lhs,rhs);
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
    local_aCache.set(lhs, rhs, result);
    

    return result;


}


vEdge vv_add_fiber2(const vEdge& lhs, const vEdge& rhs, int32_t current_var){
    if(lhs.w.isApproximatelyZero()){
        return rhs;
    }else if(rhs.w.isApproximatelyZero()){
        return lhs;
    }
    
    if(current_var == -1) {
        if(!lhs.isTerminal() || !rhs.isTerminal()){
            brp();
        
        }
        assert(lhs.isTerminal() && rhs.isTerminal());
        return {lhs.w + rhs.w, vNode::terminal};
    }
    
    vEdge result;
#ifdef CACHE
    AddCache& local_aCache = _aCache.local();
    result = local_aCache.find(lhs,rhs);
    if(result.n != nullptr){
        return result;
    }
#endif


    vEdge x, y;

    Qubit lv = lhs.getVar();
    Qubit rv = rhs.getVar();
    vNode* lnode = lhs.getNode();
    vNode* rnode = rhs.getNode();
    std::array<vEdge, 2> edges;

    std::vector<boost::fibers::future<vEdge>> sums;

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
        if(current_var>LIMIT)
        {
            boost::fibers::packaged_task<vEdge()> pt(std::bind(vv_add_fiber2,x,y, current_var - 1));
            sums.emplace_back(pt.get_future());
            boost::fibers::fiber f(std::move(pt));
            f.detach();
        }else{
            edges[i] = vv_add_fiber2(x, y, current_var - 1);
        }
    }
    if(current_var>LIMIT){
        assert(sums.size() == 2);
        edges[0] = sums[0].get();
        edges[1] = sums[1].get();
    }
    result =  makeVEdge( current_var, edges);
#ifdef CACHE
    local_aCache.set(lhs, rhs, result);
#endif
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

static void printVector_sparse2(const vEdge& edge, std::size_t row, const std_complex& w, uint64_t left, std::map<int,std_complex> &m){
        
    std_complex wp = edge.w * w;

    const double thr = 0.001;

    if(edge.isTerminal() && left == 0){
        if(norm(wp) > thr)
            m[row] = wp;
        return;
    }else if (edge.isTerminal()){
        if(norm(wp) > thr){
          row = row << left;
          for(std::size_t i = 0; i < (1<<left); i++){
              m[row|i] = wp;  
          }
        }
        return;
    }
    
    vNode* node = edge.getNode();
    printVector_sparse2(node->getEdge(0), (row<<1)|0, wp, left-1, m);
    printVector_sparse2(node->getEdge(1), (row<<1)|1, wp, left-1, m);

}

void vEdge::printVector_sparse() const {
    if(this->isTerminal()) {
        std::cout<<this->w<<std::endl;
        return;
    }
    Qubit q = this->getVar();   
    std::size_t dim = 1 << (q+1);

    std::map<int, std_complex> map;

    printVector_sparse2(*this, 0, {1.0,0.0}, q+1, map);

    for (auto itr = map.begin(); itr != map.end(); itr++){
        std::cout << std::bitset<30>(itr->first) << ": " << itr->second << std::endl;
    }
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
struct VectorGuard{
    VectorGuard(std_complex* m, std::size_t dim): _m(m), _dim(dim){} 
    ~VectorGuard(){
        delete[] _m;
    
    }

    std_complex* _m;
    const std::size_t _dim;
};

VectorXcf vEdge::getEigenVector(){
    std::size_t dim;
    auto v = getVector(&dim);
    VectorGuard g(v, dim);
    VectorXcf V(dim);

    for(auto i = 0; i < dim; i++){
       V(i) = std::complex<double>{v[i].r, v[i].i}; 
    }
    return V;
}



vEdge mv_multiply2(Worker* w, const mEdge& lhs, const vEdge& rhs, int32_t current_var){
    /*
    std::cout<<"lhs\n";
    lhs.printMatrix();
    std::cout<<"rhs\n";
    rhs.printVector();
    std::cout<<std::endl;
*/
    if(lhs.w.isApproximatelyZero() || rhs.w.isApproximatelyZero()){
        return vEdge::zero;
    }

    if(current_var == -1) {
        if(!lhs.isTerminal() || !rhs.isTerminal()){
            brp();
        }
        assert(lhs.isTerminal() && rhs.isTerminal());
        return {lhs.w * rhs.w, vNode::terminal};
    }
    
    MulCache& local_mCache = _mCache.local();
    vEdge result = local_mCache.find(lhs.n, rhs.n);
    if(result.n != nullptr){
        if(result.w.isApproximatelyZero()){
            return vEdge::zero;
        }else{
            result.w = result.w * lhs.w * rhs.w;
            if(result.w.isApproximatelyZero()) return vEdge::zero;
            else return result;
        }
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
                x = lcopy; 
            }


            if(rv == current_var && !rhs.isTerminal()){
                y = rnode->getEdge(k);
            }else{
                y = rcopy;     
            }

            
            product[k] = mv_multiply2(w, x, y, current_var - 1); 
            
        }

        edges[i] = vv_add2(w, product[0], product[1], current_var - 1);
    }

    result = makeVEdge(current_var, edges);
    local_mCache.set(lhs.n, rhs.n, result);
    result.w = result.w * lhs.w * rhs.w;
    if(result.w.isApproximatelyZero()){ 
        return vEdge::zero;
    }
    else {
        return result;    
    }
}

struct mEdgeRefGuard{
    mEdgeRefGuard(mEdge& ee): e(ee){ e.incRef();}
    ~mEdgeRefGuard(){ e.decRef();};
    mEdge e;
};
struct vEdgeRefGuard{
    vEdgeRefGuard(vEdge& ee): e(ee){ e.incRef();}
    ~vEdgeRefGuard(){ e.decRef();};
    vEdge e;
};

vEdge mv_multiply(Worker *w, mEdge lhs, vEdge rhs){
    //mEdgeRefGuard mg(lhs);
    //vEdgeRefGuard vg(rhs);
    
    if(lhs.isTerminal() && rhs.isTerminal()){
        return {lhs.w * rhs.w, vNode::terminal};
    }

    // assume lhs and rhs are the same length.
    assert(lhs.getVar() == rhs.getVar());
    vEdge v = mv_multiply2(w, lhs, rhs, lhs.getVar());
    return v;
}

vEdge mv_multiply_fiber2(const mEdge& lhs, const vEdge& rhs, int32_t current_var){


    if(lhs.w.isApproximatelyZero() || rhs.w.isApproximatelyZero()){
        return vEdge::zero;
    }

    if(current_var == -1) {

        assert(lhs.isTerminal() && rhs.isTerminal());
        return {lhs.w * rhs.w, vNode::terminal};
    }

    vEdge result;

#ifdef CACHE               
    MulCache& local_mCache = _mCache.local();
#endif
    

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

    
    std::vector<boost::fibers::future<vEdge>> products;
    std::vector<vEdge> products_nofuture;
    std::array<vEdge, 2> edges;

    std::vector<std::variant<vEdge, boost::fibers::future<vEdge>>> products_children;

    for(auto i = 0; i < 2; i++){
        std::array<vEdge, 2> product;
        for(auto k = 0; k < 2; k++){
            if(lv == current_var && !lhs.isTerminal()){
                x = lnode->getEdge((i<<1) | k);
            }else{
                x = lcopy; 
            }


            if(rv == current_var && !rhs.isTerminal()){
                y = rnode->getEdge(k);
            }else{
                y = rcopy;     
            }
    
            {
                if(current_var > LIMIT){
#ifdef CACHE
                    result = local_mCache.find(x.n, y.n);
                    if(result.n != nullptr){
                        if(result.w.isApproximatelyZero()){
                            result = vEdge::zero;
                        }else if(result.w.isApproximatelyOne()){
                            result = vEdge::one;
                        }else{
                            result.w = result.w * lhs.w * rhs.w;
                            if(result.w.isApproximatelyZero()){
                                result = vEdge::zero;
                            }else if (result.w.isApproximatelyOne()){
                                result = vEdge::one;
                            }
                        }

                        products_children.emplace_back(result);
                    }else{
                        boost::fibers::packaged_task<vEdge()> pt(std::bind(mv_multiply_fiber2,x,y, current_var - 1));
                        products_children.emplace_back(pt.get_future());
                        boost::fibers::fiber f(std::move(pt));
                        f.detach();
                    }
#else
                        boost::fibers::packaged_task<vEdge()> pt(std::bind(mv_multiply_fiber2,x,y, current_var - 1));
                        products_children.emplace_back(pt.get_future());
                        boost::fibers::fiber f(std::move(pt));
                        f.detach();
#endif
                }else{
                    products_nofuture.emplace_back(mv_multiply_fiber2(x,y, current_var - 1));
                }
            }
            
        }

    }

if(current_var > LIMIT){
    assert(products_children.size() == 4 );
    for(std::variant<vEdge, boost::fibers::future<vEdge>>& v: products_children){
            if(std::holds_alternative<vEdge>(v)){
                products_nofuture.emplace_back(std::get<vEdge>(v));
            }else{
                products_nofuture.emplace_back(std::get<boost::fibers::future<vEdge>>(v).get());
            
            }
    }
    edges[0] = vv_add_fiber2(products_nofuture[0], products_nofuture[1], current_var - 1);
    edges[1] = vv_add_fiber2(products_nofuture[2], products_nofuture[3], current_var - 1);
}else{
    assert(products_nofuture.size() == 4 );
    edges[0] = vv_add_fiber2(products_nofuture[0], products_nofuture[1], current_var - 1);
    edges[1] = vv_add_fiber2(products_nofuture[2], products_nofuture[3], current_var - 1);
}


    result = makeVEdge(current_var, edges);
#ifdef CACHE
    local_mCache.set(lhs.n, rhs.n, result);
#endif
    result.w = result.w * lhs.w * rhs.w;
    if(result.w.isApproximatelyZero()){ 
        return vEdge::zero;
    }
    if(result.w.isApproximatelyOne()){
        result.w = {1.0, 0.0};
    }
    return result;    

}

vEdge mv_multiply_fiber(mEdge lhs, vEdge rhs){
    //mEdgeRefGuard mg(lhs);
    //vEdgeRefGuard vg(rhs);
    
    if(lhs.isTerminal() && rhs.isTerminal()){
        return {lhs.w * rhs.w, vNode::terminal};
    }

    // assume lhs and rhs are the same length.
    assert(lhs.getVar() == rhs.getVar());
    LIMIT = lhs.getVar()-MINUS;


    vEdge result;
#ifdef CACHE
    MulCache& local_mCache = _mCache.local();
    
    result = local_mCache.find(lhs.n, rhs.n);
    if(result.n != nullptr){
        if(result.w.isApproximatelyZero()){
            return vEdge::zero;
        }else if(result.w.isApproximatelyOne()){
            return vEdge::one;
        }else{
            result.w = result.w * lhs.w * rhs.w;
            if(result.w.isApproximatelyZero()) return vEdge::zero;
            else if (result.w.isApproximatelyOne()) return vEdge::one;
            else return result;
        }
    }
    
#endif
    result = mv_multiply_fiber2(lhs, rhs, lhs.getVar());
    return result;
}


vEdge vv_multiply(Worker* w, const vEdge& lhs, const vEdge& rhs){
    return vEdge{}; 
}




std::string measureAll(vEdge& rootEdge, const bool collapse, std::mt19937_64& mt, double epsilon) {
    std::size_t dim; 
    std_complex* vs = rootEdge.getVector(&dim);
    std::size_t max_idx = 0;
    double max_val = vs[max_idx].mag2();
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


void mEdge::decRef() {
    if (isTerminal()) return;
    assert(n->ref > 0);     
    n->ref--;
    
    for(auto i = 0; i < 4; i++)  n->getEdge(i).decRef();
    
}
void vEdge::decRef() {
    if (isTerminal()) return;

    assert(n->ref > 0);     
    n->ref--;
    
    for(auto i = 0; i < 2; i++)  n->getEdge(i).decRef();
    
}
void mEdge::incRef() {
    if (isTerminal()) return;

    if(n->ref < MAX_REF){
        n->ref++;
    }
    
    for(auto i = 0; i < 4; i++)  n->getEdge(i).incRef();
    
}
void vEdge::incRef() {
    if (isTerminal()) return;

    if(n->ref < MAX_REF){
        n->ref++;
    }
    
    for(auto i = 0; i < 2; i++)  n->getEdge(i).incRef();
    
}
void mEdge::check(){
    if(isTerminal()) return;

    if(n->v == -2 || n->v == -3 || n->v == -4) {
        foo();
        assert(false);
    }
    for(auto i = 0; i < 4; i++) n->getEdge(i).check();
    
}

mEdge RX(QubitCount qnum, int target, float angle) {
    std::complex<float> i1 = {std::cos(angle / 2), 0};
    std::complex<float> i2 = {0, -std::sin(angle / 2)};
    return makeGate(qnum, GateMatrix{i1,i2,i2,i1}, target);
}
mEdge RY(QubitCount qnum, int target, float angle) {
    std::complex<float> i1 = {std::cos(angle / 2), 0};
    std::complex<float> i2 = {-std::sin(angle / 2), 0};
    std::complex<float> i3 = {std::sin(angle / 2), 0};
    return makeGate(qnum, GateMatrix{i1,i2,i3,i1}, target);
}
mEdge RZ(QubitCount qnum, int target, float angle) {
    std::complex<float> i1 = {std::cos(angle / 2), -std::sin(angle / 2)};
    std::complex<float> i2 = {std::cos(angle / 2), std::sin(angle / 2)};
    return makeGate(qnum, GateMatrix{i1,cf_zero,cf_zero,i2}, target);
}

mEdge CX(QubitCount qnum, int target, int control){
    std::complex<float> zero = {0, 0};
    std::complex<float> one = {1, 0};
    Controls controls;
    controls.emplace(Control{control, Control::Type::pos});
    return makeGate(qnum, GateMatrix{zero,one,one,zero}, target, controls );
}
