#include "dd.h"
#include "common.h"
#include "engine.h"
#include <algorithm>
#include "table.hpp"

#define SUBTASK_THRESHOLD 5

mNodeTable mUnique(20);
vNodeTable vUnique(20);

std::vector<mEdge> identityTable(20);
AddTable addTable;
MulTable mulTable;


mEdge mEdge::one{{1.0, 0.0}, mNode::terminal};
mEdge mEdge::zero{{0.0,0.0}, mNode::terminal};
vEdge vEdge::one{{1.0, 0.0}, vNode::terminal};
vEdge vEdge::zero{{0.0,0.0}, vNode::terminal};

mNode mNode::terminalNode{.v = -1, .children = {}, .next = nullptr};
vNode vNode::terminalNode{.v = -1, .children = {}, .next = nullptr};

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

static vEdge normalize(Worker* w,  const vEdge& e){

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

vEdge makeEdge(Worker* w, Qubit q, const std::array<vEdge, 2>& c){
    

    

    vNode* node = vUnique.getNode();
    node->v = q;
    node->children = c;

    for(int i = 0; i < 2; i++){
        assert(&node->children[i] != &vEdge::one && (&node->children[i] != &vEdge::zero));
    }

    
    vEdge e =  normalize(w, {{1.0,0.0}, node}); 

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
    if(identityTable[q].n != nullptr) {
        return identityTable[q];
    }

    mEdge e = makeEdge(w, 0, {mEdge::one, mEdge::zero, mEdge::zero, mEdge::one});
    for(Qubit i = 1; i <= q; i++){
       e = makeEdge(w, i, {{e,mEdge::zero,mEdge::zero,e}}); 
    }
    
    identityTable[q] = e;
    return e;


}

vEdge makeZeroState(Worker *w, QubitCount q){
    vEdge e = makeEdge(w, 0, {vEdge::one, vEdge::zero});
    for(Qubit i = 1; i < q; i++){
       e = makeEdge(w, i, {{e,vEdge::zero}}); 
    }
    return e;
}

vEdge makeOneState(Worker *w, QubitCount q){
    vEdge e = makeEdge(w, 0, {vEdge::zero, vEdge::one});
    for(Qubit i = 1; i < q; i++){
       e = makeEdge(w, i, {{vEdge::zero, e}}); 
    }
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

mEdge add2(Worker* w, const mEdge& lhs, const mEdge& rhs, int32_t current_var);
static mEdge sum_for_entry(Worker* w, const mEdge& lhs, const mEdge& rhs, int i, int32_t current_var){


    mEdge x, y;

    Qubit lv = lhs.getVar();
    Qubit rv = rhs.getVar();
    mNode* lnode = lhs.getNode();
    mNode* rnode = rhs.getNode();
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

    mEdge result = add2(w, x, y, current_var - 1);
    return result;

}
mEdge add2(Worker* w, const mEdge& lhs, const mEdge& rhs, int32_t current_var){
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
        return result;
    }


    std::array<mEdge, 4> edges;

    Job* jobs[2];

    for(auto i = 0; i < 4; i++){

       if(i < 2 && current_var >= SUBTASK_THRESHOLD) jobs[i] = w->submit(sum_for_entry, lhs, rhs, i, current_var ); 
       else{
        edges[i] = sum_for_entry(w, lhs, rhs, i, current_var );
       }
    }

    if(current_var >= SUBTASK_THRESHOLD){
        while(!jobs[0]->available()){
            w->run_pending();
        }
        while(!jobs[1]->available()){
            w->run_pending();
        }

        edges[0] = jobs[0]->getResult();
        edges[1] = jobs[1]->getResult();
    
    }



    result =  makeEdge(w, current_var, edges);
    w->_addCache.set(lhs, rhs, result);
    

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

mEdge multiply2(Worker* w, const mEdge& lhs, const mEdge& rhs, int32_t current_var);

static mEdge product_for_entry(Worker* w, const mEdge& lhs, const mEdge& rhs, int i, int32_t current_var){
    Qubit lv = lhs.getVar();
    Qubit rv = rhs.getVar();
    mNode* lnode = lhs.getNode();
    mNode* rnode = rhs.getNode();
    mEdge x,y;

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

        
        product[k] = multiply2(w, x, y, current_var - 1); 
        
    }

    mEdge result = add2(w, product[0], product[1], current_var - 1);
    return result;

}

mEdge multiply2(Worker* w, const mEdge& lhs, const mEdge& rhs, int32_t current_var){

    if(lhs.w.isZero() || rhs.w.isZero()){
        return mEdge::zero;
    }

    if(current_var == -1) {
        assert(lhs.isTerminal() && rhs.isTerminal());
        return {lhs.w * rhs.w, mNode::terminal};
    }
    
    mEdge result = w->_mulCache.find(lhs.n, rhs.n);
    if(result.n != nullptr){
        return {result.w * lhs.w * rhs.w, result.n};
    }
    

    mEdge lcopy = lhs;
    mEdge rcopy = rhs;
    lcopy.w = {1.0,0.0};
    rcopy.w = {1.0, 0.0};


    std::array<mEdge, 4 > edges;

    Job* jobs[2];
    for(auto i = 0; i < 4; i++){
        if(i < 2 && current_var >= SUBTASK_THRESHOLD) jobs[i] = w->submit(product_for_entry, lcopy, rcopy, i, current_var);
        else{
            edges[i] = product_for_entry(w, lcopy, rcopy, i, current_var);
        }

    }
    if(current_var >= SUBTASK_THRESHOLD){
        while(!jobs[0]->available()){
            w->run_pending();
        }
        while(!jobs[1]->available()){
            w->run_pending();
        }

        edges[0] = jobs[0]->getResult();
        edges[1] = jobs[1]->getResult();
    
    }

    result = makeEdge(w, current_var, edges);
    w->_mulCache.set(lhs.n, rhs.n, result);

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
struct kjob{
    mEdge l;
    mEdge r;
    int idx;
};
/*
mEdge kronecker2(Worker *w, const mEdge &lhs, const mEdge &rhs)
{
    std::vector<kjob> jobs;
    kjob initial_job = {lhs, rhs, 0};
    jobs.push_back(initial_job);
    std::vector<mEdge> result;

    while(jobs.size()){
        kjob& job = jobs.back();
        auto ltmp = job.l;
        auto rtmp = job.r;
        char idx = job.idx;
        if (idx == 4)
        {
            jobs.pop_back();
            Qubit lv = ltmp.getVar();
            Qubit rv = rtmp.getVar();
            std::array<mEdge, 4> last4;
            assert(result.size() > 3);
            copy(result.end() - 4, result.end(), last4.begin());
            result.erase(result.end() - 4, result.end());
            mEdge ret = makeEdge(w,  lv + rv + 1, last4);
            result.push_back(ret);
        }else{
            if (ltmp.isTerminal())
            {
                result.push_back({ltmp.w * rtmp.w, rtmp.n});
                jobs.pop_back();
            }else{
                job.idx += 1;
                auto x = ltmp.getNode()->getEdge(idx);
                x.w = ltmp.w * x.w;
                kjob next_job = {x, rhs, 0};
                jobs.push_back(next_job);
            }
        }
    }
    assert(result.size() == 1);
    return result[0];
}
*/

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

vEdge add2(Worker* w, const vEdge& lhs, const vEdge& rhs, int32_t current_var);
static vEdge sum_for_entry(Worker* w, const vEdge& lhs, const vEdge& rhs, int i, int32_t current_var){


    vEdge x, y;

    Qubit lv = lhs.getVar();
    Qubit rv = rhs.getVar();
    vNode* lnode = lhs.getNode();
    vNode* rnode = rhs.getNode();
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

    vEdge result = add2(w, x, y, current_var - 1);
    return result;

}

vEdge add2(Worker* w, const vEdge& lhs, const vEdge& rhs, int32_t current_var){
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


    std::array<vEdge, 2> edges;

    for(auto i = 0; i < 4; i++){
        edges[i] = sum_for_entry(w, lhs, rhs, i, current_var );
    }

    result =  makeEdge(w, current_var, edges);
    w->_addCache.set(lhs, rhs, result);
    

    return result;


}

vEdge add(Worker* w, const vEdge& lhs, const vEdge& rhs){
    if(lhs.isTerminal() && rhs.isTerminal()){
        return {lhs.w + rhs.w, vNode::terminal};
    }

    // assume lhs and rhs are the same size vector.
    assert(lhs.getVar()==rhs.getVar());
    return add2(w, lhs, rhs, lhs.getVar());
}

vEdge kronecker2(Worker* w, const vEdge& lhs, const vEdge& rhs){
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
        edges[i] = kronecker2(w, x, rhs);
    }

    vEdge ret = makeEdge(w,  lv + rv + 1, edges);
    return ret;
}

vEdge kronecker(Worker* w, const vEdge& lhs, const vEdge& rhs){
    if(lhs.isTerminal() && rhs.isTerminal()){
        return {lhs.w * rhs.w, vNode::terminal};
   }

   return kronecker2(w, lhs, rhs);
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

vEdge mEdge::get_column(std::size_t col)  const {

}
