#include "table.hpp"
#include "cache.hpp"
#include<iostream>
#include <random>
#ifdef isMT
#include "task.h"
#endif

std::random_device seed_gen;
std::mt19937_64 mt(seed_gen());
std::uniform_real_distribution<double> dist(0.0, 1.0L);

//int GC_SIZE = 131072*2*2*2*2;

double get_random(){
    return dist(mt);
}

// int vNode_to_vec2(vNode *node, std::vector<vContent> &table,
//                  std::unordered_map<vNode *, int> &map) {
//     /*
//     This function is to serialize vNode* recursively.
//     'table' is the outcome for serialization.
//     */

//     // If table is empty, always add a terminal node as id=0.
//     // 'map' remember the processed vNode*, so you can avoid adding the same
//     // vNode* for multiple times.
//     if (table.size() == 0) {
//         vContent terminal(-1, {0.0, 0.0}, {0.0, 0.0}, 0, 0);
//         table.push_back(terminal);
//         map[node->terminal] = 0;
//     }
//     if (map.find(node) != map.end()) {
//         return map[node];
//     }

//     // If the given vNode* is not included in 'table', new data is pushed.
//     int i0 = vNode_to_vec2(node->children[0].getNode(), table, map);
//     int i1 = vNode_to_vec2(node->children[1].getNode(), table, map);
//     vContent nodeData(node->v, node->children[0].w, node->children[1].w, i0,
//                       i1);
//     table.push_back(nodeData);
//     map[node] = table.size() - 1;
//     return table.size() - 1;
// }

// vNode *vec_to_vNode2(std::vector<vContent> &table, vNodeTable &uniqTable) {
//     /*
//     This function is to de-serialize table into vNode*.
//     You can specify which uniqueTable to be used. (usually vUnique ?)
//     */

//     // The node(0) must be terminal node.
//     std::unordered_map<int, vNode *> map;
//     map[0] = &vNode::terminalNode;

//     for (int i = 1; i < table.size(); i++) {
//         vNode *node = uniqTable.getNode();
//         node->v = table[i].v;
//         vNode *i0 = map[table[i].index[0]];
//         vNode *i1 = map[table[i].index[1]];
//         vEdge e0 = {table[i].w[0], i0};
//         vEdge e1 = {table[i].w[1], i1};
//         node->children = {e0, e1};
//         node = uniqTable.lookup(node);
//         map[i] = node;
//     }
//     return map[table.size() - 1];
// }

// vEdge gc(vEdge state, std::string message=""){
//     if(vUnique.get_allocations()<GC_SIZE){
//         return state;
//     }
//     std::cout << "vSize="<<vUnique.get_allocations() << " mSize=" << mUnique.get_allocations() << " vLimit="<<GC_SIZE;

//     std::vector<vContent> v;
//     std::unordered_map<vNode *, int> map;
//     int nNodes = vNode_to_vec2(state.n, v, map);
//     if(nNodes>GC_SIZE){
//         GC_SIZE += nNodes;
//     }

//     vNodeTable new_table(NQUBITS);
//     vUnique = std::move(new_table);
//     //mNodeTable new_table_m(NQUBITS);
//     //mUnique = std::move(new_table_m);
//     state.n = vec_to_vNode2(v, vUnique);

//     AddCache newA(NQUBITS);
//     MulCache newM(NQUBITS);
//     _aCache = std::move(newA);
//     _mCache = std::move(newM);
//     std::cout << " gc_done " << message << std::endl;
//     return state;
// }

vEdge exec(int nQubits){
    int n_SingleGates = nQubits * 28;
    int n_CXGates = nQubits * 9;
    std::cout << "nQubits=" << nQubits << " Total Gates=" << n_SingleGates + n_CXGates << std::endl;

    vEdge v = makeZeroState(nQubits);

    //first_rotation
    double angle = get_random();
    for (int target = 0; target < nQubits; target++){
        auto g = RX(nQubits, target, angle);
        v = mv_multiply(g, v);
        v = gc(v);
    }
    angle = get_random();
    for (int target = 0; target < nQubits; target++){
        auto g = RZ(nQubits, target, angle);
        v = mv_multiply(g, v);
        v = gc(v);
    }
    std::cout << "First rotation" << std::endl;

    //entangler
    for (int i = 0; i < nQubits; i++){
        int control = i;
        int target = (i + 1) % nQubits;
        auto g = CX(nQubits, target, control);
        v = mv_multiply(g, v);
        v = gc(v);
    }
    std::cout << "Entangler" << std::endl;

    for (int k=0; k<8; k++){
        // mid rotation
        angle = get_random();
        for (int target = 0; target < nQubits; target++){
            auto g = RZ(nQubits, target, angle);
            v = mv_multiply(g, v);
            v = gc(v);
        }
        std::cout << "mid rotation (1) " << k << std::endl;
        angle = get_random();
        for (int target = 0; target < nQubits; target++)
        {
            auto g = RX(nQubits, target, angle);
            v = mv_multiply(g, v);
            v = gc(v);
        }
        std::cout << "mid rotation (2) " << k << std::endl;
        angle = get_random();
        for (int target = 0; target < nQubits; target++){
            auto g = RZ(nQubits, target, angle);
            v = mv_multiply(g, v);
            v = gc(v);
        }
        std::cout << "mid rotation (3) " << k << std::endl;
        //entangler
        for (int i = 0; i < nQubits; i++){
            int control = i;
            int target = (i + 1) % nQubits;
            auto g = CX(nQubits, target, control);
            v = mv_multiply(g, v);
            v = gc(v);
        }
        std::cout << "entangler " << k << std::endl;
    }
    //last rotation
    angle = get_random();
    for (int target = 0; target < nQubits; target++){
        auto g = RZ(nQubits, target, angle);
        v = mv_multiply(g, v);
        v = gc(v);
    }
    angle = get_random();
    for (int target = 0; target < nQubits; target++){
        auto g = RX(nQubits, target, angle);
        v = mv_multiply(g, v);
        v = gc(v);
    }
    return v;
}

int main(int argc, char** argv){
    auto start = std::chrono::high_resolution_clock::now();

#ifdef isMT
    Scheduler s(8);
#endif
    
    auto v = exec(std::atoi(argv[1]));

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::micro> ms = end - start;
    std::cout << "nQubits " << std::atoi(argv[1]) << " nNodes " << get_nNodes(v) << " "<< ms.count() / 1000000 << " sec"  << std::endl;
}