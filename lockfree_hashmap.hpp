#pragma once
#include <cstdint>
#include <cassert>
#include <type_traits>

struct AddQuery;
struct MulQuery;

static constexpr int R = 8;

static constexpr int arrayLength = (1<<R);

static_assert(arrayLength == 256);

static constexpr int arrayFlag = (1 << 1);

static constexpr int contendFlag = (1);

static uint64_t double_to_bits(double val){
    return *reinterpret_cast<uint64_t*>(&val);
}

static constexpr int NLEVELS = sizeof(uint64_t)*8/R;
static_assert((sizeof(uint64_t)*8) % 8 == 0);

static uint8_t nth_byte(uint64_t key, int i ){
    assert(i < NLEVELS);
    key = key>>(i*NLEVELS);
    return key&0xFF;
}


template<typename T>
class LockFreeMap {
    public:

        struct DataNode{
            uint64_t hash;
            T* val;
            ~DataNode(){
                delete val;
            }
        };

        struct ArrayNode{

            ArrayNode() {
                for(int i = 0; i < arrayLength; i++) nodes[i] = reinterpret_cast<uintptr_t>(nullptr); 
            }
            uintptr_t nodes[arrayLength];

        };

        bool isArrayNode(uintptr_t p){
            return p&arrayFlag;
        }

        bool isContendNode(uintptr_t p ){
            return p&contendFlag;
        }

        bool isNull(uintptr_t p){
            return p == reinterpret_cast<uintptr_t>(nullptr);
        }

        uintptr_t markArray(uintptr_t p ){
            assert(!isArrayNode(p));
            return p|arrayFlag;
        }

        uintptr_t clearArray(uintptr_t p){
            assert(isArrayNode(p));
            return p&(~arrayFlag);
        }

        uintptr_t markContend(uintptr_t p ){
            assert(!isArrayNode(p));
            assert(!isContendNode(p));
            return p|contendFlag;
        }

        uintptr_t clearContend(uintptr_t p ){
            assert(isContendNode(p));
            return p&(~contendFlag);
        }

        ArrayNode* getArrayNode(uintptr_t p ){
            assert(isArrayNode(p));
            assert(!isContendNode(p));
            return reinterpret_cast<ArrayNode*>(clearArray(p));
        }

        DataNode* getDataNode(uintptr_t p){
            assert(!isArrayNode(p));
            assert(!isContendNode(p));
            return reinterpret_cast<DataNode*>(p);
        }

        uintptr_t getNode(uintptr_t arr, int i){
           assert(i < arrayLength); 
           ArrayNode* n = getArrayNode(arr);
           return n->nodes[i];
        }

        uintptr_t* getNodeP(uintptr_t arr, int i){
            assert(i < arrayLength);
            ArrayNode* n = getArrayNode(arr);
            return n->nodes + i;
        }

        void markDataNodeForContend(uintptr_t arr, int i){
            assert(i < arrayLength);
            ArrayNode* n = getArrayNode(arr);
            uintptr_t data_node = n->nodes[i];
            n->nodes[i] = markContend(data_node);
        }

        void expandTable(uintptr_t arr, int pos, uintptr_t old, int next_pos){
            ArrayNode* arr_node = new ArrayNode();
            assert(!isArrayNode(old));
            arr_node->nodes[next_pos] = old;
            uintptr_t arrp = markArray(reinterpret_cast<uintptr_t>(arr_node))    ;
            if(!cas(getNodeP(arr,pos), old, arrp)){
                delete arr_node;
            }
        }

        void clear_data(uintptr_t arr){
           if(isArrayNode(arr)){
                ArrayNode* n = getArrayNode(arr);
                for(int i = 0; i < arrayLength; i++) clear_data(n->nodes[i]);
                delete n;
           }else if(!isNull(arr)){
                DataNode* d = getDataNode(arr);
                delete d;
           }
        }

    public:
        LockFreeMap(){
            
            head = reinterpret_cast<uintptr_t>(new ArrayNode());
            head = markArray(head);
        }
        
        ~LockFreeMap();

        
        //require a perfect hashing function so that vals with equal hash are equal
        T* find_or_insert(uint64_t key, const T& val);

        //PHF is not required. vals with equal hash are overwritten, only the latest one is kept
        T* find_or_overwrite(uint64_t key, const T& val);

    private:
        const int maxFailCount = 10;
        uintptr_t head; // ArrayNode* 



};


template<typename T>
T* LockFreeMap<T>::find_or_insert(uint64_t key, const T& val){
   
    T* pt = new T(val);
    uintptr_t insertThis = reinterpret_cast<uintptr_t>(new DataNode{.hash = key, .val = pt});
    assert(!isContendNode(insertThis) && !isArrayNode(insertThis));
    uintptr_t local = head;

    for(int i = 0; i < NLEVELS; i++){
        uint8_t pos = nth_byte(key, i);
        int failCount = 0;

        while(true){
            uintptr_t node = getNode(local, pos);

            if(isArrayNode(node)){
                local = node;
                break;
            }else if(isNull(node)){
                if(cas(getNodeP(local, pos), reinterpret_cast<uintptr_t>(nullptr), insertThis)){
                    return pt;
                }else{
                    continue;
                }//cas failed
            }else{
                DataNode* d = getDataNode(node); 
                if(d->hash == key){
                    assert(val == *(d->val));
                    delete reinterpret_cast<DataNode*>(insertThis);
                    return d->val;
                }else{
                    assert(i < (NLEVELS - 1));
                    expandTable(local, pos, node, nth_byte(d->hash, i+1));
                    continue;
                }
            }
        }


    }//end for

    assert(false);
    
}


template<typename T>
T* LockFreeMap<T>::find_or_overwrite(uint64_t key, const T& val){
    
    T* pt = new T(val);
    uintptr_t insertThis = reinterpret_cast<uintptr_t>(new DataNode{.hash = key, .val = pt});
    assert(!isContendNode(insertThis) && !isArrayNode(insertThis));
    uintptr_t local = head;

    for(int i = 0; i < NLEVELS; i++){
        uint8_t pos = nth_byte(key, i);
        int failCount = 0;

        while(true){
            uintptr_t node = getNode(local, pos);

            if(isArrayNode(node)){
                local = node;
                break;
            }else if(isNull(node)){
                if(cas(getNodeP(local, pos), reinterpret_cast<uintptr_t>(nullptr), insertThis)){
                    return pt;
                }else{
                    continue;
                }//cas failed
            }else{
                DataNode* d = getDataNode(node); 
                if(d->hash == key){
                    if(val == *(d->val)){
                        delete reinterpret_cast<DataNode*>(insertThis);
                        return d->val; 
                    }else{
                        //overwrite
                        if(cas(getNodeP(local, pos), node, insertThis)){
                            //delete d;

                            return pt;
                        }else{
                           continue; 
                        }
                    }
                }else{
                    assert(i < (NLEVELS - 1));
                    expandTable(local, pos, node, nth_byte(d->hash, i+1));
                    continue;
                }
            }
        }


    }//end for

    assert(false);
    
}
template<typename T>
LockFreeMap<T>::~LockFreeMap(){
    clear_data(head);


}

