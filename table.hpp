#pragma once

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <functional>
#include <stdio.h>
#include "common.h"
#include "dd.h"

static thread_local int64_t my_region = -1; 

/*
 * CL_MASK and CL_MASK_R are for the probe sequence calculation.
 * With 64 bytes per cacheline, there are 8 64-bit values per cacheline.
 */
static const uint64_t CL_MASK     = ~(((LINE_SIZE) / 8) - 1);  //  X & CL_MASK = X/8
static const uint64_t CL_MASK_R   = ((LINE_SIZE) / 8) - 1;    // X&CL_MASK_R = X%8
                                                              //

template<typename T, typename Hash = std::hash<T>, typename ValueEqual = std::equal_to<T>>
class CHashTable{
public:

    static constexpr size_t value_size = sizeof(T);
    static constexpr size_t key_size = 8;
    static constexpr size_t region_size = 512; //512 buckets per region

    CHashTable(size_t initial_size, size_t max_size);
    ~CHashTable(){

        free(hash);
        free(data);
        free(region_bitmap);
        free(bucket_bitmap);
    }

    size_t find_or_insert( const T& value, int* created = nullptr);
    size_t find_or_overwrite( const T& value);

    T* get_data(size_t index){
/*
        if(index >= this->table_size){
            fprintf(stderr,"get data from invalid index\n");
            exit(1);
        }
        */
        return this->data + index;
    }
    
    //for debug
    void get_index_for_value(const T& value){
        
        for(size_t i = 0 ; i < this->table_size/64; i++ ){
            uint64_t v = *(this->bucket_bitmap + i);
            if(v != 0x0000000000000000LL){
                uint64_t cursor = 0x8000000000000000LL;
                for(size_t j = 0; j < 64; j ++){
                    if( v & (cursor>>j)){
                        size_t didx = i * 64 + j;
                        if(ValueEqual()(value, *(data+didx))){
                            return;
                        }
                    }
                }
            }
    }

    printf("value does not exist!\n");

    }

private:

    void set_size(size_t size){
        
        if (size > 512 && size <= this->max_size) {
            this->table_size = size;
            /* Set threshold: number of cache lines to probe before giving up on node insertion */
            this->threshold = 192 - 2 * __builtin_clzll(this->table_size);

            this->total_region = this->table_size / region_size;
        }
    }

    size_t claim_data_bucket();

    void release_data_bucket(size_t index);

    void resize(){
        assert(false);
    }

    size_t request_region(){
        const std::lock_guard<std::mutex> lock(this->region_mtx);

        if(current_region == this->total_region) exit(1);

        return current_region++;
    }

    uint64_t          *hash;
    T                 *data;
    uint64_t          *region_bitmap;
    uint64_t          *bucket_bitmap;
    size_t            max_size;
    size_t            table_size;
    size_t            total_region;
    uint16_t          threshold;       
    
    std::mutex        region_mtx;
    size_t            current_region;

};


template<typename T, typename Hash, typename ValueEqual>
CHashTable<T, Hash, ValueEqual>::CHashTable(size_t initial_size, size_t max_size){
    

    /* Check if initial_size and max_size are powers of 2 */
    if (__builtin_popcountll(initial_size) != 1) {
        fprintf(stderr, "table_create: initial_size is not a power of 2!\n");
        exit(1);
    }

    if (__builtin_popcountll(max_size) != 1) {
        fprintf(stderr, "table_create: max_size is not a power of 2!\n");
        exit(1);
    }


    if (initial_size > max_size) {
        fprintf(stderr, "table_create: initial_size > max_size!\n");
        exit(1);
    }


    uint64_t v = 0x1FFLL;
    if(initial_size & v || max_size & v){
        fprintf(stderr, "table_create: size is not a multiple of 512!\n");
        exit(1);
    }



    this->max_size = max_size; //allocate that much?
    this->set_size(max_size);

    if(posix_memalign((void**)&hash, LINE_SIZE, max_size * key_size ) != 0){
        fprintf(stderr, "table_create: allocate memory for hash failed!\n");
        exit(1);
    }
    memset(hash, 0, max_size * key_size);
    
    if(posix_memalign((void**)&data, LINE_SIZE, max_size * value_size) != 0){
        fprintf(stderr, "table_create: allocate memory for data failed!\n");
        exit(1);
    }
    memset(data, 0, max_size*value_size);
    

    if(posix_memalign((void**)&region_bitmap, LINE_SIZE, max_size/(region_size * 8)) != 0){
        fprintf(stderr, "table_create: allocate memory for region bitmap failed!\n");
        exit(1);
    }
    memset(region_bitmap,0,max_size/(region_size * 8));

    if(posix_memalign((void**)&bucket_bitmap, LINE_SIZE, max_size / 8) != 0){
        fprintf(stderr, "table_create: allocate memory for bucket bitmap failed!\n");
        exit(1);
    }
    memset(bucket_bitmap, 0, max_size / 8);

    this->bucket_bitmap[0] = 0x8000000000000000LL; // the first bucket is for the terminal node

    
    current_region = 0;
    // initialize hashtab


}




template<typename T, typename Hash, typename ValueEqual>
void CHashTable<T, Hash, ValueEqual>::release_data_bucket(size_t index){
    if(index == 0 || index >= this->table_size){
        fprintf(stderr, "Try to release a data bucket out of valid range\n");
        exit(1);
    }

    T* d = this->data + index;
    d->~T();

    uint64_t *ptr = this->bucket_bitmap + (index/64);
    uint64_t mask = 0x8000000000000000LL >> (index&63);
    *ptr &= ~mask;
}

template<typename T, typename Hash, typename ValueEqual>
size_t CHashTable<T, Hash, ValueEqual>::claim_data_bucket(){
    

    if(my_region == -1 ){
        my_region = request_region();
    }

    int64_t start_region = my_region;
    for(;;){
       uint64_t* ptr = this->bucket_bitmap + (my_region*8);
       int i = 0;
       for(;i<8;){
            uint64_t v = *ptr;
            if(v != 0xFFFFFFFFFFFFFFFFLL){
                int j = __builtin_clzll(~v);
                *ptr |=(0x8000000000000000LL>>j);
                uint64_t idx =(8*my_region + i)*64 + j; 
                return idx; 
            }
            i++;
            ptr++;
       }
find_region:
       my_region = request_region();
       if(my_region == start_region){
        /**
         * 
         * 
         * no free region? 
         * 
         * 
         * 
        */
            resize();
       } 

       ptr = this->region_bitmap + (my_region/64);
       uint64_t mask = 0x8000000000000000LL >> (my_region&63);
       uint64_t v;
reload:
       v = *ptr;
       if(v & mask) goto find_region;
       if(cas(ptr, v, v|mask)) continue;
       else goto reload;

    }
}
    
template<typename T, typename Hash, typename ValueEqual>
    size_t CHashTable<T, Hash, ValueEqual>::find_or_insert(const T& d,  int* created)
{
    uint64_t hash_rehash = Hash()(d);

    const uint64_t step = (((hash_rehash >> 20) | 1) << 3);
    const uint64_t hash = hash_rehash & HASH_MASK;
    uint64_t idx, last = 0;
    size_t cidx = 0;

    last = idx = hash_rehash % this->table_size;

    
    for (size_t i = 0;;i ++) {
        volatile uint64_t *bucket = this->hash + idx;
        uint64_t v = *bucket;

        if (v == 0) {
            if (cidx == 0) {
                // Claim data bucket and write data
                cidx = this->claim_data_bucket();
                if(cidx == 0 || cidx >= this->table_size){
                    fprintf(stderr, "claim data bucket from an invalid index: %zu, %zu\n", cidx, total_region);
                    exit(1);
                }
                
                T *d_ptr = (this->data) + cidx;
                new(d_ptr)T(d); //call the copy ctor
            }
            if (cas(bucket, 0, hash | cidx)) {
                if(created != nullptr) * created = 1;
                return cidx;
            } else {
                v = *bucket;
            }
        }

        if (hash == (v & HASH_MASK)) {
            uint64_t d_idx = v & INDEX_MASK;
            T *d_ptr = (this->data) + d_idx;
       
            if (ValueEqual()(d, *d_ptr)) {
                
                if (cidx != 0) this->release_data_bucket(cidx);
                if(created != nullptr) *created = 0;
                return d_idx;
            }
            
        }


        // find next idx on probe sequence
        idx = (idx & CL_MASK) | ((idx+1) & CL_MASK_R);
        if (idx == last) {
            /*
            if (++i == table->threshold){
                printf("failed to find in probe sequence\n");
                return 0;
            }
            */
            // go to next cache line in probe sequence
            hash_rehash += step;

            last = idx = hash_rehash % this->table_size;

        }
    }
}


template<typename T, typename Hash, typename ValueEqual>
    size_t CHashTable<T, Hash, ValueEqual>::find_or_overwrite(const T& d)
{
    uint64_t hash_rehash = Hash()(d);

    const uint64_t hash = hash_rehash & HASH_MASK;
    assert((hash&INDEX_MASK) == 0);
    uint64_t idx;
    size_t cidx = 0;

    idx = hash_rehash % this->table_size;

    
    volatile uint64_t *bucket = this->hash + idx;
    uint64_t v = *bucket;

    if (v == 0) {
        cidx = this->claim_data_bucket();
        if(cidx == 0 || cidx >= this->table_size){
            fprintf(stderr, "claim data bucket from an invalid index: %zu, %lld\n", cidx, total_region);
            exit(1);
        }

        assert((cidx & HASH_MASK) == 0);
        assert((cidx&INDEX_MASK) == cidx);
        
        T *d_ptr = (this->data) + cidx;
        new(d_ptr)T(d); //call the copy ctor
        while (!cas(bucket, v, hash | cidx)) {
            v = * bucket;
        } 
        return cidx;
    }

    if (hash == (v & HASH_MASK)) {
        uint64_t d_idx = v & INDEX_MASK;
        T *d_ptr = (this->data) + d_idx;
   
        if (ValueEqual()(d, *d_ptr)) {
            
            if (cidx != 0) this->release_data_bucket(cidx);
            return d_idx;
        }
        
    }
    if(cidx == 0) cidx = this->claim_data_bucket();
    if(cidx == 0 || cidx >= this->table_size){
        fprintf(stderr, "claim data bucket from an invalid index: %zu\n", cidx);
        exit(1);
    }
    
    T *d_ptr = (this->data) + cidx;
    new(d_ptr)T(d); //call the copy ctor
    while(!cas(bucket, v, hash | cidx)) {
        v =*bucket;
    } 
    return cidx;
    




}

using NodeTable = CHashTable<mNode, std::hash<mNode>, compare_node_ut>;
extern NodeTable uniqueTable;






