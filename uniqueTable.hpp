#pragma once

#include <functional>
#include "table.hpp"

template<typename T, typename Hash = std::hash<T>, typename ValueEqual = std::equal_to<T>>
class UniqueTable {
    public:

    private:
        std::vector<CHashTable<T, Hash, ValueEqual>> _tables;



};
