#pragma once
#include "allocator/CacheAllocator.h"
#include "allocator/CacheItem.h"
#include "allocator/memory/Slab.h"
#include "common.h"
#include "dd.h"
#include <cstdint>
#include <cstring>
#include <memory>
#include <string>

class MulOpCache {

public:
  MulOpCache(QubitCount q) {
    Cache::Config config;
    config.setCacheSize(1 * 1024 * 1024 * 1024)
        .setCacheName("MulOPCache")
        .validate();
    cache = std::make_unique<Cache>(config);
    default_pool =
        cache->addPool("default", cache->getCacheMemoryStats().ramCacheSize);
  }

  template <typename LT, typename RT,
            typename RetT =
                std::conditional_t<std::is_same_v<RT, vNode>, vEdge, mEdge>>
  RetT find(const LT *lhs, const RT *rhs) {
    auto key = std::to_string(murmur_hash(reinterpret_cast<uintptr_t>(lhs)));
    CacheReadHandle handle = cache->find(key);
    if (handle) {
      auto data = reinterpret_cast<const RetT *>(handle->getMemory());
      return *data;
    } else {
      return RetT{};
    }
  }

  template <typename LT, typename RT,
            typename ResT =
                std::conditional_t<std::is_same_v<RT, vNode>, vEdge, mEdge>>
  void set(const LT *lhs, const RT *rhs, const ResT &result) {
    auto key = std::to_string(murmur_hash(reinterpret_cast<uintptr_t>(lhs)));
    CacheWriteHandle handle = cache->allocate(default_pool, key, sizeof(ResT));
    if (handle) {
      void *pwm = handle->getMemory();
      memcpy(pwm, &result, sizeof(ResT));
      cache->insertOrReplace(handle);
    }
  }

private:
  using Cache = facebook::cachelib::LruAllocator;
  using PoolId = facebook::cachelib::PoolId;
  using CacheWriteHandle = typename Cache::WriteHandle;
  using CacheReadHandle = typename Cache::ReadHandle;

  std::unique_ptr<Cache> cache;
  PoolId default_pool;

  template <typename T> std::size_t hash(const T &lhs, const T &rhs) {
    auto h1 = std::hash<T>()(lhs);
    auto h2 = std::hash<T>()(rhs);
    return hash_combine(h1, h2);
  }
};

class AddOpCache {

public:
  AddOpCache(QubitCount q) {
    Cache::Config config;
    config.setCacheSize(1 * 1024 * 1024 * 1024)
        .setCacheName("AddOPCache")
        .validate();
    cache = std::make_unique<Cache>(config);
    default_pool =
        cache->addPool("default", cache->getCacheMemoryStats().ramCacheSize);
  }

  template <typename T> T find(T lhs, T rhs) {
    auto key = std::to_string(this->hash(lhs, rhs));
    CacheReadHandle handle = cache->find(key);
    if (handle) {
      auto data = reinterpret_cast<const T *>(handle->getMemory());
      return *data;
    } else {
      return T{};
    }
  }

  template <typename T> void set(T lhs, T rhs, const T &result) {
    auto key = this->hash(lhs, rhs);
    CacheWriteHandle handle =
        cache->allocate(default_pool, std::to_string(key), sizeof(T));
    if (handle) {
      void *pwm = handle->getMemory();
      memcpy(pwm, &result, sizeof(T));
      cache->insertOrReplace(handle);
    }
  }

private:
  using Cache = facebook::cachelib::LruAllocator;
  using PoolId = facebook::cachelib::PoolId;
  using CacheWriteHandle = typename Cache::WriteHandle;
  using CacheReadHandle = typename Cache::ReadHandle;

  std::unique_ptr<Cache> cache;
  PoolId default_pool;

  template <typename T> std::size_t hash(const T &lhs, const T &rhs) {
    auto h1 = std::hash<T>()(lhs);
    auto h2 = std::hash<T>()(rhs);
    return hash_combine(h1, h2);
  }
};
