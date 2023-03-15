#pragma once

#include <atomic>
#include <array>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/combinable.h>


struct DeleteEntry {
    std::array<void*, 32> nodes;
    uint64_t epoch;
    std::size_t nodesCount;
    DeleteEntry *next;
};


struct DeletionList{
        DeleteEntry* headDeletionList{nullptr};
        DeleteEntry* freeLabelDeletes{nullptr};
        std::size_t deletionListCount{0};

        std::atomic<uint64_t> localEpoch;
        std::size_t thresholdCounter{0};

        ~DeletionList();
        DeleteEntry* head();

        void add(void* n, uint64_t globalEpoch);
        void remove(DeleteEntry* label, DeleteEntry* prev);

        std::size_t size();

        uint64_t deleted = 0;
        uint64_t added = 0;
};

class Epoch;
class EpochGuard;

class ThreadInfo {
    friend class Epoch;
    friend class EpochGuard;
    Epoch& epoch;
    DeletionList& deletionList;


    DeletionList& getDeletionList() const;
public:

    ThreadInfo(Epoch& epoche);

    ThreadInfo(const ThreadInfo &ti) : epoch(ti.epoch), deletionList(ti.deletionList) {
    }

    ~ThreadInfo();

    Epoch& getEpoch() const;
};

class Epoch {
    friend class ThreadInfo;
    std::atomic<uint64_t> currentEpoch{0}; //global epoch

    tbb::enumerable_thread_specific<DeletionList> deletionLists;

    size_t startGCThreshhold;


public:
    Epoch(size_t startGCThreshhold) : startGCThreshhold(startGCThreshhold) { }

    ~Epoch();

    void enterEpoch(ThreadInfo &epochInfo);

    void markNodeForDeletion(void *n, ThreadInfo &epochInfo);

    void exitEpochAndCleanup(ThreadInfo &info);

    void showDeleteRatio();

};

class EpochGuard {
    ThreadInfo& threadEpochInfo;
public:

    EpochGuard(ThreadInfo& threadEpochInfo) : threadEpochInfo(threadEpochInfo) {
        threadEpochInfo.getEpoch().enterEpoch(threadEpochInfo);
    }

    ~EpochGuard() {
        threadEpochInfo.getEpoch().exitEpochAndCleanup(threadEpochInfo);
    }
};

class EpochGuardReadonly {
public:

    EpochGuardReadonly(ThreadInfo& threadEpochInfo) {
        threadEpochInfo.getEpoch().enterEpoch(threadEpochInfo);
    }

    ~EpochGuardReadonly() {
    }
};

inline ThreadInfo::~ThreadInfo() {
    deletionList.localEpoch.store(std::numeric_limits<uint64_t>::max());
}
