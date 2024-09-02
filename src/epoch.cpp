#include "epoch.h"

#include <cassert>
#include <iostream>

DeletionList::~DeletionList() {
    assert(deletionListCount == 0 && headDeletionList == nullptr);
    DeleteEntry* cur = nullptr;
    DeleteEntry* next = freeLabelDeletes;
    while (next != nullptr) {
        cur = next;
        next = cur->next;
        delete cur;
    }
    freeLabelDeletes = nullptr;
}

std::size_t DeletionList::size() { return deletionListCount; }

void DeletionList::remove(DeleteEntry* label, DeleteEntry* prev) {
    if (prev == nullptr) {
        headDeletionList = label->next;
    } else {
        prev->next = label->next;
    }
    deletionListCount -= label->nodesCount;

    label->next = freeLabelDeletes;
    freeLabelDeletes = label;
    deleted += label->nodesCount;
}
void DeletionList::add(void* n, uint64_t globalEpoch) {
    deletionListCount++;
    DeleteEntry* label;
    if (headDeletionList != nullptr &&
        headDeletionList->nodesCount < headDeletionList->nodes.size()) {
        label = headDeletionList;
    } else {
        if (freeLabelDeletes != nullptr) {
            label = freeLabelDeletes;
            freeLabelDeletes = freeLabelDeletes->next;
        } else {
            label = new DeleteEntry();
        }
        label->nodesCount = 0;
        label->next = headDeletionList;
        headDeletionList = label;
    }
    label->nodes[label->nodesCount] = n;
    label->nodesCount++;
    label->epoch = globalEpoch;

    added++;
}

DeleteEntry* DeletionList::head() { return headDeletionList; }

void Epoch::enterEpoch(ThreadInfo& epochInfo) {
    unsigned long curEpoch = currentEpoch.load(std::memory_order_relaxed);
    epochInfo.getDeletionList().localEpoch.store(curEpoch,
                                                 std::memory_order_release);
}

void Epoch::markNodeForDeletion(void* n, ThreadInfo& epochInfo) {
    epochInfo.getDeletionList().add(n, currentEpoch.load());
    epochInfo.getDeletionList().thresholdCounter++;
}

void Epoch::exitEpochAndCleanup(ThreadInfo& epochInfo) {
    DeletionList& deletionList = epochInfo.getDeletionList();
    if ((deletionList.thresholdCounter & (64 - 1)) == 1) {
        currentEpoch++;
    }
    if (deletionList.thresholdCounter > startGCThreshhold) {
        if (deletionList.size() == 0) {
            deletionList.thresholdCounter = 0;
            return;
        }
        deletionList.localEpoch.store(std::numeric_limits<uint64_t>::max());

        uint64_t oldestEpoch = std::numeric_limits<uint64_t>::max();
        for (auto& epoch : deletionLists) {
            auto e = epoch.localEpoch.load();
            if (e < oldestEpoch) {
                oldestEpoch = e;
            }
        }

        DeleteEntry *cur = deletionList.head(), *next, *prev = nullptr;
        while (cur != nullptr) {
            next = cur->next;

            if (cur->epoch < oldestEpoch) {
                for (std::size_t i = 0; i < cur->nodesCount; ++i) {
                    delete (cur->nodes[i]);
                }
                deletionList.remove(cur, prev);
            } else {
                prev = cur;
            }
            cur = next;
        }
        deletionList.thresholdCounter = 0;
    }
}

Epoch::~Epoch() {
    uint64_t oldestEpoch = std::numeric_limits<uint64_t>::max();
    for (auto& epoch : deletionLists) {
        auto e = epoch.localEpoch.load();
        if (e < oldestEpoch) {
            oldestEpoch = e;
        }
    }
    for (auto& d : deletionLists) {
        DeleteEntry *cur = d.head(), *next, *prev = nullptr;
        while (cur != nullptr) {
            next = cur->next;

            assert(cur->epoch < oldestEpoch);
            for (std::size_t i = 0; i < cur->nodesCount; ++i) {
                delete cur->nodes[i];
            }
            d.remove(cur, prev);
            cur = next;
        }
    }
}

void Epoch::showDeleteRatio() {
    for (auto& d : deletionLists) {
        std::cout << "deleted " << d.deleted << " of " << d.added << std::endl;
    }
}

ThreadInfo::ThreadInfo(Epoch& epoch)
    : epoch(epoch), deletionList(epoch.deletionLists.local()) {}

DeletionList& ThreadInfo::getDeletionList() const { return deletionList; }

Epoch& ThreadInfo::getEpoch() const { return epoch; }
