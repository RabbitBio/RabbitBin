/**
 * inv_index_embed.h  –  self-contained CSR inverted index for OPH sketch keys.
 *
 * Uses phmap::flat_hash_map (copied from RabbitSketch) for the posting-list
 * map, matching the implementation used by FastKMV / RabbitSketch exactly.
 *
 * Memory-efficient design
 * -----------------------
 * Callers never need to materialise all N × m keys at once.  Instead:
 *
 *  1. During sketch construction, call `insert(tid, key, sid)` for every key
 *     of every sketch as it is built.  Each thread keeps its own phmap
 *     accumulator; no synchronisation needed.
 *
 *  2. Call `build(nThreads)` once after construction is finished.  This:
 *       a. re-shards the per-thread maps into 64 shards (thread-parallel),
 *       b. merges shard-local lists (thread-parallel),
 *       c. removes singletons,
 *       d. flattens into a CSR array (postIdx + csrPosts).
 *
 *  3. In the traversal loop, call sketch[i].getKeys(scratch) on-the-fly to
 *     look up posting lists.  The per-sketch key list is never stored globally.
 *
 * Singleton pruning
 * -----------------
 * A (key, sid) pair where the key appears in only one sketch is silently
 * dropped during build().  Such keys cannot produce any cross-sketch hit
 * anyway (they only appear in one posting list), so they save both memory and
 * traversal work.
 */

#ifndef INV_INDEX_EMBED_H_
#define INV_INDEX_EMBED_H_

#include "phmap.h"

#include <omp.h>
#include <cstdint>
#include <cstring>
#include <vector>
#include <algorithm>

namespace fkmv_invidx {

struct PostRange { size_t off; uint32_t cnt; };

// ─────────────────────────────────────────────────────────────────────────────
// InvertedIndex
//   postIdx   : key → (offset, count) into csrPosts
//   csrPosts  : packed posting lists (sketch ids, uint32_t)
// ─────────────────────────────────────────────────────────────────────────────
struct InvertedIndex {
    phmap::flat_hash_map<uint64_t, PostRange> postIdx;
    std::vector<uint32_t>                      csrPosts;
    size_t                                     totalPostings = 0;
};

// ─────────────────────────────────────────────────────────────────────────────
// Builder – accumulate per-thread, then merge + flatten to CSR.
// ─────────────────────────────────────────────────────────────────────────────
struct InvertedIndexBuilder {
    using ThreadMap = phmap::flat_hash_map<uint64_t, std::vector<uint32_t>>;

    explicit InvertedIndexBuilder(int nThreads)
        : threadMaps_(nThreads) {}

    // Call from any thread without a lock (each thread owns its own map).
    inline void insert(int tid, uint64_t key, uint32_t sid) {
        threadMaps_[tid][key].push_back(sid);
    }

    // Merge per-thread maps → remove singletons → CSR flatten.
    // Must be called from a single-threaded context after all inserts are done.
    InvertedIndex build(int nThreads) {
        constexpr int NS = 64;
        const int T = (int)threadMaps_.size();

        // ── Phase 1: thread-local re-shard (no locks) ──────────────────────
        std::vector<std::vector<ThreadMap>> localShards(
            T, std::vector<ThreadMap>(NS));

        #pragma omp parallel for num_threads(nThreads) schedule(static)
        for (int tid = 0; tid < T; ++tid) {
            for (auto& [key, vec] : threadMaps_[tid]) {
                int s = (int)(key & (uint64_t)(NS - 1));
                auto& dst = localShards[tid][s];
                auto [it, ins] = dst.try_emplace(key, std::move(vec));
                if (!ins) {
                    auto& d = it->second;
                    d.insert(d.end(), vec.begin(), vec.end());
                }
            }
            ThreadMap().swap(threadMaps_[tid]);
        }
        threadMaps_.clear();

        // ── Phase 2: shard merge across threads ────────────────────────────
        std::vector<ThreadMap> invShards(NS);
        #pragma omp parallel for num_threads(nThreads) schedule(dynamic, 1)
        for (int s = 0; s < NS; ++s) {
            auto& merged = invShards[s];
            for (int tid = 0; tid < T; ++tid) {
                for (auto& [key, vec] : localShards[tid][s]) {
                    auto [it, ins] = merged.try_emplace(key, std::move(vec));
                    if (!ins) {
                        auto& d = it->second;
                        d.insert(d.end(), vec.begin(), vec.end());
                    }
                }
                ThreadMap().swap(localShards[tid][s]);
            }
        }

        // ── Phase 3: singleton removal + posting count ─────────────────────
        size_t totalAfter = 0, totalPostings = 0;
        #pragma omp parallel for num_threads(nThreads) \
            reduction(+:totalAfter, totalPostings)
        for (int s = 0; s < NS; ++s) {
            for (auto it = invShards[s].begin(); it != invShards[s].end(); ) {
                if (it->second.size() <= 1) it = invShards[s].erase(it);
                else { totalPostings += it->second.size(); ++totalAfter; ++it; }
            }
        }

        // ── Phase 4: CSR flatten ────────────────────────────────────────────
        InvertedIndex idx;
        idx.totalPostings = totalPostings;
        idx.postIdx.reserve(totalAfter * 2 + 16);
        idx.csrPosts.resize(totalPostings);

        // Pass 1 (serial): assign offsets, populate postIdx.
        std::vector<size_t> shardOff(NS, 0);
        {
            size_t cursor = 0;
            for (int s = 0; s < NS; ++s) {
                shardOff[s] = cursor;
                for (auto& [key, vec] : invShards[s]) {
                    idx.postIdx[key] = PostRange{cursor, (uint32_t)vec.size()};
                    cursor += vec.size();
                }
            }
        }

        // Pass 2 (parallel): memcpy posting lists into csrPosts.
        #pragma omp parallel for num_threads(nThreads) schedule(dynamic, 1)
        for (int s = 0; s < NS; ++s) {
            size_t off = shardOff[s];
            for (auto& [key, vec] : invShards[s]) {
                std::memcpy(&idx.csrPosts[off], vec.data(),
                            vec.size() * sizeof(uint32_t));
                off += vec.size();
                ThreadMap::value_type::second_type().swap(vec);
            }
        }

        return idx;
    }

private:
    std::vector<ThreadMap> threadMaps_;
};

} // namespace fkmv_invidx

#endif // INV_INDEX_EMBED_H_
