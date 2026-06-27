/**
 * rabbit_sketch.h  –  self-contained single-header k-mer sketch
 *
 * Representation: One-Permutation MinHash (OPH) with optimal densification,
 * projected to a 1-bit-per-bucket packed signature.
 *
 * Rationale (vs. the original bottom-k / KMV merge):
 *   - Sketch build  : single O(1) min-update per k-mer (no sorted insert/memmove).
 *   - Jaccard       : pure SIMD XOR + popcount over a tiny m/8-byte bitset,
 *                     branchless and extremely cache-friendly on the hot
 *                     all-pairs loop.
 *   - Centroid merge: element-wise min over the registers.
 *
 * Public API:
 *   rabbit_sketch::KmerSketch sk(num_buckets, kmer_size);
 *   sk.update(seq_cstr, length);
 *   double j = sk.jaccard(other);      // Jaccard in [0,1]
 *   KmerSketch merged = sk.merge(other);  // centroid sketch
 *
 * No external dependencies beyond the C++ standard library and SIMD headers.
 */

#ifndef RABBIT_SKETCH_H_
#define RABBIT_SKETCH_H_

#include <cstdint>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <climits>
#include <limits>
#include <memory>
#include <algorithm>
#include <vector>
#include <immintrin.h>

namespace rabbit_sketch {

// ═══════════════════════════════════════════════════════════════════════════
// murmur3_fmix  (from hash_int.h, MetaCache, André Müller, GPL3)
// ═══════════════════════════════════════════════════════════════════════════
static inline uint64_t murmur3_fmix(uint64_t x, uint64_t seed) noexcept {
    x ^= seed;
    x ^= x >> 33;
    x *= UINT64_C(0xff51afd7ed558ccd);
    x ^= x >> 33;
    x *= UINT64_C(0xc4ceb9fe1a85ec53);
    x ^= x >> 33;
    return x;
}

// ═══════════════════════════════════════════════════════════════════════════
// 2-bit lex encoding LUT  (A/a→0, C/c→1, G/g→2, T/t→3; else→255)
// ═══════════════════════════════════════════════════════════════════════════
static const uint8_t RB_BASE_ENC_LUT[256] = {
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,  0,255,  1,255,255,255,  2,255,255,255,255,255,255,255,255,
    255,255,255,255,  3,255,255,255,255,255,255,255,255,255,255,255,
    255,  0,255,  1,255,255,255,  2,255,255,255,255,255,255,255,255,
    255,255,255,255,  3,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
};
#define RB_BASE_ENC(c)   (RB_BASE_ENC_LUT[(uint8_t)(c)])
#define RB_BASE_COMP(e)  ((uint8_t)((e) <= 3 ? 3-(e) : 255))
#define RB_BASE_VALID(e) ((e) <= 3)

// ═══════════════════════════════════════════════════════════════════════════
// KmerSketch  –  one-permutation MinHash with 1-bit packed signature
// ═══════════════════════════════════════════════════════════════════════════
class KmerSketch {
public:
    static constexpr uint32_t REG_EMPTY = UINT32_MAX;

    // `k` is reinterpreted as the number of OPH buckets. `b` is the number of
    // packed bits stored per bucket (bit-planes); more bits => lower variance
    // / lower collision baseline at the cost of a larger signature.
    explicit KmerSketch(uint32_t k = 1000, int kmer_size = 21,
                        uint32_t b = 1, uint64_t seed = 42)
        : m_(k < 2 ? 2 : k), kmer_size_(kmer_size),
          b_(b < 1 ? 1 : (b > 24 ? 24 : b)), seed_(seed),
          nwords_((m_ + 63) / 64),
          reg_(new uint32_t[m_]),
          bits_(new uint64_t[(size_t)nwords_ * b_]),
          nonempty_(0),
          sig_valid_(false)
    {
        assert(kmer_size >= 1 && kmer_size <= 32);
        std::fill(reg_.get(), reg_.get() + m_, REG_EMPTY);
    }

    KmerSketch(const KmerSketch& o)
        : m_(o.m_), kmer_size_(o.kmer_size_), b_(o.b_), seed_(o.seed_),
          nwords_(o.nwords_),
          reg_(new uint32_t[o.m_]),
          bits_(new uint64_t[(size_t)o.nwords_ * o.b_]),
          nonempty_(o.nonempty_),
          sig_valid_(o.sig_valid_)
    {
        std::copy(o.reg_.get(),  o.reg_.get()  + m_, reg_.get());
        std::copy(o.bits_.get(), o.bits_.get() + (size_t)nwords_ * b_, bits_.get());
    }

    KmerSketch& operator=(KmerSketch other) {
        std::swap(m_,         other.m_);
        std::swap(kmer_size_, other.kmer_size_);
        std::swap(b_,         other.b_);
        std::swap(seed_,      other.seed_);
        std::swap(nwords_,    other.nwords_);
        std::swap(reg_,       other.reg_);
        std::swap(bits_,      other.bits_);
        std::swap(nonempty_,  other.nonempty_);
        std::swap(sig_valid_, other.sig_valid_);
        return *this;
    }

    KmerSketch(KmerSketch&&) = default;
    // operator=(KmerSketch) covers both copy and move (copy-and-swap idiom)
    ~KmerSketch() = default;

    // ── sketch building ──────────────────────────────────────────────────
    void update(const char* seq, uint64_t length);
    void finalize() {}  // no-op

    // ── similarity / distance ────────────────────────────────────────────
    double jaccard(const KmerSketch& other, double min_jaccard = 0.0) const;
    double distance(const KmerSketch& other,
                    double max_distance = std::numeric_limits<double>::infinity()) const;
    double ani(const KmerSketch& other) const;
    double containment(const KmerSketch& other) const;
    double cardinality() const;

    // ── merge (element-wise register min = union) ────────────────────────
    KmerSketch merge(const KmerSketch& other) const;

    // ── accessors ────────────────────────────────────────────────────────
    const uint64_t* getSignature() const { ensureSig(); return bits_.get(); }
    uint32_t getK()        const { return m_; }
    int      getKmerSize() const { return kmer_size_; }
    uint32_t size()        const { return nonempty_; }

    // ── memory management ────────────────────────────────────────────────
    // Free the raw register array after the b-bit signature has been built
    // (i.e. after buildSig() / getSignature() has been called). The signature
    // in bits_[] is sufficient for all subsequent jaccard() calls, so the
    // reg_[] array (m_*4 bytes) can be released to reduce peak RSS.
    // NOTE: Do NOT call this before merge() or getKeys() are finished.
    void freeRegisters() { reg_.reset(); }

    // ── inverted-index keys ──────────────────────────────────────────────
    // One hashed key per non-empty bucket: a (bucket, min-value) pair. Two
    // sketches share a key iff they have the same MinHash min in the same
    // bucket (= a bucket collision), so counting shared keys between two
    // sketches directly yields their bucket-match count (Jaccard numerator).
    // The murmur mix spreads the small (bucket<<32|val) range across 64 bits
    // so the inverted-index high-bit sharding distributes evenly.
    void getKeys(std::vector<uint64_t>& out) const {
        out.clear();
        out.reserve(nonempty_);
        const uint32_t* reg = reg_.get();
        for (uint32_t i = 0; i < m_; ++i) {
            if (reg[i] == REG_EMPTY) continue;
            const uint64_t packed = ((uint64_t)i << 32) | (uint64_t)reg[i];
            out.push_back(murmur3_fmix(packed, seed_ ^ UINT64_C(0x51A3C2B7)));
        }
    }

private:
    void     ensureSig() const { if (!sig_valid_) buildSig(); }
    void     buildSig()  const;
    uint32_t densify(uint32_t i) const;  // returns a borrowed register value

    uint32_t m_;          // number of buckets
    int      kmer_size_;
    uint32_t b_;          // bits per bucket (number of bit-planes)
    uint64_t seed_;
    uint32_t nwords_;     // (m_ + 63) / 64
    std::unique_ptr<uint32_t[]> reg_;    // OPH min registers
    mutable std::unique_ptr<uint64_t[]> bits_;  // b_ packed bit-planes, plane p at [p*nwords_]
    mutable uint32_t nonempty_;
    mutable bool     sig_valid_;
};

// ─── reduce a 64-bit hash to [0, m) via Lemire's multiply-shift ──────────
static inline uint32_t rb_sketch_reduce(uint64_t h, uint32_t m) {
    return (uint32_t)(((h >> 32) * (uint64_t)m) >> 32);
}

// ─── update  (2-bit lex rolling + murmur3_fmix + scatter min into bucket) ─
inline void KmerSketch::update(const char* seq, uint64_t length) {
    const int K = kmer_size_;
    if (length < static_cast<uint64_t>(K)) return;

    const uint64_t loc_seed = seed_;
    const uint64_t kmer_mask = (K >= 32) ? UINT64_MAX : ((1ULL << (2*K)) - 1);

    uint64_t fwd_h = 0, rc_h = 0;
    int inv = 0;
    for (int k = 0; k < K; ++k) {
        uint8_t ef = RB_BASE_ENC(seq[k]);
        if (!RB_BASE_VALID(ef)) ++inv;
        fwd_h = ((fwd_h << 2) | (RB_BASE_VALID(ef) ? (ef & 3u) : 0u)) & kmer_mask;
        uint8_t er = RB_BASE_VALID(ef) ? (RB_BASE_COMP(ef) & 3u) : 0u;
        rc_h  = (rc_h >> 2) | (static_cast<uint64_t>(er) << (2*(K-1)));
    }

    const uint64_t N_body = length - static_cast<uint64_t>(K);
    uint32_t* __restrict__ reg = reg_.get();

    for (uint64_t i = 0; ; ++i) {
        if (inv == 0) {
            uint64_t canon = (fwd_h < rc_h) ? fwd_h : rc_h;
            uint64_t h     = murmur3_fmix(canon, loc_seed);
            uint32_t b     = rb_sketch_reduce(h, m_);
            uint32_t v     = (uint32_t)h;
            if (v < reg[b]) {
                if (reg[b] == REG_EMPTY) ++nonempty_;
                reg[b] = v;
            }
        }
        if (i >= N_body) break;
        uint8_t ef_out = RB_BASE_ENC(seq[i]);
        uint8_t ef_in  = RB_BASE_ENC(seq[i + K]);
        if (!RB_BASE_VALID(ef_out)) --inv;
        if (!RB_BASE_VALID(ef_in))  ++inv;
        fwd_h = ((fwd_h << 2) | (RB_BASE_VALID(ef_in) ? (ef_in & 3u) : 0u)) & kmer_mask;
        uint8_t er_in = RB_BASE_VALID(ef_in) ? (RB_BASE_COMP(ef_in) & 3u) : 0u;
        rc_h  = (rc_h >> 2) | (static_cast<uint64_t>(er_in) << (2*(K-1)));
    }

    sig_valid_ = false;
}

// ─── densify  (optimal densification for an empty bucket) ─────────────────
// Deterministic probe sequence shared by all sketches, so two sketches that
// are both empty at bucket i borrow consistently. The hop count is mixed into
// the returned value to decorrelate distinct empty buckets.
inline uint32_t KmerSketch::densify(uint32_t i) const {
    const uint32_t* reg = reg_.get();
    uint64_t base = murmur3_fmix(i, seed_ ^ UINT64_C(0xA5A5A5A5DEADBEEF));
    uint64_t step = murmur3_fmix(i, seed_ ^ UINT64_C(0x123456789ABCDEF0)) | 1ULL;
    uint64_t cur  = base;
    for (uint32_t t = 1; t <= m_; ++t) {
        cur += step;
        uint32_t cand = rb_sketch_reduce(cur, m_);
        uint32_t v = reg[cand];
        if (v != REG_EMPTY)
            return v ^ (uint32_t)murmur3_fmix(t, UINT64_C(0x9E3779B97F4A7C15));
    }
    return 0;  // unreachable when nonempty_ > 0
}

// ─── buildSig  (densify + pack low b_ bits of each register into planes) ──
inline void KmerSketch::buildSig() const {
    uint64_t* b = bits_.get();
    std::fill(b, b + (size_t)nwords_ * b_, 0ULL);
    if (nonempty_ != 0) {
        const uint32_t* reg = reg_.get();
        for (uint32_t i = 0; i < m_; ++i) {
            uint32_t v = reg[i];
            if (v == REG_EMPTY) v = densify(i);
            const uint32_t word = i >> 6;
            const uint64_t bit  = 1ULL << (i & 63);
            for (uint32_t p = 0; p < b_; ++p)
                if ((v >> p) & 1u) b[(size_t)p * nwords_ + word] |= bit;
        }
    }
    sig_valid_ = true;
}

// ─── jaccard  (1-bit MinHash estimator via SIMD popcount) ─────────────────
inline double KmerSketch::jaccard(const KmerSketch& other, double min_jaccard) const {
    assert(m_ == other.m_ && b_ == other.b_);
    if (nonempty_ == 0 || other.nonempty_ == 0) return 0.0;
    ensureSig(); other.ensureSig();

    const uint64_t* __restrict__ a = bits_.get();
    const uint64_t* __restrict__ b = other.bits_.get();
    const uint32_t  nw   = nwords_;
    const uint32_t  np   = b_;
    const uint32_t  full = (m_ % 64 == 0) ? nw : (nw - 1);
    const uint32_t  rem  = m_ & 63;
    const uint64_t  lastmask = (rem == 0) ? ~0ULL : ((1ULL << rem) - 1);

    // A bucket matches iff all b_ bit-planes agree there: AND of XNOR(a,b).
    // SIMD: 8 × uint64 popcount per step (mirrors rabbitbin.cpp::jaccard_raw
    // and the BinDash XNOR+VPOPCNT kernel in RabbitSketch).  ternarylogic
    // imm 0xC3 = ~(A ^ B) (XNOR) computes the per-plane agreement in one op.
    uint64_t matches = 0;
    uint32_t w = 0;
#if defined(__AVX512F__) && defined(__AVX512VPOPCNTDQ__)
    {
        __m512i vsum = _mm512_setzero_si512();
        if (np == 1) {
            for (; w + 8 <= full; w += 8) {
                __m512i va = _mm512_loadu_si512((const __m512i*)(a + w));
                __m512i vb = _mm512_loadu_si512((const __m512i*)(b + w));
                __m512i vm = _mm512_ternarylogic_epi64(va, vb, va, 0xC3);
                vsum = _mm512_add_epi64(vsum, _mm512_popcnt_epi64(vm));
            }
        } else if (np == 2) {
            const uint64_t* a1 = a + nw;
            const uint64_t* b1 = b + nw;
            for (; w + 8 <= full; w += 8) {
                __m512i va0 = _mm512_loadu_si512((const __m512i*)(a  + w));
                __m512i vb0 = _mm512_loadu_si512((const __m512i*)(b  + w));
                __m512i va1 = _mm512_loadu_si512((const __m512i*)(a1 + w));
                __m512i vb1 = _mm512_loadu_si512((const __m512i*)(b1 + w));
                __m512i xn0 = _mm512_ternarylogic_epi64(va0, vb0, va0, 0xC3);
                __m512i xn1 = _mm512_ternarylogic_epi64(va1, vb1, va1, 0xC3);
                vsum = _mm512_add_epi64(vsum,
                          _mm512_popcnt_epi64(_mm512_and_si512(xn0, xn1)));
            }
        } else {
            for (; w + 8 <= full; w += 8) {
                __m512i vm = _mm512_set1_epi64(-1LL);
                for (uint32_t p = 0; p < np; ++p) {
                    __m512i va = _mm512_loadu_si512((const __m512i*)(a + (size_t)p*nw + w));
                    __m512i vb = _mm512_loadu_si512((const __m512i*)(b + (size_t)p*nw + w));
                    vm = _mm512_and_si512(vm,
                            _mm512_ternarylogic_epi64(va, vb, va, 0xC3));
                }
                vsum = _mm512_add_epi64(vsum, _mm512_popcnt_epi64(vm));
            }
        }
        matches = (uint64_t)_mm512_reduce_add_epi64(vsum);
    }
#endif
    // Scalar tail (remaining full words) + last partial word.
    for (; w < nw; ++w) {
        uint64_t m = ~0ULL;
        for (uint32_t p = 0; p < np; ++p) {
            const size_t off = (size_t)p * nw + w;
            m &= ~(a[off] ^ b[off]);
        }
        if (w >= full) m &= lastmask;
        matches += (uint64_t)__builtin_popcountll(m);
    }

    // b-bit MinHash estimator: P(all bits match) = J + (1-J) * 2^{-b}.
    const double base = std::ldexp(1.0, -(int)np);   // 2^{-b}
    double p = (double)matches / (double)m_;
    double j = (p - base) / (1.0 - base);
    // Clamp into (0, 1): downstream label propagation evaluates LOG(1 - sim)
    // (Fisher's method), so a similarity of exactly 1.0 would yield -inf.
    if (j < 0.0) j = 0.0;
    else if (j > 1.0 - 1e-6) j = 1.0 - 1e-6;
    (void)min_jaccard;
    return j;
}

// ─── cardinality  (linear counting over empty buckets) ────────────────────
inline double KmerSketch::cardinality() const {
    if (nonempty_ == 0) return 0.0;
    if (nonempty_ >= m_) return (double)m_;
    const double empty = (double)(m_ - nonempty_);
    return -(double)m_ * std::log(empty / (double)m_);
}

// ─── containment ─────────────────────────────────────────────────────────
inline double KmerSketch::containment(const KmerSketch& other) const {
    const double card_a = cardinality();
    if (card_a <= 0.0) return 0.0;
    const double j = jaccard(other);
    if (j <= 0.0) return 0.0;
    const double card_b = other.cardinality();
    return j * (card_a + card_b) / (card_a * (1.0 + j));
}

// ─── distance (Mash) ─────────────────────────────────────────────────────
inline double KmerSketch::distance(const KmerSketch& other, double max_distance) const {
    (void)max_distance;
    const double j = jaccard(other);
    if (j <= 0.0) return std::numeric_limits<double>::infinity();
    if (j >= 1.0) return 0.0;
    return -std::log(2.0 * j / (1.0 + j)) / static_cast<double>(kmer_size_);
}

// ─── ani ─────────────────────────────────────────────────────────────────
inline double KmerSketch::ani(const KmerSketch& other) const {
    const double j = jaccard(other);
    if (j <= 0.0) return 0.0;
    if (j >= 1.0) return 1.0;
    return std::pow(2.0 * j / (1.0 + j), 1.0 / static_cast<double>(kmer_size_));
}

// ─── merge (element-wise register min = set union) ───────────────────────
inline KmerSketch KmerSketch::merge(const KmerSketch& other) const {
    assert(m_ == other.m_ && kmer_size_ == other.kmer_size_);
    KmerSketch ret(m_, kmer_size_, b_, seed_);
    const uint32_t* __restrict__ x = reg_.get();
    const uint32_t* __restrict__ y = other.reg_.get();
    uint32_t* __restrict__ z = ret.reg_.get();
    uint32_t ne = 0;
    for (uint32_t i = 0; i < m_; ++i) {
        uint32_t v = (x[i] < y[i]) ? x[i] : y[i];
        z[i] = v;
        if (v != REG_EMPTY) ++ne;
    }
    ret.nonempty_ = ne;
    ret.sig_valid_ = false;
    return ret;
}

#undef RB_BASE_ENC
#undef RB_BASE_COMP
#undef RB_BASE_VALID

} // namespace rabbit_sketch

#endif // RABBIT_SKETCH_H_
