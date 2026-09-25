// Standalone streaming feature exporter for experimental representations.
// No reference labels; output is 136 canonical 4-mer counts per retained contig.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <omp.h>

struct Record {
  std::string name, sequence;
  std::string parent;
  int half = -1;
  std::array<float, 136> counts{};
  std::array<float, 4> bases{};
};

int main(int argc, char **argv) {
  if (argc != 5 && !(argc == 6 && std::string(argv[5]) == "--halves")) {
    std::cerr << "usage: k4_features ASSEMBLY PREFIX MIN_LENGTH THREADS [--halves]\n";
    return 2;
  }
  const size_t min_length = std::stoull(argv[3]);
  const int threads = std::stoi(argv[4]);
  const bool halves = argc == 6;
  const auto begin = std::chrono::steady_clock::now();
  std::ifstream input(argv[1]);
  std::ofstream names(std::string(argv[2]) + ".nodes.tsv");
  std::ofstream values(std::string(argv[2]) + ".f32", std::ios::binary);
  std::ofstream base_values(std::string(argv[2]) + ".base.f32", std::ios::binary);
  if (!input || !names || !values || !base_values) throw std::runtime_error("cannot open input/output");
  std::array<int, 256> canonical{};
  std::array<int, 256> index{};
  index.fill(-1);
  int dim = 0;
  for (unsigned code = 0; code < 256; ++code) {
    unsigned rc = 0, v = code;
    for (int k = 0; k < 4; ++k) { rc = (rc << 2) | (3u - (v & 3u)); v >>= 2; }
    canonical[code] = (int)std::min(code, rc);
    if (canonical[code] == (int)code) index[code] = dim++;
  }
  if (dim != 136) throw std::runtime_error("canonical dimension mismatch");
  for (unsigned code = 0; code < 256; ++code) index[code] = index[canonical[code]];
  std::array<int, 256> bases{};
  bases.fill(-1);
  bases['A'] = bases['a'] = 0; bases['C'] = bases['c'] = 1;
  bases['G'] = bases['g'] = 2; bases['T'] = bases['t'] = 3;
  size_t emitted = 0;
  std::vector<Record> batch;
  size_t batch_bytes = 0;
  names << "index\tcontig\tlength_bp" << (halves ? "\tparent\thalf\n" : "\n");
  auto flush = [&]() {
#pragma omp parallel for num_threads(threads) schedule(dynamic, 8)
    for (size_t i = 0; i < batch.size(); ++i) {
      auto &rec = batch[i];
      unsigned code = 0, valid = 0;
      for (unsigned char ch : rec.sequence) {
        const int base = bases[ch];
        if (base < 0) { code = valid = 0; continue; }
        rec.bases[base] += 1.0f;
        code = ((code << 2) | (unsigned)base) & 255u;
        if (++valid >= 4) rec.counts[index[code]] += 1.0f;
      }
    }
    for (const auto &rec : batch) {
      names << emitted++ << '\t' << rec.name << '\t' << rec.sequence.size();
      if (halves) names << '\t' << rec.parent << '\t' << rec.half;
      names << '\n';
      values.write(reinterpret_cast<const char *>(rec.counts.data()),
                   rec.counts.size() * sizeof(float));
      base_values.write(reinterpret_cast<const char *>(rec.bases.data()),
                        rec.bases.size() * sizeof(float));
    }
    batch.clear(); batch_bytes = 0;
  };
  Record current;
  auto emit = [&]() {
    if (halves && !current.name.empty() && current.sequence.size() / 2 >= min_length) {
      const size_t mid = current.sequence.size() / 2;
      for (int half = 0; half < 2; ++half) {
        Record rec;
        rec.parent = current.name;
        rec.half = half;
        rec.name = current.name + "/fragment" + std::to_string(half);
        rec.sequence = current.sequence.substr(half ? mid : 0, half ? std::string::npos : mid);
        batch_bytes += rec.sequence.size();
        batch.push_back(std::move(rec));
      }
      if (batch_bytes >= 64u * 1024u * 1024u) flush();
    } else if (!halves && !current.name.empty() && current.sequence.size() >= min_length) {
      batch_bytes += current.sequence.size();
      batch.push_back(std::move(current));
      if (batch_bytes >= 64u * 1024u * 1024u) flush(); // memory bound, not a model parameter
    }
    current = Record{};
  };
  std::string line;
  while (std::getline(input, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (!line.empty() && line[0] == '>') {
      emit();
      current.name = line.substr(1, line.find_first_of(" \t", 1) - 1);
    } else {
      current.sequence += line;
    }
  }
  emit(); flush();
  if (!names || !values || !base_values || input.bad()) throw std::runtime_error("I/O failure");
  std::cout << "records=" << emitted << " dims=" << dim << " seconds="
            << std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count()
            << '\n';
}
