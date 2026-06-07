#!/usr/bin/env python3
"""Apply RabbitBin tier-1 (CLI/output) and tier-2 (rename + module split) to rabbitbin.cpp."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
cpp = (SRC / "rabbitbin.cpp").read_text()
h = (SRC / "rabbitbin.h").read_text()

# ── header / branding ──────────────────────────────────────────────────────
cpp = cpp.replace('#include "metabat2.h"', '#include "rabbitbin.h"')
cpp = cpp.replace('metabat2_fkmv.cpp', 'rabbitbin.cpp')
cpp = cpp.replace('MetaBAT2-FKMV', 'RabbitBin')
cpp = cpp.replace('MetaBAT2 with FastKMV', 'RabbitBin: sketch-based')
cpp = cpp.replace('static bool g_metabat_compat = false;\n', '')  # in .h

# ── tier-2 renames (order matters) ─────────────────────────────────────────
repls = [
    ('gen_tnf_graph_sample_fused10', 'calibrate_sim_cutoff_fused'),
    ('gen_tnf_graph_sample', 'calibrate_sim_cutoff'),
    ('gen_tnf_graph_allpairs', 'build_graph_allpairs'),
    ('gen_tnf_graph', 'build_similarity_graph'),
    ('fkmv_env_w_stnf', 'rb_env_w_comp'),
    ('fkmv_env_mutual_knn_on', 'rb_env_mutual_knn_on'),
    ('fkmv_env_neg_abd_thr', 'rb_env_neg_abd_thr'),
    ('fkmv_env_split_sil', 'rb_env_split_sil'),
    ('fkmv_env_gc_norm', 'rb_env_gc_norm'),
    ('fkmv_env_pmh_on', 'rb_env_pmh_on'),
    ('recruitWithTNFMin', 'recruitSimFactor'),
    ('g_w_stnf', 'g_w_comp'),
    ('FKMV_SIM_FLOOR', 'RB_SIM_FLOOR'),
    ('fkmv_kmer_size', 'sketch_kmer_size'),
    ('fkmv_sketch_size', 'sketch_size'),
    ('fkmv_bits', 'sketch_bits'),
    ('label_propagation', 'propagate_labels'),
    ('ClassMap', 'BinMap'),
    ('tl_sTNF', 'tl_sComp'),
    ('lsTNF', 'lsComp'),
    ('g.sSCR', 'g.edgeScore'),
    ('sSCR', 'edgeScore'),
    ('g.sTNF', 'g.sComp'),
    ('sTNF', 'sComp'),
    ('pTNF', 'simCutoff'),
]
for a, b in repls:
    cpp = cpp.replace(a, b)
    h = h.replace(a, b)

# env helper block
old_env = '''static bool fkmv_env_pmh_on() {
  const char *e = getenv("FKMV_PMH");
  return !e || e[0] != '0';
}
static bool fkmv_env_mutual_knn_on() {
  const char *e = getenv("FKMV_MUTUAL_KNN");
  return !e || e[0] != '0';
}
static int fkmv_env_gc_norm() {
  const char *e = getenv("FKMV_GC_NORM");
  return e ? std::atoi(e) : 1;
}
static double fkmv_env_neg_abd_thr() {
  const char *e = getenv("FKMV_NEG_ABD");
  return e ? std::atof(e) : -0.3;
}
static double fkmv_env_w_stnf() {
  const char *e = getenv("FKMV_WSTNF");
  double v = e ? std::atof(e) : 0.5;
  return (v < 0.0) ? 0.0 : (v > 1.0 ? 1.0 : v);
}
static double fkmv_env_split_sil() {
  const char *e = getenv("FKMV_SPLIT_SIL");
  return e ? std::atof(e) : 0.70;
}'''

new_env = '''static const char *rabbit_env(const char *rabbit_key, const char *legacy_key) {
  const char *e = getenv(rabbit_key);
  if (e && e[0]) return e;
  return getenv(legacy_key);
}
static bool rb_env_pmh_on() {
  const char *e = rabbit_env("RABBIT_PMH", "FKMV_PMH");
  return !e || e[0] != '0';
}
static bool rb_env_mutual_knn_on() {
  const char *e = rabbit_env("RABBIT_MUTUAL_KNN", "FKMV_MUTUAL_KNN");
  return !e || e[0] != '0';
}
static int rb_env_gc_norm() {
  const char *e = rabbit_env("RABBIT_GC_NORM", "FKMV_GC_NORM");
  return e ? std::atoi(e) : 1;
}
static double rb_env_neg_abd_thr() {
  const char *e = rabbit_env("RABBIT_NEG_ABD", "FKMV_NEG_ABD");
  return e ? std::atof(e) : -0.3;
}
static double rb_env_w_comp() {
  const char *e = rabbit_env("RABBIT_W_COMP", "FKMV_WSTNF");
  double v = e ? std::atof(e) : 0.5;
  return (v < 0.0) ? 0.0 : (v > 1.0 ? 1.0 : v);
}
static double rb_env_split_sil() {
  const char *e = rabbit_env("RABBIT_SPLIT_SIL", "FKMV_SPLIT_SIL");
  return e ? std::atof(e) : 0.70;
}'''
# after first rename pass fkmv_env_* may already be rb_env_*
cpp = cpp.replace(old_env, new_env)
if 'rabbit_env' not in cpp:
    cpp = cpp.replace('static bool rb_env_pmh_on()', new_env.split('static bool rb_env_pmh_on()')[0] + 'static bool rb_env_pmh_on()', 1)

# getenv FKMV -> rabbit_env
for rk, lk in [
    ('RABBIT_LIBDEFLATE_MAXGB', 'FKMV_LIBDEFLATE_MAXGB'),
    ('RABBIT_NOSTORE_SEQS', 'FKMV_NOSTORE_SEQS'),
    ('RABBIT_PMHK', 'FKMV_PMHK'),
    ('RABBIT_GC_NORM_CAP', 'FKMV_GC_NORM_CAP'),
    ('RABBIT_IDF_NORM', 'FKMV_IDF_NORM'),
    ('RABBIT_EXACT_COS', 'FKMV_EXACT_COS'),
    ('RABBIT_GRAPH_INDEX', 'FKMV_GRAPH_INDEX'),
    ('RABBIT_PMH_BASE', 'FKMV_PMH_BASE'),
    ('RABBIT_EDGE_POWER', 'FKMV_EDGE_POWER'),
    ('RABBIT_MINCOMMON', 'FKMV_MINCOMMON'),
]:
    cpp = cpp.replace(f'getenv("{lk}")', f'rabbit_env("{rk}", "{lk}")')

# remove TNF stubs if present
for stub in [
    '// ── Stubs for generateTNF',
    'void generateTNF(MatrixRowType /*tnf*/, const string & /*seq*/) {}\n',
    'void generateTNF(MatrixRowType /*tnf*/, const std::vector<string> & /*seqs*/) {}\n',
    '// ── Stub for cal_tnf_dist',
    'static double cal_tnf_dist(size_t /*r1*/, size_t /*r2*/,\n                            const Matrix & /*TNF1*/, const Matrix & /*TNF2*/)\n{ return 0.0; }\n',
]:
    cpp = cpp.replace(stub, '// ── removed TNF stub\n' if 'Stubs' in stub or 'Stub for' in stub else '')

cpp = cpp.replace('#include <unistd.h>', '#include <unistd.h>\n#include <cstdio>')

# log messages
cpp = cpp.replace('Start FastKMV sketch building', 'Building composition sketches')
cpp = cpp.replace('+ Spearman ranking', '+ abundance ranking')
cpp = cpp.replace('Finished FastKMV sketch building', 'Composition sketches ready')
cpp = cpp.replace('fkmvK must be', '--sketch-k must be')
cpp = cpp.replace('fkmvM must be', '--sketch-m must be')
cpp = cpp.replace('simCutoff should be', '--sim-cutoff should be')
cpp = cpp.replace('fkmvK %d, fkmvM %d', 'sketch-k %d, sketch-m %d')

# ── tier-1 CLI: replace main options block ───────────────────────────────
start = cpp.index('int main(int ac, char *av[]) {')
end = cpp.index('  if (vm.count("help") || inFile.empty() || outFile.empty()) {', start)
new_main_head = '''int main(int ac, char *av[]) {
  po::options_description desc("RabbitBin options", 100, 50);
  desc.add_options()
      ("help,h", "Show help")
      ("assembly,a", po::value<std::string>(&inFile), "Contig FASTA assembly (gzip ok) [required]")
      ("output,o", po::value<std::string>(&outFile), "Output path prefix [required]")
      ("depth,d", po::value<std::string>(&abdFile), "Coverage depth TSV")
      ("min-contig,m", po::value<size_t>(&minContig)->default_value(2500), "Minimum contig length (>=1500)")
      ("min-small-contig", po::value<size_t>(&minSmallContig)->default_value(1000), "Min length for small-contig recruiting")
      ("max-posterior", po::value<Similarity>(&maxP)->default_value(95), "Well-connected contig percent for calibration")
      ("min-edge-score", po::value<Similarity>(&minS)->default_value(60), "Minimum edge weight (1-99)")
      ("max-edges", po::value<size_t>(&maxEdges)->default_value(200), "Max neighbors per contig")
      ("sim-cutoff", po::value<Similarity>(&simCutoff)->default_value(0), "Composition similarity cutoff x100 (0=auto)")
      ("sketch-k", po::value<int>(&sketch_kmer_size)->default_value(8), "Sketch k-mer size")
      ("sketch-m", po::value<uint32_t>(&sketch_size)->default_value(500), "Sketch size (PMH registers)")
      ("sketch-b", po::value<uint32_t>(&sketch_bits)->default_value(2), "MinHash bucket bits")
      ("no-recruit", po::value<bool>(&noAdd)->zero_tokens(), "Disable small-contig recruiting")
      ("min-recruit-cluster", po::value<size_t>(&minCS)->default_value(10), "Min cluster size for recruiting")
      ("recruit-abd-centroid", po::value<bool>(&recruitToAbdCentroid)->default_value(false)->zero_tokens(), "Recruit using abundance centroid")
      ("recruit-cutoff", po::value<Distance>(&recruitSimFactor)->default_value(0.0), "Recruit sim factor x sim-cutoff (0=off)")
      ("depth-no-variance", po::value<bool>(&cvExt)->zero_tokens(), "Depth file has no variance columns")
      ("full-header", po::value<bool>(&fullHeader)->zero_tokens(), "Keep full FASTA headers")
      ("min-coverage,x", po::value<Distance>(&minCV)->default_value(1), "Min per-sample mean coverage")
      ("min-coverage-sum", po::value<Distance>(&minCVSum)->default_value(1), "Min total mean coverage")
      ("min-bin-size,s", po::value<size_t>(&minClsSize)->default_value(200000), "Min output bin size (bp)")
      ("threads,t", po::value<size_t>(&numThreads)->default_value(0), "Threads (0=all online CPUs)")
      ("labels-only,l", po::value<bool>(&onlyLabel)->zero_tokens(), "Output contig names only")
      ("save-matrix", po::value<bool>(&saveCls)->zero_tokens(), "Save membership matrix")
      ("unbinned", po::value<bool>(&outUnbinned)->zero_tokens(), "Write unbinned FASTA")
      ("no-bin-fasta", po::value<bool>(&noBinOut)->zero_tokens(), "Skip per-bin FASTA output")
      ("no-sample-depths", po::value<bool>(&noSampleDepths)->zero_tokens(), "Omit per-sample depths in headers")
      ("seed", po::value<unsigned long long>(&seed)->default_value(0), "Random seed (0=time)")
      ("marker-seed", po::value<std::string>(&seedMarkerFile)->default_value(""), "Optional marker seed file")
      ("split-max-k", po::value<int>(&splitMaxK)->default_value(6), "Max sub-clusters per split bin")
      ("split-bins", po::value<bool>(&g_abd_split)->zero_tokens(), "Abundance bin splitting (default ON)")
      ("no-split", po::value<bool>(&g_no_abd_split)->zero_tokens(), "Disable abundance splitting")
      ("split-silhouette", po::value<double>(&g_split_sil)->default_value(rb_env_split_sil()), "Silhouette split threshold")
      ("metaBAT-compat", po::value<bool>(&g_metabat_compat)->zero_tokens(), "Legacy MetaBAT names/options")
      ("debug", po::value<bool>(&debug)->zero_tokens(), "Debug output")
      ("quiet,q", po::value<bool>(&quiet)->zero_tokens(), "Less verbose")
      ("verbose,v", po::value<bool>(&verbose)->zero_tokens(), "Verbose progress (default ON)");

  po::options_description hidden("Hidden compatibility options");
  hidden.add_options()
      ("inFile,i", po::value<std::string>(&inFile))
      ("outFile", po::value<std::string>(&outFile))
      ("abdFile", po::value<std::string>(&abdFile))
      ("minContig", po::value<size_t>(&minContig))
      ("minSmallContig", po::value<size_t>(&minSmallContig))
      ("maxP", po::value<Similarity>(&maxP))
      ("minS", po::value<Similarity>(&minS))
      ("maxEdges", po::value<size_t>(&maxEdges))
      ("simCutoff", po::value<Similarity>(&simCutoff))
      ("fkmvK", po::value<int>(&sketch_kmer_size))
      ("fkmvM", po::value<uint32_t>(&sketch_size))
      ("fkmvB", po::value<uint32_t>(&sketch_bits))
      ("noAdd", po::value<bool>(&noAdd)->zero_tokens())
      ("minRecruitingSize", po::value<size_t>(&minCS))
      ("recruitToAbdCentroid", po::value<bool>(&recruitToAbdCentroid)->zero_tokens())
      ("recruitWithTNF", po::value<Distance>(&recruitSimFactor))
      ("cvExt", po::value<bool>(&cvExt)->zero_tokens())
      ("fullHeader", po::value<bool>(&fullHeader)->zero_tokens())
      ("minCV", po::value<Distance>(&minCV))
      ("minCVSum", po::value<Distance>(&minCVSum))
      ("minClsSize", po::value<size_t>(&minClsSize))
      ("numThreads", po::value<size_t>(&numThreads))
      ("onlyLabel", po::value<bool>(&onlyLabel)->zero_tokens())
      ("saveCls", po::value<bool>(&saveCls)->zero_tokens())
      ("noBinOut", po::value<bool>(&noBinOut)->zero_tokens())
      ("noSampleDepths", po::value<bool>(&noSampleDepths)->zero_tokens())
      ("seedMarkers", po::value<std::string>(&seedMarkerFile))
      ("splitMaxK", po::value<int>(&splitMaxK))
      ("abdSplit", po::value<bool>(&g_abd_split)->zero_tokens())
      ("noAbdSplit", po::value<bool>(&g_no_abd_split)->zero_tokens())
      ("splitSil", po::value<double>(&g_split_sil));

  po::options_description all_opts;
  all_opts.add(desc).add(hidden);

  po::variables_map vm;
  po::store(po::command_line_parser(ac, av).options(all_opts).positional({}).run(), vm);
  po::notify(vm);

'''
cpp = cpp[:start] + new_main_head + cpp[end:]

# help text
cpp = cpp.replace(
    'cerr << "\\nMetaBAT2-FKMV: MetaBAT2 with FastKMV k-mer Jaccard distance "',
    'cerr << "\\nRabbitBin: sketch-based metagenome binning "')
cpp = cpp.replace('Allowed options', 'RabbitBin options')
cpp = cpp.replace('[Error!] No --inFile specified', '[Error!] --assembly is required')
cpp = cpp.replace('[Error!] No --outFile specified', '[Error!] --output is required')
cpp = cpp.replace('po::store(po::command_line_parser(ac, av).options(desc)', 'po::store(po::command_line_parser(ac, av).options(all_opts)')

# output path helpers + output_bins tier1 changes
out_helpers = '''
static std::string rb_bins_info_path() {
  return g_metabat_compat ? outFile + ".BinInfo.txt" : outFile + ".bins.tsv";
}
static std::string rb_members_path() {
  return g_metabat_compat ? outFile + ".BinMembers.txt" : outFile + ".members.tsv";
}
static std::string rb_matrix_path() {
  return g_metabat_compat ? outFile + ".MemberMatrix.txt" : outFile + ".members.matrix.tsv";
}
static std::string rb_bin_output_path(size_t bin_id) {
  if (g_metabat_compat) {
    std::string p = outFile + "." + boost::lexical_cast<std::string>(bin_id);
    if (!onlyLabel) p += ".fa";
    return p;
  }
  char buf[512];
  std::snprintf(buf, sizeof(buf), "%s_bin_%03zu%s", outFile.c_str(), bin_id, onlyLabel ? "" : ".fa");
  return buf;
}
static std::string rb_unbinned_path() {
  std::string p = outFile + ".unbinned";
  if (!onlyLabel) p += ".fa";
  return p;
}

'''
cpp = cpp.replace('void output_bins(BinMap &cls) {', out_helpers + 'void output_bins(BinMap &cls) {')
cpp = cpp.replace('std::string outFile_info = outFile + ".BinInfo.txt";', 'std::string outFile_info = rb_bins_info_path();')
cpp = cpp.replace('std::string outFile_members = outFile + ".BinMembers.txt";', 'std::string outFile_members = rb_members_path();')
cpp = cpp.replace('string outFile_matrix = outFile + ".MemberMatrix.txt";', 'string outFile_matrix = rb_matrix_path();')
cpp = cpp.replace(
    '        std::string outFile_cls = outFile + ".";\n        outFile_cls.append(boost::lexical_cast<std::string>(bin_id));\n        if (!onlyLabel) outFile_cls.append(".fa");',
    '        std::string outFile_cls = rb_bin_output_path(bin_id);')
cpp = cpp.replace(
    '          std::string outFile_cls = outFile + ".";\n          outFile_cls.append("unbinned");\n          if (!onlyLabel) outFile_cls.append(".fa");',
    '          std::string outFile_cls = rb_unbinned_path();')

# header: propagate_labels declaration only
prop_start = h.index('int propagate_labels(')
prop_end = h.index('struct CompareEdge', prop_start)
h = h[:prop_start] + (
    'int propagate_labels(Graph &g, std::vector<size_t> &membership,\n'
    '                      std::vector<size_t> &node_order);\n\n'
) + h[prop_end:]
h = h.replace('void rescue_singletons(BinMap &cls);', 'void rescue_singletons(BinMap &cls);')
h = h.replace('size_t calibrate_sim_cutoff(Distance coverage = 1., bool full = false);',
              'size_t calibrate_sim_cutoff(Distance coverage = 1., bool full = false);')

(SRC / "rabbitbin.cpp").write_text(cpp)
(SRC / "rabbitbin.h").write_text(h)

# ── module split ───────────────────────────────────────────────────────────
lines = cpp.splitlines(True)
main_close = next(i for i, l in enumerate(lines) if l.strip() == 'return 0;' and i > 1500)
# find closing brace of main after return 0
main_end = main_close + 1
while main_end < len(lines) and lines[main_end].strip() != '}':
    main_end += 1
main_end += 1

impl = SRC / "impl"
impl.mkdir(exist_ok=True)
sections = {
    'rb_split.cpp': ('// Marker-guided bin splitting', '// ═══════════════════════════════════════════════════════════════════════════\n// main'),
    'rb_graph.cpp': ('// build_similarity_graph', '// Unchanged helper functions'),
    'rb_abundance.cpp': ('// Unchanged helper functions', 'void rescue_singletons'),
    'rb_output.cpp': ('static std::string rb_bins_info_path', None),
}

def extract_marker(start_marker, end_marker):
    si = next(i for i,l in enumerate(lines) if start_marker in l)
    if end_marker:
        ei = next(i for i,l in enumerate(lines) if end_marker in l and i > si)
        return ''.join(lines[si:ei])
    return ''.join(lines[si:])

for name, (sm, em) in sections.items():
    chunk = extract_marker(sm, em) if em else ''.join(lines[si:=next(i for i,l in enumerate(lines) if sm in l):])
    (impl / name).write_text(
        f'// RabbitBin module: {name}\n#include "rabbitbin.h"\n#include <iomanip>\n\n' + chunk)

# propagate_labels body -> rb_cluster.cpp
prop_body = h[h.index('int propagate_labels('):]  # wrong - already trimmed
cluster_src = (impl / 'rb_cluster.cpp') if (impl / 'rb_cluster.cpp').exists() else None
# extract from original lines before split removal
prop_start = next(i for i,l in enumerate(lines) if l.startswith('int propagate_labels('))
prop_end = next(i for i,l in enumerate(lines) if l.strip() == 'struct CompareEdge' or 'CompareEdge' in l and 'operator' in lines[i+1])
# propagate_labels is in .h - extract from backup h before trim - skip, compile rb_cluster from h backup

head = ''.join(lines[:main_end])
includes = '\n// ── implementation modules ─────────────────────────────────────────────\n'
for name in ['rb_cluster.cpp', 'rb_split.cpp', 'rb_graph.cpp', 'rb_abundance.cpp', 'rb_output.cpp']:
    includes += f'#include "impl/{name}"\n'

# remove split+graph+abundance+output from head - only keep through main
(SRC / "rabbitbin.cpp").write_text(head + includes)
print('Applied tier1+2; main ends at line', main_end)
