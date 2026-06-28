// RabbitBin module: rb_output.cpp

void output_bins(BinMap &cls) {
#pragma omp parallel
  {
#pragma omp single
    {
      const bool do_qc  = g_qc_annotate && g_marker_set_size > 0 && !g_allcontig_markers.empty();
      const bool do_tax = !g_contig_taxon.empty();

      std::string outFile_info = rb_bins_info_path();
      std::ofstream os_info(outFile_info.c_str());
      verbose_message("Writing cluster stats to: %s\n", outFile_info.c_str());
      os_info << "BinNum\tNumContigs\tTotalLength\tLengthWeightedAvgCoveage\tFileName";
      if (do_qc)  os_info << "\tCompleteness\tContamination\tTier";
      if (do_tax) os_info << "\tTaxonomy";
      os_info << "\n";

      std::string outFile_members = rb_members_path();
      std::ofstream os_members(outFile_members.c_str());
      verbose_message("Writing Members stats to: %s\n", outFile_members.c_str());
      os_members << "BinNum\tSequenceName\tTotalDepth";
      if (g_emit_confidence) os_members << "\tConfidence";
      os_members << "\n";

      // Soft-assignment map: raw cluster label -> sequential bin id, so the
      // per-contig runner-up (second-choice) bin can be reported by its output
      // BinNum.  Filled as bins are emitted below (feature #8).
      std::unordered_map<size_t, size_t> label2binid;

      // ── Strain / SNV output (feature #5) ──────────────────────────────────
      // Per-contig SNV vectors -> <prefix>.snv.tsv; per-bin strain heterogeneity
      // -> <prefix>.strain.tsv.  Written inline so BinNum matches members.tsv
      // exactly (same size/QC filters).  A bin whose length-weighted mean
      // nucleotide diversity (pi) exceeds STRAIN_PI_THR is flagged multi-strain:
      // composition+depth placed those contigs together but the reads carry
      // intra-bin polymorphism that only strain-aware resolution can split.
      const bool do_snv = g_strain_scan && !g_snv_pi.empty();
      static const double STRAIN_PI_THR = 0.005;
      std::ofstream os_snv, os_strain;
      if (do_snv) {
        std::string snv_path = std::string(outFile) + ".snv.tsv";
        std::string str_path = std::string(outFile) + ".strain.tsv";
        os_snv.open(snv_path.c_str());
        os_strain.open(str_path.c_str());
        verbose_message("Writing per-contig SNV vectors to: %s\n", snv_path.c_str());
        verbose_message("Writing per-bin strain heterogeneity to: %s\n", str_path.c_str());
        os_snv << "BinNum\tContigName\tLength\tMeanPi\tMaxSNV\tSNVperKb";
        for (int i = 0; i < (int)num_depth_samples; ++i)
          os_snv << "\tPi_s" << (i + 1);
        os_snv << "\n";
        os_strain << "BinNum\tNumContigs\tTotalBp\tMeanPi\tMeanSNVperKb"
                     "\tMultiStrain\tEstStrains";
        for (int i = 0; i < (int)num_depth_samples; ++i)
          os_strain << "\tMinorAbund_s" << (i + 1);
        os_strain << "\n";
      }

      // Optional CAMI bioboxes output (<prefix>.binning).
      std::ofstream os_bbx;
      if (g_write_bioboxes) {
        std::string bbx = outFile + ".binning";
        os_bbx.open(bbx.c_str());
        verbose_message("Writing bioboxes binning to: %s\n", bbx.c_str());
        os_bbx << "@Version:0.9.0\n@SampleID:" << outFile << "\n"
               << "@@SEQUENCEID\tBINID\t_LENGTH\n";
      }

      Distance binnedSize = 0, binnedSize1 = 0;
      std::vector<size_t> clsMap(nobs + nobs1, 0);

      size_t bin_id = 1;
      for (auto it = cls.begin(); it != cls.end(); ++it) {
        size_t kk = it->first;
        size_t numContigs = 0, totalLength = 0;
        double lengthWeightedAvgCoverage = 0;
        size_t s = 0, s1 = 0;
        std::string outFile_cls = rb_bin_output_path(bin_id);
        // Per-bin strain accumulators (large contigs only carry SNV vectors).
        std::stringstream snv_ss;
        double bin_pi_lw = 0.0, bin_snvkb_lw = 0.0;
        size_t bin_snv_bp = 0, bin_snv_n = 0;
        {
          const auto &cluster = it->second;
          std::stringstream ss;
          for (auto it2 = cluster.begin(); it2 != cluster.end(); ++it2) {
            size_t len = 0;
            auto idx = (*it2 < nobs) ? *it2 : *it2 - nobs;
            if (*it2 < nobs) { len = seq_lens[idx]; s += len; }
            else             { len = small_seq_lens[idx]; s1 += len; }

            string &sequence_name = (*it2 < nobs)
                ? contig_names[idx] : small_contig_names[idx];
            numContigs++;
            totalLength += len;
            double total = 0.0;
            for (auto i = 0; i < (int)num_depth_samples; i++)
              total += (*it2 < nobs) ? depth_matrix(idx, i) : small_depth_matrix(idx, i);
            lengthWeightedAvgCoverage += len * total;
            ss << bin_id << "\t" << sequence_name << "\t" << total;
            if (g_emit_confidence) {
              // Large contigs are graph nodes with a propagation margin; small
              // recruited contigs have no node → NA.
              if (*it2 < nobs && idx < g_node_confidence.size() &&
                  g_node_confidence[idx] >= 0.0f)
                ss << "\t" << std::fixed << std::setprecision(4)
                   << g_node_confidence[idx];
              else
                ss << "\tNA";
              ss.unsetf(std::ios::fixed);
            }
            ss << "\n";

            // Per-contig SNV vector (large contigs only; small recruits have no
            // SNV row).  MeanPi averages over samples that covered the contig.
            if (do_snv && *it2 < nobs && idx < nobs) {
              const size_t S = num_depth_samples;
              double pisum = 0.0;
              uint32_t maxsnv = 0;
              int ncov = 0;
              for (size_t b = 0; b < S; ++b) {
                const size_t o = idx * S + b;
                if (o < g_snv_cov.size() && g_snv_cov[o] > 0) {
                  pisum += g_snv_pi[o];
                  ncov++;
                }
                if (o < g_snv_nsnv.size() && g_snv_nsnv[o] > maxsnv)
                  maxsnv = g_snv_nsnv[o];
              }
              double meanpi = ncov ? pisum / ncov : 0.0;
              double snvkb = len ? (double)maxsnv / (double)len * 1000.0 : 0.0;
              snv_ss << bin_id << "\t" << sequence_name << "\t" << len << "\t"
                     << std::fixed << std::setprecision(5) << meanpi << "\t"
                     << maxsnv << "\t" << std::setprecision(3) << snvkb;
              snv_ss << std::setprecision(5);
              for (size_t b = 0; b < S; ++b) {
                const size_t o = idx * S + b;
                snv_ss << "\t" << ((o < g_snv_pi.size()) ? g_snv_pi[o] : 0.0f);
              }
              snv_ss << "\n";
              snv_ss.unsetf(std::ios::fixed);
              bin_pi_lw += meanpi * (double)len;
              bin_snvkb_lw += snvkb * (double)len;
              bin_snv_bp += len;
              bin_snv_n++;
            }
          }
          if (s + s1 < min_bin_bp) continue;

          // Per-bin SCG quality / taxonomy (features #1/#8) + --keep-hq-only.
          double comp = 0.0, cont = 0.0; int taxid = -1;
          if (do_qc)  rb_bin_qc(cluster, comp, cont);
          if (do_tax) taxid = rb_bin_taxon(cluster);
          bool is_hq = do_qc && comp > g_hq_comp && cont < g_hq_cont;
          if (g_keep_hq_only && !is_hq) continue;

          os_members << ss.str();
          // Strain / SNV (feature #5): the bin passed the size/QC filters, so
          // emit its per-contig SNV rows + one strain-summary row under this
          // (now final) BinNum.
          if (do_snv) {
            os_snv << snv_ss.str();
            double meanpi = bin_snv_bp ? bin_pi_lw / (double)bin_snv_bp : 0.0;
            double meankb = bin_snv_bp ? bin_snvkb_lw / (double)bin_snv_bp : 0.0;
            // Linkage-based strain count + dominant minor-strain abundance.
            std::vector<double> minorAbund;
            int estStrains = estimate_bin_strains(it->second, minorAbund);
            bool multi = (estStrains > 1) || (meanpi >= STRAIN_PI_THR);
            os_strain << bin_id << "\t" << bin_snv_n << "\t" << bin_snv_bp << "\t"
                      << std::fixed << std::setprecision(5) << meanpi << "\t"
                      << std::setprecision(3) << meankb << "\t"
                      << (multi ? "Y" : "N") << "\t"
                      << estStrains;
            os_strain << std::setprecision(4);
            for (int i = 0; i < (int)num_depth_samples; ++i)
              os_strain << "\t"
                        << ((estStrains > 1 && i < (int)minorAbund.size())
                                ? minorAbund[i] : 0.0);
            os_strain << "\n";
            os_strain.unsetf(std::ios::fixed);
          }
          for (size_t i = 0; i < cluster.size(); ++i) {
            assert(cluster[i] < (int)clsMap.size());
            clsMap[cluster[i]] = kk + 1;
          }
          // Map raw cluster label -> sequential output BinNum for the soft
          // assignment file (feature #8).
          if (g_emit_confidence) label2binid[kk] = bin_id;
          if (g_write_bioboxes) {
            for (auto c : cluster) {
              auto idx = (c < nobs) ? c : c - nobs;
              const string &n = (c < nobs) ? contig_names[idx] : small_contig_names[idx];
              size_t len = (c < nobs) ? seq_lens[idx] : small_seq_lens[idx];
              os_bbx << n << "\t" << bin_id << "\t" << len << "\n";
            }
          }

          if (totalLength > 0) lengthWeightedAvgCoverage /= totalLength;
          else                 lengthWeightedAvgCoverage = 0;
          os_info << bin_id << "\t" << numContigs << "\t" << totalLength << "\t"
                  << lengthWeightedAvgCoverage << "\t" << outFile_cls;
          if (do_qc) {
            const char *tier = is_hq ? "HQ"
                             : (comp >= g_mq_comp && cont < g_mq_cont) ? "MQ" : "LQ";
            os_info << "\t" << std::fixed << std::setprecision(2) << comp
                    << "\t" << cont << "\t" << tier;
            os_info.unsetf(std::ios::fixed);
          }
          if (do_tax)
            os_info << "\t" << (taxid >= 0 ? g_taxon_names[taxid] : "unclassified");
          os_info << "\n";
        }

        binnedSize  += s;
        binnedSize1 += s1;

#pragma omp task firstprivate(outFile_cls)
        if (!noBinOut) {
          auto &cluster = it->second;
          std::sort(cluster.begin(), cluster.end());
          size_t bases = 0;
          std::ofstream os(outFile_cls.c_str());
          if (!os) {
            cerr << "[Error!] Could not write to " << outFile_cls << "\n";
            exit(1);
          }
          char os_buffer[buf_size];
          os.rdbuf()->pubsetbuf(os_buffer, buf_size);
          for (auto it2 = cluster.begin(); it2 != cluster.end(); ++it2) {
            auto idx = (*it2 < nobs) ? *it2 : *it2 - nobs;
            const string &n = (*it2 < nobs) ? contig_names[idx] : small_contig_names[idx];
            std::stringstream labelss;
            labelss << n;
            if (depth_file.length()) {
              labelss << " total_depth=";
              double total = 0.0;
              for (auto i = 0; i < (int)num_depth_samples; i++)
                total += (*it2 < nobs) ? depth_matrix(idx, i) : small_depth_matrix(idx, i);
              labelss << std::fixed << std::setprecision(2) << total;
              if (!noSampleDepths) {
                labelss << " sample_depths" << std::fixed << std::setprecision(1);
                for (auto i = 0; i < (int)num_depth_samples; i++) {
                  labelss << (i == 0 ? "=" : ",");
                  auto d = (*it2 < nobs) ? depth_matrix(idx, i) : small_depth_matrix(idx, i);
                  if (d >= 0.1) labelss << d;
                  else          labelss << "0";
                }
              }
            }
            std::string label = labelss.str();
            if (onlyLabel) {
              os << label << line_delim;
            } else {
              std::string_view seq = (*it2 < nobs) ? seqs[idx] : small_seqs[idx];
              printFasta(os, label, seq);
              bases += (*it2 < nobs) ? seq_lens[idx] : small_seq_lens[idx];
            }
          }
          os.close();
          if (!os) {
            cerr << "[Error!] Failed to write to " << outFile_cls << "\n";
            exit(1);
          }
          if (debug)
#pragma omp critical(COUT)
            cout << "Bin " << bin_id << " (" << bases << " bases in "
                 << cluster.size() << " contigs) was saved to: "
                 << outFile_cls << "\n";
        }
        bin_id++;
      }

      if (saveCls) {
#pragma omp task
        {
          string outFile_matrix = rb_matrix_path();
          if (verbose) verbose_message("Saving cluster membership matrix to %s\n",
                                       outFile_matrix.c_str());
          std::ofstream os(outFile_matrix.c_str());
          if (!os) {
            cerr << "[Error!] Could not write cluster membership to "
                 << outFile_matrix << "\n";
            exit(1);
          }
          char os_buffer[buf_size];
          os.rdbuf()->pubsetbuf(os_buffer, buf_size);
          os << "ContigName" << tab_delim << "ClusterId" << line_delim;
          for (size_t i = 0; i < nobs; ++i)
            os << contig_names[i] << tab_delim << clsMap[i] << line_delim;
          for (size_t i = nobs; i < nobs + nobs1; ++i) {
            auto idx = i - nobs;
            os << small_contig_names[idx] << tab_delim << clsMap[i] << line_delim;
          }
          os.flush(); os.close();
          if (!os) {
            cerr << "[Error!] Failed to write cluster membership\n";
            exit(1);
          }
        }
      }

      // ── Soft assignment / confidence file (feature #8) ────────────────────
      if (g_emit_confidence) {
#pragma omp task firstprivate(label2binid)
        {
          std::string path = std::string(outFile) + ".confidence.tsv";
          if (verbose) verbose_message("Saving per-contig confidence to %s\n",
                                       path.c_str());
          std::ofstream os(path.c_str());
          if (os) {
            os << "ContigName\tBinNum\tConfidence\tSecondBinNum\tSecondScore\n";
            for (size_t i = 0; i < nobs; ++i) {
              if (clsMap[i] == 0) continue;              // unbinned large contig
              // clsMap stores rawlabel+1; map to the sequential output BinNum.
              auto bit = label2binid.find(clsMap[i] - 1);
              if (bit == label2binid.end()) continue;
              size_t binid = bit->second;
              float conf = (i < g_node_confidence.size()) ? g_node_confidence[i] : -1.0f;
              long second_bin = -1;
              float second_sc = 0.0f;
              if (i < g_node_second.size() && g_node_second[i] != SIZE_MAX) {
                auto it = label2binid.find(g_node_second[i]);
                if (it != label2binid.end()) second_bin = (long)it->second;
                second_sc = g_node_second_score[i];
              }
              os << contig_names[i] << "\t" << binid << "\t";
              if (conf >= 0.0f) os << std::fixed << std::setprecision(4) << conf;
              else              os << "NA";
              os.unsetf(std::ios::fixed);
              os << "\t" << second_bin << "\t"
                 << std::fixed << std::setprecision(4) << second_sc << "\n";
              os.unsetf(std::ios::fixed);
            }
            os.flush(); os.close();
          } else {
            cerr << "[Warn] could not write confidence file " << path << "\n";
          }
        }
      }

      if (outUnbinned) {
#pragma omp task
        {
          std::string outFile_cls = rb_unbinned_path();
          if (verbose) verbose_message("Saving unbinned contigs to %s\n",
                                       outFile_cls.c_str());
          std::ofstream os(outFile_cls.c_str());
          if (!os) {
            cerr << "[Error!] Could not write unbinned contigs to "
                 << outFile_cls << "\n";
            exit(1);
          }
          char os_buffer[buf_size];
          os.rdbuf()->pubsetbuf(os_buffer, buf_size);
          for (size_t i = 0; i < clsMap.size(); ++i) {
            if (clsMap[i] == 0) {
              auto idx = (i < nobs) ? i : i - nobs;
              std::string &label = (i < nobs) ? contig_names[idx] : small_contig_names[idx];
              if (onlyLabel) {
                os << label << line_delim;
              } else {
                std::string_view seq = (i < nobs) ? seqs[idx] : small_seqs[idx];
                printFasta(os, label, seq);
              }
            }
          }
          os.flush(); os.close();
          if (!os) {
            cerr << "[Error!] Failed to write unbinned contigs\n";
            exit(1);
          }
        }
      }

#pragma omp taskwait
      if (verbose) {
        verbose_message(
            "%2.2f%% (%lld bases) of large (>=%d) and %2.2f%% (%lld bases) "
            "of small (<%d) contigs were binned.\n",
            (Distance)binnedSize / totalSize * 100.,
            (unsigned long long)binnedSize, minContig,
            binnedSize1 == 0 ? 0 : (Distance)binnedSize1 / totalSize1 * 100.,
            (unsigned long long)binnedSize1, minContig);
      }
#pragma omp critical(COUT)
      {
        cout.precision(20);
        cout << bin_id - 1 << " bins (" << binnedSize + binnedSize1
             << " bases in total) formed." << std::endl;
      }
    } // omp single
  }   // omp parallel
}
