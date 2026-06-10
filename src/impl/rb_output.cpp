// RabbitBin module: rb_output.cpp

void output_bins(BinMap &cls) {
#pragma omp parallel
  {
#pragma omp single
    {
      std::string outFile_info = rb_bins_info_path();
      std::ofstream os_info(outFile_info.c_str());
      verbose_message("Writing cluster stats to: %s\n", outFile_info.c_str());
      os_info << "BinNum\tNumContigs\tTotalLength\tLengthWeightedAvgCoveage\tFileName\n";

      std::string outFile_members = rb_members_path();
      std::ofstream os_members(outFile_members.c_str());
      verbose_message("Writing Members stats to: %s\n", outFile_members.c_str());
      os_members << "BinNum\tSequenceName\tTotalDepth\n";

      Distance binnedSize = 0, binnedSize1 = 0;
      std::vector<size_t> clsMap(nobs + nobs1, 0);

      size_t bin_id = 1;
      for (auto it = cls.begin(); it != cls.end(); ++it) {
        size_t kk = it->first;
        size_t numContigs = 0, totalLength = 0;
        double lengthWeightedAvgCoverage = 0;
        size_t s = 0, s1 = 0;
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
            ss << bin_id << "\t" << sequence_name << "\t" << total << "\n";
          }
          if (s + s1 < min_bin_bp) continue;
          os_members << ss.str();
          for (size_t i = 0; i < cluster.size(); ++i) {
            assert(cluster[i] < (int)clsMap.size());
            clsMap[cluster[i]] = kk + 1;
          }
        }

        std::string outFile_cls = rb_bin_output_path(bin_id);

        if (totalLength > 0) lengthWeightedAvgCoverage /= totalLength;
        else                 lengthWeightedAvgCoverage = 0;
        os_info << bin_id << "\t" << numContigs << "\t" << totalLength << "\t"
                << lengthWeightedAvgCoverage << "\t" << outFile_cls << "\n";

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
