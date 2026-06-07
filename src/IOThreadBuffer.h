/*
 * IOThreadBuffer
 *
 *  Generic Header to include basic OpenMP functions
 *  when OpenMP is not enabled
 *
 */

#ifndef IO_THREAD_BUFFER_H
#define IO_THREAD_BUFFER_H

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <memory>

#include "OpenMP.h"

class IOThreadBuffer {
public:
  typedef std::vector<std::shared_ptr<std::stringstream>> ThreadBuffers;
  static const int FLUSH_TRIGGER = 1024 * 1024;

private:
  static ThreadBuffers &_get(std::ostream &os) {
    static ThreadBuffers _;
    if (_.empty()) {
#pragma omp critical(ThreadBuffer_get)
      {
        if (_.empty()) {
          _.resize(omp_get_num_threads());
          for (int i = 0; i < omp_get_num_threads(); i++) {
            _[i].reset(new std::stringstream());
          }
        }
      }
    }
    assert(!_.empty());
    assert(omp_get_thread_num() < (int)_.size());
    assert(_[omp_get_thread_num()]);
    return _;
  }
  static void clear(std::ostream &os) {
#pragma omp critical(ThreadBuffer_clear)
    { _get(os).clear(); }
  }

public:
  static std::stringstream &getMyBuffer(std::ostream &os) {
    std::stringstream &ss = *(_get(os)[omp_get_thread_num()]);
    if (ss.tellp() > FLUSH_TRIGGER) {
      flush(os, ss);
    }
    return *(_get(os)[omp_get_thread_num()]);
  }
  static void close(std::ostream &os) {
    if (omp_in_parallel()) {
      flush(os);
    } else {
      ThreadBuffers &tb = _get(os);
      for (int i = 0; i < (int)tb.size(); i++) {
        flush(os, *tb[i]);
      }
    }
    clear(os);
  }
  static std::ostream &flush(std::ostream &os) {
    std::stringstream &ss = getMyBuffer(os);
    flush(os, ss);
    return os;
  }
  static std::ostream &flush(std::ostream &os, std::stringstream &ss) {
    std::string s = ss.str();
#pragma omp critical(ThreadBufferFlush)
    {
      os << s;
      os.flush();
    }
    ss.str("");
    return os;
  }
};

#endif /* IO_THREAD_BUFFER_H */
