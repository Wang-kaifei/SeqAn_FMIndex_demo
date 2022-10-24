#include <cstdlib>
#include <fstream>
#include <iostream>
#include <seqan/index.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <stdlib.h>


using namespace seqan;

bool fmindex_encode(std::string fasta_path, std::string fmindex_path) {
  SeqFileIn seqFileIn;
  std::cout << "Encoding..." << std::endl;

  // Load fasta file
  char abs_fasta_path[1024] = "";
  if (_fullpath(abs_fasta_path, fasta_path.data(), 1024) == NULL ||
      !open(seqFileIn, toCString(abs_fasta_path))) {
    std::cerr << "ERROR: Could not open the file.\n";
    return false;
  }

  StringSet<CharString> ids;
  StringSet<String<AminoAcid>> seqs;
  try {
    readRecords(ids, seqs, seqFileIn);
  } catch (Exception const &e) {
    std::cout << "ERROR: " << e.what() << std::endl;
    return false;
  }

  // Write IDs
  std::ofstream fmID;
  try {
    fmID.open((fmindex_path + ".id").data(), std::ios::out);
  } catch (Exception const &e) {
    std::cout << "ERROR: " << e.what() << std::endl;
    return false;
  }

  for (unsigned i = 0; i < length(ids); ++i) {
    fmID << ids[i] << std::endl;
  }
  fmID.close();

  // Write indexed sequences
  Index<StringSet<String<AminoAcid>>, FMIndex<>> fmIndex(seqs);
  save(fmIndex, (fmindex_path + ".fm").data());

  return true;
}

class FMIndexDecoder {
private:
  Index<StringSet<String<AminoAcid>>, FMIndex<>> fmIndex;
  std::vector<std::string> fmID;
  std::vector<std::string> fmID_short;

public:
  FMIndexDecoder(std::string fmindex_path) {
    std::cout << "Decoding..." << std::endl;
    open(fmIndex, (fmindex_path + ".fm").data());

    std::ifstream f_fmID;
    try {
      f_fmID.open((fmindex_path + ".id").data(), std::ios::in);
    } catch (Exception const &e) {
      std::cout << "ERROR: " << e.what() << std::endl;
      return;
    }

    std::string description;
    while (getline(f_fmID, description)) {
      fmID.emplace_back(description);

      std::string::size_type pos;
      pos = std::min(description.find('\t', 0), description.find(' ', 0));
      if (pos < description.size()) {
        description = description.substr(0, pos);
      }
      fmID_short.emplace_back(description);
    }
    f_fmID.close();
  }

  bool has_pattern(std::string pattern) {
    static Finder<Index<StringSet<String<AminoAcid>>, FMIndex<>>> fmFinder(
        fmIndex);

    if (find(fmFinder, pattern)) {
      clear(fmFinder);
      return true;
    }
    clear(fmFinder);
    return false;
  }

  bool has_pattern(std::string &pattern) {
    static Finder<Index<StringSet<String<AminoAcid>>, FMIndex<>>> fmFinder(
        fmIndex);

    if (find(fmFinder, pattern)) {
      clear(fmFinder);
      return true;
    }
    clear(fmFinder);
    return false;
  }

  int get_pattern_size(std::string pattern) {
    static Finder<Index<StringSet<String<AminoAcid>>, FMIndex<>>> fmFinder(
        fmIndex);

    int count = 0;
    while (find(fmFinder, pattern)) {
      count++;
    }
    clear(fmFinder);
    return count;
  }

  int get_pattern_size(std::string &pattern) {
    static Finder<Index<StringSet<String<AminoAcid>>, FMIndex<>>> fmFinder(
        fmIndex);

    int count = 0;
    while (find(fmFinder, pattern)) {
      count++;
    }
    clear(fmFinder);
    return count;
  }

  std::vector<std::pair<std::string, int>>
  get_all_pattern(std::string pattern) {
    static Finder<Index<StringSet<String<AminoAcid>>, FMIndex<>>> fmFinder(
        fmIndex);
    std::vector<std::pair<std::string, int>> res;

    while (find(fmFinder, pattern)) {
      res.emplace_back(std::pair<std::string, int>(
          fmID_short[position(fmFinder).i1], position(fmFinder).i2));
    }
    clear(fmFinder);
    return res;
  }

  std::vector<std::pair<std::string, int>>
  get_all_pattern(std::string &pattern) {
    static Finder<Index<StringSet<String<AminoAcid>>, FMIndex<>>> fmFinder(
        fmIndex);
    std::vector<std::pair<std::string, int>> res;

    while (find(fmFinder, pattern)) {
      res.emplace_back(std::pair<std::string, int>(
          fmID_short[position(fmFinder).i1], position(fmFinder).i2));
    }
    clear(fmFinder);
    return res;
  }
};

int main() {
  std::string fasta_str = "demo_data/std8.fasta";
  std::string save_file_name_str = "tmp/std8";

  fmindex_encode(fasta_str, save_file_name_str);
  FMIndexDecoder FMIndexDecoder(save_file_name_str);

  std::cout << "Checking MRSLL..." << std::endl;
  std::cout << FMIndexDecoder.has_pattern("MRSLL") << std::endl;
  std::cout << FMIndexDecoder.get_pattern_size("MRSLL") << std::endl;
  for (auto &res : FMIndexDecoder.get_all_pattern("MRSLL")) {
    std::cout << res.first << ", " << res.second << std::endl;
  }

  std::cout << "Checking K..." << std::endl;
  std::cout << FMIndexDecoder.has_pattern("K") << std::endl;
  std::cout << FMIndexDecoder.get_pattern_size("K") << std::endl;
  for (auto &res : FMIndexDecoder.get_all_pattern("K")) {
    std::cout << res.first << ", " << res.second << std::endl;
  }

  return 0;
}