#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>
#include <stdlib.h>
#include <cstdlib>

using namespace seqan;


void build1(const std::vector<std::string>& proteins_list, const char* save_file_name) {
	StringSet<String<AminoAcid> > protein_string_set;
	for (auto& str : proteins_list) {
		appendValue(protein_string_set, str);
	}

	Index<StringSet<String<AminoAcid> >, IndexEsa<> > fmIndex(protein_string_set);
	save(fmIndex, save_file_name);
}


void load1(const char* save_file_name) {
	Index<StringSet<String<AminoAcid> >, IndexEsa<> > fmIndex;
	open(fmIndex, save_file_name);

	Finder<Index<StringSet<String<AminoAcid>>, IndexEsa<> > > fmFinder(fmIndex);

	while (find(fmFinder, "XL"))
	{
		std::cout << position(fmFinder).i1 << ", " << position(fmFinder).i2 << std::endl;
	}
}


void build2(const std::vector<std::string>& proteins_list, const char* save_file_name) {
	StringSet<String<AminoAcid> > protein_string_set;
	for (auto& str : proteins_list) {
		appendValue(protein_string_set, str);
	}

	Index<StringSet<String<AminoAcid> >, FMIndex<> > fmIndex(protein_string_set);
	save(fmIndex, save_file_name);
}


void load2(const char* save_file_name) {
	Index<StringSet<String<AminoAcid> >, FMIndex<> > fmIndex;
	open(fmIndex, save_file_name);

	Finder<Index<StringSet<String<AminoAcid>>, FMIndex<> > > fmFinder(fmIndex);

	while (find(fmFinder, "XL"))
	{
		std::cout << position(fmFinder).i1 << ", " << position(fmFinder).i2 << std::endl;
	}
}


bool fmindex_encode(std::string fasta_path, std::string fmindex_path) {
	SeqFileIn seqFileIn;

	// Load fasta file
	char abs_fasta_path[1024] = "";
	if (_fullpath(abs_fasta_path, fasta_path.data(), 1024) == NULL || !open(seqFileIn, toCString(abs_fasta_path)))
	{
		std::cerr << "ERROR: Could not open the file.\n";
		return false;
	}

	StringSet<CharString> ids;
	StringSet<String<AminoAcid>> seqs;
	try
	{
		readRecords(ids, seqs, seqFileIn);
	}
	catch (Exception const& e)
	{
		std::cout << "ERROR: " << e.what() << std::endl;
		return false;
	}

	// Encode fasta file
	for (unsigned i = 0; i < length(ids); ++i)
		std::cout << ids[i] << '\t' << seqs[i] << '\n';

	Index<StringSet<String<AminoAcid> >, FMIndex<> > fmIndex(seqs);
	save(fmIndex, fmindex_path.data());

	return true;
}

class FMIndexDecoder {
private:
	Index<StringSet<String<AminoAcid> >, FMIndex<> > fmIndex;

public:
	FMIndexDecoder(std::string fmindex_path) {
		open(fmIndex, fmindex_path.data());
	}

	void decode(std::string& pattern) {
		static Finder<Index<StringSet<String<AminoAcid>>, FMIndex<> > > fmFinder(fmIndex);

		while (find(fmFinder, pattern))
		{
			std::cout << position(fmFinder).i1 << ", " << position(fmFinder).i2 << std::endl;
		}
		clear(fmFinder);
	}

	void decode(std::string pattern) {
		static Finder<Index<StringSet<String<AminoAcid>>, FMIndex<> > > fmFinder(fmIndex);

		while (find(fmFinder, pattern))
		{
			std::cout << position(fmFinder).i1 << ", " << position(fmFinder).i2 << std::endl;
		}
		clear(fmFinder);
	}
};


int main()
{
	std::string fasta_str;
	fasta_str = "demo_data/std8.fasta";
	/*
	CharString seqFileName;
	SeqFileIn seqFileIn;
	open(seqFileIn, "D:/demo_data/std8.fasta"); // input MUST be absolute path
	close(seqFileIn);
	open(seqFileIn, "D:/SeqAn_demo/demo_data/std8.fasta"); // input MUST be absolute path
	close(seqFileIn);
	open(seqFileIn, "D:/SeqAn_demo/demo_data/std8.FASTA"); // input MUST be absolute path
	close(seqFileIn);

	char abs_path[1024] = "";
	_fullpath(abs_path, "demo_data/std8.fasta", 1024);
	open(seqFileIn, toCString(abs_path));
	close(seqFileIn);

	seqFileName = getAbsolutePath("/demo_data/std8.fasta");
	open(seqFileIn, toCString(seqFileName));
	close(seqFileIn);

	seqFileName = getAbsolutePath("./demo_data/std8.fasta");
	open(seqFileIn, toCString(seqFileName));
	close(seqFileIn);

	char* fasta_path2 = "/demo_data/std8.fasta";
	CharString seqFileName2 = getAbsolutePath(fasta_path2);
	SeqFileIn seqFileIn2(toCString(fasta_path2));
	char* fasta_path3 = "./demo_data/std8.fasta";
	CharString seqFileName3 = getAbsolutePath(fasta_path3);
	SeqFileIn seqFileIn3(toCString(fasta_path3));
	*/

	std::vector<std::string> proteins({ "VXLAGZ", "GKTVXL", "XLZ" });

	std::string save_file_name_str;
	save_file_name_str = "tmp/fasta.esa";
	build1(proteins, save_file_name_str.data());
	load1(save_file_name_str.data());
	save_file_name_str = "tmp/fasta.fm";
	build2(proteins, save_file_name_str.data());
	load2(save_file_name_str.data());
	save_file_name_str = "demo_data/std8.fm";

	fmindex_encode(fasta_str, save_file_name_str);
	FMIndexDecoder FMIndexDecoder(save_file_name_str);
	FMIndexDecoder.decode("MRSLL");
	FMIndexDecoder.decode("MGSIG");
	FMIndexDecoder.decode("MKSTI");

	return 0;
}