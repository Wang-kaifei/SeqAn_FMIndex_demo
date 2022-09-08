#include <seqan/sequence.h>
#include <seqan/index.h>

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


int main()
{
	std::vector<std::string> proteins({ "VXLAGZ", "GKTVXL", "XLZ" });

	std::string save_file_name_str;
	save_file_name_str = "./tmp/fasta.esa";
	build1(proteins, save_file_name_str.data());
	load1(save_file_name_str.data());
	save_file_name_str = "./tmp/fasta.fm";
	build2(proteins, save_file_name_str.data());
	load2(save_file_name_str.data());

	return 0;
}