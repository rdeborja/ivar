#include<stdint.h>
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<algorithm>
#include<string>
#include<regex>
#include<libgen.h>

#include "allele_functions.h"
#include "vcf_writer.h"

#ifndef call_consensus_from_pileup
#define call_consensus_from_pileup

struct ret_t {
  std::string nuc;
  std::string q;
};

int call_consensus_from_plup(std::istream &cin, std::string out_file, uint8_t min_qual, double threshold, uint8_t min_depth, char gap, bool min_coverage_flag);
ret_t get_consensus_allele(std::vector<allele> &ad, uint8_t min_qual, double threshold, char gap);

#endif
