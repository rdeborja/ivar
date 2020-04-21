#include <iostream>
#include <fstream>
#include <sstream>
#include "allele_functions.h"
#include "ref_seq.h"
#include <cmath>
#include <htslib/kfunc.h>

class variants_tsv_writer{
const std::string hdr_aa = "REGION"
  "\tPOS"
  "\tREF"
  "\tALT"
  "\tREF_DP"
  "\tREF_RV"
  "\tREF_QUAL"
  "\tALT_DP"
  "\tALT_RV"
  "\tALT_QUAL"
  "\tALT_FREQ"
  "\tTOTAL_DP"
  "\tPVAL"
  "\tPASS"
  "\tGFF_FEATURE"
  "\tREF_CODON"
  "\tREF_AA"
  "\tALT_CODON"
  "\tALT_AA";
const std::string hdr = "REGION"
  "\tPOS"
  "\tREF"
  "\tALT"
  "\tREF_DP"
  "\tREF_RV"
  "\tREF_QUAL"
  "\tALT_DP"
  "\tALT_RV"
  "\tALT_QUAL"
  "\tALT_FREQ"
  "\tTOTAL_DP"
  "\tPVAL"
  "\tPASS"
  "\tGFF_FEATURE"
  "\tREF_CODON"
  "\tREF_AA"
  "\tALT_CODON"
  "\tALT_AA";
  const std::string delim = "\t";
  std::ofstream fout;
  ref_antd *refa;
public:
  variants_tsv_writer(std::string path);
  ~variants_tsv_writer();
  int write_line(allele a, allele ref, std::string region, uint64_t pos, double *freq_depth);
  int is_valid();
  void add_ref(ref_antd *ref_a);
};
