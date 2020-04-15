#include <iostream>
#include <fstream>
#include <vector>
#include "allele_functions.h"
#include "htslib/vcf.h"
#include "htslib/kstring.h"
#include "version.h"

class vcf_writer{
  vcfFile *file;
  bcf_hdr_t *hdr;
  std::string region;
  std::string sample_name;
  int init_header();
 public:
  int add_info_format_hdr();
  vcf_writer(char _mode, std::string fname, std::string region_name, std::string sample_name);
  int write_record(uint32_t pos, std::vector<allele> aalt, std::string region, char ref_nuc);
  int write_record_below_threshold(uint32_t pos, std::string region_name, char ref_nuc);
  ~vcf_writer();
};
