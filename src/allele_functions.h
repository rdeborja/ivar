#include<string>
#include<vector>
#include<stdint.h>

#ifndef allele_functions
#define allele_functions

struct allele{
  std::string nuc;
  uint32_t depth;
  uint32_t reverse;
  uint8_t mean_qual;
  uint32_t beg;
  uint32_t end;
  float tmp_mean_qual;
  std::string deleted_bases;	// For output of variants
  bool operator < (const allele& a) const{
    return (depth > a.depth) ? true : false;
  }
  bool operator == (const allele& a) const{
    return (nuc.compare(a.nuc) == 0 && deleted_bases.compare(a.deleted_bases) == 0) ? true : false;
  }
};

bool is_indel(allele a);
int check_allele_exists(std::string n, std::string deleted_bases, std::vector<allele> ad);
std::vector<allele> update_allele_depth(char ref,std::string bases, std::string qualities, uint8_t min_qual);
void print_allele_depths(std::vector<allele> ad);
int find_ref_in_allele(std::vector<allele> ad, char ref);
int find_ref_in_allele(std::vector<allele> ad, std::string deleted_bases, std::string ref);
char gt2iupac(char a, char b);
char codon2aa(char n1, char n2, char n3);

#endif
