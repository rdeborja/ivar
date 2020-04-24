#include <iostream>

#include "primer_bed.h"
#include "allele_functions.h"

#ifndef variants_by_amplicon
#define variants_by_amplicon

class var_by_amp{
  std::uint64_t pos;
  std::vector<primer*> fwd_primers;
  std::vector<primer*> rev_primers;
  std::vector<allele*> alleles;
  std::string delim = "\t";
  var_by_amp* next;
  var_by_amp* prev;
public:
  var_by_amp(uint64_t pos);
  ~var_by_amp();
  allele* get_or_add_allele(std::string nuc, std::string deleted_bases, primer *fwd, primer *rev);
  allele* get_allele(std::string nuc, std::string deleted_bases, primer *fwd, primer *rev);
  void add_allele(allele *a, primer *fwd, primer *rev);
  std::vector<primer*> get_fwd_primers();
  std::vector<primer*> get_rev_primers();
  std::vector<allele*> get_alleles();
  void add_next(var_by_amp *n);
  void add_prev(var_by_amp *p);
  var_by_amp* get_prev();
  var_by_amp* get_next();
  var_by_amp* get_node(uint64_t pos);
  var_by_amp* get_or_add_node(uint64_t pos);
  void print_graph();
  uint64_t get_pos();
};

#endif
