#include <iostream>
#include <cmath>
#include <map>
#include <numeric>

#include "primer_bed.h"
#include "allele_functions.h"
#include <htslib/kfunc.h>

#ifndef variants_by_amplicon
#define variants_by_amplicon

class var_by_amp{
  std::uint64_t pos;
  std::vector<primer*> fwd_primers;
  std::vector<primer*> rev_primers;
  std::vector<allele*> alleles;	// Every allele has an associated map of positions and variants counts
  std::vector<std::map<uint32_t, std::map<std::string, uint32_t>>> associated_variants; // map[pos] = map[allele, count_at_pos] for every amplicon covering position
  std::string delim = "\t";
  var_by_amp* next;
  var_by_amp* prev;
public:
  var_by_amp(uint64_t pos);
  ~var_by_amp();
  allele* get_or_add_allele(std::string nuc, std::string deleted_bases, primer *fwd, primer *rev);
  allele* get_allele(std::string nuc, std::string deleted_bases, primer *fwd, primer *rev, int &allele_ind);
  void add_allele(allele *a, primer *fwd, primer *rev);
  std::vector<primer*> get_fwd_primers();
  std::vector<primer*> get_rev_primers();
  std::vector<allele*> get_alleles();
  void add_next(var_by_amp *n);
  void add_prev(var_by_amp *p);
  var_by_amp* get_prev();
  var_by_amp* get_next();
  var_by_amp* get_node(uint64_t pos);
  std::vector<std::map<uint32_t, std::map<std::string, uint32_t>>> get_associated_variants();
  var_by_amp* get_or_add_node(uint64_t pos);
  void print_graph(bool recurse);
  void print_graph(double min_freq);
  uint64_t get_pos();
  uint32_t get_depth();
  std::vector<int> get_alleles_above_freq(double min_freq);
  void get_distinct_variants_amp(double min_freq, std::vector<allele*> &unique_alleles, std::vector<uint32_t> &counts, uint &unique_primers);
  int get_num_unique_primers();
  void get_linked_variants_on_amplicon(int allele_indice);
  int add_associated_variants(uint32_t pos, allele *aa, allele *a, primer *fwd, primer *rev); // Add associated variant for a particular allele
  void print_linked_variants();
};

double chisqr(int dof, double cv);
double igf(double s, double z);
double* chisqr_goodness_of_fit(uint32_t ctable[16][16], int nrow, int ncol);
double compute_critical_value(uint32_t ctable[16][16], int nrow, int ncol);

#endif
