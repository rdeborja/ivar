#include "../src/variants_by_amplicon.h"

int check_next_pos(var_by_amp *v, uint64_t pos){
  if(v->get_next() == NULL)
    return -1;
  return (v->get_next()->get_pos() == pos) == 0;
}

int check_prev_pos(var_by_amp *v, uint64_t pos){
  if(v->get_prev() == NULL)
    return -1;
  return (v->get_prev()->get_pos() == pos) == 0;
}

int main()
{
  int success = 0;
  var_by_amp *v = new var_by_amp(10);
  v = v->get_or_add_node(10);
  if(v->get_pos() != 10){
    success += 1;
  }
  // Prepend
  v = v->get_or_add_node(7);
  success += check_next_pos(v, 10);
  // Insert after
  v = v->get_or_add_node(9);
  success += check_next_pos(v, 10);
  success += check_prev_pos(v, 7);
  // Insert before
  v = v->get_or_add_node(8);
  success += check_prev_pos(v, 7);
  success += check_next_pos(v, 9);
  // Append
  v = v->get_or_add_node(12);
  success += check_prev_pos(v, 10);
  success += ((v->get_next() == NULL) == 0);
  v = v->get_or_add_node(7);
  // Test getting and adding alleles
  allele a = {
    "A",
    1,
    0,
    0,
    0,
    0,
    0,
    ""
  };
  allele a2 = {
    "T",
    1,
    0,
    0,
    0,
    0,
    0,
    ""
  };
  
  primer *fwd = new primer();
  fwd->set_name("Fwd1");
  fwd->set_start(5);
  fwd->set_end(10);
  fwd->set_indice(0);
  primer *rev = new primer();
  rev->set_name("Rev1");
  rev->set_start(10);
  rev->set_end(12);
  rev->set_indice(1);
  primer *fwd2 = new primer();
  fwd2->set_name("Fwd2");
  fwd2->set_start(5);
  fwd2->set_end(10);
  fwd2->set_indice(2);
  primer *rev2 = new primer();
  rev2->set_name("Rev2");
  rev2->set_start(10);
  rev2->set_end(12);
  rev2->set_indice(3);

  // Add to existing allele
  v->add_allele(&a, fwd, rev);
  allele *b = v->get_or_add_allele("A", "", fwd, rev);
  b->depth += 1;
  std::cout << b->nuc << ": " << b->depth << std::endl;
  success += (b->depth == 2) == 0;

  // Add new allele
  b = v->get_or_add_allele("T", "", fwd, rev);
  b->depth += 1;
  std::cout << b->nuc << ": " << b->depth << std::endl;
  success += (b->depth == 1) == 0;

  // Check with changing primers
  b = v->get_or_add_allele("T", "", fwd2, rev);
  b->depth += 1;
  std::cout << b->nuc << ": " << b->depth << std::endl;
  success += (b->depth == 1) == 0;

  b = v->get_or_add_allele("T", "", fwd, rev2);
  b->depth += 1;
  std::cout << b->nuc << ": " << b->depth << std::endl;
  success += (b->depth == 1) == 0;

  b = v->get_or_add_allele("T", "", fwd, rev);
  b->depth += 1;
  std::cout << b->nuc << ": " << b->depth << std::endl;
  success += (b->depth == 2) == 0;

  std::vector<allele *> unique_alleles;
  std::vector<uint32_t> counts;
  v->get_distinct_variants_amp(0.01, unique_alleles, counts);
  success += (unique_alleles.size() == 2) ? 0 : 1;
  success += (unique_alleles.at(0)->nuc.compare(a.nuc) == 0) ? 0 : 1;
  success += (counts.at(0) == 1) ? 0 : 1;
  success += (unique_alleles.at(1)->nuc.compare(a2.nuc) == 0) ? 0 : 1;
  success += (counts.at(1) == 3) ? 0 : 1;

  v->print_graph(true);
  return (success == 0) ? 0 : -1;
}
