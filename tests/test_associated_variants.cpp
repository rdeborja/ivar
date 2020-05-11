#include <iostream>
#include "../src/contamination.h"

int main()
{
  int success = 0;
  std::string bam_path = "../data/test_asc_var.sorted.bam";
  samFile *in = hts_open(bam_path.c_str(), "r");
  hts_idx_t *idx = sam_index_load(in, bam_path.c_str());
  bam_hdr_t *header = sam_hdr_read(in);
  hts_itr_t *iter = NULL;
  iter  = sam_itr_querys(idx, header, header->target_name[0]);
  bam1_t *aln = bam_init1();
  var_by_amp *cur= NULL;
  primer *p1 = new primer(), *p2 = new primer();
  p1->set_name("primer1");
  p2->set_name("primer2");
  p1->set_indice(0);
  p2->set_indice(1);
  uint8_t min_qual = 0;
  while(sam_itr_next(in, iter, aln) >= 0) {
    cur = get_alleles_per_read(aln, cur, p1, p2, min_qual);
  }
  cur = cur->get_node(1181);
  success += (cur->get_associated_variants().at(0)[1182]["TTGT"] == 1) ? 0 : 1;
  success += (cur->get_associated_variants().at(0)[1183]["G-NN"] == 1) ? 0 : 1;
  success += (cur->get_associated_variants().at(0)[1289]["A"] == 2) ? 0 : 1;
  cur = cur->get_node(1716);
  cur->print_graph(false);
  cur->print_linked_variants();
  delete p1;
  delete p2;
  return success;
}

