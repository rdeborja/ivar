#include "contamination.h"

void identify_amp(std::vector<primer> primers, uint64_t pos, uint64_t end_pos, std::vector<primer> &fwd_primers, std::vector<primer> &rev_primers){
  fwd_primers.clear();
  rev_primers.clear();
  primer primer_pair;
  uint16_t p2;
  for (std::vector<primer>::iterator it = primers.begin(); it != primers.end(); ++it) {
    if(it->get_pair_indice()== -1 || it->get_start() > pos)
      continue;
    p2 = it->get_pair_indice();
    if(pos >= it->get_start() && end_pos <= primers.at(p2).get_end()){
      fwd_primers.push_back(*it);
      rev_primers.push_back(primers.at(p2));
    }
  }
}

var_by_amp* get_alleles_per_read(bam1_t *aln, var_by_amp *v, primer *fwd, primer *rev){
  std::vector<allele> all;
  int ncigar = aln->core.n_cigar;
  uint32_t* cigar = bam_get_cigar(aln);
  uint32_t cig;
  int i = 0, len, ctr = 0;
  uint64_t pos = aln->core.pos, ql=0;
  uint8_t *seq = bam_get_seq(aln);
  allele *a;
  std::string indel;
  if(v == NULL)
    v = new var_by_amp(pos);
  while (i < ncigar ) {
    cig = bam_cigar_op(cigar[i]);
    len = bam_cigar_oplen(cigar[i]);
    if ((bam_cigar_type(cigar[i]) & 1) && (bam_cigar_type(cigar[i]) & 2)){ // Query and ref consuming
      ctr = 0;
      while(ctr < len -1){
	v = v->get_or_add_node(pos);
	a = v->get_or_add_allele(std::string(1, seq_nt16_str[bam_seqi(seq, ql+ctr)]), "", fwd, rev);
	a->depth += 1;
	ctr++;
	pos++;
      }
      if ( i < ncigar - 1){			// For final loop check for only query or ref consuming cigar in next position
	if((bam_cigar_type(cigar[i+1]) & 1) && !(bam_cigar_type(cigar[i+1]) & 2)){ // Only query consuming
	  indel.clear();
	  indel.assign(1, seq_nt16_str[bam_seqi(seq, ql+ctr)]);
	  ctr = 0;
	  ql += len;
	  while(ctr < bam_cigar_oplen(cigar[i+1])){
	    indel += seq_nt16_str[bam_seqi(seq, ql+ctr)];
	    ctr++;
	  }
	  v = v->get_or_add_node(pos);
	  a = v->get_or_add_allele(indel, "", fwd, rev);
	  a->depth += 1;
	  ql += bam_cigar_oplen(cigar[i+1]);
	  i++;
	} else if (!(bam_cigar_type(cigar[i+1]) & 1) && (bam_cigar_type(cigar[i+1]) & 2)){ // Only ref consuming
	  indel.clear();
	  indel.assign(1, seq_nt16_str[bam_seqi(seq, ql+ctr)]);
	  v = v->get_or_add_node(pos);
	  a = v->get_or_add_allele(indel, std::string(bam_cigar_oplen(cigar[i+1]), 'N'), fwd, rev);
	  a->depth += 1;
	  i++;
	  pos += bam_cigar_oplen(cigar[i+1]);
	} else {
	  v = v->get_or_add_node(pos);
	  a = v->get_or_add_allele(std::string(1, seq_nt16_str[bam_seqi(seq, ql+ctr)]), "", fwd, rev);
	  a->depth += 1;
	  pos++;
	}
      } else {
	v = v->get_or_add_node(pos);
	a = v->get_or_add_allele(std::string(1, seq_nt16_str[bam_seqi(seq, ql+ctr)]), "", fwd, rev);
	a->depth += 1;
	pos++;
      }
    } else if(bam_cigar_type(cigar[i]) & 1){
      ql += len;
    } else if(bam_cigar_type(cigar[i]) & 2){
      pos += len;
    }
    i++;
  }
  return v;
}

int main()
{
  std::string bed = "../data/test.bed";
  std::string bam = "../data/test.sim.sorted.bam";
  std::vector<primer> primers = populate_from_file(bed);
  populate_pair_indices(primers, "../data/pair_information.tsv");
  samFile *in = hts_open(bam.c_str(), "r");
  bam_hdr_t *header = sam_hdr_read(in);
  hts_idx_t *idx = sam_index_load(in, bam.c_str());
  hts_itr_t *iter = NULL;
  std::vector<primer> fwd, rev;
  primer fwd_max_end, rev_min_start;
  iter  = sam_itr_querys(idx, header, header->target_name[0]);
  var_by_amp *cur= NULL;
  bam1_t *aln = bam_init1();
  while(sam_itr_next(in, iter, aln) >= 0) {
    identify_amp(primers, aln->core.pos, bam_endpos(aln), fwd, rev);
    if(fwd.size() == 0 || rev.size() == 0 || ((aln->core.flag&BAM_FUNMAP) != 0))
      continue;
    std::cout << bam_get_qname(aln) << std::endl;
    fwd_max_end = get_max_end(fwd);
    rev_min_start = get_min_start(rev);
    cur = get_alleles_per_read(aln, cur, &fwd_max_end, &rev_min_start);
  }
  while(cur->get_prev() != NULL){
    cur = cur->get_prev();
  }
  cur->print_graph();
  delete cur;
  hts_itr_destroy(iter);
  hts_idx_destroy(idx);
  bam_destroy1(aln);
  bam_hdr_destroy(header);
  sam_close(in);
  return 0;
}
