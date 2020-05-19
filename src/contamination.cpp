#include "contamination.h"

void identify_amp(std::vector<primer> &primers, uint64_t pos, uint64_t end_pos, std::vector<primer*> &fwd_primers, std::vector<primer*> &rev_primers){
  fwd_primers.clear();
  rev_primers.clear();
  primer primer_pair;
  uint16_t p2;
  for (std::vector<primer>::iterator it = primers.begin(); it != primers.end(); ++it) {
    if(it->get_pair_indice()== -1 || it->get_start() > pos)
      continue;
    p2 = it->get_pair_indice();
    if(pos >= it->get_start() && end_pos <= primers.at(p2).get_end()){
      fwd_primers.push_back(&(*it));
      rev_primers.push_back(&(primers.at(p2)));
    }
  }
}

var_by_amp* get_alleles_per_read(bam1_t *aln, var_by_amp *v, primer *fwd, primer *rev, uint8_t min_qual){
  std::vector<allele> all;
  uint ncigar = aln->core.n_cigar;
  uint32_t* cigar = bam_get_cigar(aln);
  uint32_t cig;
  uint i = 0, len, ctr = 0;
  uint32_t pos = aln->core.pos, ql=0;
  uint8_t *seq = bam_get_seq(aln);
  allele *a;
  std::string indel;
  std::vector<allele*> current_read_alleles;
  std::vector<uint32_t> current_read_pos;
  uint8_t *qual = bam_get_qual(aln), q = 0;
  if(v == NULL)
    v = new var_by_amp(pos);
  while (i < ncigar ) {
    cig = bam_cigar_op(cigar[i]);
    len = bam_cigar_oplen(cigar[i]);
    if ((bam_cigar_type(cigar[i]) & 1) && (bam_cigar_type(cigar[i]) & 2)){ // Query and ref consuming
      ctr = 0;
      while(ctr < len -1){
	q = qual[ql+ctr];
	if(q >= min_qual){
	  v = v->get_or_add_node(pos);
	  a = v->get_or_add_allele(std::string(1, seq_nt16_str[bam_seqi(seq, ql+ctr)]), "", fwd, rev);
	  a->mean_qual = ((a->mean_qual * a->depth) + q)/(a->depth + 1);
	  a->depth += 1;
	  current_read_alleles.push_back(a);
	  current_read_pos.push_back(pos);
	}
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
	  if(bam_cigar_op(cigar[i+1]) != BAM_CSOFT_CLIP && q >= min_qual){
	    v = v->get_or_add_node(pos);
	    a = v->get_or_add_allele(indel, "", fwd, rev);
	    a->mean_qual = ((a->mean_qual * a->depth) + q)/(a->depth + 1);
	    a->depth += 1;
	    current_read_alleles.push_back(a);
	    current_read_pos.push_back(pos);
	  }
	  pos++;
	  ql += bam_cigar_oplen(cigar[i+1]);
	  i++;
	} else if (!(bam_cigar_type(cigar[i+1]) & 1) && (bam_cigar_type(cigar[i+1]) & 2)){ // Only ref consuming
	  indel.clear();
	  indel.assign(1, seq_nt16_str[bam_seqi(seq, ql+ctr)]);
	  if(q >= min_qual){
	    v = v->get_or_add_node(pos);
	    a = v->get_or_add_allele(indel, std::string(bam_cigar_oplen(cigar[i+1]), 'N'), fwd, rev);
	    a->mean_qual = ((a->mean_qual * a->depth) + q)/(a->depth + 1);
	    a->depth += 1;
	    current_read_alleles.push_back(a);
	    current_read_pos.push_back(pos);
	  }
	  pos += bam_cigar_oplen(cigar[i+1]) + 1; // Account for incrementing last base of previous cigar
	  ql += len;				  // Add query length for previous cigar
	  i++;
	} else {
	  q = qual[ql+ctr];
	  if(q >= min_qual){
	    v = v->get_or_add_node(pos);
	    a = v->get_or_add_allele(std::string(1, seq_nt16_str[bam_seqi(seq, ql+ctr)]), "", fwd, rev);
	    a->mean_qual = ((a->mean_qual * a->depth) + q)/(a->depth + 1);
	    a->depth += 1;
	    current_read_alleles.push_back(a);
	    current_read_pos.push_back(pos);
	  }
	  pos++;
	  ql += len;
	}
      } else {
	q = qual[ql+ctr];
	if(q >= min_qual){
	  v = v->get_or_add_node(pos);
	  a = v->get_or_add_allele(std::string(1, seq_nt16_str[bam_seqi(seq, ql+ctr)]), "", fwd, rev);
	  a->mean_qual = ((a->mean_qual * a->depth) + q)/(a->depth + 1);
	  a->depth += 1;
	  current_read_alleles.push_back(a);
	  current_read_pos.push_back(pos);
	}
	pos++;
	ql+=len;
      }
    } else if(bam_cigar_type(cigar[i]) & 1){
      ql += len;
    } else if(bam_cigar_type(cigar[i]) & 2){
      pos += len;
    }
    i++;
  }
  if(aln->core.pos <= 3988 && 3988 <= bam_endpos(aln)){
    uint32_t current_pos;
    for (uint i = 0; i < current_read_alleles.size(); ++i) {
      current_pos = current_read_pos.at(i);
      v = v->get_node(current_pos);
      for (uint j = 0; j < current_read_alleles.size(); ++j) {
	if(j==i)
	  continue;
	pos = current_read_pos.at(j);
	v->add_associated_variants(pos, current_read_alleles.at(j), current_read_alleles.at(i), fwd, rev);
      }
    }
  }
  return v;
}

int identify_contamination(std::string bed, std::string bam, std::string pairs_file, std::string prefix, uint8_t min_qual){
  prefix += ".bam";
  BGZF *out = bgzf_open(prefix.c_str(), "w");
  std::vector<primer> primers = populate_from_file(bed);
  populate_pair_indices(primers, pairs_file);
  samFile *in = hts_open(bam.c_str(), "r");
  if(in == NULL) {
    std::cout << ("Unable to open BAM file.") << std::endl;
    return -1;
  }
  //Load the index
  hts_idx_t *idx = sam_index_load(in, bam.c_str());
  if(idx == NULL) {
    std::cout << "Building BAM index" << std::endl;
    if(sam_index_build2(bam.c_str(), 0, 0)< 0){
      std::cout << ("Unable to open or build BAM index.") << std::endl;
      return -1;
    } else {
      idx = sam_index_load(in, bam.c_str());
    }
  }
  bam_hdr_t *header = sam_hdr_read(in);
  if(header == NULL) {
    sam_close(in);
    std::cout << "Unable to open BAM header." << std::endl;
  }
  if(bam_hdr_write(out, header) < 0){
    std::cout << "Unable to write BAM header to path." << std::endl;
    sam_close(in);
    return -1;
  }  
  hts_itr_t *iter = NULL;
  std::vector<primer*> fwd, rev;
  primer *fwd_max_end, *rev_min_start;
  iter  = sam_itr_querys(idx, header, header->target_name[0]);
  var_by_amp *cur= NULL;
  bam1_t *aln = bam_init1();
  double min_freq = 0.03;
  std::cout << "Reads overlapping multiple amplicons: " << std::endl;
  std::vector<std::vector<uint32_t>> overlap = get_overlap_pos(primers);
  while(sam_itr_next(in, iter, aln) >= 0) {
    identify_amp(primers, aln->core.pos, bam_endpos(aln) - 1, fwd, rev);
    if(fwd.size() == 0 || rev.size() == 0 || ((aln->core.flag&BAM_FUNMAP) != 0)){
      // std::cout  << bam_get_qname(aln) << ((bam_is_rev(aln)) ? "/2" : "/1") << "\t" << aln->core.isize << std::endl;
      continue;
    }
    if(bam_write1(out, aln) < 0){
      std::cout << "Not able to write to BAM" << std::endl;
      hts_itr_destroy(iter);
      hts_idx_destroy(idx);
      bam_destroy1(aln);
      bam_hdr_destroy(header);
      sam_close(in);
      bgzf_close(out);
      return -1;
    }
    // fwd_max_end = get_max_end_(fwd);
    // rev_min_start = get_min_start_(rev);
    fwd_max_end = fwd.at(0);
    rev_min_start = rev.at(0);
    cur = get_alleles_per_read(aln, cur, fwd_max_end, rev_min_start, min_qual);
  }
  if(cur == NULL){
    std::cout << "No reads from amplicons found" << std::endl;
  } else {
    while(cur->get_prev() != NULL){
      cur = cur->get_prev();
    }
    // Check for variants in primer positions
    var_by_amp *tmp;
    std::vector<int> allele_indices;
    for (std::vector<primer>::iterator it=primers.begin(); it != primers.end(); ++it) {
      for (uint32_t i = it->get_start(); i < it->get_end()+1; ++i) {
	tmp = cur->get_node(i);
	if (tmp== NULL)
	  continue;
	allele_indices = tmp->get_alleles_above_freq(min_freq);
	if(allele_indices.size() > 1){
	  tmp->print_graph(min_freq);
	  for (std::vector<int>::iterator ait = allele_indices.begin(); ait != allele_indices.end();  ++ait) {
	    tmp->get_linked_variants_on_amplicon(*ait);
	  }
	}
      }
    }
    std::cout << "Overlap" << std::endl << std::endl;
    // Check for variants in overlap position
    std::vector<allele*> unique_alleles;
    std::vector<uint32_t> counts;
    uint unique_primer_count = 0;
    for (std::vector<std::vector<uint32_t>>::iterator it = overlap.begin(); it!=overlap.end(); ++it) {
      for (uint32_t i = (*it).at(0); i < (*it).at(1)+1; ++i) {
	tmp = cur->get_node(i);
	if (tmp== NULL)
	  continue;
	tmp->get_distinct_variants_amp(min_freq, unique_alleles, counts, unique_primer_count);
	for (uint j = 0; j < unique_alleles.size(); ++j) {
	  if(counts.at(j) == unique_primer_count || tmp->get_fwd_primers().size() == 1)
	    continue;
	  tmp->print_graph(false);
	  tmp->print_linked_variants(min_freq);
	  std::cout << "O " << "\t"
		    << i << "\t"
		    << unique_alleles.at(j)->nuc
		    << ((!unique_alleles.at(j)->deleted_bases.empty()) ? "-"+unique_alleles.at(j)->deleted_bases : "") << "\t"
		    << (unique_alleles.at(j)->depth/(double) tmp->get_depth()) << "\t"
		    << counts.at(j) << "\t"
		    << unique_primer_count << "\t"
		    << std::endl;
	}
      }
    }
    std::cout << "---" << std::endl << std::endl;
    cur = cur->get_node(3988);
    // cur->print_graph(false);
    cur->print_linked_variants(min_freq);
    // cur->print_graph(true);
    delete cur; 
  }
  hts_itr_destroy(iter);
  hts_idx_destroy(idx);
  bam_destroy1(aln);
  bam_hdr_destroy(header);
  sam_close(in);
  bgzf_close(out);
  return 0;
}
