#include "call_variants.h"

std::vector<allele>::iterator get_ref_allele(std::vector<allele> &ad, char ref){
  for(std::vector<allele>::iterator it = ad.begin(); it != ad.end(); ++it) {
    if(it->nuc[0] == ref && !is_indel(*it))
      return it;
  }
  return ad.end();
}

double* get_frequency_depth(allele a, uint32_t pos_depth, uint32_t total_depth){ // pos_depth is ungapped depth after passing quality. total_depth is depth without quality filtering for indels.
  double *val = new double[2];
  if(is_indel(a)){	// For insertions and deletions use depth discarding quality
    val[0] = a.depth/(double)total_depth;
    val[1] = total_depth;
    return val;
  }
  val[0] = a.depth/(double)pos_depth;
  val[1] = pos_depth;
  return val;
}

int call_variants_from_plup(std::istream &cin, std::string out_file, uint8_t min_qual, double min_threshold, uint8_t min_depth, std::string ref_path, std::string gff_path){
  std::string line, cell, bases, qualities, region;
  variants_tsv_writer *tsv_writer = new variants_tsv_writer(out_file+".tsv");
  ref_antd refantd(ref_path, gff_path);
  tsv_writer->add_ref(&refantd);
  std::ostringstream out_str;
  int ctr = 0;
  int64_t pos = 0;
  uint32_t mdepth = 0, pdepth = 0; // mpdepth for mpileup depth and pdeth for ungapped depth at position
  double *freq_depth;
  std::stringstream line_stream;
  char ref;
  std::vector<allele> ad;
  std::vector<allele>::iterator ref_it;
  while (std::getline(cin, line)){
    line_stream << line;
    ctr = 0;
    ref = 'N';
    while(std::getline(line_stream,cell,'\t')){
      switch(ctr){
      case 0:
	region = cell;
	break;
      case 1:
	pos = stoi(cell);
	break;
      case 2:
	ref = refantd.get_base(pos, region);
	ref = (ref == 0) ? cell[0] : ref; // If ref does not exist then use from mpileup
	break;
      case 3:
	mdepth = stoi(cell);
	break;
      case 4:
	bases = cell;
	break;
      case 5:
	qualities = cell;
	break;
      case 6:
	break;
      }
      ctr++;
    }
    if(mdepth < min_depth) {	// Check for minimum depth
      line_stream.clear();
      continue;
    }
    ad = update_allele_depth(ref, bases, qualities, min_qual);
    if(ad.size() == 0){
      line_stream.clear();
      continue;
    }
    ref_it = get_ref_allele(ad, ref);
    if (ref_it == ad.end()) {	// If ref not present in reads.
      allele a;
      a.nuc = ref;
      a.depth = 0;
      a.reverse = 0;
      a.mean_qual = 0;
      a.deleted_bases = "";
      ad.push_back(a);
      ref_it = ad.end() - 1;
    }
    // Get ungapped coverage
    pdepth = 0;
    for(std::vector<allele>::iterator it = ad.begin(); it != ad.end(); ++it) {
      if(it->nuc[0]=='*')
	continue;
      pdepth += it->depth;
    }
    for(std::vector<allele>::iterator it = ad.begin(); it != ad.end(); ++it) {
      if((*it == *ref_it) || it->nuc[0]=='*')
	continue;
      freq_depth = get_frequency_depth(*it, pdepth, mdepth);
      if(freq_depth[0] < min_threshold)
	continue;
      tsv_writer->write_line(*it, *ref_it, region, pos, freq_depth);
    }
    line_stream.clear();
  }
  delete tsv_writer;
  return 0;
}
