#include "call_consensus_pileup.h"

ret_t get_consensus_allele(std::vector<allele> ad, uint8_t min_qual, double threshold, char gap){
  ret_t t;
  t.nuc = "";
  t.q = "";
  t.nuc=gap;
  t.q = min_qual+33;
  if(ad.size()==0)
    return t;
  uint32_t total_depth = get_total_depth(ad);
  get_alleles_by_threshold(ad, threshold, total_depth);
  total_depth = get_total_depth(ad);
  double q = 0;
  std::string nuc;
  std::string qual;
  char n;
  uint8_t max_len = get_longest_insertion(ad);
  int i =0;
  std::vector<allele>::iterator it;
  for (i = 0; i < max_len; ++i) {
    n = 0;
    q = 0;
    for (it = ad.begin(); it != ad.end(); ++it) {
      if(n == 0){
	n = it->nuc.length() > i ? it->nuc[i] : 0;
      } else {
	n = it->nuc.length() > i ? gt2iupac(n, it->nuc[i]) : n;
      }
      if(i == 0){		// Get qual from first bases
	q = q + (uint8_t)(it->mean_qual) * it->depth;
      } else {
	q = min_qual;
      }
    }
    if(n!='*'){
      q = q/total_depth;
      q += 0.5;
      qual += (((uint8_t)q) + 33);
      nuc += n;
    }
  }
  t.nuc = nuc;
  t.q = qual;
  return t;
}

int call_consensus_from_plup(std::istream &cin, std::string out_file, uint8_t min_qual, double threshold, uint8_t min_depth, char gap, bool min_coverage_flag, std::string ref_path){
  std::string line, cell;
  std::ofstream fout((out_file+".fa").c_str());
  std::ofstream tmp_qout((out_file+".qual.txt").c_str());
  char *o = new char[out_file.length() + 1];
  vcf_writer *vw = NULL;
  strcpy(o, out_file.c_str());
  fout << ">Consensus_" << basename(o) << "_threshold_" << threshold << "_quality_" << (uint16_t) min_qual  <<std::endl;
  int ctr = 0, mdepth = 0;
  uint32_t prev_pos = 0, pos = 0;
  std::stringstream lineStream;
  char ref;
  std::string bases;
  std::string qualities;
  std::vector<allele> ad;
  uint32_t bases_zero_depth = 0, bases_min_depth = 0, total_bases = 0;
  std::string region;
  while (std::getline(cin, line)){
    lineStream << line;
    ctr = 0;
    ref = 'N';
    while(std::getline(lineStream,cell,'\t')){
      switch(ctr){
      case 0:
	region = cell;
	break;
      case 1:
	pos = stoi(cell);
	break;
      case 2:
	ref = cell[0];
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
    total_bases++;
    if(prev_pos == 0)		// No -/N before alignment starts
      prev_pos = pos;
    if((pos > prev_pos && min_coverage_flag)){
      fout << std::string((pos - prev_pos) - 1, gap);
      tmp_qout << std::string((pos - prev_pos) - 1, '!'); // ! represents 0 quality score.
    }
    if(vw == NULL){ // Note: for now consensus working only for 1 region
      if(!ref_path.empty())
	vw = new vcf_writer('b', out_file+".bcf", region, basename(o), ref_path);
    }
    ret_t t;
    if(mdepth >= min_depth){
      ad = update_allele_depth(ref, bases, qualities, min_qual);
      t = get_consensus_allele(ad, min_qual, threshold, gap);
      if(vw!=NULL)
	vw->write_record(pos, ad, region, ref);
      fout << t.nuc;
      tmp_qout << t.q;
    } else{
      bases_min_depth += 1;
      if (mdepth == 0)
	bases_zero_depth += 1;
      if(min_coverage_flag){
      if(vw!=NULL)
	vw->write_record_below_threshold(pos, region, ref);
	fout << gap;
	tmp_qout << '!';
      }
    }
    lineStream.clear();
    ad.clear();
    prev_pos = pos;
  }
  fout << "\n";			// Add new line character after end of sequence
  tmp_qout << "\n";
  tmp_qout.close();
  fout.close();
  std::cout << "Reference length: " << total_bases << std::endl;
  std::cout << "Positions with 0 depth: " << bases_zero_depth << std::endl;
  std::cout << "Positions with depth below " <<(unsigned) min_depth << ": " << bases_min_depth << std::endl;
  delete vw;
  delete [] o;
  return 0;
}
