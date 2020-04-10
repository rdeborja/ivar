#include "vcf_writer.h"

int vcf_writer::init_header(){
  this->hdr = bcf_hdr_init("w");
  int res;
  res = bcf_hdr_append(this->hdr,("##source=ivar"+(std::string)VERSION).c_str());
  res = bcf_hdr_append(this->hdr,"##ALT=<ID=*,Description=\"Represents allele(s) other than observed.\">");
  res = bcf_hdr_append(this->hdr,"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">");
  res = bcf_hdr_append(this->hdr, "##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Total read depth for each allele (based on minimum quality threshold)\">");
  res = bcf_hdr_append(this->hdr, "##INFO=<ID=ADF,Number=R,Type=Integer,Description=\"Total read depth for each allele on forward strand (based on minimum quality threshold)\">");
  res = bcf_hdr_append(this->hdr, "##INFO=<ID=ADR,Number=R,Type=Integer,Description=\"Total read depth for each allele on reverse strand (based on minimum quality threshold)\">");
  res = bcf_hdr_append(this->hdr, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency for each ALT allele in the same order as listed (estimated from primary data, not called genotypes)\">");
  res = bcf_hdr_append(this->hdr, "##INFO=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
  res = bcf_hdr_append(this->hdr, ("##contig=<ID="+this->region+">").c_str());
  res = bcf_hdr_append(this->hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
  res = bcf_hdr_append(this->hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
  res = bcf_hdr_append(this->hdr, "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Total read depth for each allele (based on minimum quality threshold)\">");
  res = bcf_hdr_append(this->hdr, "##FORMAT=<ID=ADF,Number=R,Type=Integer,Description=\"Total read depth for each allele on forward strand (based on minimum quality threshold)\">");
  res = bcf_hdr_append(this->hdr, "##FORMAT=<ID=ADR,Number=R,Type=Integer,Description=\"Total read depth for each allele on reverse strand (based on minimum quality threshold)\">");
  res = bcf_hdr_append(this->hdr, "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency for each ALT allele in the same order as listed (estimated from primary data, not called genotypes)\">");
  res = bcf_hdr_append(this->hdr, "##FILTER=<ID=fobs,Description=\"Fisher's exact test to check if frequency of the iSNV is significantly higher than the mean error rate at that position\">");
  bcf_hdr_add_sample(this->hdr, this->sample_name.c_str());
  if(res != 0)
    return -1;
  return 0;
}

vcf_writer::~vcf_writer(){
  bcf_hdr_destroy(this->hdr);
  bcf_close(this->file);
}

vcf_writer::vcf_writer(char _mode, std::string fname, std::string region, std::string sample_name, std::string ref_path){
  this->ref = new ref_antd(ref_path);
  this->region = region;
  this->sample_name = sample_name;
  std::string mode = "w";
  mode += _mode;
  this->file = bcf_open(fname.c_str(), mode.c_str());
  this->init_header();
  if(this->hdr == NULL || this->file == NULL){
    std::cout << "Unable to write BCF/VCF file" << std::endl;
  }
  if(bcf_hdr_write(this->file, this->hdr) < 0){
    std::cout << "Unable to write BCF/VCF file" << std::endl;
  }
}

bool compare_allele_depth(const allele &a, const allele &b){
  return b.depth < a.depth;
}

uint32_t get_total_depth(std::vector<allele> vec){
  uint32_t d = 0;
  for (std::vector<allele>::iterator it = vec.begin(); it != vec.end(); ++it) {
    d += it->depth;
  }
  return d;
}

void get_alleles_by_threshold(std::vector<allele> &aalt, double threshold, uint32_t total_depth){
  std::sort(aalt.begin(), aalt.end(), compare_allele_depth); // Sort by depth
  std::vector<allele>::iterator it = aalt.begin();
  double cur_threshold;
  uint32_t cur_depth;
  while(it != aalt.end() && (cur_threshold < threshold || it->depth == (it+1)->depth)){
    cur_depth += it->depth;
    cur_threshold = (double)cur_depth/(double)total_depth;
    it++;
  }
  aalt.erase(it, aalt.end());
}

int* set_genotype(std::vector<allele> aalt, char ref_nuc){
  int *tmp = (int*) malloc((aalt.size())*sizeof(int));
  std::string ref(1, ref_nuc);
  int ctr = 0, ref_offset = 0;
  for(std::vector<allele>::iterator it = aalt.begin(); it != aalt.end(); ++it ){
    ctr = (it - aalt.begin()) + 1 - ref_offset;
    if(it->nuc.compare(ref) == 0){
      ctr = 0;
      ref_offset = 1;
    }
    tmp[it-aalt.begin()] = bcf_gt_unphased(ctr);
  }
  return tmp;
}

int vcf_writer::write_record(uint32_t pos, std::vector<allele> aalt, char ref_nuc, double threshold){
  bcf1_t *rec = bcf_init();
  uint32_t total_depth = get_total_depth(aalt);
  rec->rid = bcf_hdr_name2id(this->hdr, this->region.c_str());
  rec->pos  = pos - 1;		// converts to pos + 1 on write
  bcf_update_id(this->hdr, rec, ".");
  std::string allele_str, ref_str;
  int max_del_len = 0, ctr = 0;
  std::vector<allele>::iterator it;
  get_alleles_by_threshold(aalt, threshold, total_depth);
  int ref_pos = find_ref_in_allele(aalt, ref_nuc);
  if(ref_pos != -1){
    rec->qual = aalt.at(ref_pos).mean_qual;
  }
  for (it = aalt.begin(); it != aalt.end(); ++it) {
    if(it->nuc[0] == ref_nuc && it->nuc.length() == 1) // Skip ref
      continue;
    if(it->nuc[0]=='-'){
      max_del_len = (it->nuc.length() -1 > max_del_len) ? it->nuc.length() : max_del_len;
      allele_str += this->ref->get_base(pos, this->region);
    } else if (it->nuc[0] == '+') {
      allele_str += this->ref->get_base(pos, this->region) + it->nuc.substr(1);
    } else {
      allele_str += it->nuc;
    }
    if (it < aalt.end() - 1)
      allele_str += ",";
  }
  while(ctr < max_del_len + 1){	// By default add one pos from ref
    ref_str += this->ref->get_base(pos + ctr, this->region);
    ctr++;
  }
  allele_str = ref_str + "," + allele_str;
  bcf_update_alleles_str(this->hdr, rec, allele_str.c_str());
  int32_t asize = aalt.size();
  // Set genotype to all alleles required to pass threshold
  int32_t *tmp = set_genotype(aalt, ref_nuc);
  bcf_update_genotypes(this->hdr, rec, tmp, asize);
  free(tmp);
  // FORMAT: add values for all genotypes
  int32_t *tmp_depth, *tmp_qual, *tmp_depth_forward, *tmp_depth_reverse;
  tmp_qual = (int32_t*) malloc((asize) * sizeof(int));
  tmp_depth_forward = (int32_t*)malloc((asize) * sizeof(int));
  tmp_depth_reverse = (int32_t*)malloc((asize) * sizeof(int));
  tmp_depth = (int32_t*)malloc((asize) * sizeof(int));
  float *tmp_freq = (float*)malloc((asize)*sizeof(float));
  for (it = aalt.begin(); it != aalt.end(); ++it) {
    tmp_depth[it - aalt.begin()] = it->depth;
    tmp_qual[it - aalt.begin()] = it->mean_qual;
    tmp_depth_reverse[it - aalt.begin()] = it->reverse;
    tmp_depth_forward[it - aalt.begin()] = it->depth - it->reverse;
    tmp_freq[it - aalt.begin()] = (float) it->depth/(float) total_depth;
  }
  // FORMAT
  bcf_update_format_int32(this->hdr, rec, "AD", tmp_depth, asize);
  bcf_update_format_int32(this->hdr, rec, "GQ", tmp_qual, asize);
  bcf_update_format_int32(this->hdr, rec, "ADF", tmp_depth_forward, asize);
  bcf_update_format_int32(this->hdr, rec, "ADR", tmp_depth_reverse, asize);
  bcf_update_format_float(this->hdr, rec, "AF", tmp_freq, asize);
  // INFO: addvalues only for ALT alleles not REF
  if(ref_pos != -1)
    aalt.erase(aalt.begin()+ref_pos); // Erase ref
  asize = aalt.size();
  free(tmp_depth);
  free(tmp_qual);
  free(tmp_depth_forward);
  free(tmp_depth_reverse);
  free(tmp_freq);
  tmp_qual = (int32_t*) malloc((asize) * sizeof(int));
  tmp_depth_forward = (int32_t*)malloc((asize) * sizeof(int));
  tmp_depth_reverse = (int32_t*)malloc((asize) * sizeof(int));
  tmp_depth = (int32_t*)malloc((asize) * sizeof(int));
  tmp_freq = (float*)malloc((asize)*sizeof(float));
  for (it = aalt.begin(); it != aalt.end(); ++it) {
    tmp_depth[it - aalt.begin()] = it->depth;
    tmp_qual[it - aalt.begin()] = it->mean_qual;
    tmp_depth_reverse[it - aalt.begin()] = it->reverse;
    tmp_depth_forward[it - aalt.begin()] = it->depth - it->reverse;
    tmp_freq[it - aalt.begin()] = (float) it->depth/(float) total_depth;
  }
  int32_t tmpi = 1;
  // NS
  bcf_update_info_int32(this->hdr, rec, "DP", tmp_depth, asize);
  bcf_update_info_int32(this->hdr, rec, "NS", &tmpi, 1);
  tmpi = bcf_hdr_id2int(this->hdr, BCF_DT_ID, "PASS");
  bcf_update_filter(this->hdr, rec, &tmpi, 1);
  bcf_update_info_int32(this->hdr, rec, "AD", tmp_depth, asize);
  bcf_update_info_int32(this->hdr, rec, "ADF", tmp_depth_forward, asize);
  bcf_update_info_int32(this->hdr, rec, "ADR", tmp_depth_reverse, asize);
  bcf_update_info_int32(this->hdr, rec, "GQ", tmp_qual, asize);
  bcf_update_info_float(this->hdr, rec, "AF", tmp_freq, asize);
  int res;
  res = bcf_write1(this->file, this->hdr, rec);
  if(res != 0)
    std::cout << "Unable to write to VCF/BCF file!" << std::endl;
  bcf_clear1(rec);
  free(tmp_depth);
  free(tmp_qual);
  free(tmp_depth_forward);
  free(tmp_depth_reverse);
  free(tmp_freq);
  return res;
}

// int main(int argc, char *argv[]) {
//   vcf_writer *vw = new vcf_writer('b', "./test.vcf", "test", "test_sample", "../data/db/test_ref.fa");
//   allele r = {
//   depth: 8,
//   reverse: 2,
//   nuc: "C",
//   mean_qual:30
//   };
//   allele a1 = {
//   depth: 10,
//   reverse: 6,
//   nuc: "T",
//   mean_qual:25
//   };
//   allele a2 = {
//   depth: 20,
//   reverse: 6,
//   nuc: "A",
//   mean_qual:25
//   };
//   std::vector<allele> a = {r, a1, a2};
//   vw->write_record(5, a, 'C', 0.9);
//   delete vw;
//   return 0;
// }
