#include "variants_tsv_writer.h"

const float sig_level = 0.01;

variants_tsv_writer::variants_tsv_writer(std::string path){
  this->fout.open(path.c_str());
  if(this->fout.fail()){
    std::cout << "Failed to open file " << path << std::endl;
  }
}

void variants_tsv_writer::add_ref(ref_antd *ref_a){
  this->refa = ref_a;
  if(this->refa->is_empty()){
    this->fout << this->hdr << std::endl;
  } else {
    this->fout << this->hdr_aa << std::endl;
  }
}

variants_tsv_writer::~variants_tsv_writer(){
  this->fout.close();
}

int variants_tsv_writer::write_line(allele a, allele ref, std::string region, uint64_t pos, double *freq_depth){
  std::ostringstream out_str;
  double err, pval_left, pval_right, pval_twotailed;
  out_str << region << this->delim;
  out_str << pos << this->delim;
  if(!(a.deleted_bases.empty())){
    out_str << ref.nuc << a.deleted_bases << this->delim;
  } else {
    out_str << ref.nuc << this->delim;
  }
  out_str << a.nuc << this->delim;
  out_str << ref.depth << this->delim;
  out_str << ref.reverse << this->delim;
  out_str << (uint16_t) ref.mean_qual << this->delim;
  out_str << a.depth << this->delim;
  out_str << a.reverse << this->delim;
  out_str << (uint16_t) a.mean_qual << this->delim;
  out_str << freq_depth[0] << this->delim;
  out_str << freq_depth[1] << this->delim;
  err = pow(10, ( -1 * (a.mean_qual)/10));
  /*
        | Var   | Ref      |
    Exp | Error | Err free |
    Obs | AD    | RD       |
  */
  kt_fisher_exact((err * freq_depth[1]), (1-err) * freq_depth[1], a.depth, ref.depth, &pval_left, &pval_right, &pval_twotailed);
  out_str << pval_left << this->delim;
  if(pval_left <= sig_level){
    out_str << "TRUE";
  } else {
    out_str << "FALSE";
  }
  if(!this->refa->is_empty()){
    out_str << this->delim;
    if(!is_indel(a)){ // Ignore indels for aa translation
      std::vector<std::vector<std::string>> codon_aa = refa->codon_aa(region, pos, a.nuc[0]);
      for (std::vector<std::vector<std::string>>::iterator it = codon_aa.begin(); it != codon_aa.end(); ++it) {
	fout << out_str.str();
	for (std::vector<std::string>::iterator cit = (*it).begin(); cit != (*it).end(); ++cit) {
	  this->fout << *cit;
	  if (cit < (it->end() - 1))
	    this->fout << this->delim;
	}
	this->fout << std::endl;
      }
    } else {
      this->fout << out_str.str() << "NA\tNA\tNA\tNA\tNA" << std::endl;
    }
  } else {
    this->fout << out_str.str() << std::endl;
  }
  out_str.str("");
  out_str.clear();
  return 0;
}
