#include "variants_tsv_writer.h"

const float sig_level = 0.01;

variants_tsv_writer::variants_tsv_writer(std::string path){
  this->fout.open((path+".tsv").c_str());
  if(this->fout.fail()){
    std::cout << "Failed to open file " << path << std::endl;
  }
  this->fout << this->hdr << std::endl;
}

int variants_tsv_writer::write_line(allele a, allele ref, std::string region, uint64_t pos, double *freq_depth){
  std::ostringstream out_str;
  double err, pval_left, pval_right, pval_twotailed;
  out_str << region << this->delim;
  out_str << pos << this->delim;
  if(!(a.deleted_bases.empty())){
    out_str << ref.nuc << a.deleted_bases << "\t";
  } else {
    out_str << ref.nuc << "\t";
  }
  out_str << a.nuc << "\t";
  out_str << ref.depth << "\t";
  out_str << ref.reverse << "\t";
  out_str << (uint8_t) ref.mean_qual << "\t";
  out_str << a.depth << "\t";
  out_str << a.reverse << "\t";
  out_str << (uint8_t) a.mean_qual << "\t";
  out_str << freq_depth[0] << "\t";
  out_str << freq_depth[1] << "\t";
  err = pow(10, ( -1 * (a.mean_qual)/10));
  /*
        | Var   | Ref      |
    Exp | Error | Err free |
    Obs | AD    | RD       |
  */
  kt_fisher_exact((err * freq_depth[1]), (1-err) * freq_depth[1], a.depth, ref.depth, &pval_left, &pval_right, &pval_twotailed);
  out_str << pval_left << "\t";
  if(pval_left <= sig_level){
    out_str << "TRUE" << "\t";
  } else {
    out_str << "FALSE" << "\t";
  }
  // if(!is_indel(*it)){ // Ignore indels for aa translation
  //   refantd.codon_aa_stream(region, out_str, fout, pos, it->nuc[0]);
  // } else {
  //   fout << out_str.str() << "NA\tNA\tNA\tNA\tNA" << std::endl;
  // }  
  out_str.str("");
  out_str.clear();
  return 0;
}
