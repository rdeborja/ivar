#include "variants_by_amplicon.h"

var_by_amp::var_by_amp(uint64_t pos){
  this->pos = pos;
  this->prev = NULL;
  this->next = NULL;
}

var_by_amp::~var_by_amp(){

}

allele* var_by_amp::get_allele(std::string nuc, std::string deleted_bases, primer *fwd, primer*rev, int &allele_ind){
  allele_ind = -1;
  for (std::vector<allele*>::iterator it = this->alleles.begin(); it != this->alleles.end(); ++it) {
    if((*it)->nuc.compare(nuc) == 0 && (*it)->deleted_bases.compare(deleted_bases) == 0 && this->fwd_primers.at(it - this->alleles.begin())->get_indice() == fwd->get_indice() && this->rev_primers.at(it - this->alleles.begin())->get_indice() == rev->get_indice()){
      allele_ind = (it - this->alleles.begin());
      return (*it);
    }
  }
  return NULL;
}

void var_by_amp::add_allele(allele* a, primer *fwd, primer *rev){
  this->alleles.push_back(a);
  this->fwd_primers.push_back(fwd);
  this->rev_primers.push_back(rev);
  this->associated_variants.push_back(std::map<uint32_t, std::map<std::string, uint32_t>>());
}

std::vector<primer*> var_by_amp::get_fwd_primers(){
  return this->fwd_primers;
}

std::vector<primer*> var_by_amp::get_rev_primers(){
  return this->rev_primers;
}

uint64_t var_by_amp::get_pos(){
  return this->pos;
}

std::vector<allele*> var_by_amp::get_alleles(){
  return this->alleles;
}

std::vector<std::map<uint32_t, std::map<std::string, uint32_t>>> var_by_amp::get_associated_variants(){
  return this->associated_variants;
}

void var_by_amp::add_next(var_by_amp *n){
  this->next = n;
}

void var_by_amp::add_prev(var_by_amp *p){
  this->prev = p;
}

var_by_amp* var_by_amp::get_node(uint64_t p){
  if(p == this->pos)
    return this;
  var_by_amp *cur = this;
  if(p < this->pos){
    while(p != cur->get_pos()){
      if(cur->prev == NULL)
	return NULL;
      cur = cur->prev;
    }
  } else {
    while(p!=cur->get_pos()){
      if(cur->next == NULL)
	return NULL;
      cur = cur->next;
    }
  }
  return cur;
}

var_by_amp* var_by_amp::get_or_add_node(uint64_t p){
  if(p == this->pos)
    return this;
  var_by_amp *cur = this;
  if(p < this->pos){
    while(p != cur->get_pos()){
      if(cur->prev == NULL){
	var_by_amp* n = new var_by_amp(p);
	cur->add_prev(n);
	n->next = cur;
	return n;
      }
      if(cur->prev->get_pos() < p){ // If prev position missing, insert node
	var_by_amp* n = new var_by_amp(p);
	var_by_amp *tmp = cur->prev;
	cur->add_prev(n);
	n->next = cur;
	n->prev = tmp;
	tmp->add_next(n);
	return n;
      }
      cur = cur->prev;
    }
  } else {
    while(p != cur->get_pos()){
      if(cur->next == NULL){
	var_by_amp* n = new var_by_amp(p);
	cur->add_next(n);
	n->prev = cur;
	return n;
      }
      if(cur->next->get_pos() > p){ // If next position greater than pos, insert node
	var_by_amp* n = new var_by_amp(p);
	var_by_amp *tmp = cur->next;
	cur->add_next(n);
	n->next = tmp;
	n->prev = cur;
	tmp->add_prev(n);
	return n;
      }
      cur = cur->next;
    }
  }
  return cur;
}

allele* var_by_amp::get_or_add_allele(std::string nuc, std::string deleted_bases, primer *fwd, primer *rev){
  int allele_ind;
  allele *a = this->get_allele(nuc, deleted_bases, fwd, rev, allele_ind);
  if(a == NULL){
    a = new allele();
    a->nuc = nuc;
    a->deleted_bases = deleted_bases;
    a->depth = 0;
    a->reverse = 0;
    a->beg = 0;
    a->end = 0;
    a->mean_qual = 0;
    a->tmp_mean_qual = 0;
    this->add_allele(a, fwd, rev);
  }
  return a;
}

void var_by_amp::print_graph(bool recurse = true){
  for (std::vector<allele*>::iterator it = alleles.begin(); it != alleles.end(); ++it) {
    std::cout << this->pos << this->delim;
    std::cout << this->fwd_primers.at(it-alleles.begin())->get_name() << this->delim << this->rev_primers.at(it-alleles.begin())->get_name() << this->delim;
    std::cout << (*it)->nuc << this->delim;
    if((*it)->deleted_bases.size() > 0) {
      std::cout << (*it)->deleted_bases.size() << this->delim;
    }
    std::cout << (*it)->depth << std::endl;
    // std::cout << "Associated variants" << std::endl;
    // std::map<uint32_t, std::map<std::string, uint32_t>> m = this->get_associated_variants().at((it - alleles.begin()));
    // for(std::map<uint32_t, std::map<std::string, uint32_t>>::iterator m_it1 = m.begin(); m_it1 != m.end(); m_it1++){
    //   for(std::map<std::string, uint32_t>::iterator m_it2 = m_it1->second.begin(); m_it2 != m_it1->second.end(); m_it2++){
    // 	std::cout << m_it1->first << this->delim << m_it2->first << this->delim << m_it2->second << std::endl;
    //   }
    // }
    std::cout << "-----" << std::endl;
  }
  if(recurse){
    if(this->next != NULL){
      this->next->print_graph();
    }
  }
}

void var_by_amp::get_linked_variants_on_amplicon(int allele_indice){
  primer *fwd = this->get_fwd_primers().at(allele_indice);
  primer *rev = this->get_rev_primers().at(allele_indice);
  var_by_amp *tmp = this;
  std::vector<allele*> alleles;
  allele *a = this->get_alleles().at(allele_indice);
  double freq, cur_freq = a->depth/(double) this->get_depth();
  uint32_t depth;
  for (uint i = fwd->get_start(); i < rev->get_end() + 1; ++i) {
    tmp = tmp->get_node(i);
    if(tmp == NULL){
      tmp = this;
      continue;
    }
    alleles = tmp->get_alleles();
    depth = tmp->get_depth();
    for (std::vector<allele*>::iterator it = alleles.begin(); it != alleles.end(); ++it) {
      if(*fwd == *(tmp->get_fwd_primers().at(it-alleles.begin())) && *rev == *(tmp->get_rev_primers().at(it-alleles.begin()))){
	continue;
      }
      freq = (*it)->depth/(double) depth;
      if(freq < cur_freq/2 || freq > cur_freq * 2) // Check within a fold change for noq
	continue;
      std::cout << this->pos << "\t" << a->nuc << "\t" << cur_freq << "\t" << fwd->get_name() << "\t" << rev->get_name() << "\t" << tmp->pos << "\t" << (*it)->nuc << "\t" << freq << std::endl;
      }
  }
}

void var_by_amp::print_graph(double min_freq){
  uint32_t total_depth = this->get_depth();
  double freq;
  for (std::vector<allele*>::iterator it = this->alleles.begin(); it != this->alleles.end(); ++it) {
    freq = (double) (*it)->depth/(double) total_depth;
    if(freq < min_freq)
      continue;
    std::cout << this->pos << this->delim;
    std::cout << this->fwd_primers.at(it-alleles.begin())->get_name() << this->delim << this->rev_primers.at(it-alleles.begin())->get_name() << this->delim;
    std::cout << (*it)->nuc << this->delim;
    if((*it)->deleted_bases.size() > 0){
      std::cout << (*it)->deleted_bases.size() << this->delim;
    }
    std::cout << (*it)->depth << std::endl;
  }
}

var_by_amp* var_by_amp::get_prev(){
  return this->prev;
}

var_by_amp* var_by_amp::get_next(){
  return this->next;
}

uint32_t var_by_amp::get_depth(){
  uint32_t depth = 0;
  std::vector<allele*> a = this->get_alleles();
  for (std::vector<allele*>::iterator it = a.begin(); it != a.end(); ++it) {
    depth += (*it)->depth;
  }
  return depth;
}

uint get_count_of_alleles(allele *a, std::vector<allele*> _alleles){
  int count = 0;
  for (std::vector<allele*>::iterator it = _alleles.begin(); it != _alleles.end(); ++it) {
    if(**it == *a)
      count += 1;
  }
  return count;
}

int find_allele(allele *a, std::vector<allele*> _alleles){
  for (std::vector<allele*>::iterator it = _alleles.begin(); it != _alleles.end(); ++it) {
    if(**it == *a)
      return it-_alleles.begin();
  }
  return -1;
}

void var_by_amp::get_unique_primers(std::vector<primer*> &uniq_fwd, std::vector<primer*> &uniq_rev){
  std::vector<primer*> pfwd = this->get_fwd_primers();
  std::vector<primer*> prev = this->get_rev_primers();
  std::vector<primer*>::iterator fwdit;
  static uint16_t fwd_indice;
  for (uint i = 0; i< pfwd.size();i++) {
    fwd_indice = pfwd.at(i)->get_indice();
    fwdit = std::find_if(uniq_fwd.begin(), uniq_fwd.end(), [] (const primer* s) {
	return (fwd_indice == s->get_indice());
      });
    if(fwdit == uniq_fwd.end() || uniq_rev.at(fwdit - uniq_fwd.begin())->get_indice() != prev.at(i)->get_indice()){
      uniq_fwd.push_back(pfwd.at(i));
      uniq_rev.push_back(prev.at(i));
    }
  }
}

int var_by_amp::get_num_unique_primers(){
  std::vector<primer*> uniq_fwd;
  std::vector<primer*> uniq_rev;
  this->get_unique_primers(uniq_fwd, uniq_rev);
  return uniq_fwd.size();
}

void var_by_amp::get_distinct_variants_amp(double min_freq, std::vector<allele*> &unique_alleles, std::vector<uint32_t> &counts, uint &unique_primers){
  unique_alleles.clear();
  counts.clear();
  unique_primers = this->get_num_unique_primers();
  std::vector<allele*> a = this->get_alleles();
  std::vector<uint32_t> counts_per_amplicon;
  uint32_t count, total_depth = this->get_depth();
  double error_rate, freq;
  for (uint i = 0; i < a.size(); ++i) {
    error_rate = pow(10, -1*((double)(a.at(i)->mean_qual)/(double)10));
    freq = ((double)(a.at(i)->depth)/(double) total_depth);
    if((freq <= min_freq) || freq <= error_rate)
      continue;
    if(find_allele(a.at(i), unique_alleles) != -1)
      continue;
    count = get_count_of_alleles(a.at(i), a);
    counts.push_back(count);
    unique_alleles.push_back(a.at(i));
  }
}

uint32_t** var_by_amp::get_contingency_table(uint32_t pos, primer *fwd, primer *rev, double min_freq, int &nrow, int &ncol, std::vector<std::string> &row_alleles, std::vector<std::string> &col_alleles){
  row_alleles.clear();
  col_alleles.clear();
  uint32_t** ctable;
  ctable = new uint32_t*[100];
  for (int i = 0; i < 100; ++i) {
    ctable[i] = new uint32_t[100];
    for (int j = 0; j < 100; ++j) {
      ctable[i][j] = 0;
    }
  }
  std::vector<primer*> fwd_primers = this->get_fwd_primers();
  std::vector<primer*> rev_primers = this->get_rev_primers();
  std::vector<std::string> cols;
  std::vector<std::string>::iterator cit;
  std::vector<allele*> _alleles = this->get_alleles();
  int current_col = 0;
  nrow = 0;
  ncol = 0;
  /*
    Contingency table
    | Current pos\Other pos  |  |  |  |
    | Allele 1               |  |  |  |
    | Allele 2               |  |  |  |
    | Allele 3               |  |  |  |
   */
  std::map<uint32_t, std::map<std::string, uint32_t>> m;
  uint32_t total_depth = this->get_depth();
  for (std::vector<allele*>::iterator it = alleles.begin(); it != alleles.end(); ++it) {
    if(fwd_primers.at(it-alleles.begin())->get_indice() != fwd->get_indice() || rev_primers.at(it-alleles.begin())->get_indice() != rev->get_indice())
      continue;
    if(((*it)->depth/(double) total_depth) <= min_freq)
      continue;
    m = this->get_associated_variants().at((it - alleles.begin()));
    if(m.find(pos) == m.end()){
      continue;
    }
    row_alleles.push_back((*it)->nuc);
    for(std::map<std::string, uint32_t>::iterator mit = m[pos].begin(); mit != m[pos].end(); mit++){
      cit = std::find(cols.begin(), cols.end(), mit->first);
      if(cit == cols.end()){ // New col
	ncol++;
	current_col = ncol-1;
	cols.push_back(mit->first);
	col_alleles.push_back(mit->first);
      } else {
	current_col = (cit - cols.begin());
      }
      ctable[nrow][current_col] = mit->second;
    }
    nrow++;
  }
  return ctable;
}

void var_by_amp::get_pos_extent(primer *fwd, primer *rev, uint32_t *extent){
  std::vector<allele*> _alleles = this->get_alleles();
  extent[0] = std::numeric_limits<uint32_t>::max();
  extent[1] = 0;
  std::map<uint32_t, std::map<std::string, uint32_t>> m;
  for (std::vector<allele*>::iterator it = _alleles.begin(); it != _alleles.end(); ++it) {
    if(fwd_primers.at(it-_alleles.begin())->get_indice() != fwd->get_indice() || rev_primers.at(it-_alleles.begin())->get_indice() != rev->get_indice())
      continue;
    m = this->get_associated_variants().at((it - _alleles.begin()));
    for(std::map<uint32_t, std::map<std::string, uint32_t>>::iterator mit = m.begin(); mit != m.end(); mit++){
      if(extent[0] > mit->first){
	extent[0] = mit->first;
      }
      if(extent[1] < mit->first){
	extent[1] = mit->first;
      }
    }
  }
}

void var_by_amp::print_linked_variants(double min_freq){
  std::string nuc, deleted_bases;
  int nrow, ncol;
  uint32_t **ctable;
  double *res;
  std::map<uint32_t, std::map<std::string, uint32_t>> m;
  std::vector<primer*> uniq_fwd;
  std::vector<primer*> uniq_rev;
  this->get_unique_primers(uniq_fwd, uniq_rev);
  uint32_t extent[2];
  std::vector<primer*>::iterator fit = uniq_fwd.begin();
  std::vector<primer*>::iterator rit = uniq_rev.begin();
  std::vector<std::string> row_alleles, col_alleles;
  for (; fit != uniq_fwd.end(); ++fit, ++rit) {
    this->get_pos_extent(*fit, *rit, extent);
    for (uint32_t i = extent[0]; i < extent[1]+1; ++i) {
      ctable = this->get_contingency_table(i, *fit, *rit, min_freq, nrow, ncol, row_alleles, col_alleles);
      if(ctable == NULL)
	continue;
      if(ncol == 1 || nrow == 1)
	continue;
      res = chisqr_goodness_of_fit(ctable, nrow, ncol);
      if(res[0] > 0.05)
	continue;
      std::cout << this->get_pos() << this->delim << i << this->delim << res[0] << this->delim << res[1] << this->delim << (*fit)->get_name() << this->delim << (*rit)->get_name()  << std::endl;
      for (int i = 0; i < nrow+1; ++i) {
	for (int j = 0; j < ncol; ++j) {
	  if(i == 0){
	    if(j==0)
	      std::cout << "\t";
	    std::cout << col_alleles.at(j) << "\t";
	  } else {
	    if(j == 0)
	      std::cout << row_alleles.at(i-1) << "\t";
	    std::cout << ctable[i-1][j] << "\t";
	  }
	}
	std::cout << std::endl;
      }
      delete[] ctable;
    } 
  }
}

int var_by_amp::add_associated_variants(uint32_t pos, allele *aa, allele *a, primer *fwd, primer *rev){
  int ind;
  get_allele(a->nuc, a->deleted_bases, fwd, rev, ind);
  if (ind == -1){
    std::cout << "Allele not found in vector" << std::endl;
    return -1;
  }
  std::string nuc = (aa->deleted_bases.size() == 0) ? aa->nuc : aa->nuc + "-" + aa->deleted_bases;
  if(associated_variants.at(ind).find(pos) == associated_variants.at(ind).end()){
    std::map<std::string, uint32_t> m;
    m[nuc] = 1;
    this->associated_variants.at(ind)[pos] = m;
  } else {
    this->associated_variants.at(ind)[pos][nuc] += 1;
  }
  return 0;
}

std::vector<int> var_by_amp::get_alleles_above_freq(double min_freq){
  double freq;
  std::vector<int> v;
  uint32_t total_depth = this->get_depth();
  std::vector<allele*> a = this->get_alleles();
  double error_rate;
  for (std::vector<allele*>::iterator it = a.begin(); it != a.end(); ++it) {
    freq = ((double)(*it)->depth)/(double)total_depth;
    error_rate = pow(10, -1*((double)((*it)->mean_qual)/(double)10));
    if(freq >= min_freq && freq >= error_rate)
      v.push_back(it-a.begin());
  }
  return v;
}

double compute_critical_value(uint32_t** ctable, int nrow, int ncol){
  double exp[16][16];
  int i, j;
  double cv = 0;
  uint32_t row_total[16] = {0}, col_total[16] = {0}, total = 0;
  for (i = 0; i < nrow; ++i) {
    for (j = 0; j < ncol; ++j) {
      col_total[j] += ctable[i][j];
      row_total[i] += ctable[i][j];
    }
  }
  total = std::accumulate(row_total, row_total+nrow, total);
  for (i = 0; i < nrow; ++i) {
    for (j = 0; j < ncol; ++j) {
      exp[i][j] = (col_total[j] * row_total[i])/(double) total;
    }
  }
  for (i = 0; i < nrow; ++i) {
    for (j = 0; j < ncol; ++j) {
      cv += pow(ctable[i][j] - exp[i][j], 2)/exp[i][j];
    }
  }
  return cv;
}

double chisqr(uint32_t dof, double cv){
  if(cv < 0 || dof < 1){
    return 0;
  }
  double k = ((double) dof) * 0.5;
  double x = cv * 0.5;
  if(dof == 2){
    return exp(-1.0 * x);
  }
  double pvalue = igf(k, x);
  if(isnan(pvalue) || isinf(pvalue) || pvalue <= 1e-8){
    return 1e-14;
  }
  pvalue /= tgamma(k);
  return (1.0 - pvalue);
}

double igf(double s, double z){
if(z < 0.0){
  return 0;
 }
 double sc = (1.0 / s);
 sc *= pow(z, s);
 sc *= exp(-z);
 double sum = 1.0;
 double num = 1.0;
 double den = 1.0; 
 for(int i = 0; i < 200; i++){
   num *= z;
   s++;
   den *= s;
   sum += (num / den);
 }
 return sum * sc;
}

double* chisqr_goodness_of_fit(uint32_t **ctable, int nrow, int ncol){
  static double res[2] = {0};
  uint32_t dof = (nrow-1) * (ncol-1);
  double cv = compute_critical_value(ctable, nrow, ncol);
  res[0] = chisqr(dof, cv);
  res[1] = cv;
  return res;
}
