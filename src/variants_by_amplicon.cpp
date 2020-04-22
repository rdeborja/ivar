#include "variants_by_amplicon.h"

var_by_amp::var_by_amp(uint64_t pos){
  this->pos = pos;
  this->prev = NULL;
  this->next = NULL;
}

allele* var_by_amp::get_allele(std::string nuc, std::string deleted_bases, primer *fwd, primer*rev){
  for (std::vector<allele*>::iterator it = this->alleles.begin(); it != this->alleles.end(); ++it) {
    if((*it)->nuc.compare(nuc) == 0 && (*it)->deleted_bases.compare(deleted_bases) == 0 && this->fwd_primers.at(it - this->alleles.begin())->get_indice() == fwd->get_indice() && this->rev_primers.at(it - this->alleles.begin())->get_indice() == rev->get_indice())
      return (*it);
  }
  return NULL;
}

void var_by_amp::add_allele(allele* a, primer *fwd, primer *rev){
  alleles.push_back(a);
  this->fwd_primers.push_back(fwd);
  this->rev_primers.push_back(rev);
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

void var_by_amp::add_next(var_by_amp *n){
  this->next = n;
}

void var_by_amp::add_prev(var_by_amp *p){
  this->prev = p;
}

var_by_amp* var_by_amp::get_or_add_node(uint64_t p){
  if(p == this->pos)
    return this;
  var_by_amp *cur = this;
  if(p < this->pos){
    while(p != cur->get_pos()){	// Prev nodes should never be NULL
      if(cur->prev == NULL){
	var_by_amp* n = new var_by_amp(p);
	cur->add_prev(n);
	n->next = cur;
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
      cur = cur->next;
    }
  }
  return cur;
}

allele* var_by_amp::get_or_add_allele(std::string nuc, std::string deleted_bases, primer *fwd, primer *rev){
  allele *a = this->get_allele(nuc, deleted_bases, fwd, rev);
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

void var_by_amp::print_graph(){
  std::cout << this->pos << std::endl;
  for (std::vector<allele*>::iterator it = alleles.begin(); it != alleles.end(); ++it) {
    std::cout << this->fwd_primers.at(it-alleles.begin())->get_name() << ": " << this->rev_primers.at(it-alleles.begin())->get_name() << std::endl;
    std::cout << (*it)->nuc << ": " << (*it)->depth << std::endl;
  }
  std::cout << std::endl;
  if(this->next != NULL){
    this->next->print_graph();
  } else {
    std::cout << "END" << std::endl;
  }
}

var_by_amp* var_by_amp::get_prev(){
  return this->prev;
}

var_by_amp::~var_by_amp(){
  
}
