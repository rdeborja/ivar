#include<iostream>
#include "../src/call_consensus_pileup.h"
#include "../src/allele_functions.h"

int main() {
  int num_tests = 8;
  allele a1 = {
    "A",
    3,
    0,
    30,
    0,
    0,
    0,
    ""
  };
  allele ai1 = {
    "AT",
    3,
    0,
    30,
    0,
    0,
    0,
    ""
  };
  allele a2 = {
    "T",
    1,
    10,
    30,
    0,
    0,
    0,
    ""
  };
  allele a3 = {
    "C",
    1,
    10,
    30,
    0,
    0,
    0,
    ""
  };
  allele a4 = {
    "G",
    1,
    10,
    30,
    0,
    0,
    0,
    ""
  };
  allele a5 = {
    "*",
    1,
    10,
    30,
    0,
    0,
    0,
    ""
  };
  allele a6 = {
    "GA",
    2,
    10,
    30,
    0,
    0,
    0,
    ""
  };
  int success = 0;
  allele arr[] = {a1, ai1,a2,a3,a4,a5};
  ret_t s;
  std::vector<allele> ad(arr, arr+sizeof(arr)/sizeof(allele));
  int size = sizeof(arr)/sizeof(allele);
  for(int i = 0;i<size;i++){
    ad.at(i) = arr[i];
  }
  // A, AT
  s = get_consensus_allele(ad,20,.6, 'N');
  std::cout << s.nuc << ": " << s.q << std::endl;
  success += (s.nuc.compare("AT") == 0) ? 1: 0;
  success += (s.q.compare("?5") == 0) ? 1 : 0;
  // A, AT, T, C, G -> NT
  s = get_consensus_allele(ad,20,.7, 'N');
  std::cout << s.nuc << ": " << s.q << std::endl;
  success += (s.nuc.compare("NT") == 0) ? 1: 0;
  success += (s.q.compare("?5") == 0) ? 1 : 0;

  // A, AT, GA
  ad.erase(ad.begin() + 4, ad.begin()+5);
  ad.push_back(a6);
  s = get_consensus_allele(ad,20,.7, 'N');
  std::cout << s.nuc << ": " << s.q << std::endl;
  success += (s.nuc.compare("RW") == 0) ? 1: 0;
  success += (s.q.compare("?5") == 0) ? 1 : 0;

  // *
  allele a7 = {
    "*",
    10,
    0,
    30,
    0,
    0,
    0,
    ""
  };

  ad.push_back(a7);
  s = get_consensus_allele(ad,20,0, 'N');
  std::cout << s.nuc << ": " << s.q << std::endl;
  success += s.nuc.empty() ? 1: 0;
  success += s.q.empty() ? 1 : 0;  
  return (success == num_tests) ? 0 : -1;
}
