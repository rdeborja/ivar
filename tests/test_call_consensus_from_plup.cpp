#include<iostream>
#include "../src/call_consensus_pileup.h"
#include "../src/allele_functions.h"

int main() {
  int num_tests = 11;
  allele a1 = {
    "A",
    40,
    0,
    30,
    0,
    0,
    0,
    ""
  };
  allele a2 = {
    "T",
    40,
    10,
    30,
    0,
    0,
    0,
    ""
  };
  allele a3 = {
    "G",
    40,
    10,
    10,
    0,
    0,
    0,
    ""
  };
  allele a4 = {
    "AT",
    30,
    10,
    10,
    0,
    0,
    0,
    ""
  };
  allele a5 = {
    "AAG",
    30,
    10,
    10,
    0,
    0,
    0,
    ""
  };
  allele a6 = {
    "AAC",
    40,
    4,
    10,
    0,
    0,
    0,
    ""
  };
  allele a7 = {
    "AWG",
    40,
    2,
    10,
    0,
    0,
    0,
    ""
  };
  allele a8 = {
    "CCT",
    40,
    2,
    10,
    0,
    0,
    0,
    ""
  };
  int success = 0;
  std::vector<allele> ad = {a1,a2,a3,a4,a5};
  uint8_t max_len = get_longest_insertion(ad);
  success += (max_len == 3) ? 1 : 0;
  ret_t s;
  s = get_consensus_allele(ad, 0, 0, 'N');
  // A, T, G -> D
  success += (s.nuc.compare("D") == 0) ? 1: 0;
  // (40 * 30 + 40 * 10 + 40 * 10)/120
  success += (s.q.compare("8") == 0) ? 1 : 0;
  std::cout << s.nuc << std::endl;
  std::cout << s.q << std::endl;
  
  ad.push_back(a6);
  s = get_consensus_allele(ad, 0, 0, 'N');
  // A, T, G, AAC -> DAC
  success += (s.nuc.compare("DAC") == 0) ? 1: 0;
  // (40 * 30 + 40 * 30 + 40 * 10 + 40 * 10)/160
  success += (s.q.compare("5!!") == 0) ? 1 : 0;
  
  ad.push_back(a7);
  s = get_consensus_allele(ad, 0, 0, 'N');
  // A, T, G, AAC, AWG -> DWS
  success += (s.nuc.compare("DWS") == 0) ? 1: 0;
  // (40 * 30 + 40 * 30 + 40 * 10 + 40 * 10 + 40 * 10)/200
  success += (s.q.compare("3!!") == 0) ? 1 : 0;
  std::cout << s.nuc << std::endl;
  std::cout << s.q << std::endl;
  
  ad.push_back(a8);
  s = get_consensus_allele(ad, 0, 0, 'N');
  // A, T, G, AAC, AWG, CCT -> NHB
  success += (s.nuc.compare("NHB") == 0) ? 1: 0;
  // (40 * 30 + 40 * 30 + 40 * 10 + 40 * 10 + 40 * 10 + 40 * 10)/240
  success += (s.q.compare("2!!") == 0) ? 1 : 0;
  std::cout << s.nuc << std::endl;
  std::cout << s.q << std::endl;
  
  allele a9 = {
    "A",
    50,
    5,
    20,
    0,
    0,
    0,
    ""
  };
  ad.push_back(a9);
  s = get_consensus_allele(ad, 0, 0, 'N');
  // A
  success += (s.nuc.compare("A") == 0) ? 1: 0;
  // 20
  success += (s.q.compare("5") == 0) ? 1 : 0;
  std::cout << s.nuc << std::endl;
  return (success == num_tests) ? 0 : -1;
}
