#include <iostream>
#include "../src/vcf_writer.h"
#include <fstream>

void write_test_vcf(std::vector<allele> a, char ref){
  uint32_t pos = 30;
  std::string region = "test_vcf";
  vcf_writer *vw = new vcf_writer(0, "../data/test_vcf.vcf", region, "test_vcf");
  vw->write_record(pos, a, region, ref);
  delete vw;
}

int check_vcf(std::string obs){
  int comp;
  std::string l;
  std::ifstream in_file("../data/test_vcf.vcf");
  while(l.empty() || l[0] == '#')
    getline(in_file, l);
  comp = l.compare(obs) == 0 ? 0 : 1;
  in_file.close();
  return comp;
}

int main(int argc, char *argv[])
{
  int success = 0;
  char ref = 'A';
  std::vector<allele> ad = {{
      "ATG",
      1,
      2,
      3
    },{
      "T",
      4,
      5,
      6
    }, {
      "G",
      7,
      8,
      9,
      0,
      0,
      0,
      "TGC"
    }};
  write_test_vcf(ad, ref);
  write_test_vcf(ad, 'T');
  success += check_vcf("test_vcf\t30\t.\tTTGC\tATG,G\t6\tPASS\tDP=1,7;AD=1,7;ADF=-1,-1;ADR=2,8;GQ=3,9;AF=0.0833333,0.583333\tGT:AD:GQ:ADF:ADR:AF\t1|0|2:1,4,7:3,6,9:-1,-1,-1:2,5,8:0.0833333,0.333333,0.583333");

  ad.clear();
  ad = {{
      "A",
      1,
      2,
      3
    },{
      "T",
      4,
      5,
      6
    }, {
      "G",
      7,
      8,
      9,
      0,
      0,
      0,
      "TGC"
    }};
  write_test_vcf(ad, ref);
  write_test_vcf(ad, 'A');
  success += check_vcf("test_vcf\t30\t.\tATGC\tT,G\t3\tPASS\tDP=4,7;AD=4,7;ADF=-1,-1;ADR=5,8;GQ=6,9;AF=0.333333,0.583333\tGT:AD:GQ:ADF:ADR:AF\t0|1|2:1,4,7:3,6,9:-1,-1,-1:2,5,8:0.0833333,0.333333,0.583333");


  ad.clear();
  ad = {{
      "A",
      1,
      2,
      3,
      0,
      0,
      0,
      ""
    },{
      "T",
      4,
      5,
      6,
      0,
      0,
      0,
      ""
    }, {
      "G",
      7,
      8,
      9,
      0,
      0,
      0,
      ""
    }, {
      "GTA",
      7,
      8,
      9,
      0,
      0,
      0,
      ""
    }};
  write_test_vcf(ad, ref);
  write_test_vcf(ad, 'G');
  success += check_vcf("test_vcf\t30\t.\tG\tA,T,GTA\t9\tPASS\tDP=1,4,7;AD=1,4,7;ADF=-1,-1,-1;ADR=2,5,8;GQ=3,6,9;AF=0.0526316,0.210526,0.368421\tGT:AD:GQ:ADF:ADR:AF\t1|2|0|3:1,4,7,7:3,6,9,9:-1,-1,-1,-1:2,5,8,8:0.0526316,0.210526,0.368421,0.368421");
  return success;
}

