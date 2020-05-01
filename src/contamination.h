#include <iostream>
#include <string>
#include <algorithm>

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/bgzf.h"

#include "variants_by_amplicon.h"

void identify_amp(std::vector<primer> &primers, uint64_t pos, uint64_t end_pos, std::vector<primer*> &fwd_primers, std::vector<primer*> &rev_primers);
int identify_contamination(std::string bed, std::string bam, std::string pairs_file, std::string prefix, uint8_t min_qual);
