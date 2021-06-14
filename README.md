# bioinf-julia
This code for the Bioinformatics Coursera Specialization written in Julia
1) Checking for the OriC using functions Skew & FrequentMismatch 
  - Skew finds the location in the genome with the least discrepancy between G and C nucleotides, which for many bacterial genomes is the location of the origin of replication.
  - FrequentMismatch finds all of the kmers with maximum reverse complements and hamming distances <= d, all of which may be possible DnaA boxes (they don't necessarily occur in the genome)
  - For Salmonella enterica and E. coli, use the Skew function to find the location of least skew in the genome, then plug in the sequence with windows of length 500bp to the FrequentMismatch function to find possible DnaA boxes with high frequency.
