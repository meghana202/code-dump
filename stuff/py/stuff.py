import numpy as np
from Bio import SeqIO
import math

def kmers(n):
	letters = ["A", "C", "G", "T"]
	if n == 1:
		answer = letters
		return answer
	if n > 1:
		new = []
		last = kmers(n-1)
		for i in range(len(last)):
			for j in range(len(letters)):
				new.append(f"{last[i]}{letters[j]}")
		return new
		
			
#nth order, x hidden states
def train(n, x):
	all_kmers = kmers(n)
	freq = {}
	for i in range(len(all_kmers)):
		freq[all_kmers[i]] = {"A": [0]*x, "C": [0]*x, "G": [0]*x, "T": [0]*x}
		
	exons = SeqIO.parse(open("exons.fa"),'fasta')
	for fasta in exons:
  		name, seq = fasta.id, str(fasta.seq)
  		for i in range(len(seq)-n):
  			letters = seq[i:i+n]
  			bp = seq[i+n]
  			freq[letters][bp][0] += 1
	
	introns = SeqIO.parse(open("introns.fa"),'fasta')
	for fasta in introns:
  		name, seq = fasta.id, str(fasta.seq)
  		for i in range(len(seq)-n):
  			letters = seq[i:i+n]
  			bp = seq[i+n]
  			freq[letters][bp][1] += 1
  			
	for key in freq:
		exons = 0
		introns = 0
		values_ = list(freq[key].values())
		for i in range(len(values_)):
			exons += values_[i][0]
			introns += values_[i][1]
		for value in freq[key]:
		#if the kmer doesn't show up give it a value of .001 ..?
			if freq[key][value][0] == 0:
				freq[key][value][0] = .01
			if freq[key][value][1] == 0:
				freq[key][value][1] = .01
			if introns == 0:
				introns = 1
			if exons == 0:
				exons = 1
				
			freq[key][value][0] = math.log(freq[key][value][0]/exons, 10)
			freq[key][value][1] = math.log(freq[key][value][1]/introns, 10)

	return freq