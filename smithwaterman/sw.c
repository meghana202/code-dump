#include <stdio.h>
#include <string.h>

#include "genomikon.h"
#include "dpmatrix.h"
#include "align.h"

static char *usage = "\
sw - align sequences\n\n\
usage: sw <file1> <file2> [options]\n\
options:\n\
  -m <int>	match score [1]\n\
  -n <int>	mismatch score [-1]\n\
  -g <int>	gap score [-1]\n\
  -v		verbose\n";
  

int main(int argc, char **argv) {
	int M = 1;
	int N = -1;
	int G = -1;
	int verbose = 0;
	// Command Line Interface
	gkn_set_program_name(argv[0]);
	gkn_register_option("-m", 1);
	gkn_register_option("-n", 1);
	gkn_register_option("-g", 1);
	gkn_register_option("-v", 0);
	gkn_parse_options(&argc, argv);

	if (argc == 1) gkn_exit("%s", usage);

	char *file1 = argv[1];
	char *file2 = argv[2];
	if (gkn_option("-m")) M = atoi(gkn_option("-m"));
	if (gkn_option("-n")) N = atoi(gkn_option("-n"));
	if (gkn_option("-g")) G = atoi(gkn_option("-g"));
	if (gkn_option("-v")) verbose = 1;
	
	gkn_pipe fha = gkn_pipe_open(file1, "r"); 
	gkn_pipe fhb = gkn_pipe_open(file2, "r"); 
	
	gkn_fasta ffa = gkn_fasta_read(fha); 
	gkn_fasta ffb = gkn_fasta_read(fhb);
	
	printf("%s\n%s\n%s\n%s\n", ffa->def, ffa->seq, ffb->def, ffb->seq);
	dpmatrix mat = new_dpmatrix(ffa, ffb);
	align(mat, M, N, G);
	
	swalign ret = align(mat, M, N, G); 
	display_swalign(&ret);

	free_swalign(&ret);
	
	free_dpmatrix(mat);
}