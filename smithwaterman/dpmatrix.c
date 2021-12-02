#include <stdio.h>
#include "genomikon.h"
#include "dpmatrix.h"

dpmatrix new_dpmatrix(const gkn_fasta ffa, const gkn_fasta ffb) {
	dpmatrix mat = malloc(sizeof(struct DPMatrix));
	mat->seq1 = NULL;
	mat->seq2 = NULL;
	mat->len1 = ffa->length + 1;
	mat->len2 = ffb->length + 1;
	
	mat->seq1 = malloc(strlen(ffa->seq)+1);
	strcpy(mat->seq1,ffa->seq);
	mat->seq2 = malloc(strlen(ffb->seq)+1);
	strcpy(mat->seq2,ffb->seq);
	
	mat->score = malloc(sizeof(int*)*mat->len1);
	mat->trace = malloc(sizeof(char*)*mat->len1);
	for (int i = 0; i < mat->len1; i++) {
		mat->score[i] = malloc(sizeof(int)*mat->len2);
		mat->trace[i] = malloc(sizeof(char*)*mat->len2);
	}
	return mat;
}

dpmatrix new_blank(void) {
	dpmatrix mat = malloc(sizeof(struct DPMatrix));
	mat->len1 = 0;
	mat->len2 = 0;
	mat->seq1 = NULL;
	mat->seq2 = NULL;
	mat->score = NULL;
	mat->trace = NULL;
	return mat;
}

void free_dpmatrix(dpmatrix matrix) {
	if (matrix->seq1) free(matrix->seq1);
	if (matrix->seq2) free(matrix->seq2);
	
	if (matrix->score) {
		for(int i = 0; i < matrix->len1; i++) {
			free(matrix->score[i]);
			free(matrix->trace[i]);
		}
		free(matrix->score);
		free(matrix->trace);
	}
}
