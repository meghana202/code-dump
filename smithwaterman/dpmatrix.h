#ifndef DP_MATRIX_H
#define DP_MATRIX_H

struct DPMatrix {
	int    len1;
	int    len2;
	char  *seq1;
	char  *seq2;
	int  **score;
	char **trace;
};

typedef struct DPMatrix *dpmatrix;

dpmatrix new_dpmatrix(const gkn_fasta, const gkn_fasta);
dpmatrix new_blank(void);
void free_dpmatrix(dpmatrix);

#endif