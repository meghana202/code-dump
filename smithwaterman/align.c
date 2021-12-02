#include <string.h>
#include <stdio.h>
#include "genomikon.h"
#include "dpmatrix.h"
#include "align.h"


char *rev(const char *seq){
	int len = strlen(seq);
	char *reverse = (char*) malloc((len)*sizeof(char));
	for(int i = 0; i <= len-1; i++){
		reverse[i] = seq[len-i-1];
	}
	return reverse;
}

void display_swalign(swalign *output){
	printf("\n");
	printf("Score: %d\n", output->score);
	char *alignment = malloc(strlen(output->seq1));
	for (int i = 0; i < strlen(output->seq1);i++){
		if (output->seq1[i] == output->seq2[i]){
			alignment[i] = '|';
		} else {
			alignment[i] = ' ';
		}
	}
	printf("%d", output->beg1);
	for (int i = 0; i < (strlen(output->seq1)-2); i++){
		printf(" ");
	}
	printf("%d\n", output->end1);
	printf("%s\n", output->seq1);
	printf("%s\n", alignment);
	printf("%s\n", output->seq2);
	printf("%d", output->beg2);
	for (int i = 0; i < (strlen(output->seq2)-2); i++){
		printf(" ");
	}
	printf("%d\n", output->end2);
}

void free_swalign(swalign *output){
	if (output->seq1) free(output->seq1);
	if (output->seq2) free(output->seq2);
}


swalign align(dpmatrix mat, int M, int N, int G) {
	int x = -1;
	int y = -1;
	swalign *output = malloc(sizeof(struct SWAlign));
	
	//Initialization
	for (int i = 0; i < mat->len1; i++) {
		for (int j = 0; j < mat->len2; j++) {
			mat->score[i][j] = 0;
			mat->trace[i][j] = 'N';
		}
	}

	//Fill
	int mm, lgap, tgap;
	int max_score = 0;
	for (int i = 1; i < mat->len1; i++) {
		for (int j = 1; j < mat->len2; j++) {
			//Match vs Mismatch
			char char1 = mat->seq1[i-1];
			char char2 = mat->seq2[j-1];
			
			if (char1 == char2) {
				mm = mat->score[i-1][j-1] + M; 
			}	
			else {
				mm = mat->score[i-1][j-1] + N;
			}
				
			lgap = mat->score[i][j-1] + G;
			tgap = mat->score[i-1][j] + G;
			
			if (lgap > tgap && lgap > mm) {
				mat->score[i][j] = lgap;
				mat->trace[i][j] = 'L';
				if (lgap>max_score) max_score = lgap;
			} else if (tgap > mm) {
				mat->score[i][j] = tgap;
				mat->trace[i][j] = 'U';
				if (tgap>max_score) max_score = tgap;
			} else {
				mat->score[i][j] = mm;
				mat->trace[i][j] = 'D';
				if (mm>max_score) max_score = mm;
			}
				
		}
	}

	//Trace
	for (int i = 1; i < mat->len1; i++) {
		for (int j = 1; j < mat->len2; j++) {
			x = i;
			y = j;
			if (mat->score[i][j] == max_score) {
				char *seq1_align = malloc(mat->len1 + mat->len2);
				char *seq2_align = malloc(mat->len1 + mat->len2);
				output->seq1 = malloc(mat->len1 + mat->len2);
				output->seq2 = malloc(mat->len1 + mat->len2);
				output->score = max_score;
				output->end1 = x;
				output->end2 = y;
				int count = 0;
				while (mat->score[x][y]!=0) {
					if (mat->trace[x][y]=='L') {
						seq2_align[count] = mat->seq2[y-1];
						seq1_align[count] = '-';
						count++;
						y--;
					}
					else if (mat->trace[x][y]=='U') {
						seq1_align[count] = mat->seq1[x-1];
						seq2_align[count] = '-';
						count++;
						x--;
					}
					else {
						seq1_align[count] = mat->seq1[x-1];
						seq2_align[count] = mat->seq2[y-1];
						count++;
						x--;
						y--;
					}	
				} 
				if (mat->score[x][y]==0){
					output->beg1 = x;
					output->beg2 = y;
					if (mat->trace[x][y]=='L') {
						seq2_align[count] = mat->seq2[y-1];
						seq1_align[count] = '-';
					}
					else if (mat->trace[x][y]=='U') {
						seq1_align[count] = mat->seq1[x-1];
						seq2_align[count] = '-';
					}
					else {
						seq1_align[count] = mat->seq1[x-1];
						seq2_align[count] = mat->seq2[y-1];
					}	
				}
				strcpy(output->seq1, rev(seq1_align));
				strcpy(output->seq2, rev(seq2_align));
				free(seq1_align);
				free(seq2_align);
				break;
			}
		}
	}
	return *output;
}

