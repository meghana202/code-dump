#include <stdio.h>
#include <string.h>
#include <stdlib.h>

// match = 5, mismatch and indel = -1

void rev(const char *seq);


void sw(const char *seq1, const char *seq2) {
	int l1 = strlen(seq1);
	int l2 = strlen(seq2);
	int matrix[l2+1][l1+1];
	char trace[l2+1][l1+1];
	
	//initialization
	trace[0][0] = '-';
	for(int i = 1; i <= l1; i++){
		trace[0][i] = 'L';
	}
	for(int i = 1; i <= l2; i++){
		trace[i][0] = 'U';
	}
	for(int i = 1; i <= l2; i++){
		for(int j = 1; j <= l1; j++){
			trace[i][j] = 'O';
		}
	}
	for(int i = 0; i <= l2; i++){
		for(int j = 0; j <= l1; j++){
			matrix[i][j] = 0;
		}
	}
	
	//fill
	int match; int top; int left; int di;
	int max_score = 0; int max_i = 0; int max_j = 0;
	for(int i = 1; i <= l2; i++){
		for(int j = 1; j <= l1; j++){
			if (seq1[j-1] == seq2[i-1]){
				match = 5;
			} else {
				match = -1;
			}
			
			top = matrix[i-1][j] - 1;
			left = matrix[i][j-1] -1;
			di = matrix[i-1][j-1] + match;
			
			if (0 > di && 0 > top && 0 > left){
				matrix[i][j] = 0;
			} else if (di > top && di > left){
				if (di>max_score){
					max_score = di;
					max_i = i;
					max_j = j;
				}
				matrix[i][j] = di;
				trace[i][j] = 'D';
			} else if (top > left){
				matrix[i][j] = top;
				trace[i][j] = 'U';
			} else {
				matrix[i][j] = left;
				trace[i][j] = 'L';
			}
		}
	}
	
	//trace
	int max;
	if (l2>l1){
		max = l1;
	} else{
		max = l2;
	}
	char firstseq[max];
	char secondseq[max];
	int cell[2] = {max_i, max_j};
	int x = cell[0];
	int y = cell[1];
	for (int i = 0; i <= max; i++){
		if (matrix[x][y] == 0){
			break;
		}
		if (trace[x][y]=='U'){
			secondseq[i] = seq2[x-1];
			firstseq[i] = '-';
			x = x - 1;
		} else if (trace[x][y]=='L'){
			firstseq[i] = seq1[y-1];
			secondseq[i] = '-';
			y = y - 1;
		} else if (trace[x][y]=='D'){
			secondseq[i] = seq2[x-1];
			firstseq[i] = seq1[y-1];
			x = x - 1;
			y = y - 1;
		}
	}
	
	//output format
	firstseq[l1+1] = '\0';
	secondseq[l2+1] = '\0';
	rev(firstseq);
	rev(secondseq);
/*
	for(int i = 0; i <= l2; i++){
		for(int j = 0; j <= l1; j++){
			printf("%c ", trace[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	for(int i = 0; i <= l2; i++){
		for(int j = 0; j <= l1; j++){
			if (matrix[i][j] < 10){
				printf("%d  ", matrix[i][j]);
			}
			else {
				printf("%d ", matrix[i][j]);
			}
		}
		printf("\n");
	}
*/

}

void rev(const char *seq){
	int len = strlen(seq);
	char reverse[len];
	for(int i = 0; i <= len-1; i++){
		reverse[i] = seq[len-i-1];
	}
	printf("%s\n", reverse);
}


int main (int argc, char **argv) {
	char *seq1 = argv[1];
	char *seq2 = argv[2];
	sw(seq1, seq2);
}
