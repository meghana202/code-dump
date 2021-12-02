#ifndef ALIGN_H
#define ALIGN_H
struct SWAlign{
	int beg1;
	int end1;
	int beg2;
	int end2;
	int score;
	char *seq1;
	char *seq2;
};
typedef struct SWAlign swalign; 
char *rev(const char *);
swalign align(dpmatrix, int, int, int);
void display_swalign(swalign*);
void free_swalign(swalign*);
#endif