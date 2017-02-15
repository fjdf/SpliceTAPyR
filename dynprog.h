int **dpMatrix; // matrix for dynamic programming alignment scores
char **dpDirections; // matrix for dynamic programming alignment directions
int dpNumRows; // number of rows of the dynamic programming matrices
int dpNumCols; // number of columns of the dynamic programming matrices

// for debugging purposes
char *dpTrace; // dp directions for the current alignment
int dpTraceSize; // size of the alignment string
char *dpAlignedTarget; // alignment string of the text
char *dpAlignedRead; // alignment string of the pattern

int RunBandedDynamicProgramming(char *text, int n, char *pattern, int m, int maxerrors, int seedlocation, int *numgapsatend);
void GetRealErrorStatistics(char *reffilename, char *readsfilename, char *samfilename);
