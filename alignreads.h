// global variables that need to be accessible by the file format functions
char *read; // chars of the current read
char *readName; // label of the current read
FILE *readsFile; // file with the reads to be aligned
char *readsFilename; // filename with the reads to be aligned
int maxReadArraySize; // maximum size of the array that stores the read chars (is always larger or equal than maxReadSize)
int maxReadNameArraySize; // size of the longest read name in the reads file (always set to 255)
void IncreaseMaxReadArraySize(int newsize);

// global variables that need to be accessible by the file format functions to load paired end reads
int numPairs; // number of reads currently being processed at the same time (1 or 2)
int currentPair; // index of the current reads file in paired-end mode (0 or 1)
char *fwdReads[2]; // forward strand for both reads in pair
char *revReads[2]; // reverse strand for both reads in pair
char *readNames[2]; // read names for both reads in pair
int readsSizes[2]; // sizes of each one of the paired-end reads
FILE *pairFiles[2]; // both reads files in paired-end mode
char *pairFilenames[2]; // both reads filenames in paired-end mode

// global variables that need to be accessible by the dynamic programming functions
FILE *debugfile; // file to write debug results
