char *cigarStrCodes; // CIGAR string operation chars
int *cigarStrCounts; // CIGAR string operation counts
int cigarStrSize; // size of the CIGAR string

void ReportSAMStatistics(char *reffilename, char *samfilename, unsigned int realnumreads, int savebadreads);
void EvaluateSAMPositions(char *samfilename, int savewrongreads);
void ConvertSAMToCSV(char *samfilename);
