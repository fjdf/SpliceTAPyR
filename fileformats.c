#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "alignreads.h"
#include "tools.h"

#ifdef _MSC_VER
	#pragma warning(disable:4996)
#endif

#define FLX_LINKER "GTTGGAACCGAAAGGGTTTGAATTCAAACCCTTTCGGTTCCAAC"
#define TIT_LINKER1 "TCGTATAACTTCGTATAATGTATGCTATACGAAGTTATTACG"
#define TIT_LINKER2 "CGTAATAACTTCGTATAGCATACATTATACGAAGTTATACGA"
#define ION_LINKER1 "CTGCTGTACGGCCAAGGCGGATGTACGGTACAGCAG"
#define ION_LINKER2 "CTGCTGTACCGTACATCCGCCTTGGCCGTACAGCAG"

uint32_t sff_number_of_reads = 0; // global variable for total number of reads inside current sff file
uint16_t sff_number_of_flows = 0; // global variable for number of flows inside current sff file

int sffReadsCount; // SFF parsing function needs to know if it has processed all the reads in the current file
short int numLinkers;
short int linkerSize;
char linkerSequence[2][64];
unsigned long long int linkerCharMasks[2][5]; // one mask for each char: A C G T N
unsigned long long int linkerLastPosMask; // bit array for exact matching detection
unsigned long long int linkerSearchBitArray; // working bit array for exact matching
unsigned long long linkerSearchBitArrays[64]; // working bit arrays for inexact matching
int linkerStartPosInRead; // where the linker starts in the read
int linkerEndPosInRead; // where the linker ends in the read
int minLinkerMatchSize; // minimum size of the continuous exact match of the linker
int maxLinkerErrors; // maximum number of errors allowed in the linker
int numLinkerErrors; // number of errors found in the linker
int allowUnpairedReads; // if we should allow paired-end reads were one mate is missing 


// Counts the total number of reads and the size of the longest read in the reads file in the FASTA format
char AnalyzeFastaReads(int preprocess, int verbose){
	char c;
	int readSize, avgReadSize, readNameSize;
	int totalNumReads, maxReadSize, maxReadNameSize;
	/*
	FILE *readsfile;
	readsfile = fopen(readsfilename,"r");
	if( readsfile == NULL ) {
		printf("\n> ERROR: Reads file not found\n");
		exit(0);
	}
	*/
	c = fgetc(readsFile);
	if( c != '>' ) {
		printf("> ERROR: Invalid FASTA file\n");
		exit(0);
	}
	if(verbose) printf("(FASTA) ");
	fflush(stdout);
	totalNumReads = 0;
	maxReadSize = 0;
	maxReadNameSize = 0;
	avgReadSize = 0;
	while( c != EOF ) {
		if( c != '>' ) break; // start of a new read
		readNameSize = 0;
		c = fgetc(readsFile);
		while( (c != '\n') && (c != EOF) ){ // read description
			readNameSize++;
			c = fgetc(readsFile);
		}
		if( readNameSize > maxReadNameSize ) maxReadNameSize = readNameSize;
		readSize = 0;
		c = fgetc(readsFile);
		while( (c != '>') && (c != EOF) ){ // read characters
			if( c=='A' || c=='C' || c=='G' || c=='T' || c=='a' || c=='c' || c=='g' || c=='t' ) readSize++; // only count valid chars
			c = fgetc(readsFile);
		}
		if( readSize != 0 ) {
			totalNumReads++;
			if( readSize > maxReadSize ) maxReadSize = readSize;
			avgReadSize += readSize;
		}
		if(!preprocess) break; // if we do not want to preprocess the entire reads file, break after the first read
	}
	if( totalNumReads == 0 ) {
		printf("\n> ERROR: No valid data in reads file\n");
		exit(0);
	}
	rewind(readsFile); // back to beginning of file
	if(!preprocess) return 'F';
	if(verbose){
		avgReadSize = ( avgReadSize / totalNumReads );
		PrintNumber((long long int)totalNumReads);
		printf(" reads (average size = %dbp) ",avgReadSize);
		fflush(stdout);
	}
	return 'F';
}

// Loads the next read in the reads file in the FASTA format
int LoadNextReadFromFasta(){
	int n;
	char c;
	c = fgetc(readsFile); // get '>' char
	if( (c != '>') || (c == EOF) ) return 0;
	n = 0;
	c = fgetc(readsFile);
	while( (c != '\n') && (c != EOF) ){ // get read name
		if( n < maxReadNameArraySize ) readName[n++] = c;
		c = fgetc(readsFile);
	}
	readName[n] = '\0';
	n = 0;
	c = fgetc(readsFile);
	while( (c != '>') && (c != EOF) ) { // get valid read chars
		if( n == maxReadArraySize ) IncreaseMaxReadArraySize(n+1); // increase read array size if needed
		if( c=='A' || c=='C' || c=='G' || c=='T' ) read[n++] = c;
		else if( c=='a' || c=='c' || c=='g' || c=='t' )	read[n++] = (char)(c-32);
		c = fgetc(readsFile);
	}
	read[n] = '\0';
	ungetc(c,readsFile); // put back '>' char of next read
	if( n != 0 ) return n; // valid read is present
	if( c == EOF ) return 0; // end of file
	return LoadNextReadFromFasta(); // invalid read, load next one
}

// Loads both mates of a paired-end reads pair from two distinct FASTA files
int LoadNextReadPairFromFasta(){
	currentPair=0; // needed if IncreaseMaxReadArraySize is called
	readsFile=pairFiles[0];
	read=fwdReads[0];
	readName=readNames[0];
	readsSizes[0]=LoadNextReadFromFasta();
	currentPair=1;
	readsFile=pairFiles[1];
	read=fwdReads[1];
	readName=readNames[1];
	readsSizes[1]=LoadNextReadFromFasta();
	currentPair=0; // back to initial pair id
	return 1;
}

// Counts the total number of reads and the size of the longest read in the reads file in the FASTAQ format
char AnalyzeFastaQReads(int preprocess, int verbose){
	char c;
	int readSize, avgReadSize, readNameSize;
	int totalNumReads, maxReadSize, maxReadNameSize;
	/*
	FILE *readsfile;
	readsfile = fopen(readsfilename,"r");
	if( readsfile == NULL ) {
		printf("\n> ERROR: Reads file not found\n");
		exit(0);
	}
	*/
	c = fgetc(readsFile);
	if( c != '@' ) {
		printf("> ERROR: Invalid FASTQ file\n");
		exit(0);
	}
	if(verbose) printf("(FASTQ) ");
	fflush(stdout);
	totalNumReads = 0;
	maxReadSize = 0;
	maxReadNameSize = 0;
	avgReadSize = 0;
	while( c != EOF ) {
		if( c != '@' ) break; // start of a new read
		readNameSize = 0;
		c = fgetc(readsFile);
		while( (c != '\n') && (c != EOF) ){ // read description
			readNameSize++;
			c = fgetc(readsFile);
		}
		if( readNameSize > maxReadNameSize ) maxReadNameSize = readNameSize;
		readSize = 0;
		c = fgetc(readsFile);
		while( (c != '\n') && (c != EOF) ){ // read characters
			if( c=='A' || c=='C' || c=='G' || c=='T' || c=='a' || c=='c' || c=='g' || c=='t' ) readSize++; // only count valid chars
			c = fgetc(readsFile);
		}
		c = fgetc(readsFile);
		if( c != '+' ) break;
		while( (c != '\n') && (c != EOF) ) c = fgetc(readsFile); // skip line with read information
		c = fgetc(readsFile);
		while( (c != '\n') && (c != EOF) ) c = fgetc(readsFile); // skip line with phred qualities
		if( readSize != 0 ){
			totalNumReads++;
			if( readSize > maxReadSize ) maxReadSize = readSize;
			avgReadSize += readSize;
		}
		c = fgetc(readsFile); // get first char of next line
		if(!preprocess) break; // if we do not want to preprocess the entire reads file, break after the first read
	}
	if( totalNumReads == 0 ) {
		printf("\n> ERROR: No valid data in reads file\n");
		exit(0);
	}
	rewind(readsFile); // back to beginning of file
	if(!preprocess) return 'Q';
	if(verbose){
		avgReadSize = ( avgReadSize / totalNumReads );
		PrintNumber((long long int)totalNumReads);
		printf(" reads (average size = %dbp) ",avgReadSize);
		fflush(stdout);
	}
	return 'Q';
}

// Loads the next read in the reads file in the FASTQ format
int LoadNextReadFromFastaQ(){
	int n;
	char c;
	c = fgetc(readsFile); // get '@' char
	if( c == EOF ) return 0;
	while( (c != '@') && (c != EOF) ) c = fgetc(readsFile);
	n = 0;
	c = fgetc(readsFile);
	while( (c != '\n') && (c != EOF) ){ // get read name
		if( n < maxReadNameArraySize ) readName[n++] = c;
		c = fgetc(readsFile);
	}
	readName[n] = '\0';
	n = 0;
	c = fgetc(readsFile);
	while( (c != '+') && (c != EOF) ) { // get valid read chars (until description line '+')
		if( n == maxReadArraySize ) IncreaseMaxReadArraySize(n+1); // increase read array size if needed
		if( c=='A' || c=='C' || c=='G' || c=='T' ) read[n++] = c;
		else if( c=='a' || c=='c' || c=='g' || c=='t' ) read[n++] = (char)(c-32);
		c = fgetc(readsFile);
	}
	read[n] = '\0';
	//c = fgetc(readsFile);
	//if( (c != '+') || (c == EOF) ) return 0;
	while( (c != '\n') && (c != EOF) ) c = fgetc(readsFile); // skip line with read information
	c = fgetc(readsFile);
	while( (c != '\n') && (c != EOF) ) c = fgetc(readsFile); // skip line with phred qualities
	if( n != 0 ) return n; // valid read is present
	if( c == EOF ) return 0; // end of file
	return LoadNextReadFromFastaQ(); // invalid read, load next one
}

// Loads both mates of a paired-end reads pair from two distinct FASTAQ files
int LoadNextReadPairFromFastaQ(){
	currentPair=0; // needed if IncreaseMaxReadArraySize is called
	readsFile=pairFiles[0];
	read=fwdReads[0];
	readName=readNames[0];
	readsSizes[0]=LoadNextReadFromFastaQ();
	currentPair=1;
	readsFile=pairFiles[1];
	read=fwdReads[1];
	readName=readNames[1];
	readsSizes[1]=LoadNextReadFromFastaQ();
	currentPair=0; // back to initial pair id
	return 1;
}

void BuildSffLinkerMasks(){
	int i, j, k;
	char *linker, c;
	unsigned long long int *masks;
	for(j=0;j<numLinkers;j++){
		linker=linkerSequence[j];
		masks=linkerCharMasks[j];
		for(k=0;k<4;k++) masks[k] = 0ULL; // reset char masks
		linkerSize=(short int)strlen(linker);
		for(i=0;i<linkerSize;i++){
			c=linker[i];
			switch(c){ // set bit on corresponding char mask
			case 'A':
				masks[0] |= 1ULL;
				break;
			case 'C':
				masks[1] |= 1ULL;
				break;
			case 'G':
				masks[2] |= 1ULL;
				break;
			case 'T':
				masks[3] |= 1ULL;
				break;
			default:
				break;
			} // all chars of alphabet
			for(k=0;k<4;k++) masks[k] <<= 1; // shift all char masks to the left
		} // all chars of linker
		for(k=0;k<4;k++) masks[k] >>= 1; // undo the last un-needed shift
		masks[4]=0ULL; // set mask for invalid chars
	} // all linkers
	linkerLastPosMask = ( 1ULL << ( linkerSize - 1 ) );
	//linkerSearchBitArrays = (unsigned long long *)malloc((linkerSize+2)*sizeof(unsigned long long)); // bit arrays for inexact matching (no need to set, because initialized at size 64)
	linkerSearchBitArrays[0] = ~0ULL; // the 0-th position is initialized but never used
	minLinkerMatchSize=(int)(linkerSize/5); // linker should have at least 20% of exact match
	maxLinkerErrors=(int)(linkerSize/5); // 20% of errors allowed in the linker
	allowUnpairedReads=0; // only allow reads when both mates are present
}

void PrintBinary(unsigned long long number){
	unsigned long long mask;
	int i;
	mask=(1ULL<<63);
	for(i=0;i<64;i++){
		printf("%c",(number&mask)?'1':'0');
		mask>>=1;
	}
}

// Searches for an exact match of the paired-ends linker in the read
int FindExactSffLinker(int readSize){
	int i, j;
	char c;
	unsigned long long int *masks;
	for(j=0;j<numLinkers;j++){
		/*
		printf("\n");
		printf("L: %64s\n",linkerSequence[j]);
		printf("*: ");PrintBinary(linkerLastPosMask);printf("\n");
		printf("A: ");PrintBinary(linkerCharMasks[j][0]);printf("\n");
		printf("C: ");PrintBinary(linkerCharMasks[j][1]);printf("\n");
		printf("G: ");PrintBinary(linkerCharMasks[j][2]);printf("\n");
		printf("T: ");PrintBinary(linkerCharMasks[j][3]);printf("\n");
		printf("-: ");PrintBinary(linkerCharMasks[j][4]);printf("\n");
		*/
		masks=linkerCharMasks[j];
		linkerSearchBitArray = 0ULL; // reset bit array
		for(i=(readSize-1);i>=0;i--){
			linkerSearchBitArray <<= 1; // shift bit array to the left
			linkerSearchBitArray |= 1ULL; // add a bit in the first/rightmost position
			c=read[i];
			switch(c){ // bit-and the bit array with the corresponding char mask
			case 'A':
				linkerSearchBitArray &= masks[0];
				break;
			case 'C':
				linkerSearchBitArray &= masks[1];
				break;
			case 'G':
				linkerSearchBitArray &= masks[2];
				break;
			case 'T':
				linkerSearchBitArray &= masks[3];
				break;
			default:
				linkerSearchBitArray &= masks[4];
				break;
			} // all chars of alphabet
			if( linkerSearchBitArray & linkerLastPosMask ){ // check if the linker was found
				numLinkerErrors=0;
				linkerStartPosInRead=i;
				linkerEndPosInRead=(i+linkerSize-1);
				return 0; // linker with no errors
			}
			/*
			printf(" %c ",c);PrintBinary(linkerSearchBitArray);printf("\n");
			*/
		} // all chars of read
	} // all linkers
	numLinkerErrors=readSize;
	linkerStartPosInRead=readSize;
	linkerEndPosInRead=(readSize-1);
	return -1; // exact match of the linker was not found
}

// Searches for the best inexact match of the linker inside the read
// NOTE: returns the number of errors of the best match (if it is >= maxLinkerErrors, it returns -1 meaning that the linker was not found)
// NOTE: it sets the variables linkerStartPosInRead and linkerEndPosInRead to later retrieve the linker
int FindInexactSffLinker(int readSize){
	int i, j, k, n;
	char c, *linker;
	unsigned long long *allMasks, charMask, longestMatchBitArray;
	int longestMatchSize, longestMatchLinkerId, longestMatchPosInLinker, longestMatchPosInRead;
	longestMatchLinkerId=0; // id of the linker where the longest match was found
	longestMatchSize=0; // size of the longest match for the entire read
	longestMatchBitArray=0ULL; // bit array corresponding to the longest match found
	longestMatchPosInLinker=-1; // position in the linker where the longest match starts
	longestMatchPosInRead=-1; // position in the read where the longest match starts
	/*
	printf("\n");
	printf(">%s\n",readName);
	for(j=0;j<numLinkers;j++){
		printf("   ");for(i=63;i>=linkerSize;i--) printf(" ");
		for(i=0;i<=(linkerSize-1);i++) printf("%d",(i%10));printf("\n");
		printf("L: %64s\n",linkerSequence[j]);
		printf("A: ");PrintBinary(linkerCharMasks[j][0]);printf("\n");
		printf("C: ");PrintBinary(linkerCharMasks[j][1]);printf("\n");
		printf("G: ");PrintBinary(linkerCharMasks[j][2]);printf("\n");
		printf("T: ");PrintBinary(linkerCharMasks[j][3]);printf("\n");
		printf("-: ");PrintBinary(linkerCharMasks[j][4]);printf("\n");
		printf("\n");
	}
	*/
	for(j=0;j<numLinkers;j++){
		allMasks=linkerCharMasks[j];
		n=0; // size of the current longest match
		for(i=(readSize-1);i>=0;i--){ // the search for the linker is performed in reverse (from right to left)
			c=read[i];
			switch(c){ // get mask for current char
				case 'A':
					charMask=allMasks[0];
					break;
				case 'C':
					charMask=allMasks[1];
					break;
				case 'G':
					charMask=allMasks[2];
					break;
				case 'T':
					charMask=allMasks[3];
					break;
				default:
					charMask=allMasks[4];
					break;
			} // all chars of alphabet
			for( k=(n+1) ; k>1 ; k-- ) linkerSearchBitArrays[k] = ( (linkerSearchBitArrays[k-1] << 1 ) & charMask ); // the bit array at position k stores the matches with length k
			linkerSearchBitArrays[1] = charMask; // match of length 1 is the current char
			if( linkerSearchBitArrays[(n+1)] != 0ULL ){ // check if the current longest match was extended
				n++;
				if(n>longestMatchSize){
					longestMatchSize=n;
					longestMatchBitArray=linkerSearchBitArrays[n];
					longestMatchPosInRead=i;
					longestMatchLinkerId=j;
				}
			} else { // if the previous longest match ended, search back for the new one
				while( linkerSearchBitArrays[n] == 0ULL ) n--;
			}
		} // all chars of read
	} // all linkers
	/*
	printf("   ");PrintBinary(longestMatchBitArray);printf(" %d @ R:%d",longestMatchSize,longestMatchPosInRead);
	*/
	if( longestMatchSize<minLinkerMatchSize && longestMatchPosInRead!=0 && (longestMatchPosInRead+longestMatchSize)!=readSize ){ // linker not found (match is too short, and not at the ends of the read)
		linkerStartPosInRead=readSize; // set these variables like this so the function LoadNextReadPairFromSff can split the read correctly
		linkerEndPosInRead=(readSize-1);
		return -1;
	}
	longestMatchPosInLinker=0; // now get the position in the linker where the longest match ends, counting from right to left (location of the leftmost bit set to 1)
	longestMatchBitArray >>= 1; // so first position is 0, and not 1
	k=32; // positions range from 0 to 63 in the 64 bit array, so 32 is the highest power of 2 in the interval
	while( longestMatchBitArray ){ // remove decreasing powers of 2 until the number reaches 0
		if( ( longestMatchBitArray >> (k-1) ) != 0ULL ){ // when k=1, no shift is performed, so this is not 0
			longestMatchPosInLinker += k; // the number contains this power of 2
			longestMatchBitArray >>= k; // = ( pos - 2^k )
		}
		k >>= 1; // divide by 2 (decrease the power of 2)
	}
	longestMatchPosInLinker=(linkerSize-longestMatchPosInLinker-1); // it is the leftmost position but counted from right to left, and now get the position as if the 0-th position is on the left 
	/*
	printf(" L:%d\n",longestMatchPosInLinker);
	printf("\n");
	for(i=0;i<linkerSize;i++) printf("%d",(i%10));
	printf("\n");
	n=(longestMatchPosInLinker+longestMatchSize);
	for(i=0;i<longestMatchPosInLinker;i++) printf("%c",(char)(linkerSequence[longestMatchLinkerId][i]+32));
	for(i=longestMatchPosInLinker;i<n;i++) printf("%c",linkerSequence[longestMatchLinkerId][i]);
	for(i=n;i<linkerSize;i++) printf("%c",(char)(linkerSequence[longestMatchLinkerId][i]+32));
	printf("\n");
	for(i=0;i<readSize;i++) printf("%d",(i%10));
	printf("\n");
	n=(longestMatchPosInRead+longestMatchSize);
	for(i=0;i<longestMatchPosInRead;i++) printf("%c",(char)(read[i]+32));
	for(i=longestMatchPosInRead;i<n;i++) printf("%c",read[i]);
	for(i=n;i<readSize;i++) printf("%c",(char)(read[i]+32));
	printf("\n");
	fflush(stdout);
	*/
	numLinkerErrors=0; // number of errors found in the linker
	linker=linkerSequence[longestMatchLinkerId];
	i=(longestMatchPosInRead-1); // position in read
	j=(longestMatchPosInLinker-1); // position in linker
	while(i>=0 && j>=0 && numLinkerErrors<maxLinkerErrors){ // extend the found linker segment to the left
		if(read[i]==linker[j]){ // match
			i--;
			j--;
		} else if(i>0 && read[i-1]==linker[j]){ // insertion (added a char to the read)
			i-=2;
			j-=1;
			numLinkerErrors++;
		} else if(j>0 && read[i]==linker[j-1]){ // deletion (missed a char from the linker)
			i-=1;
			j-=2;
			numLinkerErrors++;
		} else { // mismatch
			i-=1;
			j-=1;
			numLinkerErrors++;
		}
	}
	if(i<0) i=0; // if it reached the left end of the read, set it to the first valid position
	else i++; // it was one position to the left
	linkerStartPosInRead=i; // save linker starting position in read for later use
	i=(longestMatchPosInRead+longestMatchSize); // position in read
	j=(longestMatchPosInLinker+longestMatchSize); // position in linker
	while(i<readSize && j<linkerSize && numLinkerErrors<maxLinkerErrors){ // extend the found linker segment to the right
		if(read[i]==linker[j]){ // match
			i++;
			j++;
		} else if(i<readSize && read[i+1]==linker[j]){ // insertion (added a char to the read)
			i+=2;
			j+=1;
			numLinkerErrors++;
		} else if(j<linkerSize && read[i]==linker[j+1]){ // deletion (missed a char from the linker)
			i+=1;
			j+=2;
			numLinkerErrors++;
		} else { // mismatch
			i+=1;
			j+=1;
			numLinkerErrors++;
		}
	}
	if(i>=readSize) i=(readSize-1); // if it reached the right end of the read, set it to the last valid position
	else i--; // it was one position to the right
	linkerEndPosInRead=i; // save linker ending position in read for later use
	/*
	printf("\n");
	printf("numLinkerErrors=%d (max=%d)\n",numLinkerErrors,maxLinkerErrors);
	if(numLinkerErrors<maxLinkerErrors){
		printf("%d [ %d - %d ]\n",(linkerEndPosInRead-linkerStartPosInRead+1),linkerStartPosInRead,linkerEndPosInRead);
		printf("\n");
		for(i=0;i<readSize;i++) printf("%d",(i%10));
		printf("\n");
		for(i=0;i<linkerStartPosInRead;i++) printf("%c",(char)(read[i]+32));
		for(i=linkerStartPosInRead;i<=linkerEndPosInRead;i++) printf("%c",read[i]);
		for(i=(linkerEndPosInRead+1);i<readSize;i++) printf("%c",(char)(read[i]+32));
		printf("\n");
		fflush(stdout);
	}
	getchar();
	*/
	//free(linkerSearchBitArrays);
	if(numLinkerErrors>=maxLinkerErrors){ // linker not found (too many errors)
		linkerStartPosInRead=readSize;
		linkerEndPosInRead=(readSize-1);
		return -1;
	}
	if( linkerStartPosInRead!=0 && linkerStartPosInRead<minLinkerMatchSize ) linkerStartPosInRead=0; // if the left mate is too short, skip it completely
	if( linkerEndPosInRead!=(readSize-1) && (readSize-linkerEndPosInRead-1)<minLinkerMatchSize ) linkerEndPosInRead=(readSize-1); // if the right mate is too short, skip it completely
	return numLinkerErrors;
}

// Loads the next read in the reads file in the SFF format
int LoadNextReadFromSff(){
	char c;
	int readSize;
	unsigned int i, n;
	unsigned int number_of_bases;
	unsigned int padding_size;
	unsigned int start_base_pos;
	unsigned int end_base_pos;
	uint16_t sff_name_length;
	uint32_t sff_number_of_bases;
	uint16_t sff_clip_qual_left;
	uint16_t sff_clip_qual_right;
	uint16_t sff_clip_adapter_left;
	uint16_t sff_clip_adapter_right;
	size_t ncount = 0; ncount = (size_t)ncount;
	if( sffReadsCount == ((int)sff_number_of_reads) ) return 0; // check if we already processed all reads inside file (to prevent reading incorrect data from the index at the end of the file)
	for(i=0;i<2;i++) c=fgetc(readsFile); // skip 2 bytes: read_header_length[2]
	fread16be(&sff_name_length,readsFile); // name_length[2]
	fread32be(&sff_number_of_bases,readsFile); // number_of_bases[4]
	fread16be(&sff_clip_qual_left,readsFile); // clip_qual_left[2]
	fread16be(&sff_clip_qual_right,readsFile); // clip_qual_right[2]
	fread16be(&sff_clip_adapter_left,readsFile); // clip_adapter_left[2]
	fread16be(&sff_clip_adapter_right,readsFile); // clip_adapter_right[2]
	if( sff_clip_qual_left > sff_clip_adapter_left ) start_base_pos = sff_clip_qual_left; // choose maximum value
	else start_base_pos = sff_clip_adapter_left;
	if( start_base_pos > 0 ) start_base_pos--; // convert to 0-based array
	if( sff_clip_qual_right == 0 ) sff_clip_qual_right = sff_number_of_bases; // highest right position is number_of_bases
	if( sff_clip_adapter_right == 0 ) sff_clip_adapter_right = sff_number_of_bases;
	if( sff_clip_qual_right < sff_clip_adapter_right) end_base_pos = sff_clip_qual_right; // choose minimum value
	else end_base_pos = sff_clip_adapter_right;
	end_base_pos--; // convert to 0-based array
	number_of_bases = ( end_base_pos - start_base_pos + 1 ); // real size of the read after clipping ends
	ncount = fread(readName,sizeof(char),(size_t)sff_name_length,readsFile); // name[name_length]
	readName[sff_name_length] = '\0'; // add terminator char to string
	padding_size = ( ( 8 - ( ( 16 + sff_name_length ) & 7 ) ) & 7 );
	n = ( start_base_pos ); // skip the left clipped bases
	n = ( padding_size + 2*sff_number_of_flows + sff_number_of_bases + n ); // skip: eight_byte_padding[padding_size] + flowgram_values[2][number_of_flows] + flow_index_per_base[number_of_bases] + (left_clip_size)
	for(i=0;i<n;i++) c=fgetc(readsFile); // end of read header section and beginning of read data section
	if( (int)number_of_bases > maxReadArraySize ) IncreaseMaxReadArraySize((int)number_of_bases); // increase read array size if needed
	n = 0;
	for(i=0;i<number_of_bases;i++){ // get valid read chars
		c=fgetc(readsFile);
		if( c=='A' || c=='C' || c=='G' || c=='T' ) read[n++] = c;
		else if( c=='a' || c=='c' || c=='g' || c=='t' ) read[n++] = (char)(c-32);
	}
	read[n] = '\0';
	readSize = ((int)n); // store read size
	padding_size = ( ( 8 - ( ( 2*sff_number_of_flows + 3*sff_number_of_bases ) & 7 ) ) & 7 );
	n = ( sff_number_of_bases - end_base_pos - 1 ); // skip the right clipped bases
	n = ( padding_size + sff_number_of_bases + n ); // skip: quality_scores[number_of_bases] + eight_byte_padding[padding_size] + (right_clip_size)
	for(i=0;i<n;i++) c=fgetc(readsFile);
	sffReadsCount++; // update reads count for this SFF file
	if( c == EOF ) return 0; // if we have reached the end of file and the last read is incomplete, do not output it
	if( readSize == 0 ) return LoadNextReadFromSff(); // the read has no valid characters, so, load the next one
	return readSize; // the read is valid
}

// Checks if the SFF file contains paired-end reads joined by a linker or not
int CheckForPairedEndsSff(){
	char *savedFwdRead, *savedRevRead, *savedReadName;
	int savedMaxReadSize, savedMaxReadNameSize, k, n;
	savedMaxReadSize=maxReadArraySize; // save old values
	savedMaxReadNameSize=maxReadNameArraySize;
	savedFwdRead=fwdReads[0];
	savedRevRead=revReads[0];
	savedReadName=readName;
	maxReadArraySize=255; // set temporary values that will be used by LoadNextReadFromSff function so it will not make any changes to the real arrays
	maxReadNameArraySize=255;
	fwdReads[0]=(char *)malloc((maxReadArraySize+1)*sizeof(char));
	revReads[0]=(char *)malloc((maxReadArraySize+1)*sizeof(char));
	readName=(char *)malloc((maxReadNameArraySize+1)*sizeof(char));
	sffReadsCount=0; // initialize counter that will be updated inside LoadNextReadFromSff
	numPairs=1; // set some variables that are used if IncreaseMaxReadArraySize is called inside LoadNextReadFromSff
	currentPair=0;
	read=fwdReads[0];
	while( sffReadsCount<100 && (k=LoadNextReadFromSff())!=0 ){ // search the first 100 reads for an exact match of the linker
		if( (n=FindExactSffLinker(k)) >= 0 ) break;
	}
	free(fwdReads[0]);
	free(revReads[0]);
	free(readName);
	maxReadArraySize=savedMaxReadSize; // restore original values
	maxReadNameArraySize=savedMaxReadNameSize;
	fwdReads[0]=savedFwdRead;
	revReads[0]=savedRevRead;
	readName=savedReadName;
	if(n>=0) return 1;
	return 0;
}

// Counts the total number of reads and the size of the longest read in the reads file in the SFF format
char AnalyzeSffReads(int preprocess, int verbose){
	fpos_t readssectionpos;
	char c, sffType;
	int readSize, avgReadSize, readNameSize;
	int totalNumReads, maxReadSize, maxReadNameSize;
	unsigned int i, n;
	unsigned int number_of_reads;
	unsigned int number_of_bases;
	unsigned int padding_size;
	unsigned int start_base_pos;
	unsigned int end_base_pos;
	uint32_t sff_magic_number;
	//uint32_t sff_number_of_reads; // moved to global variable
	uint16_t sff_key_length;
	//uint16_t sff_number_of_flows; // moved to global variable
	uint16_t sff_name_length;
	uint32_t sff_number_of_bases;
	uint16_t sff_clip_qual_left;
	uint16_t sff_clip_qual_right;
	uint16_t sff_clip_adapter_left;
	uint16_t sff_clip_adapter_right;
	readsFile = freopen(readsFilename,"rb",readsFile); // re-open file in binary mode
	if( readsFile == NULL ) {
		printf("\n> ERROR: Cannot re-open reads file\n");
		exit(0);
	}
	fread32be(&sff_magic_number,readsFile); // magic_number[4]
	if( sff_magic_number != 0x2E736666 ){ // ".sff"
		printf("\n> ERROR: Invalid SFF file\n");
		exit(0);
	}
	if(verbose) printf("(SFF) ");
	fflush(stdout);
	totalNumReads = 0;
	maxReadSize = 0;
	maxReadNameSize = 0;
	avgReadSize = 0;
	for(i=0;i<16;i++) c=fgetc(readsFile); // skip 16 bytes: version[4] + index_offset[8] + index_length[4]
	fread32be(&sff_number_of_reads,readsFile); // number_of_reads[4]
	for(i=0;i<2;i++) c=fgetc(readsFile); // skip 2 bytes: header_length[2]
	fread16be(&sff_key_length,readsFile); // key_length[2]
	fread16be(&sff_number_of_flows,readsFile); // number_of_flows[2]
	if(sff_number_of_flows==168){
		if(verbose) printf("(GS 20) ");
		numLinkers=1;
		strcpy(linkerSequence[0],FLX_LINKER);
		strcpy(linkerSequence[1],"");
	} else if(sff_number_of_flows==400){
		if(verbose) printf("(GS FLX Standard) ");
		numLinkers=1;
		strcpy(linkerSequence[0],FLX_LINKER);
		strcpy(linkerSequence[1],"");
	} else if(sff_number_of_flows==800){
		if(verbose) printf("(GS FLX Titanium / GS Junior) ");
		numLinkers=2;
		strcpy(linkerSequence[0],TIT_LINKER1);
		strcpy(linkerSequence[1],TIT_LINKER2);
	} else if(sff_number_of_flows==1600){
		if(verbose) printf("(GS FLX+) ");
		numLinkers=2;
		strcpy(linkerSequence[0],TIT_LINKER1);
		strcpy(linkerSequence[1],TIT_LINKER2);
	} else if(sff_number_of_flows==220 || sff_number_of_flows==260){
		if(verbose) printf("(Ion Torrent) ");
		numLinkers=2;
		strcpy(linkerSequence[0],ION_LINKER1);
		strcpy(linkerSequence[1],ION_LINKER2);
	} else {
		if(verbose) printf("(Unknown Instrument) ");
		numLinkers=0;
		strcpy(linkerSequence[0],"");
		strcpy(linkerSequence[1],"");
	}
	padding_size = ( ( 8 - ( ( 31 + sff_number_of_flows + sff_key_length ) & 7 ) ) & 7 );
	n = ( 1 + sff_number_of_flows + sff_key_length + padding_size ); // skip: flowgram_format_code[1] + flow_chars[number_of_flows] + key_sequence[key_length] + eight_byte_padding[padding_size]
	for(i=0;i<n;i++) c=fgetc(readsFile); // end of common header section
	fgetpos(readsFile,&readssectionpos); // save position of start of reads section
	number_of_reads=0;
	while( number_of_reads < sff_number_of_reads ){
		for(i=0;i<2;i++) c=fgetc(readsFile); // skip 2 bytes: read_header_length[2]
		fread16be(&sff_name_length,readsFile); // name_length[2]
		fread32be(&sff_number_of_bases,readsFile); // number_of_bases[4]
		fread16be(&sff_clip_qual_left,readsFile); // clip_qual_left[2]
		fread16be(&sff_clip_qual_right,readsFile); // clip_qual_right[2]
		fread16be(&sff_clip_adapter_left,readsFile); // clip_adapter_left[2]
		fread16be(&sff_clip_adapter_right,readsFile); // clip_adapter_right[2]
		if( sff_clip_qual_left > sff_clip_adapter_left ) start_base_pos = sff_clip_qual_left; // choose maximum value
		else start_base_pos = sff_clip_adapter_left;
		if( start_base_pos > 0 ) start_base_pos--; // convert to 0-based array
		if( sff_clip_qual_right == 0 ) sff_clip_qual_right = sff_number_of_bases; // highest right position is number_of_bases
		if( sff_clip_adapter_right == 0 ) sff_clip_adapter_right = sff_number_of_bases;
		if( sff_clip_qual_right < sff_clip_adapter_right) end_base_pos = sff_clip_qual_right; // choose minimum value
		else end_base_pos = sff_clip_adapter_right;
		end_base_pos--; // convert to 0-based array
		number_of_bases = ( end_base_pos - start_base_pos + 1 );
		readNameSize = ((int)sff_name_length);
		if( readNameSize > maxReadNameSize ) maxReadNameSize = readNameSize;
		padding_size = ( ( 8 - ( ( 16 + sff_name_length ) & 7 ) ) & 7 );
		n = ( start_base_pos ); // skip the left clipped bases
		n = ( sff_name_length + padding_size + 2*sff_number_of_flows + sff_number_of_bases + n ); // skip: name[name_length] + eight_byte_padding[padding_size] + flowgram_values[2][number_of_flows] + flow_index_per_base[number_of_bases] + (left_clip_size)
		for(i=0;i<n;i++) c=fgetc(readsFile); // end of read header section and beginning of read data section
		readSize = 0;
		for(i=0;i<number_of_bases;i++){
			c=fgetc(readsFile);
			if( c=='A' || c=='C' || c=='G' || c=='T' || c=='a' || c=='c' || c=='g' || c=='t' ) readSize++; // only count valid chars
		}
		padding_size = ( ( 8 - ( ( 2*sff_number_of_flows + 3*sff_number_of_bases ) & 7 ) ) & 7 );
		n = ( sff_number_of_bases - end_base_pos - 1 ); // skip the right clipped bases
		n = ( padding_size + sff_number_of_bases + n ); // skip: quality_scores[number_of_bases] + eight_byte_padding[padding_size] + (right_clip_size)
		for(i=0;i<n;i++) c=fgetc(readsFile);
		if(c==EOF) break;
		if( number_of_bases != 0 ){ // consider size including all chars (even invalid) and not readSize
			if( (int)number_of_bases > maxReadSize ) maxReadSize = (int)number_of_bases;
			avgReadSize += (int)number_of_bases;
		}
		number_of_reads++;
		if(!preprocess) break; // if we do not want to preprocess the entire reads file, break after the first read
	}
	totalNumReads = ((int)number_of_reads);
	if( totalNumReads == 0 ){
		printf("\n> ERROR: No valid data in reads file\n");
		exit(0);
	}
	fsetpos(readsFile,&readssectionpos); // back to beginning of reads section
	sffType='S';
	BuildSffLinkerMasks();
	if( CheckForPairedEndsSff() ){  // take the first reads from the SFF file and search for the linker to check if the file contains paired-end reads or not
		sffType='s';
		if(verbose) printf("(Paired-Ends) ");
	}
	fsetpos(readsFile,&readssectionpos); // back to beginning of reads section
	sffReadsCount=0; // reset reads count for the current SFF file
	if(!preprocess) return sffType;
	if(verbose){
		avgReadSize = ( avgReadSize / totalNumReads );
		PrintNumber((long long int)totalNumReads);
		printf(" reads (average size = %dbp) ",avgReadSize);
		fflush(stdout);
	}
	return sffType;
}

// Splits a paired-end read from a paired-end SFF file into both read mates
int LoadNextReadPairFromSff(){
	int i, n, readSize;
	char *charPtr;
	currentPair=0; // needed if IncreaseMaxReadArraySize is called
	readsFile=pairFiles[0];
	read=fwdReads[0];
	readName=readNames[0];
	n=-1;
	while(n==-1){ // keep searching for the next read with the linker
		readSize=LoadNextReadFromSff();
		if(readSize==0){ // no more reads in file
			readsSizes[0]=0;
			readsSizes[1]=0;
			return readSize;
		}
		//n=FindExactSffLinker(readSize);
		n=FindInexactSffLinker(readSize); // returns -1 if the linker was not found, otherwise it sets the linker start and end positions in the read
		if( linkerStartPosInRead==0 || linkerEndPosInRead==(readSize-1) ) n=-1; // check if both reads have a non-zero length
		if(allowUnpairedReads) break; // if we accept single mates, all reads are kept
	}
	fwdReads[0][linkerStartPosInRead]='\0'; // stop 1st read before linker starts
	readsSizes[0]=linkerStartPosInRead;
	charPtr=(char *)( fwdReads[0] + (linkerEndPosInRead+1) );
	i=0;
	while((*charPtr)!='\0'){ // copy all chars following the end of the linker to the 2nd read
		fwdReads[1][i++]=(*charPtr);
		charPtr++;
	}
	fwdReads[1][i]='\0';
	readsSizes[1]=i; // (readSize-(linkerEndPosInRead-linkerStartPosInRead+1))
	strcpy(readNames[1],readName); // copy name of 1st read to name of 2nd read
	strcat(readNames[0],"_left"); // append identifiers to both read names
	strcat(readNames[1],"_right");
	return readSize;
}

void SplitPairedEndSff(char *sfffilename){
	FILE *sffFile, *fastaFiles[3];
	char *fastaFilenames[3];
	unsigned int pairedCharsCount, unpairedCharsCount, pairedReadsCount, unpairedReadsCount, exactLinkerCount, partialLinkerCount;
	int i, n;
	printf("> Opening SFF file <%s> ... ",sfffilename);
	fflush(stdout);
	if((sffFile=fopen(sfffilename,"rb"))==NULL){
		printf("\n> ERROR: SFF file not found\n");
		exit(0);
	}
	readsFile=sffFile;
	readsFilename=sfffilename;
	AnalyzeSffReads(1,1);
	fastaFilenames[0]=AppendToBasename(sfffilename,"_left.fasta");
	fastaFilenames[1]=AppendToBasename(sfffilename,"_right.fasta");
	fastaFilenames[2]=AppendToBasename(sfffilename,"_unpaired.fasta");
	for(i=0;i<3;i++){
		if((fastaFiles[i]=fopen(fastaFilenames[i],"w"))==NULL){
			printf("\n> ERROR: Cannot create Fasta file <%s>\n",fastaFilenames[i]);
			exit(0);
		}
	}
	printf("OK\n");
	printf("> Splitting reads ... ");
	fflush(stdout);
	allowUnpairedReads=1; // allow reads even when one mate is missing
	numPairs=2;
	maxReadArraySize = 1;
	maxReadNameArraySize = 255;
	for(i=0;i<numPairs;i++){
		fwdReads[i] = (char *)malloc((maxReadArraySize+1)*sizeof(char));
		revReads[i] = (char *)malloc((maxReadArraySize+1)*sizeof(char)); // not used, but needed because IncreaseMaxReadArraySize will try to realloc it
		readNames[i] = (char *)malloc((maxReadNameArraySize+1)*sizeof(char));
	}
	currentPair=0;
	pairFiles[0]=sffFile;
	pairFiles[1]=NULL;
	pairedReadsCount=0;
	unpairedReadsCount=0;
	exactLinkerCount=0;
	partialLinkerCount=0;
	pairedCharsCount=0;
	unpairedCharsCount=0;
	while((n=LoadNextReadPairFromSff())!=0){ // while we do not reach the end of the file
		if(readsSizes[0]>0 && readsSizes[1]>0){ // both mates are present
			for(i=0;i<numPairs;i++){
				fprintf(fastaFiles[i],">%s\n%s\n",readNames[i],fwdReads[i]); // save both mates
				//printf("\n>%s\n%s\n",readNames[i],fwdReads[i]);
				pairedCharsCount+=readsSizes[i];
			}
			if(numLinkerErrors==0) exactLinkerCount++; // count number of reads with exact linker
			pairedReadsCount++;
		} else { // only one of the mates is present
			for(i=0;i<numPairs;i++){
				n=( (int)strlen(readNames[i]) - ((i==0)?(5):(6)) ); // remove "_left" or "_right" suffix from read name
				readNames[i][n]='\0';
				if(readsSizes[i]>0) fprintf(fastaFiles[2],">%s\n%s\n",readNames[i],fwdReads[i]); // save unpaired mate
				//if(readsSizes[i]>0) printf("\n>%s\n%s\n",readNames[i],fwdReads[i]);
				unpairedCharsCount+=readsSizes[i];
			}
			if(linkerStartPosInRead<linkerEndPosInRead) partialLinkerCount++; // count number of reads with partial linker (if not found entirely, it would be start>end)
			unpairedReadsCount++;
		}
	}
	if(pairedReadsCount==0){
		printf("\n> ERROR: No reads found with paired-ends linker\n");
		exit(0);
	}
	printf("OK\n");
	printf(":: ");
	PrintUnsignedNumber(pairedReadsCount);
	printf(" paired reads (");
	PrintUnsignedNumber(exactLinkerCount);
	printf(" with exact linker ; avg %u bp/mate)\n",(pairedCharsCount/(2*pairedReadsCount)));
	printf(":: ");
	PrintUnsignedNumber(unpairedReadsCount);
	printf(" unpaired reads (");
	PrintUnsignedNumber(partialLinkerCount);
	printf(" with partial linker ; avg %u bp/read)\n",(unpairedCharsCount/unpairedReadsCount));
	fclose(sffFile);
	for(i=0;i<numPairs;i++){
		printf("> Saving ");
		PrintUnsignedNumber(pairedReadsCount);
		printf(" %s reads to <%s> ... ",((i==0)?"left":"right"),fastaFilenames[i]);
		fclose(fastaFiles[i]);
		free(fastaFilenames[i]);
		printf("OK\n");
		free(fwdReads[i]);
		free(revReads[i]);
		free(readNames[i]);
	}
	printf("> Saving ");
	PrintUnsignedNumber(unpairedReadsCount);
	printf(" unpaired reads to <%s> ... ",fastaFilenames[2]);
	fclose(fastaFiles[2]);
	free(fastaFilenames[2]);
	printf("OK\n");
	printf("> Done!\n");
	exit(0);
}
