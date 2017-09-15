#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/timeb.h>
#include <limits.h>

#include "version.h"
#include "alignreads.h"
#include "bwtindex.h"
#include "fileformats.h"
#include "samoutput.h"
#include "dynprog.h"
#include "tools.h"
#include "console.h"
#include "graphics.h"
#include "consensus.h"

#include "splicesites.h"

//#define DEBUG 1
//#define DISABLE_SS_ARRAY 1

#ifdef _MSC_VER
	#pragma warning(disable:4996)
#endif

#define STATSINTERVAL 2.0

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

/**************/
/* STRUCTURES */
/**************/

typedef struct _SeedPos { // structure that stores the position and the corresponding seed number
	unsigned int pos; // position of the seed in the reference genome
	int readpos; // position of the seed in the read
	int seed; // corresponding seed id
	int size; // size of the seed
	int ssId; // splice site id used to create this seed, if any
} SeedPos;

typedef struct _ReadHit { // structure that stores an occurrence of the read in the genome
	unsigned int pos; // starting position of the read in the reference genome
	int score; // score of the alignment
	int numerrors; // number of errors of the alignment
	int cigarsize; // size of the CIGAR string of the alignment
	int *cigarcounts; // counts for the CIGAR string of the alignment
	char *cigarcodes; // codes for the CIGAR string of the alignment
	char strand; // strand of the read
	int matedist; // distance from this mate to the other mate (in paired-ends mode)
} ReadHit;

typedef struct _Seed { // structure that stores information about the seeds of the read in the genome
	int startPosInRead;
	int size;
	unsigned int numHits;
	unsigned int bwtTopPointer;
} Seed;

typedef struct _Strand {
	char *chars;
	int numSeeds;
	Seed *seeds;
	unsigned int *topPointers;
	unsigned int *bottomPointers;
} Strand;

typedef struct _Read {
	char *name;
	int size;
	char type;
	Strand fwdStrand;
	Strand revStrand;
} Read;

int maxNumSeeds = 0; // maximum number of seeds that the array (Read->Strand->seeds) can currently store


/********************/
/* GLOBAL VARIABLES */
/********************/

FILE *outputFile; // file to write results
FILE *unalignedReadsFile; // file to write unaligned reads
FILE *errorsFile;

int pairedEndReadsMode; // indicates if we are in paired-end mode or not
int rnaSequencingMode; // indicates if we are in RNA/Transcriptome sequencing mode or not

int reportUnalignedReads;
int reportAllHits;
int singleRefReadsOnly;

int numReadsFiles; // number of reads files
int currentReadsFile; // index of the current reads file
FILE **readsFiles; // files with the reads to be aligned
char **readsFilenames; // names of the reads files
char *readsFileTypes; // types of the reads files
int *minDists; // mininum distances for read pairs for each reads files
int *maxDists; // maximum distances for read pairs for each reads files

int minIntronLength; // minimum distance between exons
int maxIntronLength; // maximum distance between exons
unsigned int numSplicingEvents;
unsigned int numSplicedReads;
unsigned int numSeedsAddedBySpliceSites; // number of created seeds on the the opposite side of an intron, because a splice site was detected

unsigned int sizeBWT; // reference genome index size
unsigned int sizeRefGenome; // reference genome size

char *refGenome; // reference genome string
char **refsNames; // reference genome(s) name(s)
int numRefs; // number of references
unsigned int *refsEndPositions; // ending positions of each ref in the global ref string
char *genomefilename; // name of the file containing the un-indexed reference genome

unsigned int numReads; // number of reads 
unsigned int numAlignedReads; // number of aligned reads
unsigned int numUnalignedReads; // number of unaligned reads
unsigned int numAlignedHardReads; // number of hard reads aligned
unsigned int numHighRepetitionReads; // number of reads with too many hits
unsigned int maxReadHits; // maximum number of hits observed for a read

unsigned int numAlignedBothMates; // number of pairs where both reads were aligned
unsigned int numAlignedSingleMates; // number of pairs where only one read was aligned
unsigned int minObservedPairDistance; // minimum distance observed between two mate reads
unsigned int maxObservedPairDistance; // maximum distance observed between two mate reads
unsigned int avgObservedPairDistances; // aproximate average of all observed distances between the two mate reads
unsigned int weightObservedPairDistances; // weight of the aproximate average

int totalNumSeeds; // total number of seeds in all the reads
int totalNumSeedOccs; // total number of occurrences of all the found seeds
int totalNumSeedChains; // total number of best seed chains tested
unsigned long long int totalNumReadHits; // count of all hits for all reads
unsigned long long int totalNumDPErrors; // total number of DP errors in all the aligned reads
unsigned long long int totalNumReadsChars; // total number of chars in all the processed reads
unsigned long long int totalNumMappedReadsChars; // total number of chars in all the mapped reads
int maxReadSize; // size of the longest read in the reads file (always lower or equal to maxReadArraySize)

ReadHit *readsHits[2]; // stores all hits in the genome for both reads in the paired-end read set
int maxReadHitsArraySize; // maximum number of hits that the read hits array can store
int numReadsHits[2]; // number of hits for each one of the paired-end reads
int bestReadsHitIds[2]; // ids of the top scoring hits for each one of the paired-end reads
int bestReadsScores2nd[2]; // second best alignment score for each one of the paired-end reads
int minReadPairDistance; // minimum allowed distance between the read pair
int maxReadPairDistance; // maximum allowed distance between the read pair
int reqReadPairStrands; // required read pair strands (1=same or 2=opposite)

int minIdentityPct; // minimum allowed identity percentage for each read
int maxNumErrors; // maximum number of errors allowed in the banded dynamic programming among all reads
int maxSeedHits; // maximum allowed number of occurences of the seed in the reference genome

char *fwdRead; // read chars for the forward strand
char *revRead; // read chars for the reverse strand

struct timeb starttb, endtb; // variables for storing time values
long long int prevFilesData, totalFilesData, currentFileData; // variables for storing relative position inside reads file

char (*AnalyzeReads)(int,int); // function pointer for specific file format reader function
int (*LoadNextRead)(); // function pointer for specific file format reader function


/********************************/
/* FUNCTIONS FOR ALIGNMENT MODE */
/********************************/


// returns elapsed time since the given time value
double GetElapsedTime(struct timeb starttb){
	ftime( &endtb );
	return ((endtb.time) + (endtb.millitm)/1000.0) - ((starttb.time) + (starttb.millitm)/1000.0);
}

// Compare function for quicksort (sort seed positions)
int CompareSeedPositions(const void *a, const void *b){
	/*
	const SeedPos *aa = *(const SeedPos *const *) a;
	const SeedPos *bb = *(const SeedPos *const *) b;
	*/
	SeedPos *aa = *(SeedPos **)a;
	SeedPos *bb = *(SeedPos **)b;
	//return ( (aa->pos) - (bb->pos) );
	if( (aa->pos) > (bb->pos) ) return +1;
	if( (aa->pos) < (bb->pos) ) return -1;
	return ( (aa->seed) - (bb->seed) ); // if( (aa->pos) == (bb->pos) ) ; when the position is the same, sort by seed id
}

// Compare function for quicksort (sort read hit positions)
int CompareReadHitPositions(const void *a, const void *b){
	return ( ((ReadHit *)a)->pos - ((ReadHit *)b)->pos );
}

// expand the read chars arrays for both forward and reverse strands to store a read with the new larger size
void IncreaseMaxReadArraySize(int newsize){
	int i;
	if(newsize<maxReadArraySize) return;
	while(maxReadArraySize<newsize) maxReadArraySize += maxReadArraySize; // double the array's size each time
	for( i = 0 ; i < numPairs ; i++ ){
		fwdReads[i]=(char *)realloc(fwdReads[i],(maxReadArraySize+1)*sizeof(char)); // realloc arrays
		revReads[i]=(char *)realloc(revReads[i],(maxReadArraySize+1)*sizeof(char));
	}
	fwdRead = fwdReads[currentPair]; // set pointers to current reads file
	revRead = revReads[currentPair];
	read=fwdRead; // we always load the read chars from the reads file into the forward strand
}

// expand the read hits arrays to store a the new number of hits
void IncreaseMaxReadHitsArraySize(int newsize, int cigararraysize){
	int i,j,n;
	n = maxReadHitsArraySize; // old size
	maxReadHitsArraySize = newsize;
	for(i=0;i<numPairs;i++){
		readsHits[i]=(ReadHit *)realloc(readsHits[i],maxReadHitsArraySize*sizeof(ReadHit));
		for(j=n;j<maxReadHitsArraySize;j++){ // only alloc CIGAR arrays on newly created read hit entries
			readsHits[i][j].cigarcodes = (char *)malloc(cigararraysize*sizeof(char));
			readsHits[i][j].cigarcounts = (int *)malloc(cigararraysize*sizeof(int));
		}
	}
}

// loads the reverse complementary of the read
void LoadReverseStrandRead(char *read, char *revread, int readsize){
	char c;
	int i, j;
	j = 0;
	for( i = (readsize-1) ; i >= 0 ; i-- ){
		c = read[i];
		switch( c ){
			case 'A':
				revread[j] = 'T';
				break;
			case 'C':
				revread[j] = 'G';
				break;
			case 'G':
				revread[j] = 'C';
				break;
			case 'T':
				revread[j] = 'A';
				break;
			default:
				break;
		}
		j++;
	}
	revread[readsize] = '\0';
}

// updates the Cigar string with the given operation
void UpdateCigarString(char opcode, int opsize){
	#ifdef DEBUGDP
	int i;
	#endif
	if( cigarStrCodes[cigarStrSize] == opcode ){ // if the previous operation was already the same
		cigarStrCounts[cigarStrSize] += opsize;
	} else { // if not, create new operation
		cigarStrSize++; // holds the last position of the array, size is +1
		cigarStrCodes[cigarStrSize] = opcode;
		cigarStrCounts[cigarStrSize] = opsize;
	}
	#ifdef DEBUGDP
	i = 0;
	if( cigarStrCounts[0] == 0 ) i = 1;
	for( ; i <= cigarStrSize ; i++ ) fprintf(debugfile,"%d%c",cigarStrCounts[i],cigarStrCodes[i]);
	fputc('\n',debugfile);
	#endif
}

// NOTE: currently unused because the CIGAR string is now updated directly inside the main DP function while backtracking
// updates the Cigar string with the latest DP alignment
void UpdateCigarStringWithDP(){
	int i;
	for( i = 0 ; i < dpTraceSize ; i++ ){ // scan all dp trace directions (the array is already in the correct order)
		if( cigarStrCodes[cigarStrSize] == dpTrace[i] ){ // if it's the same operation as the previous one
			cigarStrCounts[cigarStrSize]++;
		} else { // if not, create new operation
			cigarStrSize++; // holds the last position of the array, size is +1
			cigarStrCodes[cigarStrSize] = dpTrace[i];
			cigarStrCounts[cigarStrSize] = 1;
		}
	}
	#ifdef DEBUGDP
	i = 0;
	if( cigarStrCounts[0] == 0 ) i = 1;
	for( ; i <= cigarStrSize ; i++ ) fprintf(debugfile,"%d%c",cigarStrCounts[i],cigarStrCodes[i]);
	fputc('\n',debugfile);
	#endif
}

// TODO: move to "samoutput.c" file
// NOTE: single reference is labeled as "R" and multiple references as "R#"
void PrintSAMHeader(){
	int i,n;
	if(numRefs==1){ // only one ref
		fprintf(outputFile,"@SQ\tSN:%s\tLN:%u\tSP:%s\tUR:%s\n","R",sizeRefGenome,refsNames[0],genomefilename);
	} else { // multiple refs
		for(i=0;i<numRefs;i++){
			if(i==0) n=(refsEndPositions[i]+1);
			else n=(refsEndPositions[i]-refsEndPositions[(i-1)]); // size of the ref
			fprintf(outputFile,"@SQ\tSN:%s%d\tLN:%u\tSP:%s\tUR:%s\n","R",(i+1),n,refsNames[i],genomefilename);
		}
	}
	fprintf(outputFile,"@PG\tID:0\tPN:TAPyR\tVN:%s\n",VERSION);
}

// Returns the reference where this read hit is mapped on
int GetReadRef(unsigned int rpos){
	int start, end, middle;
	//while( (i<numRefs) && ( readPos > refsEndPositions[i] ) ) i++;
	start = -1;
	end = ( numRefs - 1 ); // the correct ref number will be stored in here
	while( (end-start) != 1 ){ // binary search for correct ref
		middle = ( ( start + end ) / 2 );
		if( rpos > refsEndPositions[middle] ) start = middle;
		else end = middle;
	}
	return end; // current ref
}

// TODO: update the number of errors of the read if insertions are added
// When there are multiple reference sequences, finds the correct read reference and updates the read position and its CIGAR string accordingly
int FixMultiRefRead(unsigned int *readPos, int readSize){
	unsigned int n, rpos, gpos;
	int i, ref;
	char c;
	rpos = (*readPos);
	ref = GetReadRef(rpos);
	n = refsEndPositions[ref]; // last valid position in the current ref
	if( (rpos+readSize-1) <= n ) return ref; // if the read does not extend to the next ref, it is ok
	rpos = 0; // position inside the read
	gpos = (*readPos); // position inside the genome
	i = 0; // position in the CIGAR string
	c = 'M'; // set variable just so compiler doesn't complain
	while( ( (gpos-1) <= n ) && ( i <= cigarStrSize ) ){ // gpos points to the position next to the last filled one
		c = cigarStrCodes[i];
		if( c=='M' || c=='=' || c=='X' ){ // advances positions in both read and ref
			rpos += cigarStrCounts[i];
			gpos += cigarStrCounts[i];
		} else if( c=='I' || c=='S' ) rpos += cigarStrCounts[i]; // advance pos only in read
		else if( c=='D' || c=='N' ) gpos += cigarStrCounts[i]; // advance pos only in ref
		i++;
	}
	if( (gpos-1) <= n ) return ref; // if the read extends to the next ref we need to fix the CIGAR string, otherwise it is ok
	gpos--; // now both variables point to the last filled position
	rpos--;
	n = ( gpos - n ); // how many chars did we advance into the next ref
	i--; // last CIGAR operation (that crossed to the next ref)
	cigarStrCounts[i] -= n; // break operation before it enters the next ref
	if( c=='M' || c=='=' || c=='X' ) rpos -= n; // if the operation also changed the read pos, change it back
	if( rpos >= (unsigned int)(readSize/2) ){ // if the largest part of the read is in the current ref
		n = ( readSize - rpos - 1 ); // number of chars inside the next ref until the end of the read
		i++;
		cigarStrCounts[i] = n; // consider the rest of the read chars as insertions
		cigarStrCodes[i] = 'I';
		cigarStrSize = i;
	} else { // if the largest part of the read is in the next ref
		cigarStrCounts[i] = n; // set the last operation (M or D) to account only for the chars inside the new ref
		rpos++; // first read char pos inside the new ref
		if( i>0 ) i--; // previous operation
		else { // if this is the first operation, we need to create a new one before, so, shift all operations forward
			i = (cigarStrSize + 1);
			while( i>0 ){
				cigarStrCounts[i] = cigarStrCounts[(i-1)];
				cigarStrCodes[i] = cigarStrCodes[(i-1)];
				i--;
			}
			cigarStrSize++; // one more operation
		}
		cigarStrCounts[i] = rpos; // consider all the read chars behind as insertions
		cigarStrCodes[i] = 'I';
		while( i>0 ){ // reset all the other operations behind
			i--;
			cigarStrCounts[i] = 0;
			cigarStrCodes[i] = 'M';
		}
		(*readPos) = ( refsEndPositions[ref] + 1 ); // set read to start in the 1st pos of the next ref
		ref++; // next ref
	}
	return ref; // return the ref number where the read starts
}

// TODO: in case of multi refs, update pos and ref directly in readHits structure
// TODO: update the matedist variable with the size of the other mate (and use matedist field in ReadHit structure)
// TODO: better calculation of mate pair: max(p1+r1,p2+r2) - min(p1,p2)
// TODO: if the mate maps to another reference, output the correct RNEXT field and fix the position in the PNEXT field
// TODO: move to "samoutput.c" file
// NOTE: when readPos is -1 , it means that the read was not mapped
// NOTE: cigarStrSize points to the last filled entry of the CIGAR arrays (so, real size is +1)
// FORMAT: QNAME FLAG RNAME POS MAPQ CIGAR MRNM MPOS ISIZE SEQ QUAL
void PrintReadInSAMFormat(char *read, char *readName, char strand, unsigned int readPos, int numErrors, int bestScore, int bestScore2nd, int readSize){
	int i, flag, mapqual;
	int matedist, absmatedist;
	unsigned int pos, matepos;
	#ifdef DEBUG
	unsigned int pg, originalreadpos;
	int pr, j, n;
	char c;
	static int debugstringssize = 0;
	static char *cigarstring = NULL;
	static char *refstring = NULL;
	static char *readstring = NULL;
	originalreadpos = readPos;
	#endif
	matepos = 0; // initialize mate variables here to avoid compiler warnings
	matedist = 0;
	flag = 0;
	if( pairedEndReadsMode ){ // paired-end reads mode
		flag += 1; // flag = 1 means it is a paired-end read
		if( currentPair == 0 ) flag += 64; // flag = 64 means first read in a pair
		else flag += 128; // flag = 128 means second read in a pair
		i = ((currentPair + 1) & 1); // id of the mate read
		if( numReadsHits[i] > 0 ){ // if the mate read was mapped
			matepos = readsHits[i][bestReadsHitIds[i]].pos; // MPOS = mate position
			if( numReadsHits[currentPair] > 0 ){ // if this read is mapped too
				flag += 2; // flag = 2 means both reads are mapped
				pos = readsHits[currentPair][bestReadsHitIds[currentPair]].pos; // pos of this mate
				absmatedist = (int) ((matepos>pos)?(matepos-pos):(pos-matepos)); // abs( matepos - pos );
				if( (maxReadPairDistance>0) && ((absmatedist<minReadPairDistance) || (absmatedist>maxReadPairDistance)) ) // if a pair distance is defined but this distance is not valid
					absmatedist = (int)(sizeRefGenome - absmatedist); // set the extra case when the reads "wrap around" from the end to the beginning of the (circular) genome
				if( matepos > pos ) matedist = absmatedist; // fix mate distance signal according to SAM format
				else matedist = (-absmatedist);
				/*
				matedist = ( matepos - readsHits[currentPair][bestReadsHitIds[currentPair]].pos ); // TLEN = template length (other mate's pos minus this mate's pos ; distance covering the full length of both mates)
				if(matedist>0) matedist += readSize; // the other mate is to the right, add the other mate's size
				else matedist -= readSize; // the other mate is to the left, add this mate's size
				*/
				if(currentPair==0){ // update count and stats, but only one time per pair
					if( (unsigned int)absmatedist < minObservedPairDistance ) minObservedPairDistance = absmatedist;
					if( (unsigned int)absmatedist > maxObservedPairDistance ) maxObservedPairDistance = absmatedist;
					if( ((unsigned int)absmatedist > (avgObservedPairDistances/2) ) && ((unsigned int)absmatedist < (3*avgObservedPairDistances/2) ) ){ // distance within the current valid interval range (+/-50% the average)
						avgObservedPairDistances = ( (weightObservedPairDistances*avgObservedPairDistances) + absmatedist ) / (weightObservedPairDistances+1); // update current average
						weightObservedPairDistances++; // increase weight
					} else { // distance outside the current valid range
						weightObservedPairDistances--; // decrease weight
						if(weightObservedPairDistances==0){ // if the current average lost all its weight, set a new average
							avgObservedPairDistances = absmatedist;
							weightObservedPairDistances = 1;
						}
					}
					numAlignedBothMates++;
				}
			} else {
				numAlignedSingleMates++; // the mate was mapped, but this one was not
			}
			if( readsHits[i][bestReadsHitIds[i]].strand == '-' ) flag += 32; // flag = 32 means the other read mate is mapped in the reverse strand
		} else flag += 8; // flag = 8 means the other read mate is unmapped
	}
	if( readPos == (unsigned int) -1 ){ // if no positions were found for the read
		flag += 4; // flag = 4 means unmapped read
		fprintf(outputFile,"%s\t%d\t*\t0\t0\t*\t*\t0\t0\t*\t*\n",readName,flag);
		return;
	}
	if( strand == '-' ) flag += 16; // flag = 16 means reverse strand
	fprintf(outputFile,"%s\t%d\t",readName,flag);
	pos = readPos; // position in first or only ref
	if( numRefs == 1 ){
		//fprintf(outputFile,"R\t"); // single ref labeled only as "R"
		fprintf(outputFile,"%s\t",refsNames[0]);
	} else { // update read pos in case of multiple refs
		i = FixMultiRefRead((unsigned int *)(&readPos),readSize); // get the reference number
		//fprintf(outputFile,"R%d\t",(i+1));
		fprintf(outputFile,"%s\t",refsNames[i]);
		if( i != 0 ) pos = ( readPos - refsEndPositions[(i-1)] - 1 ); // change position to match the new reference
	}
	fprintf(outputFile,"%u\t",(pos+1)); // pos+1 because in SAM, genome starts at pos 1, not 0
	if( bestScore == 0 ) mapqual = 255; // unavailable mapping quality
	else {
		mapqual = 250 * ( bestScore - bestScore2nd ) / bestScore; // = 250 * c1 * c2 * ( S1 - S2 ) / S1
		if(mapqual<0) mapqual=0; // prevent negative values
	}
	fprintf(outputFile,"%d\t",mapqual); // mapping quality
	i = 0;
	while( cigarStrCounts[i] == 0 ) i++; // skip first 0 size match if it exists
	for( ; i <= cigarStrSize ; i++ ) fprintf(outputFile,"%d%c",cigarStrCounts[i],cigarStrCodes[i]); // CIGAR string
	if( !pairedEndReadsMode ) fprintf(outputFile,"\t*\t0\t0\t"); // normal read
	else fprintf(outputFile,"\t=\t%d\t%d\t",matepos,matedist); // RNEXT, PNEXT, TLEN
	fprintf(outputFile,"%s\t*",read); // read chars and base qualities
	#ifdef DEBUG
	fprintf(outputFile,"\tNM:i:%d",numErrors);
	#else
	numErrors = 0; // just so compiler does not complain about unused variable
	#endif
	fprintf(outputFile,"\n");
	#ifdef DEBUG
	n = 0;
	for( i = 0 ; i <= cigarStrSize ; i++ ) n += cigarStrCounts[i]; // get size of CIGAR strings
	if( cigarstring==NULL || refstring==NULL || readstring==NULL ){ // alloc space for the strings the first time
		debugstringssize = n;
		cigarstring = (char *)malloc((debugstringssize+1)*sizeof(char));
		refstring = (char *)malloc((debugstringssize+1)*sizeof(char));
		readstring = (char *)malloc((debugstringssize+1)*sizeof(char));
	} else if( n > debugstringssize ){ // realloc the space if needed
		debugstringssize = n;
		cigarstring = (char *)realloc(cigarstring,(debugstringssize+1)*sizeof(char));
		refstring = (char *)realloc(refstring,(debugstringssize+1)*sizeof(char));
		readstring = (char *)realloc(readstring,(debugstringssize+1)*sizeof(char));
	}
	pg = readPos; // position in genome
	pr = 0; // position in read
	j = 0; // position in debug strings
	for( i = 0 ; i <= cigarStrSize ; i++ ){
		n = cigarStrCounts[i];
		if( n == 0 ) continue;
		c = cigarStrCodes[i];
		if( c == 'M' || c == '=' || c== 'X' ) while( n > 0 ){ // match or mismatch
			cigarstring[j] = c;
			refstring[j] = refGenome[pg++];
			readstring[j] = read[pr++];
			j++;
			n--;
		}
		else if( c == 'N' ) while( n > 0 ){ // intron
			cigarstring[j] = c;
			refstring[j] = refGenome[pg++];
			readstring[j] = ' '; // empty space
			j++;
			n--;
		}
		else if( c == 'I' ) while( n > 0 ){ // insertion into reference
			cigarstring[j] = c;
			refstring[j] = '-';
			readstring[j] = read[pr++];
			j++;
			n--;
		}
		else if( c == 'D' ) while( n > 0 ){ // deletion from reference
			cigarstring[j] = c;
			refstring[j] = refGenome[pg++];
			readstring[j] = '-';
			j++;
			n--;
		}
		else if( c == 'S' ) while( n > 0 ){ // soft clipping of read
			cigarstring[j] = c;
			refstring[j] = refGenome[pg++];
			readstring[j] = (char)(read[pr++]+32); // convert to lowercase
			j++;
			n--;
		}
	}
	cigarstring[j] = '\0';
	refstring[j] = '\0';
	readstring[j] = '\0';
	fprintf(debugfile,"%s\n%s\n%s\n",cigarstring,refstring,readstring);
	if((numRefs==1) || (originalreadpos==readPos)){ // if there are multiple references and the read falls between two of them, the number of errors is not updated correctly yet
		n = 0; // number of errors retrieved from the CIGAR string alignment
		j--;
		for( ; j >= 0 ; j-- ){
			c = cigarstring[j];
			if( c == 'S' || c == 'N' ) continue; // these two operations do not count as errors
			if( refstring[j] != readstring[j] ) n++;
		}
		for( i = 0 ; i <= cigarStrSize ; i++ ){
			j = cigarStrCounts[i];
			if( (j <= 0) || ((j > readSize) && cigarStrCodes[i]!='N') ){ // the intron operation can have a length longer than the read size
				if( i==0 && j==0 ) continue; // there can be a valid 0 sized match at the first position
				break;
			}
		}
		if( (n!=numErrors) || (pr!=readSize) || (i!=(cigarStrSize+1)) ){
			if(errorsFile==NULL) errorsFile=fopen("errors.txt","w");
			fprintf(errorsFile,">%s",readName);
			if(n!=numErrors) fprintf(errorsFile," (%d/%d)",n,numErrors); // check if the found number of errors matches the one calculated before
			if(pr!=readSize) fprintf(errorsFile," [%d/%d]",pr,readSize); // check if the output read size matches the original read size
			if(i!=(cigarStrSize+1)) fprintf(errorsFile," (%d@%d)",j,i); // check if there is any incorrect cigar count value
			fprintf(errorsFile,"\n%s\n",read);
		}
	}
	#endif
}

// TODO: does it need a faster algorithm to find (all the combinations of the) pairs at valid distances ?
// NOTE: if at least one of the min/max distances is set, any single read or any read pair but with an invalid distance will not be reported
// outputs both reads of the paired-end read pair to the SAM file
void PrintReadPair(){
	int i, j, k, n, m;
	int absDist, bestScoreSum;
	char c;
	ReadHit *hit;
	if( maxReadPairDistance>0 ){ // if the limits are set, check which pairs are at a valid distance
		if( numReadsHits[0]==0 || numReadsHits[1]==0 ) return; // if at least one of the reads does not have hits, do not report any of them
		/*
		if(numReadsHits[0]>1) qsort(readsHits[0],numReadsHits[0],sizeof(ReadHit),CompareReadHitPositions); // if it has multiple hits, sort them by position
		if(numReadsHits[1]>1) qsort(readsHits[1],numReadsHits[1],sizeof(ReadHit),CompareReadHitPositions);
		*/
		m = numReadsHits[0]; // save number of hits
		n = numReadsHits[1];
		numReadsHits[0] = 0; // invalidate all previous hits
		numReadsHits[1] = 0;
		bestReadsHitIds[0] = -1;
		bestReadsHitIds[1] = -1;
		bestScoreSum = -1; // best combined score of both reads hits
		for( i = 0 ; i < m ; i++ ){
			for( j = 0 ; j < n ; j++ ){
				//absDist = abs( readsHits[1][j].pos - readsHits[0][i].pos );
				absDist = (int)( (readsHits[1][j].pos > readsHits[0][i].pos)?(readsHits[1][j].pos - readsHits[0][i].pos):(readsHits[0][i].pos - readsHits[1][j].pos) );
				if( (absDist < minReadPairDistance) || (absDist > maxReadPairDistance) ){ // less than the min dist or more than the max dist is not acceptable
					absDist = (sizeRefGenome - absDist); // check the extra case when the reads "wrap around" from the end to the beginning of the (circular) genome
					if( (absDist < minReadPairDistance) || (absDist > maxReadPairDistance) ) continue;
				}
				if( reqReadPairStrands>0 ){ // if we have requirements about the mates strands, check them too
					if(reqReadPairStrands==1 && (readsHits[1][j].strand)!=(readsHits[0][i].strand) ) continue; // same strand
					if(reqReadPairStrands==2 && (readsHits[1][j].strand)==(readsHits[0][i].strand) ) continue; // opposite strands
				}
				readsHits[0][i].matedist = absDist;
				readsHits[1][j].matedist = absDist;
				numReadsHits[0]++;
				numReadsHits[1]++;
				k = ( readsHits[0][i].score + readsHits[1][j].score ); // score of this hit pair
				if( k > bestScoreSum ){ // only set the best hit if it has a better score
					bestScoreSum = k;
					bestReadsHitIds[0] = i;
					bestReadsHitIds[1] = j;
				}
				if( reportAllHits ){ // if we want all hits at valid distances, report them now
					hit = &(readsHits[0][i]); // hit for 1st mate
					cigarStrCodes = hit->cigarcodes;
					cigarStrCounts = hit->cigarcounts;
					cigarStrSize = hit->cigarsize;
					currentPair = 0;
					c = hit->strand;
					PrintReadInSAMFormat(((c=='+')?(fwdReads[0]):(revReads[0])),readNames[0],c,hit->pos,hit->numerrors,hit->score,bestReadsScores2nd[0],readsSizes[0]);
					hit = &(readsHits[1][j]); // hit for 2nd mate
					cigarStrCodes = hit->cigarcodes;
					cigarStrCounts = hit->cigarcounts;
					cigarStrSize = hit->cigarsize;
					currentPair = 1;
					c = hit->strand;
					PrintReadInSAMFormat(((c=='+')?(fwdReads[1]):(revReads[1])),readNames[1],c,hit->pos,hit->numerrors,hit->score,bestReadsScores2nd[1],readsSizes[1]);
				}
			} // loop for 2nd mate
		} // loop for 1st mate
		if( numReadsHits[0] == 0 ) return; // if no pairs were found at a valid distance, do not report any hits
		if( reportAllHits ) return; // if we wanted all hits, they were all reported above already
	} else if( reportAllHits ){ // if no limits have been set but we want all the hits
		for( i = 0 ; i < numPairs ; i++ ){ // report all hits for both reads
			currentPair = i;
			n = numReadsHits[i];
			for( j = 0 ; j < n ; j++ ){ // all hits of this read
				hit = &(readsHits[i][j]);
				cigarStrCodes = hit->cigarcodes;
				cigarStrCounts = hit->cigarcounts;
				cigarStrSize = hit->cigarsize;
				c = hit->strand;
				PrintReadInSAMFormat(((c=='+')?(fwdReads[i]):(revReads[i])),readNames[i],c,hit->pos,hit->numerrors,hit->score,bestReadsScores2nd[i],readsSizes[i]);
			}
		}
		currentPair = 1;
		return;
	}
	for( i = 0 ; i < numPairs ; i++ ){ // report the best hit only
		if( numReadsHits[i] == 0 ) continue; // cannot report the best hit if none exists
		hit = &(readsHits[i][bestReadsHitIds[i]]);
		cigarStrCodes = hit->cigarcodes;
		cigarStrCounts = hit->cigarcounts;
		cigarStrSize = hit->cigarsize;
		currentPair = i; // set this because the PrintReadInSAMFormat function will use it
		c = hit->strand;
		PrintReadInSAMFormat(((c=='+')?(fwdReads[i]):(revReads[i])),readNames[i],c,hit->pos,hit->numerrors,hit->score,bestReadsScores2nd[i],readsSizes[i]);
	}
	currentPair = 1; // set to 1, so the next read will be loaded from the correct reads file
}

void PrintSeeds(Read *read, char strand, int detailedmode){
	char *readChars;
	int numSeeds, seedId, seedPos, i;
	Seed *seedsArray, *seed;
	if(strand=='+'){
		strand='+';
		readChars=(read->fwdStrand.chars);
		numSeeds=(read->fwdStrand.numSeeds);
		seedsArray=(read->fwdStrand.seeds);
	} else {
		strand='-';
		readChars=(read->revStrand.chars);
		numSeeds=(read->revStrand.numSeeds);
		seedsArray=(read->revStrand.seeds);
	}
	if(numSeeds==0) return;
	if(detailedmode){ // details mode
		fprintf(debugfile,"%c\n%s\n",strand,readChars);
		for(i=0;i<numSeeds;i++){
			seed=&(seedsArray[i]);
			fprintf(debugfile,"%2d [%-3d-%3d] %3d %u \t",(i+1),(seed->startPosInRead),((seed->startPosInRead)+(seed->size)-1),(seed->size),(seed->numHits));
			if(i!=(numSeeds-1)) seedPos=((seedsArray[i+1].startPosInRead)+(seedsArray[i+1].size));
			else seedPos=0;
			while(seedPos<(seed->startPosInRead)) fprintf(debugfile,"%c",(char)(readChars[seedPos++]+32));
			fprintf(debugfile,"%.*s\n",(seed->size),(char *)(readChars+(seed->startPosInRead)));
		}
		if((seedsArray[(numSeeds-1)].startPosInRead)>1) fprintf(debugfile,"...\n");
	} else { // overview mode
		seedPos=0;
		for(seedId=(numSeeds-1);seedId>=0;seedId--){
			seed=&(seedsArray[seedId]);
			while(seedPos<(seed->startPosInRead)){
				fputc(' ',debugfile);
				seedPos++;
			}
			for(i=0;i<(seed->size);i++){
				if(i==0) fputc('<',debugfile);
				else if(i==((seed->size)-1)) fputc('|',debugfile);
				else fputc('=',debugfile);
				seedPos++;
			}
		}
		fputc('\n',debugfile);
	}
}

void GetSeeds(Read *read){
	unsigned int topPtr, bottomPtr;
	unsigned int *fwdTopPointers, *fwdBottomPointers;
	unsigned int *revTopPointers, *revBottomPointers;
	int seedSize, numFwdSeeds, numRevSeeds;
	int readSize, fwdPos, revPos;
	unsigned int numOccs;
	char *readChars, c;
	Seed *seed;
	readSize=(read->size);
	fwdTopPointers=(read->fwdStrand.topPointers);
	fwdBottomPointers=(read->fwdStrand.bottomPointers);
	revTopPointers=(read->revStrand.topPointers);
	revBottomPointers=(read->revStrand.bottomPointers);
	numFwdSeeds=0;
	numRevSeeds=0;
	fwdPos=(readSize-1); // current position in fwd read
	revPos=(readSize-1); // current position in rev read
	while(1){ // get one seed at a time for both fwd and rev strands
		readChars=(read->fwdStrand.chars); // process fwd chars
		topPtr=0; // initialize pointers
		bottomPtr=sizeBWT;
		seedSize=0;
		for(;fwdPos>=0;fwdPos--){ // get longest seed starting at current position of the read (from right to left)
			c=readChars[fwdPos]; // process all read chars until no occurrences of the substring exist
			numOccs=FMI_FollowLetter(c,&topPtr,&bottomPtr);
			if(numOccs==0){ // if there are no matches by this letter, stop seed here
				fwdTopPointers[fwdPos]=0; // set pointers of failed position to zero
				fwdBottomPointers[fwdPos]=0;
				break;
			}
			fwdTopPointers[fwdPos]=topPtr; // save index pointers for this position
			fwdBottomPointers[fwdPos]=bottomPtr;
			seedSize++;
		}
		if(seedSize!=0){ // save information about the newly created seed
			if(numFwdSeeds==maxNumSeeds){ // realloc seed arrays if needed
				maxNumSeeds=(numFwdSeeds+1);
				(read->fwdStrand.seeds)=(Seed *)realloc((read->fwdStrand.seeds),maxNumSeeds*sizeof(Seed));
				(read->revStrand.seeds)=(Seed *)realloc((read->revStrand.seeds),maxNumSeeds*sizeof(Seed));
			}
			seed=&(read->fwdStrand.seeds[numFwdSeeds]);
			(seed->size)=seedSize; // save seed size
			fwdPos++; // pos of the last valid char
			(seed->startPosInRead)=fwdPos; // save pos where this seed ends
			(seed->bwtTopPointer)=fwdTopPointers[fwdPos];
			(seed->numHits)=(fwdBottomPointers[fwdPos]-fwdTopPointers[fwdPos]+1);
			numFwdSeeds++; // one more seed created
			fwdPos-=2; // skip the erroneous position and start next seed on next position (to the left)
		}
		readChars=(read->revStrand.chars); // process rev chars
		topPtr=0; // initialize pointers
		bottomPtr=sizeBWT;
		seedSize=0;
		for(;revPos>=0;revPos--){ // get longest seed starting at current position of the read (from right to left)
			c=readChars[revPos]; // process all read chars until no occurrences of the substring exist
			numOccs=FMI_FollowLetter(c,&topPtr,&bottomPtr);
			if(numOccs==0){ // if there are no matches by this letter, stop seed here
				revTopPointers[revPos]=0; // set pointers of failed position to zero
				revBottomPointers[revPos]=0;
				break;
			}
			revTopPointers[revPos]=topPtr; // save index pointers for this position
			revBottomPointers[revPos]=bottomPtr;
			seedSize++;
		}
		if(seedSize!=0){ // save information about the newly created seed
			if(numRevSeeds==maxNumSeeds){ // realloc seed arrays if needed
				maxNumSeeds=(numRevSeeds+1);
				(read->fwdStrand.seeds)=(Seed *)realloc((read->fwdStrand.seeds),maxNumSeeds*sizeof(Seed));
				(read->revStrand.seeds)=(Seed *)realloc((read->revStrand.seeds),maxNumSeeds*sizeof(Seed));
			}
			seed=&(read->revStrand.seeds[numRevSeeds]);
			(seed->size)=seedSize; // save seed size
			revPos++; // pos of the last valid char
			(seed->startPosInRead)=revPos; // save pos where this seed ends
			(seed->bwtTopPointer)=revTopPointers[revPos];
			(seed->numHits)=(revBottomPointers[revPos]-revTopPointers[revPos]+1);
			numRevSeeds++; // one more seed created
			revPos-=2; // skip the erroneous position and start next seed on next position (to the left)
		}
		if(fwdPos<0 && revPos<0){ // stop if we have reached the end of the read in both strands
			(read->fwdStrand.numSeeds)=numFwdSeeds;
			(read->revStrand.numSeeds)=numRevSeeds;
			break;
		}
		if(numFwdSeeds!=numRevSeeds){ // stop if there will be a difference of more than 1 seed between both strands
			if(numFwdSeeds<numRevSeeds){
				(read->fwdStrand.numSeeds)=numFwdSeeds;
				(read->revStrand.numSeeds)=0; // discard strand seeds if they didn't get to the beggining of the read
			} else {
				(read->fwdStrand.numSeeds)=0;
				(read->revStrand.numSeeds)=numRevSeeds;
			}
			break;
		}
	}
	#ifdef DEBUG
	PrintSeeds(read,'+',1);
	PrintSeeds(read,'-',1);
	#endif
}

void GetImprovedSeeds(Read *read){
	unsigned int topPtr, bottomPtr;
	unsigned int *topPointers, *bottomPointers;
	int seedExtraSize, numSeeds, seedId;
	int seedPos, seedOldLeftPos, seedNewLeftPos, seedOldRightPos, seedNewRightPos, rightmostPos;
	unsigned int numOccs;
	char *readChars, strandcount, strandchar, c;
	Seed *seed;
	Strand *strand;
	for(strandcount=0;strandcount<2;strandcount++){
		if(strandcount==0){ // forward strand
			strand=&(read->fwdStrand);
			strandchar='+';
		} else { // reverse strand
			strand=&(read->revStrand);
			strandchar='-';
		}
		numSeeds=(strand->numSeeds);
		if(numSeeds==0) continue;
		#ifdef DEBUG
		// print old seeds
		fprintf(debugfile,"%c\n%s\n",strandchar,(strand->chars));
		PrintSeeds(read,strandchar,0);
		#endif
		readChars=(strand->chars);
		topPointers=(strand->topPointers);
		bottomPointers=(strand->bottomPointers);
		seedExtraSize=0;
		for(seedId=(numSeeds-1);seedId>=0;seedId--){ // process all seeds from the leftmost one to the (before the) rightmost one
			seed=&(strand->seeds[seedId]);
			seedOldLeftPos=(seed->startPosInRead); // position where the original seed ends (on the left)
			seedNewLeftPos=seedOldLeftPos; // new left pos of the seed after (possibly) being changed by the prev seed (at the left)
			if(seedExtraSize!=0){ // if the prev seed (to the left) was extended, it will change this seed too
				seedExtraSize--; // the mismatched char (that was jumped over) between the 2 seeds does not count
				seedNewLeftPos+=seedExtraSize; // update leftmost seed pos
				(seed->startPosInRead)=seedNewLeftPos; // the starting (left) pos of the seed was advanced to the right
				(seed->size)-=seedExtraSize; // and its size was shortened
				(seed->bwtTopPointer)=topPointers[seedNewLeftPos]; // set the correct pointers in case we never get to set them below
				(seed->numHits)=(bottomPointers[seedNewLeftPos]-topPointers[seedNewLeftPos]+1);
			}
			if(seedId==0) break; // the 1st (rightmost) seed cannot be further extended to the right
			seedNewRightPos=(seedNewLeftPos+(seed->size)); // start one position to the right of the current rightmost end of the seed
			seedOldRightPos=(seedNewRightPos-1); // last valid position on the left of the original seed
			rightmostPos=((strand->seeds[seedId-1].startPosInRead)+(strand->seeds[seedId-1].size)-1); // rightmost pos of the following (to the right) seed
			seedExtraSize=0; // extended seed size (only the extra length)
			for(;seedNewRightPos<rightmostPos;seedNewRightPos++){ // try starting the seed one char at a time to the right, but do not let the seed start beyond the right end of the next seed
				topPtr=0; // initialize pointers
				bottomPtr=sizeBWT;
				for(seedPos=seedNewRightPos;seedPos>=seedOldLeftPos;seedPos--){ // set new starting position of the seed and check how far to the left the seed can reach
					c=readChars[seedPos];
					numOccs=FMI_FollowLetter(c,&topPtr,&bottomPtr);
					if(numOccs==0){ // the seed ended
						seedPos++; // the last valid char was the one before (to the right)
						break;
					}
					if(seedPos<=seedOldRightPos){ // only compare and pointers if the pos is over the original seed pointers
						if((topPtr==topPointers[seedPos]) && (bottomPtr==bottomPointers[seedPos])){ // if the pointers are the same as previous found pointers
							seedPos=seedOldLeftPos; // the seed will end at a pos lower or equal at the original end pos, so we can continue extending it
							break;
						}
						topPointers[seedPos]=topPtr; // save pointers
						bottomPointers[seedPos]=bottomPtr;
					}
				}
				if(seedPos>seedOldLeftPos) break; // continue to extend the seed while the end is equal to (or -1) the original end
				seedExtraSize++; // seed was successfuly extended one char to the right
				(seed->bwtTopPointer)=topPointers[seedNewLeftPos]; // save pointers now because in the last (incomplete) run, they can have a wrong value
				(seed->numHits)=(bottomPointers[seedNewLeftPos]-topPointers[seedNewLeftPos]+1);
			} // next, start the seed one more char to the right
			(seed->size)+=seedExtraSize; // save the new size of this seed
		}
		#ifdef DEBUG
		PrintSeeds(read,strandchar,0);
		PrintSeeds(read,strandchar,1);
		#endif
	}
}

void CheckDepthSeeds(Read *read){
	unsigned int topPtr, bottomPtr;
	unsigned int *fwdTopPointers, *fwdBottomPointers;
	unsigned int *revTopPointers, *revBottomPointers;
	int *fwdSeedsEndPos, *revSeedsEndPos;
	//int numFwdSeeds, numRevSeeds;
	int readSize, i, j;
	unsigned int n;
	char *readChars, c;
	
	/*
	#ifdef DEBUG
	fprintf(debugfile,">%s\n%d\n",(read->name),(read->size));
	#endif
	*/

	readSize=(read->size);
	fwdTopPointers=(unsigned int *)malloc(readSize*sizeof(unsigned int));
	revTopPointers=(unsigned int *)malloc(readSize*sizeof(unsigned int));
	fwdBottomPointers=(unsigned int *)malloc(readSize*sizeof(unsigned int));
	revBottomPointers=(unsigned int *)malloc(readSize*sizeof(unsigned int));
	fwdSeedsEndPos=(int *)malloc(readSize*sizeof(int));
	revSeedsEndPos=(int *)malloc(readSize*sizeof(int));

	for(i=0;i<readSize;i++){ // reset index pointers
		fwdTopPointers[i]=0;
		revTopPointers[i]=0;
		fwdBottomPointers[i]=0;
		revBottomPointers[i]=0;
		fwdSeedsEndPos[i]=-1;
		revSeedsEndPos[i]=-1;
	}

	readChars=(read->fwdStrand.chars); // process fwd chars
	for(i=(readSize-1);i>=0;i--){ // get longest seed starting on each position of the read (from right to left)
		topPtr=0;
		bottomPtr=sizeBWT;
		for(j=i;j>=0;j--){ // process all read chars from position (i) to position 0
			c=readChars[j];
			n=FMI_FollowLetter(c,&topPtr,&bottomPtr);
			if(n==0) break;// if there are no matches by this letter, stop seed here
			if( topPtr==fwdTopPointers[j] && bottomPtr==fwdBottomPointers[j] ){ // if a prev start pos (i>j) already set this same index interval at this pos (j), the end pos will be the same
				j=fwdSeedsEndPos[(i+1)]; // set the same end pos
				j--; // decrease pos because it will be increased later
				break;
			}
			fwdTopPointers[j]=topPtr; // save index pointers for this position
			fwdBottomPointers[j]=bottomPtr;
		}
		j++; // last pos with matches
		fwdSeedsEndPos[i]=j; // save pos where this seed ended (last valid char)
		fwdTopPointers[i]=fwdTopPointers[j]; // save last valid index pointers for the seed starting at this pos
		fwdBottomPointers[i]=fwdBottomPointers[j];
	}

	#ifdef DEBUG
	fprintf(debugfile,"+\n%s\n",(read->fwdStrand.chars));
	for(i=(readSize-1);i>=0;i--){
		for(j=0;j<fwdSeedsEndPos[i];j++) fputc(' ',debugfile);
		if(j!=i){
			fputc('<',debugfile);
			j++;
		}
		for(;j<i;j++) fputc('=',debugfile);
		fputc('|',debugfile);
		j++;
		for(;j<readSize;j++) fputc(' ',debugfile);
		fprintf(debugfile,"[%d-%d] %d %u\n",fwdSeedsEndPos[i],i,(i-fwdSeedsEndPos[i]+1),(fwdBottomPointers[i]-fwdTopPointers[i]+1));
	}
	#endif

	readChars=(read->revStrand.chars); // process rev chars
	for(i=(readSize-1);i>=0;i--){ // get longest seed starting on each position of the read (from right to left)
		topPtr=0;
		bottomPtr=sizeBWT;
		for(j=i;j>=0;j--){ // process all read chars from position (i) to position 0
			c=readChars[j];
			n=FMI_FollowLetter(c,&topPtr,&bottomPtr);
			if(n==0) break;// if there are no matches by this letter, stop seed here
			if( topPtr==revTopPointers[j] && bottomPtr==revBottomPointers[j] ){ // if a prev start pos (i) already set this same index interval at this pos (j), the end pos will be the same
				j=revSeedsEndPos[(i+1)]; // set the same end pos
				j--; // decrease pos because it will be increased later
				break;
			}
			revTopPointers[j]=topPtr; // save index pointers for this position
			revBottomPointers[j]=bottomPtr;
		}
		j++; // last pos with matches
		revSeedsEndPos[i]=j; // save pos where this seed ended (last valid char)
		revTopPointers[i]=revTopPointers[j]; // save last valid index pointers for the seed starting at this pos
		revBottomPointers[i]=revBottomPointers[j];
	}

	#ifdef DEBUG
	fprintf(debugfile,"-\n%s\n",(read->revStrand.chars));
	for(i=(readSize-1);i>=0;i--){
		for(j=0;j<revSeedsEndPos[i];j++) fputc(' ',debugfile);
		if(j!=i){
			fputc('<',debugfile);
			j++;
		}
		for(;j<i;j++) fputc('=',debugfile);
		fputc('|',debugfile);
		j++;
		for(;j<readSize;j++) fputc(' ',debugfile);
		fprintf(debugfile,"[%d-%d] %d %u\n",revSeedsEndPos[i],i,(i-revSeedsEndPos[i]+1),(revBottomPointers[i]-revTopPointers[i]+1));
	}
	#endif

	/*
	if(maxNumSeeds==0){
		maxNumSeeds=1;
		(read->fwdSeeds)=(Seed *)malloc(maxNumSeeds*sizeof(Seed));
		(read->revSeeds)=(Seed *)malloc(maxNumSeeds*sizeof(Seed));
	}

	numFwdSeeds=0;
	i=(readSize+1);
	while(1){
		j=(i+2); // start searching 2 positions to the right of the end of the prev seed
		while(j<readSize && fwdSeedsEndPos[j]>=i) j--; // get the rightmost seed that ends to the left of the prev seed (at the right)
		if(j<0) break;
		i=fwdSeedsEndPos[j]; // end pos of the current seed
		j--;
		while(j>=0 && fwdSeedsEndPos[j]==i) j--; // get the rightmost seed that ends at this pos
		j++;
	}
	*/

	free(fwdTopPointers);
	free(revTopPointers);
	free(fwdBottomPointers);
	free(revBottomPointers);
	free(fwdSeedsEndPos);
	free(revSeedsEndPos);
}
/**/

// TODO: check missaligned reads
// TODO: (?) when the gap is larger than the maxnumerrors, it can still be a valid hit (if it is all matches from one side)
// TODO: (!) when numseeds>(maxnumerrors+1) don't get positions for seeds and get new seeds
// TODO: add Seed and SeedChain structs
void RunSeededAlignment(){
	int readSize; // size of the read
	char strand; // indicates if we are dealing with the forward ('+') or the reverse ('-') strand
	int maxReadErrors; // maximum number of errors allowed for the current read
	int tryHarder; // flag to indicate wether to use more aggressive read mapping method or not
	
	int i, j, d, m, n, s, prevs; // general purpose counters
	unsigned int p, k, prevk; // general purpose counters
	char c; // current character in seed
	int numSeeds; // number of seeds in the current read
	int *seedSizes; // array of the sizes of each seed
	int *numSeedHits; // array of the number of positions of each seed in the reference genome
	int numTotalSeedHits; // number of positions in the reference genome of all the seeds
	unsigned int *seedLeftPtrs; // array of the index left pointer of each seed
	int *seedStartPosInRead; // array of the position in the read where each seed starts

	int fwdNumSeeds;
	int revNumSeeds;
	int fwdNumTotalSeedHits;
	int revNumTotalSeedHits;
	int *fwdSeedSizes;
	int *revSeedSizes;
	int *fwdNumSeedHits;
	int *revNumSeedHits;
	unsigned int *fwdSeedLeftPtrs;
	unsigned int *revSeedLeftPtrs;
	int *fwdSeedStartPosInRead;
	int *revSeedStartPosInRead;

	SeedPos **fwdSeedPositions;
	SeedPos **revSeedPositions;
	int *fwdSeedChainSizes;
	int *revSeedChainSizes;
	int *fwdSeedChainScores;
	int *revSeedChainScores;
	int *fwdSeedChainStarts;
	int *revSeedChainStarts;

	SeedPos **seedPositions; // array of the positions in the reference genome of all the seeds
	int maxSeedArraysSize; // maximum number of positions that each seed or chain related array can currently store
	
	int numSeedChains; // number of groups of consecutive seeds with consistent positions
	int *seedChainSizes; // array of the number of seeds contained in each seed chain
	int *seedChainScores; // array of the sum of sizes of all the seeds in each seed chain
	int *seedChainStarts; // array of positions in the seedPositions array where the first seed of each seed chain is located
	int numBestSeedChains; // number of top scoring seed chains
	int *bestSeedChains; // array of identifiers of the top scoring seed chains
	char *bestSeedChainsStrands; // array of strands of the top scoring seed chains
	int bestSeedChainsScore; // score of the top scoring seed chains
	int distanceRelaxation; // allowed distance relaxation on end segments of chains
	
	double progress, timeleft; // counters for progress status
	struct timeb laststatus; // time of last status information

	char *charPtrAtRead; // pointer to a temporary position in the read
	char *charPtrAtGenome; // pointer to a temporary position in the genome
	
	int gapSizeInRead;
	int gapSizeInGenome;
	int gapStartPosInRead;
	unsigned int gapStartPosInGenome;
	int numSeedCharsToSkip; // seed size to skip in overlapping seeds when introducing seed size in cigar string
	int numIndelsAtEnd; // number of insertions and deletions in the tops of the read returned by DP when the extreme-most seed is missing

	unsigned int readPosInGenome; // position in the reference genome where the read aligns with less number of errors
	int numDPErrors; // sum of the errors of the dynamic programming alignment between the consistent seeds
	int numReadErrors; // number of errors of the alignment of the read position
	int bestReadNumErrors; // minimum number of errors among all the tested read positions
	int testOtherStrandToo; // variable that defines wether or not it will search for seed chains in the opposite strand

	ReadHit *readHits; // list of hits in the genome for the current read being processed
	int numReadHits; // number of hits for the current read
	int bestReadHitId; // index in the readHits array of the top scoring hit
	int numErrors2ndBestHit; // number of errors of the hit with the second best score
	
	int splicedRead; // indicates if a splicing event occurred in the read or not
	int checkGap;
	int numExtraSeedHits;
	unsigned int ssp;

	Read tempRead;
	Seed *seed, *seeds;

	maxNumSeeds = 1;
	tempRead.fwdStrand.seeds = (Seed *)malloc(maxNumSeeds*sizeof(Seed));
	tempRead.revStrand.seeds = (Seed *)malloc(maxNumSeeds*sizeof(Seed));

	maxReadSize = 1;
	tempRead.fwdStrand.topPointers = (unsigned int *)malloc(maxReadSize*sizeof(unsigned int));
	tempRead.revStrand.topPointers = (unsigned int *)malloc(maxReadSize*sizeof(unsigned int));
	tempRead.fwdStrand.bottomPointers = (unsigned int *)malloc(maxReadSize*sizeof(unsigned int));
	tempRead.revStrand.bottomPointers = (unsigned int *)malloc(maxReadSize*sizeof(unsigned int));

	// initialize these variables here to avoid compiler warnings
	numErrors2ndBestHit = 0;
	bestReadHitId = 0;
	bestReadNumErrors = 0;
	numReadErrors = 0;
	numBestSeedChains = 0;
	fwdNumTotalSeedHits = 0;
	fwdNumSeeds = 0;
	readHits = NULL;

	if( minIdentityPct > 0 ){ // if an identity percentage was defined, set the maximum number of errors now
		maxNumErrors = ( ( maxReadSize * ( 100 - minIdentityPct ) ) / 100 );
	} // else, the maxNumErrors variable as already been set with the command line arguments

	dpNumRows = ( 2*maxNumErrors + 1); // for banded DP, allow this number of errors above and below the matrix diagonal
	dpNumCols = ( maxReadSize + maxNumErrors + 1 ); // maximum segment size (plus errors if they were all indels at the ends) that we are going to align by DP (+1 for 0-th column)
	dpMatrix = (int **)malloc(dpNumRows*sizeof(int *));
	dpDirections = (char **)malloc(dpNumRows*sizeof(char *));
	for( i = 0 ; i < dpNumRows ; i++ ){
		dpMatrix[i] = (int *)malloc(dpNumCols*sizeof(int));
		dpDirections[i] = (char *)malloc(dpNumCols*sizeof(char));
	}
	#ifdef DEBUGDP
	dpTrace = (char *)malloc((dpNumCols)*sizeof(char));
	dpAlignedTarget = (char *)malloc((dpNumCols)*sizeof(char));
	dpAlignedRead = (char *)malloc((dpNumCols)*sizeof(char));
	#endif

	maxSeedArraysSize = 1;
	fwdSeedPositions = (SeedPos **)malloc(maxSeedArraysSize*sizeof(SeedPos *));
	for(i=0;i<maxSeedArraysSize;i++) fwdSeedPositions[i] = (SeedPos *)malloc(sizeof(SeedPos));
	fwdSeedChainSizes = (int *)malloc(maxSeedArraysSize*sizeof(int));
	fwdSeedChainScores = (int *)malloc(maxSeedArraysSize*sizeof(int));
	fwdSeedChainStarts = (int *)malloc(maxSeedArraysSize*sizeof(int));
	revSeedPositions = (SeedPos **)malloc(maxSeedArraysSize*sizeof(SeedPos *));
	for(i=0;i<maxSeedArraysSize;i++) revSeedPositions[i] = (SeedPos *)malloc(sizeof(SeedPos));
	revSeedChainSizes = (int *)malloc(maxSeedArraysSize*sizeof(int));
	revSeedChainScores = (int *)malloc(maxSeedArraysSize*sizeof(int));
	revSeedChainStarts = (int *)malloc(maxSeedArraysSize*sizeof(int));

	bestSeedChains = (int *)malloc(maxSeedArraysSize*sizeof(int));
	bestSeedChainsStrands = (char *)malloc(maxSeedArraysSize*sizeof(char));

	// TODO: clean arrays at end
	fwdSeedSizes = (int *)malloc(maxReadSize*sizeof(int));
	fwdNumSeedHits = (int *)malloc(maxReadSize*sizeof(int));
	fwdSeedStartPosInRead = (int *)malloc(maxReadSize*sizeof(int));
	fwdSeedLeftPtrs = (unsigned int *)malloc(maxReadSize*sizeof(unsigned int));
	revSeedSizes = (int *)malloc(maxReadSize*sizeof(int));
	revNumSeedHits = (int *)malloc(maxReadSize*sizeof(int));
	revSeedStartPosInRead = (int *)malloc(maxReadSize*sizeof(int));
	revSeedLeftPtrs = (unsigned int *)malloc(maxReadSize*sizeof(unsigned int));

	numReads = 0;
	numAlignedReads = 0;
	numUnalignedReads = 0;
	numAlignedHardReads = 0;
	numHighRepetitionReads = 0;
	totalNumSeeds = 0;
	totalNumSeedOccs = 0;
	totalNumSeedChains = 0;
	totalNumDPErrors = 0;
	totalNumReadsChars = 0;
	totalNumMappedReadsChars = 0;
	totalNumReadHits = 0;
	maxReadHits = 0;
	
	if( rnaSequencingMode ){
		numSplicingEvents = 0;
		numSplicedReads = 0;
		numSeedsAddedBySpliceSites = 0;
		InitializeSpliceSitesArray(sizeRefGenome);
	}

	if( pairedEndReadsMode ){
		numAlignedBothMates = 0;
		numAlignedSingleMates = 0;
		minObservedPairDistance = sizeRefGenome;
		maxObservedPairDistance = 0;
		avgObservedPairDistances = 0;
		weightObservedPairDistances = 0;
		numPairs = 2; // if at least one of the reads file is of paired-ends type, set this variable to correctly alloc all arrays
	} else {
		numPairs = 1;
	}

	maxReadArraySize = 1;
	maxReadNameArraySize = 255;
	for( i = 0 ; i < numPairs ; i++ ){ // initialize read chars and read name arrays
		fwdReads[i] = (char *)malloc((maxReadArraySize+1)*sizeof(char)); // forward strand read
		revReads[i] = (char *)malloc((maxReadArraySize+1)*sizeof(char)); // reverse strand read
		readNames[i] = (char *)malloc((maxReadNameArraySize+1)*sizeof(char)); // description of the read
	}

	maxReadHitsArraySize = 1; // initialize read hits lists
	for(i=0;i<numPairs;i++){ // no need to malloc (ReadHit **)(readsHits) because it was initialized as (ReadHit *)(readsHits[2])
		readsHits[i]=(ReadHit *)malloc(maxReadHitsArraySize*sizeof(ReadHit));
		readHits=readsHits[i];
		for(j=0;j<maxReadHitsArraySize;j++){
			readHits[j].cigarcodes = (char *)malloc((maxReadSize+maxNumErrors)*sizeof(char));
			readHits[j].cigarcounts = (int *)malloc((maxReadSize+maxNumErrors)*sizeof(int));
		}
	}

	if( !pairedEndReadsMode ){
		fwdReads[1] = NULL;
		revReads[1] = NULL;
		readNames[1] = NULL;
		readsHits[1] = NULL;
	}

	PrintSAMHeader();

	ftime(&laststatus);
	ftime(&starttb); // start time counter
	currentReadsFile=-1; // it will be incremented to 0 next
	goto _next_file;

	while(1){ // loop for all the reads
		if( pairedEndReadsMode ){ // if we are in paired-end reads mode, switch to the other reads file
			//currentPair = ( (currentPair+1) & 1); // next index in array (mod 2)
			if(currentPair==1){ // pair was 1, next will be 0
				currentPair=0;
				LoadNextRead(); // load a new set of 2 reads
			} else currentPair=1; // pair was 0, now is 1
			readsFile = pairFiles[currentPair]; // switch to the variables of the other mate read
			readHits = readsHits[currentPair];
			fwdRead = fwdReads[currentPair];
			revRead = revReads[currentPair];
			readName = readNames[currentPair];
			read = fwdRead;
			readSize = readsSizes[currentPair];
		} else { // single reads mode
			read = fwdRead; // load read in forward read array
			readSize = LoadNextRead(); // load next read from the reads file
		}
		if( readSize == 0 ){ // no more reads in this reads file, go to next one
			_next_file:
			currentReadsFile++;
			if(currentReadsFile==numReadsFiles) break; // last reads file
			if(currentReadsFile>0) for(i=0;i<numPairs;i++) if(pairFiles[i]!=NULL) prevFilesData += (long long int)ftell(pairFiles[i]); // store sizes of last processed files (not NULL because 2nd file in P.E. SFF does not exist)
			readsFile=readsFiles[currentReadsFile];
			readsFilename=readsFilenames[currentReadsFile];
			c=readsFileTypes[currentReadsFile];
			if(c=='f' || c=='q' || c=='s'){ // lowercase filetype means paired-end reads file
				pairedEndReadsMode=1;
				maxReadPairDistance=maxDists[currentReadsFile];
				minReadPairDistance=minDists[currentReadsFile];
				if(c=='f' || c=='q'){ // two files, one for each mate reads
					pairFiles[0]=readsFiles[(currentReadsFile)];
					pairFiles[1]=readsFiles[(currentReadsFile+1)];
					rewind(pairFiles[0]);
					rewind(pairFiles[1]);
					if(c=='f') LoadNextRead=LoadNextReadPairFromFasta;
					else if(c=='q') LoadNextRead=LoadNextReadPairFromFastaQ;
					currentReadsFile++; // advance to 2nd pair, so next time it will load a new file
				} else { // one single SFF file for both mate reads
					pairFiles[0]=readsFiles[currentReadsFile];
					pairFiles[1]=NULL;
					rewind(pairFiles[0]);
					AnalyzeSffReads(0,0); // load in memory the correct SFF variables for this SFF file
					LoadNextRead=LoadNextReadPairFromSff;
				}
				numPairs=2;
				currentPair=1; // this will change to 0 the next time we enter the loop
			} else { // uppercase filetype means single reads file
				pairedEndReadsMode=0;
				rewind(readsFile);
				if(c=='F'){ // set correct read loading functions according to file type
					LoadNextRead=LoadNextReadFromFasta;
				} else if(c=='Q'){
					LoadNextRead=LoadNextReadFromFastaQ;
				} else if(c=='S'){
					AnalyzeSffReads(0,0); // load in memory the correct SFF variables for this SFF file
					LoadNextRead=LoadNextReadFromSff;
				}
				numPairs=1;
				currentPair=0;
				pairFiles[0]=readsFile;
				pairFiles[1]=NULL;
				fwdRead=fwdReads[0]; // set pointers to the arrays of the first/only reads file
				revRead=revReads[0];
				readName=readNames[0];
				readHits=readsHits[0];
			}
			continue; // back to beginning of loop to start loading reads from new file
		}
		if( ((numReads ^ 0x0400) & 0x07FF) == 0 ){ // only get time each 1024 reads
			if( GetElapsedTime(laststatus) > STATSINTERVAL ) { // print progress from time to time
				currentFileData = prevFilesData; // progress of previous already completed files
				for(i=0;i<numPairs;i++) if(pairFiles[i]!=NULL) currentFileData += (long long int)ftell(pairFiles[i]); // combined file pos of both pairs (or only one) (not NULL because 2nd file in P.E. SFF does not exist)
				progress = (( (double) currentFileData ) / ( (double) totalFilesData )) * 100.0;
				timeleft = ( ( (double) ( totalFilesData - currentFileData ) * GetElapsedTime(starttb) ) / ( (double) currentFileData ) );
				#ifdef __unix__
				PrintProgressBar(progress,1);
				#endif
				printf("\n> Aligning reads... %3.2lf%% (%d processed: %d aligned + %d unaligned) ",progress,numReads,numAlignedReads,numUnalignedReads);
				PrintTime(timeleft);
				printf("left");
				fflush(stdout);
				ftime(&laststatus);
			}
		}
		numReads++;
		totalNumReadsChars += (unsigned int)readSize;
		if( readSize > maxReadSize ){ // realloc arrays to support the larger read size, if needed
			maxReadSize = readSize;
			n = dpNumRows; // old value for dpNumRows
			if( minIdentityPct > 0 ){ // if an identity percentage was defined, the maxNumErrors variable changes too
				maxNumErrors = ( ( maxReadSize * ( 100 - minIdentityPct ) ) / 100 );
				dpNumRows = ( 2*maxNumErrors + 1);
				dpMatrix = (int **)realloc(dpMatrix,dpNumRows*sizeof(int *));
				dpDirections = (char **)realloc(dpDirections,dpNumRows*sizeof(char *));
			}
			dpNumCols = ( maxReadSize + maxNumErrors + 1 );
			for( i = 0 ; i < n ; i++ ){ // realloc for old existing rows
				dpMatrix[i] = (int *)realloc(dpMatrix[i],dpNumCols*sizeof(int));
				dpDirections[i] = (char *)realloc(dpDirections[i],dpNumCols*sizeof(char));
			}
			for( i = n ; i < dpNumRows ; i++ ){ // simple initial alloc for new rows only
				dpMatrix[i] = (int *)calloc(dpNumCols,sizeof(int));
				dpDirections[i] = (char *)calloc(dpNumCols,sizeof(char));
			}
			#ifdef DEBUGDP
			dpTrace = (char *)realloc(dpTrace,(dpNumCols)*sizeof(char));
			dpAlignedTarget = (char *)realloc(dpAlignedTarget,(dpNumCols)*sizeof(char));
			dpAlignedRead = (char *)realloc(dpAlignedRead,(dpNumCols)*sizeof(char));
			#endif
			for(i=0;i<numPairs;i++){
				for(j=0;j<maxReadHitsArraySize;j++){
					readsHits[i][j].cigarcodes = (char *)realloc((readsHits[i][j].cigarcodes),(maxReadSize+maxNumErrors)*sizeof(char));
					readsHits[i][j].cigarcounts = (int *)realloc((readsHits[i][j].cigarcounts),(maxReadSize+maxNumErrors)*sizeof(int));
				}
			}
			fwdSeedSizes = (int *)realloc(fwdSeedSizes,maxReadSize*sizeof(int));
			fwdNumSeedHits = (int *)realloc(fwdNumSeedHits,maxReadSize*sizeof(int));
			fwdSeedStartPosInRead = (int *)realloc(fwdSeedStartPosInRead,maxReadSize*sizeof(int));
			fwdSeedLeftPtrs = (unsigned int *)realloc(fwdSeedLeftPtrs,maxReadSize*sizeof(unsigned int));
			revSeedSizes = (int *)realloc(revSeedSizes,maxReadSize*sizeof(int));
			revNumSeedHits = (int *)realloc(revNumSeedHits,maxReadSize*sizeof(int));
			revSeedStartPosInRead = (int *)realloc(revSeedStartPosInRead,maxReadSize*sizeof(int));
			revSeedLeftPtrs = (unsigned int *)realloc(revSeedLeftPtrs,maxReadSize*sizeof(unsigned int));

			tempRead.fwdStrand.topPointers=(unsigned int *)realloc((tempRead.fwdStrand.topPointers),maxReadSize*sizeof(unsigned int));
			tempRead.revStrand.topPointers=(unsigned int *)realloc((tempRead.revStrand.topPointers),maxReadSize*sizeof(unsigned int));
			tempRead.fwdStrand.bottomPointers=(unsigned int *)realloc((tempRead.fwdStrand.bottomPointers),maxReadSize*sizeof(unsigned int));
			tempRead.revStrand.bottomPointers=(unsigned int *)realloc((tempRead.revStrand.bottomPointers),maxReadSize*sizeof(unsigned int));

		} // end of arrays realloc for new max read size
		
		
		if( minIdentityPct > 0 ) maxReadErrors = ( ( readSize * ( 100 - minIdentityPct ) ) / 100 ); // define number of errors based on the non-identity percentage
		else maxReadErrors = maxNumErrors;
		#ifdef DEBUG
		fprintf(debugfile,">%s\n%u %d %d\n",readName,numReads,readSize,maxReadErrors);
		#endif
		splicedRead=0;

		LoadReverseStrandRead(fwdRead,revRead,readSize);
		tempRead.name=readName;
		tempRead.size=readSize;
		tempRead.fwdStrand.chars=fwdRead;
		tempRead.revStrand.chars=revRead;
		GetSeeds(&tempRead);
		
		/*
		CheckDepthSeeds(&tempRead);
		continue;
		*/

		fwdNumSeeds = (tempRead.fwdStrand.numSeeds);
		seeds = (tempRead.fwdStrand.seeds);
		fwdNumTotalSeedHits = 0;
		for( s = 0 ; s < fwdNumSeeds ; s++ ){
			seed = &(seeds[s]);
			fwdSeedSizes[s] = (seed->size);
			fwdNumSeedHits[s] = (seed->numHits);
			fwdSeedStartPosInRead[s] = (seed->startPosInRead);
			fwdSeedLeftPtrs[s] = (seed->bwtTopPointer);
			if( fwdNumSeedHits[s] <= maxSeedHits ) fwdNumTotalSeedHits += fwdNumSeedHits[s];
		}
		revNumSeeds = (tempRead.revStrand.numSeeds);
		seeds = (tempRead.revStrand.seeds);
		revNumTotalSeedHits = 0;
		for( s = 0 ; s < revNumSeeds ; s++ ){
			seed = &(seeds[s]);
			revSeedSizes[s] = (seed->size);
			revNumSeedHits[s] = (seed->numHits);
			revSeedStartPosInRead[s] = (seed->startPosInRead);
			revSeedLeftPtrs[s] = (seed->bwtTopPointer);
			if( revNumSeedHits[s] <= maxSeedHits ) revNumTotalSeedHits += revNumSeedHits[s];
		}

		if( ( fwdNumTotalSeedHits == 0 ) && ( revNumTotalSeedHits == 0 ) ){ // if there are too many hits for both strands, there is nothing to do
			numHighRepetitionReads++;
			numReadHits = 0;
			goto _unmapped_read; // report that read was not aligned and go to next read
		}

		if( ( fwdNumSeeds == 1 ) || ( revNumSeeds == 1 ) ){ // if the read occurs in the genome with no errors, skip the seeds chains part and output it directly
			numReadHits = 0;
			numTotalSeedHits = 0; // to count the total number of hits between the two strands
			strand = '+';
			read = fwdRead;
			numSeeds = fwdNumSeeds;
			numSeedHits = fwdNumSeedHits;
			seedLeftPtrs = fwdSeedLeftPtrs;
			seedStartPosInRead = fwdSeedStartPosInRead;
			while(1){ // iterate one time for the fwd strand and one time for the rev strand
				if( ( numSeeds == 1 ) && ( numSeedHits[0] <= maxSeedHits ) ){ // check if the fwd strand has a single seed but not too many hits
					numTotalSeedHits += numSeedHits[0];
					if( numTotalSeedHits > maxReadHitsArraySize ){ // resize the read hit lists if needed
						IncreaseMaxReadHitsArraySize(numTotalSeedHits,(maxReadSize+maxNumErrors));
						readHits=readsHits[currentPair];
					}
					numReadErrors = 0;
					// TODO: if the read has a mismatch on the 1st pos, instead of changing all these variables, set (seedPositions[*]->readpos=1) (?)
					if( seedStartPosInRead[0] == 1 ){ // mismatch in first char of the read
						numReadErrors++;
						if(maxReadErrors==0) numSeedHits[0]=0; // if no errors all allowed, reads with 1 mismatch at the end are not valid, and we will not enter the loop bellow
					}
					k = seedLeftPtrs[0];
					for( i = 0 ; i < numSeedHits[0] ; i++ ){ // iterate through all the hits
						p = FMI_PositionInText(k);
						if( seedStartPosInRead[0] == 1 && p!=0 ) p--; // if there was a mismatch in the first char of the read, shift the ref position accordingly
						readHits[numReadHits].pos = p;
						readHits[numReadHits].strand = strand;
						readHits[numReadHits].numerrors = numReadErrors;
						readHits[numReadHits].score = ( readSize - numReadErrors );
						readHits[numReadHits].cigarsize = 0; // it is not the size, but the last filled position
						readHits[numReadHits].cigarcodes[0] = 'M';
						readHits[numReadHits].cigarcounts[0] = readSize;
						cigarStrCodes = readHits[numReadHits].cigarcodes;
						cigarStrCounts = readHits[numReadHits].cigarcounts;
						cigarStrSize = readHits[numReadHits].cigarsize;
						numReadHits++;
						k++;
					}
				}
				if( strand == '-' ) break;
				strand = '-';
				read = revRead;
				numSeeds = revNumSeeds;
				numSeedHits = revNumSeedHits;
				seedLeftPtrs = revSeedLeftPtrs;
				seedStartPosInRead = revSeedStartPosInRead;
			} // end of loop for both strands
			numSeeds = 1; // set these vars to update stats after the jump
			numBestSeedChains = numTotalSeedHits;
			bestReadHitId = 0;
			bestReadNumErrors = numReadErrors;
			numErrors2ndBestHit = readSize;
			if( numReadHits == 0 ) goto _unmapped_read; // can happen if we want 0 errors but both strands have a mismatch at the end
			if( numReadHits > 1 ) numErrors2ndBestHit = numReadErrors;
			#ifdef DEBUG
			fprintf(debugfile, "%c%d %d\n", readHits[0].strand, readHits[0].pos, readHits[0].numerrors);
			#endif
			goto _mapped_read; // report aligned reads and continue to next read
		}

		testOtherStrandToo = 0; // test only one strand (for now)
		if( fwdNumSeeds <= revNumSeeds ){ // select the strand with lower number of seeds (has larger seeds which means less errors)
			strand = '+';
			if( ( revNumSeeds - fwdNumSeeds ) < 2 ) testOtherStrandToo = 1; // if the difference is just 0 or 1, test the other strand too later
		} else {
			strand = '-';
			if( ( fwdNumSeeds - revNumSeeds ) < 2 ) testOtherStrandToo = 1;
		}
		
		tryHarder = 0;
		while(1){ // iterate one time with normal mapping method and, if needed, another time using more agressive mapping method

			// TODO: when get improved seeds, for each seed hit, extend it to the left too by exact matching
			if( tryHarder ){ // if we need to try harder

				GetImprovedSeeds(&tempRead);
				fwdNumSeeds = (tempRead.fwdStrand.numSeeds);
				seeds = (tempRead.fwdStrand.seeds);
				fwdNumTotalSeedHits = 0;
				for( s = 0 ; s < fwdNumSeeds ; s++ ){
					seed = &(seeds[s]);
					fwdSeedSizes[s] = (seed->size);
					fwdNumSeedHits[s] = (seed->numHits);
					fwdSeedStartPosInRead[s] = (seed->startPosInRead);
					fwdSeedLeftPtrs[s] = (seed->bwtTopPointer);
					if( fwdNumSeedHits[s] <= maxSeedHits ) fwdNumTotalSeedHits += fwdNumSeedHits[s];
				}
				revNumSeeds = (tempRead.revStrand.numSeeds);
				seeds = (tempRead.revStrand.seeds);
				revNumTotalSeedHits = 0;
				for( s = 0 ; s < revNumSeeds ; s++ ){
					seed = &(seeds[s]);
					revSeedSizes[s] = (seed->size);
					revNumSeedHits[s] = (seed->numHits);
					revSeedStartPosInRead[s] = (seed->startPosInRead);
					revSeedLeftPtrs[s] = (seed->bwtTopPointer);
					if( revNumSeedHits[s] <= maxSeedHits ) revNumTotalSeedHits += revNumSeedHits[s];
				}

				if( ( fwdNumTotalSeedHits == 0 ) && ( revNumTotalSeedHits == 0 ) ){ // if there are too many hits for both strands, there is nothing to do
					numHighRepetitionReads++;
					numReadHits = 0;
					goto _unmapped_read; // report that read was not aligned and go to next read
				}

			}

			if( ( fwdNumTotalSeedHits == 0 ) || ( revNumTotalSeedHits == 0 ) ){ // check if there are more than maxSeedHits occurrences of all seeds
				testOtherStrandToo = 0; // do not test the other strand
				if( fwdNumTotalSeedHits != 0 ) strand = '+';
				else strand = '-';
			}
		
			// TODO: for *SeedPositions arrays, use array of structs (SeedPos*) instead of array of pointers to structs (SeedPos**)
			n = ( fwdNumTotalSeedHits + revNumTotalSeedHits );
			if( n > maxSeedArraysSize ){ // if there could be more chains or positions than the current arrays can store, then reallocate the arrays
				fwdSeedPositions = (SeedPos **)realloc(fwdSeedPositions,n*sizeof(SeedPos *));
				for( i = maxSeedArraysSize ; i < n ; i++ ){
					fwdSeedPositions[i] = (SeedPos *)malloc(sizeof(SeedPos));
				}
				revSeedPositions = (SeedPos **)realloc(revSeedPositions,n*sizeof(SeedPos *));
				for( i = maxSeedArraysSize ; i < n ; i++ ){
					revSeedPositions[i] = (SeedPos *)malloc(sizeof(SeedPos));
				}
				maxSeedArraysSize = n;
				fwdSeedChainSizes = (int *)realloc(fwdSeedChainSizes,maxSeedArraysSize*sizeof(int));
				fwdSeedChainScores = (int *)realloc(fwdSeedChainScores,maxSeedArraysSize*sizeof(int));
				fwdSeedChainStarts = (int *)realloc(fwdSeedChainStarts,maxSeedArraysSize*sizeof(int));
				revSeedChainSizes = (int *)realloc(revSeedChainSizes,maxSeedArraysSize*sizeof(int));
				revSeedChainScores = (int *)realloc(revSeedChainScores,maxSeedArraysSize*sizeof(int));
				revSeedChainStarts = (int *)realloc(revSeedChainStarts,maxSeedArraysSize*sizeof(int));
				bestSeedChains = (int *)realloc(bestSeedChains,maxSeedArraysSize*sizeof(int));
				bestSeedChainsStrands = (char *)realloc(bestSeedChainsStrands,maxSeedArraysSize*sizeof(char));
			}

			// TODO: when (testOtherStrandToo==1), test if it is worth keeping separate best chain scores for each strand or one global score for both strands
			numBestSeedChains = 0; // keep best seed chains for both strands (if not reseted to zero inside loop)
			//bestSeedChainsScore = 0; // if initialized here (and not in beggining of inner loop), will only keep best chains among the two strands
			while(1){ // process one or both strands and find the best seed chains
				if( strand == '+' ){ // load variables for forward strand
					read = fwdRead;
					numSeeds = fwdNumSeeds;
					numTotalSeedHits = fwdNumTotalSeedHits;
					seedSizes = fwdSeedSizes;
					numSeedHits = fwdNumSeedHits;
					seedLeftPtrs = fwdSeedLeftPtrs;
					seedStartPosInRead = fwdSeedStartPosInRead;
					seedPositions = fwdSeedPositions;
					seedChainSizes = fwdSeedChainSizes;
					seedChainScores = fwdSeedChainScores;
					seedChainStarts = fwdSeedChainStarts;
				} else { // load variables for reverse strand
					read = revRead;
					numSeeds = revNumSeeds;
					numTotalSeedHits = revNumTotalSeedHits;
					seedSizes = revSeedSizes;
					numSeedHits = revNumSeedHits;
					seedLeftPtrs = revSeedLeftPtrs;
					seedStartPosInRead = revSeedStartPosInRead;
					seedPositions = revSeedPositions;
					seedChainSizes = revSeedChainSizes;
					seedChainScores = revSeedChainScores;
					seedChainStarts = revSeedChainStarts;
				}
				#ifdef DEBUG
				fprintf(debugfile,"%c\n",strand);
				#endif
		
				n = 0; // current position in the seedPositions array
				for( s = 0 ; s < numSeeds ; s++ ){ // fill the positions array with all the positions of all the seeds
					if( numSeedHits[s] > maxSeedHits ) continue; // skip seed if it has too many hits
					k = seedLeftPtrs[s];
					for( i = 0 ; i < numSeedHits[s] ; i++ ){
						seedPositions[n]->seed = s;
						seedPositions[n]->size = seedSizes[s];
						seedPositions[n]->pos = FMI_PositionInText(k);
						seedPositions[n]->readpos = seedStartPosInRead[s]; // needed because individual seed ocorrences can shift their position inside the read when in RNA-Seq mode
						seedPositions[n]->ssId = 0;
						n++;
						k++;
					}
				}
				qsort(seedPositions,numTotalSeedHits,sizeof(SeedPos *),CompareSeedPositions); // sort positions

				/*
				#ifdef DEBUG
				for( n = 0 ; n < numTotalSeedHits ; n++ ){
				s = seedPositions[n]->seed;
				p = seedPositions[n]->pos;
				m = seedSizes[s];
				charPtrAtRead = (char *)(read + seedStartPosInRead[s]);
				charPtrAtGenome = (char *)(refGenome + p);
				fprintf(debugfile,"(%d,%u,'%.*s','%.*s',%d)\n",s,p,m,charPtrAtRead,m,charPtrAtGenome,m);
				}
				#endif
				*/

				// TODO: check if this makes sense to do now in the new hard mode or not
				// extend the seeds to the right by exact matching (but not beyond the start of the seed at the right)
				if (!tryHarder) { // disable this loop in hard mode because the seeds are already extended to the right
					for (n = 0; n < numTotalSeedHits; n++) { // extend all seeds to the right (because they were only extended maximally to the left)
						s = seedPositions[n]->seed;
						if (s == 0) continue; // if it is not the right-most seed, try to extend it as much as possible to the right by exact matching (the previous index matching only extended it to the left)
						if (n != (numTotalSeedHits - 1)) k = seedPositions[(n + 1)]->pos; // start pos of the next sorted seed to the right in the genome
						else k = sizeRefGenome;
						p = (seedPositions[n]->pos + seedSizes[s]);
						charPtrAtRead = (char *)(read + seedStartPosInRead[s] + seedSizes[s]); // pointer to the char next to the right-most char of the seed
						charPtrAtGenome = (char *)(refGenome + p);
						m = 0; // number of additionally matched characters of the seed
						while ((p < k) && (charPtrAtRead[0] == charPtrAtGenome[0])) { // extend while chars match, but not beyond the start of the next seed
							charPtrAtRead++; // advance to next char in read
							charPtrAtGenome++; // advance to next char in genome
							m++;
							p++;
						}
						seedPositions[n]->size += m; // update seed size
					}
				}

				#ifndef DISABLE_SS_ARRAY
				/**********************/
				/* Begin RNA-Seq Code */
				/**********************/
				// TODO: merge this part with the seed chaining part, and add new seeds dynamically by linked lists (and allowing to discover multiple new seeds in a row)
				// TODO: when looking for seeds chains later on, instead of only checking the seeds immediatelly before and after, check for an interval consistent with the maximum allowed intron length (non-linear seed chainning)
				// TODO: also check the GetImprovedSeeds code to see if anything can be improved
				if (rnaSequencingMode) {
					numExtraSeedHits = 0; // number of extra seeds added for this read
					prevs = numSeeds;
					prevk = 0;
					gapSizeInGenome = 0;
					for (i = 0; i < numTotalSeedHits; i++) {
						if (numReads == 39877 && i == 3) {
							numReads = 39877;
						}
						s = seedPositions[i]->seed;
						k = seedPositions[i]->pos;
						checkGap = (s != (numSeeds - 1)) && ( // is not the leftmost/last seed id, and
							(i == 0) // is the first seed hit
							|| (s != (prevs-1)) // or the previous seed hit does not have the previous seed id
							// TODO: the gap size here might be consistent but still exists a better seed in the middle (but check first if it has the same position as the existing seed or not)
							|| (gapSizeInGenome > maxIntronLength) // or the gap length with the other seed does not already correspond to an intron
							// TODO: when non-linear chains are available, re-add the check "||(gapSizeInGenome<minIntronLength)" here too;
							// now it's not useful because the new extra seed would be added before the previous seed and would not be checked/chained (and would cause errors too!)
							);
						if (checkGap) { // check if there is an existing splice site to the left
							ssp = k; // leftmost position of the seed hit
							p = GetOppositePositionIfSpliceSiteExists(&ssp); // get real splice site position (ssp) and opposite position (p)
							if ((p != 0) && (p < ssp)) { // if a splice site does exist to the left
								d = (int)(ssp - k); // difference between seed starting position (possibly to the left) and splice site position (possibly to the right)
								n = (seedPositions[i]->readpos) + d; // account for that difference in the seed's starting position in the read
								n--; // next position on the left to match in the read
								if (n >= 0) {
									m = 0; // number of matched chars
									charPtrAtRead = (char *)(read + n); // pointer to the char at the left of the beginning of the seed (corrected for the difference)
									charPtrAtGenome = (char *)(refGenome + p); // pointer to the char at the rightmost end of the left exon
									while ((*charPtrAtRead) == (*charPtrAtGenome)) { // while the characters in the read and in the genome match,
										m++;
										n--;
										p--;
										if ((charPtrAtRead == read) || (charPtrAtGenome == refGenome)) break; // ensure that the new seed does not overlap the existing one on the left
										charPtrAtRead--; // keep advancing both pointer to the left
										charPtrAtGenome--;
									}
									if (m > 0) { // if some chars were matched, create a new seed to the left
										seedPositions[i]->pos += d; // but first change the current seed to begin exactly at the splice site's (right) position
										seedPositions[i]->readpos += d;
										seedPositions[i]->size -= d;
										numSeedsAddedBySpliceSites++;
										numExtraSeedHits++;
										d = (numTotalSeedHits + numExtraSeedHits); // new seed hit id (+1)
										if (d > maxSeedArraysSize) { // reallocate the arrays if needed
											fwdSeedPositions = (SeedPos **)realloc(fwdSeedPositions, d * sizeof(SeedPos *));
											for (j = maxSeedArraysSize; j < d; j++) {
												fwdSeedPositions[j] = (SeedPos *)malloc(sizeof(SeedPos));
											}
											revSeedPositions = (SeedPos **)realloc(revSeedPositions, d * sizeof(SeedPos *));
											for (j = maxSeedArraysSize; j < d; j++) {
												revSeedPositions[j] = (SeedPos *)malloc(sizeof(SeedPos));
											}
											maxSeedArraysSize = d;
											fwdSeedChainSizes = (int *)realloc(fwdSeedChainSizes, maxSeedArraysSize * sizeof(int));
											fwdSeedChainScores = (int *)realloc(fwdSeedChainScores, maxSeedArraysSize * sizeof(int));
											fwdSeedChainStarts = (int *)realloc(fwdSeedChainStarts, maxSeedArraysSize * sizeof(int));
											revSeedChainSizes = (int *)realloc(revSeedChainSizes, maxSeedArraysSize * sizeof(int));
											revSeedChainScores = (int *)realloc(revSeedChainScores, maxSeedArraysSize * sizeof(int));
											revSeedChainStarts = (int *)realloc(revSeedChainStarts, maxSeedArraysSize * sizeof(int));
											bestSeedChains = (int *)realloc(bestSeedChains, maxSeedArraysSize * sizeof(int));
											bestSeedChainsStrands = (char *)realloc(bestSeedChainsStrands, maxSeedArraysSize * sizeof(char));
											if (strand == '+') { // reload modified array variables for the current strand
												seedPositions = fwdSeedPositions;
												seedChainSizes = fwdSeedChainSizes;
												seedChainScores = fwdSeedChainScores;
												seedChainStarts = fwdSeedChainStarts;
											}
											else {
												seedPositions = revSeedPositions;
												seedChainSizes = revSeedChainSizes;
												seedChainScores = revSeedChainScores;
												seedChainStarts = revSeedChainStarts;
											}
										} // realloc arrays
										d--; // d-th seed goes in the (d-1)-th position
										seedPositions[d]->seed = (s + 1); // seed id of the seed to the left
										seedPositions[d]->size = m;
										seedPositions[d]->pos = (p + 1); // +1 because the last position was a mismatch
										seedPositions[d]->readpos = (n + 1);
										// TODO: correctly fill this with splice site id to indicate that the seed was created from that splice site, and don't have to go get the id later on
										seedPositions[d]->ssId = 1;
									} // create new seed
								}
							} // splice site check
						} // left gap check
						if (i != (numTotalSeedHits - 1)) { // get size of gap between this seed and the next one
							prevk = (seedPositions[i + 1]->pos); // starting position of the existing seed to the right (here it should be not "prev" but "next")
							gapSizeInGenome = prevk - ((seedPositions[i]->pos) + (seedPositions[i]->size));
						}
						else prevk = sizeRefGenome;
						checkGap = ((s != 0) && ( // is not the rightmost/first seed id, and
										(i == (numTotalSeedHits - 1)) // is the last seed hit
										|| ((seedPositions[i + 1]->seed) != (s - 1)) // or the next seed hit does not have the next seed id
										|| (gapSizeInGenome > maxIntronLength) // or the gap length with the other seed does not already correspond to an intron
										// TODO: when non-linear chains are available, re-add the check "||(gapSizeInGenome<minIntronLength)" here too;
										// now it's not useful because the new extra seed would be added after the next seed and would not be checked/chained (and would cause errors too!)
										));
						if (checkGap) { // check if there is an existing splice site to the right
							k += ((seedPositions[i]->size) - 1); // rightmost position of the seed hit
							ssp = k;
							p = GetOppositePositionIfSpliceSiteExists(&ssp); // get real splice site position (ssp) and opposite position (p)
							if ((p != 0) && (p > ssp)) { // if a splice site does exist to the right
								d = (int)(k - ssp); // difference between seed ending position (possibly to the right) and splice site position (possibly to the left)
								n = ((seedPositions[i]->readpos) + (seedPositions[i]->size)); // next position on the right to match in the read
								n -= d; // account for that difference in the seed's ending position in the read
								if (n < readSize) {
									m = 0; // number of matched chars
									charPtrAtRead = (char *)(read + n); // pointer to the char at the right of the ending of the seed (corrected for the difference)
									charPtrAtGenome = (char *)(refGenome + p); // pointer to the char at the leftmost end of the right exon
									while ((*charPtrAtRead) == (*charPtrAtGenome)) { // while the characters in the read and in the genome match,
										m++;
										n++;
										p++;
										if (((*charPtrAtRead) == '\0') || ((*charPtrAtGenome) == '\0')) break; // ensure that the new seed does not overlap the existing one on the right
										charPtrAtRead++; // keep advancing both pointer to the right
										charPtrAtGenome++;
									}
									if (m > 0) { // if some chars were matched, create a new seed to the right
										seedPositions[i]->size -= d; // but first change the current seed to end exactly at the splice site's (left) position
										numSeedsAddedBySpliceSites++;
										numExtraSeedHits++;
										d = (numTotalSeedHits + numExtraSeedHits); // new seed hit id (+1)
										if (d > maxSeedArraysSize) { // reallocate the arrays if needed
											fwdSeedPositions = (SeedPos **)realloc(fwdSeedPositions, d * sizeof(SeedPos *));
											for (j = maxSeedArraysSize; j < d; j++) {
												fwdSeedPositions[j] = (SeedPos *)malloc(sizeof(SeedPos));
											}
											revSeedPositions = (SeedPos **)realloc(revSeedPositions, d * sizeof(SeedPos *));
											for (j = maxSeedArraysSize; j < d; j++) {
												revSeedPositions[j] = (SeedPos *)malloc(sizeof(SeedPos));
											}
											maxSeedArraysSize = d;
											fwdSeedChainSizes = (int *)realloc(fwdSeedChainSizes, maxSeedArraysSize * sizeof(int));
											fwdSeedChainScores = (int *)realloc(fwdSeedChainScores, maxSeedArraysSize * sizeof(int));
											fwdSeedChainStarts = (int *)realloc(fwdSeedChainStarts, maxSeedArraysSize * sizeof(int));
											revSeedChainSizes = (int *)realloc(revSeedChainSizes, maxSeedArraysSize * sizeof(int));
											revSeedChainScores = (int *)realloc(revSeedChainScores, maxSeedArraysSize * sizeof(int));
											revSeedChainStarts = (int *)realloc(revSeedChainStarts, maxSeedArraysSize * sizeof(int));
											bestSeedChains = (int *)realloc(bestSeedChains, maxSeedArraysSize * sizeof(int));
											bestSeedChainsStrands = (char *)realloc(bestSeedChainsStrands, maxSeedArraysSize * sizeof(char));
											if (strand == '+') { // reload modified array variables for the current strand
												seedPositions = fwdSeedPositions;
												seedChainSizes = fwdSeedChainSizes;
												seedChainScores = fwdSeedChainScores;
												seedChainStarts = fwdSeedChainStarts;
											}
											else {
												seedPositions = revSeedPositions;
												seedChainSizes = revSeedChainSizes;
												seedChainScores = revSeedChainScores;
												seedChainStarts = revSeedChainStarts;
											}
										} // realloc arrays
										d--; // d-th seed goes in the (d-1)-th position
										seedPositions[d]->seed = (s - 1); // seed id of the seed to the right
										seedPositions[d]->size = m;
										seedPositions[d]->pos = (p - m); // position of first matched char (the current position is a mismatch)
										seedPositions[d]->readpos = (n - m);
										// TODO: correctly fill this with splice site id to indicate that the seed was created from that splice site, and don't have to go get the id later on
										seedPositions[d]->ssId = 1;
									} // create new seed
								}
							} // splice site check
						} // left gap check
						prevs = s;
						prevk = (seedPositions[i]->pos) + (seedPositions[i]->size) - 1; // rightmost char position of this seed
					} // end of loop for all seed hits
				} // RNA-Seq mode
				// TODO: add extra linked list pointers in the seeds' struct, so they do not need to be sorted again
				if (numExtraSeedHits != 0) { // sort the seeds again, in case some new ones were added
					numTotalSeedHits += numExtraSeedHits;
					qsort(seedPositions, numTotalSeedHits, sizeof(SeedPos *), CompareSeedPositions);
				}
				/********************/
				/* End RNA-Seq Code */
				/********************/
				#endif

				
				#ifdef DEBUG
				fprintf(debugfile, "%d", numTotalSeedHits);
				for (n = 0; n < numTotalSeedHits; n++) fprintf(debugfile, " %s[%d]%d:%u(%u)", ((seedPositions[n]->ssId) ? "!" : ""), (n + 1), ((seedPositions[n]->seed) + 1), (seedPositions[n]->pos), (seedPositions[n]->size));
				fprintf(debugfile, "\n");
				#endif

				// TODO: when the SeedChain struct exists, add a field for "gapSizeInGenome"
				bestSeedChainsScore = 0; // if initialized once for each strand, will keep best chains for each strand separately, if not, it will only keep best chains among both strands
				numSeedChains = 0;
				s = seedPositions[0]->seed;
				k = seedPositions[0]->pos;
				seedChainSizes[numSeedChains] = 1;
				seedChainScores[numSeedChains] = seedPositions[0]->size;
				seedChainStarts[numSeedChains] = 0;
				prevs = s; // values from the previous seed in the sorted array
				prevk = k;
				for( i = 1 ; i < numTotalSeedHits ; i++ ){ // find all groups of consecutive seeds (chains) in the sorted array that have consistent positions among them
					s = seedPositions[i]->seed;
					k = seedPositions[i]->pos;
					// check if the seeds appear in the same order in the read and in the genome and have approximate relative positions
					if( s < prevs ){ // lower than because the seed identifier grows from right to left (from the end to the beginning)
						m = (seedPositions[i]->readpos - seedPositions[i-1]->readpos); // distance between both seeds in the read
						n = (int)( k - prevk ); // distance between both seeds in the reference genome
						if( n < 0 ) goto _stop_chain; // prevent overflow if distance is larger than 2^31
						m = abs( n - m ); // absolute difference between both distances				
						if( m > maxReadErrors ){ // if the difference is more than the errors we have to spare, stop the chain
							// TODO: if RNA-Seq mode, also check if seeds are followed (e.g.: s and s+1)
							// TODO: check case where the read is very long, and the maxNumErrors > minIntronSize, and there could still be a splicing event
							/**********************/
							/* Begin RNA-Seq Code */
							/**********************/
							// TODO: even if we are in RNA-Seq mode and the intron lengths are valid, it currently only processes the seeds if they have consecutive ids,
							//       but it should check for simple errors/indels which are the likely cause of 1 or at most 2 seeds being missing in the middle
							if( (!rnaSequencingMode) // in case we are in RNA-Seq mode, check if the distance does not correspond to an intron in the genome before quitting
								|| (n < minIntronLength) || (n > maxIntronLength) ) goto _stop_chain;
							/********************/
							/* End RNA-Seq Code */
							/********************/
						}
						// TODO: check if this only happens in RNA-Seq mode or in regular mode too
						m = ((seedPositions[i - 1]->readpos) + (seedPositions[i - 1]->size)); // end (+1) of the previous seed
						if ((seedPositions[i]->readpos) < m) { // if, on the read, the beggining of this seed overlaps the end of the previous one,
							// TODO: ensure that this is not <=0 ?
							(seedPositions[i - 1]->size) -= (m - (seedPositions[i]->readpos)); // cut the overlapping size of the previous seed
						}
						if( (k + (seedPositions[i]->size)) <= (prevk + (seedPositions[i-1]->size)) ) goto _stop_chain; // if this seed occurs in the genome completely inside the previous (left) seed, discard it from the chain
						seedChainSizes[numSeedChains]++; // one more seed in the current chain
						seedChainScores[numSeedChains] += seedPositions[i]->size; // add this seed size of the chain score
						prevs = s; // store current values to check against next seed
						prevk = k;
						continue; // go to next seed
					} // if the current seed cannot be incorporated in the current chain, begin a new chain starting at this seed
					_stop_chain:
					if( seedChainScores[numSeedChains] > bestSeedChainsScore ){ // first check if this is a new top scoring chain
						bestSeedChainsScore = seedChainScores[numSeedChains];
						//numBestSeedChains = 0;
					}
					numSeedChains++; // new chain
					seedChainSizes[numSeedChains] = 1; // still with just one seed
					seedChainScores[numSeedChains] = seedPositions[i]->size; // of this size
					seedChainStarts[numSeedChains] = i; // starting at current position
					prevs = s; // store current values to check against next seed
					prevk = k;
				}
				if( seedChainScores[numSeedChains] > bestSeedChainsScore ){ // check if last chain is a new top scoring chain
					bestSeedChainsScore = seedChainScores[numSeedChains];
					//numBestSeedChains = 0;
				}
				numSeedChains++; // update counters: chain at pos i is the (i+1)-th chain
				#ifdef DEBUG
				fprintf(debugfile,"%d",numSeedChains);
				for( n = 0 ; n < numSeedChains ; n++ ) fprintf(debugfile," [%d](%d)%d,",(seedChainStarts[n]+1),seedChainSizes[n],seedChainScores[n]);
				fprintf(debugfile,"\n");
				#endif

				// TODO: test if we should keep this part or only the below
				/*
				if( tryHarder ){ // if we need to try harder, keep all seed chains, no matter their score
					for( i = 0 ; i < numSeedChains ; i++ ){
						bestSeedChains[numBestSeedChains] = i;
						bestSeedChainsStrands[numBestSeedChains] = strand;
						numBestSeedChains++;
					}
					
				} else { // if we do not need to try harder, keep only the top scoring seed chains
					for( i = 0 ; i < numSeedChains ; i++ ){
						if( seedChainScores[i] == bestSeedChainsScore ){ // check if it has the same score as the top scoring chain
							bestSeedChains[numBestSeedChains] = i;
							bestSeedChainsStrands[numBestSeedChains] = strand;
							numBestSeedChains++;
						}
					}
				}
				*/
				/**/
				// TODO: when SeedChain struct exists, add linked pointer to next best chain
				for( i = 0 ; i < numSeedChains ; i++ ){ // test only the top scoring chains for each strand, instead of all chains
					if( seedChainScores[i] == bestSeedChainsScore ){ // check if it has the same score as the top scoring chain
						bestSeedChains[numBestSeedChains] = i;
						bestSeedChainsStrands[numBestSeedChains] = strand;
						numBestSeedChains++;
					}
				}
				/**/
				if( testOtherStrandToo ){ // process the other strand too if needed
					testOtherStrandToo = 0;
					if( strand == '+' ) strand = '-';
					else strand = '+';
					continue;
				}
				break;
			} // end of finding best chains for both strands loop
		
			#ifdef DEBUG
			fprintf(debugfile,"%d",numBestSeedChains);
			for( n = 0 ; n < numBestSeedChains ; n++ ){
				strand=bestSeedChainsStrands[n];
				if(strand=='+'){
					seedChainStarts=fwdSeedChainStarts;
					seedChainSizes=fwdSeedChainSizes;
					seedChainScores=fwdSeedChainScores;
				} else {
					seedChainStarts=revSeedChainStarts;
					seedChainSizes=revSeedChainSizes;
					seedChainScores=revSeedChainScores;
				}
				fprintf(debugfile," %c[%d](%d)%d,",strand,(seedChainStarts[bestSeedChains[n]]+1),seedChainSizes[bestSeedChains[n]],seedChainScores[bestSeedChains[n]]);
			}
			fprintf(debugfile,"\n");
			fflush(debugfile);
			#endif

			/*
			if( tryHarder ){ // increase number of seeds, because the last/first seed is incomplete to the left, so when processing the (numSeed-1)-th seed, it will do DP to the left
				fwdNumSeeds++;
				revNumSeeds++;
			}
			*/

			if( numBestSeedChains > maxReadHitsArraySize ){ // resize the read hit lists if needed
				IncreaseMaxReadHitsArraySize(numBestSeedChains,(maxReadSize+maxNumErrors));
				readHits=readsHits[currentPair]; // set pointer again because memory address probably changed
			}

			numReadHits = 0;
			bestReadHitId = -1;
			numErrors2ndBestHit = readSize; // 2nd lower number of alignment errors among all the hits/positions of the read
			bestReadNumErrors = readSize; // minimum number of alignment errors among all the hits/positions of the read
			
			//TODO: (done?) use each seed's updated size (seedPositions[i]->size) instead of original one (seedSizes[s]) which is the same for all seeds
			for( j = 0 ; j < numBestSeedChains ; j++ ){ // process all top scoring seed chains
				cigarStrCodes = readHits[numReadHits].cigarcodes; // set pointers to arrays in hits list
				cigarStrCounts = readHits[numReadHits].cigarcounts;
				cigarStrSize = 0; // reset cigar string for use with this read position
				cigarStrCodes[0] = 'M'; // start with a 0-size match
				cigarStrCounts[0] = 0;
				numReadErrors = 0; // start error count for this hit/position
				numIndelsAtEnd = 0;
				strand = bestSeedChainsStrands[j];
				if( strand == '+' ){ // load variables for corresponding chain strand
					seedChainStarts = fwdSeedChainStarts;
					seedChainSizes = fwdSeedChainSizes;
					seedPositions = fwdSeedPositions;
					read = fwdRead;
					numSeeds = fwdNumSeeds;
					seedSizes = fwdSeedSizes;
					seedStartPosInRead = fwdSeedStartPosInRead;
				} else {
					seedChainStarts = revSeedChainStarts;
					seedChainSizes = revSeedChainSizes;
					seedPositions = revSeedPositions;
					read = revRead;
					numSeeds = revNumSeeds;
					seedSizes = revSeedSizes;
					seedStartPosInRead = revSeedStartPosInRead;
				}
				k = bestSeedChains[j]; // current chain
				i = seedChainStarts[k]; // first seed of chain
				n = ( i + seedChainSizes[k] ); // last (+1) seed of chain
				m = 0; // number of seeds already processed
				s = seedPositions[i]->seed; // information of first seed
				p = seedPositions[i]->pos;
				if( s == (numSeeds-1) && seedStartPosInRead[s] <= 1 ){ // if it's the left-most seed (and it starts at the 0-th or 1-st positions), no need for DP, just store seed size and go to next seed
					if( seedStartPosInRead[s] == 1 ){ // if the left-most seed ended with a mismatch (in the first char of the read), automatically consider it a mismatch (without DP)
						seedStartPosInRead[s]--;
						(seedPositions[i]->size)++;
						(seedPositions[i]->readpos)--;
						if(p!=0) p--;
						numReadErrors++;
					}
					gapStartPosInRead = (seedPositions[i]->size); // position in read next to end of chain's first seed
					gapStartPosInGenome = ( p + (seedPositions[i]->size) ); // position in genome next to end of chain's first seed
					/*
					if( i != (n-1) ){ // if it's not the last seed of the chain, limit its extended size if needed
						d = ( gapStartPosInRead - seedStartPosInRead[(seedPositions[(i+1)]->seed)] ); // distance between the end of this extended seed and the start of the next seed
						if( d > 0 ){ // if the current seed was extended (to the right) over the next seed (on the right), limit its size to the beginning of the next seed
							(seedPositions[i]->size) -= d;
							gapStartPosInGenome -= d;
							gapStartPosInRead -= d;
						}
					}
					*/
					UpdateCigarString('M',(seedPositions[i]->size));
					#ifdef DEBUG
					fprintf(debugfile,"[Seed:%02d/%02d](%02d/%02d) %d\n",(s+1),numSeeds,(m+1),seedChainSizes[bestSeedChains[j]],(seedPositions[i]->size));
					fprintf(debugfile,"%.*s\n",(seedPositions[i]->size),(char *)(refGenome+p));
					fprintf(debugfile,"%.*s\n",(seedPositions[i]->size),(char *)(read+seedStartPosInRead[s]));
					#endif
					if( numReadErrors == 1 ){ // if the left-most seed ended with a mismatch (the number of errors here would otherwise be 0)
						seedStartPosInRead[s]++; // change the variables back to their original values (because we need to decrease the position again in all the occurrences of this same left-most seed)
						(seedPositions[i]->size)--;
						(seedPositions[i]->readpos)++;
						if( maxReadErrors == 0 ){ // if no errors are allowed, we already have an error too much
							#ifdef DEBUG
							fprintf(debugfile,"%c%u ~%d X\n",strand,readPosInGenome,numReadErrors);
							#endif
							continue; // go to next seed chain
						}
					}
					distanceRelaxation = 0;
					readPosInGenome = p; // get read position from position of first (left-most) seed
					m = 1; // one seed processed
					i++; // go to next seed
				} else { // if not the left-most seed, define start positions so the starting gap can be aligned by DP in the next step
					distanceRelaxation = (2*numSeeds); // distance relaxation for DP on left end segment
					d = ( maxReadErrors + distanceRelaxation ); // number of errors that will be considered in next DP call
					if( d > maxNumErrors ) distanceRelaxation = ( maxNumErrors - maxReadErrors ); // check if we will not exceed the available maximum columns in the DP matrix
					gapStartPosInRead = 0;
					gapStartPosInGenome = ( p - (seedPositions[i]->readpos) - distanceRelaxation ); // allow some relaxation to allow indels in this segment
					if( gapStartPosInGenome > p ){ // if segment extends beyond beggining of genome ( < 0 if it was signed int)
						gapStartPosInGenome = 0; // limit gap to beginning of genome
						gapSizeInRead = (seedPositions[i]->readpos);
						gapSizeInGenome = (int)p;
						d = abs(gapSizeInGenome-gapSizeInRead); // difference in the sizes of the segments
						if( d > maxReadErrors ){ // if the alignment would have more errors than we have to spare, do not do it
							#ifdef DEBUG
							fprintf(debugfile,"%c%u ~%d X\n",strand,readPosInGenome,numReadErrors);
							#endif
							continue; // go to next seed chain
						}
						if( d < distanceRelaxation ) distanceRelaxation = d; // only allow this distance relaxation
					}
					readPosInGenome = gapStartPosInGenome; // store read start position (temporary because it will be updated by numIndelsAtEnd)
				} // end of processing first segment (not necessarilly the left-most)
				for( ; i < n ; i++ ){ // process all (middle) seeds from this seed chain (consecutive positions in the array with all the seeds), and check alignment of the gap between the current seed and the previous one
					s = seedPositions[i]->seed;
					p = seedPositions[i]->pos;
					gapSizeInRead = ( (seedPositions[i]->readpos) - gapStartPosInRead );
					gapSizeInGenome = ( p >= gapStartPosInGenome ) ? (int)( p - gapStartPosInGenome ) : -(int)( gapStartPosInGenome - p ) ; // get signed difference of unsigned number
					#ifdef DEBUG
					fprintf(debugfile,"{%u-%u} (%d)\n",gapStartPosInGenome,(p-1),gapSizeInGenome);
					fprintf(debugfile,"{%d-%d} (%d)\n",gapStartPosInRead,(seedStartPosInRead[s]-1),gapSizeInRead);
					#endif
					if (m == 0) { // if we are dealing with the left most gap
						if (gapSizeInGenome >= (2 * gapSizeInRead)) { // at most it will be all insertions followed by all deletions, which is the same that all mismatches
							gapSizeInGenome = (2 * gapSizeInRead - 1); // reduce gap size
							gapStartPosInGenome = (p - gapSizeInGenome); // update already set positions to the left accordingly
							readPosInGenome = gapStartPosInGenome;
						}
					} else { // if it's not the left most gap
						// TODO: if a seed was added because of splicing ((seedPositions[i]->ssId)!=0) skip some of the code bellow and directly set an intron
						// TODO: merge some of the code bellow (e.g. call to "IsSpliceSite") with the RNA-Seq code in the seed chaining section
						// TODO: check case where the read is very long, and the maxNumErrors > minIntronSize, and there could still be a splicing event
						if( abs(gapSizeInGenome-gapSizeInRead) > (maxReadErrors-numReadErrors) ){ // if the difference in the segment sizes is larger than the number of errors that we can still use, stop checking this chain
							/**********************/
							/* Begin RNA-Seq Code */
							/**********************/
							if(!rnaSequencingMode) break; // leave this seed chain and go to next one
							if (gapSizeInRead > 0) break; // if the seeds are not right next to each other on the read, we will not be able to check for splicing accurately
							if (gapSizeInRead < 0) { // if the left seed was extended over the right seed (in the read, not in the genome), fix it
								gapSizeInRead = (-gapSizeInRead); // set positive value
								(seedPositions[(i-1)]->size) -= gapSizeInRead; // update (shorten) the size of the previous/left seed
								UpdateCigarString('M', (-gapSizeInRead)); // update with negative value
								gapStartPosInRead -= gapSizeInRead;
								gapStartPosInGenome -= gapSizeInRead;
								gapSizeInGenome += gapSizeInRead;
								gapSizeInRead = 0; // we will only check for splicing if the seeds are right next to each other in the read
							}
							// TODO: check if this is needed, because if we get here, the intron lengths have already been checked when building the seed chains
							if( (gapSizeInGenome<minIntronLength) || (gapSizeInGenome>maxIntronLength)) break; // check if the intron length is valid
							gapSizeInRead = IsSpliceSite(refGenome, gapStartPosInGenome, gapSizeInGenome); // check for splice-site signals
							if( gapSizeInRead == INT_MIN ) break; // no splicing was found
							if (gapSizeInRead != 0) { // if the seeds shifted to the left or to the right, update variables
								gapStartPosInRead += gapSizeInRead;
								gapStartPosInGenome += gapSizeInRead;
								(seedPositions[(i - 1)]->size) += gapSizeInRead; // change size of left (previous) seed
								UpdateCigarString('M', gapSizeInRead);
								(seedPositions[i]->size) -= gapSizeInRead; // change size of right (current) seed
								(seedPositions[i]->pos) += gapSizeInRead; // move starting genome position of right (current) seed
								p += gapSizeInRead; // update the value in this variable too because it will be used ahead
								(seedPositions[i]->readpos) += gapSizeInRead; // move starting read position of right (current) seed
							}
							#ifndef DISABLE_SS_ARRAY
							AddOrUpdateSpliceSite((gapStartPosInGenome - 1), (gapStartPosInGenome + gapSizeInGenome)); // save information about this splice site
							#endif
							UpdateCigarString('N',gapSizeInGenome); // update cigar string with intron length
							// TODO: (?) if there's still a part of the read left to be aligned, align it with DP
							gapSizeInRead = 0; // reset used variable to real gap size in read, which is zero
							gapSizeInGenome = 0; // reset the space in the genome, so it will no be considered a deletion next
							numSplicingEvents++;
							splicedRead = 1;
							/********************/
							/* End RNA-Seq Code */
							/********************/
						} // end of max num errors check
					} // end of leftmost gap check
					// TODO: see if we can use (seedPositions[i]->readpos) to do this instead of using the extra variable next
					numSeedCharsToSkip = 0;
					if( ( gapSizeInRead <= 0 ) || ( gapSizeInGenome <= 0 ) ){ // directly deal with simple cases of small insertions and deletions
						if( gapSizeInRead < gapSizeInGenome ){ // deletion
							gapSizeInRead = (-gapSizeInRead); // is <= 0, so get absolute value
							numSeedCharsToSkip = gapSizeInRead; // seed size to skip
							numDPErrors = ( gapSizeInRead + gapSizeInGenome ); // error size
							UpdateCigarString('D',numDPErrors); // update cigar string
						} else if( gapSizeInGenome < gapSizeInRead ){ // insertion
							gapSizeInGenome = (-gapSizeInGenome); // is <= 0, so get absolute value
							numSeedCharsToSkip = gapSizeInGenome; // seed size to skip
							numDPErrors = ( gapSizeInRead + gapSizeInGenome ); // error size
							UpdateCigarString('I',numDPErrors); // update cigar string
						} else { // both gaps are zero (both seeds are connected)
							// TODO: check if the only case we enter here is in fact both gaps = 0
							numDPErrors = 0;
						}
					} else { // use dynamic programming
						//TODO: check if this part ever happens and remove it if not needed
						/*
						if( gapSizeInGenome < 0 ){ // fix problem with overlapping seeds in genome
							numSeedCharsToSkip = (-gapSizeInGenome); // seed size to skip
							gapSizeInRead += numSeedCharsToSkip; // add size of overlap
							gapSizeInGenome = 0;
						}
						*/
						numDPErrors = RunBandedDynamicProgramming( (char *)( refGenome + gapStartPosInGenome ) , gapSizeInGenome, (char *)( read + gapStartPosInRead ) , gapSizeInRead , ( maxReadErrors - numReadErrors + distanceRelaxation ) , m , &numIndelsAtEnd );
						//UpdateCigarStringWithDP(m);
					}
					numReadErrors += numDPErrors;
					if( m == 0 ){ // if this is the first seed to be processed (but not the left-most one), update read position
						readPosInGenome += numIndelsAtEnd;
						distanceRelaxation = 0; // remove distance relaxation for next segments
					}
					#ifdef DEBUG
					fprintf(debugfile,"%d %d\n",numDPErrors,dpTraceSize);
					fprintf(debugfile,"[Seed:%02d/%02d](%02d/%02d) %d\n",(s+1),numSeeds,(m+1),seedChainSizes[bestSeedChains[j]],(seedPositions[i]->size));
					fprintf(debugfile,"%.*s\n",(seedPositions[i]->size),(char *)(refGenome+p));
					fprintf(debugfile,"%.*s\n",(seedPositions[i]->size),(char *)(read+(seedPositions[i]->readpos)));
					#endif
					if( numReadErrors > maxReadErrors ) break; // if we have too many errors already, stop checking this chain
					gapStartPosInRead = ( (seedPositions[i]->readpos) + (seedPositions[i]->size) ); // update gap start positions for next alignment
					gapStartPosInGenome = ( p + (seedPositions[i]->size) );
					/*
					if( i != (n-1) ){ // if it's not the last seed of the chain, limit its extended size if needed
						d = ( gapStartPosInRead - seedStartPosInRead[(seedPositions[(i+1)]->seed)] ); // distance between the end of this extended seed and the start of the next seed
						if( d > 0 ){ // if the current seed was extended (to the right) over the next seed (on the right), limit its size to the beginning of the next seed
							(seedPositions[i]->size) -= d;
							gapStartPosInGenome -= d;
							gapStartPosInRead -= d;
						}
					}
					*/
					UpdateCigarString('M',((seedPositions[i]->size)-numSeedCharsToSkip)); // add seed information to cigar string
					m++; // one more processed seed
				} // end seeds loop
				if( i != n ){ // if we did not get to the end of the loop (process all seeds of this chain) because of too many errors
					#ifdef DEBUG
					fprintf(debugfile,"%c%u ~%d X\n",strand,readPosInGenome,numReadErrors);
					#endif
					continue; // go to next seed chain
				}
				if( (s != 0) || (gapStartPosInRead != readSize) ){ // if the last processed seed was not the right-most one (or if it was an "incomplete" seed created by a Splice Site check), fill the remaining right ending gap with DP
					distanceRelaxation = (2*numSeeds); // distance relaxation for DP on right end segment
					d = ( maxReadErrors - numReadErrors + distanceRelaxation ); // number of errors to consider in next DP call
					if( d > maxNumErrors ) distanceRelaxation = ( maxNumErrors - ( maxReadErrors - numReadErrors ) );
					gapSizeInRead = ( readSize - gapStartPosInRead );
					gapSizeInGenome = ( gapSizeInRead + distanceRelaxation ); // allow small size relaxation at end too
					if( gapSizeInGenome >= (2*gapSizeInRead) ) gapSizeInGenome = ( 2*gapSizeInRead - 1 ); // at most it will be all insertions/deletions
					if( ( gapStartPosInGenome + gapSizeInGenome ) >= sizeRefGenome ){ // if segment extends beyhond end of genome
						gapSizeInGenome = (int)( sizeRefGenome - gapStartPosInGenome ); // limit gap to end of genome
						d = abs(gapSizeInGenome-gapSizeInRead); // difference in the sizes of the segments
						if( d > (maxReadErrors-numReadErrors) ){ // if the alignment would have more errors than we have to spare, do not do it
							#ifdef DEBUG
							fprintf(debugfile,"%c%u ~%d X\n",strand,readPosInGenome,numReadErrors);
							#endif
							continue; // go to next seed chain
						}
						if( d < distanceRelaxation ) distanceRelaxation = d; // only allow this distance relaxation
					}
					#ifdef DEBUG
					fprintf(debugfile,"{%u-%u} (%d)\n",gapStartPosInGenome,(gapStartPosInGenome+gapSizeInGenome-1),gapSizeInGenome);
					fprintf(debugfile,"{%d-%d} (%d)\n",gapStartPosInRead,(readSize-1),gapSizeInRead);
					#endif
					m = -1; // segment at the right end
					numDPErrors = RunBandedDynamicProgramming( (char *)( refGenome + gapStartPosInGenome ) , gapSizeInGenome, (char *)( read + gapStartPosInRead ) , gapSizeInRead , ( maxReadErrors - numReadErrors + distanceRelaxation ) , m , &numIndelsAtEnd );
					numReadErrors += numDPErrors;
					//UpdateCigarStringWithDP(m);
					#ifdef DEBUG
					fprintf(debugfile,"%d %d\n",numDPErrors,dpTraceSize);
					#endif
					if( numReadErrors > maxReadErrors ){ // check if we exceeded the maximum number of errors in this last segment
						#ifdef DEBUG
						fprintf(debugfile,"%c%u ~%d X\n",strand,readPosInGenome,numReadErrors);
						#endif
						continue; // go to next seed chain
					}
				} // end of processing the right-most gap if it exists
				#ifdef DEBUG
				fprintf(debugfile,"%c%d %d\n",strand,readPosInGenome,numReadErrors);
				#endif
				// if we got to here, it means we have less than the maximum number of errors
				if( numReadErrors < bestReadNumErrors ){ // check if this is a new top scoring hit
					numErrors2ndBestHit = bestReadNumErrors; // the best is now the 2nd best
					bestReadNumErrors = numReadErrors;
					bestReadHitId = numReadHits;
				} else { // if it is not the best hit, check if it is the 2nd best hit
					if( numReadErrors < numErrors2ndBestHit ) numErrors2ndBestHit = numReadErrors;
				}
				readHits[numReadHits].pos = readPosInGenome;
				readHits[numReadHits].strand = strand;
				readHits[numReadHits].numerrors = numReadErrors;
				readHits[numReadHits].score = ( readSize - numReadErrors );
				readHits[numReadHits].cigarsize = cigarStrSize; // the CIGAR arrays for the next hit will be correctly set when we go back to the beginning of the loop
				numReadHits++;
				if( tryHarder ){ // if the occurrence was found in hard mode, stop checking the other candidates
					numAlignedHardReads++;
					break;
				}
			} // end best seed chains loop

			// NOTE: these stats updates need to be here and not at _report_read because there they would not be updated when in paired-end reads mode
			if( numReadHits > 0 ){ // if we got hits for this read, the read was mapped
				_mapped_read:
				numAlignedReads++; // update statistics
				totalNumMappedReadsChars += readSize;
				totalNumReadHits += numReadHits;
				if( numReadHits > (int)maxReadHits ) maxReadHits = numReadHits;
				totalNumDPErrors += bestReadNumErrors;
				totalNumSeeds += numSeeds;
				totalNumSeedOccs += numTotalSeedHits;
				totalNumSeedChains += numBestSeedChains;
				if( rnaSequencingMode && splicedRead ) numSplicedReads++;
				break; // report read and go to next one
			}
			if( tryHarder ){ // if we already did try harder but it did not work, the read was not mapped
				_unmapped_read:
				numUnalignedReads++; // update statistics
				#ifdef DEBUG
				fprintf(debugfile,"X\n");
				#endif
				break; // report read and go to next one
			} else { // if the read was not aligned but we did not try harder yet
				tryHarder = 1; // try aligning read in hard mode
				#ifdef DEBUG
				fprintf(debugfile,"*\n");
				#endif
			}

		} // end of normal/agressive method loop

		if( pairedEndReadsMode ){ // when finished mapping, save variables for the current read pair
			readsSizes[currentPair] = readSize;
			numReadsHits[currentPair] = numReadHits; // even if we get here skipping the best chains part, this var has already been initialized before
			bestReadsHitIds[currentPair] = bestReadHitId;
			bestReadsScores2nd[currentPair] = ( readSize - numErrors2ndBestHit );
			if( currentPair == 1 ) PrintReadPair(); // only output the reads if we finished processing both of them
		} else { // single read mode
			if( numReadHits > 0 ){ // if we got hits for this read, the read was correctly mapped
				// TODO: add this check to paired-ends mode too
				if( singleRefReadsOnly ){ // if we only wants reads that map to one ref, check all hits refs
					m = GetReadRef(readHits[0].pos); // ref for the first hit
					for( i = 1 ; i < numReadHits ; i++ ){
						if( GetReadRef(readHits[i].pos) != m ){ // if a hit has a different ref, discard read
							numReadHits = 0;
							numAlignedReads--;
							goto _unmapped_read;
						}
					}
				}
				i = 0;
				n = numReadHits;
				if( !reportAllHits ){ // if we only want the best hit
					i = bestReadHitId;
					n = (i+1);
				}
				for( ; i < n ; i++ ){
					cigarStrCodes = readHits[i].cigarcodes;
					cigarStrCounts = readHits[i].cigarcounts;
					cigarStrSize = readHits[i].cigarsize;
					//PrintReadInSAMFormat(read,readName,strand,readPosInGenome,numReadErrors,(readSize-numReadErrors),(readSize-numErrors2ndBestHit),readSize);
					PrintReadInSAMFormat(((readHits[i].strand=='+')?(fwdRead):(revRead)),readName,readHits[i].strand,readHits[i].pos,readHits[i].numerrors,readHits[i].score,(readSize-numErrors2ndBestHit),readSize);
				}
			} else { // if the read was not mapped
				PrintReadInSAMFormat(fwdRead,readName,'+',-1,0,0,0,0); // report in SAM that the read was not aligned
				if( reportUnalignedReads ) fprintf(unalignedReadsFile,">%s\n%s\n",readName,fwdRead);
			}
		} // report read in single read mode

	} // end of loop for all the reads

	free(tempRead.fwdStrand.seeds);
	free(tempRead.revStrand.seeds);
	free(tempRead.fwdStrand.topPointers);
	free(tempRead.revStrand.topPointers);
	free(tempRead.fwdStrand.bottomPointers);
	free(tempRead.revStrand.bottomPointers);
	for( i = 0 ; i < dpNumRows ; i++ ){
		free(dpMatrix[i]);
		free(dpDirections[i]);
	}
	free(dpMatrix);
	free(dpDirections);
	#ifdef DEBUGDP
	free(dpTrace);
	free(dpAlignedTarget);
	free(dpAlignedRead);
	#endif
	for(i=0;i<maxSeedArraysSize;i++){
		free(fwdSeedPositions[i]);
		free(revSeedPositions[i]);
	}
	free(fwdSeedPositions);
	free(revSeedPositions);
	free(fwdSeedChainSizes);
	free(fwdSeedChainScores);
	free(fwdSeedChainStarts);
	free(revSeedChainSizes);
	free(revSeedChainScores);
	free(revSeedChainStarts);
	free(bestSeedChains);
	free(bestSeedChainsStrands);
	free(fwdSeedSizes);
	free(fwdNumSeedHits);
	free(fwdSeedStartPosInRead);
	free(fwdSeedLeftPtrs);
	free(revSeedSizes);
	free(revNumSeedHits);
	free(revSeedStartPosInRead);
	free(revSeedLeftPtrs);
	if( pairedEndReadsMode ) numPairs=2;
	else numPairs = 1;
	for(i=0;i<numPairs;i++){
		free(fwdReads[i]);
		free(revReads[i]);
		free(readNames[i]);		
		for(j=0;j<maxReadHitsArraySize;j++){
			free(readsHits[i][j].cigarcodes);
			free(readsHits[i][j].cigarcounts);
		}
		free(readsHits[i]);
	}

}


void LoadReferenceGenome(char *basefilename){
	FILE *genomefile;
	unsigned int n, i, k;
	char c;
	printf("> Opening reference genome file ");
	n=(int)strlen(FMI_GetTextFilename());
	genomefilename=(char *)calloc((n+1),sizeof(char));
	strcpy(genomefilename,FMI_GetTextFilename()); // get filename from internal index information
	genomefile=fopen(genomefilename,"r");
	if(genomefile==NULL){ // try ".fasta" extention
		free(genomefilename);
		genomefilename=AppendToBasename(basefilename,".fasta"); // get filename from input parameters
		genomefile=fopen(genomefilename,"r");
	}
	if(genomefile==NULL){ // try ".fas" extention
		n=0;
		while(genomefilename[n]!='\0') n++; // get number of chars of filename
		genomefilename[(n-2)]='\0';
		genomefile=fopen(genomefilename,"r");
	}
	if(genomefile==NULL){ // try ".fa" extention
		genomefilename[(n-3)]='\0';
		genomefile=fopen(genomefilename,"r");
	}
	if(genomefile==NULL){ // try ".fna" extention
		genomefilename[(n-4)]='n';
		genomefilename[(n-3)]='a';
		genomefilename[(n-2)]='\0';
		genomefile=fopen(genomefilename,"r");
	}
	if(genomefile==NULL){ // try ".fsa" extention
		genomefilename[(n-4)]='s';
		genomefile=fopen(genomefilename,"r");
	}
	if(genomefile==NULL){ // try ".txt" extention
		genomefilename[(n-5)]='t';
		genomefilename[(n-4)]='x';
		genomefilename[(n-3)]='t';
		genomefile=fopen(genomefilename,"r");
	}
	if(genomefile==NULL){ // report error
		genomefilename[(n-5)]='f';
		genomefilename[(n-4)]='a';
		genomefilename[(n-3)]='s';
		genomefilename[(n-2)]='t'; // back to ".fasta" extention
		printf("<%s> ... ",genomefilename);
		printf("\n> ERROR: Reference genome file not found\n");
		exit(-1);
	}
	printf("<%s> ... ",genomefilename);
	fflush(stdout);
	refGenome=(char *)malloc((sizeRefGenome+1)*sizeof(char));
	refGenome[sizeRefGenome]='\0';
	numRefs=1;
	refsEndPositions=(unsigned int *)calloc(numRefs,sizeof(unsigned int));
	refsEndPositions[0]=(sizeRefGenome-1);
	refsNames=(char **)calloc(numRefs,sizeof(char *));
	refsNames[0]=(char *)calloc((255+1),sizeof(char));
	c=fgetc(genomefile);
	if(c=='>'){ // get reference description (for 1st sequence)
		i=0;
		while(c!='\n' && c!=EOF){
			c=fgetc(genomefile);
			if(c=='\t') c=' '; // convert tabs to spaces
			if(i<255) refsNames[0][i++]=c;
		}
		if(i) i--;
		refsNames[0][i]='\0';
	} else { // no reference description found
		ungetc(c,genomefile); // put back valid character
		strcpy(refsNames[0],genomefilename); // load filename as genome description
	}
	k=0; // number of invalid chars (converted to N's)
	n=0; // current count of all usable chars inside file (including N's)
	while((c=fgetc(genomefile))!=EOF){
		if(c=='>'){ // new sequence inside reference file
			if(n) refsEndPositions[(numRefs-1)]=(n-1); // store position (in global seq) of the last valid char of the previous seq
			numRefs++; // one more seq
			refsEndPositions=(unsigned int *)realloc(refsEndPositions,numRefs*sizeof(unsigned int)); // realloc arrays
			refsNames=(char **)realloc(refsNames,numRefs*sizeof(char *));
			refsNames[(numRefs-1)]=(char *)calloc(256,sizeof(char));
			i=0;
			while(c!='\n' && c!=EOF){ // get label of the new sequence
				c=fgetc(genomefile);
				if(c=='\t') c=' '; // convert tabs to spaces
				if(i<255) refsNames[(numRefs-1)][i++]=c;
			}
			if(i) i--;
			refsNames[(numRefs-1)][i]='\0';
			if(c==EOF){ // invalidate last seq if it ends in the description
				free(refsNames[(numRefs-1)]);
				numRefs--;
			}
			continue; // get next char
		}
		if(c>=97 && c<=122) c=(char)(c-32); // lowercase alphabet letter, convert to uppercase
		if(c>=65 && c<=90){ // uppercase alphabet letter
			if(c!='A' && c!='C' && c!='G' && c!='T'){
				c='N'; // all other non ACGT alphabet letters are converted to N's
				k++;
			}
			if(n<sizeRefGenome) refGenome[n]=c; // assure that we do not exceed the ref size that we got from the ref index
			n++;
		}
	}
	if(n) refsEndPositions[(numRefs-1)]=(n-1); // store position of the last valid char of the last seq
	if(n!=sizeRefGenome){
		printf("\n> ERROR: Genome size in index (%u) and in original file (%u) do not match",sizeRefGenome,n);
		printf("\n>        Please use the correct reference file or (re)build the index for this reference\n");
		exit(-1);
	}
	fclose(genomefile);
	if(numRefs>1){ printf("(");PrintUnsignedNumber((unsigned int)numRefs);printf(" sequences) "); }
	printf("(");PrintUnsignedNumber(sizeRefGenome);printf(" basepairs) ");
	if(k!=0){ printf("(");PrintUnsignedNumber(k);printf(" N's) "); }
	printf("OK\n");
}

// NOTE: uppercase file type means single-read mode and lowercase file types means paired-read mode
void LoadReadsFiles(){
	int n;
	char c;
	long long int filesize;
	readsFiles=(FILE **)calloc(numReadsFiles,sizeof(FILE *)); // allocate arrays for all reads files
	readsFileTypes=(char *)calloc(numReadsFiles,sizeof(char));
	minDists=(int *)calloc(numReadsFiles,sizeof(int));
	maxDists=(int *)calloc(numReadsFiles,sizeof(int));
	prevFilesData=0; // set variables to keep track of the progress while mapping
	totalFilesData=0;
	for(n=0;n<numReadsFiles;n++){
		printf("> Opening reads file <%s> ... ",readsFilenames[n]);
		if((readsFiles[n]=fopen(readsFilenames[n],"r"))==NULL){
			printf("\n> ERROR: Reads file not found\n");
			exit(-1);
		}
		printf("(");
		fseek(readsFiles[n],0L,SEEK_END);
		filesize=(long long int)ftell(readsFiles[n]);
		fseek(readsFiles[n],0L,SEEK_SET);
		PrintNumber(filesize);
		printf(" bytes) ");
		totalFilesData+=filesize;
		c=fgetc(readsFiles[n]); // decide file format based on the file's first char
		if(c=='>'){ // FASTA format
			AnalyzeReads=AnalyzeFastaReads;
			LoadNextRead=LoadNextReadFromFasta;
		} else if(c=='@'){ // FASTQ format
			AnalyzeReads=AnalyzeFastaQReads;
			LoadNextRead=LoadNextReadFromFastaQ;
		} else if(c=='.'){ // SFF format
			AnalyzeReads=AnalyzeSffReads;
			LoadNextRead=LoadNextReadFromSff;
		} else { // invalid format
			printf("\n> ERROR: Invalid reads file format\n");
			exit(-1);
		}
		rewind(readsFiles[n]); // back to beginning of file
		readsFile=readsFiles[n]; // set both these two varibles because they are going to be needed by the AnalyzeReads function
		readsFilename=readsFilenames[n];
		readsFileTypes[n]=AnalyzeReads(0,1); // analyze reads file and store file type
		if(pairedEndReadsMode && readsFileTypes[n]!='s') readsFileTypes[n]=(char)(readsFileTypes[n]+32); // if we are in (forced) paired-end reads mode, set file types accordingly (lowercase)
		printf("OK\n");
	}
	if(!pairedEndReadsMode && maxReadPairDistance>0){ // automatically enable paired-end reads mode if a read pair distance is set and there are no SFF files
		for(n=0;n<numReadsFiles;n++){ // check if we go through all files without finding any SFF
			c=readsFileTypes[n];
			if(c=='S' || c=='s') break;
		}
		if(n==numReadsFiles){
			pairedEndReadsMode=1;
			for(n=0;n<numReadsFiles;n++) readsFileTypes[n]=(char)(readsFileTypes[n]+32); // change all file types for paired-end reads mode
		}
	}
	if(!pairedEndReadsMode && numReadsFiles==2){ // automatically enable paired-end reads mode if there are only 2 files and both filenames end in "_1" and "_2"
		n=0;
		while( (c=(readsFilenames[0][n]))==(readsFilenames[1][n]) && c!='\0' ) n++;
		if( n!=0 && readsFileTypes[0]!='S' && readsFileTypes[1]!='S'
			&& (readsFilenames[0][n-1])=='_' && (readsFilenames[1][n-1])=='_'
			&& (readsFilenames[0][n])=='1' && (readsFilenames[1][n])=='2'
			&& (readsFilenames[0][n+1])=='.' && (readsFilenames[1][n+1])=='.' ){
				n++;
				while( (c=(readsFilenames[0][n]))==(readsFilenames[1][n]) && c!='\0' ) n++;
				if( (readsFilenames[0][n])=='\0' && (readsFilenames[1][n])=='\0' ){
					pairedEndReadsMode=1;
					readsFileTypes[0]=(char)(readsFileTypes[0]+32);
					readsFileTypes[1]=(char)(readsFileTypes[1]+32);
				}
		}
	}
	if(pairedEndReadsMode){ // check if read file pairs are correctly formed
		for(n=0;n<numReadsFiles;n++){
			c=readsFileTypes[n];
			if(c=='S'){ // SFF with no paired-ends
				printf("> ERROR: SFF reads file <%s> was not detected as containing paired-end reads\n",readsFilenames[n]);
				exit(-1);
			}
			if(c!='f' && c!='q') continue; // no need to check SFF or single read
			if(n==(numReadsFiles-1)){ // the 2nd pair does not exist
				printf("> ERROR: Read file <%s> does not have a 2nd pair\n",readsFilenames[n]);
				exit(-1);
			}
			if(c!=readsFileTypes[(n+1)]){ // the 2nd pair is not of the same format
				printf("> ERROR: Reads files <%s> and <%s> have distinct file formats\n",readsFilenames[n],readsFilenames[(n+1)]);
				exit(-1);
			}
			n++; // skip the next read file (correct pair)
		}
	} else { // if we are not in paired ends reads mode, but at least one of the SFF files is paired-end, set paired ends mode for next functions to set and print paired-end info
		for(n=0;n<numReadsFiles;n++){
			c=readsFileTypes[n];
			if(c=='f' || c=='q' || c=='s') pairedEndReadsMode=1;
		}
	}
}

/********/
/* MAIN */
/********/

int main(int argc, char *argv[]) {
	struct timeb starttb,endtb;
	double timetb;
	int n;
	unsigned int k;
	char *indexfilename, *outputfilename, *unalignedreadsfilename, *debugfilename, c;
	int *distsarray, *stddevsarray, nd, ns, id, is;
	int argMinIdentity,argMaxErrors,argMaxHits;
	//for(n=0;n<=255;n++){printf("[%.3d:%.2X:'%c']",n,n,n);}printf("\n");
	ConsoleDrawBoxChar("TL",0);
	n=(27+(int)strlen(VERSION));
	ConsoleDrawBoxChar("T",n);
	ConsoleDrawBoxChar("TR",0);
	printf("\n");
	ConsoleDrawBoxChar("L",0);
	ConsoleSetTextColor(COLOR_BRIGHT, COLOR_RED, COLOR_WHITE); // bright, red font, white background
	printf(" [Splice][T][A][P][y][R] v%s ",VERSION);
	ConsoleResetTextColor(); // reset colors
	ConsoleDrawBoxChar("R",0);
	#ifdef DEBUG
	printf(" DEBUG");
	#endif
	printf(" (%s %s)\n",__DATE__,__TIME__);
	ConsoleDrawBoxChar("BL",0);
	n=(27+(int)strlen(VERSION));
	ConsoleDrawBoxChar("B",n);
	ConsoleDrawBoxChar("BR",0);
	printf("\n");
	if(argc<2){
		printf("> USAGE:\n");
		printf("\nBuild index:\n");
		ConsoleSetTextColor(COLOR_BRIGHT, COLOR_WHITE, COLOR_BLACK);
		printf("\t%s I <ref-file>\n",argv[0]);
		ConsoleResetTextColor();
		printf("\t<ref-file>\treference genome to index ( required ; format=FASTA )\n");
		printf("\tOutput:\t\t'<ref-file>.fmi' (FMI format)\n");
		printf("\nMap reads:\n");
		ConsoleSetTextColor(COLOR_BRIGHT, COLOR_WHITE, COLOR_BLACK);
		printf("\t%s <index-file> <reads-files> <options>\n",argv[0]);
		ConsoleResetTextColor();
		printf("\t<index-file>\tindexed reference genome file ( required ; format=FMI )\n");
		printf("\t<reads-files>\tfile(s) with reads to map ( required ; format=FASTA,FASTQ,SFF )\n");
		printf("\t-i\t\tminimum Identity percentage allowed per read ( optional ; default=90 )\n");
		printf("\t-e\t\tmaximum number of Errors allowed per read ( optional ; default=none )\n");
		printf("\t-h\t\tmaximum number of Hits allowed per seed ( optional ; default=64 )\n");
		printf("\t-b\t\treport only the Best read hit instead of all hits ( optional ; default=off )\n");
		printf("\t-o\t\tuse this Output file name instead of default ( optional ; default=none )\n");
		printf("\t-u\t\tsave Unmapped reads to file ( optional ; default=off )\n");
		printf("\tOutput:\t\t'<reads-file>.sam' (SAM format)\n");
		printf("\nOptions for Mate-Pair / Paired-End reads:\n");
		//printf("\t<reads-file1>\tfile with 1st/left mates of the read pairs ( required ; format=FASTA,FASTQ,SFF )\n");
		//printf("\t<reads-file2>\tfile with 2nd/right mates of the read pairs ( optional ; format=FASTA,FASTQ,SFF )\n");
		printf("\t-p\t\tforce mapping in Paired-end reads mode if not auto-detected ( optional ; default=off )\n");
		printf("\t-d\t\taverage Distance allowed between both mate reads ( optional ; default=none )\n");
		printf("\t-sd\t\tdistance Standard Deviation between both mate reads ( optional ; default=none )\n");
		printf("\t-ss\t\tboth mate reads must be on the Same Strand ( optional ; default=none )\n");
		printf("\t-os\t\tboth mate reads must be on Opposite Strands ( optional ; default=none )\n");
		printf("\nOptions for Transcriptomics:\n");
		//printf("\t-r\t\tenable splice-aware mapping for Rna-seq reads ( required ; default=on )\n");
		printf("\t-min\t\tMINimum intron length ( optional ; default=50 )\n");
		printf("\t-max\t\tMAXimum intron length ( optional ; default=500000 )\n");
		printf("\nOptions for Metagenomics:\n");
		printf("\t-s\t\tonly output reads that are mapped to one Single reference ( optional ; default=off )\n");
		printf("\n( Press any key to continue )");
		getchar();
		printf("\nGenerate consensus:\n");
		ConsoleSetTextColor(COLOR_BRIGHT, COLOR_WHITE, COLOR_BLACK);
		printf("\t%s C <ref-file> <sam-file> <options>\n",argv[0]);
		ConsoleResetTextColor();
		printf("\t<ref-file>\treference genome file ( required ; format=FASTA )\n");
		printf("\t<sam-file>\tmapped reads file ( required ; format=SAM )\n");
		printf("\tOptions:\tRun '%s C' without arguments for full list of options\n",argv[0]);
		printf("\tOutput: \tFASTA , AGP , VCF , ACE , CSV\n");
		printf("\nVisualize alignment:\n");
		ConsoleSetTextColor(COLOR_BRIGHT, COLOR_WHITE, COLOR_BLACK);
		printf("\t%s V <ref-file> <sam-file> (<gff-file>) (<start-pos> <end-pos>)\n",argv[0]);
		ConsoleResetTextColor();
		printf("\t<ref-file>\treference genome file ( required ; format=FASTA )\n");
		printf("\t<sam-file>\tmapped reads file ( required ; format=SAM )\n");
		printf("\t<gff-file>\tannotations file ( optional ; format=GFF,GTF )\n");
		printf("\t<start-pos>\tstarting reference position ( optional ; default=0 )\n");
		printf("\t<end-pos>\tending reference position ( optional ; default=inf )\n");
		printf("\tOutput:\t'<sam-file>.bmp' (BMP format)\n");
		printf("\n> EXAMPLE:\n");
		ConsoleSetTextColor(COLOR_BRIGHT, COLOR_WHITE, COLOR_BLACK);
		printf("\t%s i genome.fasta\n",argv[0]);
		printf("\t%s genome.fmi reads.fasta\n",argv[0]);
		printf("\t%s c genome.fasta reads.sam\n",argv[0]);
		printf("\t%s v genome.fasta reads.sam\n",argv[0]);
		ConsoleResetTextColor();
		printf("\n> TOOLS:\n");
		printf("  Report SAM mapping statistics     : ");
		ConsoleSetTextColor(COLOR_BRIGHT, COLOR_WHITE, COLOR_BLACK);
		printf("\t%s S <reference-file> <sam-file> (<num-reads>)\n", argv[0]);
		ConsoleResetTextColor();
		printf("  Convert SFF file to FASTQ file    : ");
		ConsoleSetTextColor(COLOR_BRIGHT, COLOR_WHITE, COLOR_BLACK);
		printf("\t%s F <sff-file>\n",argv[0]);
		ConsoleResetTextColor();
		printf("  Split paired-ends SFF file        : ");
		ConsoleSetTextColor(COLOR_BRIGHT, COLOR_WHITE, COLOR_BLACK);
		printf("\t%s P <sff-file>\n",argv[0]);
		ConsoleResetTextColor();
		printf("  Convert FASTQ file to FASTA file  : ");		
		ConsoleSetTextColor(COLOR_BRIGHT, COLOR_WHITE, COLOR_BLACK);
		printf("\t%s Q <fastq-file>\n", argv[0]);
		ConsoleResetTextColor();
		/*
		printf("  Clean and join multiple sequences : ");
		printf("\t%s A <fasta-file>\n",argv[0]);
		printf("  Convert SAM to sorted CSV format  : ");
		printf("\t%s T <sam-file>\n",argv[0]);
		printf("  Evaluate 'wgsim' reads in SAM     : ");
		printf("\t%s E <sam-file>\n",argv[0]);
		printf("  Calculate rough coverage percentage of mapped reads: \n");
		printf("\t%s X <positions-file> <genome-size>\n",argv[0]);
		printf("  Convert Newbler file ('454ReadStatus.txt') to simple positions file: \n");
		printf("\t%s N <newbler-file>\n",argv[0]);
		printf("  Convert SAM file (BWA) to simple positions file: \n");
		printf("\t%s S <sam-file> <total-num-reads>\n",argv[0]);
		printf("  Convert MAP file (Segemehl) to simple positions file: \n");
		printf("\t%s M <map-file> <total-num-reads>\n",argv[0]);
		*/
		printf("\n> Done!\n");
		exit(0);
	}

	if(( !strcmp(argv[1],"Q") || !strcmp(argv[1],"q") ) && argc==3 ) ConvertFastaQ(argv[2]);
	if(( !strcmp(argv[1],"F") || !strcmp(argv[1],"f") ) && argc==3 ) ConvertSFF(argv[2]);
	if(( !strcmp(argv[1],"P") || !strcmp(argv[1],"p") ) && argc==3 ) SplitPairedEndSff(argv[2]);
	if(( !strcmp(argv[1],"A") || !strcmp(argv[1],"a") ) && argc==3 ) CleanFasta(argv[2]);
	/*
	if(( !strcmp(argv[1],"T") || !strcmp(argv[1],"t") ) && argc==3 ) ConvertSAMToCSV(argv[2]);
	if(( !strcmp(argv[1],"X") || !strcmp(argv[1],"x") ) && argc==4 ) CalculateCoverage(argv[2],atoi(argv[3]));
	if(( !strcmp(argv[1],"N") || !strcmp(argv[1],"n") ) && argc==3 ) Process454File(argv[2]);
	if(( !strcmp(argv[1],"S") || !strcmp(argv[1],"s") ) && argc==4 ) ProcessSamFile(argv[2],atoi(argv[3]));
	if(( !strcmp(argv[1],"M") || !strcmp(argv[1],"m") ) && argc==4 ) ProcessMapFile(argv[2],atoi(argv[3]));
	if(( !strcmp(argv[1],"O") || !strcmp(argv[1],"o") ) && argc==4 ) ProcessSoapFile(argv[2],atoi(argv[3]));
	if(( !strcmp(argv[1],"B") || !strcmp(argv[1],"b") ) && argc==4 ) ProcessBowtieFile(argv[2],atoi(argv[3]));
	*/

	if(( !strcmp(argv[1],"C") || !strcmp(argv[1],"c") ) && argc>=2 ) GenerateConsensus(argv,argc);
	if(( !strcmp(argv[1],"V") || !strcmp(argv[1],"v") ) && (argc==4 || argc==6 || argc==7) ){
		if (argc < 4 || argc > 7) {
			printf("> ERROR: Invalid number of arguments to create visualization\n");
			exit(-1);
		}
		if ((argc == 6 || argc == 7) && (atoi(argv[argc - 1]) < 0 || atoi(argv[argc - 2]) <= 0)) {
			printf("> ERROR: Invalid range positions for visualization\n");
			exit(-1);
		}
		if(argc==4) DrawMappingPlot(argv[2],argv[3],NULL,0,0); // plot for entire sequence
		else if(argc==5) DrawMappingPlot(argv[2],argv[3],argv[4],0,0); // also use GFF
		else if(argc==6) DrawMappingPlot(argv[2],argv[3],NULL,atoi(argv[4]),atoi(argv[5])); // plot for specified range
		else if(argc==7) DrawMappingPlot(argv[2],argv[3],argv[4],atoi(argv[5]),atoi(argv[6])); // draw GFF annotations too
		exit(0);
	}
	if(( !strcmp(argv[1],"S") || !strcmp(argv[1],"s") ) && (argc==4 || argc==5 || argc==6) ){
		if(argc==4) ReportSAMStatistics(argv[2],argv[3],0,0); // ref + sam
		else {
			k=(unsigned int)strtoul(argv[4],NULL,10);
			if(argc==6) ReportSAMStatistics(argv[2],argv[3],k,atoi(argv[5])); // ref + sam + numreads + savebadreads
			else { // if(argc==5)
				if(k>1) ReportSAMStatistics(argv[2],argv[3],k,0); // ref + sam + numreads
				else ReportSAMStatistics(argv[2],argv[3],0,1); // ref + sam + savebadreads
			}
		}
	}
	if(( !strcmp(argv[1],"E") || !strcmp(argv[1],"e") ) && (argc==3 || argc==4 ) ){
		if(argc==3) EvaluateSAMPositions(argv[2],0);
		else EvaluateSAMPositions(argv[2],1);
	}
	if(( !strcmp(argv[1],"R") || !strcmp(argv[1],"r") ) && argc==5 ) GetRealErrorStatistics(argv[2],argv[3],argv[4]);

	if( !strcmp(argv[1],"I") || !strcmp(argv[1],"i") ){ // build index
		indexfilename=AppendToBasename(argv[2],".fmi");
		FMI_BuildIndex(argv[2],0);
		FMI_SaveIndex(indexfilename);
		FMI_FreeIndex();
		free(indexfilename);
		printf("> Done!\n");
		#ifdef PAUSEATEXIT
		getchar();
		#endif
		return 0;
	}

	// parse command line arguments
	if( ( argMinIdentity = ParseArgument(argc,argv,"I",1) ) < 0 ) argMinIdentity=90; // default is 90% identity
	if( ( argMaxErrors = ParseArgument(argc,argv,"E",1) ) < 0 ) argMaxErrors=(-1); // this can be set to zero
	if( ( argMaxHits = ParseArgument(argc,argv,"H",1) ) < 0 ) argMaxHits=64;

	if(argMinIdentity>100) argMinIdentity=100;
	maxReadPairDistance = ParseArgument(argc,argv,"D",0); // initialize this here because it will be needed by LoadReadsFiles, but it will be set correctly later
	reportUnalignedReads = ParseArgument(argc,argv,"U",0);
	reportAllHits = ( 1 - ParseArgument(argc,argv,"B",0) ); // opposite of report the best hit only
	singleRefReadsOnly = ParseArgument(argc,argv,"S",0);
	
	pairedEndReadsMode = ParseArgument(argc,argv,"P",0);
	//rnaSequencingMode = ParseArgument(argc,argv,"R",0);
	// Enable RNA Sequencing Mode by default
	rnaSequencingMode = 1;
	
	if(rnaSequencingMode){ // get rna-seq options and load defaults if needed
		if( ( minIntronLength = ParseArgument(argc,argv,"MI",1) ) < 0 ) minIntronLength=50;
		if( ( maxIntronLength = ParseArgument(argc,argv,"MA",1) ) < 0 ) maxIntronLength=500000;
	}

	indexfilename=NULL;
	genomefilename=NULL;
	unalignedreadsfilename=NULL;
	errorsFile=NULL;
	debugfile=NULL;
	debugfilename=NULL;
	readsFiles=NULL;
	readsFilenames=NULL;
	readsFile=NULL;
	readsFilename=NULL;

	n=1;
	while(n<argc){ // get index filename
		if(argv[n][0]=='-'){
			c=argv[n][1];
			if(c=='p'||c=='r'||c=='s'||c=='b'||c=='u'||c=='P'||c=='R'||c=='S'||c=='B'||c=='U') n++; // skip options with no content (TODO: missing "SS" and "OS")
			else n+=2; // skip options name and its content
			continue;
		}
		indexfilename=argv[n];
		break;
	}
	if(indexfilename==NULL){
		printf("> ERROR: Missing index filename argument\n");
		exit(-1);
	}
	n++;
	numReadsFiles=0;
	readsFilenames=(char **)malloc(1*sizeof(char *));
	while(n<argc){ // get all reads filenames
		if(argv[n][0]=='-'){
			c=argv[n][1];
			if(c=='p'||c=='r'||c=='s'||c=='b'||c=='u'||c=='P'||c=='R'||c=='S'||c=='B'||c=='U') n++; // skip options with no content (TODO: missing "SS" and "OS")
			else n+=2; // skip options name and its content
			continue;
		}
		if(numReadsFiles>0) readsFilenames=(char **)realloc(readsFilenames,(numReadsFiles+1)*sizeof(char *));
		readsFilenames[numReadsFiles]=argv[n];
		numReadsFiles++;
		n++;
	}
	if(numReadsFiles==0){
		printf("> ERROR: Missing reads filename argument\n");
		exit(-1);
	}
	
	FMI_LoadIndex(indexfilename);
	sizeRefGenome = FMI_GetTextSize(); // get size from the index to later compare it to the size of the un-indexed reference
	sizeBWT = FMI_GetBWTSize(); // size of the BWT to later initialize the pointers in the index search
	
	LoadReferenceGenome(indexfilename);
	LoadReadsFiles();

	nd=0;
	ns=0;
	distsarray=ParseMultiArgument(argc,argv,"D",&nd); // get arrays with all input values
	stddevsarray=ParseMultiArgument(argc,argv,"SD",&ns);
	if(pairedEndReadsMode){ // distribute distances through paired-end files
		id=0;
		is=0;
		for(n=0;n<numReadsFiles;n++){
			c=readsFileTypes[n];
			minDists[n]=-1; // set default for no distance or single reads files
			maxDists[n]=-1;
			if((c=='f' || c=='q' || c=='s') && (id<nd)){ // paired-end files
				minDists[n]=distsarray[id]; // set average distance
				maxDists[n]=distsarray[id];
				if(is<ns){ // allow standard deviation for both sides
					minDists[n]-=stddevsarray[is];
					maxDists[n]+=stddevsarray[is];
				} else { // set stddev default to 25% of the distance
					minDists[n]-=(distsarray[id]/4);
					maxDists[n]+=(distsarray[id]/4);
				}
				if(minDists[n]<0) minDists[n]=0;
				if(c=='f' || c=='q'){ // set distances for 2nd pair too
					minDists[(n+1)]=minDists[n];
					maxDists[(n+1)]=maxDists[n];
					n++; // skip 2nd pair
				}
				id++; // advance in arrays
				is++;
			}
		}
		reqReadPairStrands = -1; // get strand requirements too
		if(ParseArgument(argc,argv,"SS",0)) reqReadPairStrands = 1; // same strand
		if(ParseArgument(argc,argv,"OS",0)) reqReadPairStrands = 2; // opposite strands
	}

	if( pairedEndReadsMode && numReadsFiles==2 && readsFileTypes[0]!='s' ){ // remove '1'/'2' chars from output filenames
		n=(int)strlen(readsFilenames[0]); // size of 1st read file name
		while(n>0 && readsFilenames[0][n]!='1') n--; // rightmost position of char '1'
		if(n!=0 && n<(int)strlen(readsFilenames[1]) && readsFilenames[1][n]=='2'){ // if a char '2' exists in the same position of 2nd reads file
			c=readsFilenames[0][(n-1)]; // char in the position before
			if( (c>=48 && c<=57) || (c>=65 && c<=90) || (c>=97 && c<=122) ) // if it has an alphanumeric char before
				while(readsFilenames[0][n]!='\0'){
					readsFilenames[0][n]=readsFilenames[0][(n+1)]; // delete '1' char by shifting all right chars 1 pos to the left
					n++;
				}
			else // non-alphanumeric char before number
				while(readsFilenames[0][n]!='\0'){
					readsFilenames[0][(n-1)]=readsFilenames[0][(n+1)]; // delete both non-alphanumeric and '1' chars
					n++;
			}
		}
	}
	n=ParseArgument(argc,argv,"O",2); // index in the arguments array of the alternative output filename if it exists
	if(n!=(-1)) outputfilename=AppendToBasename(argv[n],".sam"); // override default
	else outputfilename=AppendToBasename(readsFilenames[0],".sam"); // default is the reads filename
	outputFile=fopen(outputfilename,"w");
	if(outputFile==NULL) {
		printf("> ERROR: Can't write output file\n");
		exit(-1);
	}

	if( reportUnalignedReads ){
		unalignedreadsfilename=AppendToBasename(outputfilename,"-unaligned.txt");
		unalignedReadsFile=fopen(unalignedreadsfilename,"w");
		if(unalignedReadsFile==NULL) {
			printf("> ERROR: Can't write unaligned reads file\n");
			exit(-1);
		}
	}

	#if defined DEBUG || defined DEBUGDP
	debugfilename=AppendToBasename(readsFilenames[0],"-debug.txt");
	debugfile=fopen(debugfilename,"w");
	if(debugfile==NULL) {
		printf("> ERROR: Can't write debug file\n");
		exit(-1);
	}
	fprintf(debugfile,">ReadName\nNumReads ReadSize MaxReadErrors\n");
	fprintf(debugfile,"Strand\nReadChars\n");
	fprintf(debugfile,"SeedId+1 [SeedStartPosInRead-SeedEndPosInRead] SeedSize NumSeedHits \tSeedChars\n...\n");
	fprintf(debugfile,"Strand\n");
	fprintf(debugfile,"NumSeedHits [SortedHitId+1] SeedId+1 : SeedPosition (SeedSize) ...\n");
	fprintf(debugfile,"NumSeedChains [ChainStart+1] (ChainSize) ChainScore, ...\n...\n");
	fprintf(debugfile,"NumBestSeedChains Strand [ChainStart+1] (ChainSize) ChainScore, ...\n");
	fprintf(debugfile,"[Seed:SeedId/NumSeeds](ChainSeedId/ChainSize) SeedSize\n");
	fprintf(debugfile,"SeedCharsInGenome\nSeedCharsInRead\n");
	fprintf(debugfile,"{gapStartPosInGenome-GapEndPosInGenome} (GapSizeInGenome)\n");
	fprintf(debugfile,"{gapStartPosInRead-GapEndPosInRead} (GapSizeInRead)\n");
	fprintf(debugfile,"NumDPErrors DPTraceSize\n...\n");
	fprintf(debugfile,"Strand ReadPosInGenome NumReadErrors\n\n");
	#endif
	
	if( pairedEndReadsMode || rnaSequencingMode){
		printf("> Using mode: ");
		if(rnaSequencingMode) printf("RNASequencingMode");
		else printf("DNASequencingMode");
		if(pairedEndReadsMode) printf(" + PairedEndReadsMode");
		printf("\n");
	}
	printf("> Using parameters: ");
	if(argMaxErrors==-1){ // number of errors not defined, so, use identity percentage
		printf("MinimumIdentityPercentage=%d",argMinIdentity);
		minIdentityPct=argMinIdentity;
		maxNumErrors=0;
		maxSeedHits=argMaxHits;
	} else {
		printf("MaximumNumberErrors=%d",argMaxErrors);
		minIdentityPct=0;
		maxNumErrors=argMaxErrors;
		maxSeedHits=argMaxHits;
	}
	printf(" ; MaximumSeedHits=%d",argMaxHits);
	if(rnaSequencingMode){
		printf(" ; MinimumIntronLength=%d",minIntronLength);
		printf(" ; MaximumIntronLength=%d",maxIntronLength);
	}
	if(pairedEndReadsMode){
		printf(" ; MinimumReadPairDistance=");
		for(n=0;n<numReadsFiles;n++){
			if(n!=0) printf("|");
			c=readsFileTypes[n];
			if(c=='F'||c=='Q'||c=='S') printf("n/a");
			else if(minDists[n]==-1) printf("none");
			else printf("%d",minDists[n]);
		}
		printf(" ; MaximumReadPairDistance=");
		for(n=0;n<numReadsFiles;n++){
			if(n!=0) printf("|");
			c=readsFileTypes[n];
			if(c=='F'||c=='Q'||c=='S') printf("n/a");
			else if(maxDists[n]==-1) printf("none");
			else printf("%d",maxDists[n]);
		}
		printf(" ; ReadPairStrand=");
		if(reqReadPairStrands==1) printf("same");
		else if(reqReadPairStrands==1) printf("opposite");
		else printf("none");
	}
	printf("\n");

	#ifdef __unix__
	PrintProgressBar(0.0,0); // initialize progress bar
	#endif
	printf("> Aligning reads... 0.00%%");
	fflush(stdout);
	ftime(&starttb);
	RunSeededAlignment();
	ftime(&endtb);
	timetb=((endtb.time) + (endtb.millitm)/1000.0) - ((starttb.time) + (starttb.millitm)/1000.0);
	FMI_FreeIndex();
	#ifdef __unix__
	PrintProgressBar(100.0,1); // print 100% progress
	#endif
	printf("\n> Aligning reads... 100.00%% (");
	PrintUnsignedNumber(numReads);printf(" processed : ");
	PrintUnsignedNumber(numAlignedReads);printf(" aligned + ");
	PrintUnsignedNumber(numUnalignedReads);printf(" unaligned) completed!\n");
	printf("> Writing mapped reads file <%s>... (",outputfilename);
	PrintNumber((long long int)ftell(outputFile));
	printf(" bytes) ");
	fflush(stdout);
	for(n=0;n<numReadsFiles;n++) fclose(readsFiles[n]);
	fclose(outputFile);
	printf("OK\n");

	if(unalignedReadsFile!=NULL){
		printf("> Writing un-mapped reads file <%s>... (",unalignedreadsfilename);
		PrintNumber((long long int)ftell(unalignedReadsFile));
		printf(" bytes) ");
		fflush(stdout);
		fclose(unalignedReadsFile);
		free(unalignedreadsfilename);
		printf("OK\n");
	}

	if(errorsFile!=NULL){
		printf("> WARNING: Errors in CIGAR strings alignments were found and written to file <%s>... (","errors.txt");
		PrintNumber((long long int)ftell(errorsFile));
		printf(" bytes) ");
		fflush(stdout);
		fclose(errorsFile);
		printf("OK\n");
	}

	#if defined DEBUG || defined DEBUGDP
	printf("> Writing debug file <%s>... (",debugfilename);
	PrintNumber((long long int)ftell(debugfile));
	printf(" bytes) ");
	fclose(debugfile);
	free(debugfilename);
	printf("OK\n");
	#endif

	for(n=0;n<numReadsFiles;n++){ // check if at least one of the reads files has paired ends to print info next (if last file was not paired-ends, the variable was set to 0)
		c=readsFileTypes[n];
		if(c=='f' || c=='q' || c=='s') pairedEndReadsMode=1;
	}

	printf("> Running time  : ");PrintTime(timetb);printf("(%.3lf reads/s)\n",(double)((double)numReads/timetb));
	printf("> Aligned reads : ");PrintUnsignedNumber(numAlignedReads);printf(" of ");PrintUnsignedNumber(numReads);
	printf(" (%.3lf%%)\n",(numReads==0)?(0):(double)((((double)numAlignedReads)*100.0)/((double)numReads)));
	printf("> Read counts   : ");
	PrintUnsignedNumber(numReads);printf(" total reads (%.2lfx coverage) , ",(double)((double)totalNumReadsChars/(double)sizeRefGenome));
	PrintUnsignedNumber(numUnalignedReads);printf(" unaligned");
	//#ifdef DEBUG
	if (numAlignedHardReads != 0) { printf(" ("); PrintUnsignedNumber(numAlignedHardReads); printf(" aligned in 'hard mode')"); }
	if (numHighRepetitionReads != 0) { printf(" ("); PrintUnsignedNumber(numHighRepetitionReads); printf(" with too many hits)"); }
	//#endif
	printf("\n");
	printf("> Read sizes    : ");
	PrintUnsignedNumber((numReads == 0) ? (0) : (unsigned int)(totalNumReadsChars / (long long unsigned int)numReads)); printf(" bp avg size , ");
	PrintUnsignedNumber(maxReadSize); printf(" bp maximum size\n");
	printf("> Read hits     : %.3lf hits/read , %u maximum hits\n",(numAlignedReads==0)?(0.0):(double)(((double)totalNumReadHits)/((double)numAlignedReads)),maxReadHits);
	#ifdef DEBUG
	printf("> Seed stats    : %.3f seeds/read , %.3f occurrences/seed , %.3f chains/read\n",(numAlignedReads==0)?(0.0):(float)(((float)totalNumSeeds)/((float)numAlignedReads)),(float)(((float)totalNumSeedOccs)/((float)totalNumSeeds)),(numAlignedReads==0)?(0.0):(float)(((float)totalNumSeedChains)/((float)numAlignedReads)));	
	#endif
	printf("> Error stats   : ");
	printf("%.3lf errors/read , ",(numAlignedReads == 0) ? (0.0) : (double)(((double)totalNumDPErrors) / ((double)numAlignedReads)));
	printf("%.3lf errors/bp ", (totalNumMappedReadsChars == 0) ? (0.0) : (double)(((double)totalNumDPErrors) / ((double)totalNumMappedReadsChars)));
	printf("(1 error per "); PrintUnsignedNumber((totalNumDPErrors == 0) ? (0) : (unsigned int)(totalNumMappedReadsChars / totalNumDPErrors)); printf(" bp)\n");
	if( rnaSequencingMode ){
		printf("> RNA-Seq stats : ");
		PrintUnsignedNumber(numSplicedReads); printf(" spliced reads ; ");
		PrintUnsignedNumber(numSplicingEvents); printf(" splicing events detected\n");
		#ifdef DEBUG
		printf("> SS Seeds      : "); PrintUnsignedNumber(numSeedsAddedBySpliceSites); printf(" extra seeds added because of detected splice sites\n");
		PrintSSBlocksStats();
		PrintSpliceSites(readsFilenames[0]);
		#endif
		FreeSpliceSitesArray();
	}
	if (pairedEndReadsMode) {
		printf("> Pair counts   : ");
		PrintUnsignedNumber((numReads / 2)); printf(" total pairs , ");
		PrintUnsignedNumber(numAlignedBothMates); printf(" with both mates aligned , ");
		PrintUnsignedNumber(numAlignedSingleMates); printf(" with only one mate aligned\n");
		printf("> Pair distances: ");
		PrintUnsignedNumber(avgObservedPairDistances); printf(" bp avg , ");
		PrintUnsignedNumber(minObservedPairDistance); printf(" bp minimum , ");
		PrintUnsignedNumber(maxObservedPairDistance); printf(" bp maximum\n");
	}
	printf("> Done!\n");
	fflush(stdout);
	/*
	#ifdef DEBUG
	ReportSAMStatistics(genomefilename,outputfilename,0,0);
	#endif
	*/
	free(readsFiles);
	free(readsFileTypes);
	free(readsFilenames);
	free(minDists);
	free(maxDists);
	free(genomefilename);
	free(outputfilename);
	free(refGenome);
	for(n=0;n<numRefs;n++) free(refsNames[n]);
	free(refsNames);
	free(refsEndPositions);
	#ifdef PAUSEATEXIT
	getchar();
	#endif
	return 0;
}
