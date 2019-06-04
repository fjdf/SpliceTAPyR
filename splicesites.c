#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "splicesites.h"
#ifdef DEBUG
#include "alignreads.h"
#include "tools.h"
#endif

typedef struct _SpliceSite {
	unsigned int position;		// last (or first) position of the exon before (or after) the intron
	unsigned int oppositePosition;	// first (or last) position of the exon after (or before) the intron
	unsigned int nextIdInBlock; // next splice site in this reference block in the array of all splice sites
	unsigned int count; // number of times that this splice site was detected
} SpliceSite;


static SpliceSite *spliceSitesArray = NULL; // array of all splice sites, indexed by id; position 0 of this array is never used
static unsigned int *genomeSpliceSiteBlocks = NULL; // contains the id of the first splice site found on that block's range of positions
static unsigned int lastGenomePos = 0;
static unsigned int numSpliceSites = 0;
static unsigned int maxNumSpliceSites = 0;
static unsigned int numUniqueSpliceSites = 0;
static const unsigned int ssBlocksShift = 15; // blocks of 2^15 = 32.768 positions each
static const unsigned int numMoreSSIdsToAdd = 1024; // number of new splice site ids added to the array each time it is expanded
static const unsigned int checkLengthAroundSS = 5; // when looking for splice site positions, a position is accepted if it is this many position to the right or to the left of the actual position

#ifdef DEBUG
static unsigned int numExtraSSIds = 0;
static unsigned int numSSBlocks = 0;
static unsigned int sizeSSBlocks = 0;
static unsigned int numFilledSSBlocks = 0;
static unsigned int numExtendedSSBlocks = 0;
static unsigned int maxSSBlockExtention = 0;
static unsigned int numAlternativeSplicingCollisions = 0;
#endif
  
void InitializeSpliceSitesArray(unsigned int genomeSize){
	unsigned int numBlocks, blocksMask;
	blocksMask = (unsigned int)((2U << ssBlocksShift) - 1U); // to get the remainder of the division by (2^ssBlocksShift)
	numBlocks = (genomeSize >> ssBlocksShift); // number of blocks of size (2^ssBlocksShift)
	if ((genomeSize & blocksMask) != 0U) numBlocks++; // if there are positions left, add an extra final block
	genomeSpliceSiteBlocks = (unsigned int *)calloc(numBlocks, sizeof(unsigned int));
	spliceSitesArray = (SpliceSite *)malloc((1 + numMoreSSIdsToAdd) * sizeof(SpliceSite));
	maxNumSpliceSites = numMoreSSIdsToAdd;
	lastGenomePos = (genomeSize - 1);
	numSpliceSites = 0;
	numUniqueSpliceSites = 0;
	#ifdef DEBUG
	numSSBlocks = numBlocks;
	sizeSSBlocks = (genomeSize / numBlocks);
	#endif
}

void FreeSpliceSitesArray() {
	if (spliceSitesArray != NULL) free(spliceSitesArray);
	if (genomeSpliceSiteBlocks != NULL) free(genomeSpliceSiteBlocks);
}

// TODO: if alternative splicing exists, one position can have multiple opposite positions
// Checks if a Splice Site exists in the interval of [(position - checkLengthAroundSS),(position + checkLengthAroundSS)]
// and if yes, returns the opposite position of that Splice Site
// NOTE: The input argument "position" is modified to contain the real Splice Site position
unsigned int GetOppositePositionIfSpliceSiteExists(unsigned int *position) {
	unsigned int blockId, nextBlockId, ssId, ssPos, lowPos, highPos;
	lowPos = ((*position) - checkLengthAroundSS);
	highPos = ((*position) + checkLengthAroundSS);
	if (lowPos >= lastGenomePos) lowPos = 0; // prevent underflow
	if (highPos >= lastGenomePos) highPos = lastGenomePos; // prevent going further that the size of the genome
	blockId = (lowPos >> ssBlocksShift); // by summing or subtracting positions, we can fall on a different block
	nextBlockId = (highPos >> ssBlocksShift);
	ssId = genomeSpliceSiteBlocks[blockId];
	while (1) {
		if (ssId == 0) { // empty block or we reached the end of an extended block
			if (blockId == nextBlockId) break; // no new block to check
			blockId = nextBlockId; // check next block
			ssId = genomeSpliceSiteBlocks[blockId];
			continue;
		}
		ssPos = spliceSitesArray[ssId].position;
		if (ssPos > highPos) break; // if the pos we are looking for is already behind, stop checking
		if (ssPos < lowPos) { // if the pos we are looking for is still ahead, continue checking
			ssId = spliceSitesArray[ssId].nextIdInBlock;
			continue;
		} // if ((ssPos >= lowPos) && (ssPos <= highPos))
		(*position) = ssPos; // position found in interval
		return spliceSitesArray[ssId].oppositePosition;
	}
	return 0;
}

unsigned int GetSpliceSiteId(unsigned int position) {
	unsigned int blockId, ssId;
	blockId = (position >> ssBlocksShift);
	ssId = genomeSpliceSiteBlocks[blockId];
	while (ssId != 0) {
		if (spliceSitesArray[ssId].position == position) return ssId;
		ssId = spliceSitesArray[ssId].nextIdInBlock;
	}
	return 0;
}

// TODO: if alternative splicing exists, one position can have multiple opposite positions
// Adds the Splice Site to the array of Splices Sites if it is new, or updates its count if it already exists
unsigned int AddOrUpdateSpliceSite(unsigned int leftPosition, unsigned int rightPosition) {
	unsigned int pos, opPos, ssPos, blockId, prevSsId, nextSsId, newSsId;
	#ifdef DEBUG
	unsigned int extCount = 0;
	#endif
	newSsId = 0;
	pos = leftPosition; // first add the left position
	opPos = rightPosition;
	while (1) {
		ssPos = 0; // position of the current splice site being checked
		blockId = (pos >> ssBlocksShift);
		prevSsId = 0;
		nextSsId = genomeSpliceSiteBlocks[blockId];
		while ((nextSsId != 0) && ((ssPos = spliceSitesArray[nextSsId].position) < pos)) {
			prevSsId = nextSsId;
			nextSsId = spliceSitesArray[nextSsId].nextIdInBlock;
			#ifdef DEBUG
			extCount++;
			#endif
		}
		if (nextSsId != 0 && ssPos == pos) { // splice site already exists
			spliceSitesArray[nextSsId].count++;
			#ifdef DEBUG
			if (spliceSitesArray[nextSsId].oppositePosition != opPos) numAlternativeSplicingCollisions++;
			#endif
		} else { // create new splice site
			#ifdef DEBUG
			if (extCount > maxSSBlockExtention) maxSSBlockExtention = extCount;
			#endif
			numSpliceSites++;
			if (numSpliceSites == maxNumSpliceSites) {
				maxNumSpliceSites += numMoreSSIdsToAdd;
				spliceSitesArray = (SpliceSite *)realloc(spliceSitesArray, maxNumSpliceSites * sizeof(SpliceSite));
			}
			newSsId = numSpliceSites;
			spliceSitesArray[newSsId].position = pos;
			spliceSitesArray[newSsId].oppositePosition = opPos;
			spliceSitesArray[newSsId].nextIdInBlock = nextSsId;
			spliceSitesArray[newSsId].count = 1;
			if (prevSsId == 0) {
				genomeSpliceSiteBlocks[blockId] = newSsId;
				#ifdef DEBUG
				numFilledSSBlocks++;
				#endif
			}
			else {
				spliceSitesArray[prevSsId].nextIdInBlock = newSsId;
				#ifdef DEBUG
				numExtraSSIds++;
				if (genomeSpliceSiteBlocks[blockId] == prevSsId) numExtendedSSBlocks++;
				#endif
			}
		} // new splice site
		if (pos == rightPosition) break;
		else { // now do the same to add the right position
			pos = rightPosition;
			opPos = leftPosition;
		}
	} // end of loop for both positions
	if (newSsId != 0) {
		numUniqueSpliceSites++;
		#ifdef DEBUG
		fprintf(debugfile,"![ %u - %u ]\n",leftPosition,rightPosition);
		#endif
	}
	return newSsId;
}

#ifdef DEBUG
void PrintSSBlocksStats() {
	unsigned int i, id, pos, prevpos, dist, mindist, maxdist, numss, numsspblock;
	long long unsigned int avgdist;
	unsigned int avgsspblock, maxsspblock;
	avgdist = 0;
	mindist = UINT_MAX;
	maxdist = 0;
	avgsspblock = 0;
	maxsspblock = 0;
	numss = 0;
	prevpos = 0;
	for (i = 0; i < numSSBlocks; i++) {
		if ((id = genomeSpliceSiteBlocks[i]) != 0) {
			numsspblock = 0;
			while (id != 0) {
				pos = spliceSitesArray[id].position;
				if (prevpos != 0) {
					dist = (pos - prevpos);
					if (dist < mindist) mindist = dist;
					if (dist > maxdist) maxdist = dist;
					avgdist += dist;
				}
				id = spliceSitesArray[id].nextIdInBlock;
				prevpos = pos;
				numsspblock++;
				numss++;
			}
			if (numsspblock > maxsspblock) maxsspblock = numsspblock;
			avgsspblock += numsspblock;
		}
	}
	avgsspblock /= numSSBlocks;
	if(numss!=0) avgdist /= numss;
	else avgdist = 0;
	printf("> SS Blocks     : %u blocks ; %u filled ; %u extended ; ", numSSBlocks, numFilledSSBlocks, numExtendedSSBlocks);
	PrintUnsignedNumber((unsigned int)(numSpliceSites * sizeof(spliceSitesArray[0]) + numSSBlocks * sizeof(genomeSpliceSiteBlocks[0]))); printf(" bytes\n");
	printf("> SS Blocks +   : %u bp/block ; avg %u ss/block (max=%u) ; avg %u bp between ss's (min=%u, max=%u)\n", sizeSSBlocks, avgsspblock, maxsspblock, (unsigned int)avgdist, mindist, maxdist);
	printf("> SS IDs        : %u total IDs ; %u unique ; %u extra ; max in block = %u\n", numSpliceSites, numUniqueSpliceSites, numExtraSSIds, maxSSBlockExtention);
	printf("> A.S. coll.s   : %u alternative splicing collisions\n", numAlternativeSplicingCollisions);
}

void PrintSpliceSites(char *readsfilename) {
	FILE *ssfile;
	char *ssfilename;
	unsigned int i, n, id, p, op;
	ssfilename = AppendToBasename(readsfilename, "-splicesites.txt");
	ssfile = fopen(ssfilename, "w");
	if (ssfile == NULL) {
		printf("> ERROR: Can't write splice sites file\n");
		return;
	}
	for (i = 0; i < numSSBlocks; i++) {
		if ((id = genomeSpliceSiteBlocks[i]) != 0) {
			fprintf(ssfile, "{%u}\n", i);
			n = 0;
			while (id != 0) {
				p = spliceSitesArray[id].position;
				op = spliceSitesArray[id].oppositePosition;
				fprintf(ssfile, "[%u] %u %c-%c %u (%u)\n", id, p, (p > op) ? '<' : '-', (p < op) ? '>' : '-', op, spliceSitesArray[id].count);
				id = spliceSitesArray[id].nextIdInBlock;
				n++;
			}
			fprintf(ssfile, "(%u)\n", n);
		}
	}
	printf("> Writing splice sites file <%s>... (", ssfilename);
	PrintNumber((long long int)ftell(ssfile));
	printf(" bytes) ");
	fclose(ssfile);
	free(ssfilename);
	printf("OK\n");
}
#endif


#if ! defined(DEBUG) && ! defined(GDB)
inline
#endif
// Checks for the single canonical splice sites: GT-AG
int CheckSite_1GTAG(char *genome, unsigned int leftPos, unsigned int rightPos) {
	if ( (genome[leftPos]=='G') && (genome[leftPos+1]=='T')
			&& (genome[rightPos]=='A') && (genome[rightPos+1]=='G') ) return 1;
	else return 0;
}

#if ! defined(DEBUG) && ! defined(GDB)
inline
#endif
// Checks for the 3 canonical splice sites: GT-AG , GC-AG , AT-AC
int CheckSite_3Canonical(char *genome, unsigned int leftPos, unsigned int rightPos) {
	if (genome[rightPos] == 'A') { // gt-[A]g , gc-[A]g , at-[A]c
		if (genome[leftPos] == 'G') { // [G]t-Ag , [G]c-Ag
			if (genome[rightPos + 1] == 'G') { // Gt-A[G] , Gc-A[G]
				if ((genome[leftPos + 1] == 'T') || (genome[leftPos + 1] == 'C')) // G[T]-AG , G[C]-AG
					return 1;
			}
			return 0;
		}
		if ((genome[leftPos] == 'A') && (genome[leftPos+1] == 'T') && (genome[rightPos+1] == 'C')) // [A][T]-A[C]
			return 1;
	}
	return 0;
}

#if ! defined(DEBUG) && ! defined(GDB)
inline
#endif
// Checks for the 3 canonical splice sites and their reverse complementaries: GT-AG , GC-AG , AT-AC and CT-AC , CT-GC , GT-AT
int CheckSite_3CanonicalAndRC(char *genome, unsigned int leftPos, unsigned int rightPos) {
	if (genome[rightPos] == 'A') { // gt-[A]g , gc-[A]g , at-[A]c ; ct-[A]c , gt-[A]t
		if (genome[leftPos] == 'G') { // [G]t-Ag , [G]c-Ag ; [G]t-At
			if (genome[rightPos + 1] == 'G') { // Gt-A[G] , Gc-A[G]
				if ((genome[leftPos + 1] == 'T') || (genome[leftPos + 1] == 'C')) // G[T]-AG , G[C]-AG
					return 1;
			} // Gx-A[~G]
			else if ((genome[rightPos + 1] == 'T') && (genome[leftPos + 1] == 'T')) // ; G[T]-A[T]
				return 1;
		} // [~G]x-Ax
		else if ((genome[leftPos + 1] == 'T') && (genome[rightPos + 1] == 'C')) { // a[T]-A[C] ; c[T]-A[C]
			if ((genome[leftPos] == 'A') || (genome[leftPos] == 'C')) // [A]T-AC ; [C]T-AC
				return 1;
		}
	} // xx-[~A]x
	else if ((genome[rightPos] == 'G') && (genome[rightPos + 1] == 'C') && (genome[leftPos] == 'C') && (genome[leftPos + 1] == 'T')) // ; [C][T]-[G][C]
		return 1;
	return 0;
}

#if ! defined(DEBUG) && ! defined(GDB)
inline
#endif
// Checks for the 11 canonical and non-canonical splice sites mentioned in the RUM paper
int CheckSite_11RUM(char *genome, unsigned int leftPos, unsigned int rightPos) {
	if ((genome[leftPos] == 'G') && (genome[leftPos + 1] == 'T')) {
		if ((genome[rightPos] == 'A') && (genome[rightPos + 1] == 'G')) return 1;	// GT-AG
		if ((genome[rightPos] == 'T') && (genome[rightPos + 1] == 'G')) return 1;	// GT-TG
		if ((genome[rightPos] == 'A') && (genome[rightPos + 1] == 'A')) return 1;	// GT-AA
		return 0;
	}
	if ((genome[leftPos] == 'G') && (genome[leftPos + 1] == 'C')) {
		if ((genome[rightPos] == 'A') && (genome[rightPos + 1] == 'G')) return 1;	// GC-AG
		if ((genome[rightPos] == 'T') && (genome[rightPos + 1] == 'G')) return 1;	// GC-TG
		if ((genome[rightPos] == 'A') && (genome[rightPos + 1] == 'A')) return 1;	// GC-AA
		if ((genome[rightPos] == 'G') && (genome[rightPos + 1] == 'G')) return 1;	// GC-GG
		return 0;
	}
	if ((genome[leftPos] == 'A') && (genome[leftPos + 1] == 'T')) {
		if ((genome[rightPos] == 'A') && (genome[rightPos + 1] == 'C')) return 1;	// AT-AC
		if ((genome[rightPos] == 'A') && (genome[rightPos + 1] == 'A')) return 1;	// AT-AA
		if ((genome[rightPos] == 'A') && (genome[rightPos + 1] == 'G')) return 1;	// AT-AG
		if ((genome[rightPos] == 'A') && (genome[rightPos + 1] == 'T')) return 1;	// AT-AT
		return 0;
	}
	return 0;
}


// TODO: rename to "GapHasSpliceSignals()"
// Check if there are splice signals on both ends of the gap between two consecutive seeds
// NOTE: At the beginning, the right seed is fully extended to the left (the next char on the left is a mismatch)
// NOTE: At the beginning, the gap size inside the read is 0, i.e. the seeds are right next to each other on the read
// NOTE: The right end of the left seed and the left end of the right seed may be shifted to the left or right when looking for splice signals
// NOTE: Returns INT_MIN if no splicing is found, and the (signed) number of shifted positions (might be 0) if splicing was found
int IsSpliceSite(char *genome, unsigned int gapStartInGenome, int gapSizeInGenome) {
	int seedsShift;
	unsigned int leftPos, rightPos;
	
	// first, try to shift the gap to the left, if possible
	// both variables bellow have (+1) because they will be decremented next at the beginning of the loop
	leftPos = (gapStartInGenome - 1) + 1; // rightmost position of the left seed
	rightPos = (gapStartInGenome + gapSizeInGenome - 1) + 1; // rightmost position of the gap (one position before the right seed)
	seedsShift = +1; // count how many positions to the left the ends of the seeds moved (+1 so it starts at 0 next on the loop)
	do { // look for splice signals on both ends of the gap, right next to the seeds
		leftPos--; // move both pointers one position to the left
		rightPos--;
		seedsShift--;
		// both pointers next have (+1) and (-1) so they are placed at the beginning of the dinucleotide donor/acceptor splice signals
		if ( CheckSite_3CanonicalAndRC(genome, (leftPos + 1), (rightPos - 1)) ) { // check splice site
			return seedsShift;
		}
		// keep looking for signals and shift the gap to the left,
		// while the rightmost chars of the left seed match the rightmost chars of the gap (immediately at the left of the right seed)
	} while (genome[leftPos] == genome[rightPos]);

	// second, try to shift the gap to the right, if possible
	leftPos = (gapStartInGenome); // leftmost position of the gap (one position after the left seed)
	rightPos = (gapStartInGenome + gapSizeInGenome);  // leftmost position of the right seed
	seedsShift = 0; // count how many positions to the right the ends of the seeds moved (it will start at +1 next)
	// keep looking for signals and shift the gap to the right,
	// while the leftmost chars of the gap (immediately at the right of the left seed) match the leftmost chars of the right seed
	while (genome[leftPos] == genome[rightPos]) { // look for splice signals on both ends of the gap, right next to the seeds
		leftPos++; // move both pointers one position to the right
		rightPos++;
		seedsShift++;
		// both pointers next have (+0) and (-2) so they are placed at the beginning of the dinucleotide donor/acceptor splice signals
		if ( CheckSite_3CanonicalAndRC(genome, (leftPos), (rightPos - 2)) ) { // check splice site
			return seedsShift;
		}
	}
	
	// if shifting to either direction did not work, report failure
	return INT_MIN;
}

/*
int IsSpliceSite(char *genome, unsigned int gapStartInGenome, int gapSizeInGenome) {
	int seedsShift;
	unsigned int leftPos, rightPos;

	leftPos = (gapStartInGenome - 1); // rightmost position of the left seed
	rightPos = (gapStartInGenome + gapSizeInGenome - 1); // rightmost position of the gap (one position before the right seed)
	
	// check for splice sites in the original position
	// both pointers next have (+1) and (-1) respectively so they are placed at the beginning of the dinucleotide donor/acceptor splice signals
	if (CheckSite_3CanonicalAndRC(genome, (leftPos + 1) + 0, (rightPos - 1) + 0)) return (+0); // check for splice site on current position
	if (CheckSite_3CanonicalAndRC(genome, (leftPos + 1) - 1, (rightPos - 1) - 1)) return (-1); // check for splice site one position to the left
	if (CheckSite_3CanonicalAndRC(genome, (leftPos + 1) + 1, (rightPos - 1) + 1)) return (+1); // check for splice site one position to the right
	
	// first, shift the gap all the way to the left, if possible
	seedsShift = 0; // count how many positions to the left the ends of the seeds moved
	while (genome[leftPos] == genome[rightPos]) {
		leftPos--; // move both pointers to the left
		rightPos--;
		seedsShift--;
	}
	if (seedsShift != 0) { // if the seeds shifted, check for splice sites in the new position
		if (CheckSite_3CanonicalAndRC(genome, (leftPos + 1) + 0, (rightPos - 1) + 0)) return (seedsShift + 0); // check for splice site on current position
		if (CheckSite_3CanonicalAndRC(genome, (leftPos + 1) - 1, (rightPos - 1) - 1)) return (seedsShift - 1); // check for splice site one position to the left
		if (CheckSite_3CanonicalAndRC(genome, (leftPos + 1) + 1, (rightPos - 1) + 1)) return (seedsShift + 1); // check for splice site one position to the right
	}

	leftPos = (gapStartInGenome); // leftmost position of the gap (one position after the left seed)
	rightPos = (gapStartInGenome + gapSizeInGenome); // leftmost position of the right seed
	
	// second, shift the gap all the way to the right, if possible
	seedsShift = 0; // count how many positions to the right the ends of the seeds moved
	while (genome[leftPos] == genome[rightPos]) {
		leftPos++; // move both pointers to the right
		rightPos++;
		seedsShift++;
	}
	// both pointers next have (+0) and (-2) so they are placed at the beginning of the dinucleotide donor/acceptor splice signals
	if (seedsShift != 0) { // if the seeds shifted, check for splice sites in the new position
		if (CheckSite_3CanonicalAndRC(genome, (leftPos) + 0, (rightPos - 2) + 0)) return (seedsShift + 0); // check for splice site on current position
		if (CheckSite_3CanonicalAndRC(genome, (leftPos) - 1, (rightPos - 2) - 1)) return (seedsShift - 1); // check for splice site one position to the left
		if (CheckSite_3CanonicalAndRC(genome, (leftPos) + 1, (rightPos - 2) + 1)) return (seedsShift + 1); // check for splice site one position to the right
	}

	// if no splice sites were found, report failure
	return INT_MIN;
}
*/
