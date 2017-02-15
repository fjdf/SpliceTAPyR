#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "bwtindex.h"
#include "tools.h"

#ifdef _MSC_VER
	#pragma warning(disable:4996)
#endif

#define DEBUG_BWTIS 1

#define FILEHEADER "IDX0"

typedef struct _IndexBlock {
	unsigned int bwtLowBits;
	unsigned int bwtHighBits;
	unsigned int specialLettersMask;
	unsigned int letterJumpsSample[5];
	unsigned int textPositionSample;
} IndexBlock;

unsigned int offsetMasks[32] = { // = (1UL<<offset)
	0x00000001, // 1st bit
	0x00000002, // 2nd bit
	0x00000004, // 3rd bit
	0x00000008, // 4th bit
	0x00000010, // 5th bit
	0x00000020, // 6th bit
	0x00000040, // 7th bit
	0x00000080, // 8th bit
	0x00000100, // 9th bit
	0x00000200, // 10th bit
	0x00000400, // 11th bit
	0x00000800, // 12th bit
	0x00001000, // 13th bit
	0x00002000, // 14th bit
	0x00004000, // 15th bit
	0x00008000, // 16th bit
	0x00010000, // 17th bit
	0x00020000, // 18th bit
	0x00040000, // 19th bit
	0x00080000, // 20th bit
	0x00100000, // 21st bit
	0x00200000, // 22nd bit
	0x00400000, // 23rd bit
	0x00800000, // 24th bit
	0x01000000, // 25th bit
	0x02000000, // 26th bit
	0x04000000, // 27th bit
	0x08000000, // 28th bit
	0x10000000, // 29th bit
	0x20000000, // 30th bit
	0x40000000, // 31st bit
	0x80000000  // 32nd bit
};
unsigned int inverseOffsetMasks[32] = { // = ~(1UL<<offset)
	0xFFFFFFFE, // 1st bit
	0xFFFFFFFD, // 2nd bit
	0xFFFFFFFB, // 3rd bit
	0xFFFFFFF7, // 4th bit
	0xFFFFFFEF, // 5th bit
	0xFFFFFFDF, // 6th bit
	0xFFFFFFBF, // 7th bit
	0xFFFFFF7F, // 8th bit
	0xFFFFFEFF, // 9th bit
	0xFFFFFDFF, // 10th bit
	0xFFFFFBFF, // 11th bit
	0xFFFFF7FF, // 12th bit
	0xFFFFEFFF, // 13th bit
	0xFFFFDFFF, // 14th bit
	0xFFFFBFFF, // 15th bit
	0xFFFF7FFF, // 16th bit
	0xFFFEFFFF, // 17th bit
	0xFFFDFFFF, // 18th bit
	0xFFFBFFFF, // 19th bit
	0xFFF7FFFF, // 20th bit
	0xFFEFFFFF, // 21st bit
	0xFFDFFFFF, // 22nd bit
	0xFFBFFFFF, // 23rd bit
	0xFF7FFFFF, // 24th bit
	0xFEFFFFFF, // 25th bit
	0xFDFFFFFF, // 26th bit
	0xFBFFFFFF, // 27th bit
	0xF7FFFFFF, // 28th bit
	0xEFFFFFFF, // 29th bit
	0xDFFFFFFF, // 30th bit
	0xBFFFFFFF, // 31st bit
	0x7FFFFFFF  // 32nd bit
};
unsigned int firstLettersMasks[33] = { // all bits before the offset = ((1UL<<offset)-1UL)
	0x00000000, // lower 0 letters
	0x00000001, // lower 1 letters
	0x00000003, // lower 2 letters
	0x00000007, // lower 3 letters
	0x0000000F, // lower 4 letters
	0x0000001F, // lower 5 letters
	0x0000003F, // lower 6 letters
	0x0000007F, // lower 7 letters
	0x000000FF, // lower 8 letters
	0x000001FF, // lower 9 letters
	0x000003FF, // lower 10 letters
	0x000007FF, // lower 11 letters
	0x00000FFF, // lower 12 letters
	0x00001FFF, // lower 13 letters
	0x00003FFF, // lower 14 letters
	0x00007FFF, // lower 15 letters
	0x0000FFFF, // lower 16 letters
	0x0001FFFF, // lower 17 letters
	0x0003FFFF, // lower 18 letters
	0x0007FFFF, // lower 19 letters
	0x000FFFFF, // lower 20 letters
	0x001FFFFF, // lower 21 letters
	0x003FFFFF, // lower 22 letters
	0x007FFFFF, // lower 23 letters
	0x00FFFFFF, // lower 24 letters
	0x01FFFFFF, // lower 25 letters
	0x03FFFFFF, // lower 26 letters
	0x07FFFFFF, // lower 27 letters
	0x0FFFFFFF, // lower 28 letters
	0x1FFFFFFF, // lower 29 letters
	0x3FFFFFFF, // lower 30 letters
	0x7FFFFFFF, // lower 31 letters
	0xFFFFFFFF  // lower 32 letters
};
unsigned int lastLettersMasks[33] = { // all bits at or after offset = ~((1UL<<offset)-1UL)
	0xFFFFFFFF, // higher 32 letters
	0xFFFFFFFE, // higher 31 letters
	0xFFFFFFFC, // higher 30 letters
	0xFFFFFFF8, // higher 29 letters
	0xFFFFFFF0, // higher 28 letters
	0xFFFFFFE0, // higher 27 letters
	0xFFFFFFC0, // higher 26 letters
	0xFFFFFF80, // higher 25 letters
	0xFFFFFF00, // higher 24 letters
	0xFFFFFE00, // higher 23 letters
	0xFFFFFC00, // higher 22 letters
	0xFFFFF800, // higher 21 letters
	0xFFFFF000, // higher 20 letters
	0xFFFFE000, // higher 19 letters
	0xFFFFC000, // higher 18 letters
	0xFFFF8000, // higher 17 letters
	0xFFFF0000, // higher 16 letters
	0xFFFE0000, // higher 15 letters
	0xFFFC0000, // higher 14 letters
	0xFFF80000, // higher 13 letters
	0xFFF00000, // higher 12 letters
	0xFFE00000, // higher 11 letters
	0xFFC00000, // higher 10 letters
	0xFF800000, // higher 9 letters
	0xFF000000, // higher 8 letters
	0xFE000000, // higher 7 letters
	0xFC000000, // higher 6 letters
	0xF8000000, // higher 5 letters
	0xF0000000, // higher 4 letters
	0xE0000000, // higher 3 letters
	0xC0000000, // higher 2 letters
	0x80000000, // higher 1 letters
	0x00000000  // higher 0 letters
};
unsigned int searchOffsetMasks[32] = { // all bits before or at the offset, except 1st one = (((1UL<<(offset+1))-1UL)&(~1UL))
	0x00000000, // lower 1 letters
	0x00000002, // lower 2 letters
	0x00000006, // lower 3 letters
	0x0000000E, // lower 4 letters
	0x0000001E, // lower 5 letters
	0x0000003E, // lower 6 letters
	0x0000007E, // lower 7 letters
	0x000000FE, // lower 8 letters
	0x000001FE, // lower 9 letters
	0x000003FE, // lower 10 letters
	0x000007FE, // lower 11 letters
	0x00000FFE, // lower 12 letters
	0x00001FFE, // lower 13 letters
	0x00003FFE, // lower 14 letters
	0x00007FFE, // lower 15 letters
	0x0000FFFE, // lower 16 letters
	0x0001FFFE, // lower 17 letters
	0x0003FFFE, // lower 18 letters
	0x0007FFFE, // lower 19 letters
	0x000FFFFE, // lower 20 letters
	0x001FFFFE, // lower 21 letters
	0x003FFFFE, // lower 22 letters
	0x007FFFFE, // lower 23 letters
	0x00FFFFFE, // lower 24 letters
	0x01FFFFFE, // lower 25 letters
	0x03FFFFFE, // lower 26 letters
	0x07FFFFFE, // lower 27 letters
	0x0FFFFFFE, // lower 28 letters
	0x1FFFFFFE, // lower 29 letters
	0x3FFFFFFE, // lower 30 letters
	0x7FFFFFFE, // lower 31 letters
	0xFFFFFFFE  // lower 32 letters
};
unsigned int firstLetterMask = 0x00000001; // lowest bit
unsigned int lastLetterMask  = 0x80000000; // highest bit
unsigned int secondLetterMask = 0x00000002; // 2nd lowest bit
unsigned int firstHalfLettersMask = 0x0000FFFF; // lowest 16 bits
unsigned int lastHalfLettersMask  = 0xFFFF0000; // highest 16 bits
unsigned int letterLowBitMasks[6] = {
	0x00000000, // unused (for '$' (00))
	0xFFFFFFFF, // unused (for 'N' (01))
	0x00000000, // 1st bit mask for 'A' (00): 0...0
	0xFFFFFFFF, // 1st bit mask for 'C' (01): 1...1
	0x00000000, // 1st bit mask for 'G' (10): 0...0
	0xFFFFFFFF  // 1st bit mask for 'T' (11): 1...1
};
unsigned int letterHighBitMasks[6] = {
	0x00000000, // unused (for '$' (00))
	0x00000000, // unused (for 'N' (01))
	0x00000000, // 2nd bit mask for 'A' (00): 0...0
	0x00000000, // 2nd bit mask for 'C' (01): 0...0
	0xFFFFFFFF, // 2nd bit mask for 'G' (10): 1...1
	0xFFFFFFFF  // 2nd bit mask for 'T' (11): 1...1
};
unsigned int inverseLetterLowBitMasks[6] = {
	0xFFFFFFFF, // unused (for '$' (00))
	0x00000000, // unused (for 'N' (01))
	0xFFFFFFFF, // 1st bit mask for 'A' (00): ~0...0 = 1...1
	0x00000000, // 1st bit mask for 'C' (01): ~1...1 = 0...0
	0xFFFFFFFF, // 1st bit mask for 'G' (10): ~0...0 = 1...1
	0x00000000  // 1st bit mask for 'T' (11): ~1...1 = 0...0
};
unsigned int inverseLetterHighBitMasks[6] = {
	0xFFFFFFFF, // unused (for '$' (00))
	0xFFFFFFFF, // unused (for 'N' (01))
	0xFFFFFFFF, // 2nd bit mask for 'A' (00): ~0...0 = 1...1
	0xFFFFFFFF, // 2nd bit mask for 'C' (01): ~0...0 = 1...1
	0x00000000, // 2nd bit mask for 'G' (10): ~1...1 = 0...0
	0x00000000  // 2nd bit mask for 'T' (11): ~1...1 = 0...0
};
char letterChars[6] = { '$' , 'N' , 'A' , 'C' , 'G' , 'T' }; // get letter char from letter id
unsigned int sampleIntervalShift = 5; // sample interval of 32 positions (2^5=32)
unsigned int sampleIntervalMask = 0x0000001F; // = ((1<<sampleIntervalShift)-1) = (32-1)
unsigned int sampleIntervalSize = 32; // = (1<<sampleIntervalShift) = 32
unsigned int sampleIntervalHalfSize = 16; // = (1<<(sampleIntervalShift-1)) = 16

// varibles needed to load/store index
struct _IndexBlock *Index = NULL;
unsigned int textSize = 0; // inside the index functions, textSize always counts with the terminator char
unsigned int numSamples = 0;
unsigned int lastBwtPos = 0; // position in the BWT of the 0-th entry of the last sample (multiple of 32 to speed up search of 1st char of pattern)
char *textFilename = NULL;
unsigned char *letterIds = NULL;

void SetCharAtBWTPos(unsigned int charid, unsigned int bwtpos){
	unsigned int sample = ( bwtpos >> sampleIntervalShift );
	IndexBlock *block = &(Index[sample]); // get sample block
	unsigned int offset = ( bwtpos & sampleIntervalMask );
	unsigned int mask = inverseOffsetMasks[offset]; // get all bits except the one at the offset
	(block->specialLettersMask) &= mask; // reset special letter mark
	(block->bwtLowBits) &= mask; // reset low bit
	(block->bwtHighBits) &= mask; // reset high bit
	mask = offsetMasks[offset]; // get only the bit at the offset
	if( charid < 2 ) (block->specialLettersMask) |= mask; // set special letter mark if needed (for ids 0 and 1)
	(block->bwtLowBits) |= ( letterLowBitMasks[charid] & mask ); // set low bit
	(block->bwtHighBits) |= ( letterHighBitMasks[charid] & mask ); // set high bit
}

unsigned int GetCharIdAtBWTPos(unsigned int bwtpos){
	unsigned int sample, offset, mask, charid;
	IndexBlock *block;
	sample = ( bwtpos >> sampleIntervalShift );
	offset = ( bwtpos & sampleIntervalMask );
	block = &(Index[sample]); // get sample block
	mask = offsetMasks[offset]; // get only the bit at the offset
	charid = ( ( (block->bwtLowBits) >> offset ) & firstLetterMask ); // get low bit
	charid |= ( ( ( (block->bwtHighBits) >> offset ) & firstLetterMask ) << 1 ); // get high bit
	if( !( (block->specialLettersMask) & mask ) ) charid += 2; // check if special letter bit is set, because regular ids start at id=2
	return charid;
}

char GetCharAtBWTPos(unsigned int bwtpos){
	return letterChars[ GetCharIdAtBWTPos(bwtpos) ]; // get letter char
}

/*
unsigned int GetRightCharIdAtBWTPos(unsigned int bwtpos){
	unsigned int charid;
	charid = 5;
	while( (charid != 0) && (bwtpos < letterStartPos[charid]) ) charid--;
	return charid;
}
*/

void InitializeLetterIdsArray(){
	int i;
	letterIds=(unsigned char *)malloc(256*sizeof(unsigned char));
	for(i=0;i<256;i++) letterIds[i]=(unsigned char)1; // ACGT, N and $
	letterIds[(int)'$']=(unsigned char)0;
	letterIds[(int)'N']=(unsigned char)1;
	letterIds[(int)'A']=(unsigned char)2;
	letterIds[(int)'C']=(unsigned char)3;
	letterIds[(int)'G']=(unsigned char)4;
	letterIds[(int)'T']=(unsigned char)5;
	letterIds[(int)'n']=(unsigned char)1;
	letterIds[(int)'a']=(unsigned char)2;
	letterIds[(int)'c']=(unsigned char)3;
	letterIds[(int)'g']=(unsigned char)4;
	letterIds[(int)'t']=(unsigned char)5;
}

unsigned int FMI_GetTextSize(){
	return (textSize-1); // inside the index functions, the textSize variable counts the terminator char too
}

unsigned int FMI_GetBWTSize(){
	return lastBwtPos; // last valid filled position of the BWT (aligned to a multiple of the sample interval)
}

char *FMI_GetTextFilename(){
	return textFilename;
}

// TODO: find way to remove the decrement in letterJumpsSample[(letterId-1)]
// TODO: move code with the check for special letters to another function, because on normal search we never use jumps by $'s or N's (but on position serch, we do by N's)
unsigned int FMI_LetterJump( unsigned int letterId , unsigned int bwtPos ){
	unsigned int offset, bitArray, letterJump;
	IndexBlock *block;
	block = &(Index[( bwtPos >> sampleIntervalShift )]);
	offset = ( bwtPos & sampleIntervalMask );
	bitArray = searchOffsetMasks[offset]; // all bits bellow offset, except 1st one
	if( letterId < 2 ) bitArray &= (block->specialLettersMask); // if we are looking for '$' or 'N', keep only special chars
	else bitArray &= ( ~ (block->specialLettersMask) ); // remove special chars
	bitArray &= ( (block->bwtLowBits) ^ inverseLetterLowBitMasks[letterId] ); // keep only ones with the same low bit
	bitArray &= ( (block->bwtHighBits) ^ inverseLetterHighBitMasks[letterId] ); // keep only ones with the same high bit
	letterJump = (block->letterJumpsSample[(letterId-1)]); // get top letter jump
	while( bitArray ){
		letterJump++; // of other equal letter exist, increase letter jump
		bitArray &= ( bitArray - 1U ); // clear lower occurrence bit
	}
	return letterJump;
}

// TODO: find better way to deal with original toppointer value
// TODO: move code of called functions inside and reuse variables for speed up
// NOTE: returns the size of the BWT interval if a match exists, and 0 otherwise
unsigned int FMI_FollowLetter( char c , unsigned int *topPointer , unsigned int *bottomPointer ){
	unsigned int charid, originaltoppointer;
	originaltoppointer = (*topPointer);
	charid = letterIds[(int)c];
	(*topPointer) = FMI_LetterJump( charid , originaltoppointer );
	if( GetCharIdAtBWTPos(originaltoppointer) != charid ) (*topPointer)++; // if the letter is not in the topPointer position (the initial, not the updated one), its next occurrence is after that
	(*bottomPointer) = FMI_LetterJump( charid , (*bottomPointer) );
	if( (*topPointer) > (*bottomPointer) ) return 0;
	return ( (*bottomPointer) - (*topPointer) + 1 );
	/*
	unsigned int bwtPos, letterId, offset, bitArray, letterJump;
	IndexBlock *block;
	letterId = letterIds[c];
	bwtPos = (*topPointer); // process top pointer
	block = &(Index[( bwtPos >> sampleIntervalShift )]); // get sample block
	bitArray = ( ~ (block->specialLettersMask) ); // remove special chars
	bitArray &= ( (block->bwtLowBits) ^ inverseLetterLowBitMasks[letterId] ); // keep only ones with the same low bit
	bitArray &= ( (block->bwtHighBits) ^ inverseLetterHighBitMasks[letterId] ); // keep only ones with the same high bit
	letterJump = (block->letterJumpsSample[letterId]); // get top letter jump
	offset = ( bwtPos & sampleIntervalMask ); // get offset inside sample
	if( !( bitArray & offsetMasks[offset] ) ) letterJump++; // if the letter is not in the topPointer position, its next occurrence is bellow that
	bitArray &= searchOffsetMasks[offset]; // keep all bits bellow offset, except 1st one
	while( bitArray ){ // while there are occurrences of this letter in the interval
		letterJump++; // if other equal letter exist, increase letter jump
		bitArray &= ( bitArray - 1U ); // clear lower occurrence bit
	}
	(*topPointer) = letterJump;
	bwtPos = (*bottomPointer); // process bottom pointer
	block = &(Index[( bwtPos >> sampleIntervalShift )]);
	offset = ( bwtPos & sampleIntervalMask );
	bitArray = searchOffsetMasks[offset]; // all bits bellow offset, except 1st one
	bitArray &= ( ~ (block->specialLettersMask) ); // remove special chars
	bitArray &= ( (block->bwtLowBits) ^ inverseLetterLowBitMasks[letterId] ); // keep only ones with the same low bit
	bitArray &= ( (block->bwtHighBits) ^ inverseLetterHighBitMasks[letterId] ); // keep only ones with the same high bit
	letterJump = (block->letterJumpsSample[letterId]); // get top letter jump
	while( bitArray ){
		letterJump++; // if other equal letter exist, increase letter jump
		bitArray &= ( bitArray - 1U ); // clear lower occurrence bit
	}
	(*bottomPointer) = letterJump;
	if( (*topPointer) > (*bottomPointer) ) return 0;
	return ( (*bottomPointer) - (*topPointer) + 1 );
	*/
}

// TODO: move code of called functions inside and reuse variables for speed up
unsigned int FMI_PositionInText( unsigned int bwtpos ){
	unsigned int charid, addpos;
	addpos = 0;
	while( bwtpos & sampleIntervalMask ){ // move backwards until we land on a position with a sample
		charid = GetCharIdAtBWTPos(bwtpos);
		if( charid == 0 ) return addpos; // check if this is the terminator char
		bwtpos = FMI_LetterJump( charid , bwtpos ); // follow the left letter backwards
		addpos++; // one more position away from our original position
	}
	return ( (Index[( bwtpos >> sampleIntervalShift )].textPositionSample) + addpos );
}

void FMI_LoadIndex(char *indexfilename){
	FILE *indexfile;
	size_t readcount;
	fpos_t filepos;
	unsigned int i;
	long long int seqstart, seqend, seqsize;
	char c;
	char fileHeader[5] = FILEHEADER;
	printf("> Loading index from file <%s> ... ",indexfilename);
	fflush(stdout);
	indexfile = fopen(indexfilename,"rb");
	if( indexfile == NULL ){
		printf("\n> ERROR: Cannot open file\n");
		exit(0);
	}
	seqstart=(long long int)ftell(indexfile);
	fseek(indexfile,0L,SEEK_END);
	seqend=(long long int)ftell(indexfile);
	seqsize=(seqend-seqstart);
	rewind(indexfile);
	printf("(");PrintNumber(seqsize);printf(" bytes)");
	fflush(stdout);
	for(i=0;i<4;i++){ // check if header is "IDX0"
		fread( &c , sizeof(char) , (size_t)1 , indexfile );
		if( c != fileHeader[i] ) break;
	}
	if( i != 4 ){
		printf("\n> ERROR: Invalid index file\n");
		exit(0);
	}
	fgetpos(indexfile,&filepos);
	i=0;
	while( c!='\0' && c!=EOF ){ // get text filename size
		fread( &c , sizeof(char) , (size_t)1 , indexfile );
		i++;
	}
	textFilename=(char *)calloc((i),sizeof(char));
	fsetpos(indexfile,&filepos);
	fread( textFilename , sizeof(char) , (size_t)i , indexfile ); // get text filename chars
	fread( &textSize , sizeof(unsigned int) , (size_t)1 , indexfile ); // counts the terminator char too, so it is actually the BWT size
	fread( &numSamples , sizeof(unsigned int) , (size_t)1 , indexfile );
	#ifdef DEBUG
	printf("\n  [textFilename=\"%s\";textSize=%u;numSamples=%u]",textFilename,textSize,numSamples);
	fflush(stdout);
	#endif
	i = ( ( ( textSize-1 ) >> sampleIntervalShift ) + 1 );
	if( ((textSize-1) & sampleIntervalMask) != 0 ) i++; // extra sample
	if( numSamples != i ){ // check if number of samples is correct based on text size
		printf("\n> ERROR: Invalid index data\n");
		exit(0);
	}
	Index = (IndexBlock *)malloc( numSamples * sizeof(IndexBlock) ); // allocate memory for all index blocks
	if( Index == NULL ){
		printf("\n> ERROR: Not enough memory\n");
		exit(0);
	}
	readcount = fread( Index , sizeof(IndexBlock) , (size_t)(numSamples) , indexfile ); // read all index blocks
	if( readcount != (size_t)numSamples ){
		printf("\n> ERROR: Incomplete index data\n");
		exit(0);
	}
	fclose(indexfile);
	InitializeLetterIdsArray(); // initialize arrays used by index functions
	lastBwtPos = ((numSamples-1) << sampleIntervalShift); // set last BWT pos for pattern matching initialization
	printf(" OK\n");
	fflush(stdout);
}

// TODO: use indexBlocksStartInFile, distanceToNextBlockBits
void FMI_SaveIndex(char *indexfilename){
	FILE *indexfile;
	size_t writecount;
	long long int filesize;
	int i;
	printf("> Saving index to file <%s> ... ",indexfilename);
	fflush(stdout);
	indexfile = fopen(indexfilename,"wb");
	if( indexfile == NULL ){
		printf("\n> ERROR: Cannot create file\n");
		exit(0);
	}
	i=0;
	while(textFilename[i]!='\0') i++; // get text filename size
	i++;
	fwrite( (FILEHEADER) , sizeof(char) , (size_t)4 , indexfile );
	fwrite( textFilename , sizeof(char) , (size_t)i , indexfile );
	fwrite( &textSize , sizeof(unsigned int) , (size_t)1 , indexfile );
	fwrite( &numSamples , sizeof(unsigned int) , (size_t)1 , indexfile );
	writecount = fwrite( Index , sizeof(IndexBlock) , (size_t)(numSamples) , indexfile );
	if( writecount != (size_t)(numSamples) ){
		printf("\n> ERROR: Cannot write file\n");
		exit(0);
	}
	filesize=(long long int)ftell(indexfile);
	fclose(indexfile);
	printf("(");PrintNumber(filesize);printf(" bytes) OK\n");
	fflush(stdout);
}

void FMI_FreeIndex(){
	if(textFilename!=NULL){
		free(textFilename);
		textFilename=NULL;
	}
	if(Index!=NULL){
		free(Index);
		Index=NULL;
	}
	if(letterIds!=NULL){
		free(letterIds);
		letterIds=NULL;
	}
}

void PrintBWT(char *text, unsigned int *letterStartPos){
	unsigned int i, n, p;
	printf("%u {", textSize );
	for( n = 1 ; n < 6 ; n++ ){
		printf(" %c [%02u-%02u] %c", letterChars[n] , letterStartPos[n] , (n==5)?(textSize-1):(letterStartPos[(n+1)]-1) , (n==5)?'}':',' );
	}
	printf(" 2^%u=%u %u %#.8X\n", sampleIntervalShift , sampleIntervalSize , numSamples , sampleIntervalMask );
	printf("[  ] (  ) {");
	for( n = 1 ; n < 6 ; n++ ){
		printf(" %c%c", letterChars[n] , (n==5)?'}':',' );
	}
	printf("\n");
	for( i = 0 ; i < textSize ; i++ ){ // position in BWT
		p = FMI_PositionInText( i );
		printf("[%02u]%c(%2u) {", i , (i & sampleIntervalMask)?' ':'*' , p );
		for( n = 1 ; n < 6 ; n++ ){
			printf("%02u%c", FMI_LetterJump( n , i ) , (n==5)?'}':',' );
		}
		printf(" %c ", GetCharAtBWTPos(i) );
		n = 0;
		while( ( n < 5 ) && ( i >= letterStartPos[(n+1)] ) ) n++;
		printf(" %c ", letterChars[n] );
		if( text != NULL ) printf("%s", (char *)(text+p+1) );
		printf("\n");
	}
}


typedef struct _PackedNumberArray {
	unsigned long long *bitsArray;
	unsigned char bitsPerInt;
	unsigned int numWords;
	//unsigned char bitsPerWord; // = 64
} PackedNumberArray;

PackedNumberArray *NewPackedNumberArray(unsigned int numInts, unsigned int maxInt){
	PackedNumberArray *intArray;
	unsigned long long numBits;
	unsigned int n;
	intArray = (PackedNumberArray *)malloc(sizeof(PackedNumberArray));
	//(intArray->bitsPerWord) = (unsigned char)(sizeof(unsigned long long)*8); // use 64 bit words (8 bytes * 8 bits/byte)
	n = 1; // number of bits needed to store one number
	while( ( (1UL << n) - 1 ) < maxInt ) n++; // n bits per number (stores 2^n numbers, but the last one is (2^n-1))
	(intArray->bitsPerInt) = (unsigned char)n;
	numBits = ( ((unsigned long long)numInts) * ((unsigned long long)n) ); // total number of bits occupied by all the numbers
	if( numBits != 0) numBits--; // if it was a multiple of 64 , it would create an extra unused word
	n = (unsigned int)( (numBits/64ULL) + 1ULL); // number of 64 bit words required to store (numInts) numbers of (bitsPerInt) bits each
	(intArray->numWords) = n;
	(intArray->bitsArray) = (unsigned long long *)calloc(n,sizeof(unsigned long long)); // bit array that will store the numbers
	return intArray;
}

void FreePackedNumberArray(PackedNumberArray *intArray){
	free(intArray->bitsArray);
	free(intArray);
}

void ResetPackedNumberArray(PackedNumberArray *intArray){
	unsigned int n;
	n = (intArray->numWords);
	while( n != 0 ) (intArray->bitsArray)[(--n)] = 0ULL;
}

unsigned int GetPackedNumber(PackedNumberArray *intArray, unsigned int pos){
	unsigned char offset, bitsPerInt;
	unsigned long long temp;
	unsigned int number;
	bitsPerInt = (intArray->bitsPerInt);
	temp = ((unsigned long long)(pos)) * ((unsigned long long)(bitsPerInt)); // number of bits
	offset = (unsigned char)( temp & 63ULL ); // mask by (64-1) to get position of 1st bit inside the word
	pos = (unsigned int)( temp >> 6 ); // divide by 64 to get word position inside bit array
	temp = ( ( 1ULL << bitsPerInt ) - 1ULL ); // mask for the n bits of each number
	number = (unsigned int)( ( (intArray->bitsArray)[pos] >> offset ) & temp );
	offset = ( 64 - offset ); // check if the bits of the number extend to the next word (get number of bits left until the end of this word)
	if( offset < bitsPerInt ) number |= (unsigned int)( ( (intArray->bitsArray)[(++pos)] << offset ) & temp );
	return number;
}

void SetPackedNumber(PackedNumberArray *intArray, unsigned int pos, unsigned int num){
	unsigned char offset, bitsPerInt;
	unsigned long long temp;
	bitsPerInt = (intArray->bitsPerInt);
	temp = ((unsigned long long)(pos)) * ((unsigned long long)(bitsPerInt)); // number of bits
	offset = (unsigned char)( temp & 63ULL ); // mask by (64-1) to get position of 1st bit inside the word
	pos = (unsigned int)( temp >> 6 ); // divide by 64 to get word position inside bit array
	(intArray->bitsArray)[pos] |= ( ((unsigned long long)num) << offset );
	offset = ( 64 - offset ); // check if the bits of the number extend to the next word (get number of bits left until the end of this word)
	if( offset < bitsPerInt ) (intArray->bitsArray)[(++pos)] |= ( ((unsigned long long)num) >> offset );
}

void ReplacePackedNumber(PackedNumberArray *intArray, unsigned int pos, unsigned int num){
	unsigned char offset, bitsPerInt;
	unsigned long long temp;
	bitsPerInt = (intArray->bitsPerInt);
	temp = ((unsigned long long)(pos)) * ((unsigned long long)(bitsPerInt)); // number of bits
	offset = (unsigned char)( temp & 63ULL ); // mask by (64-1) to get position of 1st bit inside the word
	pos = (unsigned int)( temp >> 6 ); // divide by 64 to get word position inside bit array
	temp = ( ((unsigned long long)num) << offset ); // rightmost bits to set
	(intArray->bitsArray)[pos] &= temp; // only keep common 1 bits
	(intArray->bitsArray)[pos] |= temp; // set the missing 1 bits (same as reset and then set)
	offset = ( 64 - offset ); // check if the bits of the number extend to the next word (get number of bits left until the end of this word)
	if( offset < bitsPerInt ){
		pos++; // next word
		temp = ( ((unsigned long long)num) >> offset ); // leftmost bits to set
		(intArray->bitsArray)[pos] &= temp; // only keep common 1 bits
		(intArray->bitsArray)[pos] |= temp; // set the missing 1 bits (same as reset and then set)
	}
}

typedef struct _PackedIncreasingNumberArray {
	unsigned char numTotalBits;
	unsigned char numHighBits;
	unsigned char numLowBits;
	unsigned int highBitsMask;
	unsigned int lowBitsMask;
	unsigned int *highLimits;
	PackedNumberArray *lowBits;
} PackedIncreasingNumberArray;

PackedIncreasingNumberArray *NewPackedIncreasingNumberArray(unsigned int numInts, unsigned int maxInt){
	PackedIncreasingNumberArray *incIntArray;
	unsigned long long numBits;
	unsigned int k;
	unsigned char n;
	incIntArray = (PackedIncreasingNumberArray *)malloc(sizeof(PackedIncreasingNumberArray));
	//k = (unsigned int)(sizeof(unsigned int)*8); // use 32 bits words
	numBits = (unsigned long long)(numInts*32); // total number of bits to store all the words regularly
	n = 0; // get best number of high bits to compact
	while( (((1ULL << n)-1ULL)*32ULL) < numBits ) n++; // for compacting n high bits we need (2^n-1) extra numbers
	if( n != 0 ) n--; // the current value already exceeded the regular size
	(incIntArray->numHighBits) = n;
	n = 0; // number of bits per number
	k = 1; // power of two (2^n)
	while( k && ( (k-1) < maxInt ) ){ // (2^n-1) is the largest number that can be represented by n bits
		n++;
		k <<= 1;
	}
	if( k != 0 ) k--; // mask for all bits
	else k = (~k); // n = 32, whole word (~0UL)
	(incIntArray->numTotalBits) = n; // total bits of each number
	n -= (incIntArray->numHighBits); // number of lower bits
	(incIntArray->numLowBits) = n;
	(incIntArray->lowBitsMask) = ( (1UL << n) - 1UL ); // mask for lower bits
	(incIntArray->highBitsMask) = ( k ^ (incIntArray->lowBitsMask) ); // all bits of number except the low ones
	k = ( ( 1UL << (incIntArray->numHighBits) ) - 1UL ); // number of limits needed for compacting the high bits
	(incIntArray->highLimits) = (unsigned int *)calloc(k,sizeof(unsigned int));
	k = ( ( 1UL << (incIntArray->numLowBits) ) - 1UL ); // highest number representable by the low bits
	(incIntArray->lowBits) = NewPackedNumberArray(numInts,k);
	return incIntArray;
}

void FreePackedIncreasingNumberArray(PackedIncreasingNumberArray *incIntArray){
	FreePackedNumberArray(incIntArray->lowBits);
	free(incIntArray->highLimits);
	free(incIntArray);
}

void ResetPackedIncreasingNumberArray(PackedIncreasingNumberArray *incIntArray){
	unsigned int n;
	n = ( ( 1UL << (incIntArray->numHighBits) ) - 1UL ); // number of limits used
	while( n != 0 ) (incIntArray->highLimits)[(--n)] = 0UL;
	ResetPackedNumberArray(incIntArray->lowBits);
}

unsigned int GetPackedIncreasingNumber(PackedIncreasingNumberArray *incIntArray, unsigned int pos){
	unsigned char n;
	unsigned int number, limit, i, k;
	number = 0UL;
	i = pos; // position in the increasing numbers array
	k = 0; // current position in the limits array
	n = (incIntArray->numHighBits); // current high bit
	while( n && i ){
		limit = (incIntArray->highLimits)[k];
		if( i < limit ){ // bit 0, and advance to limit (2*k+1)
			k <<= 1;
			k++;
		} else { // bit 1, and advance to limit (2*(k+1))
			k++;
			k <<= 1;
			i -= limit; // fix array position relative to the current limit
			number++; // the number has a bit 1 at this position
		}
		n--;
		number <<= 1; // go to next bit of number
	}
	number <<= n; // if we didn't check all high bits, shift the number correctly
	return ( ( number << (incIntArray->numLowBits) ) | GetPackedNumber((incIntArray->lowBits),pos) );
}

void SetPackedIncreasingNumber(PackedIncreasingNumberArray *incIntArray, unsigned int pos, unsigned int num){
	unsigned char n;
	unsigned int mask, k;
	mask = (1UL << ((incIntArray->numTotalBits)-1)); // mask for the highest/leftmost bit of the number
	k = 0; // current position in the limits array
	n = (incIntArray->numHighBits); // current high bits count
	while( n > 0 ){
		if( num & mask ){ // bit 1, and advance to limit (2*(k+1))
			k++;
			k <<= 1;
		} else { // bit 0, and advance to limit (2*k+1)
			((incIntArray->highLimits)[k])++; // increase limit count
			k <<= 1;
			k++;
		}
		n--;
		mask >>= 1; // process next bit at the right
	}
	SetPackedNumber((incIntArray->lowBits),pos,(num & (incIntArray->lowBitsMask)));
}


// variables used by induced sort functions
int *nextId;
int *topLMS, *bottomLML;
int *firstId, *lastId;
unsigned int *suffixArray;
int maxAlphabetSize;

// variable used to store final BWT of the initial text
PackedNumberArray *bwtArray;

//#define GETCHAR(text,n,level) ( (level) ? ( ((int *)text)[n] ) : ( (int)letterIds[ (int)((char *)text)[n] ] ) )

#define SETBIT(bitarray,n) ( bitarray[(n>>5)] |= ( 1UL << (n&31) ) )
#define GETBIT(bitarray,n) ( bitarray[(n>>5)] & ( 1UL << (n&31) ) )
#define GETCHARTYPE(bitarray,n) ( (n==0) ? ( GETBIT(bitarray,n)?'s':'L' ) : ( GETBIT(bitarray,n) ? (GETBIT(bitarray,(n-1))?'s':'S') : (GETBIT(bitarray,(n-1))?'L':'l') ) )

// S-type: bit 1 at pos not 0 and with bit 0 at pos (n-1)
#define IS_S_TYPE(bitarray,n) ( GETBIT(bitarray,n) && ( (n!=0) && !GETBIT(bitarray,(n-1)) ) )
/*
// s-type: bit 1 at pos 0 or with bit 1 at pos (n-1)
#define IS_s_TYPE(bitarray,n) ( GETBIT(bitarray,n) && ( (n==0) || GETBIT(bitarray,(n-1)) ) )
// L-type: bit 0 at pos 0 or with bit 1 at pos (n-1)
#define IS_L_TYPE(bitarray,n) ( !GETBIT(bitarray,n) && ( (n==0) || GETBIT(bitarray,(n-1)) ) )
// l-type: bit 0 at pos not 0 and with bit 0 at pos (n-1)
#define IS_l_TYPE(bitarray,n) ( !GETBIT(bitarray,n) && ( (n!=0) && !GETBIT(bitarray,(n-1)) ) )
*/

void PrintISState( char *header, PackedNumberArray *text , unsigned int *posLMS , unsigned int *charTypes , unsigned int textSize , int alphabetSize , unsigned int *usedLMS, int level ){
	unsigned int n;
	int i,j,k;
	if(topLMS[0]!=(-1)){ // S-type suffixes
		printf("\n%s\n",header);
		for(i=0;i<16;i++) putchar('-'); putchar('\n');
		k=0;
		for(j=0;j<alphabetSize;j++){
			i=topLMS[j];
			while(i!=(-1)){
				n=posLMS[i];
				printf("[%02d]{%02d} %c %02u %c",k,i,GETCHARTYPE(charTypes,n),n,((usedLMS && GETBIT(usedLMS,i))?'*':' '));
				if(level){
					printf("(%02u)",GetPackedNumber(text,n)); n++;
					if(n==textSize) n=0;
					while(GETCHARTYPE(charTypes,n)!='S'){ printf("(%02u)",GetPackedNumber(text,n)); n++; }
					printf("(%02u)",GetPackedNumber(text,n));
				} else {
					printf("%c",letterChars[GetPackedNumber(text,n)]); n++;
					if(n==textSize) n=0;
					while(GETCHARTYPE(charTypes,n)!='S'){ printf("%c",letterChars[GetPackedNumber(text,n)]); n++; }
					printf("%c",letterChars[GetPackedNumber(text,n)]);
				}
				putchar('\n');
				i=nextId[i];
				k++;
			}
			for(i=0;i<16;i++) putchar('-'); putchar('\n');
		}
	} else { // L-type suffixes
		printf("\n%s\n",header);
		for(i=0;i<16;i++) putchar('-'); putchar('\n');
		k=0;
		for(j=alphabetSize;j!=0;){
			j--;
			i=bottomLML[j];
			while(i!=(-1)){
				n=posLMS[i];
				printf("[%02d]{%02d} %c %02u %c",k,i,GETCHARTYPE(charTypes,n),n,((usedLMS && GETBIT(usedLMS,i))?'*':' '));
				if(level){
					printf("(%02u)",GetPackedNumber(text,n)); n++;
					if(n==textSize) n=0;
					while(GETCHARTYPE(charTypes,n)!='S'){ printf("(%02u)",GetPackedNumber(text,n)); n++; }
					printf("(%02u)",GetPackedNumber(text,n));
				} else {
					printf("%c",letterChars[GetPackedNumber(text,n)]); n++;
					if(n==textSize) n=0;
					while(GETCHARTYPE(charTypes,n)!='S'){ printf("%c",letterChars[GetPackedNumber(text,n)]); n++; }
					printf("%c",letterChars[GetPackedNumber(text,n)]);
				}
				putchar('\n');
				i=nextId[i];
				k++;
			}
			for(i=0;i<16;i++) putchar('-'); putchar('\n');
		}
	}
}

// NOTE: if fillSA is set, the suffixArray variable is filled too
// NOTE: if fillSA is set and the level is 0, the BWT (of the initial text) is created instead
void InducedSort( PackedNumberArray *text , unsigned int *posLMS , unsigned int *charTypes , unsigned int textSize , int alphabetSize , int fillSA, int level ){
	int i, j, k, m;
	unsigned int n;
	unsigned int *bucketSize, *bucketPointer;
	char processingType, currentType;
	k=0; // just so compiler does not complain
	bucketSize=NULL;
	bucketPointer=NULL;
	if( fillSA ){ // if we want to fill the suffix array while doing the induced sorting, set pointers to the beginning of the buckets
		bucketSize = (unsigned int *)malloc(alphabetSize*sizeof(unsigned int));
		for( j = 0 ; j < alphabetSize ; j++ ) bucketSize[j] = 0; // reset bucket size
		for( n = 0 ; n < textSize ; n++ ) bucketSize[ GetPackedNumber(text,n) ]++; // count number of each alphabet letter in the text (size of each bucket)
		bucketPointer = (unsigned int *)malloc(alphabetSize*sizeof(unsigned int)); // pointer to the location in the suffixArray that will be filled up next
		bucketPointer[0] = 0;
		for( j = 1 ; j < alphabetSize ; j++ ) bucketPointer[j] = ( bucketPointer[(j-1)] + bucketSize[(j-1)] ); // points to the first position in each bucket
	}
	for( j = 0 ; j < alphabetSize ; j++ ){ // reset pointers that will be used
		bottomLML[j] = (-1);
		firstId[j] = (-1);
		lastId[j] = (-1);
	}
	#ifdef DEBUG_BWTIS
	if( textSize<100 && level<5 && !fillSA ) PrintISState("[->S]",text,posLMS,charTypes,textSize,alphabetSize,NULL,level);
	#endif
	for( j = 0 ; j < alphabetSize ; j++ ){ // induce sort L-type chars from S-type chars
		processingType = 'l';
		i = firstId[j]; // process l-type chars at the top of the bucket
		L_from_S:
		while( i != (-1) ){
			posLMS[i]--; // go backwards
			n = posLMS[i]; // position in the text
			currentType = GETCHARTYPE(charTypes,n); // char type at this position
			m = GetPackedNumber(text,n); // char at this position
			if( fillSA ){ // set this suffix position at the correct position in the top of its bucket inside the suffix array
				if(level) suffixArray[ bucketPointer[m] ] = n;
				else { // at level 0, fill the BWT array
					if( n == 0 ) n = textSize;
					n--; // position behind the current one
					n = GetPackedNumber(text,n);
					SetPackedNumber(bwtArray,(bucketPointer[m]),n);
				}
				bucketPointer[m]++;
			}
			if( currentType == 'l' ){ // add to bottom of l-type list at the top of the bucket
				if( firstId[m] == (-1) ) firstId[m] = i;
				else nextId[ lastId[m] ] = i;
				k = nextId[i]; // save next id here because if n==j and it was the last id, the next ptr was updated in the line above
				nextId[i] = (-1);
				lastId[m] = i;
			} else if( currentType == 'L' ){ // add to bottom of L-type list at the top of the bucket
				k = nextId[i]; // save next id
				nextId[i] = bottomLML[m];
				bottomLML[m] = i;
			}
			i = k; // next id
		}
		if( processingType == 'l' ){
			processingType = 'S';
			i = topLMS[j]; // process S-type chars at the bottom of the bucket
			goto L_from_S;
		}
	}
	if( fillSA ){ // if we want to fill the suffix array while doing the induced sorting, set pointers to the end of the buckets
		bucketPointer[0] = bucketSize[0];
		for( j = 1 ; j < alphabetSize ; j++ ) bucketPointer[j] = ( bucketPointer[(j-1)] + bucketSize[j] ); // total number of positions at and before this bucket
		for( j = 0 ; j < alphabetSize ; j++ ) bucketPointer[j]--; // points to the last position in each bucket
	}
	for( j = 0 ; j < alphabetSize ; j++ ){ // reset pointers that will be used
		topLMS[j] = (-1);
		firstId[j] = (-1);
		lastId[j] = (-1);
	}
	#ifdef DEBUG_BWTIS
	if( textSize<100 && level<5 && !fillSA ) PrintISState("[S->L]",text,posLMS,charTypes,textSize,alphabetSize,NULL,level);
	#endif
	for( j = alphabetSize ; j != 0 ; ){ // induce sort S-type chars from L-type chars
		j--; // process from bottom to top
		processingType = 's';
		i = firstId[j]; // process s-type chars at the bottom of the bucket
		S_from_L:
		while( i != (-1) ){
			if( posLMS[i] == 0 ){ // first position in text
				posLMS[i] = (textSize-1); // last position in text (it is the '$' symbol, which is an S-type char)
				if( fillSA ){ // it corresponds to the position of the terminal symbol
					if(level) suffixArray[ bucketPointer[0] ] = posLMS[i];
					else { // at level 0, fill the BWT array
						n = ( posLMS[i] - 1 ); // position (textSize-2)
						n = GetPackedNumber(text,n);
						SetPackedNumber(bwtArray,(bucketPointer[0]),n);
					}
				}
				topLMS[0] = i; // it goes to the first and only position in the bucket of the 0-th alphabet char
				k = nextId[i]; // save next id
				nextId[i] = (-1);
				i = k;
				continue;
			}
			posLMS[i]--; // go backwards
			n = posLMS[i]; // position in the text
			currentType = GETCHARTYPE(charTypes,n); // char type at this position
			m = GetPackedNumber(text,n); // char at this position
			if( fillSA ){ // set this suffix position at the correct position in the bottom of its bucket inside the suffix array
				if(level) suffixArray[ bucketPointer[m] ] = n;
				else { // at level 0, fill the BWT array
					if( n == 0 ) n = textSize;
					n--; // position behind the current one
					n = GetPackedNumber(text,n);
					SetPackedNumber(bwtArray,(bucketPointer[m]),n);
				}
				bucketPointer[m]--;
			}
			if( currentType == 's' ){ // add to top of s-type list at the bottom of the bucket
				if( firstId[m] == (-1) ) firstId[m] = i;
				else nextId[ lastId[m] ] = i;
				k = nextId[i]; // save next id here because if n==j and it was the last id, the next ptr was updated in the line above
				nextId[i] = (-1);
				lastId[m] = i;
			} else if( currentType == 'S' ){ // add to top of S-type list at the bottom of the bucket
				k = nextId[i]; // save next id
				nextId[i] = topLMS[m];
				topLMS[m] = i;
			}
			i = k; // next id
		}
		if( processingType == 's' ){
			processingType = 'L';
			i = bottomLML[j]; // process L-type chars at the top of the bucket
			goto S_from_L;
		}
	}
	if( fillSA ){
		free(bucketSize);
		free(bucketPointer);
	}
	/*
	#ifdef DEBUG_BWTIS
	if( textSize<100 && level<5 && !fillSA ) PrintISState("[L->S]",text,posLMS,charType,textSize,alphabetSize,NULL,level);
	#endif
	*/
}

// TODO: re-enable check until next S-type char inclusive
// TODO: call recursion with used LMS's only
// TODO: use suffixArray of size (numLMS) instead of keeping pointers and directly update as the recursion output
void BWTIS( PackedNumberArray *text , unsigned int textSize , int alphabetSize , int level ){
	unsigned int *charTypes;
	unsigned int *posLMS;
	PackedNumberArray *orderArray;
	int numLMS;
	int numUniqueLMS;
	int i, j, k, m;
	unsigned int n;
	unsigned int posOne, posTwo;
	unsigned int mask;
	/*
	int numUsedLMS;
	unsigned int *repeatedLMS, *usedLMS, *sortedUsedLMS;
	*/


	printf("\n(%u:%d)",textSize,alphabetSize);fflush(stdout);
	n = ( (textSize >> 5) + 1); // number of 32bit words that fit inside the textSize
	charTypes = (unsigned int *)calloc(n,sizeof(unsigned int)); // bit array with 0 for L-type chars and 1 for S-type chars
	n = ( textSize - 1); // last position of the text
	m = ( n >> 5 ); // word position inside charTypes array for current text position
	mask = ( 1UL << ( n & 31 ) ); // bit mask for current position inside word
	charTypes[m] = mask; // last char of the text is S-type
	j = GetPackedNumber(text,n); // last char
	k = 1; // type of last char (1 for S-type and 0 for L-type)
	numLMS = 0; // current number of LMS (last char is S*-type but it will only be counted next)
	while( n != 0 ){ // process text in reverse from end to start
		n--; // previous position
		mask >>= 1; // shift mask
		if( mask == 0 ){ // multiple of 32
			m--; // previous word
			//charTypes[m] = 0UL; // reset word bits
			mask = ( 1UL << 31 ); // re-initialize mask
		}
		i = GetPackedNumber(text,n); // current char
		if( i < j ) k = 1; // lower suffix, S-type char
		else if ( i > j ){ // higher suffix, L-type char
			if( k ) numLMS++; // if last char was S-type and this one is L-type, last one was S*-type
			k = 0;
		} // else (i==j), so use last used type
		if( k ) charTypes[m] |= mask; // set bit for S-type char
		j = i; // current char will be previous char on next step
	}

	/*
	charType = (char *)malloc(textSize*sizeof(char));
	n = (textSize-1);
	currentType = 's';
	charType[n] = currentType;
	while( n != 0 ){ // mark each text position as s-type or l-type
		n--;
		if( GetPackedNumber(text,n) < GetPackedNumber(text,(n+1)) ) currentType = 's';
		else if( GetPackedNumber(text,n) > GetPackedNumber(text,(n+1)) ) currentType = 'l';
		charType[n] = currentType; // if( text[n] == text[n+1] )
	}
	charType[(textSize-1)] = 'S';
	numLMS = 1;
	for( n = (textSize-2) ; n >= 1; n-- ){ // mark S-type and L-type chars and count LMS's
		if( charType[n] == 's' ){
			if( charType[n-1] == 'l' ){ // if it's s-type and char at the left is l-type, it's S-type
				charType[n] = 'S';
				numLMS++;
			}
		} else { // charType[i] == 'l'
			if( charType[n-1] == 's' ) charType[n] = 'L'; // if it's l-type and char at the left is s-type, it's L-type
		}
	}
	if( charType[0] == 'l' ) charType[0] = 'L'; // if the first/leftmost position is l-type, it's L-type because the last position is S-type
	*/

	#ifdef DEBUG_BWTIS
	if( textSize<100 && level<5 ){
		printf("\nlevel = %d\n",level);
		printf("textSize = %u\n",textSize);
		printf("alphabetSize = %d\n",alphabetSize);
		if(level){
			for(n=0;n<textSize;n++) printf(" %02u ",n);
			printf("\n");
			for(n=0;n<textSize;n++) printf("(%02u)",GetPackedNumber(text,n));
			printf("\n");
			for(n=0;n<textSize;n++) printf(" %c  ",GETCHARTYPE(charTypes,n));
			printf("\n");
			for(n=0;n<textSize;n++) printf(" %c  ",((GETCHARTYPE(charTypes,n)=='S')?'*':' '));
			printf("\n");
		} else {
			for(n=0;n<textSize;n++) printf("%u",(n%10));
			printf("\n");
			for(n=0;n<textSize;n++) printf("%c",letterChars[GetPackedNumber(text,n)]);
			printf("\n");
			for(n=0;n<textSize;n++) printf("%c",GETCHARTYPE(charTypes,n));
			printf("\n");
			for(n=0;n<textSize;n++) printf("%c",((GETCHARTYPE(charTypes,n)=='S')?'*':' '));
			printf("\n");
		}
	}
	#endif
	posLMS = (unsigned int *)malloc(numLMS*sizeof(unsigned int)); // position in the text of each S-type char (sorted by text position) (actually it is initialized with the position of the last char of this LMS, which is the position of the next S-type char to the right)
	if( level == 0 ){ // initialize these arrays only once with their maximum size (at level 0)
		//suffixArray = (unsigned int *)malloc(textSize*sizeof(unsigned int)); // positions in the text of all the sorted suffixes
		suffixArray = (unsigned int *)malloc(numLMS*sizeof(unsigned int)); // positions in the text of all the sorted suffixes (the s.a. for the initial text will never be used, so the largest size is numLMS for level 1)
		nextId = (int *)malloc(numLMS*sizeof(int));
		maxAlphabetSize = alphabetSize;
		topLMS = (int *)malloc(alphabetSize*sizeof(int));
		bottomLML = (int *)malloc(alphabetSize*sizeof(int));
		firstId = (int *)malloc(alphabetSize*sizeof(int));
		lastId = (int *)malloc(alphabetSize*sizeof(int));
	} else if( alphabetSize > maxAlphabetSize ){ // in level 1, the alphabetSize will be equal to numUniqueLMS, possibly the maximum value it will take
		maxAlphabetSize = alphabetSize;
		topLMS = (int *)realloc(topLMS,alphabetSize*sizeof(int));
		bottomLML = (int *)realloc(bottomLML,alphabetSize*sizeof(int));
		firstId = (int *)realloc(firstId,alphabetSize*sizeof(int));
		lastId = (int *)realloc(lastId,alphabetSize*sizeof(int));
	}
	for( j = 0 ; j < alphabetSize ; j++ ) topLMS[j] = (-1); // initialize bucket pointers
	i = (numLMS-1); // current LMS id (goes from (numLMS-1) , 0 , ... , (numLMS-2) )
	for( n = 0 ; n < textSize ; n++ ){ // add S-type chars to the bottom of their corresponding buckets
		if( IS_S_TYPE(charTypes,n) ){ // charType[n] == 'S'
			j = GetPackedNumber(text,n); // alphabet char (bucket)
			nextId[i] = topLMS[j]; // new top item links to old top item bellow
			topLMS[j] = i; // this is the new top item
			posLMS[i] = n; // set LMS position in text
			i++; // next LMS
			if( i == numLMS ) i = 0; // the position of the i-th LMS is stored in id (i-1) because next we walk backwards until we get to the LMS behind ((i-1)-th)
		}
	}
	#ifdef DEBUG_BWTIS
	if( textSize<100 && level<5 ){
		printf("\nnumLMS = %d\n",numLMS);
		for(i=0;i<numLMS;i++) printf("%02d ",i);
		printf("\n");
		for(i=0;i<numLMS;i++) printf("%02u,",posLMS[i]);
		printf("\n");
	}
	#endif

	InducedSort(text,posLMS,charTypes,textSize,alphabetSize,0,level); // sort the LMS's lexicographically using induced sort
	
	/*
	repeatedLMS = NULL;
	usedLMS = NULL;
	sortedUsedLMS = NULL;
	numUsedLMS = 0;
	m = ( (numLMS >> 5) + 1); // number of 32bit words needed to store numLMS bits
	repeatedLMS = (unsigned int *)calloc(m,sizeof(unsigned int)); // marks the LMS's that are not unique
	usedLMS = (unsigned int *)calloc(m,sizeof(unsigned int)); // marks the LMS's that will be used in the recursive call (repeated or unique after a repeated one)
	SETBIT(usedLMS,topLMS[0]); // the last LMS (for the terminator char) is always used
	//sortedUsedLMS = (unsigned int *)calloc(m,sizeof(unsigned int)); // marks the used LMS's but when they are in their sorted order (as in the suffix array)
	*/

	orderArray = NewPackedNumberArray(numLMS,numLMS); // orders of sorted LMS's
	numUniqueLMS = 0; // current id for unique LMS's
	for( j = 0 ; j < alphabetSize ; j++ ){ // get order of each sorted LMS (and detect repetitions)
		i = topLMS[j]; // id of the first LMS for this letter bucket
		while( i != (-1) ){
			SetPackedNumber(orderArray,i,numUniqueLMS); // set the order of this LMS in the array (it's sorted by text position)
			k = nextId[i]; // next id bellow
			if( k == (-1) ){ // if there are no LMS bellow, there is nothing to compare
				numUniqueLMS++; // increase id count so next letter bucket starts with a new id
				break;
			}
			posOne = (posLMS[i] + 1); // compare this LMS with the one bellow (the 1st char we know that is the same)
			posTwo = (posLMS[k] + 1);
			if( ( posLMS[i+1] - posOne ) != ( posLMS[k+1] - posTwo ) ){ // if the two LMS's have different lengths, they can't be equal
				numUniqueLMS++;
				i = k;
				continue;
			}
			while( ( GetPackedNumber(text,posOne) == GetPackedNumber(text,posTwo) ) && ( !IS_S_TYPE(charTypes,posOne) ) ){ // compare the two LMS's chars until they differ or we reach the end
				posOne++;
				posTwo++;
			}
			if( IS_S_TYPE(charTypes,posOne) && IS_S_TYPE(charTypes,posTwo) ){ // if we reached the end of the two LMS's, they are equal (compared until their last S-type char inclusive)
				/*
			if( GetPackedNumber(text,posOne) == GetPackedNumber(text,posTwo) ){ // if we reached the end of the two LMS's, they are equal (compared until their last S-type char inclusive)
				// TODO: remove this line bellow
				SETBIT(repeatedLMS,i); // just here to count real number of repeated LMS
				SETBIT(repeatedLMS,k); // mark only the later LMS as repeated
				SETBIT(usedLMS,i); // mark both LMS's as used
				SETBIT(usedLMS,k);
				if( (++i) != numLMS ) SETBIT(usedLMS,i); // mark the LMS's in front of these ones as used too
				if( (++k) != numLMS ) SETBIT(usedLMS,k);
				i--; // back to their previous values
				k--;
				*/
			} else { // if the two LMS are different
				numUniqueLMS++; // increase count of unique LMS's
			}
			i = k; // next LMS
		}
	}
	#ifdef DEBUG_BWTIS
	if( textSize<100 && level<5 ){
		//PrintISState("[L->S]",text,posLMS,charTypes,textSize,alphabetSize,usedLMS,level);
		PrintISState("[L->S]",text,posLMS,charTypes,textSize,alphabetSize,NULL,level);
		printf("\nnumUniqueLMS = %u\n",numUniqueLMS);
		for(i=0;i<numLMS;i++) printf("%02u ",i);
		printf("\n");
		for(i=0;i<numLMS;i++) printf("%02u,",GetPackedNumber(orderArray,i));
		printf("\n");
	}
	#endif
	
	/*
	m = 0; // used
	n = 0; // repeated
	for(i=0;i<numLMS;i++){
		if( GETBIT(usedLMS,i) ) m++;
		if( GETBIT(repeatedLMS,i) ) n++;
	}
	printf("\t[#LMS=%d,#uLMS=%d,#rLMS=%d]",numLMS,m,n);fflush(stdout);
	*/
	
	if( numUniqueLMS != numLMS ){ // if there are still repeated LMS's, do a recursive call here
		/*
		numUsedLMS = 0;
		m = 0; // current order of used LMS's (repetitions are counted as having the same order)
		n = 0; // current global order for all the LMS's
		for( j = 0 ; j < alphabetSize ; j++ ){ // get the order but now only inside the used LMS's
			i = topLMS[j];
			while( i != (-1) ){
				if( GETBIT(usedLMS,i) ){
					SetPackedNumber(orderArray,i,m); // set the order of this used LMS
					if( !GETBIT(repeatedLMS,i) ) m++; // the first LMS of a group of identical LMS's is always marked as non repeated, and the following ones as repeated
					SETBIT(sortedUsedLMS,n);
					numUsedLMS++;
				} else SetPackedNumber(orderArray,i,n); // set the order of this not used LMS
				n++; // one more LMS
				i = nextId[i]; // next id bellow
			}
		}
		*/
		BWTIS( orderArray , numLMS , numUniqueLMS , (level+1) ); // get the suffix array for the text formed by these LMS's
		if( level == 0 ){ // the alphabet size at level 0 is only 6, almost certainly much smaller than the maximum size currently allocated
			if( alphabetSize < maxAlphabetSize ){
				maxAlphabetSize = alphabetSize;
				topLMS = (int *)realloc(topLMS,alphabetSize*sizeof(int));
				bottomLML = (int *)realloc(bottomLML,alphabetSize*sizeof(int));
				firstId = (int *)realloc(firstId,alphabetSize*sizeof(int));
				lastId = (int *)realloc(lastId,alphabetSize*sizeof(int));
			}
		}
		for( j = 0 ; j < alphabetSize ; j++ ) topLMS[j] = (-1); // reset the bucket pointers
		for( m = numLMS ; m != 0 ; ){ // use the suffix array to add the sorted LMS's (S-type) to the bottom of their corresponding buckets
			m--; // process the suffix array from bottom to top because the S-type suffixes are added this way
			i = suffixArray[m]; // current LMS id
			n = posLMS[i]; // position in the text corresponding to this LMS
			j = GetPackedNumber(text,n); // first char of this LMS (its bucket)
			nextId[i] = topLMS[j]; // the id bellow this one will be the current top
			topLMS[j] = i; // this is the new top
		}
		#ifdef DEBUG_BWTIS
		if( textSize<100 && level<5 ){
			printf("\nlevel = %d\n",level);
			for(i=0;i<numLMS;i++) printf("%02d ",i);
			printf("\n");
			for(i=0;i<numLMS;i++) printf("%02u,",suffixArray[i]);
			printf("\n");
			PrintISState("[(SA)->S]",text,posLMS,charTypes,textSize,alphabetSize,NULL,level);
		}
		#endif
	} // bellow this point, all the S-type suffixes (LMS's) are all correctly sorted and linked to each other
	
	/*
	free(repeatedLMS);
	free(usedLMS);
	//free(sortedUsedLMS);
	*/

	FreePackedNumberArray(orderArray); // it was the text of the level bellow
	
	if( level == 0 ){ // when finishing level 0, we will no longer need the suffix array, because the final result will be the BWT
		free(suffixArray);
		bwtArray = NewPackedNumberArray(textSize,alphabetSize); // bit array that will store the BWT of the text
	}

	InducedSort(text,posLMS,charTypes,textSize,alphabetSize,1,level); // induce sort to fill the suffix array for this text

	#ifdef DEBUG_BWTIS
	if( textSize<100 && level<5 ){
		if(level){
			printf("\n[S->(SA)]\n");
			for(i=0;i<(int)textSize;i++){
				n=suffixArray[i];
				printf("[%02u] %c %02u ",i,GETCHARTYPE(charTypes,n),n);
				if(level) while(n!=textSize){ printf("(%02u)",GetPackedNumber(text,n)); n++; }
				else while(n!=textSize){ printf("%c",letterChars[GetPackedNumber(text,n)]); n++; }
				putchar('\n');
			}
		} else {
			printf("\n[S->(BWT)]\n");
			for(i=0;i<(int)textSize;i++) printf("[%02u] %c\n",i,letterChars[GetPackedNumber(bwtArray,i)]);
		}
	}
	#endif

	free(posLMS); // the positions are all stored in the suffix array, so this is not needed anymore
	free(charTypes);

	if( level == 0 ){
		free(nextId);
		free(topLMS);
		free(bottomLML);
		free(firstId);
		free(lastId);
	}
	printf(".");fflush(stdout);
	//return suffixArray;
}

// TODO: change LoadReference function to use filename stored in index
// TODO: on multiple references, add special char to separate them
// TODO: use type double/float for variables used in statistics at the end
// TODO: use inverse of speciallettersmask to prevent extra "~" operation
// TODO: replace original letter*BitsMask arrays by ~inverseLetter*BitsMask
// TODO: move global variables that are not needed for searching, to inside this function
// TODO: initialize mask arrays instead of storing them explicitly
// TODO: in main functions, replace array fetches by defines with the shifts to see if it's faster
// TODO: check speed of FMI_LetterJump with only half 16bits counts
void FMI_BuildIndex(char *filename, int silentmode){
	FILE *textFile;
	char *text, *textPtr, c;
	unsigned int letterId, i, k, n;
	unsigned int bwtPos, samplePos, textPos;
	unsigned int progressCounter, progressStep;
	unsigned int letterCounts[6], letterStartPos[6];
	PackedNumberArray *textArray;
	IndexBlock *block;
	if(silentmode) silentmode=0;
	printf("> Opening reference genome file <%s> ... ",filename);
	fflush(stdout);
	textFile=fopen(filename,"r");
	if(textFile==NULL){
		printf("\n> ERROR: Cannot open file\n");
		exit(0);
	}
	c=fgetc(textFile);
	if(c!='>'){
		printf("\n> ERROR: Invalid FASTA file\n");
		exit(0);
	}
	ungetc(c,textFile);
	n=0; // number of sequences
	textSize=0; // counts all letters plus one sequence terminator symbol
	InitializeLetterIdsArray(); // initialize letter ids array
	for(i=0;i<6;i++) letterCounts[i]=0; // reset counts for each letter (ACGTN)
	c='\0';
	while(c!=EOF){
		c=fgetc(textFile);
		if(c=='>'){ // new sequence
			while(c!='\n' && c!=EOF) c=fgetc(textFile); // skip label
			/* // NOTE: enable to set one special char at the beginning of each new sequence
			letterCounts[0]++; // terminator char symbol
			textSize++;
			if(textSize==UINT_MAX) break; // prevent overflow
			*/
			n++;
		}
		if(c=='\n') continue;
		if( (c>=65 && c<=90) || (c>=97 && c<=122) ){ // all alphabet letters (uppercase and lowercase)
			letterId = (unsigned int)letterIds[(int)c]; // get id of this letter
			letterCounts[letterId]++; // increase count of this letter
			textSize++; // one more text char
			if(textSize==UINT_MAX) break; // prevent overflow
		}
	}
	if(n==0){
		printf("\n> ERROR: No sequences found\n");
		exit(0);
	}
	if(textSize==0){
		printf("\n> ERROR: No valid characters found\n");
		exit(0);
	}
	if(textSize==UINT_MAX){
		printf("\n> ERROR: Maximum reference size exceeded\n");
		exit(0);
	}
	if(n>1) { printf("(");PrintUnsignedNumber(n);printf(" sequences) "); }
	/*
	k=(letterCounts[0]-n); // number of invalid chars, not counting sequence terminators
	if(k>0) printf("(%u N's) ",k);
	printf("(%u basepairs) OK\n",(textSize-n)); // number of valid chars, not counting sequence terminators
	*/
	/**/
	k=letterCounts[1]; // number of invalid chars ('N's)
	if(k>0) { printf("(");PrintUnsignedNumber(k);printf(" N's) "); }
	printf("(");PrintUnsignedNumber(textSize);printf(" basepairs) OK\n"); // number of valid chars (not including sequence terminator)
	letterCounts[0]++; // terminator char
	textSize++; // count terminator char
	/**/
	printf("> Reading reference chars ... ");
	fflush(stdout);
	#ifdef DEBUG
	text=(char *)malloc((textSize+1)*sizeof(char)); // +1 for the string '\0'
	if(text==NULL){
		printf("\n> ERROR: Not enough memory\n");
		exit(0);
	}
	#else
	text=NULL;
	#endif
	textArray = NewPackedNumberArray(textSize,6); // bit array that stores all the chars of the text in packed bits form
	if(textArray==NULL){
		printf("\n> ERROR: Not enough memory\n");
		exit(0);
	}
	rewind(textFile);
	i=0;
	while((c=fgetc(textFile))!=EOF){ // get text chars
		if( (c>=65 && c<=90) || (c>=97 && c<=122) ){ // all alphabet letters
			letterId=(unsigned int)letterIds[(int)c];
			SetPackedNumber(textArray,i,letterId);
			#ifdef DEBUG
			if(c>=97 && c<=122) c=(char)(c-32); // convert to uppercase
			if(c!='A' && c!='C' && c!='G' && c!='T') c='N'; // invalid letter
			text[i]=c;
			if((c=letterChars[GetPackedNumber(textArray,i)])!=text[i]){
				printf("\n> ERROR: textArray[%u](%u)('%c')=!=text[%u](%u)('%c')\n",i,GetPackedNumber(textArray,i),c,i,letterId,text[i]);
				exit(0);
			}
			#endif
			i++;
		} else if(c=='>'){ // new sequence
			while(c!='\n' && c!=EOF) c=fgetc(textFile); // skip label
			/*
			text[i++]='$'; // terminator char at the beggining of each sequence
			*/
		}
	}
	SetPackedNumber(textArray,i,0); // set terminator char (id=0) at the end
	#ifdef DEBUG
	text[i++]='$';
	text[i]='\0';
	#endif
	printf("OK\n");
	/* // NOTE: code to process text chars in the forward direction from the file without storing them in memory
	rewind(textFile);
	c='\0';
	while(c!=EOF){ // process all text chars
		c=fgetc(textFile); // get next char
		if( (c>=65 && c<=90) || (c>=97 && c<=122) ){ // all letters
			if(c>=97 && c<=122) c=(char)(c-32); // convert to uppercase
			letterId=(unsigned int)letterIds[(int)(*textPtr)]; // get letter id
		} else if(c=='>'){ // new sequence
			while(c!='\n' && c!=EOF) c=fgetc(textFile); // skip label
			letterId=0; // end of sequence char
		} else continue; // new line or EOF or non-alpha char
	*/
	printf("> Building BWT ");
	fflush(stdout);
	BWTIS( textArray , textSize , 6 , 0 );
	FreePackedNumberArray(textArray); // the packed text array is not needed anymore, only the BWT array
	printf(" OK\n");
	printf("> Allocating memory space for index ... ");
	fflush(stdout);
	/*
	sampleIntervalShift = 2; // samples with 4 positions each
	sampleIntervalSize = ( 1UL << sampleIntervalShift );
	sampleIntervalMask = ( sampleIntervalSize - 1UL );
	sampleIntervalHalfSize = ( sampleIntervalSize >> 1 );
	firstLetterMask = 1UL;
	secondLetterMask = 2UL;
	lastLetterMask  = ( 1UL << (sampleIntervalSize-1) );
	firstHalfLettersMask = ( ( 1UL << (sampleIntervalSize/2) ) - 1UL );
	lastHalfLettersMask  = ( firstHalfLettersMask << (sampleIntervalSize/2) );
	*/
	numSamples = ( ( ( textSize - 1 ) >> sampleIntervalShift ) + 1 ); // blocks of 32 chars (if textSize was a multiple of 32 it needed +1 additional unused sample, because last used pos is (textSize-1))
	if( ((textSize-1) & sampleIntervalMask) != 0 ) numSamples++; // if the last used position does not fall over a sample (0-th position), add an extra sample after the end (to speed up letter counts on 1st char of pattern search)
	lastBwtPos = ((numSamples-1) << sampleIntervalShift); // so that when getting the sample for this position, if falls in the 0-th pos of the last sample (numSamples-1)
	Index=(IndexBlock *)malloc(numSamples*sizeof(IndexBlock)); // 
	if(Index==NULL){
		printf("\n> ERROR: Not enough memory\n");
		exit(0);
	}
	printf("(%u MB) ",( (unsigned int)((numSamples)*sizeof(IndexBlock))/1000000));
	printf("OK\n");
	i=1;
	while(filename[(i-1)]!='\0') i++; // get filename size
	textFilename=(char *)malloc(i*sizeof(char));
	i=0;
	while(filename[i]!='\0'){ // copy filename chars
		textFilename[i]=filename[i];
		i++;
	}
	textFilename[i]='\0';
	printf("> Collecting letter jump samples ");
	fflush(stdout);	
	letterStartPos[0]=0; // the terminator char is at the top (0-th) position of the BWT (but on the right)
	for(i=1;i<6;i++) letterStartPos[i]=(letterStartPos[(i-1)]+letterCounts[(i-1)]); // where previous letter starts plus number of previous letter occurrences
	for(i=1;i<6;i++) letterCounts[i] = (letterStartPos[i]-1); // initialize all letter jumps with the position before the start of the letter
	progressStep=(textSize/10);
	progressCounter=0;
	samplePos = 0; // start in top position of the BWT and go down
	for( n = 0 ; n < textSize ; n++ ){
		progressCounter++;
		if(progressCounter==progressStep){ // print progress dots
			printf(".");
			fflush(stdout);
			progressCounter=0;
		}
		letterId = GetPackedNumber(bwtArray,n);
		letterCounts[letterId]++;
		if( ( n & sampleIntervalMask ) == 0 ){ // if we are over a sample, store here the current letter counts
			block = &(Index[samplePos]);
			(block->bwtLowBits) = 0U; // reset block
			(block->bwtHighBits) = 0U;
			(block->specialLettersMask) = 0U;
			for(i=1;i<6;i++) (block->letterJumpsSample[(i-1)]) = letterCounts[i]; // the i-th letter here is the (i-1)-th letter in the index
			samplePos++;
		}
		SetCharAtBWTPos(letterId,n); // copy the current letter from the packed BWT to the BWT in the index
	}
	if( ((textSize-1) & sampleIntervalMask) != 0 ){ // if one last extra sample was added at the end, fill it too
		block = &(Index[samplePos]);
		(block->bwtLowBits) = 0U; // reset block
		(block->bwtHighBits) = 0U;
		(block->specialLettersMask) = (~0U); // set special letters block to prevent confusing with real letters
		for(i=1;i<6;i++) (block->letterJumpsSample[(i-1)]) = letterCounts[i];
	}
	FreePackedNumberArray(bwtArray); // the packed BWT array is not needed anymore
	printf(" OK\n");
	printf("> Collecting position samples ");
	fflush(stdout);
	textPos=(textSize-1); // start at the terminator char (first text position is 0)
	n=0; // start at the first/topmost BWT position
	progressStep=(textSize/10);
	progressCounter=0;
	while(1){
		progressCounter++;
		if(progressCounter==progressStep){ // print progress dots
			printf(".");
			fflush(stdout);
			progressCounter=0;
		}
		if( ( n & sampleIntervalMask ) == 0 ){ // if we are over a sample, store here the current position of the text
			samplePos = ( n >> sampleIntervalShift );
			(Index[samplePos].textPositionSample) = textPos;
		}
		if(textPos==0) break;
		i = GetCharIdAtBWTPos(n); // get char at this BWT position (in the left)
		n = FMI_LetterJump(i,n); // follow the letter backwards to go to next position in the BWT
		textPos--;
	}
	printf(" OK\n");
	fflush(stdout);
	bwtPos=0; // just so compiler does not complain
	textPtr=NULL;
	#ifdef DEBUG
	if(textSize<100) PrintBWT(text,letterStartPos);
	textPtr = AppendToBasename(filename,".fmi"); // temporarily store index filename
	FMI_SaveIndex(textPtr); // save and load index
	FMI_FreeIndex();
	FMI_LoadIndex(textPtr);
	free(textPtr);
	printf("> Checking BWT sort ");
	fflush(stdout);
	progressCounter=0;
	k = FMI_PositionInText(0); // position in the text of the top BWT position (should be equal to (textSize-1))
	for(bwtPos=1;bwtPos<textSize;bwtPos++){ // compare current position with position above
		progressCounter++;
		if(progressCounter==progressStep){ // print progress dots
			printf(".");
			fflush(stdout);
			progressCounter=0;
		}
		i = FMI_PositionInText(bwtPos); // position in the text of the suffix in this row
		n = 0; // current suffix depth
		while( text[(k+n)] == text[(i+n)] ) n++; // keep following suffix chars to the right while their letters are equal
		if( letterIds[(int)(text[(k+n)])] > letterIds[(int)(text[(i+n)])] ) break; // if the top letter is larger than the bottom letter, it is incorrectly sorted
		k = i; // current pos will be prev pos in next step
	}
	if(bwtPos!=textSize) printf(" FAILED (error at BWT position %u)\n",bwtPos);
	else printf(" OK\n");
	printf("> Checking jump samples ");
	fflush(stdout);
	for(i=0;i<6;i++){ // initialize counters
		if(letterStartPos[i]!=0) letterCounts[i]=(letterStartPos[i]-1);
		else letterCounts[i]=0;
	}
	progressCounter=0;
	for(n=0;n<textSize;n++){
		progressCounter++;
		if(progressCounter==progressStep){ // print progress dots
			printf(".");
			fflush(stdout);
			progressCounter=0;
		}
		i=GetCharIdAtBWTPos(n); // get letter and update count
		if(i!=0){ // check if it is the terminator char because we cannot jump by it
			letterCounts[i]++;
			if( FMI_LetterJump(i,n) != letterCounts[i] ) break;
		}
		if( ( n & sampleIntervalMask ) == 0 ){ // if there is a sample at this position, check counts
			samplePos = ( n >> sampleIntervalShift );
			for(i=1;i<6;i++) if( (Index[samplePos].letterJumpsSample[(i-1)]) != (letterCounts[i]) ) break;
			if(i!=6) break;
		}
	}
	if(n!=textSize) printf(" FAILED (error at BWT position %u)\n",n);
	else printf(" OK\n");
	printf("> Checking position samples ");
	fflush(stdout);
	progressCounter=0;
	textPos=(textSize-1); // start at the terminator char
	n=0; // start at the first/topmost BWT position
	while(1){
		progressCounter++;
		if(progressCounter==progressStep){ // print progress dots
			printf(".");
			fflush(stdout);
			progressCounter=0;
		}
		if( FMI_PositionInText(n) != textPos ) break;
		if( ( n & sampleIntervalMask ) == 0 ){ // if there is a sample at this BWT position, check text position
			samplePos = ( n >> sampleIntervalShift );
			if( (Index[samplePos].textPositionSample) != textPos ) break;
		}
		if(textPos==0) break;
		textPos--;
		c=text[textPos];
		letterId=letterIds[(int)c];
		i=GetCharIdAtBWTPos(n);
		if(letterId!=i) break; // check if it is the same letter at the BWT and at the text
		if(i==0) break; // check if it is the terminator char because it should not be here
		n = FMI_LetterJump(i,n); // follow the letter backwards to go to next position in the BWT
	}
	if( textPos!=0 || FMI_PositionInText(n)!=0 || letterId!=i || i==0 || GetCharIdAtBWTPos(n)!=0 ) printf(" FAILED (error at text position %u)\n",textPos);
	else printf(" OK\n");
	/*
	printf("> Calculating BWT statistics ");
	fflush(stdout);
	k = 0; // number of runs of equal letters
	i = 0; // length of the current run
	prevLetterId = 0; // letter in the previous position
	progressCounter=0;
	for(bwtPos=0;bwtPos<textSize;bwtPos++){
		progressCounter++;
		if(progressCounter==progressStep){
			printf(".");
			fflush(stdout);
			progressCounter=0;
		}
		if( !(bwtPos & sampleIntervalMask) ) text[ Index[(bwtPos >> sampleIntervalShift)].textPositionSample ] = '.'; // mark sampled text positions with a dot
		letterId = GetCharIdAtBWTPos(bwtPos);
		if( (letterId==prevLetterId) && (i<255) ) i++;
		else {
			k++;
			i=1;
			prevLetterId=letterId;
		}
	}
	printf(" (BWT=");PrintUnsignedNumber((textSize*3)/8+(textSize/32)*5*4);
	printf("b;RLE=");PrintUnsignedNumber((k*(8+3))/8+(k/32)*5*4);printf("b)");
	fflush(stdout);
	bwtPos=0; // maximum distance (between 2 consecutive sampled text positions)
	n=0; // current average distance
	k=0; // number of points measured
	i=0; // distance from the last observed point
	for(textPos=0;textPos<textSize;textPos++){
		if( text[textPos]=='.' ){
			if(i>bwtPos) bwtPos=i; // check if this is a new maximum
			n = ( n*k + i/( ((i>n)?(i-n):(n-i)) + 1 ) ) / (k+1); // update weighted average: (current_average*number_of_prev_points + this_point/++distance_of_this_point_from_average)/++number_of_prev_points
			k++; // one more observed sample
			i=0; // reset distance
		}
		i++; // increase distance from last observed sample
	}
	printf(" (max=");PrintUnsignedNumber(bwtPos);
	printf(";wavg=");PrintUnsignedNumber(n);printf(")");
	printf(" OK\n");
	fflush(stdout);
	*/
	free(text);
	#endif
	fclose(textFile); // and the index memory will be freed by the function calling this one
}
