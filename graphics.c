#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <limits.h>
#include <string.h>
#include "bitmap.h"
#include "tools.h"

#ifdef _MSC_VER
	#pragma warning(disable:4996)
#endif

#define CHAR_WIDTH 5
#define CHAR_HEIGHT 6
#define CHAR_SPACE 2
#define ALPHABET_SIZE 95
const unsigned int alphabet[] = { // least significant bit is the upper left pixel, and so on
	0x00000000, // (space)
	0x08021084, // !
	0x0000014A, // "
	0x15F57D40, // #
	0x1F4717C4, // $
	0x22221111, // %
	0x2C9A88A2, // &
	0x00000084, // '
	0x1821084C, // (
	0x0C842106, // )
	0x2AEFBAA0, // *
	0x084F9080, // +
	0x04400000, // ,
	0x00007C00, // -
	0x08000000, // .
	0x02221110, // /
	0x1D3AD72E, // 0
	0x3E4214C4, // 1
	0x3E17420F, // 2
	0x1F083A1F, // 3
	0x108FA988, // 4
	0x1F083C3F, // 5
	0x1D18BC3E, // 6
	0x0222221F, // 7
	0x1D18BA2E, // 8
	0x1F0F462E, // 9
	0x00400080, // :
	0x04401000, // ;
	0x30609B00, // <
	0x01F07C00, // =
	0x06C83060, // >
	0x0802322E, // ?
	0x05DAF62E, // @
	0x23153944, // A
	0x1F18BE2F, // B
	0x3C10843E, // C
	0x1F18C62F, // D
	0x3E109C3F, // E
	0x02109C3F, // F
	0x3D1C843E, // G
	0x2318FE31, // H
	0x3E42109F, // I
	0x0C94A11E, // J
	0x22928CB9, // K
	0x3E108421, // L
	0x2318D771, // M
	0x239AD671, // N
	0x1D18C62E, // O
	0x0217C62F, // P
	0x2C9AC62E, // Q
	0x2297C62F, // R
	0x1F08383E, // S
	0x0842109F, // T
	0x1D18C631, // U
	0x08A54631, // V
	0x23BAC631, // W
	0x23151151, // X
	0x08422A31, // Y
	0x3E22111F, // Z
	0x1C21084E, // [
	0x20821041, // (backslash)
	0x1C84210E, // ]
	0x00004544, // ^
	0x3E000000, // _
	0x00000082, // `
	0x1CA721C0, // a
	0x1CA53842, // b
	0x1C2109C0, // c
	0x1CA53908, // d
	0x1C2729C0, // e
	0x0842388C, // f
	0x1C8729C0, // g
	0x14A53842, // h
	0x08421004, // i
	0x0C421004, // j
	0x14A32842, // k
	0x08421084, // l
	0x2B5AD560, // m
	0x14A528C0, // n
	0x1CA529C0, // o
	0x042729C0, // p
	0x108729C0, // q
	0x04253840, // r
	0x1C8709C0, // s
	0x0C4211C4, // t
	0x1CA52940, // u
	0x08452940, // v
	0x14AAD6A0, // w
	0x14A22940, // x
	0x1C872940, // y
	0x1C2221C0, // z
	0x18210C4C, // {
	0x08421084, // |
	0x0C846106, // }
	0x009AC800, // ~
	0x3FFFFFFF // invalid
};

/*
const unsigned int dna_alphabet[] = {
	0x23153944, // A
	0x3C10843E, // C
	0x3D1C843E, // G
	0x0842109F, // T
	0x00007C00, // -
	0xFFFFFFFF  // invalid
};

const unsigned int numbers_alphabet[] = {
	0x1D3AD72E, // 0
	0x3E4214C4, // 1
	0x3E17420F, // 2
	0x1F083A1F, // 3
	0x108FA988, // 4
	0x1F083C3F, // 5
	0x1D18BC3E, // 6
	0x0222221F, // 7
	0x1D18BA2E, // 8
	0x1F0F462E  // 9
};
*/

#define WHITE getColorFromPalette(255,255,255)
#define BLACK getColorFromPalette(0,0,0)
#define GREY getColorFromPalette(128,128,128)
#define RED getColorFromPalette(255,0,0)
#define GREEN getColorFromPalette(0,255,0)
#define BLUE getColorFromPalette(0,0,255)

void DrawChar(char c, int x, int y, uint8_t color){
	unsigned int charpixels;
	unsigned int i, j;
	c-=32; // convert char code to position in alphabet array
	if(c>ALPHABET_SIZE) c=ALPHABET_SIZE; // invalid char
	charpixels=alphabet[(int)c];
	j=CHAR_HEIGHT;
	while(j--){ // all rows
		i=CHAR_WIDTH;
		while(i--){ // all columns of each row
			if(charpixels & 1u) drawPoint(x,y,color);
			charpixels>>=1;
			x++;
		}
		x-=CHAR_WIDTH;
		y++;
	} // draw from left to right, and top to bottom
}

void DrawString(char *s, int x, int y, uint8_t color){
	while((*s)!='\0'){
		DrawChar((*s),x,y,color);
		x+=(CHAR_WIDTH+CHAR_SPACE);
		s++;
	}
}

void TestAlphabet(){
	char alphabetstring[(ALPHABET_SIZE+2)];
	int i;
	for(i=0;i<=ALPHABET_SIZE;i++) alphabetstring[i]=(char)(i+32); // all printable ASCII chars
	alphabetstring[i]='\0';
	DrawString(alphabetstring,1,1,BLACK);
}

int DigitsCount(int number){
	int count;
	count=1;
	while(number>=10){
		number /= 10;
		count++;
	}
	return count;
}

void DrawNumber(int number, int x, int y, uint8_t color){
	int n,d,m;
	n=DigitsCount(number); // print number centered
	n=( n*CHAR_WIDTH + (n-1)*CHAR_SPACE )/2;
	if(x>n) x-=n;
	else x=1;
	n=number;
	d=1;
	while(n>=10){
		n /= 10;
		d *= 10;
	}
	n=number;
	while(d>0){
		m = (n/d); // from 0 to 9
		DrawChar((char)(48+m),x,y,color); // draw char corresponding to this digit
		x+=(CHAR_WIDTH+CHAR_SPACE);
		n -= (m*d);
		d /= 10;
	}
}

void DrawHorizontalLine(int x, int y, int length, uint8_t color){
	while(length--) drawPoint(x++,y,color);
}

void DrawVerticalLine(int x, int y, int length, uint8_t color){
	while(length--) drawPoint(x,y++,color);
}

void DrawEmptyRectangle(int x, int y, int xsize, int ysize, uint8_t color){
	DrawHorizontalLine(x,y,xsize,color);
	DrawHorizontalLine(x,y+ysize,xsize,color);
	DrawVerticalLine(x,y,ysize,color);
	DrawVerticalLine(x+xsize,y,ysize,color);
}

void DrawFilledRectangle(int x, int y, int xsize, int ysize, uint8_t color){
	while(ysize--) DrawHorizontalLine(x,y++,xsize,color);
}


void CreateGlobalMap(unsigned int startpos, unsigned int endpos, unsigned int numpos, char *readsfilename){
	FILE *file;
	unsigned char *poscount;
	unsigned int *rowcount;
	unsigned int i, j, k, n, pos, numrows, count, totalnumchars, numfilledpos;
	unsigned int imgheight, imgwidth, refwidth, refheight, rowheight, margindim;
	unsigned int readstartpos, readendpos, readsize, numreads, totalnumreads;
	double posperpixel;
	unsigned int x, y, startx, starty;
	uint8_t color;
	char c;
	int nret = 0; nret = (int)nret;
	printf("> Allocating space for positions from %u to %u ... ",(startpos+1),(endpos+1)); // +1 for 1-based positions
	fflush(stdout);
	poscount=(unsigned char *)calloc(numpos,sizeof(unsigned char));
	if(poscount==NULL){
		printf("\n> ERROR: Not enought memory\n");
		exit(0);
	}
	printf("OK (%u positions)\n",numpos);
	printf("> Processing mapped reads file <%s> ... ",readsfilename);
	fflush(stdout);
	if((file=fopen(readsfilename,"r"))==NULL){
		printf("\n> ERROR: File not found\n");
		exit(0);
	}
	c=fgetc(file);
	while(c=='@'){ // skip header lines
		while(c!='\n' && c!=EOF) c=fgetc(file);
		c=fgetc(file);
	}
	ungetc(c,file);
	numreads=0;
	totalnumreads=0;
	totalnumchars=0;
	while(c!=EOF){
		for(i=0;i<2;i++){ // advance to 2nd tab
			c=fgetc(file);
			while(c!='\t' && c!=EOF) c=fgetc(file);
		}
		c=fgetc(file); // get first letter of reference
		if(c=='*'){ // if the read was not mapped, advance to the next one
			while(c!='\n' && c!=EOF) c=fgetc(file);
			continue;
		}
		for(i=0;i<1;i++){ // advance one more tab
			c=fgetc(file);
			while(c!='\t' && c!=EOF) c=fgetc(file);
		}
		readstartpos=0;
		nret = fscanf(file,"%u",&readstartpos); // get read mapped position
		readstartpos--; // fix 1-based location
		for(i=0;i<6;i++){ // advance 6 more tabs
			c=fgetc(file);
			while(c!='\t' && c!=EOF) c=fgetc(file);
		}
		readsize=0;
		while((c=fgetc(file))!='\t' && c!=EOF) readsize++; // get size of the read
		if(c==EOF) break;
		readendpos=(readstartpos+readsize-1);
		if( (readstartpos>=startpos && readstartpos<=endpos) || (readendpos>=startpos && readendpos<=endpos) ){ // if the read falls in the valid range of positions
			if(readstartpos>=startpos){ // if the read starts in the middle of the valid interval
				i=0; // current position in the read
				pos=(readstartpos-startpos); // current position in the array of counts
			} else { // if the read starts before the first valid position, advance to it
				i=(startpos-readstartpos);
				pos=0;
			}
			while(i<readsize && pos<numpos){
				if(poscount[pos]<255) poscount[pos]++; // do not exceed the maximum capacity of the counters
				pos++;
				i++;
			}
			numreads++; // count valid reads
			totalnumchars+=readsize; // count total number of mapped bases
		}
		totalnumreads++; // count total mapped reads
		while(c!='\n' && c!=EOF) c=fgetc(file); // advance to next read
	}
	fclose(file);
	if(numreads==0){
		printf("\n> ERROR: No valid reads found\n");
		exit(0);
	}
	printf("OK (%u reads - %u inside range)\n",totalnumreads,numreads);
	printf("> Creating global reference coverage plot ... ");
	fflush(stdout);
	numfilledpos=0;
	numrows=0;
	for(i=0;i<numpos;i++){
		if(poscount[i]>numrows) numrows=poscount[i]; // find maximum number of rows
		if(poscount[i]!=0) numfilledpos++; // count number of filled positions
	}
	if(numrows>50) numrows=50; // limit maximum number of rows to 50
	rowcount=(unsigned int *)calloc(numrows,sizeof(unsigned int));
	rowheight=4;
	margindim=(2*rowheight); // side and top margins with size 2 times the row thickness
	refheight=(2*rowheight); // height of the reference in pixels
	imgwidth=1024; // fixed image width
	imgheight=( 2*margindim + 2*CHAR_HEIGHT + refheight + 2*numrows*rowheight ); // top margins, position markers, reference, and all rows with a spacing before
	initializeBitmap((int)imgwidth,(int)imgheight,2);
	k=DigitsCount(endpos+1); // make room to print half of a position number on both ends
	refwidth=( imgwidth - 2*( margindim + k*(CHAR_WIDTH+CHAR_SPACE)/2 ) ); // width of the reference in pixels
	posperpixel=(((double)numpos)/((double)refwidth)); // number of ref positions represented by each pixel
	startx=( margindim + (k*(CHAR_WIDTH+CHAR_SPACE))/2 ); // where to start drawing from
	starty=margindim;
	color=BLACK;
	x=startx;
	y=starty;
	DrawNumber((startpos+1),(x-1),y,color); // print left pos number
	DrawVerticalLine((x-1),(y+CHAR_HEIGHT+CHAR_HEIGHT/2),(CHAR_HEIGHT/2),color);
	x=(imgwidth-startx);
	DrawNumber((endpos+1),x,y,color); // print right pos number
	DrawVerticalLine(x,(y+CHAR_HEIGHT+CHAR_HEIGHT/2),(CHAR_HEIGHT/2),color);
	j=(unsigned int)floor(posperpixel*(double)(refwidth-k*(CHAR_WIDTH+CHAR_SPACE))); // last middle position that can have a marked number
	k=( k*CHAR_WIDTH + (k-1)*CHAR_SPACE ); // space occupied by each number
	n=( refwidth - k ); // space available to print numbers without the numbers already at both ends
	k+=CHAR_WIDTH; // each number needs a space after it
	k=(n/k); // how many numbers can be printed in that space
	n=(endpos-startpos+1); // total number of positions
	n=(n/(k+1)); // how many positions can go between each mark (minimum)
	k=1; // how many positions between marks
	while(k<n){
		k=2*k; // x2
		if(k>=n) break;
		k=(k/2)*5; // x5
		if(k>=n) break;
		k=2*k; // x10
	}
	n=(unsigned int)floor(((double)k)/posperpixel); // number of pixels between marks
	x=(startx+n);
	y=starty;
	i=k;
	while(i<j){ // print marks
		DrawNumber((startpos+1+i),x,y,color);
		DrawVerticalLine(x,(y+CHAR_HEIGHT+CHAR_HEIGHT/2),(CHAR_HEIGHT/2),color);
		x+=n;
		i+=k;
	}
	k=DigitsCount(endpos+1);
	startx=( margindim + (k*(CHAR_WIDTH+CHAR_SPACE))/2 ); // margin + half the space occupied by the number of the 1st position
	starty=( margindim + 2*CHAR_HEIGHT ); // bellow position numbers and marks
	color=BLACK;
	DrawHorizontalLine((startx-1),(starty-1),(refwidth+2),color); // draw black box frame around reference
	DrawHorizontalLine((startx-1),(starty+refheight),(refwidth+2),color);
	DrawVerticalLine((startx-1),(starty-1),(refheight+2),color);
	DrawVerticalLine((startx+refwidth),(starty-1),(refheight+2),color);
	for(i=0;i<refwidth;i++){ // process all horizontal pixels along the reference
		startpos=(unsigned int)floor(((double)i)*posperpixel);
		endpos=(unsigned int)floor(((double)(i+1))*posperpixel);
		if(endpos>=numpos) endpos=(numpos-1); // prevent overflow in last pixel
		n=(endpos-startpos+1); // number of positions included in this pixel
		count=0; // get reference coverage
		for(pos=startpos;pos<=endpos;pos++){
			//count+=poscount[pos];
			if(poscount[pos]>0) count++; // count how many positions are filled in this pixel
			j=0; // fill counts for each row
			while(poscount[pos]>j && j<numrows) rowcount[j++]++; // if the position has chars on that row (fill from row 0 to row (poscount[pos]-1))
		}
		//color=(uint8_t)((255*count)/(n*numrows)); // get color intensity from filled ratio (full intensity = 255)
		color=(uint8_t)((255*count)/(n)); // get color intensity from occupancy ratio (full intensity = 255)
		color=getColorFromPalette(255,(255-color),(255-color)); // paint reference in shades of red
		for(k=0;k<refheight;k++) drawPoint((int)(startx+i),(int)(starty+k),color); // draw reference line thickness
		for(j=0;j<numrows;j++){ // process all vertical rows
			if(rowcount[j]==0) break; // no more filled rows next
			color=(uint8_t)((255*rowcount[j])/n); // get black color intensity from filled ratio (full intensity = 255)
			color=getColorFromPalette((255-color),(255-color),255); // paint in shades of blue
			x=(startx+i); // horizontal position of the pixel
			y=(starty+refheight+(2*j+1)*rowheight); // vertical position of the pixel in the current column
			for(k=0;k<rowheight;k++) drawPoint((int)x,(int)(y+k),color); // draw row thickness
			rowcount[j]=0; // reset row count for next pixel column
		}
	}
	printf("OK\n");
	printf("> Coverage: depth = %.2lfx ; fill = %.2lf%%\n", (((double)totalnumchars)/((double)numpos)) , ((((double)numfilledpos)*100.0)/((double)numpos)) );
	free(rowcount);
	free(poscount);
}

// TODO: add support for padding operation (concatenate strings)
void CreateLocalMap(unsigned int startpos, unsigned int endpos, unsigned int numpos, unsigned int refsize, char *reffilename, char *readsfilename, char*gfffilename){
	FILE *file, *gfffile;
	char *refchars, **chars, ***inschars, *string;
	unsigned short int **freeahead, **charcount;
	unsigned short int *inscount, *inssize, **insrows;
	unsigned int i, j, k, n, m, refpos, readpos, numrows, insidscount, *insids, count, numcharcountpos;
	unsigned int imgheight, imgwidth, rowheight, colwidth, margindim;
	unsigned int readstartpos, readsize, numreads, totalnumreads, maxreadsize, tempreadstartpos, tempreadendpos, tempreadsize;
	unsigned int x, y, startx, starty;
	unsigned int redvalue, greenvalue, addvalue;
	uint8_t color;
	char c, cigarcode, *read;
	int cigarcount, numcigarops, extrareadsize;
	fpos_t cigarstringpos, sequencepos;
	int nret = 0; nret = (int)nret;
	if(numpos>(unsigned int)USHRT_MAX){
		printf("\n> ERROR: Sequence range is too large to draw local mapping plot\n");
		exit(0);
	}
	gfffile=NULL;
	if(gfffilename!=NULL){
		printf("> Processing annotations file <%s> ... ",gfffilename);
		fflush(stdout);
		if((gfffile=fopen(gfffilename,"r"))==NULL){
			printf("\n> ERROR: File not found\n");
			exit(0);
		}
		printf("OK\n");
	}
	printf("> Allocating space for positions from %u to %u ... ",(startpos+1),(endpos+1)); // +1 because SAM positions start at 1 and these passed arguments start at 0
	fflush(stdout);
	numpos++; // +1 to account for insertions after the last position
	endpos++;
	refchars=(char *)calloc((numpos+1),sizeof(char));
	numrows=1;
	chars=(char **)calloc(numrows,sizeof(char *)); // stores pointers to rows of chars (reallocable with numrows)
	chars[0]=(char *)calloc(numpos,sizeof(char)); // each row stores numpos chars (initialized as '\0')
	freeahead=(unsigned short int **)calloc(numrows,sizeof(unsigned short int *)); // stores pointers to rows of counters (reallocable with numrows)
	freeahead[0]=(unsigned short int *)calloc(numpos,sizeof(unsigned short int)); // each row stores numpos counters (number of available positions to the right)
	insidscount=1; // current number of positions in the seq with insertions there
	insids=(unsigned int *)calloc((refsize+1),sizeof(unsigned int)); // ids of the insertion chars between (before) each char/position of the sequence (0 means it has no insertions)
	inscount=(unsigned short *)calloc(insidscount,sizeof(unsigned short)); // number of insertions for this id (array reallocable with insidscount ; position 0 is never used)
	inssize=(unsigned short *)calloc(insidscount,sizeof(unsigned short)); // maximum size of the insertions for this id (array reallocable with insidscount ; position 0 is never used)
	insrows=(unsigned short **)calloc(insidscount,sizeof(unsigned short *)); // which rows have insertions in this id (insrows[i] reallocable with inscount[i])
	inschars=(char ***)calloc(insidscount,sizeof(char **)); // strings of the inserted chars for each row of each id (inschars[i] reallocable with inscount[i])
	maxreadsize=1;
	read=(char *)calloc((maxreadsize+1),sizeof(char));
	if( refchars==NULL || chars==NULL || chars[0]==NULL || freeahead==NULL || freeahead[0]==NULL
		|| insids==NULL || inscount==NULL || inssize==NULL || insrows==NULL || inschars==NULL || read==NULL ){
		printf("\n> ERROR: Not enought memory\n");
		exit(0);
	}
	for(i=0;i<numpos;i++) freeahead[0][i]=USHRT_MAX; // initialize all positions of first row as free
	printf("OK (%u positions)\n",(numpos-1)); // -1 because at the beginning we added +1 to account for insertions at the end
	
	printf("> Reading reference genome positions ... ");
	fflush(stdout);
	if((file=fopen(reffilename,"r"))==NULL){
		printf("\n> ERROR: File not found\n");
		exit(0);
	}
	i=0; // position in stored chars vector
	k=0; // position in reference
	while((c=fgetc(file))!=EOF){ // get reference chars
		if(c=='>') while(c!='\n' && c!=EOF) c=fgetc(file); // skip description in case of multi-fasta file
		if(c>=97 && c<=122) c=(char)(c-32); // lowercase char, convert to uppercase
		if(c>=65 && c<=90){ // uppercase char
			if(c!='A' && c!='C' && c!='G' && c!='T') c='N'; // not ACGT char
			if(k>=startpos) refchars[i++]=c; // if the position is inside the valid range, store this char
			if(i==(numpos-1)) break; // stop when we have all valid position chars (-1 because we added +1 to numpos before)
			k++; // count valid chars (first valid char is at the 0-th position)
		}
	}
	fclose(file);
	printf("OK (%u chars)\n",i);
	printf("> Processing mapped reads file <%s> ... ",readsfilename);
	fflush(stdout);
	if((file=fopen(readsfilename,"r"))==NULL){
		printf("\n> ERROR: File not found\n");
		exit(0);
	}
	c=fgetc(file);
	while(c=='@'){ // skip header lines
		while(c!='\n' && c!=EOF) c=fgetc(file);
		c=fgetc(file);
	}
	ungetc(c,file);
	totalnumreads=0;
	numreads=0;
	while(c!=EOF){ // loop for all reads inside SAM file
		for(i=0;i<2;i++){ // advance to 2nd tab (reference label)
			c=fgetc(file);
			while(c!='\t' && c!=EOF) c=fgetc(file);
		}
		c=fgetc(file); // get first letter of reference
		if(c=='*' || c==EOF){ // if the read was not mapped, advance to the next one
			while(c!='\n' && c!=EOF) c=fgetc(file);
			continue;
		}
		for(i=0;i<1;i++){ // advance one more tab (position)
			c=fgetc(file);
			while(c!='\t' && c!=EOF) c=fgetc(file);
		}
		readstartpos=0;
		nret = fscanf(file,"%u",&readstartpos); // get read mapped position
		if(readstartpos!=0) readstartpos--; // fix 1-based location
		for(i=0;i<2;i++){ // advance 2 more tabs (CIGAR string)
			c=fgetc(file);
			while(c!='\t' && c!=EOF) c=fgetc(file);
		}
		c=fgetc(file); // get first letter of the CIGAR string
		if(c=='*' || c==EOF){ // if the CIGAR string is not available, advance to the next read
			while(c!='\n' && c!=EOF) c=fgetc(file);
			continue;
		}
		ungetc(c,file);
		fgetpos(file,&cigarstringpos); // save location
		tempreadstartpos=readstartpos; // position where characters will effectively start being written
		extrareadsize=0;
		numcigarops=0;
		while(fscanf(file,"%d%c",&cigarcount,&cigarcode)>0){ // process CIGAR string to track alterations to the read's mapped size in the arrays
			numcigarops++;
			if(cigarcode=='I') extrareadsize-=cigarcount; // insertion: chars added to special array, not to regular one
			else if(cigarcode=='D' || cigarcode=='N') extrareadsize+=cigarcount; // deletion: dash characters are added
			else if(cigarcode=='S'){ // soft clip: chars are present but dots are written instead
				if(numcigarops==1){ // if the chars are clipped at the beginning, there has to be space for them
					if(tempreadstartpos>=(unsigned int)cigarcount) tempreadstartpos-=cigarcount;
					else tempreadstartpos=0; // prevent underflow
				}
				numcigarops--; // do not count this operation
			} else if(cigarcode=='H'){ // hard clip: chars are not present
				numcigarops--;
			}
		}
		for(i=0;i<3;i++){ // advance 3 more tabs (sequence) (should be 4 but fscanf already consumes one '\t')
			c=fgetc(file);
			while(c!='\t' && c!=EOF) c=fgetc(file);
		}
		fgetpos(file,&sequencepos); // save location
		readsize=0;
		while((c=fgetc(file))!='\t' && c!=EOF) readsize++; // get size of the read
		if(c==EOF) break;
		if(extrareadsize<0 && (-extrareadsize)>((int)readsize)) tempreadsize=0; // prevent underflow
		else tempreadsize=(unsigned int)(((int)readsize)+extrareadsize); // number of positions that will be effectively written to the arrays
		tempreadendpos=(tempreadstartpos+tempreadsize-1);
		totalnumreads++; // count total number of mapped reads
		if( ( tempreadstartpos>endpos || tempreadendpos<startpos ) || tempreadsize==0 ){ // if the read does not fall in the valid range of positions, continue to next one
			while(c!='\n' && c!=EOF) c=fgetc(file);
			continue;
		}
		if(tempreadstartpos<startpos){ // if the read starts before the first valid position in the arrays, update read size and start position
			tempreadsize-=(startpos-tempreadstartpos);
			tempreadstartpos=startpos;
		}
		if(tempreadendpos>endpos){ // if the read ends after the last valid position in the arrays, update read size and end position
			tempreadsize-=(tempreadendpos-endpos);
			tempreadendpos=endpos;
		}
		if(tempreadstartpos!=startpos){ // require an empty space on both ends of the read (so it does not appear joined with another read)
			tempreadstartpos--;
			tempreadsize++;
		}
		if(tempreadendpos!=endpos){
			tempreadendpos++;
			tempreadsize++;
		}
		if(tempreadsize>numpos) tempreadsize=numpos;
		j=(tempreadstartpos-startpos); // starting position (in the arrays) were we will start looking for free space
		k=0; // row where the chars will be stored
		while( k<numrows && freeahead[k][j]<((unsigned short int)tempreadsize)) k++; // find the first row that has enought free space to store the read chars
		if(k==numrows){ // if none of the existing rows has the necessary free space, create a new row
			if(numrows==USHRT_MAX){ // more rows than we can store
				/*
				printf("\n> ERROR: More stacked reads than the arrays can store\n");
				exit(0);
				*/
				numreads++; // it is a valid read too
				while(c!='\n' && c!=EOF) c=fgetc(file); // go to next read
				continue;
			}
			numrows++;
			freeahead=(unsigned short int **)realloc(freeahead,numrows*sizeof(unsigned short int *));
			freeahead[(numrows-1)]=(unsigned short int *)calloc(numpos,sizeof(unsigned short int));
			for(i=0;i<numpos;i++) freeahead[(numrows-1)][i]=USHRT_MAX;
			chars=(char **)realloc(chars,numrows*sizeof(char *));
			chars[(numrows-1)]=(char *)calloc(numpos,sizeof(char));
		}
		if(readsize>maxreadsize){ // realloc space for read chars if needed
			maxreadsize=readsize;
			read=(char *)realloc(read,(maxreadsize+1)*sizeof(char));
		}
		fsetpos(file,&sequencepos); // restore location to beginning of sequence chars
		readsize=0;
		while((c=fgetc(file))!='\t' && c!=EOF){
			if( c>=97 && c<=122 ) c=(char)(c-32); // convert to uppercase if needed
			read[readsize++]=c; // store read chars
		}
		read[readsize]='\0';
		fsetpos(file,&cigarstringpos); // restore location to beginning of CIGAR string
		numcigarops=0;
		refpos=readstartpos; // current position relative to the full reference
		readpos=0; // current position in the read
		while(fscanf(file,"%d%c",&cigarcount,&cigarcode)>0){ // process CIGAR string
			numcigarops++;
			if(numcigarops==1 && cigarcode=='S'){ // if the first operation is a soft clipping, make room for the clipped chars behind
				if(refpos>=(unsigned int)cigarcount) refpos-=cigarcount;
				else { // prevent underflow, and discard chars that go beyond the beginning of the ref
					refpos=0;
					n=(cigarcount-refpos); // number of discarded positions
					readpos+=n;
					cigarcount-=n;
				}
				numcigarops--;
			}
			if(refpos<startpos){ // operation starts before the start of the valid interval
				if((refpos+cigarcount)<startpos || cigarcode=='I' || cigarcode=='P' || cigarcode=='H'){ // completely out of the interval, or insertion in the genome before the start position
					if(cigarcode!='I' && cigarcode!='P' && cigarcode!='H') refpos+=cigarcount; // advance position in reference
					if(cigarcode!='D' && cigarcode!='N' && cigarcode!='H' && cigarcode!='P') readpos+=cigarcount; // advance position in read
					continue; // next CIGAR operation
				}
				n=(startpos-refpos); // difference until beginning of valid interval
				cigarcount-=n; // skip invalid positions
				if(cigarcode!='I' && cigarcode!='P' && cigarcode!='H') refpos=startpos; // advance to beginning of interval
				if(cigarcode!='D' && cigarcode!='N' && cigarcode!='H' && cigarcode!='P') readpos+=n; // skip the read positions out of the interval
			}
			if((refpos+cigarcount)>endpos){ // operation ends after the end of the valid interval
				if(refpos>endpos || cigarcode=='H') break; // completely out of the interval or hard clipping (go to next read)
				n=(endpos-refpos+1); // difference until end of valid interval
				if(cigarcode!='I' && cigarcode!='P') cigarcount=n; // only count valid positions
			}
			j=(refpos-startpos); // first position/column inside the arrays
			if(cigarcode=='M' || cigarcode=='=' || cigarcode=='X'){ // match or mismatch
				i=j;
				while(cigarcount>0){ // copy chars from read to array
					chars[k][i]=read[readpos++]; // advance in read
					freeahead[k][i++]=0;
					refpos++; // advance in reference
					cigarcount--;
				}
				i=1;
				while(j>0 && freeahead[k][(j-1)]!=0){ // update number of free positions behind
					freeahead[k][(j-1)]=(i++);
					j--;
				}
			} else if(cigarcode=='I'){ // insertion into reference
				n=insids[j]; // insertion id in the current column
				if(n==0){ // still no insertions in this column
					n=insidscount;
					insids[j]=n; // set insertion id of this position
					insidscount++; // one more insertion id (the position in the arrays is always (insidscount-1) )
					inssize=(unsigned short *)realloc(inssize,insidscount*sizeof(unsigned short)); // alloc arrays for this new insertion
					inscount=(unsigned short *)realloc(inscount,insidscount*sizeof(unsigned short));
					insrows=(unsigned short **)realloc(insrows,insidscount*sizeof(unsigned short *));
					inschars=(char ***)realloc(inschars,insidscount*sizeof(char **));
					inssize[n]=cigarcount;
					inscount[n]=1;
					insrows[n]=(unsigned short *)malloc(inscount[n]*sizeof(unsigned short));
					inschars[n]=(char **)malloc(inscount[n]*sizeof(char *));
				} else { // this column already had insertions
					if(cigarcount>inssize[n]) inssize[n]=cigarcount; // update size of the largest insertion in this column
					inscount[n]++; // one more insertion in this column
					insrows[n]=(unsigned short *)realloc(insrows[n],inscount[n]*sizeof(unsigned short)); // realloc arrays to store one more column number and chars
					inschars[n]=(char **)realloc(inschars[n],inscount[n]*sizeof(char *));
				} // TODO: if 2 insertions are linked (or separated by padding), concatenate the strings
				insrows[n][(inscount[n]-1)]=(unsigned short)k; // store this row in the list of rows with insertions for this column
				inschars[n][(inscount[n]-1)]=(char *)malloc((cigarcount+1)*sizeof(char)); // alloc space for the insertion chars
				string=inschars[n][(inscount[n]-1)];
				i=0;
				while(cigarcount>0){ // copy chars from read to insertion string
					string[i++]=read[readpos++]; // advance in read
					cigarcount--;
				}
				string[i]='\0';
			} else if(cigarcode=='D'){ // deletion from reference
				i=j;
				while(cigarcount>0){ // put gap chars in array
					chars[k][i]='-';
					freeahead[k][i++]=0;
					refpos++;
					cigarcount--;
				}
				i=1;
				while(j>0 && freeahead[k][(j-1)]!=0){
					freeahead[k][(j-1)]=(i++);
					j--;
				}
			} else if(cigarcode=='S' || cigarcode=='N'){ // soft clipping or skipped
				refpos+=cigarcount; // advance in reference
				if(cigarcode=='S') readpos+=cigarcount; // advance in read if it is soft clipping
				i=j;
				while(cigarcount>0){ // put empty chars (dots) in array
					chars[k][i]='.';
					freeahead[k][i++]=0; // do not let other chars be put in this location
					cigarcount--;
				}
				i=1;
				while(j>0 && freeahead[k][(j-1)]!=0){
					freeahead[k][(j-1)]=(i++);
					j--;
				}
			} else if(cigarcode=='H'){ // hard clipping
				continue;
			} else if(cigarcode=='P'){ // padding
				if((n=insids[j])!=0 && insrows[n][(inscount[n]-1)]==((short int)k)) continue; // TODO: add padding to already inserted chars (concatenate)
			}
		} // end of processing CIGAR string
		numreads++; // count valid reads
		while(c!='\n' && c!=EOF) c=fgetc(file); // advance to next read
	} // end of loop for all reads inside SAM file
	fclose(file);
	if(numreads==0){
		printf("\n> ERROR: No reads found inside range\n");
		exit(0);
	}
	printf("OK (%u mapped reads - %u inside range)\n",totalnumreads,numreads);
	//if(numrows==USHRT_MAX) printf("> WARNING: More than %lu stacked reads detected but were ignored\n",USHRT_MAX);
	if(numrows>1024) printf("> WARNING: More than %d stacked reads were detected bellow at least one position but will not be shown in image\n",1024); // limit size of image
	printf("> Creating local reference coverage plot ... ");
	fflush(stdout);
	numpos--; // back to original number of positions
	endpos--;
	rowheight=( CHAR_HEIGHT + CHAR_SPACE );
	colwidth=( CHAR_WIDTH + CHAR_SPACE );
	margindim=colwidth;
	n=numpos;
	for(i=0;i<=numpos;i++){ // count total number of columns, including insertions between chars (+1 to account insertions at end)
		if((j=insids[i])!=0) n+=inssize[j];
	}
	k=DigitsCount(endpos+1); // make room to print half of a position number on both ends
	imgwidth=( 2*margindim + k*(CHAR_WIDTH+CHAR_SPACE) + n*colwidth ); // margin, half number, all ref chars, half number, margin
	if(numrows>1024) imgheight=( 2*margindim + (5+1024+4)*rowheight ); // limit size of image
	else imgheight=( 2*margindim + (5+numrows+4)*rowheight ); // top margins, numbers, markers, reference, consensus, dividing line, all rows, empty space and color labels
	if(gfffilename!=NULL) imgheight+=(rowheight+CHAR_SPACE); // gene names at the top
	initializeBitmap((int)imgwidth,(int)imgheight,2);
	if(gfffilename!=NULL){ // draw GFF annotations if present
		k=DigitsCount(endpos+1);
		startx=( margindim + (k*(CHAR_WIDTH+CHAR_SPACE))/2 ); // where to start drawing from
		starty=margindim;
		while(1){
			c=fgetc(gfffile);
			if(c==EOF) break;
			if(c=='#'){ // skip comment lines
				while(c!='\n' && c!=EOF) c=fgetc(gfffile);
				continue;
			}
			for(k=0;k<3;k++){ // skip: seqname, source, feature
				c=fgetc(gfffile);
				while(c!='\t' && c!='\n' && c!=EOF) c=fgetc(gfffile);
			}
			if(c!='\t') continue;
			nret = fscanf(gfffile,"%u\t%u",&i,&j); // get: start, end
			if(i!=0) i--;
			if(j!=0) j--;
			if(j<startpos || i>endpos){ // positions not inside the range of the image
				c=fgetc(gfffile);
				while(c!='\n' && c!=EOF) c=fgetc(gfffile);
				continue;
			}
			for(k=0;k<4;k++){ // skip: score, strand, frame
				c=fgetc(gfffile);
				while(c!='\t' && c!='\n' && c!=EOF) c=fgetc(gfffile);
			}
			if(c!='\t') continue;
			if(i<startpos) i=startpos; // clip positions to valid interval
			if(j>endpos) j=endpos;
			i-=startpos; // change positions in the genome to match positions in the arrays
			j-=startpos;
			x=startx; // horizontal pos in image
			m=0; // start position in the arrays
			for(k=0;k<=i;k++){ // get the position (m) in the arrays corresponding to the shifted position (i) in the genome 
				if(insids[k]!=0) m+=inssize[(insids[k])];
				m++;
			}
			m--; // it's one position ahead
			n=1; // number of available chars to write annotation
			for(k=(i+1);k<=j;k++){ // how many chars (n) in the arrays until the shifted position (j) in the genome
				if(insids[k]!=0) n+=inssize[(insids[k])];
				n++;
			}
			x=(startx+m*(CHAR_WIDTH+CHAR_SPACE));
			y=(starty-1); // draw one more pixel line above and bellow
			DrawFilledRectangle(x,y,(n*(CHAR_WIDTH+CHAR_SPACE)-CHAR_SPACE),(CHAR_HEIGHT+2),BLUE); // draw blue background
			x+=(CHAR_WIDTH+CHAR_SPACE); // leave one unfilled char on the right and one on the left
			y++;
			color=WHITE;
			if(n>2) n-=2;
			else n=0;
			c=fgetc(gfffile);
			while(n!=0 && c!='\t' && c!='\n' && c!=EOF){ // get: attributes
				if(c==';') c=' ';
				DrawChar(c,x,y,color);
				c=fgetc(gfffile);
				x+=(CHAR_WIDTH+CHAR_SPACE);
				n--;
			}
			while(c!='\n' && c!=EOF) c=fgetc(gfffile); // go to next line
		}
		fclose(gfffile);
	}
	if( imgwidth >= (2*margindim+24*colwidth) ){ // if there is enought horizontal space, draw color legend box
		startx=( (imgwidth/2) - (24*colwidth/2) ); // draw in the bottom center of the image
		starty=( imgheight - (margindim+3*rowheight) );
		x=startx;
		y=starty;
		DrawFilledRectangle(x,y,CHAR_HEIGHT,CHAR_HEIGHT,BLACK); // draw colored square
		x+=(CHAR_HEIGHT+CHAR_SPACE);
		DrawString("Reference",x,y,BLACK);
		x=startx;
		y+=rowheight;
		DrawFilledRectangle(x,y,CHAR_HEIGHT,CHAR_HEIGHT,GREY);
		x+=(CHAR_HEIGHT+CHAR_SPACE);
		DrawString("Consensus",x,y,BLACK);
		x=startx;
		y+=rowheight;
		DrawFilledRectangle(x,y,CHAR_HEIGHT,CHAR_HEIGHT,BLUE);
		x+=(CHAR_HEIGHT+CHAR_SPACE);
		DrawString("SNP/InDel",x,y,BLACK);
		startx+=( (CHAR_HEIGHT+CHAR_SPACE) + 10*colwidth );
		starty+=( rowheight/2 );
		x=startx;
		y=starty;
		DrawString("Conservation",x,y,BLACK);
		y+=rowheight;
		n=(12*(CHAR_WIDTH+CHAR_SPACE)-CHAR_SPACE);
		k=(256/n);
		for(i=0;i<n;i++){ // draw color gradient from red to green
			redvalue=(255-i*k);
			greenvalue=(255-redvalue);
			if(redvalue>greenvalue) addvalue=(255-redvalue);
			else addvalue=(255-greenvalue);
			redvalue+=addvalue;
			greenvalue+=addvalue;
			color=getColorFromPalette((uint8_t)redvalue,(uint8_t)greenvalue,0);
			DrawVerticalLine(x,y,(CHAR_HEIGHT+1),color);
			x++;
		}
		x-=n;
		DrawChar('+',(x+2),y,BLACK);
		x+=(11*colwidth);
		DrawChar('-',(x-2),y,BLACK);
	}
	k=DigitsCount(endpos+1);
	startx=( margindim + (k*(CHAR_WIDTH+CHAR_SPACE))/2 ); // where to start drawing from
	starty=margindim;
	if(gfffilename!=NULL) starty+=(rowheight+CHAR_SPACE); // space for annotations
	color=BLACK;
	x=( startx + CHAR_WIDTH/2 );
	y=starty;
	if((m=insids[0])!=0) x+=(inssize[m]*colwidth); // add extra space if insertions are present before the first position
	//DrawHorizontalLine((x-1),(starty+rowheight+CHAR_HEIGHT/2-1),((n-1)*colwidth+1),color); // draw horizontal line
	DrawNumber((startpos+1),x,y,color); // print left pos number
	DrawVerticalLine(x,(y+rowheight),(CHAR_HEIGHT/2+1),color);
	x=( imgwidth - startx - CHAR_SPACE - CHAR_WIDTH/2 );
	if((m=insids[numpos])!=0) x-=(inssize[m]*colwidth); // remove extra space if insertions are present after the last position
	DrawNumber((endpos+1),x,y,color); // print right pos number
	DrawVerticalLine(x,(y+rowheight),(CHAR_HEIGHT/2+1),color);
	n=(endpos-startpos+1); // number of available char positions
	n-=k; // take the space occupied by both halves of the numbers at the ends
	k=(k+1); // number of chars/columns occupied by each number, and a space following it
	m=1; // interval between mark numbers
	while( (n/m)*k > n ){ // find the smallest interval that can fit all the numbers in the available char positions
		m=(2*m); // double the spacing between numbers
		if( (n/m)*k <= n ) break;
		m=(5*(m/2)); // 5x spacing
		if( (n/m)*k <= n ) break;
		m=(2*m); // 10x spacing
	}
	if(n>=m){ // draw all position markers
		i=((startpos+1)+k); // position corresponding to the first column where a number can be printed
		n=((endpos+1)-k); // position corresponding to the last column where a number can be printed
		j=(((startpos+1)/m)*m); // leftmost position that is a multiple of the interval
		while(j<i) j+=m; // advance position to the first valid/printable position on the left
		i=0; // position inside the arrays
		k=(startpos+1);
		x=(startx+(j-k)*colwidth); // advance column to the corresponding position of j
		y=starty;
		while(j<=n){
			while(i<=(j-k)){ // if there are insertions between the previous mark and this one, advance in the columns
				if(insids[i]!=0) x+=((inssize[(insids[i])])*colwidth);
				i++;
			}
			DrawNumber(j,(x+CHAR_WIDTH/2+1),starty,color);
			DrawVerticalLine((x+CHAR_WIDTH/2),(starty+rowheight),(CHAR_HEIGHT/2+1),color);
			j+=m;
			x+=(m*colwidth); // advance to the next mark
		}
	}
	if(insids[0]!=0){ // if there are insertions before the first position, align all insertion chars to the right
		k=insids[0];
		n=inscount[k];
		for(i=0;i<n;i++){
			m=inssize[k];
			inschars[k][i]=(char *)realloc((inschars[k][i]),(m+1)*sizeof(char)); // realloc all strings to the same max size
			string=(inschars[k][i]);
			j=0;
			while(string[j]!='\0') j++; // find the end of each string
			if(j==m) continue;
			while(j!=0) string[m--]=string[j--]; // move chars to end of string
			string[m--]=string[j]; // move last/first char
			while(m!=0) string[m--]=' '; // fill remaining positions with spaces
			string[m]=' '; // fill last/first char
		}
	}
	starty+=(2*rowheight);
	x=startx;
	y=starty;
	for(i=0;i<=numpos;i++){ // draw reference chars at top (+1 for insertions at end)
		if((j=insids[i])!=0){
			n=(unsigned int)inssize[j];
			for(j=0;j<n;j++){ // draw dash chars where there are insertions
				DrawChar('-',x,y,color);
				x+=colwidth;
			}
		}
		if(i==numpos) break; // quit after insertions after last char
		DrawChar(refchars[i],x,y,color);
		x+=colwidth;
	}
	starty+=(2*rowheight);
	DrawHorizontalLine(startx,(starty+(rowheight/2)),(imgwidth-2*colwidth),color); // draw line between consensus and reads
	starty+=rowheight;
	numcharcountpos=1;
	charcount=(unsigned short int **)calloc(numcharcountpos,sizeof(unsigned short int *)); // char counters for one column (reallocable with insertion size)
	for(i=0;i<numcharcountpos;i++) charcount[i]=(unsigned short int *)calloc(5,sizeof(unsigned short int)); // counters for A,C,G,T,- chars
	for(j=0;j<=numpos;j++){ // process all columns/positions (+1 for insertions at end)
		m=(unsigned int)insids[j]; // id of the insertions at (before) this position
		if(m!=0){ // process insertions
			n=(unsigned int)(inssize[m]); // size of the maximum insertion in this position
			if(n>numcharcountpos){ // reallocate space for more column counters if needed
				charcount=(unsigned short int **)realloc(charcount,(n*sizeof(unsigned short int *)));
				for(i=numcharcountpos;i<n;i++) charcount[i]=(unsigned short int *)calloc(5,sizeof(unsigned short int));
				numcharcountpos=n;
			}
			for(i=0;i<n;i++){ // reset counters
				for(k=0;k<5;k++) charcount[i][k]=0;
			}
			count=0; // total count of rows that cross this column
			for(i=0;i<numrows;i++){ // count rows that have a char before and after the insertion point
				if( ( (j==0) || (j>0 && chars[i][(j-1)]!='\0') ) && chars[i][j]!='\0' ) count++;
			}
			n=(unsigned int)(inscount[m]); // number of rows with insertions
			if(count<n) count=n; // if all reads start with an insertion, they were not accounted above
			for(k=0;k<n;k++){ // update counts of all inserted chars in their respective columns
				string=inschars[m][k]; // get string of chars from insertion id
				i=0;
				while(string[i]!='\0'){
					c=string[i]; // update char count for column
					if(c=='A') charcount[i][0]++;
					else if(c=='C') charcount[i][1]++;
					else if(c=='G') charcount[i][2]++;
					else if(c=='T') charcount[i][3]++;
					else if(c=='-') charcount[i][4]++;
					i++;
				}
			}
			for(k=0;k<n;k++){ // draw all inserted chars in their respective rows
				i=(unsigned int)(insrows[m][k]); // get insertion row
				if(i>1024) continue; // limit size of image
				x=startx;
				y=(starty+i*rowheight);
				string=inschars[m][k]; // get string of chars from insertion id
				i=0;
				while(string[i]!='\0'){
					c=string[i];
					redvalue=0;
					if(c=='A') redvalue=(unsigned int)charcount[i][0]; // get red color intensity from char frequency in column
					else if(c=='C') redvalue=(unsigned int)charcount[i][1];
					else if(c=='G') redvalue=(unsigned int)charcount[i][2];
					else if(c=='T') redvalue=(unsigned int)charcount[i][3];
					else if(c=='-') redvalue=(unsigned int)charcount[i][4];
					redvalue=((255*redvalue)/count);
					greenvalue=(255-redvalue);
					if(redvalue>greenvalue) addvalue=(255-redvalue);
					else addvalue=(255-greenvalue);
					redvalue+=addvalue;
					greenvalue+=addvalue;
					color=getColorFromPalette((uint8_t)redvalue,(uint8_t)greenvalue,0); // draw char in tones from green to red (low to high frequency)
					DrawChar(c,x,y,color);
					x+=colwidth;
					i++;
				}
			}
			x=startx;
			y=(starty-2*rowheight); // write consensus sequence
			color=BLUE;
			n=(unsigned int)(inssize[m]); // size of the maximum insertion in this position
			for(i=0;i<n;i++){ // check if any of the inserted chars is frequent enought (>50%) to be added to the consensus
				for(k=0;k<5;k++){
					if(((unsigned int)charcount[i][k])>(count/2)){
						c='\0';
						if(k==0) c='A';
						else if(k==1) c='C';
						else if(k==2) c='G';
						else if(k==3) c='T';
						DrawChar(c,x,y,color);
						break;
					}
				}
				x+=colwidth;
			}
		} // end of dealing with insertions at/before this position
		if(j==numpos) break; // if we are one position ahead the last position, we just wanted to print the insertions at the right
		for(k=0;k<5;k++) charcount[0][k]=0; // reset char counts for this column
		count=0; // total count of chars in this column
		for(i=0;i<numrows;i++){ // count individual chars in this column
			c=chars[i][j];
			if(c=='A') charcount[0][0]++;
			else if(c=='C') charcount[0][1]++;
			else if(c=='G') charcount[0][2]++;
			else if(c=='T') charcount[0][3]++;
			else if(c=='-') charcount[0][4]++;
			if(c!='\0') count++;
		}
		n=0;
		c='\0';
		for(k=0;k<5;k++){ // get the most frequent char for this column (does *not* require to be >50%)
			if(((unsigned int)charcount[0][k])>n){
				if(k==0) c='A';
				else if(k==1) c='C';
				else if(k==2) c='G';
				else if(k==3) c='T';
				else if(k==4) c='-';
				n=((unsigned int)charcount[0][k]); // save highest count and corresponding char
			}
		}
		if(refchars[j]=='A' && charcount[0][0]==n) c='A'; // if the ref char has the same (most frequent) count as another char, prefer the ref char
		else if(refchars[j]=='C' && charcount[0][1]==n) c='C';
		else if(refchars[j]=='G' && charcount[0][2]==n) c='G';
		else if(refchars[j]=='T' && charcount[0][3]==n) c='T';
		if(c==refchars[j] || n==0) color=GREY; // if the char is the same as in the reference or if there are no chars bellow, draw it in grey
		else color=BLUE; // if it is a different char, draw it in blue
		x=(startx+(inssize[m]*colwidth)); // draw char in column next to the inserted chars
		y=starty;
		if(c!='\0' && count>0) DrawChar(c,x,(y-2*rowheight),color); // draw consensus char (only if it has reads bellow it)
		for(i=0;i<numrows;i++){ // draw vertically all chars of this column
			if(i>1024) break; // limit size of image
			c=chars[i][j];
			if(c!='\0'){
				redvalue=0;
				if(c=='A') redvalue=(unsigned int)charcount[0][0]; // get red color intensity from char frequency in column
				else if(c=='C') redvalue=(unsigned int)charcount[0][1];
				else if(c=='G') redvalue=(unsigned int)charcount[0][2];
				else if(c=='T') redvalue=(unsigned int)charcount[0][3];
				else if(c=='-') redvalue=(unsigned int)charcount[0][4];
				redvalue=((255*redvalue)/count);
				greenvalue=(255-redvalue);
				if(redvalue>greenvalue) addvalue=(255-redvalue);
				else addvalue=(255-greenvalue);
				redvalue+=addvalue;
				greenvalue+=addvalue;
				color=getColorFromPalette((uint8_t)redvalue,(uint8_t)greenvalue,0); // draw char in tones from green to red (low to high frequency)
				DrawChar(c,x,y,color);
			}
			y+=rowheight;
		}
		startx+=((1+inssize[m])*colwidth); // next drawn char will be after all inserted chars in this position
	}
	free(refchars);
	for(i=0;i<numrows;i++){
		free(chars[i]);
		free(freeahead[i]);
	}
	free(chars);
	free(freeahead);
	free(insids);
	free(inssize);
	for(i=1;i<insidscount;i++){
		n=inscount[i];
		for(j=0;j<n;j++) free(inschars[i][j]);
		free(insrows[i]);
		free(inschars[i]);
	}
	free(inscount);
	free(insrows);
	free(inschars);
	free(read);
	for(i=0;i<numcharcountpos;i++) free(charcount[i]);
	free(charcount);
	printf("OK\n");
}


void DrawMappingPlot(char *reffilename, char *readsfilename, char *gfffilename, unsigned int startpos, unsigned int endpos){
	FILE *file;
	unsigned int refsize, numpos;
	unsigned int i, k;
	char c, *reflabel;
	char *imagefilename, filenamesuffix[28];
	printf("> Processing reference genome file <%s> ... ",reffilename);
	fflush(stdout);
	if((file=fopen(reffilename,"r"))==NULL){
		printf("\n> ERROR: File not found\n");
		exit(0);
	}
	if((c=fgetc(file))!='>'){
		printf("\n> ERROR: Not a file in the FASTA format\n");
		exit(0);
	}
	k=0;
	while((c=fgetc(file))!='\n' && c!=EOF) k++; // get label size
	reflabel=(char *)calloc((k+1),sizeof(char));
	rewind(file);
	k=0;
	c=fgetc(file);
	while((c=fgetc(file))!='\n' && c!=EOF) reflabel[k++]=c; // get label
	reflabel[k]='\0';
	refsize=0;
	i=1;
	k=0; // number of non ACGT chars, but still letters
	while((c=fgetc(file))!=EOF){ // get reference size
		if(c=='>'){ // skip description in case of multi-fasta file
			while(c!='\n' && c!=EOF) c=fgetc(file);
			i++;
		}
		if(c=='A' || c=='C' || c=='G' || c=='T' || c=='a' || c=='c' || c=='g' || c=='t') refsize++; // count valid chars
		else if( (c>=65 && c<=90) || (c>=97 && c<=122) ) k++; // other letters are invalid and all other chars are not counted
	}
	if(refsize==0){
		printf("\n> ERROR: No DNA characters in file\n");
		exit(0);
	}
	if(i>1){
		printf("\n> ERROR: Found %d sequences but visualization is only implemented for one single reference sequence\n",i);
		exit(0);
	}
	if(k==0) printf("OK (%u basepairs)\n",refsize);
	else {
		printf("OK (%u basepairs ; %d invalid chars)\n",refsize,k);
		//printf("> WARNING: Non ACGT characters were found and may lead to wrong alignment visualization results\n");
		refsize+=k;
	}
	fclose(file);
	if( startpos>0 && endpos>0 ){ // convert from SAM 1-based position to arrays 0-based position
		startpos--;
		endpos--;
	}
	if( startpos>=refsize ) startpos=0;
	if( endpos<=startpos || endpos>refsize ) endpos=(refsize-1);
	numpos=(endpos-startpos+1); // size of the array of counts
	if (numpos > 10000) { // decide whether to create global or local map depending on image size for all positions
		CreateGlobalMap(startpos, endpos, numpos, readsfilename); // if there are more than 10000 positions to draw, create global map
	} else {
		CreateLocalMap(startpos,endpos,numpos,refsize,reffilename,readsfilename,gfffilename);
	}
	if (numpos == refsize) strcpy(filenamesuffix, ".bmp");
	else sprintf(filenamesuffix, "(%u-%u).bmp", (startpos + 1), (endpos + 1));
	/*
	k=0; // replace reads filename extention with ".bmp"
	while(readsfilename[k]!='\0') k++; // get reads filename string size
	i=k;
	while(i>0 && readsfilename[i]!='.') i--; // get position of rightmost dot
	if(i==0) i=k; // set whole filename if no dot exists
	imagefilename=(char *)calloc((i+4+1),sizeof(char));
	for(k=0;k<i;k++) imagefilename[k]=readsfilename[k]; // copy chars until dot
	imagefilename[k++]='.'; // append extention ".bmp"
	imagefilename[k++]='b';
	imagefilename[k++]='m';
	imagefilename[k++]='p';
	imagefilename[k++]='\0';
	*/
	imagefilename=AppendToBasename(readsfilename,filenamesuffix);
	printf("> Saving image to <%s> ... ",imagefilename);
	fflush(stdout);
	if(saveBitmap(imagefilename)==0){
		printf("\n> ERROR: Cannot write file\n");
		exit(0);
	}
	printf("OK\n");
	free(reflabel);
	free(imagefilename);
	printf("> Done!\n");
	exit(0);
}
