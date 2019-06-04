#include <stdio.h>
#include <stdlib.h>
#include "alignreads.h"
#include "dynprog.h"
#include "samoutput.h"

#ifdef _MSC_VER
	#pragma warning(disable:4996)
#endif

#define MATCH 0
#define MISMATCH 1
#define GAP 1

#if defined DEBUG || defined DEBUGDP
void PrintBandedDPMatrix(char *pattern, int m, char *text, int n, int maxerrors){
	int i, j, ii, mm;
	char c;
	mm = ( 2 * maxerrors ); // diameter of diagonal (central cell plus maxerrors up and maxerrors down)
	fprintf(debugfile,"\n");
	fprintf(debugfile,"       |"); // empty space for row number and char below
	for(j=0;j<=n;j++) fprintf(debugfile,"%3d |",j); // column number
	fprintf(debugfile,"\n");
	fprintf(debugfile,"     %c |",('*')); // empty space for row number and char below
	fprintf(debugfile,"  - |"); // 0-th char column
	for(j=1;j<=n;j++) fprintf(debugfile,"  %c |",(text[(n-j)])); // text chars: text[(j-1)]
	fprintf(debugfile,"\n");
	for(i=0;i<=m;i++){ // process all rows
		fprintf(debugfile,"%3d",i); // row number
		if(i==0) fprintf(debugfile,"| - |"); // 0-th char row
		else fprintf(debugfile,"| %c |",(pattern[(m-i)])); // pattern chars: pattern[(i-1)]
		ii=(maxerrors+i); // row in real matrix
		for(j=0;j<=n;j++){ // process all columns
			if(ii<0 || ii>mm) fprintf(debugfile,"----|"); // cells outside banded diagonal
			else { // cells inside banded diagonal
				switch(dpDirections[ii][j]){
				case 'd':
					c='`';	break;
				case 'D':
					c='\\';	break;
				case 'U':
					c='^'; break;
				case 'L':
					c='<'; break;
				default:
					c='?';
				}
				fprintf(debugfile,"%3d%c|",dpMatrix[ii][j],c);
			}
			ii--; // the cell to the left is one position above in the real matrix
		}
		fprintf(debugfile,"\n");
	}
	fprintf(debugfile,"\n");
}
#endif

// Performs banded dynamic programming on the read segments between the seeds
// NOTE: seedlocation is 0 on the leftmost end of the read and -1 on the rightmost end of the read
// NOTE: all strings are processed in reverse, so all homopolymer indel errors fall on the left side (beginning) of the homopolymer
// NOTE: the CIGAR strings are filled in the correct direction by backtracking the DP matrix directions
// NOTE: the returned number of errors (edit distance) is equal to (matrix[m][n]) only if the default scores are kept: MATCH=0, MISMATCH=1, GAP=1
int RunBandedDynamicProgramming(char *text, int n, char *pattern, int m, int maxerrors, int seedlocation, int *numgapsatend){
	int i, ii, j, mm, k;
	int matchscore, up, left, diag;
	char dir;
	#ifdef DEBUGDP
	static int h, inscount, delcount, miscount;
	#endif

	if(maxerrors==0){ // if no errors are allowed
		if( cigarStrCodes[cigarStrSize] == 'M' ){ // if the previous operation was also a matching
			cigarStrCounts[cigarStrSize] += m;
		} else {
			cigarStrSize++;
			cigarStrCodes[cigarStrSize] = 'M';
			cigarStrCounts[cigarStrSize] = m;
		}
		#ifdef DEBUGDP
		fprintf( debugfile , "%.*s\n" , m , text );
		fprintf( debugfile , "%.*s\n" , m , pattern );
		for( i = 0 ; i <= m ; i++ ) fprintf( debugfile, "M" );
		fprintf( debugfile , "\n" );
		#endif
		(*numgapsatend)=0; // both sizes m and n should have the same value, else something is wrong
		k=0; // number of errors
		for(i=0;i<m;i++) if(pattern[i]!=text[i]) k++;
		return k;
	}
	
	k = n; // maximum between text size and pattern size
	if( m > n ) k = m;
	if( maxerrors > k ) maxerrors = k; // the number of errors can be limited by the maximum of both lengths

	// need to initialize the first row and column everytime, when using a different number of errors each time
	ii = maxerrors; // 0-th virtual row is ME-th real row
	k = maxerrors;
	if( m < maxerrors ) k = m; // if the vertical string is shorter than ME, no need to initialize more cells than those
	for( i = 0 ; i <= k ; i++ ){ // initialize partial first column (part that falls inside banded diagonal)
		dpMatrix[ii][0] = (i*GAP);
		dpDirections[ii][0] = 'U';
		ii++; // from real row ii=ME to row ii=2*ME (virtual row i=0 to i=ME)
	}
	ii = maxerrors; // 0-th virtual row is ME-th real row
	k = maxerrors;
	if( n < maxerrors ) k = n; // if the horizontal string is shorter than ME, no need to initialize more cells than those
	if( seedlocation != -1 ){ // not rightmost segment
		for( j = 0 ; j <= k ; j++ ){ // initialize partial first row (part that falls inside banded diagonal)
			dpMatrix[ii][j] = (j*GAP);
			dpDirections[ii][j] = 'L';
			ii--; // from real row ii=ME to row ii=0 (virtual row i=0)
		}
	} else { // in the rightmost segment, the top row is initialized to all zeros
		for( j = 0 ; j <= k ; j++ ){
			dpMatrix[ii][j] = 0;
			dpDirections[ii][j] = 'L';
			ii--;
		}
	}
	dpDirections[maxerrors][0] = 'D';

	mm = ( 2 * maxerrors ); // total number of available rows (size of each column) to fill (around the diagonal)
	// the pattern/text are matched in reverse, i.e. from end/right to beginning/left, so indels appear at the beginning of matches
	for( j = 1 ; j <= n ; j++ ){ // iterate through all columns, one at a time ; the (real) rows are processed in 3 parts: first row, all middle rows, and last row
		i = ( j - maxerrors ); // i = (j-NE) : (j+NE) ; from NE positions above the diagonal cell [i,j] of this column to NE positions bellow it
		ii = 0; // ii = 0 : (2*NE)
		if( i >= 1 ){ // first cell of (real) column, ii=0, can only be from L or D, but only fill if (i > 0), because if (i = 0) it is the top (virtual) row and is already a L
			if( pattern[(m-i)] == text[(n-j)] ) matchscore = MATCH; // the string chars are fetched in reverse order
			else matchscore = MISMATCH;
			diag = (dpMatrix[ii][j-1] + matchscore);
			left = (dpMatrix[ii+1][j-1] + GAP);
			if( (diag <= left) ) {
				dpMatrix[ii][j] = diag;
				dpDirections[ii][j] = ( (matchscore==MATCH) ? 'D' : 'd' );
			} else {
				dpMatrix[ii][j] = left;
				dpDirections[ii][j] = 'L';
			}
			i++; // next row (cell bellow)
			ii++; // ii = 1
		} else { // if we are before the first virtual row (i<=0), advance to it
			ii = (1-i);
			i = 1;
		}
		for( ; ii < mm ; ii++ ){ // middle cells of column, but only fill if (i >= 1) and (i <= m)
			if( pattern[(m-i)] == text[(n-j)] ) matchscore = MATCH;
			else matchscore = MISMATCH;
			diag = (dpMatrix[ii][j-1] + matchscore); // diagonal: (i-1,j-1) = (ii,j-1)
			left = (dpMatrix[ii+1][j-1] + GAP); // left: (i,j-1) = (ii+1,j-1)
			up = (dpMatrix[ii-1][j] + GAP); // up: (i-1,j) = (ii-1,j)
			if( (up <= left) && (up <= diag) ) { // on error prefer indels over mismatches
				dpMatrix[ii][j] = up;
				dpDirections[ii][j] = 'U';
			} else if( left <= diag ) { // on indel prefer consuming chars of the read (up) than consuming chars of the genome (left)
				dpMatrix[ii][j] = left;
				dpDirections[ii][j] = 'L';
			} else {
				dpMatrix[ii][j] = diag;
				dpDirections[ii][j] = ( (matchscore==MATCH) ? 'D' : 'd' ); // uppercase or lowercase in case of match or mismatch
			}
			i++; // next row (cell bellow)
			if( i > m ) break; // if this was the last valid virtual row of the column, do no process any more real but unnecessary rows
		} // at the end of the cycle, both i and ii variable are pointing to the last row
		if( i <= m ){ // last cell of column, ii=(2*NE), can only be from U or D, but only fill if (i <= m)
			if( pattern[(m-i)] == text[(n-j)] ) matchscore = MATCH;
			else matchscore = MISMATCH;
			diag = (dpMatrix[ii][j-1] + matchscore);
			up = (dpMatrix[ii-1][j] + GAP);
			if( diag <= up ) {
				dpMatrix[ii][j] = diag;
				dpDirections[ii][j] = ( (matchscore==MATCH) ? 'D' : 'd' );
			} else {
				dpMatrix[ii][j] = up;
				dpDirections[ii][j] = 'U';
			}
		}
	} // end of columns loop
	/*
	#ifdef DEBUGDP
	PrintBandedDPMatrix(pattern,m,text,n,maxerrors);
	#endif
	*/
	k=n; // in which column to start the backtracking process
	if(seedlocation==0){ // if it is the leftmost segment, the un-needed part of the reference chars to the left (right in the matrix) will be cut off
		j=n; // last column
		ii=(maxerrors+m-n); // real row of the bottom-right matrix cell
		mm=(2*maxerrors); // last real row in each column
		i=dpMatrix[ii][j]; // stores the best score found so far
		while(ii<=mm && j>0){ // process entire (virtual) bottom row
			if(dpMatrix[ii][j]<=i){ // and search for the right-most column (but still inside the band) with the best score
				i=dpMatrix[ii][j]; // track current best score
				k=j; // track current best column
			}
			j--; // previous column
			ii++;
		}
	}

	i = m; // start at bottom-right end of matrix and backtrace directions
	j = k; // starting column
	ii = ( maxerrors + m - k ); // if we had (|n-m|>NE) we needed another check here: (ii=>0 && ii<=2*NE)
	while( i > 0 && j > 0 ){
		dir = dpDirections[ii][j];
		switch( dir ){
		case 'd': // diagonal with mismatch, and continue with case 'D'
		case 'D': // diagonal with match: i-- , j-- , ii=ii
			dir = 'M';
			i--;
			j--;
			break;
		case 'U': // up: i-- , j=j , ii--
			dir = 'I';
			i--;
			ii--;
			break;
		case 'L': // left: i=i , j-- , ii++
			dir = 'D';
			j--;
			ii++;
			break;
		default:
			break;
		}
		// fill the CIGAR strings
		if( cigarStrCodes[cigarStrSize] == dir ){ // if it's the same operation as the previous one
			cigarStrCounts[cigarStrSize]++;
		} else { // if not, create new operation
			cigarStrSize++; // holds the last position of the array, size is +1
			cigarStrCodes[cigarStrSize] = dir;
			cigarStrCounts[cigarStrSize] = 1;
		}
	}
	if( i > 0 ){ // 'U', remaining gaps in the beginning of the pattern
		cigarStrSize++;
		cigarStrCodes[cigarStrSize] = 'I';
		cigarStrCounts[cigarStrSize] = i;
	}
	if( seedlocation != -1 ){ // on the rightmost segment, the first row is filled with zeros and we do not need to consume text chars
		if( j > 0 ){ // 'L', remaining gaps in the beginning of the text
			cigarStrSize++;
			cigarStrCodes[cigarStrSize] = 'D';
			cigarStrCounts[cigarStrSize] = j;
		}
	}
	#ifdef DEBUGDP
	inscount = 0;
	delcount = 0;
	miscount = 0;
	i = m;
	j = k;
	ii = ( maxerrors + m - k );
	h = 0;
	dpTraceSize = 0; // current size of dpTrace array
	while( i > 0 && j > 0 ){
		dir = dpDirections[ii][j];
		switch( dir ){
		case 'd':
			miscount++; // update mismatches count and continue with case 'D'
		case 'D':
			dpTrace[dpTraceSize] = 'M';
			dpAlignedTarget[h] = text[(n-j)]; // normal direction would be: text[(j-1)]
			dpAlignedRead[h] = pattern[(m-i)]; // normal direction would be: pattern[(i-1)];
			i--;
			j--;
			break;
		case 'U':
			dpTrace[dpTraceSize] = 'I';
			dpAlignedTarget[h] = '-';
			dpAlignedRead[h] = pattern[(m-i)];
			inscount++;
			i--;
			ii--;
			break;
		case 'L':
			dpTrace[dpTraceSize] = 'D';
			dpAlignedTarget[h] = text[(n-j)];
			dpAlignedRead[h] = '-';
			delcount++;
			j--;
			ii++;
			break;
		default:
			break;
		}
		dpTraceSize++;
		h++;
	}
	while( i > 0 ){ // 'U', follow remaining gaps in the beginning of the pattern
		dpTrace[dpTraceSize] = 'I';
		dpTraceSize++;
		dpAlignedTarget[h] = '-';
		dpAlignedRead[h] = pattern[(m-i)];
		h++;
		inscount++;
		i--;
	}
	if( seedlocation != -1 ){
		while( j > 0 ){ // 'L', follow remaining gaps in the beginning of the text
			dpTrace[dpTraceSize] = 'D';
			dpTraceSize++;
			dpAlignedTarget[h] = text[(n-j)];
			dpAlignedRead[h] = '-';
			h++;
			delcount++;
			j--;
		}
	}
	h--;
	for( i = 0 ; i <= h ; i++ ) fprintf( debugfile , "%c" , dpAlignedTarget[i] );
	fprintf( debugfile , "\n" );
	for( i = 0 ; i <= h ; i++ ) fprintf( debugfile , "%c" , dpAlignedRead[i] );
	fprintf( debugfile , "\n" );
	for( i = 0 ; i <= h ; i++ ) fprintf( debugfile , "%c" , dpTrace[i] );
	fprintf( debugfile , "\n" );
	//numerrors = (miscount + inscount + delcount);
	#endif

	(*numgapsatend) = (n-k); // number of gaps at the left end of the pattern (always zero except on leftmost segment)
	ii = ( maxerrors + m - k );
	return (dpMatrix[ii][k]);
}

// performs basic semi-global full dynamic programming to find the number of errors in the best position where the pattern matches the text
int BasicDynamicProgramming(char *text, int txtsize, char *pattern, int patsize){
	int i, j;
	int matchscore, up, left, diag, bestscore, bestcol;
	static int **dpmatrix = NULL;
	static int maxrows = 0;
	static int maxcols = 0;
	bestcol = 0; // to prevent un-initialized variable warning
	bestcol = (int)bestcol;
	/*
	static char **dpdirections, *alignedtxt, *alignedpat, dir;
	static int inscount, delcount, miscount, alignsize;
	printf("(%d)%.*s\n(%d)%.*s\n",txtsize,txtsize,text,patsize,patsize,pattern);
	*/
	if( dpmatrix == NULL ){ // simple alloc the first time
		maxrows = patsize;
		maxcols = txtsize;
		dpmatrix = (int **)malloc((maxrows+1)*sizeof(int *));
		for( i = 0 ; i <= maxrows ; i++ ){
			dpmatrix[i] = (int *)malloc((maxcols+1)*sizeof(int));
			dpmatrix[i][0] = i; // 0-th column: score=(i*gapscore)
		}
		for( j = 0 ; j <= maxcols ; j++ ){
			dpmatrix[0][j] = 0; // 0-th row: score=0 (to allow the pattern to start at any position in the text)
		}
		/*
		dpdirections = (char **)malloc((maxrows+1)*sizeof(char *));
		for( i = 0 ; i <= maxrows ; i++ ){
			dpdirections[i] = (char *)malloc((maxcols+1)*sizeof(char));
			dpdirections[i][0] = 'U';
		}
		for( j = 0 ; j <= maxcols ; j++ ){
			dpdirections[0][j] = 'L';
		}
		dpdirections[0][0] = 'D';
		alignsize = ( maxrows + maxcols );
		alignedtxt = (char *)malloc((alignsize+1)*sizeof(char));
		alignedpat = (char *)malloc((alignsize+1)*sizeof(char));
		*/
	}
	if( (txtsize > maxcols) || (patsize > maxrows) ){ // realloc matrix space if needed
		if( patsize > maxrows ){
			dpmatrix = (int **)realloc(dpmatrix,(patsize+1)*sizeof(int *)); // realloc matrix to have more rows
			for( i = (maxrows+1) ; i <= patsize ; i++ ){ // alloc new rows only
				dpmatrix[i] = (int *)malloc((maxcols+1)*sizeof(int));
				dpmatrix[i][0] = i; // 0-th column
			}
			maxrows = patsize;
		}
		if( txtsize > maxcols ){
			maxcols = txtsize;
			for( i = 0 ; i <= maxrows ; i++ ){ // realloc all rows to larger size
				dpmatrix[i] = (int *)realloc(dpmatrix[i],(maxcols+1)*sizeof(int));
				dpmatrix[i][0] = i; // 0-th column
			}
			for( j = 0 ; j <= maxcols ; j++ ){ // 0-th row
				dpmatrix[0][j] = 0;
			}
		}
		/*
		dpdirections = (char **)realloc(dpdirections,(maxrows+1)*sizeof(char *));
		for( i = 0 ; i <= maxrows ; i++ ){
			dpdirections[i] = (char *)realloc(dpdirections[i],(maxcols+1)*sizeof(char));
			dpdirections[i][0] = 'U';
		}
		for( j = 0 ; j <= maxcols ; j++ ){
			dpdirections[0][j] = 'L';
		}
		dpdirections[0][0] = 'D';
		alignsize = ( maxrows + maxcols );
		alignedtxt = (char *)realloc(alignedtxt,(alignsize+1)*sizeof(char));
		alignedpat = (char *)realloc(alignedpat,(alignsize+1)*sizeof(char));
		*/
	}
	for( i = 1 ; i <= patsize ; i++ ){ // loop for all rows
		for( j = 1 ; j <= txtsize ; j++ ){ // loop for all columns
			if( pattern[(i-1)] == text[(j-1)] ) matchscore = 0; // match
			else matchscore = 1; // mismatch
			diag = (dpmatrix[i-1][j-1] + matchscore);
			left = (dpmatrix[i][j-1] + 1); // gap (deletion from text)
			up = (dpmatrix[i-1][j] + 1); // gap (insertion into text)
			if( (diag <= left) && (diag <= up) ) dpmatrix[i][j] = diag;
			else if( up <= left ) dpmatrix[i][j] = up;
			else dpmatrix[i][j] = left;
			/*
			if( (diag <= left) && (diag <= up) ) dpdirections[i][j] = 'D';
			else if( up <= left ) dpdirections[i][j] = 'U';
			else dpdirections[i][j] = 'L';
			*/
		} // end of cols loop
	} // end of rows loop
	i = patsize; // bottom row
	j = txtsize; // right-most column
	bestscore = dpmatrix[i][j];
	bestcol = txtsize;
	while( j > 0 ){ // find best (lowest) score on bottom row
		if( dpmatrix[i][j] <= bestscore ){
			bestscore = dpmatrix[i][j];
			bestcol = j;
		}
		j--;
	}
	/*
	inscount = 0;
	delcount = 0;
	miscount = 0;
	i = patsize;
	j = bestcol;
	alignsize = 0;
	while( i > 0 && j > 0 ){
		dir = dpdirections[i][j];
		switch( dir ){
		case 'D':
			alignedtxt[alignsize] = text[(j-1)];
			alignedpat[alignsize] = pattern[(i-1)];
			if( pattern[(i-1)] != text[(j-1)] ){
				miscount++;
			}
			i--;
			j--;
			break;
		case 'U':
			alignedtxt[alignsize] = '-';
			alignedpat[alignsize] = pattern[(i-1)];
			inscount++;
			i--;
			break;
		case 'L':
			alignedtxt[alignsize] = text[(j-1)];
			alignedpat[alignsize] = '-';
			delcount++;
			j--;
			break;
		default:
			break;
		}
		alignsize++;
	}
	while( i > 0 ){
		alignedtxt[alignsize] = '-';
		alignedpat[alignsize] = pattern[(i-1)];
		alignsize++;
		inscount++;
		i--;
	}
	for( i = (alignsize-1) ; i >= 0 ; i-- ) putchar( alignedpat[i] );
	putchar( '\n' );
	for( i = (alignsize-1) ; i >= 0 ; i-- ) putchar( alignedtxt[i] );
	putchar( '\n' );
	printf( "%d %d\n" , ( miscount + inscount + delcount ) , bestscore );
	*/
	return bestscore; // edit distance (number of errors)
}

void GetRealErrorStatistics(char *reffilename, char *readsfilename, char *samfilename){
	FILE *samfile, *reffile, *readsfile;
	char c, *refchars, *readchars, *readname;
	int i, n, m, refsize, readpos, readsize, sameread, maxreadsize;
	int realnumreads, numreads, numalignedreads, numerrors, prevnumerrors, strand;
	unsigned totalnumchars, totalnumerrors;
	int nret = 0; // number of items returned by fscanf
	nret = (int)nret; // prevent unused variable warning
	printf("> Processing reference file <%s> ... ",reffilename);
	fflush(stdout);
	if((reffile=fopen(reffilename,"r"))==NULL){
		printf("\n> ERROR: File not found\n");
		exit(0);
	}
	c=fgetc(reffile);
	if(c!='>'){
		printf("\n> ERROR: Not a valid FASTA file\n");
		exit(0);
	}
	ungetc(c,reffile);
	n=0; // number of reference sequences in file
	refsize=0;
	while((c=fgetc(reffile))!=EOF){ // get reference chars
		if(c=='>'){
			while(c!='\n' && c!=EOF) c=fgetc(reffile); // skip reference label
			n++;
			continue;
		}
		if( (c>=65 && c<=90) || (c>=97 && c<=122) ) refsize++; // count only alphanumeric chars
	}
	if(refsize==0 || n!=1){
		printf("\n> ERROR: No valid chars or not single reference\n");
		exit(0);
	}
	refchars=(char *)malloc((refsize+1)*sizeof(char));
	rewind(reffile);
	n=0; // number of invalid chars
	i=0; // position of current char
	while((c=fgetc(reffile))!=EOF){ // get reference chars
		if(c=='>'){
			while(c!='\n' && c!=EOF) c=fgetc(reffile); // skip description
			continue;
		}
		if( c>=97 && c<=122 ) c=(char)(c-32); // convert to uppercase if needed
		if( c>=65 && c<=90 ){
			refchars[i++]=c;
			if( c!='A' && c!='C' && c!='G' && c!='T' ) n++; // alphabet letter but not ACGT
		}
	}
	refchars[i]='\0';
	fclose(reffile);
	printf("(%d chars) ",refsize);
	if(n!=0) printf("(%d invalid chars) ",n);
	printf("OK\n");
	printf("> Processing reads file <%s> ... ",readsfilename);
	fflush(stdout);
	if((readsfile=fopen(readsfilename,"r"))==NULL){
		printf("\n> ERROR: File not found\n");
		exit(0);
	}
	c=fgetc(readsfile);
	if(c!='>'){
		printf("\n> ERROR: Not a valid FASTA file\n");
		exit(0);
	}
	ungetc(c,readsfile);
	realnumreads=0; // number of reads inside file
	while((c=fgetc(readsfile))!=EOF) if(c=='>') realnumreads++;
	rewind(readsfile);
	printf("(%d reads) OK\n",realnumreads);
	printf("> Processing mapped reads file <%s> ... ",samfilename);
	fflush(stdout);
	if((samfile=fopen(samfilename,"r"))==NULL){
		printf("\n> ERROR: File not found\n");
		exit(0);
	}
	maxreadsize=1;
	readchars=(char *)malloc((maxreadsize+1)*sizeof(char));
	readname=(char *)malloc((255+1)*sizeof(char));
	sameread=0;
	readsize=0;
	strand=0;
	prevnumerrors=0;
	totalnumchars=0;
	totalnumerrors=0;
	numreads=0;
	numalignedreads=0;
	c=fgetc(samfile);
	while(c=='@'){ // header
		while(c!='\n' && c!=EOF) c=fgetc(samfile);
		c=fgetc(samfile);
	}
	ungetc(c,samfile);
	while((c=fgetc(samfile))!=EOF){
		sameread=1;
		i=0;
		while(c!='\t' && c!=EOF){ // get QNAME
			if(i<255){
				sameread &= (c==readname[i]);
				readname[i++]=c;
			}
			c=fgetc(samfile);
		}
		readname[i]='\0';
		if(!sameread){ // if it is a new read
			numreads++; // one more read
			n=i; // size of the read name
			while(1){ // get the read chars for this read name in the reads file
				c=fgetc(readsfile);
				if(c!='>' && c!=EOF){
					printf("\n> ERROR: Invalid read #%d in reads file\n",numreads);
					exit(0);
				}
				i=0;
				c=fgetc(readsfile); // compare read name in the reads file with the read name from the sam file
				while(i<n && c!='\n' && c!=EOF && c==readname[i]){
					c=fgetc(readsfile);
					i++;
				}
				if(c==EOF){ // the read name was not found
					printf("\n> ERROR: Read name '%s' not found in reads file\n",readname);
					exit(0);
				}
				if(i!=n){ // not the same read, check the next one
					while(c!='>' && c!=EOF) c=fgetc(readsfile); // skip this read
					ungetc(c,readsfile);
					continue;
				}
				while(c!='\n' && c!=EOF) c=fgetc(readsfile); // skip rest of read name, if it exists
				readsize=0;
				c=fgetc(readsfile);
				while(c!='\n' && c!=EOF){ // get read chars
					if(c>=97 && c<=122) readchars[readsize]=(char)(c-32);
					else readchars[readsize]=c;
					readsize++;
					if(readsize==maxreadsize){ // increase read array size if needed
						maxreadsize+=readsize;
						readchars=(char *)realloc(readchars,(maxreadsize+1)*sizeof(char));
					}
					c=fgetc(readsfile);
				}
				readchars[readsize]='\0';
				strand=0; // forward strand
				break;
			}
		}
		/*
		c=fgetc(samfile);
		while(c!='\t' && c!=EOF) c=fgetc(samfile); // skip FLAG
		*/
		nret = fscanf(samfile,"%d\t",&n); // get FLAG
		if((n & 16) != strand){ // different strand from the read chars
			strand=(n & 16); // set new strand
			for(i=0;i<readsize;i++){ // complement all the read chars
				c=readchars[i];
				if(c=='A') readchars[i]='T';
				else if(c=='C') readchars[i]='G';
				else if(c=='G') readchars[i]='C';
				else if(c=='T') readchars[i]='A';
			}
			n=(readsize/2); // middle of read
			for(i=0;i<=(n-1);i++){ // reverse the read
				m=(readsize-1-i); // position in the oposite side of the read
				readchars[i] ^= readchars[m]; // bitwise XOR swap (left position with right position)
				readchars[m] ^= readchars[i];
				readchars[i] ^= readchars[m];
			}
		}
		c=fgetc(samfile);
		if(c=='*' || c==EOF){
			while(c!='\n' && c!=EOF) c=fgetc(samfile);
			continue;
		}
		while(c!='\t' && c!=EOF) c=fgetc(samfile); // skip RNAME
		readpos=0;
		nret = fscanf(samfile,"%d\t",&readpos); // get POS
		if(readpos==0 || readpos>=refsize){
			while(c!='\n' && c!=EOF) c=fgetc(samfile);
			continue;
		}
		readpos--; // fix 1-based location
		c=fgetc(samfile);
		while(c!='\t' && c!=EOF) c=fgetc(samfile); // skip MAPQ
		m=0; // size of the region in the reference
		i=0; // number of cigar operations
		while(fscanf(samfile,"%d%c",&n,&c)>0){ // get CIGAR
			if(c=='I') m-=n; // insertion: ref chars are not consumed
			else if(c=='D') m+=n; // deletion: more ref chars ahead are used
			else if(c=='S' || c=='H'){ // soft or hard clipping
				if(i==0) readpos-=n; // at the beginning of the read: more ref chars are needed behind
				else m+=n; // at the end of the read: more ref chars are needed ahead
				i--; // do not count these operations
			}
			i++;
		}
		for(i=0;i<4;i++){ // skip RNEXT, PNEXT, TLEN , SEQ
			c=fgetc(samfile);
			while(c!='\t' && c!=EOF) c=fgetc(samfile);
		}
		/*
		readsize=0;
		c=fgetc(samfile);
		while(c!='\t' && c!=EOF){ // get SEQ
			if(c>=97 && c<=122) readchars[readsize]=(char)(c-32);
			else readchars[readsize]=c;
			readsize++;
			if(readsize==maxreadsize){ // increase read array size if needed
				maxreadsize+=readsize;
				readchars=(char *)realloc(readchars,(maxreadsize+1)*sizeof(char));
			}
			c=fgetc(samfile);
		}
		readchars[readsize]='\0';
		*/
		while(c!='\n' && c!=EOF) c=fgetc(samfile); // skip QUAL
		readpos-=16; // allow a relaxation of 16 chars behind
		if(readpos<0) readpos=0;
		m+=readsize; // size of the region in the ref
		m+=32; // relaxation of 16 chars on both ends
		if((readpos+m)>refsize) m=(refsize-readpos);
		numerrors=BasicDynamicProgramming((char *)(refchars+readpos),m,readchars,readsize);
		if(!sameread){ // new read
			numalignedreads++;
			totalnumchars+=readsize;
			totalnumerrors+=numerrors;
			prevnumerrors=numerrors;
		} else if(numerrors<prevnumerrors){ // same read but with a lowest number of errors
			totalnumerrors-=prevnumerrors;
			totalnumerrors+=numerrors;
			prevnumerrors=numerrors;
		}
	}
	printf("(%d reads) OK\n",numreads);
	fclose(readsfile);
	fclose(samfile);
	free(refchars);
	free(readchars);
	free(readname);
	printf(":: %.2lf%% mapped reads (%d of %d)\n",(((double)numalignedreads)/((double)realnumreads)*100.0),numalignedreads,realnumreads);
	printf(":: %.2lf errors/read (%.0lf bp/error)\n",((double)totalnumerrors)/((double)numalignedreads),((double)totalnumchars)/((double)totalnumerrors));
	printf("> Done!\n");
	exit(0);
}
