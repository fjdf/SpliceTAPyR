#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <math.h>
#include "version.h"
#include "tools.h"

#ifdef _MSC_VER
	#pragma warning(disable:4996)
#endif

typedef struct _MappedRead {
	char readName[256];
	char strand;
	char refName[256];
	unsigned int pos;
	unsigned int mapQual;
	unsigned int readSize;
	char *readChars;
	unsigned char *baseQuals;
	unsigned int extReadSize;
	char *extReadChars;
	char *extCigarCodes;
	unsigned char *extBaseQuals;
	unsigned char numHits;
} MappedRead;

static unsigned int maxReadSize; // "static" keyword to prevent compilation warning about a variable with the same name on another file
unsigned int maxExtReadSize;

typedef struct _VCFLineRecord {
	char *chrom;
	int pos;
	char *id;
	char *ref;
	char *alt;
	int qual;
	char filter_r50;
	char filter_q50;
	char filter_rs1;
	int info_dp;
	int info_bq;
	int info_mq;
	int info_mq0;
} VCFLineRecord;

typedef struct _GFFLineRecord {
	char seqid[20];
	char source[20];
	char type[20];
	int start;
	int end;
	char score[10];
	char strand;
	char phase;
	char attributes[256];
} GFFLineRecord;

// variables needed by main function and ACE format
FILE *samfile;
unsigned int *readsIdsByPos; // array the same size as the reference; if at least one reads starts at a position, it has an not zero id there
unsigned short *readsStartCountsInPos; // array the same size as the number of positions with reads starting there; stores how many reads start at that position id
fpos_t **readsStartFilePos; // array the same size as the number of positions with reads starting there; each position stores an array of fpos structs for the location of the read inside the sam file

// variables only needed for generating ACE format
FILE *acefile;
char *acefilename;
unsigned char *contigchars, *contigrefcs, *contigcodes, *contigquals;
unsigned int extcontigsize, maxextcontigsize, numcontigreads, contigstartpos, contigendpos;
fpos_t aceheaderfpos;

// variables needed by main function and VCF format
unsigned int numrows, *rowpos;
MappedRead **rows;
char *normalizedrefname;

// variables only needed for generating VCF format
FILE *vcffile;
char *vcffilename;
char *altchars;
int numaltchars, maxnumaltchars;
fpos_t refcharfpos;
VCFLineRecord vcfline;

// variables only needed for generating GFF format
FILE *gfffile;
FILE *genesgfffile;
FILE *genestxtfile;
FILE *genessnpsfile;
char *genesgfffilename;
char *genestxtfilename;
char *genessnpsfilename;
GFFLineRecord **gfflines, *nextgffline;
int *genestartpos, *geneendpos;
int *numgenereads;
int *numgenechars;
int numfoundgenes;
int numtotalgenes;
int numgenes;
int maxnumgenes;

MappedRead *NewMappedRead(){
	MappedRead *read;
	read=(MappedRead *)calloc(1,sizeof(MappedRead));
	read->readChars=(char *)calloc((maxReadSize+1),sizeof(char));
	read->baseQuals=(unsigned char *)calloc((maxReadSize+1),sizeof(unsigned char));
	read->extReadChars=(char *)calloc((maxExtReadSize+1),sizeof(char));
	read->extCigarCodes=(char *)calloc((maxExtReadSize+1),sizeof(char));
	read->extBaseQuals=(unsigned char *)calloc((maxExtReadSize+1),sizeof(unsigned char));
	return read;
}

void FreeMappedRead(MappedRead *read){
	free(read->readChars);
	free(read->baseQuals);
	free(read->extReadChars);
	free(read->extCigarCodes);
	free(read->extBaseQuals);
	free(read);
}

void PrintMappedRead(FILE *file, MappedRead *read){
	fprintf(file,"%s %c %u %s %u %s\n",read->readName,read->strand,read->pos,read->refName,read->readSize,read->readChars);
}

// TODO: add (read->cigarSize) , (read->cigarChars) and (read->cigarCodes)
// TODO: add check to verify if malformed CIGAR string and QUAL chars do not consume more chars than those that exist in the read
void LoadMappedRead(MappedRead *read, FILE *readsfile, fpos_t readfpos){
	char c;
	unsigned int flag, i, j, n;
	fpos_t cigarfpos;
	fsetpos(readsfile,&readfpos); // get: QNAME, FLAG, RNAME, POS, MAPQ
	fscanf(readsfile,"%[^\t]\t%d\t%255s\t%d\t%d\t",(read->readName),&flag,(read->refName),&(read->pos),&(read->mapQual));
	if((read->mapQual)>255) (read->mapQual)=0;
	if(flag & 16) read->strand='-';
	else read->strand='+';
	fgetpos(readsfile,&cigarfpos);
	for(i=0;i<4;i++){ // skip: CIGAR, RNEXT, PNEXT, TLEN
		c=fgetc(readsfile);
		while(c!='\t' && c!=EOF) c=fgetc(readsfile);
	}
	n=0;
	c=fgetc(readsfile);
	while(c!='\t' && c!=EOF){
		read->readChars[n++]=c; // get: SEQ
		c=fgetc(readsfile);
	}
	read->readChars[n]='\0';
	read->readSize=n;
	i=0;
	c=fgetc(readsfile);
	if(c!='*'){
		while(c!='\t' && c!='\n' && c!=EOF){
			read->baseQuals[i++]=(unsigned char)(c-33); // get: QUAL ; convert from char to phred ( 33 = '!' )
			c=fgetc(readsfile);
		}
	}
	while(i<n) read->baseQuals[i++]=(unsigned char)93; // set default phred quality as max ( 33 + 93 = 196 = '~' )
	read->baseQuals[i]='\0';
	fsetpos(readsfile,&cigarfpos);
	i=0; // position in extended arrays
	j=0; // position in normal arrays
	while(fscanf(readsfile,"%d%c",&n,&c)>0){ // get: CIGAR
		if(c=='M' || c=='I' || c=='=' || c=='X'){ // match, mismatch or insertion
			while(n>0){
				read->extReadChars[i]=read->readChars[j];
				read->extBaseQuals[i]=read->baseQuals[j];
				read->extCigarCodes[i]=c;
				j++;
				i++;
				n--;
			}
		} else if(c=='D' || c=='N'){ // deletion or skipped
			while(n>0){
				read->extReadChars[i]='-';
				read->extBaseQuals[i]=0;
				read->extCigarCodes[i]=c;
				i++;
				n--;
			}
		} else if(c=='S'){ // soft clipping
			while(n>0){
				j++;
				n--;
			}
		} else if(c=='H' || c=='P'){ // hard clipping or padding
			while(n>0){
				n--;
			}
		}
	}
	read->extReadSize=i;
	read->extReadChars[i]='\0';
	read->extBaseQuals[i]='\0';
	read->extCigarCodes[i]='\0';
}

// TODO: write "contig" tag for each reference sequence
void PrintVCFHeader(char *reffilename, char *samfilename){
	time_t rawtime;
	struct tm *timeinfo;
	char timestring[9];
	time(&rawtime);
	timeinfo=localtime(&rawtime);
	strftime(timestring,9,"%Y%m%d",timeinfo);
	fprintf(vcffile,"##fileformat=VCFv4.1\n");
	fprintf(vcffile,"##fileDate=%s\n",timestring);
	fprintf(vcffile,"##source=TAPyR%s\n",VERSION);
	fprintf(vcffile,"##reference=file://%s\n",reffilename);
	//fprintf(vcffile,"##contig=<id=%s,length=%u>\n",refname,refsize);
	fprintf(vcffile,"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n");
	fprintf(vcffile,"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read depth at this position for this sample\">\n");
	fprintf(vcffile,"##INFO=<ID=BQ,Number=1,Type=Integer,Description=\"RMS base quality at this position\">\n");
	fprintf(vcffile,"##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"RMS mapping quality\">\n");
	fprintf(vcffile,"##INFO=<ID=MQ0,Number=1,Type=Integer,Description=\"Number of MAPQ == 0 reads covering this record\">\n");
	fprintf(vcffile,"##FILTER=<ID=r50,Description=\"Less than 50%% of average read coverage\">\n");
	fprintf(vcffile,"##FILTER=<ID=q50,Description=\"Less than 50%% of average mapping quality\">\n");
	fprintf(vcffile,"##FILTER=<ID=rs1,Description=\"Reads from only one of the strands\">\n");
	fprintf(vcffile,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n",samfilename);
	fprintf(vcffile,"%s\t0\t.\t.",normalizedrefname); // initialize first line
	fgetpos(vcffile,&refcharfpos);
	vcfline.filter_q50=(char)0; // initialize VCF variables
	vcfline.filter_r50=(char)0;
	vcfline.filter_rs1=(char)0;
	vcfline.info_bq=0;
	vcfline.info_dp=0;
	vcfline.info_mq=0;
	vcfline.info_mq0=0;
}

// TODO: print INFO fields for SNP/indel
// TODO: print polyploid data when present
void PrintVCFPosition(unsigned int refpos, char type, char refchar, char calledchar, unsigned int posqual){
	unsigned int i, k, mq;
	unsigned char bq;
	MappedRead *temprow;
	refpos++; // convert to 1-based position
	if(type=='M' || type=='X'){ // regular match or SNP
		fprintf(vcffile,"%s\t%u\t.\t%c",normalizedrefname,refpos,refchar); // CHROM POS ID REF
		fgetpos(vcffile,&refcharfpos); // save position after printing the ref char
		altchars[0]=calledchar; // save last used alt char
		numaltchars=1;
		vcfline.qual=(int)posqual;
		vcfline.info_dp=0;
		vcfline.info_bq=0;
		vcfline.info_mq=0;
		vcfline.info_mq0=0;
		for(i=0;i<numrows;i++){ // get INFO values for this position
			k=rowpos[i];
			temprow=rows[i];
			if( (temprow->extReadChars)[k] != calledchar ) continue;
			bq=(temprow->extBaseQuals[k]);
			mq=(temprow->mapQual);
			(vcfline.info_dp)++; // read depth
			(vcfline.info_bq)+=(int)(bq*bq); // rms of base quality
			if(mq!=0) (vcfline.info_mq)+=(int)(mq*mq); // rms of mapping quality
			else (vcfline.info_mq0)++; // number of reads with 0 mapping quality
		}
		(vcfline.info_bq)= (int) floor( sqrt( (double)((vcfline.info_bq)/(vcfline.info_dp)) ) ); // RMS (square root of the sum of squares divided by the number of values )
		(vcfline.info_mq)= (int) floor( sqrt( (double)((vcfline.info_mq)/(vcfline.info_dp)) ) );
		if(type=='X') fprintf(vcffile,"\t%c\t%d\tPASS\t",calledchar,(vcfline.qual)); // ALT QUAL FILTER
		else fprintf(vcffile,"\t.\t%d\t.\t",(vcfline.qual)); // ALT QUAL FILTER
	} else if(type=='D'){ // deletion
		fsetpos(vcffile,&refcharfpos); // go back in line
		fprintf(vcffile,"%c",refchar); // append to the ref chars the char that was deleted
		fgetpos(vcffile,&refcharfpos); // save file position after this new last written ref char
		fprintf(vcffile,"\t%.*s\t%d\tPASS\t",numaltchars,altchars,(vcfline.qual)); // ALT QUAL FILTER
	} else if(type=='I'){ // insertion
		altchars[numaltchars]=calledchar; // append to the alt chars the char that was inserted
		numaltchars++;
		if(numaltchars>maxnumaltchars){ // alloc more space for the alt chars array if needed
			maxnumaltchars=numaltchars;
			altchars=(char *)realloc(altchars,(maxnumaltchars+1)*sizeof(char));
		}
		fsetpos(vcffile,&refcharfpos); // go back in line
		fprintf(vcffile,"\t%.*s\t%d\tPASS\t",numaltchars,altchars,(vcfline.qual)); // ALT QUAL FILTER
	} else { // other operation
		fprintf(vcffile,"%s\t%u\t.\t%c\t.\t.\t.\t.\t.\t.\n",normalizedrefname,refpos,refchar); // no reads covering this position
		return;
	}
	fprintf(vcffile,"DP=%d;BQ=%d;MQ=%d;MQ0=%d\t.\t.\n",(vcfline.info_dp),(vcfline.info_bq),(vcfline.info_mq),(vcfline.info_mq0)); // INFO FORMAT SAMPLE
}

void PrintACEHeader(unsigned int totalnumreads){
	fprintf(acefile,"AS "); // AS <number of contigs> <total number of reads in ace file>
	fgetpos(acefile,&aceheaderfpos); // save this position to later print the correct number of contigs here
	fprintf(acefile,"%10u %u\n\n",1,totalnumreads); // unsigned ints have at most 10 digits
}

void UpdateACEHeader(unsigned int numcontigs){
	fsetpos(acefile,&aceheaderfpos);
	fprintf(acefile,"%10u",numcontigs); // update the ACE header at the beggining with the number of contigs in the whole file
}

// TODO: write .PHD file with read base qualities and paired-ends/mate-pairs info
void PrintACEContig(unsigned int contigid){
	MappedRead *read;
	unsigned int id, refpos, contigpos, paddedreadsize, paddedreadpos, maxinssize;
	unsigned int i, j, k, n;
	char c, contigname[17], timestring[25];
	fpos_t readsizefpos, currentfpos;
	time_t timet;
	struct tm *timeinfo;
	time(&timet);
	timeinfo=localtime(&timet);
	strftime(timestring,25,"%a %b %d %H:%M:%S %Y",timeinfo);
	sprintf(contigname,"contig%u",contigid);
	fprintf(acefile,"CO %s %u %u 1 U\n",contigname,extcontigsize,(numcontigreads+1)); // CO <contig name> <# of bases> <# of reads in contig> <# of base segments in contig> <U or C>
	i=0;
	for(j=0;j<extcontigsize;j++){ // print consensus chars
		c=contigchars[j];
		if(c=='-') c='*'; // change deletion chars
		fputc(c,acefile);
		i++;
		if(i==50){ // lines with 50 chars
			fputc('\n',acefile);
			i=0;
		}
	}
	if(i!=0) fprintf(acefile,"\n");
	fprintf(acefile,"\n");
	fprintf(acefile,"BQ\n");
	i=0;
	for(j=0;j<extcontigsize;j++){ // print consensus quals
		c=contigchars[j];
		if(c=='-' || c=='*') continue; // only print quals for unpadded chars
		fprintf(acefile," %u",(unsigned int)((contigquals[j]>99)?(99):(contigquals[j]))); // limit quality values to 99 max
		i++;
		if(i==50){ // lines with 50 quals
			fputc('\n',acefile);
			i=0;
		}
	}
	if(i!=0) fprintf(acefile,"\n");
	fprintf(acefile,"\n");
	fprintf(acefile,"AF %s U 1\n",normalizedrefname); // fake read with the padded original reference portion
	/*
	while((c=(*refname))!='\0'){
		if(c!=' ' && c!='\t') fputc(c,acefile);
		refname++;
	}
	*/
	read=NewMappedRead();
	contigpos=0; // position in the padded contig
	for(refpos=contigstartpos;refpos<=contigendpos;refpos++){ // search positions of the reference where this contig belongs, and look for reads there
		while(contigpos<extcontigsize && (contigcodes[contigpos]=='I' || contigcodes[contigpos]=='P') ) contigpos++; // skip insertions and padding to get to the pos in the contig corresponding to the pos in the ref
		if((id=readsIdsByPos[refpos])!=0){ // check if there are reads starting at this pos
			n=readsStartCountsInPos[id];
			for(i=0;i<n;i++){ // process all reads starting at this position
				LoadMappedRead(read,samfile,readsStartFilePos[id][i]);
				paddedreadpos=contigpos;
				k=0;
				while( (read->extCigarCodes)[k]=='I' && paddedreadpos>0 ){ // if the read starts with insertions, we need to go back on the padded contig
					paddedreadpos--;
					k++;
				}
				fprintf(acefile,"AF %s %c %u\n",(read->readName),(((read->strand)=='+')?'U':'C'),(paddedreadpos+1)); // AF <read name> <C or U> <padded start consensus position>
			}
		}
		contigpos++;
	}
	fprintf(acefile,"BS 1 %u %s\n",extcontigsize,normalizedrefname); // BS <padded start consensus position> <padded end consensus position> <read name>
	fprintf(acefile,"\n"); // print fake first read corresponding to the original reference portion (but padded)
	fprintf(acefile,"RD %s %u 0 0\n",normalizedrefname,extcontigsize); // RD <read name> <# of padded bases> <# of whole read info items> <# of read tags>
	i=0;
	for(j=0;j<extcontigsize;j++){ // print reference chars
		c=contigrefcs[j];
		if(c=='-') c='*'; // change insertion chars
		fputc(c,acefile);
		i++;
		if(i==50){ // lines with 50 chars
			fputc('\n',acefile);
			i=0;
		}
	}
	if(i!=0) fprintf(acefile,"\n");
	fprintf(acefile,"\n");
	fprintf(acefile,"QA 1 %u 1 %u\n",extcontigsize,extcontigsize); // QA <qual clipping start> <qual clipping end> <align clipping start> <align clipping end>
	fprintf(acefile,"DS VERSION: 1 TIME: %s CHEM: fasta\n\n",timestring);
	contigpos=0; // position in the padded contig
	for(refpos=contigstartpos;refpos<=contigendpos;refpos++){ // search positions of the reference where this contig belongs, and look for reads there
		while(contigpos<extcontigsize && (contigcodes[contigpos]=='I' || contigcodes[contigpos]=='P') ) contigpos++; // skip padding to get to the pos in the contig corresponding to the pos in the ref
		if((id=readsIdsByPos[refpos])==0){ // if there are no reads starting at this pos, go to next pos
			contigpos++;
			continue;
		}
		n=readsStartCountsInPos[id];
		for(i=0;i<n;i++){ // process all reads starting at this position
			LoadMappedRead(read,samfile,readsStartFilePos[id][i]);
			paddedreadpos=contigpos;
			k=0;
			while( (read->extCigarCodes)[k]=='I' && paddedreadpos>0 ){ // if the read starts with insertions, we need to go back on the padded contig
				paddedreadpos--;
				k++;
			}
			fprintf(acefile,"RD %s ",(read->readName)); // RD <read name> <# of padded bases> <# of whole read info items> <# of read tags>
			fgetpos(acefile,&readsizefpos); // save file position to later print the correct number of padded bases
			fprintf(acefile,"%10u 0 0\n",(read->readSize));
			paddedreadsize=0;
			k=0;
			for(j=0;j<(read->extReadSize);j++){
				maxinssize=0; // get the size of the maximum insertion at this pos, from the padded consensus info
				while( paddedreadpos<extcontigsize && (contigcodes[paddedreadpos]=='I' || contigcodes[paddedreadpos]=='P') ){
					maxinssize++;
					paddedreadpos++;
				}
				while((read->extCigarCodes[j])=='I'){ // print all contiguous inserted chars in the read
					fputc((read->extReadChars[j]),acefile);
					paddedreadsize++;
					k++;
					if(k==50){
						fputc('\n',acefile);
						k=0;
					}
					j++;
					if(maxinssize!=0) maxinssize--;
				}
				while(maxinssize>0){ // if we did not fill all the space for inserted chars, print padding chars
					fputc('*',acefile);
					paddedreadsize++;
					k++;
					if(k==50){
						fputc('\n',acefile);
						k=0;
					}
					maxinssize--;
				}
				c=(read->extReadChars[j]);
				if(c=='\0') break; // reached the end of the read
				if(c=='-') c='*'; // change deletion char
				fputc(c,acefile);
				paddedreadsize++;
				k++;
				if(k==50){
					fputc('\n',acefile);
					k=0;
				}
				paddedreadpos++; // advance position in consensus too
			} // end of processing all chars of read
			if(k!=0) fputc('\n',acefile);
			fputc('\n',acefile);
			fgetpos(acefile,&currentfpos); // save current file pos
			fsetpos(acefile,&readsizefpos);
			fprintf(acefile,"%10u",paddedreadsize); // write padded read size at its correct file pos
			fsetpos(acefile,&currentfpos); // restore current file pos
			fprintf(acefile,"QA 1 %u 1 %u\n",paddedreadsize,paddedreadsize); // QA <qual clipping start> <qual clipping end> <align clipping start> <align clipping end>
			fprintf(acefile,"DS VERSION: 1 TIME: %s CHEM: sam\n\n",timestring);
		} // end of processing all reads
		contigpos++;
	} // end of processing all chars of reference segment
	FreeMappedRead(read);
}

void LoadGFFLine(){
	char c;
	int i;
	c=fgetc(gfffile);
	while(c=='#'){ // skip header lines
		while(c!='\n' && c!=EOF) c=fgetc(gfffile);
		c=fgetc(gfffile);
	}
	ungetc(c,gfffile);
	i=0;
	while((c=fgetc(gfffile))!='\t' && c!=EOF) if(i<19) nextgffline->seqid[i++]=c;
	nextgffline->seqid[i]='\0';
	i=0;
	while((c=fgetc(gfffile))!='\t' && c!=EOF) if(i<19) nextgffline->source[i++]=c;
	nextgffline->source[i]='\0';
	i=0;
	while((c=fgetc(gfffile))!='\t' && c!=EOF) if(i<19) nextgffline->type[i++]=c;
	nextgffline->type[i]='\0';
	fscanf(gfffile,"%d\t%d\t",&(nextgffline->start),&(nextgffline->end));
	i=0;
	while((c=fgetc(gfffile))!='\t' && c!=EOF) if(i<9) nextgffline->score[i++]=c;
	nextgffline->score[i]='\0';
	fscanf(gfffile,"%c\t%c\t",&(nextgffline->strand),&(nextgffline->phase));
	i=0;
	while((c=fgetc(gfffile))!='\n' && c!=EOF) if(i<255) nextgffline->attributes[i++]=c;
	nextgffline->attributes[i]='\0';
}

// TODO: remove exclusive ref names "R?" for TAPyR
// TODO: enable "50% avg" calculation for minmapq and mincoverage options
// TODO: add GFF split
// TODO: all references at once!
// TODO: output genome "filling" percentage as well as coverage
// TODO: select which reference to process (prefix of name)
// TODO: free refs* arrays at end
void GenerateConsensus(char **argsvector, int argscount){
	FILE *reffile, *consensusfile, *contigsfile, *agpfile, *countsfile, *strandcountsfile[2];
	char *reffilename, *samfilename;
	char *consensusfilename, *contigsfilename, *agpfilename, *countsfilename, *strandcountsfilename[2], *outputfilename, *outfilearg;
	long long int filesize;
	unsigned int numrefs, *refsizes, refsize;
	fpos_t *refsfpos, samfilepos;
	char c, **refnames, *refname, *targetrefname, refchar, cigarcode;
	//char *refchars;
	unsigned int i, j, k, n, m, r, cid;
	unsigned int numreads, numuniquereads, numids, maxids;
	unsigned int maxrefnamesize, maxrefsize, cigarcount;
	unsigned int mincoverage, minmapq;
	unsigned int counts[6], mapqsum[6], fwdcount[6], revcount[6];
	long long totalnumchars, totalmapqsum, totalnumreads;
	unsigned int consensussize, gapsize, numagpcomps, numsnps, numunmappedpos;
	unsigned int maxrows;
	MappedRead *temprow;
	char *string;
	char charFromId[5], idFromChar[256];
	unsigned int numcontigs, contigsize, mincontigsize, maxcontigsize;
	long long int progresscounter, progressstep;
	unsigned char outputfasta, outputace, outputvcf, outputgff, outputcounts, outputstrandcounts;
	unsigned char fromtapyr, onefileperref, allrefsatonce, bothstrands;
	int ploidy, mincovarg, minmapqarg;
	char *prevreadname;
	unsigned char numreadhits, numreadhitsinallrefs, maxnumreadhits, sameread;
	unsigned int *prevreadhitsids;
	unsigned char **numreadhitsbyid;
	float normcount;

	//fpos_t test;
	
	/*
	char *readrefname;
	unsigned int *numreadsperref, *maxidsperref, *numidsperref, **readsIdsByPosPerRef;
	unsigned int prevref, *numreadhitsperref, *maxnumreadhitsperref, **prevreadhitsidsperref;
	long long int *totalnumcharsperref, *totalmapqsumperref, *totalnumreadsperref;
	unsigned short **readsStartCountsInPosPerRef;
	fpos_t ***readsStartFilePosPerRef;
	unsigned char ***numreadhitsbyid;
	*/
	if(argscount<4){
		printf("\nConsensus calling:\n");
		#ifdef __unix__
			printf("\x1B[1m");
		#endif
		printf(" %s C <fasta-file> <sam-file> <options>\n",argsvector[0]);
		#ifdef __unix__
			printf("\x1B[0m");
		#endif
		printf("\n");
		printf(" <fasta-file>\treference genome file (FASTA)\n");
		printf(" <sam-file>\tmapped reads file (SAM)\n");
		printf("\n");
		printf(" -F\toutput consensus and contigs in FASTA + AGP format [on]\n");
		printf(" -A\toutput consensus and contigs in ACE format [off]\n");
		printf(" -V\toutput consensus and contigs in VCF format [off]\n");
		printf(" -N\toutput number of reads aligned under each position [off]\n");
		printf(" -NS\toutput number of reads aligned under each position for each strand [off]\n");
		printf("\n");
		printf(" -C\tminimum coverage required to call SNP/indel [50%% of average coverage]\n");
		printf(" -Q\tminimum MAPQ required to call SNP/indel [50%% of average MAPQ]\n");
		printf(" -S\treads from both strands are required to call SNP/indel [off]\n");
		//printf(" -P\tploidy of the reference genome (for VCF output) [1]\n");
		printf("\n");
		printf(" -O\tchoose output base filename [mapped reads filename]\n");
		printf(" -R\tprocess only the reference with this name [all]\n");
		printf("\n");
		printf(" -G\tautomatic gene annotation using this reference GFF file [none]\n");
		printf("\n");
		printf("> Done!\n");
		exit(0);
	}
	reffilename=NULL;
	samfilename=NULL;
	outputfilename=NULL;
	for(i=2;i<(unsigned int)argscount;i++){ // get ref filename
		if(argsvector[i][0]!='-'){ // skip options args
			reffilename=argsvector[i];
			break;
		}
	}
	if(i==(unsigned int)argscount){
		printf("> ERROR: Missing FASTA reference filename argument\n");
		exit(-1);
	}
	for(i=(i+1);i<(unsigned int)argscount;i++){ // get sam filename
		if(argsvector[i][0]!='-'){ // skip options args
			samfilename=argsvector[i];
			break;
		}
	}
	if(i==(unsigned int)argscount){
		printf("> ERROR: Missing SAM mapped reads filename argument\n");
		exit(-1);
	}
	outputfasta = ParseArgument(argscount,argsvector,"F",0);
	outputace = ParseArgument(argscount,argsvector,"A",0);
	outputvcf = ParseArgument(argscount,argsvector,"V",0);
	outputgff = ParseArgument(argscount,argsvector,"G",0);
	outputcounts = ParseArgument(argscount,argsvector,"N",0);
	outputstrandcounts = ParseArgument(argscount,argsvector,"NS",0);
	bothstrands = ParseArgument(argscount,argsvector,"S",0);
	mincovarg = ParseArgument(argscount,argsvector,"C",1);
	minmapqarg = ParseArgument(argscount,argsvector,"Q",1);
	ploidy = ParseArgument(argscount,argsvector,"P",1);
	if( (outputfasta+outputace+outputvcf+outputcounts+outputstrandcounts+outputgff) == 0 ) outputfasta = 1;
	if(ploidy<1) ploidy=1;
	else if(ploidy>12) ploidy=12;
	outfilearg=NULL;
	k=ParseArgument(argscount,argsvector,"O",2);
	if(k!=(unsigned int)(-1)) outfilearg=argsvector[k]; // get user defined output filename
	targetrefname=NULL;
	k=ParseArgument(argscount,argsvector,"R",2);
	if(k!=(unsigned int)(-1)){ // get chosen reference name
		n=0; // get size of the string
		for(i=k;i<(unsigned int)argscount;i++){ // process arguments because name can be split by spaces
			j=0; // position inside argument
			if(i==k && argsvector[i][0]=='"') j=1; // skip left most quote if it exists
			else if(argsvector[i][0]=='-') break; // stop when reached next valid argument
			while(argsvector[i][j]!='\0'){ n++; j++; }
			if(j!=0 && argsvector[i][(j-1)]=='"') { n--; break; } // discard right most quote if it exists
		}
		string=(char *)malloc((n+1)*sizeof(char));
		n=0; // get string chars
		for(i=k;i<(unsigned int)argscount;i++){
			j=0;
			c=argsvector[i][0];
			if(i==k && c=='"') j=1;
			else if(c=='-') break;
			while((c=argsvector[i][j])!='\0'){ string[n++]=c; j++; }
			if(j!=0 && argsvector[i][(j-1)]=='"') { n--; break; }
		}
		string[n]='\0';
		targetrefname=NormalizeSeqName(string,1); // normalize name
		free(string);
		string=NULL;
	}
	allrefsatonce=0;
	printf("> Loading reference genome file <%s> ... ",reffilename);
	fflush(stdout);
	if((reffile=fopen(reffilename,"r"))==NULL){
		printf("\n> ERROR: Reference file not found\n");
		exit(-1);
	}
	if((c=fgetc(reffile))!='>'){
		printf("\n> ERROR: Not a file in the FASTA format\n");
		exit(-1);
	}
	refsizes=NULL;
	refsfpos=NULL;
	refnames=NULL;
	refname=NULL;
	numrefs=0;
	refsize=0;
	maxrefnamesize=255;
	while(c!=EOF){ // get number of references and their sizes
		if(c=='>'){ // new sequence
			if(numrefs==0){ // first sequence
				refsizes=(unsigned int *)malloc(1*sizeof(unsigned int));
				refnames=(char **)malloc(1*sizeof(char *));
				refsfpos=(fpos_t *)malloc(1*sizeof(fpos_t));
			} else {
				refsizes=(unsigned int *)realloc(refsizes,(numrefs+1)*sizeof(unsigned int));
				refnames=(char **)realloc(refnames,(numrefs+1)*sizeof(char *));
				refsfpos=(fpos_t *)realloc(refsfpos,(numrefs+1)*sizeof(fpos_t));
			}
			refnames[numrefs]=(char *)malloc((maxrefnamesize+1)*sizeof(char));
			refname=refnames[numrefs];
			fgetpos(reffile,&(refsfpos[numrefs])); // store location of start of sequence inside file
			i=0; // ref label size
			while((c=fgetc(reffile))!='\n' && c!=EOF) if(i<maxrefnamesize) refname[i++]=c; // store ref label
			refname[i]='\0';
			if(numrefs!=0) refsizes[(numrefs-1)]=refsize; // store size of previous sequence
			refsize=0; // reset sequence size
			numrefs++;
		} else if( (c>=65 && c<=90) || (c>=97 && c<=122) ) refsize++; // count alphabetic chars
		c=fgetc(reffile);
	}
	if(numrefs!=0) refsizes[(numrefs-1)]=refsize; // store size of last sequence
	maxrefsize=0;
	n=0; // combined size of all sequences
	for(i=0;i<numrefs;i++){
		n+=refsizes[i];
		if(refsizes[i]>maxrefsize) maxrefsize=refsizes[i];
	}
	if(n==0){
		printf("\n> ERROR: No valid characters in file\n");
		exit(-1);
	}
	printf("(");
	if(numrefs>1){ PrintUnsignedNumber(numrefs);printf(" sequences ; "); }
	PrintUnsignedNumber(n);printf(" basepairs) OK\n");
	onefileperref=0; // use global output file name for all refs
	outputfilename=samfilename;
	if(numrefs>1 && targetrefname==NULL) outfilearg=NULL; // do not allow choosing output file when multiple refs and no ref choosing
	if(outfilearg!=NULL) outputfilename=outfilearg; // set user defined output filename if needed
	if(targetrefname!=NULL && outfilearg==NULL) onefileperref=1; // if choosing a ref, but not choosing out file, use the ref name as the out file
	if(numrefs>1 && targetrefname==NULL) onefileperref=1; // if multiple refs, but not choosing any of them, output one distinct file per ref
	printf("> Loading mapped reads file <%s> ... ",samfilename);
	fflush(stdout);
	if((samfile=fopen(samfilename,"r"))==NULL){
		printf("\n> ERROR: Reads file not found\n");
		exit(-1);
	}
	fromtapyr=0; // not a sam file from TAPyR
	c=fgetc(samfile);
	while(c=='@'){ // get header lines to check for TAPYR in program name tag
		for(i=0;i<2;i++) if( (c=fgetc(samfile)) != ("PG"[i]) ) break; // find line starting with "@PG"
		if(i!=2){ // if it is not this line, check next one
			while(c!='\n' && c!=EOF) c=fgetc(samfile);
			c=fgetc(samfile);
			continue;
		}
		for(i=0;i<2;i++){ // go to 2nd tab
			c=fgetc(samfile);
			while(c!='\t' && c!=EOF) c=fgetc(samfile);
		}
		for(i=0;i<8;i++) if( (c=fgetc(samfile)) != ("PN:TAPyR"[i]) ) break; // check tag "PN"
		if(i==8) fromtapyr=1; // TAPYR tag was found
		while(c!='\n' && c!=EOF) c=fgetc(samfile); // check next line
		c=fgetc(samfile);
	}
	rewind(samfile);
	printf("(");
	fseek(samfile,0L,SEEK_END);
	filesize=(long long int)ftell(samfile);
	fseek(samfile,0L,SEEK_SET);
	PrintNumber(filesize);
	printf(" bytes) OK\n");
	k=ParseArgument(argscount,argsvector,"G",2);
	if(k!=(unsigned int)(-1)){
		printf("> Loading annotations file <%s> ... ",argsvector[k]);
		fflush(stdout);
		if((gfffile=fopen(argsvector[k],"r"))==NULL){
			printf("\n> ERROR: Annotations file not found\n");
			exit(-1);
		}
		for(i=0;i<5;i++) if( (c=fgetc(gfffile)) != ("##gff"[i]) ) break; // check tag "##gff"
		if(i!=5){
			printf("\n> ERROR: Not a file in the GFF or GTF formats\n");
			exit(-1);
		}
		nextgffline=(GFFLineRecord *)malloc(1*sizeof(GFFLineRecord)); // this variable is kept between refs when all GFFs are in a single file
		LoadGFFLine();
		printf("(");
		fseek(gfffile,0L,SEEK_END);
		PrintNumber((long long int)ftell(gfffile));
		fseek(gfffile,0L,SEEK_SET);
		printf(" bytes) OK\n");
	}
	charFromId[0]='A'; // initialize char lookup arrays
	charFromId[1]='C';
	charFromId[2]='G';
	charFromId[3]='T';
	charFromId[4]='-';
	for(i=0;i<256;i++) idFromChar[i]=(char)5;
	idFromChar[(int)'A']=(char)0;
	idFromChar[(int)'C']=(char)1;
	idFromChar[(int)'G']=(char)2;
	idFromChar[(int)'T']=(char)3;
	idFromChar[(int)'-']=(char)4;
	idFromChar[(int)'N']=(char)4;
	consensusfilename=NULL;
	contigsfilename=NULL;
	agpfilename=NULL;
	acefilename=NULL;
	vcffilename=NULL;
	genesgfffilename=NULL;
	genestxtfilename=NULL;
	genessnpsfilename=NULL;
	countsfilename=NULL;
	strandcountsfilename[0]=NULL;
	strandcountsfilename[1]=NULL;
	consensusfile=NULL;
	contigsfile=NULL;
	agpfile=NULL;
	acefile=NULL;
	vcffile=NULL;
	genesgfffile=NULL;
	genestxtfile=NULL;
	genessnpsfile=NULL;
	countsfile=NULL;
	strandcountsfile[0]=NULL;
	strandcountsfile[1]=NULL;
	numreadhitsbyid=NULL;
	prevreadhitsids=NULL;
	prevreadname=NULL;
	numreads=0;
	numreadhits=0;
	numreadhitsinallrefs=0;
	maxnumreadhits=0;
	totalnumreads=0;
	sameread=0;
	numids=0;
	r=0;
	if(!onefileperref){ // one global file for all refs
_create_output_files:
		if(outputfasta){
			consensusfilename=AppendToBasename(outputfilename,"-consensus.fasta");
			if((consensusfile=fopen(consensusfilename,"w"))==NULL){
				printf("> ERROR: Cannot create consensus file\n");
				exit(-1);
			}
			contigsfilename=AppendToBasename(outputfilename,"-contigs.fasta");
			if((contigsfile=fopen(contigsfilename,"w"))==NULL){
				printf("> ERROR: Cannot create contigs file\n");
				exit(-1);
			}
			agpfilename=AppendToBasename(outputfilename,"-contigs.agp");
			if((agpfile=fopen(agpfilename,"w"))==NULL){
				printf("> ERROR: Cannot create contigs file\n");
				exit(-1);
			}
			fprintf(agpfile,"##agp-version 2.0\n"); // print AGP header
			fprintf(agpfile,"#%s\n",reffilename);
			fprintf(agpfile,"#%s\n",samfilename);
		}
		if(outputcounts){
			countsfilename=AppendToBasename(outputfilename,"-counts.txt");
			if((countsfile=fopen(countsfilename,"w"))==NULL){
				printf("> ERROR: Cannot create CSV file\n");
				exit(-1);
			}
		}
		if(outputstrandcounts){
			strandcountsfilename[0]=AppendToBasename(outputfilename,"-fwd_counts.txt");
			if((strandcountsfile[0]=fopen(strandcountsfilename[0],"w"))==NULL){
				printf("> ERROR: Cannot create CSV file\n");
				exit(-1);
			}
			strandcountsfilename[1]=AppendToBasename(outputfilename,"-rev_counts.txt");
			if((strandcountsfile[1]=fopen(strandcountsfilename[1],"w"))==NULL){
				printf("> ERROR: Cannot create CSV file\n");
				exit(-1);
			}
		}
		if(outputace){
			acefilename=AppendToBasename(outputfilename,".ace");
			if((acefile=fopen(acefilename,"w"))==NULL){
				printf("> ERROR: Cannot create ACE file\n");
				exit(-1);
			}
		}
		if(outputvcf){
			vcffilename=AppendToBasename(outputfilename,".vcf");
			if((vcffile=fopen(vcffilename,"w"))==NULL){
				printf("> ERROR: Cannot create VCF file\n");
				exit(-1);
			}
			PrintVCFHeader(reffilename,samfilename);
		}
		if(outputgff){
			genesgfffilename=AppendToBasename(outputfilename,"-genes.gff");
			if((genesgfffile=fopen(genesgfffilename,"w"))==NULL){
				printf("> ERROR: Cannot create genes GFF file\n");
				exit(-1);
			}
			fprintf(genesgfffile,"##gff-version 3\n"); // print GFF header
			fprintf(genesgfffile,"#!processor TAPyR automatic annotation tool\n");
			genestxtfilename=AppendToBasename(outputfilename,"-genes.txt");
			if((genestxtfile=fopen(genestxtfilename,"w"))==NULL){
				printf("> ERROR: Cannot create genes CSV file\n");
				exit(-1);
			}
			genessnpsfilename=AppendToBasename(outputfilename,"-snps.txt");
			if((genessnpsfile=fopen(genessnpsfilename,"w"))==NULL){
				printf("> ERROR: Cannot create SNPs CSV file\n");
				exit(-1);
			}
		}
		if(onefileperref) goto _process_ref; // we could have got here through the _create_output_files tag
	}
	printf("> Using parameters: Ploidy=%d ; StrandEvidence=%s\n",ploidy,((bothstrands)?("both"):("single")));	
	for(r=0;r<numrefs;r++){
		refsize=refsizes[r]; // set ref size
		fsetpos(reffile,&(refsfpos[r])); // restore ref file pos
		maxrefnamesize=255;
		refname=(char *)malloc((maxrefnamesize+1)*sizeof(char));
		i=0; // ref label size
		while((c=fgetc(reffile))!='\n' && c!='\r' && c!=EOF) if(i<maxrefnamesize) refname[i++]=c; // store ref label
		refname[i]='\0';
		normalizedrefname=NormalizeSeqName(refname,1);
		if(targetrefname!=NULL){ // if we are looking for a specific reference, check its name
			i=0;
			while((c=targetrefname[i])!='\0'){ // check if target name is a prefix of current name 
				if(c!=normalizedrefname[i]) break;
				i++;
			}
			if(c!='\0') continue; // if the name did not match, check next reference
		}
		fgetpos(reffile,&(refsfpos[r])); // save file pos for beginning of ref chars
		printf("> Processing sequence \"%s\"",refname);
		if(fromtapyr){ // if it's a sam file created by TAPYR, update the ref names to ref ids
			if(numrefs>1) sprintf(refname,"R%d",(r+1));
			else sprintf(refname,"R");
			printf(" (%s)",refname);
		}
		printf(" ... ");
		fflush(stdout);
		//refchars=(char *)malloc((refsize+1)*sizeof(char));
		readsIdsByPos=(unsigned int *)calloc((refsize+1),sizeof(unsigned int)); // if a ref pos has one or more reads starting there, that pos has an non-zero id ; +1 to allow insertions after last pos
		//if( refchars==NULL || readsIdsByPos==NULL ){
		if( readsIdsByPos==NULL ){
			printf("> ERROR: Not enough memory to allocate arrays\n");
			exit(-1);
		}
		//for(i=0;i<=refsize;i++) readsIdsByPos[i]=0; // reset read start ids ; +1 to allow insertions after last pos
		maxids=1;
		numids=0; // last used id (0 is never used for storing reads)
		readsStartCountsInPos=(unsigned short *)malloc(maxids*sizeof(unsigned short)); // number of reads starting at the ref pos with that id ; id=0 is never used
		readsStartFilePos=(fpos_t **)malloc(maxids*sizeof(fpos_t *)); // for each id, the file pos in the reads file of the reads starting on that id ; id=0 is never used
		if(outputcounts || outputstrandcounts){ // both these arrays are needed while processing the reads file
			prevreadname=(char *)malloc((255+1)*sizeof(char));
			for(i=0;i<=255;i++) prevreadname[i]=0;
			numreadhits=0; // number of current hits for the current read
			numreadhitsinallrefs=0;
			maxnumreadhits=1;
			numreadhitsbyid=(unsigned char **)malloc(maxids*sizeof(unsigned char *)); // array of number of read hits for each of the reads starting at the ref pos with that id
			prevreadhitsids=(unsigned int *)malloc(maxnumreadhits*sizeof(unsigned int)); // ids where the other occurrences of the same read were placed
			if( numreadhitsbyid==NULL ){
				printf("> ERROR: Not enough memory to allocate arrays\n");
				exit(-1);
			}
		}
		if(outputace){
			maxextcontigsize=refsize; // the extended arrays need extra space for possible insertion, deletion or padding chars
			contigchars=(unsigned char *)malloc((maxextcontigsize+1)*sizeof(unsigned char));
			contigrefcs=(unsigned char *)malloc((maxextcontigsize+1)*sizeof(unsigned char));
			contigcodes=(unsigned char *)malloc((maxextcontigsize+1)*sizeof(unsigned char));
			contigquals=(unsigned char *)malloc((maxextcontigsize+1)*sizeof(unsigned char)); // unsigned to support quals up to 255
			if( contigchars==NULL || contigrefcs==NULL || contigcodes==NULL || contigquals==NULL ){
				printf("> ERROR: Not enough memory to allocate arrays\n");
				exit(-1);
			}
		}
		if(outputvcf){
			maxnumaltchars=1;
			numaltchars=0;
			altchars=(char *)malloc((maxnumaltchars+1)*sizeof(char)); // variable to save the chars of the (possibly growing) ALT string
			altchars[0]='\0';
		}
		if(outputgff){
			maxnumgenes=1;
			genestartpos=(int *)malloc(maxnumgenes*sizeof(int));
			geneendpos=(int *)malloc(maxnumgenes*sizeof(int));
			numgenereads=(int *)malloc(maxnumgenes*sizeof(int));
			gfflines=(GFFLineRecord **)malloc(maxnumgenes*sizeof(GFFLineRecord *));
			gfflines[(maxnumgenes-1)]=(GFFLineRecord *)malloc(1*sizeof(GFFLineRecord));
			numgenes=0;
			numfoundgenes=0;
			numtotalgenes=0;
			/*
			GFFLineRecord **gfflines, *nextgffline;
			int *genestartpos, *geneendpos;
			int *numgenereads;
			int *numgenechars;
			int numfoundgenes;
			int numtotalgenes;
			int numgenes;
			int maxnumgenes;
			*/
		}
		i=0; // ref size
		k=0; // number of invalid chars
		while((c=fgetc(reffile))!='>' && c!=EOF){ // get ref chars until next sequence or end of file
			if(c>=97 && c<=122) c=(char)(c-32); // convert to uppercase
			if(c>=65 && c<=90){ // valid alphabet uppercase character
				if(c!='A' && c!='C' && c!='G' && c!='T'){ // invalid char
					c='N';
					k++;
				}
				//refchars[i++]=c;
				i++;
			}
		}
		//refchars[i]='\0';
		printf("(");PrintUnsignedNumber(i);printf(" bp");
		if(k!=0){ printf(" ; ");PrintUnsignedNumber(k);printf(" N's"); }
		printf(") OK\n");
		printf("> Processing reads");
		fflush(stdout);
		progressstep=(filesize/10);
		progresscounter=0;
		maxReadSize=0;
		maxExtReadSize=0;
		numreads=0; // number of mapped reads (including multiple hits) for this ref
		numuniquereads=0; // number of unique mapped reads (excluding multiple hits) for this ref
		totalnumchars=0;
		totalmapqsum=0;
		totalnumreads=0; // number of mapped reads (including multiple hits) inside SAM file for all refs
		rewind(samfile);
		c=fgetc(samfile);
		while(c=='@'){ // skip header lines
			while(c!='\n' && c!=EOF) c=fgetc(samfile);
			c=fgetc(samfile);
		}
		ungetc(c,samfile);
		while(c!=EOF){ // loop for all reads inside SAM file
			progresscounter=(long long int)ftell(samfile);
			if(progresscounter>=progressstep){ // print progress dots
				printf(".");
				fflush(stdout);
				progressstep+=progressstep;
			}
			fgetpos(samfile,&samfilepos); // save position inside file of the beginning of this read
			if(outputcounts || outputstrandcounts){ // check if this is another occurrence of the previous read or not
				i=0;
				sameread=(unsigned char)1;
				while((c=fgetc(samfile))!='\t' && c!=EOF){ // compare read name with previous read name
					sameread &= (unsigned char)( prevreadname[i] == c );
					prevreadname[i]=c;
					if(i!=255) i++;
				}
				if(!sameread){ // this check is necessary because if the first hits are in other refs, when we get to a hit in the current ref, sameread is true but numreadhits still has the value for the previous read
					numreadhits=0;
					numreadhitsinallrefs=1; // first hit of a new read, no matter what reference (numreadhits will only be incremented next if we detect its the current seq)
				} else {
					if(numreadhitsinallrefs!=UCHAR_MAX) numreadhitsinallrefs++; // one more hit for an existing read, no matter what reference
				}
				c=fgetc(samfile);
				while(c!='\t' && c!=EOF) c=fgetc(samfile); // skip flag
			} else for(i=0;i<2;i++){ // advance to 2nd tab (reference label)
				c=fgetc(samfile);
				while(c!='\t' && c!=EOF) c=fgetc(samfile);
			}
			c=fgetc(samfile); // get first letter of reference
			if(c=='*' || c==EOF){ // if the read was not mapped, advance to the next one
				while(c!='\n' && c!=EOF) c=fgetc(samfile);
				continue;
			}
			totalnumreads++; // count total number of mapped reads on all refs
			/*
			if(allrefsatonce){ // if we are processing all refs at once, check which one this read belongs to
				i=0;
				ungetc(c,samfile);
				while((c=fgetc(samfile))!='\t' && c!=EOF) if(i<maxrefnamesize) readrefname[i++]=c; // store ref name where read was mapped
				readrefname[i]='\0';
				for(n=0;n<numrefs;n++){
					refname=refnames[n];
					i=0;
					while(refname[i]==readrefname[i] && refname[i]!='\0') i++; // compare ref names
					if(refname[i]=='\0' && readrefname[i]=='\0') break; // if both names are equal, this is the ref
				}
				if(n==numrefs){ // not matching ref were found
					while(c!='\n' && c!=EOF) c=fgetc(samfile);
					continue;
				}
				// load variables for this ref

			} else
			*/
			if(numrefs>1){ // if there is more than one ref, check if this read is mapped in the current ref
				i=0; // current ref label position
				ungetc(c,samfile);
				while((c=fgetc(samfile))!='\t' && c!=EOF) if(c!=refname[i++]) break; // compare ref label chars
				if(c!='\t' || refname[i]!='\0'){ // if the refs do not match, advance to next read
					while(c!='\n' && c!=EOF) c=fgetc(samfile);
					continue;
				}
			} else for(i=0;i<1;i++){ // advance one more tab (position)
				c=fgetc(samfile);
				while(c!='\t' && c!=EOF) c=fgetc(samfile); // skip ref label
			}
			j=0;
			fscanf(samfile,"%u",&j); // get read mapped position
			if(j==0 || j>refsize){ // if the position is invalid, go to next read
				while(c!='\n' && c!=EOF) c=fgetc(samfile);
				continue;
			}
			j--; // fix SAM's 1-based location
			for(i=0;i<1;i++){ // advance 1 more tab (MAPQ)
				c=fgetc(samfile);
				while(c!='\t' && c!=EOF) c=fgetc(samfile);
			}
			n=0;
			fscanf(samfile,"%u",&n); // get mapping quality
			if(n>255) n=0; // limit value
			totalmapqsum+=n;
			for(i=0;i<1;i++){ // advance 1 more tab (CIGAR string)
				c=fgetc(samfile);
				while(c!='\t' && c!=EOF) c=fgetc(samfile);
			}
			c=fgetc(samfile); // get first letter of the CIGAR string
			if(c=='*' || c==EOF){ // if the CIGAR string is not available, advance to the next read
				while(c!='\n' && c!=EOF) c=fgetc(samfile);
				continue;
			}
			ungetc(c,samfile);
			n=0;
			while(fscanf(samfile,"%d%c",&cigarcount,&cigarcode)>0) n+=cigarcount; // get size of the whole CIGAR string
			if(n>maxExtReadSize) maxExtReadSize=n;
			for(i=0;i<3;i++){ // advance 3 more tabs (sequence) (should be 4 but fscanf already consumed one '\t')
				c=fgetc(samfile);
				while(c!='\t' && c!=EOF) c=fgetc(samfile);
			}
			n=0;
			while((c=fgetc(samfile))!='\t' && c!=EOF) n++; // get size of the read
			if(n>maxReadSize) maxReadSize=n;
			totalnumchars+=n;
			if(c==EOF) continue;
			numreads++; // one more read for this ref
			if(!sameread) numuniquereads++;
			if(readsIdsByPos[j]==0){ // still no reads starting at this position
				numids++; // create new id
				if(numids==maxids){ // realloc more space for arrays if needed
					maxids+=1024; // increase size in steps of 1024
					//maxids++;
					readsStartCountsInPos=(unsigned short *)realloc(readsStartCountsInPos,maxids*sizeof(unsigned short));
					readsStartFilePos=(fpos_t **)realloc(readsStartFilePos,maxids*sizeof(fpos_t *));
				} // maxids should always be at least equal to (numids+1) because it has the 0-th position
				i=numids;
				readsIdsByPos[j]=i;
				readsStartCountsInPos[i]=1;
				readsStartFilePos[i]=(fpos_t *)malloc(1*sizeof(fpos_t));
			} else { // there are already some reads starting at this position
				i=readsIdsByPos[j]; // existent id
				if(readsStartCountsInPos[i]!=USHRT_MAX) readsStartCountsInPos[i]++;
				readsStartFilePos[i]=(fpos_t *)realloc(readsStartFilePos[i],(readsStartCountsInPos[i])*sizeof(fpos_t));
			}
			readsStartFilePos[i][(readsStartCountsInPos[i]-1)]=samfilepos; // save read location inside sam file
			if(outputcounts || outputstrandcounts){
				if(readsStartCountsInPos[i]==1){ // new id was created
					if(maxids==(numids+1024)) numreadhitsbyid=(unsigned char **)realloc(numreadhitsbyid,maxids*sizeof(unsigned char *));
					numreadhitsbyid[i]=(unsigned char *)malloc(1*sizeof(unsigned char));
				} else { // there is one more read for an existing id
					numreadhitsbyid[i]=(unsigned char *)realloc(numreadhitsbyid[i],(readsStartCountsInPos[i])*sizeof(unsigned char));
				}
				if(sameread){ // another hit for the same read
					for(n=0;n<numreadhits;n++){ // update count on prev hits
						j=(prevreadhitsids[n]); // id of prev hit
						//numreadhitsbyid[j][(readsStartCountsInPos[j]-1)]++; // our read was the last one added to that id
						numreadhitsbyid[j][(readsStartCountsInPos[j]-1)] = numreadhitsinallrefs;
					}
					if(numreadhits!=UCHAR_MAX) numreadhits++; // one more hit
					if(numreadhits>maxnumreadhits){ // realloc array of prev hits' ids
						maxnumreadhits=numreadhits;
						prevreadhitsids=(unsigned int *)realloc(prevreadhitsids,maxnumreadhits*sizeof(unsigned int));
					}
				} else numreadhits=1; // new read
				prevreadhitsids[(numreadhits-1)]=i; // save id used for this hit
				//numreadhitsbyid[i][(readsStartCountsInPos[i]-1)]=numreadhits; // set read count
				numreadhitsbyid[i][(readsStartCountsInPos[i]-1)] = numreadhitsinallrefs;
			}
			while(c!='\n' && c!=EOF) c=fgetc(samfile); // advance to next read
		} // end of loop for all reads inside SAM file
		if(numreads==0){
			printf("\n> WARNING: No valid reads found for this reference\n");
			goto _free_and_next_ref;
		}
		printf("(");PrintUnsignedNumber(numuniquereads);
		printf(" distinct reads (%.2lf%%) ; ",((double)(numuniquereads*100)/(double)totalnumreads));PrintUnsignedNumber(numreads);
		printf(" read hits ; %.2lfx coverage ; %u avg MAPQ) OK\n",((double)totalnumchars/(double)refsize),((unsigned int)totalmapqsum/(unsigned int)numreads));
		if(onefileperref){
			outputfilename=normalizedrefname; // prefix output file names with the name of the sequence
			goto _create_output_files;
		}
_process_ref:
		if(mincovarg>0) mincoverage=mincovarg;
		else mincoverage=1;
		//else mincoverage=(unsigned int)(totalnumchars/(long long)((ploidy+1)*refsize)); // half the coverage
		if(minmapqarg>0) minmapq=minmapqarg;
		else minmapq=0;
		//else minmapq=(unsigned int)(totalmapqsum/(long long)(2*numreads)); // half the avg mapq
		printf("> Using extra parameters: MinimumReadCoverage=%u ; MinimumMappingQuality=%u\n",mincoverage,minmapq);
		printf("> Building consensus");
		fflush(stdout);
		progressstep=(long long int)(refsize/10);
		progresscounter=0;
		if(outputfasta) fprintf(consensusfile,">%s\n",refname);
		if(outputace) PrintACEHeader(numreads); // print ACE file header
		maxrows=1;
		rowpos=(unsigned int *)malloc(maxrows*sizeof(unsigned int));
		rows=(MappedRead **)malloc(maxrows*sizeof(MappedRead *));
		for(i=0;i<maxrows;i++){
			rows[i]=NewMappedRead();
			rows[i]->extReadSize=0;
			rowpos[i]=0;
		}
		consensussize=0;
		gapsize=0;
		numagpcomps=0;
		numsnps=0;
		numunmappedpos=0;
		numcontigs=0;
		contigsize=0;
		mincontigsize=refsize;
		maxcontigsize=0;
		extcontigsize=0;
		numcontigreads=0;
		contigstartpos=0;
		contigendpos=0;
		numrows=0;
		fsetpos(reffile,&(refsfpos[r])); // restore file pos for beggining of ref chars
		for(j=0;j<=refsize;j++){ // process all positions of the reference ; +1 to allow insertions after last pos
			while(1){ // get next ref char
				refchar=fgetc(reffile);
				if(refchar>=97 && refchar<=122) refchar=(char)(refchar-32); // convert to uppercase
				if(refchar>=65 && refchar<=90){ // valid alphabet uppercase character
					if(refchar!='A' && refchar!='C' && refchar!='G' && refchar!='T') refchar='N'; // invalid char
					break;
				} else if(refchar=='>' || refchar==EOF) break; // would get here if checking for insertions after last pos
			}
			if(progresscounter==progressstep){ // print progress dots
				printf(".");
				fflush(stdout);
				progresscounter=0;
			} else progresscounter++;
			i=0;
			while(i<numrows){ // check if any of the reads has ended
				if(rowpos[i]>=(rows[i]->extReadSize)){
					temprow=rows[i]; // switch this row with last row
					rows[i]=rows[(numrows-1)];
					rows[(numrows-1)]=temprow;
					rowpos[i]=(rowpos[(numrows-1)]);
					numrows--;
					/*
					// NOTE: only needed if we want to keep the relative order of the reads, e.g. in the PileUp format
					k=i;
					temprow=rows[k]; // this row will be moved to the last position
					while((k+1)<numrows){ // shift all next rows up one position
						rows[k]=rows[(k+1)];
						rowpos[k]=rowpos[(k+1)];
						k++;
					}
					rows[(numrows-1)]=temprow;
					numrows--;
					*/
				} else i++; // if the row at this array pos ended, check the pos again (because it is a different one now)
			}
			if((k=readsIdsByPos[j])!=0){ // check if there are any reads starting at this ref pos
				n=readsStartCountsInPos[k];
				for(i=0;i<n;i++){
					numrows++;
					if(numrows>maxrows){ // alloc space for one more row if needed
						maxrows=numrows;
						rowpos=(unsigned int *)realloc(rowpos,maxrows*sizeof(unsigned int));
						rows=(MappedRead **)realloc(rows,maxrows*sizeof(MappedRead *));
						rows[(numrows-1)]=NewMappedRead();
					}
					LoadMappedRead(rows[(numrows-1)],samfile,readsStartFilePos[k][i]); // place read after the last filled row
					rowpos[(numrows-1)]=0; // reset pos inside row arrays
				}
				numcontigreads+=n;
			}
			if( (outputcounts || outputstrandcounts) && (j!=refsize) ){
				if(k!=0){ // if new reads were added
					n=readsStartCountsInPos[k]; // number of new reads added
					m=(numrows-n); // first new read
					for(i=0;i<n;i++) (rows[(m+i)]->numHits) = (numreadhitsbyid[k][i]); // set number of hits on newly added reads
				}
				if(outputcounts){
					normcount=0.0; // if a read occurs N times, it's contribution is (1/N)
					for(i=0;i<numrows;i++) normcount += ((float)1.0)/((float)(rows[i]->numHits));
					fprintf(countsfile,"%u,%u,%.3f\n",(j+1),(numrows),(normcount)); // position , number of reads , normalized number of reads
				}
				if(outputstrandcounts){
					n=0; // counts for forward strand
					normcount=0.0;
					for(i=0;i<numrows;i++){
						if((rows[i]->strand)=='+'){
							n++;
							normcount += ((float)1.0)/((float)(rows[i]->numHits));
						}
					}
					fprintf(strandcountsfile[0],"%u,%u,%.3f\n",(j+1),(n),(normcount));
					n=0; // counts for reverse strand
					normcount=0.0;
					for(i=0;i<numrows;i++){
						if((rows[i]->strand)=='-'){
							n++;
							normcount += ((float)1.0)/((float)(rows[i]->numHits));
						}
					}
					fprintf(strandcountsfile[1],"%u,%u,%.3f\n",(j+1),(n),(normcount));
				}
			}
			if(numrows==0){ // no reads covering this ref position
				if(j==refsize) continue; // last pos only counted if there were insertions
				if(outputfasta) fputc('N',consensusfile);
				numunmappedpos++;
				if(outputvcf) PrintVCFPosition(j,'N',refchar,'.',0);
				if(contigsize>0){ // save the already started contig
					consensussize+=contigsize;
					if(contigsize<mincontigsize) mincontigsize=contigsize;
					if(contigsize>maxcontigsize) maxcontigsize=contigsize;
					contigendpos=(j-1);
					if(outputfasta){
						fputc('\n',contigsfile);
						fprintf(agpfile,"%s\t%u\t%u\t%u\tD\tcontig%u\t1\t%u\t+\n",refname,(contigstartpos+1),(contigendpos+1),numagpcomps,numcontigs,contigsize);
					}
					if(outputace) PrintACEContig(numcontigs);
					contigsize=0; // reset contig counters
					extcontigsize=0;
					numcontigreads=0;
					numagpcomps++; // new gap
					gapsize=0; // reset size of new gap between contigs (this last one and the next one)
				}
				gapsize++;
				continue; // go to next ref pos
			}
			m=0; // size of the largest insertion
			for(i=0;i<numrows;i++){ // check if there are any insertions starting at this ref position
				k=rowpos[i];
				string=(rows[i]->extCigarCodes);
				n=0;
				while(string[k]=='I'){
					n++;
					k++;
				}
				if(n>m) m=n;
			}
			while(m>0){ // if insertions were found
				for(i=0;i<5;i++){ // reset char counters
					counts[i]=0;
					mapqsum[i]=0;
					fwdcount[i]=0;
					revcount[i]=0;
				}
				n=0; // number of columns with insertions
				for(i=0;i<numrows;i++){ // get count of all chars from rows with insertions in this column
					k=rowpos[i];
					temprow=rows[i];
					cigarcode=(temprow->extCigarCodes)[k];
					if(cigarcode!='I') continue;
					c=(temprow->extReadChars)[k];
					/*
					if(c=='A') cid=0;
					else if(c=='C') cid=1;
					else if(c=='G') cid=2;
					else if(c=='T') cid=3;
					else if(c=='-') cid=4;
					*/
					cid=(unsigned int)idFromChar[(int)c];
					counts[cid]++;
					mapqsum[cid]+=(temprow->mapQual);
					if((temprow->strand)=='+') fwdcount[cid]++;
					else revcount[cid]++;
					rowpos[i]++; // advance position in this row
					n++;
				}
				if(n<=(numrows/2)){ // too few columns (<=50%) with insertions to create an indel
					for(i=0;i<numrows;i++){ // advance counters in all rows with insertions to next non-insertion
						string=(rows[i]->extCigarCodes);
						while(string[(rowpos[i])]=='I') rowpos[i]++;
					}
					break;
				}
				cid=0; // get the index of the most frequent char
				n=counts[0];
				for(i=1;i<5;i++){
					if(counts[i]>n){
						n=counts[i];
						cid=i;
					}
				}
				/*
				if(cid==0) c='A'; // get the letter of the most frequent char
				else if(cid==1) c='C';
				else if(cid==2) c='G';
				else if(cid==3) c='T';
				else if(cid==4) c='-';
				*/
				c=charFromId[cid]; // get the letter of the most frequent char
				mapqsum[cid]=(mapqsum[cid]/counts[cid]); // average MAPQ sum over all rows with this insertion
				if( counts[cid]<mincoverage || mapqsum[cid]<minmapq || (bothstrands && (fwdcount[cid]==0 || revcount[cid]==0) ) ){ // if the requirements to create an indel are not satisfied
					for(i=0;i<numrows;i++){ // advance counters in all rows with insertions to next non-insertion
						string=(rows[i]->extCigarCodes);
						while(string[(rowpos[i])]=='I') rowpos[i]++;
					}
					break;
				}
				if(outputfasta) fputc(c,consensusfile);
				if(contigsize==0){ // if starting a new contig, print contig name
					numcontigs++; // update contig count
					contigstartpos=j; // 1st position of the contig in the reference
					if(outputfasta){
						if(gapsize>0) fprintf(agpfile,"%s\t%u\t%u\t%u\tN\t%u\tcontig\tno\tna\n",refname,(j-gapsize+1),(j),numagpcomps,gapsize); // save previous gap info
						fprintf(contigsfile,">contig%u\n",numcontigs);
					}
					numagpcomps++; // new contig
				}
				if(outputfasta) fputc(c,contigsfile);
				contigsize++;
				numsnps++; // one more SNP (insertion)
				if(outputvcf) PrintVCFPosition(j,'I',refchar,c,(unsigned int)mapqsum[cid]);
				if(outputace){
					contigchars[extcontigsize]=c;
					contigrefcs[extcontigsize]='-';
					contigcodes[extcontigsize]='I';
					contigquals[extcontigsize]=(unsigned char)mapqsum[cid];
					extcontigsize++;
					if(extcontigsize==maxextcontigsize){ // alloc more space for the arrays if needed
						maxextcontigsize+=1024;
						contigchars=(unsigned char *)realloc(contigchars,(maxextcontigsize+1)*sizeof(unsigned char));
						contigrefcs=(unsigned char *)realloc(contigrefcs,(maxextcontigsize+1)*sizeof(unsigned char));
						contigcodes=(unsigned char *)realloc(contigcodes,(maxextcontigsize+1)*sizeof(unsigned char));
						contigquals=(unsigned char *)realloc(contigquals,(maxextcontigsize+1)*sizeof(unsigned char));
					}
				}
				m--; // next column with insertions
			} // end of processing insertions
			if(outputace){
				while(m>0){ // store information about the longest insertion size at this pos among all reads in the form of padding operations
					contigchars[extcontigsize]='*';
					contigrefcs[extcontigsize]='*';
					contigcodes[extcontigsize]='P';
					contigquals[extcontigsize]=(unsigned char)0;
					extcontigsize++;
					if(extcontigsize==maxextcontigsize){
						maxextcontigsize+=1024;
						contigchars=(unsigned char *)realloc(contigchars,(maxextcontigsize+1)*sizeof(unsigned char));
						contigrefcs=(unsigned char *)realloc(contigrefcs,(maxextcontigsize+1)*sizeof(unsigned char));
						contigcodes=(unsigned char *)realloc(contigcodes,(maxextcontigsize+1)*sizeof(unsigned char));
						contigquals=(unsigned char *)realloc(contigquals,(maxextcontigsize+1)*sizeof(unsigned char));
					}
					m--;
				}
			}
			if(j==refsize) break; // if we are over the last ref pos (to check for insertions in the last pos) break here
			for(i=0;i<5;i++){ // reset char counters
				counts[i]=0;
				mapqsum[i]=0;
				fwdcount[i]=0;
				revcount[i]=0;
			}
			for(i=0;i<numrows;i++){ // count chars on all rows of this column
				k=rowpos[i];
				temprow=rows[i];
				c=(temprow->extReadChars)[k];
				/*
				if(c=='A') cid=0;
				else if(c=='C') cid=1;
				else if(c=='G') cid=2;
				else if(c=='T') cid=3;
				else if(c=='-') cid=4;
				*/
				cid=(unsigned int)idFromChar[(int)c];
				counts[cid]++;
				mapqsum[cid]+=(temprow->mapQual);
				if((temprow->strand)=='+') fwdcount[cid]++;
				else revcount[cid]++;
			}
			cid=0;
			n=counts[0]; // current highest count
			for(i=1;i<5;i++){ // get the char that appears in the most rows
				if(counts[i]>n){
					n=counts[i]; // count of the most frequent char
					cid=i; // index of the most frequent char
				}
			}
			if(n==0) continue; // if all existent reads ended with insertions, they were processed above and we have reached the end of the reads but they were not deleted, and as no chars remain, this count is zero
			//c=refchars[j];
			/*
			if(c=='A') k=0; // index of the char at the reference in this column
			else if(c=='C') k=1;
			else if(c=='G') k=2;
			else if(c=='T') k=3;
			else k=4; // 'N' char (but not a gap symbol)
			*/
			k=(unsigned int)idFromChar[(int)refchar]; // index of the char at the reference in this column
			if( (k!=cid) && (k!=4) && (counts[k]==counts[cid]) ){ // when k=4 at the ref, it is a 'N' char, not a gap symbol
				cid=k; // if the ref char has the same count as the most freq char, prefer the ref char
			}
			/*
			if(cid==0) c='A'; // get most freq char letter
			else if(cid==1) c='C';
			else if(cid==2) c='G';
			else if(cid==3) c='T';
			else if(cid==4) c='-';
			*/
			c=charFromId[cid]; // get most freq char letter
			mapqsum[cid]=(mapqsum[cid]/counts[cid]); // average MAPQ sum over all rows with this char
			//if(c!=refchars[j]){
			if(c!=refchar){ // possible SNP or deletion
				if( counts[cid]>=mincoverage && mapqsum[cid]>=minmapq && ( (!bothstrands) || (fwdcount[cid]>0 && revcount[cid]>0) ) ){
					numsnps++; // one more SNP
				} else { // not enough conditions to call a SNP
					//c=refchars[j];
					c=refchar; // back to original ref char
					cid=k;
				}
			} // end of SNP/del confirmation
			if(c!='-'){ // deletions do not increase contig size
				if(outputfasta) fputc(c,consensusfile);
				if(contigsize==0){ // if starting a new contig, print contig name
					numcontigs++; // update contig count
					contigstartpos=j; // 1st position of the contig in the reference
					if(outputfasta){
						if(gapsize>0) fprintf(agpfile,"%s\t%u\t%u\t%u\tN\t%u\tcontig\tno\tna\n",refname,(j-gapsize+1),(j),numagpcomps,gapsize); // save previous gap info
						fprintf(contigsfile,">contig%u\n",numcontigs);
					}
					numagpcomps++; // new contig
				}
				if(outputfasta) fputc(c,contigsfile);
				contigsize++;
			}
			if(outputvcf){
				if(c==refchar) PrintVCFPosition(j,'M',refchar,c,(unsigned int)mapqsum[cid]); // match
				else if(c=='-') PrintVCFPosition(j,'D',refchar,c,(unsigned int)mapqsum[cid]); // deletion
				else PrintVCFPosition(j,'X',refchar,c,(unsigned int)mapqsum[cid]); // mismatch
			}
			if(outputace){
				contigchars[extcontigsize]=c;
				contigrefcs[extcontigsize]=refchar;
				//if(c==refchars[j])
				if(c==refchar) contigcodes[extcontigsize]='M'; // match
				else if(c=='-') contigcodes[extcontigsize]='D'; // deletion
				else contigcodes[extcontigsize]='X'; // mismatch
				contigquals[extcontigsize]=(unsigned char)mapqsum[cid];
				extcontigsize++;
				if(extcontigsize==maxextcontigsize){ // alloc more space for the arrays if needed
					maxextcontigsize+=1024;
					contigchars=(unsigned char *)realloc(contigchars,(maxextcontigsize+1)*sizeof(unsigned char));
					contigrefcs=(unsigned char *)realloc(contigrefcs,(maxextcontigsize+1)*sizeof(unsigned char));
					contigcodes=(unsigned char *)realloc(contigcodes,(maxextcontigsize+1)*sizeof(unsigned char));
					contigquals=(unsigned char *)realloc(contigquals,(maxextcontigsize+1)*sizeof(unsigned char));
				}
			}
			for(i=0;i<numrows;i++) rowpos[i]++; // go to next pos inside each row's arrays
		} // end of processing all positions of the reference
		j=refsize; // limit j to refsize because we could have advanced to (refsize+1)
		if(contigsize>0){ // save the last started contig
			consensussize+=contigsize;
			if(contigsize<mincontigsize) mincontigsize=contigsize;
			if(contigsize>maxcontigsize) maxcontigsize=contigsize;
			contigendpos=(j-1);
			if(outputfasta){
				fputc('\n',contigsfile);
				fprintf(agpfile,"%s\t%u\t%u\t%u\tD\tcontig%u\t1\t%u\t+\n",refname,(contigstartpos+1),(contigendpos+1),numagpcomps,numcontigs,contigsize);
			}
			if(outputace) PrintACEContig(numcontigs);
			gapsize=0;
		}
		if(gapsize>0){ // if consensus ended with a gap, save it
			if(outputfasta) fprintf(agpfile,"%s\t%u\t%u\t%u\tN\t%u\tcontig\tyes\n",refname,(j-gapsize+1),(j),numagpcomps,gapsize);
		}
		if(outputace) UpdateACEHeader(numcontigs);
		if(outputfasta) fputc('\n',consensusfile);
		printf(" (");PrintUnsignedNumber(consensussize);printf(" bp in ");PrintUnsignedNumber(numcontigs);printf(" contigs ; ");
		PrintUnsignedNumber(numunmappedpos);printf(" bp unmapped ; ");PrintUnsignedNumber(numsnps);printf(" SNPs) OK\n");
		printf("> Contigs sizes: max = ");PrintUnsignedNumber(maxcontigsize);
		printf(" bp ; avg = ");PrintUnsignedNumber((numcontigs!=0)?(consensussize/numcontigs):(0));
		printf(" bp ; min = ");PrintUnsignedNumber(mincontigsize);printf(" bp\n");
		fflush(stdout);
		free(rowpos);
		for(i=0;i<maxrows;i++) FreeMappedRead(rows[i]);
		free(rows);
		if(onefileperref) goto _close_output_files;
_free_and_next_ref: // free variables that depend on size of current ref
		free(refname);
		//free(refchars);
		free(readsIdsByPos);
		for(i=1;i<=numids;i++) free(readsStartFilePos[i]);
		free(readsStartFilePos);
		free(readsStartCountsInPos);
		free(normalizedrefname);
		if(outputcounts || outputstrandcounts){
			free(prevreadname);
			free(prevreadhitsids);
			for(i=1;i<=numids;i++) free(numreadhitsbyid[i]);
			free(numreadhitsbyid);
		}
		if(outputace){
			free(contigchars);
			free(contigrefcs);
			free(contigcodes);
			free(contigquals);
		}
		if(outputvcf){
			free(altchars);
		}
		if(targetrefname!=NULL) break; // if we only wanted a specific ref, we already did it, so stop here
	} // end of processing all refs
	if(targetrefname!=NULL && r==numrefs){ // if we only wanted a specific ref, but we checked all names but didn't find it, report error
		printf("> ERROR: Reference \"%s\" not found\n",targetrefname);
		exit(-1);
	}
	if(!onefileperref){ // one global file for all refs
_close_output_files:
		if(outputfasta){
			printf("> Saving consensus to <%s> ... ",consensusfilename);
			fflush(stdout);
			fclose(consensusfile);
			free(consensusfilename);
			printf("OK\n");
			printf("> Saving contigs to <%s> ... ",contigsfilename);
			fflush(stdout);
			fclose(contigsfile);
			free(contigsfilename);
			printf("OK\n");
			printf("> Saving contigs information to <%s> ... ",agpfilename);
			fflush(stdout);
			fclose(agpfile);
			free(agpfilename);
			printf("OK\n");
		}
		if(outputace){
			printf("> Saving ACE file to <%s> ... ",acefilename);
			fflush(stdout);
			fclose(acefile);
			free(acefilename);
			printf("OK\n");
		}
		if(outputvcf){
			printf("> Saving VCF file to <%s> ... ",vcffilename);
			fflush(stdout);
			fclose(vcffile);
			free(vcffilename);
			printf("OK\n");
		}
		if(outputgff){
			printf("> Saving annotations GFF file to <%s> ... ",genesgfffilename);
			fflush(stdout);
			fclose(genesgfffile);
			free(genesgfffilename);
			printf("OK\n");
			printf("> Saving annotations CSV file to <%s> ... ",genestxtfilename);
			fflush(stdout);
			fclose(genestxtfile);
			free(genestxtfilename);
			printf("OK\n");
			//printf("  FORMAT: <gene_id>,<number_of_reads>,<coverage>,<gene_consensus_sequence>\n");
			printf("  FORMAT: <gene_name>,<start_position>,<end_position>,<number_of_reads>,<coverage>\n");
			printf("> Saving SNPs CSV file to <%s> ... ",genessnpsfilename);
			fflush(stdout);
			fclose(genessnpsfile);
			free(genessnpsfilename);
			printf("OK\n");
			printf("  FORMAT: <gene_name>,<position>,<char_in_reference>,<char_in_consensus>\n");
		}
		if(outputcounts){
			printf("> Saving counts to <%s> ... ",countsfilename);
			fflush(stdout);
			fclose(countsfile);
			free(countsfilename);
			printf("OK\n");
			if(!outputstrandcounts) printf("  FORMAT: <position>,<number_of_reads>,<normalized_number_of_reads>\n");
		}
		if(outputstrandcounts){
			printf("> Saving forward strand counts to <%s> ... ",strandcountsfilename[0]);
			fflush(stdout);
			fclose(strandcountsfile[0]);
			free(strandcountsfilename[0]);
			printf("OK\n");
			printf("> Saving reverse strand counts to <%s> ... ",strandcountsfilename[1]);
			fflush(stdout);
			fclose(strandcountsfile[1]);
			free(strandcountsfilename[1]);
			printf("OK\n");
			printf("  FORMAT: <position>,<number_of_reads>,<normalized_number_of_reads>\n");
		}
		if(onefileperref) goto _free_and_next_ref; // we could have got here through the _close_output_files tag
	}
	printf("> Done!\n");
	fclose(reffile);
	fclose(samfile);
	free(refsizes);
	free(refnames);
	free(refsfpos);
	if(targetrefname!=NULL) free(targetrefname);
	#ifdef PAUSEATEXIT
	getchar();
	#endif
	exit(0);
}
