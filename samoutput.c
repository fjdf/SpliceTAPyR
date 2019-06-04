#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "tools.h"
#include "samoutput.h"
#include "dynprog.h"

#ifdef _MSC_VER
	#pragma warning(disable:4996)
#endif

// TODO: sort splice signals by frequency before outputting them (create struct with signal string and count, create compare function and use "qsort")
// TODO: Use binary search to find correct SAM reference name inside array of references names (use "bsearch" function)
// TODO: Output Paired-End / Mate-Pair statistics too
void ReportSAMStatistics(char *reffilename, char *samfilename, unsigned int realnumreads, int savebadreads){
	FILE *reffile, *samfile, *outfile, *maxerrorsfile;
	long long int filesize, progresscounter, progressstep, progressmark;
	char *refname, *reference, **refsnames, **references, *samrefname;
	char *read, *readname, *prevreadname, *outputfilename, *maxerrorsfilename, *stringptr;
	char c, cigarcode;
	int i, readpos, cigarcount, numcigarops, newread;
	int readsize, readnamesize, maxreadsize, maxreadnamesize;
	unsigned int n, numrefs, refsize, totalrefssize, *refssizes;
	fpos_t filepos, cigarstringpos;
	char *cigaralign, *refalign, *readalign;
	unsigned int alignsize, maxalignsize;
	int prevbestnumerrors, numerrors, maxerrors, isclippedread, numsoftclips, numhardclips, hascigarerrors, numsamerrors, numclippedchars;
	unsigned int totalnumreads, totalnumuniquereads, totalnummappedreads, totalnumuniquemappedreads;
	unsigned int numunmappedreads, numbadreads, numclippedreads, nummissingcigars, numcigarerrors;
	unsigned long long int totalnumchars, totalnumerrors;
	double errorpct, maxerrorpct, mappedreadspct, avgreadsize, avgerrorsperread, avgbppererror, coverage;
	double limitpct = 50.0;
	unsigned char splicesignalid;
	int lut_char2id[8];
	char lut_id2char[6];
	int issplicedread;
	unsigned int numsplicedreads, avgintronsize, minintronsize, maxintronsize;
	unsigned long long int totalnumsplicingevents, totalnumintronschars, *splicesignalscount;
	double splicesignalpct, avgsplicingperread;
	char splicesignal[5];
	splicesignal[0] = '\0'; splicesignal[1] = splicesignal[0]; // to prevent unused variable warnings
	int nret = 0; nret = (int)nret;
	printf("> Processing reference file <%s> ... ",reffilename);
	fflush(stdout);
	if((reffile=fopen(reffilename,"r"))==NULL){
		printf("\n> ERROR: File not found\n");
		exit(-1);
	}
	c=fgetc(reffile);
	if(c!='>'){
		printf("\n> ERROR: Not a valid FASTA file\n");
		exit(-1);
	}
	while(c!='\n' && c!=EOF) c=fgetc(reffile); // skip reference label
	if(c==EOF){
		printf("\n> ERROR: Not a valid FASTA file\n");
		exit(-1);
	}
	rewind(reffile);
	refssizes=NULL;
	numrefs=0;
	refsize=0;
	totalrefssize=0;
	while((c=fgetc(reffile))!=EOF){ // get reference chars
		if(c=='>'){
			if(numrefs!=0) refssizes[(numrefs-1)]=refsize; // save size of previous ref
			numrefs++; // number of reference sequences in file
			refssizes=(unsigned int *)realloc(refssizes,numrefs*sizeof(unsigned int));
			refsize=0;
			while(c!='\n' && c!=EOF) c=fgetc(reffile); // skip description in case of multi-fasta file
		} else if( (c>=65 && c<=90) || (c>=97 && c<=122) ){ // count only alphabet chars
			refsize++;
			totalrefssize++;
		}
	}
	refssizes[(numrefs-1)]=refsize; // save size of last ref
	if(totalrefssize==0){
		printf("\n> ERROR: Not a valid FASTA file\n");
		exit(-1);
	}
	printf("(%u sequence%s) ",numrefs,((numrefs!=1)?"s":""));
	fflush(stdout);
	rewind(reffile);
	references=(char **)calloc(numrefs,sizeof(char *));
	refsnames=(char **)calloc(numrefs,sizeof(char *));
	reference=NULL;
	refname=NULL;
	numrefs=0; // id of current ref
	refsize=0; // position of current char
	n=0; // number of invalid chars
	while((c=fgetc(reffile))!=EOF){ // get reference chars
		if(c=='>'){
			if(numrefs!=0) reference[refsize]='\0'; // finish chars string of previous ref
			refsize=refssizes[numrefs]; // load variables for current ref
			references[numrefs]=(char *)malloc((refsize+1)*sizeof(char));
			reference=references[numrefs];
			refname=(char *)malloc((255+1)*sizeof(char));
			i=0; // size of ref name
			while(c!='\n' && c!=EOF){ // get ref description
				c=fgetc(reffile);
				if(i==255) continue;
				if( c>='A' && c<='Z' ) c=(char)(c+32); // convert to lowercase
				if( (c>='a' && c<='z') || (c>='0' && c<='9') ) refname[i++]=c; // get only lowercase letters and numbers for ref name
			}
			refname[i]='\0';
			refname=(char *)realloc(refname,(i+1)*sizeof(char));
			refsnames[numrefs]=refname;
			numrefs++;
			refsize=0;
			continue;
		}
		if( c>='a' && c<='z' ) c=(char)(c-32); // convert to uppercase
		if( c>='A' && c<='Z' ){
			if( c!='A' && c!='C' && c!='G' && c!='T' && c!='N' ){ // alphabet letter but not ACGTN
				c='N';
				n++;
			}
			reference[refsize++]=c;
		}
	}
	reference[refsize]='\0';
	fclose(reffile);
	printf("(%u%s chars) ",refsize,((numrefs!=1)?" total":""));
	printf("OK\n");
	if(n!=0){
		printf("> WARNING: %u invalid bases (non ACGTN letters) were found in the reference\n",n);
		printf("           If the read mapping tool discarded them, this could lead to wrong alignment statistics\n");
	}
	printf("> Processing mapped reads file <%s> ",samfilename);
	fflush(stdout);
	if((samfile=fopen(samfilename,"r"))==NULL){
		printf("\n> ERROR: File not found\n");
		exit(-1);
	}
	fseek(samfile,0L,SEEK_END);
	filesize=(long long int)ftell(samfile);
	fseek(samfile,0L,SEEK_SET);
	outputfilename=NULL;
	outfile=NULL;
	maxerrorsfilename=NULL;
	maxerrorsfile=NULL;
	if(savebadreads){
		outputfilename=AppendToBasename(samfilename,"-badalignments.txt");
		maxerrorsfilename=AppendToBasename(samfilename,"-badreads.fasta");
		outfile=fopen(outputfilename,"w");
		maxerrorsfile=fopen(maxerrorsfilename,"w");
		if( outfile==NULL || maxerrorsfile==NULL ){
			printf("\n> ERROR: Unable to write output file\n");
			exit(-1);
		}
	}
	for(i=0;i<8;i++) lut_char2id[i]=0; // initialize lookup table to convert letter char to nucleotide id
	lut_char2id[(int)('A' & 0x07)]=0; // 'A' = 65 = ...001b & 0x07 = [1]
	lut_char2id[(int)('C' & 0x07)]=1; // 'C' = 67 = ...011b & 0x07 = [3]
	lut_char2id[(int)('G' & 0x07)]=2; // 'G' = 71 = ...111b & 0x07 = [7]
	lut_char2id[(int)('T' & 0x07)]=3; // 'T' = 84 = ...100b & 0x07 = [4]
	//lut_char2id[(int)('N' & 0x07)] = 4; // 'N' = 78 = ...110b & 0x07 = [6]
	lut_id2char[0]='A'; // initialize lookup table to convert nucleotide id to letter char
	lut_id2char[1]='C';
	lut_id2char[2]='G';
	lut_id2char[3]='T';
	lut_id2char[4]='N';
	lut_id2char[5]='\0';
	splicesignalscount=(long long unsigned int *)calloc(256,sizeof(long long unsigned int)); // array with 4x4x4x4=256 positions for all 4 letters splice signals
	splicesignal[0]='\0';
	splicesignal[1]='\0';
	splicesignal[2]='\0';
	splicesignal[3]='\0';
	splicesignal[4]='\0';
	maxreadsize=1;
	maxreadnamesize=1;
	read=(char *)malloc((maxreadsize+1)*sizeof(char));
	readname=(char *)malloc((maxreadnamesize+1)*sizeof(char));
	prevreadname=(char *)malloc((maxreadnamesize+1)*sizeof(char));
	strcpy(prevreadname,"\n");
	samrefname=(char *)malloc((255+1)*sizeof(char));
	maxalignsize=1;
	cigaralign=(char *)malloc((maxalignsize+1)*sizeof(char));
	refalign=(char *)malloc((maxalignsize+1)*sizeof(char));
	readalign=(char *)malloc((maxalignsize+1)*sizeof(char));
	c=fgetc(samfile);
	while(c=='@'){ // skip header lines
		while(c!='\n' && c!=EOF) c=fgetc(samfile);
		c=fgetc(samfile);
	}
	ungetc(c,samfile);
	maxerrors=0;
	maxerrorpct=0.0;
	numsoftclips=0;
	numhardclips=0;
	nummissingcigars=0;
	numcigarerrors=0;
	numsamerrors=0;
	numbadreads=0;
	numclippedreads=0;
	numunmappedreads=0;
	prevbestnumerrors=0;
	totalnumreads=0;
	totalnummappedreads=0;
	totalnumuniquereads=0;
	totalnumuniquemappedreads=0;
	totalnumchars=0;
	totalnumerrors=0;
	totalnumsplicingevents=0;
	totalnumintronschars=0;
	numsplicedreads=0;
	minintronsize=UINT_MAX;
	maxintronsize=0;
	progressstep=(filesize/10);
	progressmark=progressstep;
	progresscounter=0;
	while(c!=EOF){ // loop for all reads inside SAM file
		progresscounter=(long long int)ftell(samfile);
		if(progresscounter>=progressmark){ // print progress dots
			printf(".");
			fflush(stdout);
			progressmark+=progressstep;
		}
		fgetpos(samfile,&filepos);
		readnamesize=0;
		while((c=fgetc(samfile))!='\t' && c!=EOF) readnamesize++; // get read name size
		if(c==EOF) break;
		if(readnamesize>maxreadnamesize){
			maxreadnamesize=readnamesize;
			readname=(char *)realloc(readname,(maxreadnamesize+1)*sizeof(char));
			prevreadname=(char *)realloc(prevreadname,(maxreadnamesize+1)*sizeof(char));
		}
		fsetpos(samfile,&filepos);
		i=0;
		while((c=fgetc(samfile))!='\t' && c!=EOF) readname[i++]=c; // get read name
		readname[i]='\0';
		totalnumreads++; // one more read
		newread=1; // if it is a different read from the previous one
		if(strcmp(prevreadname,readname)==0) newread=0; // or the same read
		if(newread) totalnumuniquereads++; // one more (distinct) read
		for(i=0;i<1;i++){ // skip flag and advance to 2nd tab (reference label)
			c=fgetc(samfile);
			while(c!='\t' && c!=EOF) c=fgetc(samfile);
		}
		c=fgetc(samfile); // get first letter of reference
		if(c=='*' || c==EOF){ // if the read was not mapped, advance to the next one
			while(c!='\n' && c!=EOF) c=fgetc(samfile);
			numunmappedreads++;
			continue;
		}
		if(numrefs!=1){ // if there are multiple references, load the correct one
			ungetc(c,samfile);
			i=0;
			while(c!='\t' && c!=EOF){ // get ref name in SAM
				c=fgetc(samfile);
				if(i==255) continue;
				if( c>='A' && c<='Z' ) c=(char)(c+32); // convert to lowercase
				if( (c>='a' && c<='z') || (c>='0' && c<='9') ) samrefname[i++]=c; // get only lowercase letters and numbers for ref name
			}
			samrefname[i]='\0';
			for(n=0;n<numrefs;n++){ // search through all the refs names for the correct one
				refname=refsnames[n];
				i=0;
				while((samrefname[i]==refname[i]) && (samrefname[i]!='\0')) i++;
				if(samrefname[i]=='\0' && refname[i]=='\0') break; // match found
			}
			if(n==numrefs){ // ref name was not found
				printf("\n> ERROR: Reference named \"%s\" from read #%u \"%s\" was not found in the references file\n",samrefname,totalnumreads,readname);
				printf("         Please ensure that this SAM file was created using the same references file given\n");
				exit(-1);
			}
			reference=references[n]; // load correct reference variables
			refsize=refssizes[n];
		} else for(i=0;i<1;i++){ // otherwise, advance one more tab (position)
			c=fgetc(samfile);
			while(c!='\t' && c!=EOF) c=fgetc(samfile);
		}
		readpos=0;
		nret = fscanf(samfile,"%d",&readpos); // get read mapped position
		if(readpos>0) readpos--; // fix 1-based location
		for(i=0;i<2;i++){ // advance 2 more tabs (CIGAR string)
			c=fgetc(samfile);
			while(c!='\t' && c!=EOF) c=fgetc(samfile);
		}
		c=fgetc(samfile); // get first letter of the CIGAR string
		if(c=='*' || c==EOF){ // if the CIGAR string is not available, advance to the next read
			while(c!='\n' && c!=EOF) c=fgetc(samfile);
			nummissingcigars++;
			continue;
		}
		ungetc(c,samfile);
		fgetpos(samfile,&cigarstringpos); // save location of CIGAR string
		for(i=0;i<4;i++){ // advance 4 more tabs (sequence)
			c=fgetc(samfile);
			c=fgetc(samfile); // consume 2 chars to prevent two tabs in a row
			while(c!='\t' && c!=EOF) c=fgetc(samfile);
		}
		fgetpos(samfile,&filepos); // save location of read chars
		readsize=0;
		while((c=fgetc(samfile))!='\t' && c!=EOF) readsize++; // get size of the read
		if(c==EOF) break;
		if(readsize>maxreadsize){
			maxreadsize=readsize;
			read=(char *)realloc(read,(maxreadsize+1)*sizeof(char));
		}
		fsetpos(samfile,&filepos); // restore location to beginning of read chars
		i=0;
		while((c=fgetc(samfile))!='\t' && c!=EOF){ // get read chars
			if(c>=97 && c<=122) read[i++]=(char)(c-32); // convert to uppercase if needed
			else read[i++]=c;
		}
		read[i]='\0';
		fsetpos(samfile,&cigarstringpos); // restore location to beginning of CIGAR string
		isclippedread=0;
		numclippedchars=0;
		alignsize=0;
		numcigarops=0;
		while(fscanf(samfile,"%d%c",&cigarcount,&cigarcode)>0){ // get size of alignment strings
			numcigarops++;
			if(cigarcode=='P') continue;
			if(cigarcode=='S' || cigarcode=='H'){ // if it is a clipping operation
				if(cigarcount<=readsize) numclippedchars+=cigarcount; // used to later exclude from the read size the size of the clipped region
				if(numcigarops==1){ // if clipping is present at the beginning of the read, we need to make room for some positions before the starting position
					if(cigarcount>readpos) cigarcount=readpos; // prevent passing beyond the start of the reference
					numcigarops--;
				}
				if(cigarcode=='S') numsoftclips++;
				else if(cigarcode=='H') numhardclips++;
				isclippedread=1;
			}
			alignsize+=cigarcount;
		}
		if(alignsize>maxalignsize){ // realloc alignment strings if needed
			maxalignsize=alignsize;
			cigaralign=(char *)realloc(cigaralign,(maxalignsize+1)*sizeof(char));
			refalign=(char *)realloc(refalign,(maxalignsize+1)*sizeof(char));
			readalign=(char *)realloc(readalign,(maxalignsize+1)*sizeof(char));
		}
		fsetpos(samfile,&cigarstringpos); // restore location to beginning of CIGAR string
		i=0;
		numcigarops=0;
		while(fscanf(samfile,"%d%c",&cigarcount,&cigarcode)>0){ // create string with CIGAR code letters
			numcigarops++;
			if(cigarcode=='P') continue;
			if(numcigarops==1 && (cigarcode=='S' || cigarcode=='H')){ // if clipping is present at the beginning of the read
				if(cigarcount>readpos) cigarcount=readpos; // do not use more space behind than we have left
				numcigarops--;
			}
			while(cigarcount>0){
				cigaralign[i++]=cigarcode;
				cigarcount--;
			}
		}
		cigaralign[i]='\0';
		issplicedread=0;
		splicesignal[0]='\0';
		hascigarerrors=0;
		stringptr=(char *)(reference+readpos); // pointer to reference string location
		fsetpos(samfile,&cigarstringpos); // restore location to beginning of CIGAR string
		i=0; // current position inside reference alignment string
		n=(unsigned int)readpos; // current position inside reference
		numcigarops=0;
		while(fscanf(samfile,"%d%c",&cigarcount,&cigarcode)>0){ // create string with reference alignment
			numcigarops++;
			if(cigarcode=='M' || cigarcode=='D' || cigarcode=='N' || cigarcode=='=' || cigarcode=='X'){
				n+=cigarcount;
				if(n>refsize){ // check if it does not go beyond the end of the reference
					hascigarerrors=1;
					numcigarerrors++;
					cigarcount-=(n-refsize);
					n=refsize;
				}
				if(cigarcode=='N' && n!=refsize){ // get splice junction statistics
					issplicedread=1;
					totalnumsplicingevents++;
					totalnumintronschars+=cigarcount;
					if((unsigned int)cigarcount<minintronsize) minintronsize=cigarcount;
					if((unsigned int)cigarcount>maxintronsize) maxintronsize=cigarcount;
					c=reference[(n-cigarcount)]; // 1st position of the intron
					splicesignal[0]=c;
					splicesignalid=(unsigned char)lut_char2id[(int)(c & 0x07)];
					splicesignalid<<=2;
					c=reference[(n-cigarcount+1)]; // 2nd position of the intron
					splicesignal[1]=c;
					splicesignalid|=(unsigned char)lut_char2id[(int)(c & 0x07)];
					splicesignalid<<=2;
					c=reference[(n-2)]; // penultimate position of the intron
					splicesignal[2]=c;
					splicesignalid|=(unsigned char)lut_char2id[(int)(c & 0x07)];
					splicesignalid<<=2;
					c=reference[(n-1)]; // last position of the intron
					splicesignal[3]=c;
					splicesignalid|=(unsigned char)lut_char2id[(int)(c & 0x07)];
					splicesignalscount[(int)splicesignalid]++; // increase count of the 4 letters splice signal in the 4x4x4x4=256 positions array
				}
				while(cigarcount>0){
					refalign[i++]=(*stringptr);
					stringptr++; // advance in reference
					cigarcount--;
				}
			} else if(cigarcode=='I'){ // insertion into reference
				while(cigarcount>0){
					refalign[i++]='-';
					cigarcount--;
				}
			} else if(cigarcode=='S' || cigarcode=='H'){ // soft clipping or hard clipping
				if(numcigarops==1){ // if it occurs at the beginning of the read
					if(cigarcount>readpos) cigarcount=readpos; // do not use more positions than the number of reference positions behind
					stringptr=(char *)(stringptr-cigarcount); // go back in the reference
					n-=cigarcount;
					numcigarops--;
				}
				n+=cigarcount;
				if(n>refsize){ // check if it does not go beyond the end of the reference
					hascigarerrors=1;
					numcigarerrors++;
					cigarcount-=(n-refsize);
					n=refsize;
				}
				while(cigarcount>0){ // reference chars will be aligned against fake read chars, but they have to be there
					refalign[i++]=(*stringptr);
					stringptr++; // advance in reference
					cigarcount--;
				}
			} else if(cigarcode=='P'){ // padding
				cigarcount=0; // print nothing
			}
		} // end of processing CIGAR string
		refalign[i]='\0';
		if(issplicedread==1) numsplicedreads++;
		stringptr=(char *)(read); // pointer to beginning of read string
		fsetpos(samfile,&cigarstringpos); // restore location to beginning of CIGAR string
		i=0; // current position inside read alignment string
		n=0; // current position inside read
		numcigarops=0;
		while(fscanf(samfile,"%d%c",&cigarcount,&cigarcode)>0){ // create string with read alignment
			numcigarops++;
			if(cigarcode=='M' || cigarcode=='I' || cigarcode=='=' || cigarcode=='X'){ // match or mismatch or insertion into reference
				n+=cigarcount;
				if((int)n>readsize){ // check if it does not go beyond the end of the read
					hascigarerrors=1;
					numcigarerrors++;
					cigarcount-=(n-readsize);
					n=readsize;
				}
				while(cigarcount>0){
					readalign[i++]=(*stringptr);
					stringptr++; // advance in read
					cigarcount--;
				}
			} else if(cigarcode=='D'){ // deletion from reference
				while(cigarcount>0){
					readalign[i++]='-'; // print gap char
					cigarcount--;
				}
			} else if(cigarcode=='N'){ // intron
				while(cigarcount>0){
					readalign[i++]='_'; // print empty space
					cigarcount--;
				}
			} else if(cigarcode=='S'){ // soft clipping
				if(numcigarops==1){ // if it occurs at the beginning of the read
					if(cigarcount>readpos){
						stringptr=(char *)(stringptr+(cigarcount-readpos)); // skip chars that would go beyond the start of the reference
						n+=(cigarcount-readpos);
						cigarcount=readpos; // do not use more positions than the number of reference positions behind
					}
					numcigarops--;
				}
				n+=cigarcount;
				if((int)n>readsize){ // check if it does not go beyond the end of the read
					hascigarerrors=1;
					numcigarerrors++;
					cigarcount-=(n-readsize);
					n=readsize;
				}
				while(cigarcount>0){
					readalign[i++]=(char)((*stringptr)+32); // print char in lowercase
					stringptr++; // advance in read
					cigarcount--;
				}
			} else if(cigarcode=='H'){ // hard clipping
				if(numcigarops==1){ // if it occurs at the beginning of the read
					if(cigarcount>readpos) cigarcount=readpos; // do not use more positions than the number of reference positions behind
					numcigarops--;
				}
				while(cigarcount>0){
					readalign[i++]=' '; // print empty space
					cigarcount--;
				}
			} else if(cigarcode=='P'){ // padding
				cigarcount=0; // print nothing
			}
		} // end of processing CIGAR string
		if((int)n!=readsize){ // check if we got to the end of the read string
			hascigarerrors=1;
			numcigarerrors++;
			while((int)n<readsize){ // print the rest of the read chars
				if(cigaralign[i]=='S') readalign[i++]=(char)((*stringptr)+32); // print char in lowercase
				else readalign[i++]=(*stringptr);
				stringptr++; // advance in read
				n++;
			}
			readalign[i]='\0';
		}
		while(i<=(int)alignsize) readalign[i++]='\0';
		numerrors=0;
		for(n=0;n<alignsize;n++) if(refalign[n]!=readalign[n]){ // get number of errors
			c=cigaralign[n];
			if(c=='N' || c=='S' || c=='H' || c=='P') continue; // these operations do not count as errors
			if(c=='M') cigaralign[n]='X'; // change symbol in case of mismatch
			numerrors++;
		}
		if(numerrors>maxerrors) maxerrors=numerrors;
		errorpct=(((double)numerrors)*100.0)/((double)(readsize-numclippedchars));
		if(errorpct>maxerrorpct) maxerrorpct=errorpct;
		if(savebadreads){
			numsamerrors=numerrors; // get number of errors stored in SAM's optional tag, if available
			for(i=0;i<5;i++){ // advance 5 more tabs (optional fields)
				c=fgetc(samfile);
				c=fgetc(samfile); // consume 2 chars to prevent two tabs in a row
				while(c!='\t' && c!='\n' && c!=EOF) c=fgetc(samfile);
			}
			if(c!='\n' && c!=EOF){ // find tag "NM:i:%d"
				while(1){
					c=fgetc(samfile);
					if(c=='\n' || c==EOF) break;
					if(c!='N') continue;
					c=fgetc(samfile);
					if(c!='M') continue;
					c=fgetc(samfile);
					if(c!=':') continue;
					c=fgetc(samfile);
					if(c!='i') continue;
					c=fgetc(samfile);
					if(c!=':') continue;
					nret = fscanf(samfile,"%d",&numsamerrors);
					break;
				}
			}
			if((errorpct>limitpct) || (numsamerrors!=numerrors) || hascigarerrors){
				// || strcmp(splicesignal,"GTAG")!=0
				fprintf(outfile,">%s\n%s\n%s\n%s\n",readname,cigaralign,refalign,readalign);
				if(errorpct>limitpct) fprintf(outfile,"%d errors (%.0lf%% > %.0lf%%) at %d\n",numerrors,errorpct,limitpct,readpos);
				if(numsamerrors!=numerrors) fprintf(outfile,"%d errors found but %d reported in SAM\n",numerrors,numsamerrors);
				if(isclippedread) fprintf(outfile,"Has clipping\n");
				if(hascigarerrors) fprintf(outfile, "Has CIGAR errors\n");
				//if(strcmp(splicesignal,"GTAG")!=0) fprintf(outfile, "Splice signal: %s\n",splicesignal);
				fprintf(maxerrorsfile,">%s\n%s\n",readname,read);
				numbadreads++;
			}
		}
		if(newread){ // only count chars and errors one time per read
			if(!isclippedread){ // clipped reads to not count for the statistics
				totalnumuniquemappedreads++; // one more (distinct) mapped read
				totalnumchars+=readsize;
				totalnumerrors+=numerrors;
				prevbestnumerrors=numerrors;
			} else {
				numclippedreads++;
				prevbestnumerrors=0;
			}
		} else { // if it is just another hit for the same read
			if(!isclippedread){
				if(numerrors<prevbestnumerrors){ // if this new hit has a new better number of errors
					totalnumerrors-=prevbestnumerrors;
					totalnumerrors+=numerrors;
					prevbestnumerrors=numerrors; // only keep the number of errors of the best hit
				}
			}
		}
		totalnummappedreads++; // one more mapped read
		strcpy(prevreadname,readname);
		while(c!='\n' && c!=EOF) c=fgetc(samfile); // advance to next read
	} // end of loop for all reads inside SAM file
	fclose(samfile);
	printf(" (%u reads) OK\n",totalnumreads);
	if(numunmappedreads==0){
		printf("> WARNING: No un-mapped read records were found.\n");
		if (realnumreads<totalnumuniquereads) {
			printf("           If the mapping tool does not save the un-mapped reads in the SAM file, the mapping percentages will be incorrect.\n");
			printf("           Please provide the real total number of reads as an extra parameter, for correct statistics.\n");
		}
	}
	if(numcigarerrors!=0){
		printf("> WARNING: %u CIGAR records have inconsistencies related to the size of the read or the reference.\n",numcigarerrors);
		printf("           This prevents the correct calculation of the read error statistics.\n");
	}
	if(nummissingcigars!=0){
		printf("> WARNING: %u read records have the CIGAR string information missing.\n",nummissingcigars);
		printf("           This prevents the correct calculation of the read error statistics.\n");
	}
	if(numclippedreads!=0){
		printf("> WARNING: %u reads were detected with ",numclippedreads);
		if(numsoftclips>0 && numhardclips>0) printf("%d soft + hard",(numsoftclips+numhardclips));
		else if(numsoftclips>0) printf("%d soft",numsoftclips);
		else if(numhardclips>0) printf("%d hard",numhardclips);
		printf(" clipping operations\n");
		printf("           These reads are not considered as mapped and therefore will not be counted for the statistics.\n");
	}
	if(realnumreads>totalnumuniquereads) totalnumuniquereads=realnumreads; // set the real number of reads if defined
	totalnumuniquemappedreads-=numclippedreads;
	mappedreadspct=(((double)totalnumuniquemappedreads)*100.0)/((double)totalnumuniquereads);
	avgreadsize=((double)totalnumchars)/((double)totalnumuniquemappedreads);
	avgerrorsperread=((double)totalnumerrors)/((double)totalnumuniquemappedreads);
	if(totalnumerrors!=0) avgbppererror=((double)totalnumchars)/((double)totalnumerrors);
	else avgbppererror=0;
	coverage=((double)totalnumchars)/((double)totalrefssize);
	printf("> Statistics\n");
	printf("  Mapping    : %.2lf%% mapped distinct reads (%u of %u) ; %.2lf%% unmapped\n",mappedreadspct,totalnumuniquemappedreads,totalnumuniquereads,(100.0-mappedreadspct));
	printf("  Read size  : avg = %.0lf bp ; max = %d bp\n",avgreadsize,maxreadsize);
	printf("  Errors     : max = %d bp ; maxpct = %.2lf%% ; avg = %.2lf errors/read (%.0lf bp/error)\n",maxerrors,maxerrorpct,avgerrorsperread,avgbppererror);
	printf("  Coverage   : %.2lfx\n",coverage);
	if(totalnumsplicingevents!=0){
		avgintronsize=(unsigned int)(totalnumintronschars/totalnumsplicingevents);
		avgsplicingperread=((double)totalnumsplicingevents)/((double)numsplicedreads);
		printf("  RNA stats  : %u spliced reads ; %llu splicing events (%.2lf splices/read)\n",numsplicedreads,totalnumsplicingevents,avgsplicingperread);
		printf("  Intron size: avg = %u bp ; min = %u bp ; max = %u bp\n",avgintronsize,minintronsize,maxintronsize);
		for(i=0;i<256;i++){
			if(splicesignalscount[i]==0) continue;
			printf("\t\t[%c%c-%c%c]",lut_id2char[(i>>6)&(0x03)],lut_id2char[(i>>4)&(0x03)],lut_id2char[(i>>2)&(0x03)],lut_id2char[(i)&(0x03)]);
			splicesignalpct=(((double)splicesignalscount[i])*100.0)/((double)totalnumsplicingevents);
			printf(": %llu (%.2lf%%)\n",splicesignalscount[i],splicesignalpct);
		}
	}
	if(savebadreads){
		if(numbadreads==0){
			printf("> No reads found with any type of inconsistency\n");
		} else {
			printf("> Saving %u reads with some type of inconsistency to <%s> ... ",numbadreads,maxerrorsfilename);
			fflush(stdout);
			fclose(maxerrorsfile);
			free(maxerrorsfilename);
			printf("OK\n");
			printf("> Saving bad read alignments to <%s> ... ",outputfilename);
			fflush(stdout);
			fclose(outfile);
			free(outputfilename);
			printf("OK\n");
		}
	}
	printf("> Done!\n");
	free(read);
	free(readname);
	free(prevreadname);
	free(cigaralign);
	free(refalign);
	free(readalign);
	free(splicesignalscount);
	free(samrefname);
	for(n=0;n<numrefs;n++){
		free(references[n]);
		free(refsnames[n]);
	}
	free(references);
	free(refsnames);
	free(refssizes);
	#ifdef PAUSEATEXIT
	getchar();
	#endif
	exit(0);
}

void EvaluateSAMPositions(char *samfilename, int savewrongreads){
	FILE *samfile, *outfile;
	char c, *readname, *prevreadname, *tempptr, *outputfilename;
	int i, readpos, relax, correctreadpos[2];
	char strand, refchar, newread, correct;
	int readsize, maxreadnamesize;
	unsigned int numreads, numcorrectreads, numwrongreads, nummappedreads, numunmappedreads;
	int numreadhits, maxreadhits, totalreadhits;
	fpos_t readnamefpos, readcharsfpos;
	int nret = 0; nret = (int)nret;
	strand = '+'; strand = (char)strand;
	printf("> Processing SAM file <%s> ... ",samfilename);
	fflush(stdout);
	if((samfile=fopen(samfilename,"r"))==NULL){
		printf("\n> ERROR: File not found\n");
		exit(0);
	}
	outputfilename=NULL;
	outfile=NULL;
	if(savewrongreads){
		outputfilename=AppendToBasename(samfilename,"-incorrect.fasta");
		outfile=fopen(outputfilename,"w");
		if(outfile==NULL){
			printf("\n> ERROR: Unable to write output file\n");
			exit(0);
		}
	}
	maxreadnamesize=255;
	readname=(char *)calloc((maxreadnamesize+1),sizeof(char));
	prevreadname=(char *)calloc((maxreadnamesize+1),sizeof(char));
	c=fgetc(samfile);
	while(c=='@'){ // skip header lines
		while(c!='\n' && c!=EOF) c=fgetc(samfile);
		c=fgetc(samfile);
	}
	ungetc(c,samfile);
	relax=32; // allow this distance from the correct position
	numreads=0;
	nummappedreads=0;
	numunmappedreads=0;
	numcorrectreads=0;
	numwrongreads=0;
	numreadhits=0;
	maxreadhits=0;
	totalreadhits=0;
	correct=0; // set status for the first read
	while(c!=EOF){ // loop for all reads inside SAM file
		newread=0; // assume this is not a new read (it will be checked next)
		fgetpos(samfile,&readnamefpos); // save location of read name
		i=0;
		while((c=fgetc(samfile))!='\t' && c!=EOF){ // get read name
			newread |= (c!=prevreadname[i]); // if at least one char differs from the prev read name, it's a new read
			if(i<maxreadnamesize) readname[i++]=c;
		}
		readname[i]='\0';
		if(newread){ // if it is a new read
			numreads++;
			if(numreads>1){ // if it's not the first read, save results from the previous read
				if(correct) numcorrectreads++;
				else {
					numwrongreads++;
					if(savewrongreads){ // save unmapped or incorrectly mapped reads if needed
						fprintf(outfile,">%s\n",prevreadname);
						fsetpos(samfile,&readcharsfpos); // restore location of chars of prev read
						while((c=fgetc(samfile))!='\t' && c!=EOF) fputc(c,outfile); // write chars of prev read
						fputc('\n',outfile);
					}
				}
				if(numreadhits>maxreadhits) maxreadhits=numreadhits;
				totalreadhits+=numreadhits;
			}
			correct=0; // reset status of the current read
			numreadhits=1;
			tempptr=prevreadname;
			prevreadname=readname; // switch pointers so that the next read name can be compared against the current one
			readname=tempptr;
			fsetpos(samfile,&readnamefpos); // restore location of read name
			while(1){ // advance to the 1st underscore followed by a number
				c=fgetc(samfile);
				while(c!='_' && c!=EOF) c=fgetc(samfile);
				c=fgetc(samfile);
				if(c==EOF || (c>='0' && c<='9')){
					ungetc(c,samfile);
					break;
				}
			}
			nret = fscanf(samfile,"%d",&(correctreadpos[0])); // get correct (left) position
			for(i=0;i<1;i++){ // advance to the next underscore
				c=fgetc(samfile);
				while(c!='_' && c!=EOF) c=fgetc(samfile);
			}
			nret = fscanf(samfile,"%d",&(correctreadpos[1])); // get correct (right) position
			while(c!='\t' && c!=EOF) c=fgetc(samfile); // advance to next tab (flag)
		} else { // not a new read
			numreadhits++; // another hit for the same read
		}
		nret = fscanf(samfile,"%d",&i); // get flag
		if( i & 16 ) strand='-'; // get strand
		else strand='+';
		c=fgetc(samfile); // get '\t'
		refchar=fgetc(samfile); // get first letter of reference
		c=fgetc(samfile);
		while(c!='\t' && c!=EOF) c=fgetc(samfile); // advance to next tab (position)
		readpos=0;
		nret = fscanf(samfile,"%d",&readpos); // get read mapped position
		for(i=0;i<6;i++){ // advance 6 more tabs (sequence)
			c=fgetc(samfile);
			while(c!='\t' && c!=EOF) c=fgetc(samfile);
		}
		fgetpos(samfile,&readcharsfpos); // save location of read chars
		readsize=0;
		while((c=fgetc(samfile))!='\t' && c!=EOF) readsize++; // get size of the read
		if(refchar!='*' && readpos!=0){ // check if the read was mapped
			if(newread) nummappedreads++;
			if( abs(correctreadpos[0]-readpos)<=relax || abs(correctreadpos[1]-(readpos+readsize-1))<=relax )
				correct=1; // check if the read was mapped in the correct pos
		} else { // the read was not mapped
			if(newread) numunmappedreads++;
		}
		while(c!='\n' && c!=EOF) c=fgetc(samfile); // advance to next read
	} // end of loop for all reads inside SAM file
	if(numreads>1){ // save results for the last read
		if(correct) numcorrectreads++;
		else {
			numwrongreads++;
			if(savewrongreads){ // save unmapped or incorrectly mapped reads if needed
				fprintf(outfile,">%s\n",prevreadname);
				fsetpos(samfile,&readcharsfpos); // restore location of chars of prev read
				while((c=fgetc(samfile))!='\t' && c!=EOF) fputc(c,samfile); // write chars of prev read
				fputc('\n',samfile);
			}
		}
		if(numreadhits>maxreadhits) maxreadhits=numreadhits;
		totalreadhits+=numreadhits;
	}
	fclose(samfile);
	printf("OK\n");
	printf(":: %u total reads (%u mapped, %u unmapped)\n",numreads,nummappedreads,numunmappedreads);
	printf(":: %u correctly mapped reads (%.2lf%%)\n",numcorrectreads,(100.0*(double)numcorrectreads)/((double)numreads));
	printf(":: %u incorrectly mapped or unmapped reads (%.2lf%%)\n",numwrongreads,(100.0*(double)numwrongreads)/((double)numreads));
	printf(":: %.2lf avg hits/read (%u max hits)\n",((double)totalreadhits)/((double)nummappedreads),maxreadhits);
	if(savewrongreads){
		printf("> Saving wrong reads to <%s> ... ",outputfilename);
		fflush(stdout);
		fclose(outfile);
		free(outputfilename);
		printf("OK\n");
	}
	free(readname);
	free(prevreadname);
	printf("> Done!\n");
	#ifdef PAUSEATEXIT
	getchar();
	#endif
	exit(0);
}


/*
void SortTest(unsigned int *array, unsigned int arraysize){
	unsigned int mask, equalPositions;
	unsigned int top, bottom, savedTop, savedBottom[31];
	unsigned int depth;
	unsigned int swaptemp;
	unsigned int i;
	if( arraysize < 2 ) return; // already sorted
	mask = ~(0U); // bit 1 in all columns means all columns have 0s and 1s
	equalPositions = 0U; // with bit 1 if at least one row differs in that column, 0 otherwise
	i = 0;
	while( ( i <= (arraysize-2) ) && ( equalPositions != mask ) )
		equalPositions |= ( array[i] ^ array[++i] );
	equalPositions = ~(equalPositions); // with bit 1 if all rows are equal in that column, 0 otherwise
	mask = ( 1U << 31 ); // leftmost bit
	while( mask & equalPositions ) mask >>= 1; // find the first column from the left that differs
	if( mask == 0 ) return; // already sorted
	depth = 0;
	top = 0;
	bottom = (arraysize-1);
	while(1){ // loop for columns
		savedTop = top;
		savedBottom[depth] = bottom;
		while(1){ // loop for rows
			while( ( top < bottom ) && !( array[top] & mask ) ) top++; // descend in array while there are 0s at the beginning
			while( ( bottom > top ) && ( array[bottom] & mask ) ) bottom--; // ascend in array while there are 1s at the end
			if( top >= bottom ) break; // pointer reached each other, else there is a 1 where there should be a 0 and vice-versa
			swaptemp = array[top]; // swap top and bottom positions
			array[top] = array[bottom];
			array[bottom] = swaptemp;
			top++; // advance to next position bellow
			bottom--; // advance to next position above
		}
		if( bottom > top ){ // if the loop above ended with a swap
			top--; // top was pointing at the first 1 bellow the 0s
			bottom++; // bottom was pointing at the first 0 above the 1s
		}
		if( ( top - savedTop ) >= 2 ){ // at least 2 items with 0s to sort
			bottom = top; // position of last 0
			top = savedTop; // position of first 0
		} else if( ( savedBottom[depth] - bottom ) >= 2 ) { // at least 2 items with 1s to sort
			top = bottom; // position of first 1
			bottom = savedBottom[depth]; // position of last 1
		} else { // none or only one item to sort

		}
		depth++;
		mask >>= 1;
		while( mask & equalPositions ) mask >>= 1; // skip columns to the left with all equal rows
		if( mask == 0 ) mask = 1;
	}
}
*/

typedef struct _ReadPositionAndFPos {
	int pos;
	fpos_t filepos;
} ReadPositionAndFPos;

int CompareReadPositionAndFPos(const void *a, const void *b){
	return ( ((ReadPositionAndFPos *)a)->pos - ((ReadPositionAndFPos *)b)->pos );
}

void ConvertSAMToCSV(char *samfilename){
	FILE *samfile, *txtfile;
	int readpos, flag, i;
	unsigned int numreads;
	char *txtfilename, *readname, *refname, strand, c;
	unsigned int readsarraysize;
	ReadPositionAndFPos *readsarray;
	fpos_t readlinestart;
	int nret = 0; nret = (int)nret;
	printf("> Opening SAM file <%s> ... ",samfilename);
	fflush(stdout);
	if((samfile=fopen(samfilename,"r"))==NULL){
		printf("\n> ERROR: SAM file not found\n");
		exit(0);
	}
	numreads=0;
	readsarraysize=1024;
	readsarray=(ReadPositionAndFPos *)malloc(readsarraysize*sizeof(ReadPositionAndFPos));
	c=fgetc(samfile);
	while(c=='@'){ // header
		while(c!='\n') c=fgetc(samfile);
		c=fgetc(samfile);
	}
	ungetc(c,samfile);
	while((c=fgetc(samfile))!=EOF){
		ungetc(c,samfile);
		fgetpos(samfile,&readlinestart);
		for(i=0;i<2;i++){ // skip read name and flag
			c=fgetc(samfile);
			while(c!='\t' && c!=EOF) c=fgetc(samfile);
		}
		c=fgetc(samfile);
		if(c=='*'){ // no ref
			while(c!='\n' && c!=EOF) c=fgetc(samfile);
			continue;
		}
		while(c!='\t' && c!=EOF) c=fgetc(samfile); // skip ref name
		readpos=0;
		nret = fscanf(samfile,"%d\t",&readpos); // read position
		if(readpos==0){ // no pos
			while(c!='\n' && c!=EOF) c=fgetc(samfile);
			continue;
		}
		if(numreads==readsarraysize){ // realloc array with more space
			readsarraysize+=1024;
			readsarray=(ReadPositionAndFPos *)realloc(readsarray,readsarraysize*sizeof(ReadPositionAndFPos));
		}
		readsarray[numreads].pos = readpos;
		readsarray[numreads].filepos = readlinestart;
		numreads++;
		while(c!='\n' && c!=EOF) c=fgetc(samfile); // skip rest of line
	}
	readsarraysize=numreads;
	printf("(%u reads) OK\n",numreads);
	printf("> Sorting positions ... ");
	fflush(stdout);
	qsort(readsarray,readsarraysize,sizeof(ReadPositionAndFPos),CompareReadPositionAndFPos);
	printf("OK\n");
	txtfilename=AppendToBasename(samfilename,".txt");
	printf("> Creating CSV file <%s> ... ",txtfilename);
	fflush(stdout);
	if((txtfile=fopen(txtfilename,"w"))==NULL){
		printf("\n> ERROR: Can't write CSV file\n");
		exit(0);
	}
	refname=(char *)calloc((255+1),sizeof(char));
	readname=(char *)calloc((255+1),sizeof(char));
	for(numreads=0;numreads<readsarraysize;numreads++){
		readlinestart=(readsarray[numreads].filepos);
		fsetpos(samfile,&readlinestart);
		i=0;
		c=fgetc(samfile);
		while(c!='\t' && c!=EOF){ // read name
			if(i<255) readname[i++]=c;
			c=fgetc(samfile);
		}
		readname[i]='\0';
		nret = fscanf(samfile,"%d\t",&flag); // flag
		if( flag & 16 ) strand='-';
		else strand='+';
		i=0;
		c=fgetc(samfile);
		while(c!='\t' && c!=EOF){ // reference name
			if(i<255) refname[i++]=c;
			c=fgetc(samfile);
		}
		refname[i]='\0';
		readpos=0;
		nret = fscanf(samfile,"%d\t",&readpos); // read position
		fprintf(txtfile,"%d,%c,%s,%s,",readpos,strand,readname,refname); // write to file
		for(i=0;i<5;i++){ // skip MAPQ, CIGAR, RNEXT, PNEXT, TLEN
			c=fgetc(samfile);
			while(c!='\t' && c!=EOF) c=fgetc(samfile);
		}
		c=fgetc(samfile);
		while(c!='\t' && c!=EOF){ // write read sequence
			fputc(c,txtfile);
			c=fgetc(samfile);
		}
		fputc('\n',txtfile);
	}
	printf(" OK\n");
	printf("> FORMAT: <position>,<strand>,<read_name>,<reference_name>,<read_sequence>\n");
	fflush(stdout);
	fclose(samfile);
	fclose(txtfile);
	free(txtfilename);
	free(refname);
	free(readname);
	free(readsarray);
	printf("> Done!\n");
	exit(0);
}
