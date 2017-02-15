#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#ifdef _MSC_VER
	#pragma warning(disable:4996)
#endif

// Checks and parses program arguments sent using the command line parameters
// NOTE: if parse = 0, returns 1 if the argument exists and 0 otherwise
// NOTE: if parse = 1, returns its value if the argument exists and -1 otherwise
// NOTE: if parse = 2, returns its index in the arglist array if the argument exists and -1 otherwise
int ParseArgument(int numargs, char** arglist, char *optionchars, int parse){
	int i;
	char uc[2], lc[2];
	for(i=0;i<2;i++){ // only process the first two chars of the option name
		if(optionchars[i]=='\0'){ // no char
			uc[i] = '\0';
			lc[i] = '\0';
			continue;
		}
		if( ( ((int)optionchars[i]) >= 65 ) && ( ((int)optionchars[i]) <= 90 ) ){ // uppercase
			uc[i] = optionchars[i];
			lc[i] = (char)(optionchars[i]+32);
		} else { // lowercase
			uc[i] = (char)(optionchars[i]-32);
			lc[i] = optionchars[i];
		}
	}
	for( i = 1 ; i < numargs ; i++ ){
		if( arglist[i][0]=='-' &&
				( arglist[i][1]==lc[0] || arglist[i][1]==uc[0] ) &&
				( arglist[i][2]==lc[1] || arglist[i][2]==uc[1] ) ){
			if(parse){
				if(i==(numargs-1)) return (-1); // if it's the last argument, it has nothing in front
				if(parse==1) return atoi(arglist[(i+1)]); // parse argument in front
				else return (i+1);
			}
			return 1; // if we only want to check if it exits, report 1
		}
	}
	if(parse) return (-1); // we wanted a value, but the argument was not found
	return 0; // argument not present
}

// Parses program arguments that appear multiple times in the command line parameters
int *ParseMultiArgument(int numargs, char** arglist, char *optionchars, int *count){
	int i, n, *outarray;
	char uc[2], lc[2];
	for(i=0;i<2;i++){ // only process the first two chars of the option name
		if(optionchars[i]=='\0'){ // no char
			uc[i] = '\0';
			lc[i] = '\0';
			continue;
		}
		if( ( ((int)optionchars[i]) >= 65 ) && ( ((int)optionchars[i]) <= 90 ) ){ // uppercase
			uc[i] = optionchars[i];
			lc[i] = (char)(optionchars[i]+32);
		} else { // lowercase
			uc[i] = (char)(optionchars[i]-32);
			lc[i] = optionchars[i];
		}
	}
	n = 0; // number of argument hits found
	outarray = NULL;
	for( i = 1 ; i < (numargs-1) ; i++ ){ // (numargs-1) because last one cannot be an option char
		if( arglist[i][0]=='-' &&
				( arglist[i][1]==lc[0] || arglist[i][1]==uc[0] ) &&
				( arglist[i][2]==lc[1] || arglist[i][2]==uc[1] ) ){
			n++;
			if(outarray==NULL) outarray=(int *)malloc(n*sizeof(int));
			else outarray=(int *)realloc(outarray,n*sizeof(int));
			outarray[(n-1)]=atoi(arglist[(i+1)]); // parse argument in front
		}
	}
	(*count)=n;
	return outarray;
}

// joins the basename of a name of a file with another string
char* AppendToBasename(char *filename, char *extra){
	char *resultfilename;
	int i, n;
	n=(int)strlen(filename);
	for(i=(n-1);i>0;i--){
		if(filename[i]=='.') break;
	}
	if(i==0) i=n;
	n=(int)strlen(extra);
	resultfilename=(char *)calloc((i+n+1),sizeof(char));
	strncpy(resultfilename,filename,i);
	resultfilename[i]='\0';
	strcat(resultfilename,extra);
	return resultfilename;
}

// Returns a new string with all special characters removed or converted
// NOTE: mode=0 removes all non-alphanumeric chars
// NOTE: mode=1 converts all non-alphanumeric chars to underscores
// NOTE: mode=2 keeps all chars except spaces, tabs, newlines and non-printable chars
char *NormalizeSeqName(char *name, int mode){
	char c, *newname;
	int i, n;
	n=0;
	while(name[n]!='\0') n++;
	newname=(char *)malloc((n+1)*sizeof(char));
	i=0;
	n=0;
	while((c=name[n++])!='\0'){
		if( (c>='0' && c<='9') || (c>='A' && c<='Z') || (c>='a' && c<='z') ) newname[i++]=c; // alpha-numeric chars
		else if(mode!=0) { // convert all other chars to underscore
			if(mode==1){ if( i==0 || newname[(i-1)]!='_' ) newname[i++]='_'; } // prevent writing of two followed underscores
			else if( c>='!' && c<='~' ) newname[i++]=c; // non-blank chars
		}
	}
	newname[i]='\0';
	return newname;
}

// prints an integer value in a formated way
void PrintNumber( long long int number ){
	long long num, denom, quot;
	//div_t divres;
	if( number < 1000 ){
		printf( "%d" , (int)number );
		return;
	}
	denom = 1;
	num = number;
	while( num >= 1000 ){
		//divres = div( num , 1000 );
		//num = divres.quot;
		num /= 1000;
		denom *= 1000;
	}
	printf( "%d," , (int)num );
	num = ( number - num * denom );
	denom /= 1000;
	while( denom > 1 ){
		//divres = div( num , denom );
		//printf( "%03d." , divres.quot );
		//num = divres.rem;
		quot = ( num / denom );
		printf( "%03d," , (int)quot );
		num -= ( quot * denom );
		denom /= 1000;
	}
	printf( "%03d" , (int)num );
}

// prints an unsigned integer value in a formated way
void PrintUnsignedNumber( unsigned int number ){
	unsigned int num, denom, quot, rem;
	if( number < 1000 ){
		printf( "%u" , number );
		return;
	}
	denom = 1;
	num = number;
	while( num >= 1000 ){
		num /= 1000;
		denom *= 1000;
	}
	printf( "%u," , num );
	num = ( number - num * denom );
	denom /= 1000;
	while( denom > 1 ){
		quot = ( num / denom );
		rem = ( num - quot * denom );
		printf( "%03u," , quot );
		num = rem;
		denom /= 1000;
	}
	printf( "%03u" , num );
}

// prints a time value in a formated way
void PrintTime( double timeval ){
	int time, days, hours, minutes, seconds;
	div_t divres;
	time = (int) timeval;
	if( time == 0 ){
		printf( "0s " );
		return;
	}
	divres = div( time , 60 );
	seconds = divres.rem;
	time = divres.quot;
	divres = div( time , 60 );
	minutes = divres.rem;
	time = divres.quot;
	divres = div( time , 24 );
	hours = divres.rem;
	days = divres.quot;
	if( days != 0 ) printf( "%dD " , days );
	if( hours != 0 ) printf( "%dH " , hours );
	if( minutes != 0 ) printf( "%dm " , minutes );
	if( seconds != 0 ) printf( "%ds " , seconds );
}

// Shows and updates a progress bar (in the current line or in the line above the current one)
void PrintProgressBar(double percentage, int lineabove){
	static int prevbarpos = -1;
	//const char symbols[] = "\x20\xB0\xB1\xB2\xDB\x2D";
	const char symbols[] = "\x2D\x2D\x2D\x2D\x61\x2D";
	const int barsize = 50;
	const int numchars = 5;
	const int leftspaces = 2;
	int barpos, charpos, i;
	if(percentage<0.0) percentage=0.0;
	if(percentage>100.0) percentage=100.0;
	if(percentage==0.0 && prevbarpos==-1){ // initialize progress bar
		for(i=0;i<(leftspaces-1);i++) putchar(' ');
		printf("[");
		for(i=0;i<barsize;i++) putchar(symbols[5]); // fill with gap symbols
		printf("]\n");
		prevbarpos = 0;
		return;
	}
	barpos=(int)floor(percentage*((double)barsize)/100.0); // which pos
	charpos=(int)floor(percentage*(double)(barsize*numchars)/100.0)-(barpos*numchars); // which char (percentage mod numchars)
	if(barpos==barsize){ // even if we are at 100%, the last bar pos is at (barsize-1)
		barpos=(barsize-1);
		charpos=(numchars-1);
	}
	printf("\x1B(0"); // printf("\x1B[11m"); // alternate font
	if(lineabove) printf("\x1B[99D\x1B[K"); // move cursor to beginning of line and delete until the end of the line
	if(lineabove) printf("\x1B[1F"); // previous line (1A = line up)
	i=prevbarpos; // previous filled position in bar
	printf("\x1B[%dC",(leftspaces+i)); // move forward (#G = go to column #)
	while((i++)<barpos) putchar(symbols[4]); // fill incomplete spaces behind if needed
	putchar(symbols[charpos]);
	printf("\x1B[%dD",(leftspaces+i)); // move back (1G = go to column 1)
	//if(lineabove) printf("\x1B[1E"); // next line (1B = line down)
	printf("\x1B(B"); // printf("\x1B[10m"); // back to normal font
	prevbarpos=barpos; // save last bar position
	//if(lineabove) printf("\x1B[99D\x1B[K"); // move cursor to beginning of line and delete until the end of the line
}
/*
void PrintProgressBar(long int current, long int total){
	char *symbols = "\x20\xB0\xB1\xB2\xDB\x2D";
	int barsize = 20;
	int steppct = 5;
	//static int prevbarpos = 0;
	int percentage, barpos, charpos, i;
	if(current==0){ // initialize progress bar
		printf(" [");
		for(i=0;i<barsize;i++) putchar(symbols[5]);
		printf("]\n");
		return;
	}
	percentage=(int)(current*100/total);
	if(percentage<0) percentage=0;
	if(percentage>100) percentage=100;
	barpos=(percentage/steppct); // 0=[0-4], 1=[5-9], 2=[10-14], ... , 19=[95-99]
	charpos=(percentage-(barpos*steppct)); // percentage mod 5 (0-4)
	if(barpos==20){ // even if we are at 100%, the last one is at 19
		barpos=19;
		charpos=4;
	}
	printf("\x1B[11m"); // alternate font
	printf("\x1B[1F"); // previous line (1A = line up)
	//i=prevbarpos; // starting column
	i=barpos; // starting column
	printf("\x1B[%dC",(2+i)); // move forward (#G = go to column #)
	//while(i<barpos) putchar(symbols[4]); // fill incomplete spaces behind if needed
	putchar(symbols[charpos]);
	printf("\x1B[%dD",i); // move back (1G = go to column 1)
	printf("\x1B[1E"); // next line (1B = line down)
	printf("\x1B[10m"); // back to normal font
	//prevbarpos=barpos;
}
*/

// Cleans all invalid characters from a FASTA file, leaving only A, C, G and T characters; and joins Multi-FASTA files
void CleanFasta(char *fastafilename){
	FILE *fastafile,*cleanfile;
	char *cleanfilename;
	int i,n,charcount,invalidcharcount,sequencecount;
	char c;
	printf("> Opening Fasta file <%s> ... ",fastafilename);
	fflush(stdout);
	if((fastafile=fopen(fastafilename,"r"))==NULL){
		printf("\n> ERROR: Fasta file not found\n");
		exit(0);
	}
	c=fgetc(fastafile);
	if(c!='>'){
		printf("\n> ERROR: Invalid Fasta file\n");
		exit(0);
	}
	printf("OK\n");
	n=(int)strlen(fastafilename);
	for(i=(n-1);i>0;i--){
		if(fastafilename[i]=='.') break;
	}
	if(i==0) i=n;
	cleanfilename=(char *)calloc((n+6+1),sizeof(char));
	strncpy(cleanfilename,fastafilename,i);
	strcat(cleanfilename,"-clean");
	strcat(cleanfilename,(char *)(fastafilename+i));
	printf("> Creating Clean Fasta file <%s> ... ",cleanfilename);
	fflush(stdout);
	if((cleanfile=fopen(cleanfilename,"w"))==NULL){
		printf("\n> ERROR: Can't write Clean Fasta file\n");
		exit(0);
	}
	fprintf(cleanfile,">%s\n",fastafilename); // filename as label
	charcount=0;
	invalidcharcount=0;
	sequencecount=0;
	while(c!=EOF){
		if(c=='A' || c=='C' || c=='G' || c=='T'){
			fputc((int)c,cleanfile);
			charcount++;
		} else if(c=='a' || c=='c' || c=='g' || c=='t'){
			fputc((int)(c-32),cleanfile);
			charcount++;
		} else if(c=='>'){ // new sequence
			while(c!='\n' && c!=EOF) c=fgetc(fastafile); // skip description of sequence
			sequencecount++;
		} else if(c>32 && c<127) invalidcharcount++; // invalid alphanumeric character
		c=fgetc(fastafile);
	}
	fputc('\n',cleanfile);
	fclose(fastafile);
	fclose(cleanfile);
	free(cleanfilename);
	printf("(%d valid chars ; %d invalid chars ; %d sequences) OK\n",charcount,invalidcharcount,sequencecount);
	printf("> Done!\n");
	exit(0);
}

// Converts a FASTAQ file to the FASTA format
void ConvertFastaQ(char *fastaqfilename){
	FILE *fastaqfile,*fastafile;
	char *fastafilename;
	int i,n;
	unsigned int charcount,invalidcharcount,readscount;
	char c;
	printf("> Opening FastaQ file <%s> ... ",fastaqfilename);
	fflush(stdout);
	if((fastaqfile=fopen(fastaqfilename,"r"))==NULL){
		printf("\n> ERROR: FastaQ file not found\n");
		exit(0);
	}
	c=fgetc(fastaqfile);
	if(c!='@'){
		printf("\n> ERROR: Invalid FastaQ file\n");
		exit(0);
	}
	printf("OK\n");
	n=(int)strlen(fastaqfilename);
	for(i=(n-1);i>0;i--){
		if(fastaqfilename[i]=='.') break;
	}
	if(i==0) i=n;
	fastafilename=(char *)calloc((i+6+1),sizeof(char));
	strncpy(fastafilename,fastaqfilename,i);
	strcat(fastafilename,".fasta");
	printf("> Creating Fasta file <%s> ... ",fastafilename);
	fflush(stdout);
	if((fastafile=fopen(fastafilename,"w"))==NULL){
		printf("\n> ERROR: Can't write Fasta file\n");
		exit(0);
	}
	readscount=0;
	charcount=0;
	invalidcharcount=0;
	while(c!=EOF){ // new read starting with c='@'
		fputc('>',fastafile);
		c=fgetc(fastaqfile);
		while(c!='\n' && c!=EOF){ // line with read description
			fputc(c,fastafile);
			c=fgetc(fastaqfile);
		}
		fputc('\n',fastafile);
		c=fgetc(fastaqfile);
		n=0; // read size
		while(c!='+' && c!=EOF){ // line(s) with read sequence
			if( (c>=55 && c<=90) || (c>=97 && c<=122) ){ // uppercase or lowercase letters
				fputc(c,fastafile);
				charcount++;
				n++;
			}
			if(c!='\n' && c!='A' && c!='C' && c!='G' && c!='T' && c!='a' && c!='c' && c!='g' && c!='t') invalidcharcount++;
			c=fgetc(fastaqfile);
		}
		fputc('\n',fastafile);
		c=fgetc(fastaqfile);
		while(c!='\n' && c!=EOF) c=fgetc(fastaqfile); // skip line starting with '+'
		i=0;
		while(i<n && c!=EOF){ // skip line(s) with ascii quality values
			c=fgetc(fastaqfile);
			if(c>32 && c<127) i++; // printable ASCII chars
		}
		c=fgetc(fastaqfile);
		while(c!='@' && c!=EOF) c=fgetc(fastaqfile); // go to beginning of next read
		readscount++;
	}
	fclose(fastaqfile);
	fclose(fastafile);
	free(fastafilename);
	printf("(%u reads; avg %u bp/read; %u invalid chars) OK\n",readscount,(charcount/readscount),invalidcharcount);
	printf("> Done!\n");
	exit(0);
}

/*
// reads an array of 16bit big-endian integers from file and converts them to little-endian
void freada16be(uint16_t *dest, int count, FILE *file){
	int i;
	for(i=0;i<count;i++){
		fread(dest,sizeof(uint16_t),(size_t)1,file);
		(*dest) = ((*dest) << 8) | ((*dest) >> 8);
		dest++;
	}
}

// reads an array of 32bit big-endian integers from file and converts them to little-endian
void freada32be(uint32_t *dest, int count, FILE *file){
	int i;
	for(i=0;i<count;i++){
		fread(dest,sizeof(uint32_t),(size_t)1,file);
		(*dest) = ((*dest) << 24) | (((*dest) & 0x0000FF00) << 8) | (((*dest) & 0x00FF0000) >> 8) | ((*dest) >> 24);
		dest++;
	}
}

// reads an array of 64bit big-endian integers from file and converts them to little-endian
void freada64be(uint64_t *dest, int count, FILE *file){
	int i;
	unsigned char *destarray;
	for(i=0;i<count;i++){
		destarray=(unsigned char *)dest;
		destarray[7]=fgetc(file);
		destarray[6]=fgetc(file);
		destarray[5]=fgetc(file);
		destarray[4]=fgetc(file);
		destarray[3]=fgetc(file);
		destarray[2]=fgetc(file);
		destarray[1]=fgetc(file);
		destarray[0]=fgetc(file);
		dest++;
	}
}
*/

// reads a 16bit big-endian integer from a file and converts it to little-endian
void fread16be(uint16_t *dest, FILE *file){
	fread(dest,sizeof(uint16_t),(size_t)1,file);
	(*dest) = ((*dest) << 8) | ((*dest) >> 8);
}

// reads a 32bit big-endian integer from a file and converts it to little-endian
void fread32be(uint32_t *dest, FILE *file){
	fread(dest,sizeof(uint32_t),(size_t)1,file);
	(*dest) = ((*dest) << 24) | (((*dest) & 0x0000FF00) << 8) | (((*dest) & 0x00FF0000) >> 8) | ((*dest) >> 24);
}

// Converts a SFF file to the FASTQ format
void ConvertSFF(char *sfffilename){
	FILE *sfffile, *fastafile;
	char *fastafilename;
	unsigned int i,n,charscount,readscount;
	char c;
	/*
	// SFF Format: Common Header Section
	uint32_t magic_number;
	char version[4]; // char[4]
	uint64_t index_offset;
	uint32_t index_length;
	uint32_t number_of_reads;
	uint16_t header_length;
	uint16_t key_length;
	uint16_t number_of_flows_per_read;
	uint8_t flowgram_format_code;
	char *flow_chars; // char[number_of_flows_per_read]
	char *key_sequence; // char[key_length]
	uint8_t *eight_byte_padding; // uint8_t[*]
	// SFF Format: Read Header Section
	uint16_t read_header_length;
	uint16_t name_length;
	uint32_t number_of_bases;
	uint16_t clip_qual_left;
	uint16_t clip_qual_right;
	uint16_t clip_adapter_left;
	uint16_t clip_adapter_right;
	char *name; // char[name_length]
	//uint8_t *eight_byte_padding; // uint8_t[*]
	// SFF Format: Read Data Section
	uint16_t *flowgram_values; // uint16_t[number_of_flows_per_read]
	uint8_t *flow_index_per_base; // uint8_t[number_of_bases]
	char *bases; // char[number_of_bases]
	uint8_t *quality_scores; // uint8_t[number_of_bases]
	//uint8_t *eight_byte_padding; // uint8_t[*]
	unsigned int k;
	*/
	unsigned int number_of_reads;
	unsigned int number_of_bases;
	unsigned int padding_size;
	unsigned int start_base_pos;
	unsigned int end_base_pos;
	unsigned int max_name_length;
	unsigned int max_number_of_bases;
	uint32_t sff_magic_number;
	uint32_t sff_number_of_reads;
	uint16_t sff_key_length;
	uint16_t sff_number_of_flows;
	uint16_t sff_name_length;
	uint32_t sff_number_of_bases;
	uint16_t sff_clip_qual_left;
	uint16_t sff_clip_qual_right;
	uint16_t sff_clip_adapter_left;
	uint16_t sff_clip_adapter_right;
	char *sff_name;
	char *sff_bases;
	uint8_t *sff_quals;
	readscount=0;
	charscount=0;
	printf("> Opening SFF file <%s> ... ",sfffilename);
	fflush(stdout);
	if((sfffile=fopen(sfffilename,"rb"))==NULL){
		printf("\n> ERROR: SFF file not found\n");
		exit(0);
	}
	/*
	printf("\nCOMMON HEADER:\n");
	freada32be(&magic_number,1,sfffile);
	printf("  magic_number\t\t\t= %#X\n",magic_number);
	fread(version,sizeof(char),(size_t)4,sfffile);
	printf("  version\t\t\t= ");
	for(k=0;k<4;k++) printf("%hu",(unsigned short int)version[k]);
	printf("\n");
	freada64be(&index_offset,1,sfffile);
	printf("  index_offset\t\t\t= %lu\n",index_offset);
	freada32be(&index_length,1,sfffile);
	printf("  index_length\t\t\t= %u\n",index_length);
	freada32be(&number_of_reads,1,sfffile);
	printf("  number_of_reads\t\t= %u\n",number_of_reads);
	freada16be(&header_length,1,sfffile);
	printf("  header_length\t\t\t= %hu\n",header_length);
	freada16be(&key_length,1,sfffile);
	printf("  key_length\t\t\t= %hu\n",key_length);
	freada16be(&number_of_flows_per_read,1,sfffile);
	printf("  number_of_flows_per_read\t= %hu\n",number_of_flows_per_read);
	fread(&flowgram_format_code,sizeof(uint8_t),(size_t)1,sfffile);
	printf("  flowgram_format_code\t\t= %hu\n",flowgram_format_code);
	flow_chars=(char *)calloc(number_of_flows_per_read,sizeof(char));
	fread(flow_chars,sizeof(char),(size_t)number_of_flows_per_read,sfffile);
	printf("  flow_chars\t\t\t= %.*s\n",number_of_flows_per_read,flow_chars);
	key_sequence=(char *)calloc(key_length,sizeof(char));
	fread(key_sequence,sizeof(char),(size_t)key_length,sfffile);
	printf("  key_sequence\t\t\t= %.*s\n",key_length,key_sequence);
	eight_byte_padding=(uint8_t *)calloc(8,sizeof(uint8_t));
	padding_size = ( ( 8 - ( ( 31 + number_of_flows_per_read + key_length ) & 7 ) ) & 7 );
	printf("  eight_byte_padding\t\t= ");
	if(padding_size!=0){
		fread(eight_byte_padding,sizeof(uint8_t),(size_t)padding_size,sfffile);
		for(k=0;k<padding_size;k++) printf("%hu",eight_byte_padding[k]);
	}
	printf("\n");
	printf("READ HEADER:\n");
	freada16be(&read_header_length,1,sfffile);
	printf("  read_header_length\t\t= %hu\n",read_header_length);
	freada16be(&name_length,1,sfffile);
	printf("  name_length\t\t\t= %hu\n",name_length);
	freada32be(&number_of_bases,1,sfffile);
	printf("  number_of_bases\t\t= %u\n",number_of_bases);
	freada16be(&clip_qual_left,1,sfffile);
	printf("  clip_qual_left\t\t= %hu\n",clip_qual_left);
	freada16be(&clip_qual_right,1,sfffile);
	printf("  clip_qual_right\t\t= %hu\n",clip_qual_right);
	freada16be(&clip_adapter_left,1,sfffile);
	printf("  clip_adapter_left\t\t= %hu\n",clip_adapter_left);
	freada16be(&clip_adapter_right,1,sfffile);
	printf("  clip_adapter_right\t\t= %hu\n",clip_adapter_right);
	name=(char *)calloc(name_length,sizeof(char));
	fread(name,sizeof(char),(size_t)name_length,sfffile);
	printf("  name\t\t\t\t= %.*s\n",name_length,name);
	padding_size = ( ( 8 - ( ( 16 + name_length ) & 7 ) ) & 7 );
	printf("  eight_byte_padding\t\t= ");
	if(padding_size!=0){
		fread(eight_byte_padding,sizeof(uint8_t),(size_t)padding_size,sfffile);
		for(k=0;k<padding_size;k++) printf("%hu",eight_byte_padding[k]);
	}
	printf("\n");
	if( clip_qual_left == 0 ) clip_qual_left = 1; // lowest left position is 1 (1-based array)
	if( clip_adapter_left == 0 ) clip_adapter_left = 1;
	if( clip_qual_left > clip_adapter_left ) start_base_pos = clip_qual_left; // choose maximum value
	else start_base_pos = clip_adapter_left;
	start_base_pos--; // convert to 0-based array
	if( clip_qual_right == 0 ) clip_qual_right = number_of_bases; // highest right position is number_of_bases
	if( clip_adapter_right == 0 ) clip_adapter_right = number_of_bases;
	if( clip_qual_right < clip_adapter_right) end_base_pos = clip_qual_right; // choose minimum value
	else end_base_pos = clip_adapter_right;
	end_base_pos--; // convert to 0-based array
	printf("READ DATA:\n");
	printf("  flowgram_values\t\t= ");
	flowgram_values=(uint16_t *)calloc(number_of_flows_per_read,sizeof(uint16_t));
	freada16be(flowgram_values,number_of_flows_per_read,sfffile);
	for(k=0;k<number_of_flows_per_read;k++) printf("%03hu,",flowgram_values[k]);
	printf("\n");
	printf("  flow_index_per_base\t\t= ");
	flow_index_per_base=(uint8_t *)calloc(number_of_bases,sizeof(uint8_t));
	fread(flow_index_per_base,sizeof(uint8_t),(size_t)number_of_bases,sfffile);
	for(k=0;k<number_of_bases;k++) printf("%hu,",flow_index_per_base[k]);
	printf("\n");
	bases=(char *)calloc(number_of_bases,sizeof(char));
	fread(bases,sizeof(char),(size_t)number_of_bases,sfffile);
	printf("  bases\t\t\t\t= %.*s\n",number_of_bases,bases);
	printf("  quality_scores\t\t= ");
	quality_scores=(uint8_t *)calloc(number_of_bases,sizeof(uint8_t));
	fread(quality_scores,sizeof(uint8_t),(size_t)number_of_bases,sfffile);
	for(k=0;k<number_of_bases;k++) printf("%02hu,",quality_scores[k]);
	printf("\n");
	padding_size = ( ( 8 - ( ( 2*number_of_flows_per_read + 3*number_of_bases ) & 7 ) ) & 7 );
	printf("  eight_byte_padding\t\t= ");
	if(padding_size!=0){
		fread(eight_byte_padding,sizeof(uint8_t),(size_t)padding_size,sfffile);
		for(k=0;k<padding_size;k++) printf("%hu",eight_byte_padding[k]);
	}
	printf("\n");
	free(flow_chars);
	free(key_sequence);
	free(eight_byte_padding);
	free(name);
	free(flowgram_values);
	free(flow_index_per_base);
	free(bases);
	free(quality_scores);
	getchar();
	rewind(sfffile);
	*/
	fread32be(&sff_magic_number,sfffile); // magic_number[4]
	if(sff_magic_number!=0x2E736666){ // ".sff"
		printf("\n> ERROR: Invalid SFF file\n");
		exit(0);
	}
	max_name_length=16;
	sff_name=(char *)calloc((max_name_length+1),sizeof(char));
	max_number_of_bases=255;
	sff_bases=(char *)calloc((max_number_of_bases+1),sizeof(char));
	sff_quals=(uint8_t *)calloc((max_number_of_bases+1),sizeof(uint8_t));
	for(i=0;i<16;i++) c=fgetc(sfffile); // skip 16 bytes: version[4] + index_offset[8] + index_length[4]
	fread32be(&sff_number_of_reads,sfffile); // number_of_reads[4]
	for(i=0;i<2;i++) c=fgetc(sfffile); // skip 2 bytes: header_length[2]
	fread16be(&sff_key_length,sfffile); // key_length[2]
	fread16be(&sff_number_of_flows,sfffile); // number_of_flows[2]
	padding_size = ( ( 8 - ( ( 31 + sff_number_of_flows + sff_key_length ) & 7 ) ) & 7 );
	n=(1+sff_number_of_flows+sff_key_length+padding_size); // skip: flowgram_format_code[1] + flow_chars[number_of_flows] + key_sequence[key_length] + eight_byte_padding[padding_size]
	for(i=0;i<n;i++) c=fgetc(sfffile); // end of common header section
	printf("(%u reads) OK\n",sff_number_of_reads);
	n=(int)strlen(sfffilename);
	for(i=(n-1);i>0;i--){
		if(sfffilename[i]=='.') break;
	}
	if(i==0) i=n;
	fastafilename=(char *)calloc((i+6+1),sizeof(char));
	strncpy(fastafilename,sfffilename,i);
	strcat(fastafilename,".fastq");
	printf("> Creating FastQ file <%s> ... ",fastafilename);
	fflush(stdout);
	if((fastafile=fopen(fastafilename,"w"))==NULL){
		printf("\n> ERROR: Can't write FastQ file\n");
		exit(0);
	}
	number_of_reads=0;
	while( number_of_reads < sff_number_of_reads ){
		for(i=0;i<2;i++) c=fgetc(sfffile); // skip 2 bytes: read_header_length[2]
		fread16be(&sff_name_length,sfffile); // name_length[2]
		fread32be(&sff_number_of_bases,sfffile); // number_of_bases[4]
		fread16be(&sff_clip_qual_left,sfffile); // clip_qual_left[2]
		fread16be(&sff_clip_qual_right,sfffile); // clip_qual_right[2]
		fread16be(&sff_clip_adapter_left,sfffile); // clip_adapter_left[2]
		fread16be(&sff_clip_adapter_right,sfffile); // clip_adapter_right[2]
		if( sff_clip_qual_left > sff_clip_adapter_left ) start_base_pos = sff_clip_qual_left; // choose maximum value
		else start_base_pos = sff_clip_adapter_left;
		if( start_base_pos > 0 ) start_base_pos--; // convert to 0-based array
		if( sff_clip_qual_right == 0 ) sff_clip_qual_right = sff_number_of_bases; // highest right position is number_of_bases
		if( sff_clip_adapter_right == 0 ) sff_clip_adapter_right = sff_number_of_bases;
		if( sff_clip_qual_right < sff_clip_adapter_right) end_base_pos = sff_clip_qual_right; // choose minimum value
		else end_base_pos = sff_clip_adapter_right;
		end_base_pos--; // convert to 0-based array
		if( sff_name_length > max_name_length ){ // allocate more memory for array if needed
			max_name_length = sff_name_length;
			sff_name=(char *)realloc(sff_name,(max_name_length+1)*sizeof(char));
		}
		fread(sff_name,sizeof(char),(size_t)sff_name_length,sfffile); // name[name_length]
		sff_name[sff_name_length] = '\0'; // add terminator char to string
		padding_size = ( ( 8 - ( ( 16 + sff_name_length ) & 7 ) ) & 7 );
		n = ( padding_size + 2*sff_number_of_flows + sff_number_of_bases ); // skip: eight_byte_padding[padding_size] + flowgram_values[2][number_of_flows] + flow_index_per_base[number_of_bases]
		for(i=0;i<n;i++) c=fgetc(sfffile); // end of read header section and beginning of read data section
		number_of_bases = ( end_base_pos - start_base_pos + 1 );
		if( number_of_bases > max_number_of_bases ){ // allocate more memory for array if needed
			max_number_of_bases = number_of_bases;
			sff_bases=(char *)realloc(sff_bases,(max_number_of_bases+1)*sizeof(char));
			sff_quals=(uint8_t *)realloc(sff_quals,(max_number_of_bases+1)*sizeof(uint8_t));
		}
		n = ( start_base_pos ); // skip the left clipped bases
		for(i=0;i<n;i++) c=fgetc(sfffile);
		fread(sff_bases,sizeof(char),(size_t)number_of_bases,sfffile); // bases[number_of_bases]
		sff_bases[number_of_bases] = '\0'; // add terminator char to string
		n = ( sff_number_of_bases - end_base_pos - 1 ); // skip the right clipped bases
		for(i=0;i<n;i++) c=fgetc(sfffile);
		n = ( start_base_pos ); // skip the left clipped qualities
		for(i=0;i<n;i++) c=fgetc(sfffile);
		fread(sff_quals,sizeof(uint8_t),(size_t)number_of_bases,sfffile); // quality_scores[number_of_bases]
		for(i=0;i<number_of_bases;i++) sff_quals[i] += 33; // convert quality scores to Sanger format (Phred+33)
		sff_quals[number_of_bases] = '\0'; // add terminator char to string
		n = ( sff_number_of_bases - end_base_pos - 1 ); // skip the right clipped qualities
		for(i=0;i<n;i++) c=fgetc(sfffile);
		padding_size = ( ( 8 - ( ( 2*sff_number_of_flows + 3*sff_number_of_bases ) & 7 ) ) & 7 );
		n = ( padding_size ); // skip: eight_byte_padding[padding_size]
		for(i=0;i<n;i++) c=fgetc(sfffile);
		if(c==EOF) break;
		number_of_reads++;
		fprintf(fastafile,"@%s\n%s\n+\n%s\n",sff_name,sff_bases,sff_quals);
		readscount++;
		charscount+=number_of_bases;
	}
	fclose(sfffile);
	fclose(fastafile);
	free(fastafilename);
	printf("(%u reads; avg %u bp/read) OK\n",readscount,(charscount/readscount));
	printf("> Done!\n");
	exit(0);
}

// Calculates the coverage percentage using the genome size and a file with the positions of the aligned reads
void CalculateCoverage(char *positionsfilename, int genomesize){
	FILE *positionsfile;
	unsigned short int *poscount;
	int readcount,pos,length,numerrors,filledposcount,totalreadsize,maxreadsize;
	double coverage;
	char c,*readdesc;
	printf("> Opening Aligned Positions file <%s> ... ",positionsfilename);
	fflush(stdout);
	if((positionsfile=fopen(positionsfilename,"r"))==NULL){
		printf("\n> ERROR: Aligned Positions file not found\n");
		exit(0);
	}
	printf("OK\n");
	printf("> Allocating space for %d position counters ... ",genomesize);
	fflush(stdout);
	poscount=(unsigned short int *)calloc((genomesize+1),sizeof(unsigned short int));
	if(poscount==NULL){
		printf("\n> ERROR: Not enough memory\n");
		exit(0);
	}
	printf("OK\n");
	printf("> Marking positions ... ");
	fflush(stdout);
	readdesc=(char *)malloc(255*sizeof(char));
	readcount=0;
	totalreadsize=0;
	maxreadsize=0;
	c=fgetc(positionsfile);
	while(c=='>'){
		pos=0;
		while(c!='\n'){ // get read description
			readdesc[pos++]=c;
			c=fgetc(positionsfile);
		}
		readdesc[pos]='\0';
		c=fgetc(positionsfile);
		if(c=='>' || c==EOF) continue; // skip empty reads
		ungetc(c,positionsfile);
		fscanf(positionsfile,"%d %d %d\n",&pos,&length,&numerrors);
		if(length>maxreadsize) maxreadsize=length;
		if(pos<0) pos=(-pos);
		totalreadsize+=length;
		while(pos<=genomesize && length>0){
			poscount[pos]++;
			pos++;
			length--;
		}
		readcount++;
		c=fgetc(positionsfile); // get '>' from next read or EOF if at end of positions file
	}
	printf("(%d reads with a total of %d bp ; max read size = %d bp) OK\n",readcount,totalreadsize,maxreadsize);
	printf("> Counting positions ... ");
	fflush(stdout);
	filledposcount=0;
	for(pos=0;pos<=genomesize;pos++){
		if(poscount[pos]!=0) filledposcount++;
	}
	printf("(%d filled positions) OK\n",filledposcount);
	fflush(stdout);
	fclose(positionsfile);
	free(poscount);
	free(readdesc);
	coverage=((double)filledposcount)/((double)genomesize)*100.0;
	printf("> Coverage = %.2lf%%\n",coverage);
	printf("> Done!\n");
	exit(0);
}

// TODO: how to count mismatches?
// TODO: detection of hard clipping or soft clipping when calculating read length?
// Retrieves the positions of the aligned reads from a file in the SAM format (from MAQ/BWA output)
void ProcessSamFile(char *samfilename, int totalnumreads){
	FILE *samfile,*coveragefile;
	char *coveragefilename,*prevreaddesc,*readdesc,*auxreaddesc;
	int i,n,pos,length,sameread,numerrors;
	int numreads,numalignedreads;
	unsigned int flag;
	char c, strand;
	printf("> Opening SAM file <%s> ... ",samfilename);
	fflush(stdout);
	if((samfile=fopen(samfilename,"r"))==NULL){
		printf("\n> ERROR: SAM file not found\n");
		exit(0);
	}
	printf("OK\n");
	i=(int)strlen(samfilename);
	coveragefilename=(char *)calloc((i+8+1),sizeof(char));
	strcpy(coveragefilename,samfilename);
	strcat(coveragefilename,".cov.txt");
	printf("> Creating Coverage file <%s> ... ",coveragefilename);
	fflush(stdout);
	if((coveragefile=fopen(coveragefilename,"w"))==NULL){
		printf("\n> ERROR: Can't write Coverage file\n");
		exit(0);
	}
	prevreaddesc=(char *)calloc(255,sizeof(char));
	readdesc=(char *)calloc(255,sizeof(char));
	numreads=0;
	numalignedreads=0;
	c=fgetc(samfile);
	while(c=='@'){ // header
		while(c!='\n') c=fgetc(samfile);
		c=fgetc(samfile);
	}
	ungetc(c,samfile);
	while((c=fgetc(samfile))!=EOF){
		auxreaddesc=prevreaddesc;
		prevreaddesc=readdesc;
		readdesc=auxreaddesc;
		sameread=1;
		i=0;
		while(c!='\t'){ // read description
			readdesc[i]=c;
			sameread=(sameread & (readdesc[i]==prevreaddesc[i]));
			c=fgetc(samfile);
			i++;
		}
		numreads++;
		fputc('>',coveragefile);
		n=i;
		for(i=0;i<n;i++) fputc(readdesc[i],coveragefile);
		fputc('\n',coveragefile);
		fscanf(samfile,"%u\t%c",&flag,&c); // flag and reference
		if( flag & 16 ) strand='-';
		else strand='+';
		if(c!='*'){
			if(!sameread) numalignedreads++;
			while(c!='\t') c=fgetc(samfile); // rest of genome description
			fscanf(samfile,"%d\t",&pos); // position
			c=fgetc(samfile);
			while(c!='\t') c=fgetc(samfile); // mapping quality
			numerrors=0;
			c=fgetc(samfile);
			while(c!='\t'){ // CIGAR string
				ungetc(c,samfile);
				fscanf(samfile,"%d%c",&n,&c);
				if( c=='I' || c=='D' ) numerrors+=n;
				c=fgetc(samfile);
			}
			for(i=0;i<3;i++){ // mate name, mate pos, insert size
				c=fgetc(samfile);
				while(c!='\t') c=fgetc(samfile);
			}
			length=0;
			c=fgetc(samfile);
			while(c!='\t'){ // read sequence
				length++;
				c=fgetc(samfile);
			}
			fprintf(coveragefile,"%c%d %d %d\n",strand,pos,length,numerrors);
		}
		while(c!='\n') c=fgetc(samfile);
	}
	fclose(samfile);
	fclose(coveragefile);
	free(coveragefilename);
	free(prevreaddesc);
	free(readdesc);
	numreads=totalnumreads;
	printf("(%.2lf%% aligned reads - %d of %d) OK\n",(((double)numalignedreads)/((double)numreads)*100.0),numalignedreads,numreads);
	printf("> Done!\n");
	exit(0);
}

// Retrieves the positions of the aligned reads from a file in the 454 format (454ReadStatus.txt file from Newbler output)
void Process454File(char *rsfilename){
	FILE *rsfile,*coveragefile;
	char *coveragefilename;
	int i,startpos,endpos,size,noerrorpct,readpct,numerrors;
	int numreads,numalignedreads;
	char c,strand;
	printf("> Opening 454 file <%s> ... ",rsfilename);
	fflush(stdout);
	if((rsfile=fopen(rsfilename,"r"))==NULL){
		printf("\n> ERROR: 454 file not found\n");
		exit(0);
	}
	printf("OK\n");
	i=(int)strlen(rsfilename);
	coveragefilename=(char *)calloc((i+8+1),sizeof(char));
	strcpy(coveragefilename,rsfilename);
	strcat(coveragefilename,".cov.txt");
	printf("> Creating Coverage file <%s> ... ",coveragefilename);
	fflush(stdout);
	if((coveragefile=fopen(coveragefilename,"w"))==NULL){
		printf("\n> ERROR: Can't write Coverage file\n");
		exit(0);
	}
	c=fgetc(rsfile);
	if(c!='R'){ // "Read"
		printf("\n> ERROR: Invalid 454 file\n");
		exit(0);
	}
	while(c!='\n') c=fgetc(rsfile); // header line 1
	c=fgetc(rsfile);
	while(c!='\n') c=fgetc(rsfile); // header line 2
	numreads=0;
	numalignedreads=0;
	while((c=fgetc(rsfile))!=EOF){
		numreads++;
		fputc('>',coveragefile);
		while(c!='\t'){ // read description
			fputc(c,coveragefile);
			c=fgetc(rsfile);
		}
		fputc('\n',coveragefile);
		c=fgetc(rsfile);
		if(c!='F' && c!='P'){ // not "Full" and not "Partial"
			if(c=='R' || c=='C') numalignedreads++; // "Repeat" or "Chimeric"
			while(c!='\n') c=fgetc(rsfile);
			continue;
		}
		numalignedreads++;
		while(c!='\t') c=fgetc(rsfile); // mapping status
		fscanf(rsfile,"%d\t",&noerrorpct); // mapped accuracy percentage
		fscanf(rsfile,"%d\t",&readpct); // percentage of read mapped
		c=fgetc(rsfile);
		while(c!='\t') c=fgetc(rsfile); // reference description
		fscanf(rsfile,"%d\t%d\t%c",&startpos,&endpos,&strand); // start position, end position, strand
		size=(endpos-startpos+1);
		readpct=((size*readpct)/100);
		noerrorpct=((readpct*noerrorpct)/100);
		numerrors=(size-noerrorpct); // calculate number of errors
		fprintf(coveragefile,"%c%d %d %d\n",strand,startpos,size,numerrors);
		while(c!='\n') c=fgetc(rsfile);
	}
	fclose(rsfile);
	fclose(coveragefile);
	free(coveragefilename);
	printf("(%.2lf%% aligned reads - %d of %d) OK\n",(((double)numalignedreads)/((double)numreads)*100.0),numalignedreads,numreads);
	printf("> Done!\n");
	exit(0);
}

// Retrieves the positions of the aligned reads from a file in the MAP format (from SEGEMEHL output)
void ProcessMapFile(char *mapfilename, int totalnumreads){
	FILE *mapfile,*coveragefile;
	char *coveragefilename,*prevreaddesc,*readdesc,*auxreaddesc;
	int i,n,startpos,endpos,length,sameread,numerrors;
	int numalignedreads;
	char c, strand;
	printf("> Opening MAP file <%s> ... ",mapfilename);
	fflush(stdout);
	if((mapfile=fopen(mapfilename,"r"))==NULL){
		printf("\n> ERROR: MAP file not found\n");
		exit(0);
	}
	c=fgetc(mapfile);
	if(c!='#'){
		printf("\n> ERROR: Invalid MAP file\n");
		exit(0);
	}
	printf("OK\n");
	i=(int)strlen(mapfilename);
	coveragefilename=(char *)calloc((i+8+1),sizeof(char));
	strcpy(coveragefilename,mapfilename);
	strcat(coveragefilename,".cov.txt");
	printf("> Creating Coverage file <%s> ... ",coveragefilename);
	fflush(stdout);
	if((coveragefile=fopen(coveragefilename,"w"))==NULL){
		printf("\n> ERROR: Can't write Coverage file\n");
		exit(0);
	}
	prevreaddesc=(char *)calloc(255,sizeof(char));
	readdesc=(char *)calloc(255,sizeof(char));
	numalignedreads=0;
	while(c!='>') c=fgetc(mapfile); // header
	ungetc(c,mapfile);
	while((c=fgetc(mapfile))!=EOF){
		auxreaddesc=prevreaddesc;
		prevreaddesc=readdesc;
		readdesc=auxreaddesc;
		sameread=1;
		i=0;
		while(c!='\t'){ // read description
			c=fgetc(mapfile);
			readdesc[i]=c;
			sameread=(sameread & (readdesc[i]==prevreaddesc[i]));
			i++;
		}
		if(!sameread) numalignedreads++;
		fputc('>',coveragefile);
		n=i;
		for(i=0;i<(n-1);i++) fputc(readdesc[i],coveragefile);
		fputc('\n',coveragefile);
		fscanf(mapfile,"%d\t",&numerrors); // edit distance
		for(i=0;i<8;i++){ // other information
			c=fgetc(mapfile);
			while(c!='\t') c=fgetc(mapfile);
		}
		fscanf(mapfile,"%c\t%d\t%d",&strand,&startpos,&endpos);
		length=(endpos-startpos+1);
		fprintf(coveragefile,"%c%d %d %d\n",strand,startpos,length,numerrors);
		while(c!='\n') c=fgetc(mapfile);
	}
	fclose(mapfile);
	fclose(coveragefile);
	free(coveragefilename);
	free(prevreaddesc);
	free(readdesc);
	printf("(%.2lf%% aligned reads - %d of %d) OK\n",(((double)numalignedreads)/((double)totalnumreads)*100.0),numalignedreads,totalnumreads);
	printf("> Done!\n");
	exit(0);
}

// TODO: output the number of errors too
// Retrieves the positions of the aligned reads from a file in the SOAP output format
void ProcessSoapFile(char *soapfilename, int totalnumreads){
	FILE *soapfile,*coveragefile;
	char *coveragefilename;
	int i,pos,length,numtabs;
	int numalignedreads;
	char c,strand;
	printf("> Opening SOAP file <%s> ... ",soapfilename);
	fflush(stdout);
	if((soapfile=fopen(soapfilename,"r"))==NULL){
		printf("\n> ERROR: SOAP file not found\n");
		exit(0);
	}
	printf("OK\n");
	i=(int)strlen(soapfilename);
	coveragefilename=(char *)calloc((i+8+1),sizeof(char));
	strcpy(coveragefilename,soapfilename);
	strcat(coveragefilename,".cov.txt");
	printf("> Creating Coverage file <%s> ... ",coveragefilename);
	fflush(stdout);
	if((coveragefile=fopen(coveragefilename,"w"))==NULL){
		printf("\n> ERROR: Can't write Coverage file\n");
		exit(0);
	}
	numalignedreads=0;
	while((c=fgetc(soapfile))!=EOF){
		numtabs=0;
		fputc('>',coveragefile);
		while(c!='\t'){
			fputc(c,coveragefile);
			c=fgetc(soapfile);
		}
		fputc('\n',coveragefile);
		numtabs=1;
		while(numtabs!=5){
			c=fgetc(soapfile);
			if(c=='\t') numtabs++;
		}
		fscanf(soapfile,"%d\t%c\t",&length,&strand);
		c='\n';
		while(c!='\t') c=fgetc(soapfile);
		fscanf(soapfile,"%d\t",&pos);
		if(strand=='-') pos=-pos;
		fprintf(coveragefile,"%d %d\n",pos,length);
		while(c!='\n') c=fgetc(soapfile);
		numalignedreads++;
	}
	fclose(soapfile);
	fclose(coveragefile);
	free(coveragefilename);
	printf("(%.2lf%% aligned reads - %d of %d) OK\n",(((double)numalignedreads)/((double)totalnumreads)*100.0),numalignedreads,totalnumreads);
	printf("> Done!\n");
	exit(0);
}

/*
// TODO: output the number of errors too
// Retrieves the positions of the aligned reads from a file in the BOWTIE output format
void ProcessBowtieFile(char *bowtiefilename, int totalnumreads){
	FILE *bowtiefile,*coveragefile;
	char *coveragefilename;
	int i,pos,length;
	int numalignedreads;
	char c,strand;
	printf("> Opening BOWTIE file <%s> ... ",bowtiefilename);
	fflush(stdout);
	if((bowtiefile=fopen(bowtiefilename,"r"))==NULL){
		printf("\n> ERROR: BOWTIE file not found\n");
		exit(0);
	}
	printf("OK\n");
	i=(int)strlen(bowtiefilename);
	coveragefilename=(char *)calloc((i+8+1),sizeof(char));
	strcpy(coveragefilename,bowtiefilename);
	strcat(coveragefilename,".cov.txt");
	printf("> Creating Coverage file <%s> ... ",coveragefilename);
	fflush(stdout);
	if((coveragefile=fopen(coveragefilename,"w"))==NULL){
		printf("\n> ERROR: Can't write Coverage file\n");
		exit(0);
	}
	numalignedreads=0;
	while((c=fgetc(bowtiefile))!=EOF){
		fputc('>',coveragefile);
		while(c!='\t'){
			fputc(c,coveragefile);
			c=fgetc(bowtiefile);
		}
		fputc('\n',coveragefile);
		fscanf(bowtiefile,"%c\t",&strand);
		c='\n';
		while(c!='\t') c=fgetc(bowtiefile);
		fscanf(bowtiefile,"%d\t",&pos);
		if(strand=='-') pos=-pos;
		length=0;
		while((c=fgetc(bowtiefile))!='\t') length++;
		fprintf(coveragefile,"%d %d\n",pos,length);
		while(c!='\n') c=fgetc(bowtiefile);
		numalignedreads++;
	}
	fclose(bowtiefile);
	fclose(coveragefile);
	free(coveragefilename);
	printf("(%.2lf%% aligned reads - %d of %d) OK\n",(((double)numalignedreads)/((double)totalnumreads)*100.0),numalignedreads,totalnumreads);
	printf("> Done!\n");
	exit(0);
}
*/
