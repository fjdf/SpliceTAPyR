unsigned int FMI_PositionInText( unsigned int bwtpos );
unsigned int FMI_FollowLetter( char c , unsigned int *topPointer , unsigned int *bottomPointer );
void FMI_LoadIndex(char *indexFilename);
void FMI_SaveIndex(char *indexFilename);
void FMI_FreeIndex();
void FMI_BuildIndex(char *filename, int silentmode);
unsigned int FMI_GetTextSize();
unsigned int FMI_GetBWTSize();
char *FMI_GetTextFilename();
