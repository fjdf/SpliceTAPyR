void InitializeSpliceSitesArray(unsigned int genomeSize);
void FreeSpliceSitesArray();
unsigned int GetOppositePositionIfSpliceSiteExists(unsigned int *position);
unsigned int AddOrUpdateSpliceSite(unsigned int leftPosition, unsigned int rightPosition);
int IsSpliceSite(char *genome, unsigned int gapStartInGenome, int gapSizeInGenome);
#ifdef DEBUG
void PrintSSBlocksStats();
void PrintSpliceSites(char *readsfilename);
#endif
