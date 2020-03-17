
#ifndef COMMONH_INCLUDED
#define COMMONH_INCLUDED

// read file names in Results directory
int findLastStep(const char *path);
void get_domain_size(unsigned int *rown, unsigned int *coln, std::ofstream& logFLUXOSfile);

#endif