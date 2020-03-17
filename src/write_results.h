
#ifndef WRITE_RESULTSH_INCLUDED
#define WRITE_RESULTSH_INCLUDED

#include "GlobVar.h"

bool write_results(GlobVar& ds, int print_tag, unsigned int print_step, std::chrono::duration<double> elapsed_seconds);

#endif // !WRITE_RESULTSH_INCLUDED
