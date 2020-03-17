
#include <vector> 
#include <dirent.h> 
#include <memory> 
#include <string.h>
#include <armadillo>
#include <iostream>

#include "common.h"

std::string SplitFilename (const std::string& str)
{
  size_t found;
  found=str.find_last_of("/\\");
  return str.substr(0,found);

}

// read file names in Results directory
int findLastStep(const char *path) 
{

   struct dirent *entry;
   int i, timestart, filenum = 0, simnum;
   std::vector<char*> filenames; //stringvec filenames, filename_i;
   const char *filename_i;
   char *simnum_str_i;
   DIR *dir = opendir(path);
   
   if (dir != NULL) {
        while ((entry = readdir(dir)) != NULL) {
        filenames.push_back(entry->d_name); // storing the file names
        filenum = filenum + 1;
        }
   }
   closedir(dir);
   
   timestart = 0;
   for(i=2;i<filenum;i++){
       filename_i = filenames[i]; //.assign(filenames[i]); //strcpy(filename_i,(char *)(&filenames[i]));
        simnum_str_i = (char*) malloc(sizeof(filename_i)-2);
        strncpy (simnum_str_i, filename_i, sizeof(filename_i)-2);
        simnum = atoi(simnum_str_i);
        timestart = std::max(timestart,simnum);
        free(simnum_str_i);
   }
   
   //free(filename_i);
   free(entry);
   //free(dir);
   //std::cout << "Start time (s): " << timestart << " (initial conditions available)" << std::endl;
   return timestart;
   
}


// get size of the domain
void get_domain_size(unsigned int *rown, unsigned int *coln, 
                    const std::string& filename, const std::string& pathfile,
                    std::ofstream& logFLUXOSfile)
{

    arma::mat filedata; 
    
    std::ifstream file(filename);
    std::string dem_file_temp, msg;
    std::getline(file, dem_file_temp);
    std::getline(file, dem_file_temp);
    file.close();
    
    bool flstatus =  filedata.load(pathfile + "/" + dem_file_temp,arma::raw_ascii);
   
    *rown = 0;
    *coln = 0;
    
    if(flstatus == true) {
        *rown = filedata.col(1).n_elem;
        *coln = filedata.row(1).n_elem;
    }
 
}

