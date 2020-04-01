
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

int getIntNumberFromString(std::string s){
   std::stringstream str_strm;
   str_strm << s; //convert the string s into stringstream
   std::string temp_str;
   int temp_int;
   while(!str_strm.eof()) {
      str_strm >> temp_str; //take words into temp_str one by one
      if(std::stringstream(temp_str) >> temp_int) { //try to convert string to int
         return temp_int;
      }
      temp_str = ""; //clear temp string
   }
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

    // Read DEM file fom modset
    std::ifstream file(filename);
    std::string dem_file_temp, msg;
    std::getline(file, dem_file_temp);
    std::getline(file, dem_file_temp);
    file.close();

    unsigned int linei,numi;
    std::string line, stri; 
    std::ifstream myfile (dem_file_temp); //opening the file.

    std::string str_nrows = "NROWS";
    std::string str_ncols = "NCOLS";

    if (myfile.is_open()) //if the file is open
    {

        for(linei=1;linei<=2;linei++) // Just need to read the first 2 lines
        {
            getline (myfile,line); //get one line from the file
            stri = line.substr(0,5);
            numi = getIntNumberFromString(line);

            if (stri.compare(str_nrows)==0){
                *rown = numi;
            }else if(stri.compare(str_ncols)==0){
                *coln = numi;
            }
        }

        myfile.close(); //closing the file
    }
    else std::cout << "Unable to open file: " + filename << std::endl; //if the file is not open output
    
}

