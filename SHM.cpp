#include <iostream>
#include <fstream>
#include <map>
#include <algorithm> // for replace function
#include "SHM.h"
#include "cell.h"
#include<vector>
#include<string>
#include <random>
#include <regex>
#include <boost/math/distributions/inverse_gamma.hpp>
#include <boost/algorithm/string.hpp>
#include <time.h>
#include <experimental/filesystem>
#include <sstream>
#include <iterator>
#include <experimental/filesystem>

using namespace std;
///R-SHM
#include "random.h"
#include "SHM.h"
///R-SHM


std::map<std::string, std::vector<float>> Seq_affinity;
std::map<std::string, std::map<std::string, std::map<int, std::map<std::string, std::string>>>> Restricted_db;
std::map<std::vector<std::string>, std::vector<float>> CDRs_affinity;
std::map<std::vector<std::string>, std::string> FWRs_db;

static std::mt19937 generator(global_seed+3);
static std::mt19937 generator2(global_seed +1532);
std::default_random_engine generator_fasta(global_seed +1);
static std::mt19937 gen(global_seed +2);


void confirmBCR(vector<float> &BCReceptor, std::string sequence, std::map<std::string, std::string> AA_REGIONS, int mother_ID, int cell_ID) {

    char filename[ ] = "bcinflow09/Sequence_db.csv";
    fstream database;
    std::ifstream database2 ("bcinflow09/Sequence_db.csv");

    char filename_restricted[ ] = "bcinflow09/Restricted_db.csv";
    fstream database_restricted;
    std::ifstream database2_restricted ("bcinflow09/Restricted_db.csv");

    char filename_CDRs[ ] = "bcinflow09/CDRs_db.csv";
    fstream database_CDRs;
    std::ifstream database2_CDRs ("bcinflow09/CDRs_db.csv");

    char filename_FWRs[ ] = "bcinflow09/FWRs_db.csv";
    fstream database_FWRs;
    std::ifstream database2_FWRs ("bcinflow09/FWRs_db.csv");

         if (Seq_affinity.size() == 0) {

             database.open(filename,fstream::app);
              // If file does not exist, Create new file
              if (!database )
              {
                cout << "Cannot open file, file does not exist. Creating new file.." << endl;
                database.close();
                database.open(filename,  fstream::out);
                //database.close();

              }
              database.close();

              database_restricted.open(filename_restricted,fstream::app);
               // If file does not exist, Create new file
               if (!database_restricted )
               {
                 cout << "Cannot open file, file does not exist. Creating new file.." << endl;
                 database_restricted.close();
                 database_restricted.open(filename_restricted,  fstream::out);
                 //database.close();

               }
               database_restricted.close();


               database_CDRs.open(filename_CDRs,fstream::app);
                // If file does not exist, Create new file
                if (!database_CDRs )
                {
                  cout << "Cannot open file, file does not exist. Creating new file.." << endl;
                  database_CDRs.close();
                  database_CDRs.open(filename_CDRs,  fstream::out);
                  //database.close();

                }
                database_CDRs.close();


                database_FWRs.open(filename_FWRs,fstream::app);
                 // If file does not exist, Create new file
                 if (!database_FWRs )
                 {
                   cout << "Cannot open file, file does not exist. Creating new file.." << endl;
                   database_FWRs.close();
                   database_FWRs.open(filename_FWRs,  fstream::out);
                   //database.close();

                 }
                 database_FWRs.close();

              int iy = 0;
              vector<string> row;
              string line, word, temp;


              if (database2.is_open())
              {
                  while (database2.good()) {
                      iy +=1;
                      row.clear();

                      if (getline(database2, line, '\n' )) //test the read for success
                          {
                             if (line.length() != 0) {


                                 istringstream s(line);

                                 while (std::getline(s, word, ',')) {

                                     row.push_back(word);

                                 }
                                 Seq_affinity[row[0]] = {std::stof(row[1]), std::stof(row[2]), std::stof(row[3]),std::stof(row[4])};
                              }
                            else
                              {
                                  cout << "failed to read file" << endl;
                              }


                      }


                  }
                  database2.close();
              }



              int iy_restricted = 0;
              vector<string> row_restricted;
              string line_restricted, word_restricted, temp_restricted;


              if (database2_restricted.is_open())
              {
                  while (database2_restricted.good()) {
                      iy_restricted +=1;
                      row_restricted.clear();

                      if (getline(database2_restricted, line_restricted, '\n' )) //test the read for success
                          {
                             if (line_restricted.length() != 0) {


                                 istringstream s_restricted(line_restricted);

                                 while (std::getline(s_restricted, word_restricted, ',')) {

                                     row_restricted.push_back(word_restricted);

                                 }
                                 Restricted_db[row_restricted[0]][row_restricted[1]][std::stof(row_restricted[2])][row_restricted[3]] = row_restricted[4];
                              }
                            else
                              {
                                  cout << "failed to read file" << endl;
                              }


                      }


                  }
                  database2_restricted.close();
              }


              vector<string> row_CDR;
              string lineCDR, wordCDR;
              if (database2_CDRs.is_open())
              {
                while (database2_CDRs.good()) {
                    row_CDR.clear();
                    if (getline(database2_CDRs, lineCDR, '\n' )) //test the read for success
                        {
                           if (lineCDR.length() != 0) {
                               istringstream sCDR(lineCDR);

                               while (std::getline(sCDR, wordCDR, ',')) {

                                   row_CDR.push_back(wordCDR);
                               }
                               CDRs_affinity[{row_CDR[0], row_CDR[1],row_CDR[2]}] = {std::stof(row_CDR[3]), std::stof(row_CDR[4]), std::stof(row_CDR[5]),std::stof(row_CDR[6])};
                            }
                          else
                            {
                                cout << "failed to read file" << endl;
                            }

                    }

                }
                database2_CDRs.close();
              }

             vector<string> row_FWR;
             string lineFWR, wordFWR;
             if (database2_FWRs.is_open())
             {
               while (database2_FWRs.good()) {
                   row_FWR.clear();
                   if (getline(database2_FWRs, lineFWR, '\n' )) //test the read for success
                       {
                          if (lineFWR.length() != 0) {

                              istringstream sFWR(lineFWR);

                              while (std::getline(sFWR, wordFWR, ',')) {

                                  row_FWR.push_back(wordFWR);
                              }
                              if (row_FWR[1]!= "") {
                                  FWRs_db[{row_FWR[0], row_FWR[1]}] = row_FWR[2];
                              } else {
                                  //cout << lineFWR <<endl;
                              }

                           }
                         else
                           {
                               cout << "failed to read file" << endl;
                           }

                   }

               }
               database2_FWRs.close();
             }
         }

         vector<float> results;




         if (Seq_affinity.find(sequence)!= Seq_affinity.end()) {
             results=Seq_affinity.find(sequence)->second;
             BCReceptor[0]=results[0];
             BCReceptor[1]=results[1];
             BCReceptor[2]=results[2];
             BCReceptor[3]=results[3];
             if(mother_ID==-1) {

                 if (FWRs_db.find({"FWR1", AA_REGIONS["FWR1"]})== FWRs_db.end() ) {
                    FWRs_db[{"FWR1", AA_REGIONS["FWR1"]}]="N";
                 }
                 if (FWRs_db.find({"FWR2", AA_REGIONS["FWR2"]})== FWRs_db.end() ) {
                    FWRs_db[{"FWR2", AA_REGIONS["FWR2"]}]="N";
                 }
                 if (FWRs_db.find({"FWR3", AA_REGIONS["FWR3"]})== FWRs_db.end() ) {
                    FWRs_db[{"FWR3", AA_REGIONS["FWR3"]}]="N";
                 }
                 if (FWRs_db.find({"FWR4", AA_REGIONS["FWR4"]})== FWRs_db.end() ) {
                    FWRs_db[{"FWR4", AA_REGIONS["FWR4"]}]="N";
                 }
                 if (CDRs_affinity.find({AA_REGIONS["CDR1"],AA_REGIONS["CDR2"],AA_REGIONS["CDR3"]}) == CDRs_affinity.end() && FWRs_db[{"FWR1", AA_REGIONS["FWR1"]}]=="N" && FWRs_db[{"FWR2", AA_REGIONS["FWR2"]}]=="N" && FWRs_db[{"FWR3", AA_REGIONS["FWR3"]}]=="N" && FWRs_db[{"FWR4", AA_REGIONS["FWR4"]}]=="N") {
                        CDRs_affinity[{AA_REGIONS["CDR1"],AA_REGIONS["CDR2"],AA_REGIONS["CDR3"]}] =  {BCReceptor[0],BCReceptor[1],BCReceptor[2],BCReceptor[3]};
                 }

                 if (CDRs_affinity.find({AA_REGIONS["CDR1"],AA_REGIONS["CDR2"],AA_REGIONS["CDR3"]}) != CDRs_affinity.end() && FWRs_db[{"FWR1", AA_REGIONS["FWR1"]}]=="N" && FWRs_db[{"FWR2", AA_REGIONS["FWR2"]}]=="N" && FWRs_db[{"FWR3", AA_REGIONS["FWR3"]}]=="N" && FWRs_db[{"FWR4", AA_REGIONS["FWR4"]}]=="N") {
                        if(BCReceptor[0]!=CDRs_affinity[{AA_REGIONS["CDR1"],AA_REGIONS["CDR2"],AA_REGIONS["CDR3"]}][0] ||
                        BCReceptor[1]!=CDRs_affinity[{AA_REGIONS["CDR1"],AA_REGIONS["CDR2"],AA_REGIONS["CDR3"]}][1] ||
                        BCReceptor[2]!=CDRs_affinity[{AA_REGIONS["CDR1"],AA_REGIONS["CDR2"],AA_REGIONS["CDR3"]}][2] ||
                        BCReceptor[3]!=CDRs_affinity[{AA_REGIONS["CDR1"],AA_REGIONS["CDR2"],AA_REGIONS["CDR3"]}][3]){
                            cout<<"Same CDRs, different affinity! Check CDRs "<<AA_REGIONS["CDR1"]<<" "<<AA_REGIONS["CDR2"]<<" "<<AA_REGIONS["CDR3"]<<endl;
                        }
                  }
             }

         } else if (CDRs_affinity.find({AA_REGIONS["CDR1"],AA_REGIONS["CDR2"],AA_REGIONS["CDR3"]}) != CDRs_affinity.end() && FWRs_db[{"FWR1", AA_REGIONS["FWR1"]}]=="N" && FWRs_db[{"FWR2", AA_REGIONS["FWR2"]}]=="N" && FWRs_db[{"FWR3", AA_REGIONS["FWR3"]}]=="N"  && FWRs_db[{"FWR4", AA_REGIONS["FWR4"]}]=="N"){
             //cout<<"CDR_aff"<<endl;
             results=CDRs_affinity.find({AA_REGIONS["CDR1"],AA_REGIONS["CDR2"],AA_REGIONS["CDR3"]})->second;
             BCReceptor[0]=results[0];
             BCReceptor[1]=results[1];
             BCReceptor[2]=results[2];
             BCReceptor[3]=results[3];
             Seq_affinity[sequence] = {BCReceptor[0],BCReceptor[1],BCReceptor[2],BCReceptor[3]};

             if(mother_ID==-1) {
                 cout<<"Founder cell " << cell_ID << "not found in Seq db" <<endl;
             }
         } else {
             //cout<<"Strangely last"<<endl;

             if(FWRs_db[{"FWR1", AA_REGIONS["FWR1"]}]=="L" || FWRs_db[{"FWR2", AA_REGIONS["FWR2"]}]=="L" || FWRs_db[{"FWR3", AA_REGIONS["FWR3"]}]=="L" || FWRs_db[{"FWR4", AA_REGIONS["FWR4"]}]=="L" ){
                  //cout<<"Lethal"<<endl;
                  Seq_affinity[sequence] = {9,9,9,9};

              } else {
                  Seq_affinity[sequence] = {BCReceptor[0],BCReceptor[1],BCReceptor[2],BCReceptor[3]};
              }
             if(mother_ID==-1) {
                 cout<<"Founder cell " << cell_ID << "not found in Seq db" <<endl;
             }
         }

         if (CDRs_affinity.find({AA_REGIONS["CDR1"],AA_REGIONS["CDR2"],AA_REGIONS["CDR3"]}) == CDRs_affinity.end() && FWRs_db[{"FWR1", AA_REGIONS["FWR1"]}]=="N" && FWRs_db[{"FWR2", AA_REGIONS["FWR2"]}]=="N" && FWRs_db[{"FWR3", AA_REGIONS["FWR3"]}]=="N" && FWRs_db[{"FWR4", AA_REGIONS["FWR4"]}]=="N") {
             //cout<<"Added to CDRs_aff"<<endl;
                CDRs_affinity[{AA_REGIONS["CDR1"],AA_REGIONS["CDR2"],AA_REGIONS["CDR3"]}] =  {BCReceptor[0],BCReceptor[1],BCReceptor[2],BCReceptor[3]};
         }




}

void setGermBCR(vector<float> &BCReceptor, vector<float> &Germ_BCR,  int mother_ID) {
    if (mother_ID == -1) {
        Germ_BCR=BCReceptor;
    }
}
///1. Load FASTA files
void add_seq_to_Restricted (std::map<std::string, std::string> & AA_REGIONS, std::string sequence_name) {

    if (Restricted_db.find(sequence_name)== Restricted_db.end()) {
      std::vector<string> regions= {"FWR1","FWR2","FWR3","FWR4"};
      for(int reg =0; reg <  regions.size(); reg++) {

          for(int aas=0; aas < AA_REGIONS[regions[reg]].length() ; aas++) {
              Restricted_db[sequence_name][regions[reg]][aas][string(1,AA_REGIONS[regions[reg]][aas])] = "Neutral";
              if (string(1,AA_REGIONS[regions[reg]][aas]).compare("")==0 || Restricted_db[sequence_name][regions[reg]][aas][string(1,AA_REGIONS[regions[reg]][aas])].compare("")==0) {
                  cout << "add_seq_to_Restricted"<<endl;
              }
          }
      }
    }
}
////1.1. Check DNA from fasta is DNA
bool checkDNA(std::string DNA){
    bool Unvalid = false;
    for (int i = 0; i < DNA.length(); i++) {
        if (DNA[i] != 'C' && DNA[i] != 'A' && DNA[i] != 'G' && DNA[i] != 'T'){
                  Unvalid = true;
                  break;
        }
    }
    return Unvalid;
}


////1.2.DNA to RNA
void findAndReplaceAll(std::string & data, std::string toSearch, std::string replaceStr)
{
	// Get the first occurrence
	size_t pos = data.find(toSearch);
 
	// Repeat till end is reached
	while( pos != std::string::npos)
	{
		// Replace this occurrence of Sub String
		data.replace(pos, toSearch.size(), replaceStr);
		// Get the next occurrence from the current position
		pos =data.find(toSearch, pos + toSearch.size());
	}
}


///FUTURE: 
///1.3. Find paths of fasta files

void findFastas(std::string path) { //Fastas have to end in .fasta


    for (const auto & entry : std::experimental::filesystem::directory_iterator(path)) {
        if (!entry.path().extension().compare(".fasta")) {
             fastas.push_back(entry.path());
        }

    }
}


///1.4. mapToString

std::string mapToString(std::map<std::string, std::vector<int>>  &m) {
    std::stringstream result2;
    int ida = 1;

    for (auto it = m.cbegin(); it != m.cend(); it++) {
        if (ida != 1) {
                result2 << "/";
        }
                result2 << it->first << ":";
                std::copy(it->second.begin(), it->second.end(), std::ostream_iterator<int>(result2, "-"));
        ida +=1;
    }


  return result2.str();
}


std::string mapToString2(std::map<std::string, std::string>  &m) {
    std::stringstream result2;
    int ida = 1;
    for(auto& kv : m) {
        if (ida != 1) {
                result2 << "/";
        }
          result2 << kv.first <<  ':'  << kv.second;

        ida +=1;


    }
    std::string resultados = result2.str();
    if (!resultados.empty() && resultados[resultados.length()-1] == '\n') {
        resultados.erase(resultados.length()-1);
    }

  return result2.str();
}


///1.5 Load sequence, regions and mutations in B-cell

void getSequence(std::map<std::string, std::string> & REGIONS, std::map<std::string, std::vector<int>> & MUTATIONS, std::string &sequence, std::string &germline_name) {
    int argc = 2;

    bool All_not_OK= true;

    do {
        if (fastas.size () == 0) {
             findFastas(fasta_path);
        }
        if (fastas_remaining.size() == 0) {
            fastas_remaining=fastas;
            std::random_device rd;

            static auto rng = std::default_random_engine {};
            std::shuffle(std::begin(fastas_remaining), std::end(fastas_remaining), generator_fasta);
        }

        std::string input = fastas_remaining.back();
        fastas_remaining.pop_back();

        germline_name=std::experimental::filesystem::path(input).filename();

        //const char *argv[] = {"/home/rgarcia/Escritorio/NGly_scripts/Fastas/seq.fasta"};  ///R cambiar en algun momento a parametro
        if( argc <= 1 ){
            std::cerr << "Usage: "<<input<<" [infile]" << std::endl;
        }
        std::ifstream input3(input, std::ios::binary);
        if(!input3.good()){
                std::cerr << "Error opening '"<<input<<"'. Bailing out." << std::endl;
        }
            std::string line, name, content;

            while( std::getline( input3, line )){

                if( line[0] == '>' ){ // Identifier marker
                    if( !line.empty() ){
                        name = line.substr(1);
                    }
                    content.clear();

                } else if( !name.empty() ){
                    boost::to_upper(line);
                    if( line.find(' ') != std::string::npos ){ // Invalid sequence--no spaces allowed
                        std::cerr << "Error, check FASTA file. Bailing out." << std::endl;
                    } else {
                        if ( checkDNA(line)) {
                            std::cerr << "Error, check strange characters in FASTA line :" << line << ". Bailing out." << std::endl;
                        };
                        findAndReplaceAll(line, "T", "U");
                        content += line;
                        REGIONS[name]=content;
                    }
                } else if ( line.empty()) {
                    std::cerr << "Error, check FASTA file. Bailing out." << std::endl;
                }
            }

            if (REGIONS.size() != 7) {
                std::cerr << "Error, check FASTA file FWR and CDR nomenclature. Bailing out." << std::endl;
            }

        sequence= REGIONS["FWR1"] +REGIONS["CDR1"] +REGIONS["FWR2"] +REGIONS["CDR2"] +REGIONS["FWR3"]  +REGIONS["CDR3"] +REGIONS["FWR4"] ;

        All_not_OK= false;
        if(isThereAnStop(DNAtoprotein(sequence))) {
            cout  << "A seq "<< input << " has no good 1st ORF" << endl;
            std::map<std::string, std::string> tmporal =REGIONS;
            REGIONS["FWR1"].erase(0,1);
            sequence= REGIONS["FWR1"] +REGIONS["CDR1"] +REGIONS["FWR2"] +REGIONS["CDR2"] +REGIONS["FWR3"]  +REGIONS["CDR3"] +REGIONS["FWR4"] ;
            if(isThereAnStop(DNAtoprotein(sequence))) {

                REGIONS["FWR1"].erase(0,1);
                sequence= REGIONS["FWR1"] +REGIONS["CDR1"] +REGIONS["FWR2"] +REGIONS["CDR2"] +REGIONS["FWR3"]  +REGIONS["CDR3"] +REGIONS["FWR4"] ;

            } if(isThereAnStop(DNAtoprotein(sequence))) {
                cout << "A seq "<< input << " has no valid ORF" << endl;
                All_not_OK= true;

            }
        }

    } while (All_not_OK);

    MUTATIONS = {
        { "FWR1", vector<int>(REGIONS["FWR1"].length()) },
        { "FWR2", vector<int>(REGIONS["FWR2"].length()) },
        { "FWR3", vector<int>(REGIONS["FWR3"].length()) },
        { "FWR4", vector<int>(REGIONS["FWR4"].length()) },
        { "CDR1", vector<int>(REGIONS["CDR1"].length()) }, 
        { "CDR2", vector<int>(REGIONS["CDR2"].length()) }, 
        { "CDR3", vector<int>(REGIONS["CDR3"].length()) }, 
    };


}


///2. Conversion to AA, NGly detection, STOP detection

////2.1. Conversion to AA

std::string DNAtoprotein (std::string newSequence) {
    int firstBase, secondBase, thirdBase;
    std::string protein ="";
    char aminoAcid[4][4][4];        ////CONFIRMAR QUE ESTA BIEN! :D
    //A = 0, C = 1, G = 2, U = 3

    //phenylalanine - F
    aminoAcid[3][3][3] = 'F';
        aminoAcid[3][3][1] = 'F';
    //Leucine - L
        aminoAcid[3][3][0] = 'L';
        aminoAcid[3][3][2] = 'L';
    //Serine - S
        aminoAcid[3][1][3] = 'S';
        aminoAcid[3][1][1] = 'S';
        aminoAcid[3][1][0] = 'S';
        aminoAcid[3][1][2] = 'S';
    //tyrosine - Y
        aminoAcid[3][0][3] = 'Y';
        aminoAcid[3][0][1] = 'Y';
    //stop codon
        aminoAcid[3][0][0] = '.';
        aminoAcid[3][0][2] = '.';
    //cysteine - C
        aminoAcid[3][2][3] = 'C';
        aminoAcid[3][2][1] = 'C';
    //stop codon
        aminoAcid[3][2][0] = '.';
    //tryptophan - W
        aminoAcid[3][2][2] = 'W';
    //leucine - L
        aminoAcid[1][3][3] = 'L';
        aminoAcid[1][3][1] = 'L';
        aminoAcid[1][3][0] = 'L';
        aminoAcid[1][3][2] = 'L';
    //proline - P
        aminoAcid[1][1][3] = 'P';
        aminoAcid[1][1][1] = 'P';
        aminoAcid[1][1][0] = 'P';
        aminoAcid[1][1][2] = 'P';
    //histidine - H
        aminoAcid[1][0][3] = 'H';
        aminoAcid[1][0][1] = 'H';
    //glutamine - Q
        aminoAcid[1][0][0] = 'Q';
        aminoAcid[1][0][2] = 'Q';
    //arginine - R
        aminoAcid[1][2][3] = 'R';
        aminoAcid[1][2][1] = 'R';
        aminoAcid[1][2][0] = 'R';
        aminoAcid[1][2][2] = 'R';
    //isoleucine - I
        aminoAcid[0][3][3] = 'I';
        aminoAcid[0][3][1] = 'I';
        aminoAcid[0][3][0] = 'I';
    //methionine(start codon) - M
        aminoAcid[0][3][2] = 'M';
    //threonine -T
        aminoAcid[0][1][3] = 'T';
        aminoAcid[0][1][1] = 'T';
        aminoAcid[0][1][0] = 'T';
        aminoAcid[0][1][2] = 'T';
    //asparagine - N
        aminoAcid[0][0][3] = 'N';
        aminoAcid[0][0][1] = 'N';
    //lysine - K
        aminoAcid[0][0][0] = 'K';
        aminoAcid[0][0][2] = 'K';
    //serine - S
        aminoAcid[0][2][3] = 'S';
        aminoAcid[0][2][1] = 'S';
    //arginine - R
        aminoAcid[0][2][0] = 'R';
        aminoAcid[0][2][2] = 'R';
    //valine - V
        aminoAcid[2][3][3] = 'V';
        aminoAcid[2][3][1] = 'V';
        aminoAcid[2][3][0] = 'V';
        aminoAcid[2][3][2] = 'V';
    //alanine - A
        aminoAcid[2][1][3] = 'A';
        aminoAcid[2][1][1] = 'A';
        aminoAcid[2][1][0] = 'A';
        aminoAcid[2][1][2] = 'A';
    //aspartic acid - D
        aminoAcid[2][0][3] = 'D';
        aminoAcid[2][0][1] = 'D';
    //glutamic acid - E
        aminoAcid[2][0][0] = 'E';
        aminoAcid[2][0][2] = 'E';
    //glycine - G
        aminoAcid[2][2][3] = 'G';
        aminoAcid[2][2][1] = 'G';
        aminoAcid[2][2][0] = 'G';
        aminoAcid[2][2][2] = 'G';

    bool isthereadot = false;
    for(int i = 0; i < newSequence.length() - 2; i += 3)
    {
        if(newSequence[i] == 'A')
        {
            firstBase = 0;
        }
        else if(newSequence[i] == 'C')
        {
            firstBase = 1;
        }
        else if(newSequence[i] == 'G')
        {
            firstBase = 2;
        }
        else if(newSequence[i] == 'U')
        {
            firstBase = 3;
        } else {
            isthereadot = true;
        }

        if(newSequence[i+1] == 'A')
        {
            secondBase = 0;
        }
        else if(newSequence[i+1] == 'C')
        {
            secondBase = 1;
        }
        else if(newSequence[i+1] == 'G')
        {
            secondBase = 2;
        }
        else if(newSequence[i+1] == 'U')
        {
            secondBase = 3;
        } else {
            isthereadot = true;
        }


        if(newSequence[i+2] == 'A')
        {
            thirdBase = 0;
        }
        else if(newSequence[i+2] == 'C')
        {
            thirdBase = 1;
        }
        else if(newSequence[i+2] == 'G')
        {
            thirdBase = 2;
        }
        else if(newSequence[i+2] == 'U')
        {
            thirdBase = 3;
        } else {
            isthereadot = true;
        }


        bool readSequence = true;

        if (aminoAcid[firstBase][secondBase][thirdBase] == aminoAcid[0][3][2])
        {
            readSequence = true; //hace falta hacer AUG translation?
        }

        if (isthereadot) {
            protein= protein + ".";
            ///break;
        } else {
            if(readSequence)
            {
                protein = protein + aminoAcid[firstBase][secondBase][thirdBase];
            }
            else
            {
                continue;
            }
        }

    }
    return protein;
}

////2.2. NGly detection

std::map<std::string, std::string> Separating(std::map<std::string, vector<int>> MUTATIONS, std::string protein_sequence) {

    std::map<std::string, std::string> AA_REGIONS;

    int FWR1=((MUTATIONS["FWR1"].size()+1)/3) ;
    int FWR2=((MUTATIONS["FWR2"].size()+1)/3) ;
    int FWR3=((MUTATIONS["FWR3"].size()+1)/3) ;
    int FWR4=((MUTATIONS["FWR4"].size()+1)/3) ;
    int CDR1=((MUTATIONS["CDR1"].size()+1)/3) ;
    int CDR2=((MUTATIONS["CDR2"].size()+1)/3) ;
    int CDR3=((MUTATIONS["CDR3"].size()+1)/3) ;


    AA_REGIONS["FWR1"]=protein_sequence.substr(0,FWR1);
    AA_REGIONS["CDR1"]=protein_sequence.substr(FWR1,CDR1);
    AA_REGIONS["FWR2"]=protein_sequence.substr(FWR1+CDR1,FWR2);
    AA_REGIONS["CDR2"]=protein_sequence.substr(FWR1+CDR1+FWR2,CDR2);
    AA_REGIONS["FWR3"]=protein_sequence.substr(FWR1+CDR1+FWR2+CDR2,FWR3);
    AA_REGIONS["CDR3"]=protein_sequence.substr(FWR1+CDR1+FWR2+CDR2+FWR3,CDR3);
    AA_REGIONS["FWR4"]=protein_sequence.substr(FWR1+CDR1+FWR2+CDR2+FWR3+CDR3,FWR4);

    if((AA_REGIONS["FWR1"] + AA_REGIONS["CDR1"]+ AA_REGIONS["FWR2"] +AA_REGIONS["CDR2"] +AA_REGIONS["FWR3"] + AA_REGIONS["CDR3"] + AA_REGIONS["FWR4"]) != protein_sequence) {
        std::cout<<"Strangely " << AA_REGIONS["FWR1"] + AA_REGIONS["CDR1"]+ AA_REGIONS["FWR2"] +AA_REGIONS["CDR2"] +AA_REGIONS["FWR3"] + AA_REGIONS["CDR3"] + AA_REGIONS["FWR4"] << "not the same as "<< protein_sequence << std::endl;
    }
    return AA_REGIONS;
}


std::vector <int> Blacklisting(std::map<std::string, vector<int>> MUTATIONS) {
    vector <int> blacklist;
    int FWR1=((MUTATIONS["FWR1"].size()+1)/3) ;
    int FWR2=((MUTATIONS["FWR2"].size()+1)/3) ;
    int FWR3=((MUTATIONS["FWR3"].size()+1)/3) ;
    int FWR4=((MUTATIONS["FWR4"].size()+1)/3) ;
    int CDR1=((MUTATIONS["CDR1"].size()+1)/3) ;
    int CDR2=((MUTATIONS["CDR2"].size()+1)/3) ;
    int CDR3=((MUTATIONS["CDR3"].size()+1)/3) ;


    blacklist.resize(FWR1);
    std::iota(blacklist.begin(), blacklist.end(), 0);

    blacklist.resize(blacklist.size()+FWR2);
    std::iota(blacklist.begin()+FWR1, blacklist.end(), CDR1+FWR1);

    blacklist.resize(blacklist.size()+FWR3);
    std::iota(blacklist.begin()+FWR1+FWR2, blacklist.end(), FWR1+CDR1+FWR2+CDR2);

    blacklist.resize(blacklist.size()+FWR4);
    std::iota(blacklist.begin()+FWR1+FWR2+FWR3, blacklist.end(), FWR1+CDR1+FWR2+CDR2+FWR3+CDR3);

    return blacklist;
}
std::vector<int> NGly_positions (std::string sequence, std::vector<int> blacklist, std::vector<int> &FULL_NGly_sites) {
    std::string pattern("N[^P](S|T)");         // Regex expression
    std::regex rx(pattern);             // Getting the regex object

    /////Get number matches (no needed, can get length indexes)
    //std::ptrdiff_t number_of_matches = std::distance(std::sregex_iterator(sequence.begin(), sequence.end(), rx),std::sregex_iterator());

    //std::cout << number_of_matches << std::endl;  // Displaying results



    ////Get position indexes EMPIEZA EN 0
    vector<int> index_matches; // results saved here
                               // (should be {2, 8}, but always get only {2})

    for(auto it = std::sregex_iterator(sequence.begin(), sequence.end(), rx);
        it != std::sregex_iterator();      ++it)
    {
        index_matches.push_back(it->position());
    }
    //std::cout << "Numero de N-Gly sites: "<< index_matches.size()<< std::endl;

    FULL_NGly_sites = index_matches;


    std::vector<int> filtered;
    //not need to sort since it already sorted
    std::set_difference(index_matches.begin(), index_matches.end(), blacklist.begin(), blacklist.end(), std::inserter(filtered, filtered.begin()));

    //for (int x = 0; x != index_matches.size(); ++x)
    //{
    //     std::cout << index_matches[x] << std::endl;
   //      //std::cout << dar.at(x) << "- calling at member" << std::endl;
    //}
    return filtered;
}




////2.3. STOP detection
bool isThereAnStop (std::string sequence) {
    bool thereIsAStop = false;

    if (sequence.find('.') != std::string::npos)
    {
        thereIsAStop = true;
    }
    return thereIsAStop;
}



///3. Mutation and affinities





////3.1. Choose which cells will mutate using poisson distribution
std::vector<int> whichCellsToMutate(int new_n_cells, double lambda) {
    std::vector<int> p;
    static std::random_device rd;


    std::poisson_distribution<int> pd(lambda);

    for (int i = 0; i < new_n_cells; ++i){
        p.push_back(pd(gen));
    }

    return p;
}

////3.2 SHM tree

void randomMutWhere(int num_of_mutations, std::map<std::string, std::string> & REGIONS, std::map<std::string, std::vector<int>> & MUTATIONS, std::map<std::string, std::vector<int>> & SILENT_MUTATIONS, std::string & sequence, std::string & protein_sequence, std::map<std::string, double> RATES_OF_MUTATIONS,std::string & sequence_name) {
    //std::cout<<sequence<<std::endl;
    std::vector<string> regions= {"FWR1","CDR1","FWR2","CDR2","FWR3","CDR3","FWR4"};

    std::random_device gen;

    std::uniform_real_distribution<double> distribution(0.0,std::nextafter(1.0, std::numeric_limits<double>::max()));


    float probRegion;
    bool Replacement;
    bool Lethal;
    int reg;
    struct timespec ts;



    std::string originalSeq;
    std::string originalSeqAA;
    std::string mutatedSeq;
    int mutationPlace;
    std::random_device gen2;

    std::uniform_int_distribution<int> point(0,1);
    std::map<std::string, std::string> REGIONS_TMP= REGIONS;
    int prueba = 0;
    Lethal = false;
    bool Keep_searching=true;
    bool Restricted_is_OK=false;
    vector<string> regiones_F={"FWR1","FWR2", "FWR3","FWR4"};
    std::map<std::string, std::string> AA_REGIONS_TMP;
    std::map<std::string, std::string> AA_REGIONS;
    int MutationPlace;
    string AA;
    //cout << "Num of mutations: " << num_of_mutations;
    for (int i = 0; i < num_of_mutations; ++i){

        Replacement = false;
        Lethal=false;
            if (distribution(generator) <= RATES_OF_MUTATIONS["alpha_FWR"] && !Lethal) {
                probRegion= distribution(generator);
                if (probRegion <= RATES_OF_MUTATIONS["alpha_FWR1"]) {
                    reg=0;
                    if (distribution(generator) <= RATES_OF_MUTATIONS["beta_FWR1"]) {
                        Replacement=true;
                    }
                } else if (probRegion <= (RATES_OF_MUTATIONS["alpha_FWR1"]+RATES_OF_MUTATIONS["alpha_FWR2"])) {
                    reg=2;
                    if (distribution(generator) <= RATES_OF_MUTATIONS["beta_FWR2"]) {
                        Replacement=true;
                    }
                } else if (probRegion <= (RATES_OF_MUTATIONS["alpha_FWR1"]+RATES_OF_MUTATIONS["alpha_FWR2"]+RATES_OF_MUTATIONS["alpha_FWR3"])) {
                    reg=4;
                    if (distribution(generator) <= RATES_OF_MUTATIONS["beta_FWR3"]) {
                        Replacement=true;
                    }
                } else {
                    reg=6;
                    if (distribution(generator) <= RATES_OF_MUTATIONS["beta_FWR4"]) {
                        Replacement=true;
                    }
                }



            } else {

                probRegion= distribution(generator);
                if (probRegion <= RATES_OF_MUTATIONS["alpha_CDR1"]) {
                    reg=1;
                    if (distribution(generator) <= RATES_OF_MUTATIONS["gamma_CDR1"]) {
                        Replacement=true;
                    }
                } else if (probRegion <= (RATES_OF_MUTATIONS["alpha_CDR1"]+RATES_OF_MUTATIONS["alpha_CDR2"])) {
                    reg=3;
                    if (distribution(generator) <= RATES_OF_MUTATIONS["gamma_CDR2"]) {
                        Replacement=true;
                    }
                } else {
                    reg=5;
                    if (distribution(generator) <= RATES_OF_MUTATIONS["gamma_CDR3"]) {
                        Replacement=true;
                    }
                }

            }



            if (Replacement==true && (reg == 0 || reg == 2 || reg == 4|| reg == 6) && distribution(generator) <= RATES_OF_MUTATIONS["delta"]){
                Lethal = true;
            }
              //std::cout<<"igual es a proteina"<<std::endl;

            REGIONS_TMP= REGIONS;
            Keep_searching=true;
            if (Replacement && !Lethal) { ////////|| (reg == 1 || reg == 3 || reg == 5), tener en cuenta que puede haber stops en los CDR... y ver que se esta adicionando a restricted db neutral
                //cout << "R no L "  << regions[reg] << endl;
                point=std::uniform_int_distribution<int>(0,REGIONS[regions[reg]].length()-1);
                mutationPlace = point(generator2);
                originalSeq=REGIONS["FWR1"] +REGIONS["CDR1"] +REGIONS["FWR2"] +REGIONS["CDR2"] +REGIONS["FWR3"]  +REGIONS["CDR3"] +REGIONS["FWR4"];
                originalSeqAA=DNAtoprotein(originalSeq);
                mutatedSeq=originalSeq;
                prueba=0;
                //std::cout<<"Replacement"<<std::endl;
                while (Keep_searching) {
                    Restricted_is_OK=false;
                    if (prueba > 25) {  //hay aminoacidos que no se mutan pase lo que pase en una posicion
                        //std::cout << "mutacion previa en" << mutationPlace << std::endl;
                        mutationPlace = point(generator2);
                        prueba=0;

                        //std::cout << "cambio mutacion a " << mutationPlace << std::endl;
                    }
                    REGIONS_TMP= REGIONS;
                    prueba++;
                    //clock_gettime(CLOCK_MONOTONIC, &ts);

                    /* using nano-seconds instead of seconds */
                    //srand((time_t)ts.tv_nsec); ////RRRR-remove
                      //srand ( time(NULL) ); ///RRR Se puede optimizar mas
                      const char arrayNucleotides[4] = {'A', 'U', 'G', 'C'};
                      //int RandIndex = rand() % 4;
                      int RandIndex = random::randomInteger(4);
                      REGIONS_TMP[regions[reg]][mutationPlace] = arrayNucleotides[RandIndex];
                      mutatedSeq=REGIONS_TMP["FWR1"] +REGIONS_TMP["CDR1"] +REGIONS_TMP["FWR2"] +REGIONS_TMP["CDR2"] +REGIONS_TMP["FWR3"]  +REGIONS_TMP["CDR3"] +REGIONS_TMP["FWR4"];
                      AA_REGIONS_TMP= Separating(MUTATIONS, DNAtoprotein(mutatedSeq));
                      AA_REGIONS= Separating(MUTATIONS, originalSeqAA);

                      if(AA_REGIONS_TMP[regions[reg]].compare(AA_REGIONS[regions[reg]])!= 0) {
                          if (std::find(regiones_F.begin(), regiones_F.end(), regions[reg]) != regiones_F.end()) {

                              int confirmar=0;
                              int position_AA;
                              for (unsigned li = 0; li < AA_REGIONS[regions[reg]].size(); li++ ) {
                                  if (AA_REGIONS[regions[reg]][li] != AA_REGIONS_TMP[regions[reg]][li]) {
                                      MutationPlace=li;
                                      AA=AA_REGIONS_TMP[regions[reg]][li];
                                      position_AA=li;
                                      confirmar++;
                                  }
                              }

                              if (confirmar > 1) {
                              cout << "Shady" << endl;

                              } else if (confirmar == 0) {
                                  cout << "Shady 2" << endl;

                              }

                              if ( Restricted_db[sequence_name][regions[reg]][MutationPlace].find(AA) != Restricted_db[sequence_name][regions[reg]][MutationPlace].end() ) {
                                  if (Restricted_db[sequence_name][regions[reg]][MutationPlace][AA].compare("Lethal") != 0) {

                                              if (FWRs_db.find({regions[reg], AA_REGIONS_TMP[regions[reg]]})!= FWRs_db.end() ) {
                                                  if (FWRs_db[{regions[reg], AA_REGIONS_TMP[regions[reg]]}].compare("N")!= 0 && FWRs_db[{regions[reg], AA_REGIONS[regions[reg]]}].compare("N")== 0 && num_of_mutations == 1) {
                                                      cout << "Error with FWRS! 3 "  << regions[reg] << " " << AA_REGIONS_TMP[regions[reg]] << " " << FWRs_db[{regions[reg], AA_REGIONS_TMP[regions[reg]]}] << " " << AA_REGIONS[regions[reg]] << " " << FWRs_db[{regions[reg], AA_REGIONS[regions[reg]]}] << endl;
                                                  }

                                              } else { /////Y si es R no L en una zona L de un FWR?
                                                  if (FWRs_db[{regions[reg], AA_REGIONS[regions[reg]]}].compare("N")==0) {
                                                      FWRs_db[{regions[reg], AA_REGIONS_TMP[regions[reg]]}] ="N";
                                                  } else {
                                                      int num_L=0;
                                                      string AA_mother;
                                                      for (unsigned li = 0; li < AA_REGIONS[regions[reg]].size(); li++ ) {
                                                          AA_mother=AA_REGIONS[regions[reg]][li];
                                                          if (Restricted_db[sequence_name][regions[reg]][li][AA_mother].compare("Lethal") == 0) {
                                                              num_L++;
                                                          }
                                                      }
                                                      AA_mother=AA_REGIONS[regions[reg]][position_AA];
                                                      if (num_L==1 && Restricted_db[sequence_name][regions[reg]][position_AA][AA_mother].compare("Lethal") == 0){
                                                          FWRs_db[{regions[reg], AA_REGIONS_TMP[regions[reg]]}] ="N";
                                                      } else {
                                                          FWRs_db[{regions[reg], AA_REGIONS_TMP[regions[reg]]}] ="L";
                                                      }

                                                  }

                                              }
                                              Restricted_is_OK=true;

                                  } else {
                                      Restricted_is_OK=false;

                                  }
                              } else {
                                  if(AA.compare(".")!= 0){
                                      if (FWRs_db.find({regions[reg], AA_REGIONS_TMP[regions[reg]]})!= FWRs_db.end() ) {
                                          if (FWRs_db[{regions[reg], AA_REGIONS_TMP[regions[reg]]}].compare("N")!= 0) {
                                              cout << "Error with FWRS! 2 " << regions[reg] << " " << AA_REGIONS_TMP[regions[reg]] << " " << FWRs_db[{regions[reg], AA_REGIONS_TMP[regions[reg]]}] << endl;
                                          }
                                      } else { /////Y si es R no L en una zona L de un FWR?
                                          if (FWRs_db[{regions[reg], AA_REGIONS[regions[reg]]}].compare("N")==0) {
                                              FWRs_db[{regions[reg], AA_REGIONS_TMP[regions[reg]]}] ="N";
                                          } else {
                                              int num_L=0;
                                              string AA_mother;
                                              for (unsigned li = 0; li < AA_REGIONS[regions[reg]].size(); li++ ) {
                                                  AA_mother=AA_REGIONS[regions[reg]][li];
                                                  if (Restricted_db[sequence_name][regions[reg]][li][AA_mother].compare("Lethal") == 0) {
                                                      num_L++;
                                                  }
                                              }
                                              AA_mother=AA_REGIONS[regions[reg]][position_AA];
                                              if (num_L==1 && Restricted_db[sequence_name][regions[reg]][position_AA][AA_mother].compare("Lethal") == 0){
                                                  FWRs_db[{regions[reg], AA_REGIONS_TMP[regions[reg]]}] ="N";
                                              } else {
                                                  FWRs_db[{regions[reg], AA_REGIONS_TMP[regions[reg]]}] ="L";
                                              }

                                          }
                                      }
                                      Restricted_db[sequence_name][regions[reg]][MutationPlace][AA]= "Neutral";
                                      Restricted_is_OK=true;
                                      if (Restricted_db[sequence_name][regions[reg]][MutationPlace][AA].compare("")==0) {
                                          cout << "FALLO2"<<endl;
                                      }
                                  } else {
                                      if (FWRs_db.find({regions[reg], AA_REGIONS_TMP[regions[reg]]})!= FWRs_db.end() ) {
                                          if (FWRs_db[{regions[reg], AA_REGIONS_TMP[regions[reg]]}].compare("L")!= 0) {
                                              cout << "Error with FWRS!"  << regions[reg] << " " << AA_REGIONS_TMP[regions[reg]] << " " << FWRs_db[{regions[reg], AA_REGIONS_TMP[regions[reg]]}] << endl;
                                          }
                                      } else {
                                          FWRs_db[{regions[reg], AA_REGIONS_TMP[regions[reg]]}] ="L";
                                      }
                                      Restricted_db[sequence_name][regions[reg]][MutationPlace][AA]= "Lethal";
                                      Restricted_is_OK=false;
                                      if (Restricted_db[sequence_name][regions[reg]][MutationPlace][AA].compare("")==0) {
                                          cout << "FALLO1"<<endl;
                                      }
                                  }

                              }

                          } else{
                              if(!isThereAnStop(originalSeqAA)) {
                                  if(!isThereAnStop(DNAtoprotein(mutatedSeq))) {
                                      Restricted_is_OK=true;
                                  }
                              } else {
                                Restricted_is_OK=true;
                              }

                          }

                      }


                      if(DNAtoprotein(mutatedSeq).compare(originalSeqAA)!= 0) {
                          string mutatedSeqAA=DNAtoprotein(mutatedSeq);

                          if (std::count(mutatedSeqAA.begin(), mutatedSeqAA.end(), '.') == std::count(originalSeqAA.begin(), originalSeqAA.end(), '.')){ ///////y si es en CDRs?
                              if (Restricted_is_OK){
                                  Keep_searching = false;
                              }
                          }

                      }


               }
                REGIONS=REGIONS_TMP;
                MUTATIONS[regions[reg]][mutationPlace] = MUTATIONS[regions[reg]][mutationPlace] + 1;
           }

            if (!Replacement && !Lethal) {
                //cout << "Silenciosa " << regions[reg] << endl;
                point=std::uniform_int_distribution<int>(0,REGIONS[regions[reg]].length()-1);
                mutationPlace = point(generator2);
                originalSeq=REGIONS["FWR1"] +REGIONS["CDR1"] +REGIONS["FWR2"] +REGIONS["CDR2"] +REGIONS["FWR3"]  +REGIONS["CDR3"] +REGIONS["FWR4"];
                originalSeqAA=DNAtoprotein(originalSeq);
                mutatedSeq=originalSeq;
                prueba=0;
                //std::cout<<"Replacement"<<std::endl;

                do {
                    if (prueba > 25) {  //hay aminoacidos que no se mutan pase lo que pase en una posicion
                        //std::cout << "mutacion previa en" << mutationPlace << std::endl;
                        mutationPlace = point(generator2);
                        prueba=0;

                        //std::cout << "cambio mutacion a " << mutationPlace << std::endl;
                    }
                    REGIONS_TMP= REGIONS;
                    //clock_gettime(CLOCK_MONOTONIC, &ts);
                    /* using nano-seconds instead of seconds */
                    //srand((time_t)ts.tv_nsec);
                      //srand ( time(NULL) ); ///RRR Se puede optimizar mas
                      const char arrayNucleotides[4] = {'A', 'U', 'G', 'C'};
                      //int RandIndex = rand() % 4;
                      int RandIndex = random::randomInteger(4);
                      if (arrayNucleotides[RandIndex] != REGIONS[regions[reg]][mutationPlace]) {
                          REGIONS_TMP[regions[reg]][mutationPlace] = arrayNucleotides[RandIndex];
                          mutatedSeq=REGIONS_TMP["FWR1"] +REGIONS["CDR1"] +REGIONS_TMP["FWR2"] +REGIONS_TMP["CDR2"] +REGIONS_TMP["FWR3"]  +REGIONS_TMP["CDR3"] +REGIONS_TMP["FWR4"];
                      }
                     prueba++;
                }
                while (DNAtoprotein(mutatedSeq).compare(originalSeqAA)!= 0);
                REGIONS=REGIONS_TMP;
                SILENT_MUTATIONS[regions[reg]][mutationPlace] = SILENT_MUTATIONS[regions[reg]][mutationPlace] + 1;

           }

            if (Replacement && Lethal) {
                //cout << "R  L "  << regions[reg] << endl;
                point=std::uniform_int_distribution<int>(0,REGIONS[regions[reg]].length()-1);
                 mutationPlace = point(generator2);
                 originalSeq=REGIONS["FWR1"] +REGIONS["CDR1"] +REGIONS["FWR2"] +REGIONS["CDR2"] +REGIONS["FWR3"]  +REGIONS["CDR3"] +REGIONS["FWR4"];
                 originalSeqAA=DNAtoprotein(originalSeq);
                 mutatedSeq=originalSeq;
                 prueba=0;

                 while (Keep_searching) {
                     Restricted_is_OK=false;
                     if (prueba > 25) {  //hay aminoacidos que no se mutan pase lo que pase en una posicion
                         //std::cout << "mutacion previa en" << mutationPlace << std::endl;
                         mutationPlace = point(generator2);
                         prueba=0;

                         //std::cout << "cambio mutacion a " << mutationPlace << std::endl;
                     }
                     REGIONS_TMP= REGIONS;
                     prueba++;

                       const char arrayNucleotides[4] = {'A', 'U', 'G', 'C'};
                       int RandIndex = random::randomInteger(4);
                       REGIONS_TMP[regions[reg]][mutationPlace] = arrayNucleotides[RandIndex];
                       mutatedSeq=REGIONS_TMP["FWR1"] +REGIONS_TMP["CDR1"] +REGIONS_TMP["FWR2"] +REGIONS_TMP["CDR2"] +REGIONS_TMP["FWR3"]  +REGIONS_TMP["CDR3"] +REGIONS_TMP["FWR4"];
                       AA_REGIONS_TMP= Separating(MUTATIONS, DNAtoprotein(mutatedSeq));
                       AA_REGIONS= Separating(MUTATIONS, originalSeqAA);


                       if(AA_REGIONS_TMP[regions[reg]].compare(AA_REGIONS[regions[reg]])!= 0) {
                           if (std::find(regiones_F.begin(), regiones_F.end(), regions[reg]) != regiones_F.end()) {


                               int confirmar=0;
                               for (unsigned li = 0; li < AA_REGIONS[regions[reg]].size(); li++ ) {
                                   if (string(1,AA_REGIONS[regions[reg]][li]).compare(string(1,AA_REGIONS_TMP[regions[reg]][li])) != 0) {
                                       MutationPlace=li;
                                       AA=AA_REGIONS_TMP[regions[reg]][li];
                                       confirmar++;
                                   }
                               }

                               if (confirmar > 1) {
                               cout << "Shady" << endl;

                               } else if (confirmar == 0) {
                                   cout << "Shady2" << endl;

                               }

                               Restricted_is_OK=false;
                               if ( Restricted_db[sequence_name][regions[reg]][MutationPlace].find(AA) != Restricted_db[sequence_name][regions[reg]][MutationPlace].end() ) {

                                   if (Restricted_db[sequence_name][regions[reg]][MutationPlace][AA].compare("Lethal")== 0) {
                                               Restricted_is_OK=true;
                                               if (FWRs_db.find({regions[reg], AA_REGIONS_TMP[regions[reg]]})!= FWRs_db.end() ) {
                                                   if (FWRs_db[{regions[reg], AA_REGIONS_TMP[regions[reg]]}].compare("L")== 0) {
                                                       Restricted_is_OK=true;
                                                   } else {
                                                       cout <<"Error incongruency FWRs "  << regions[reg] << " " << AA_REGIONS_TMP[regions[reg]] << AA_REGIONS[regions[reg]] << " " << FWRs_db[{regions[reg], AA_REGIONS_TMP[regions[reg]]}] << endl;
                                                   }
                                               } else {
                                                   FWRs_db[{regions[reg], AA_REGIONS_TMP[regions[reg]]}] ="L";
                                                   Restricted_is_OK=true;
                                               }
                                       }

                               } else {
                                   if (FWRs_db.find({regions[reg], AA_REGIONS_TMP[regions[reg]]})!= FWRs_db.end() ) {
                                       if (FWRs_db[{regions[reg], AA_REGIONS_TMP[regions[reg]]}].compare("L")== 0) {
                                           Restricted_db[sequence_name][regions[reg]][MutationPlace][AA]= "Lethal";
                                           Restricted_is_OK=true;
                                       }
                                   } else {
                                       FWRs_db[{regions[reg], AA_REGIONS_TMP[regions[reg]]}] ="L";
                                       Restricted_db[sequence_name][regions[reg]][MutationPlace][AA]= "Lethal";
                                       Restricted_is_OK=true;
                                   }


                               }

                           }
                       }


                       if(DNAtoprotein(mutatedSeq).compare(originalSeqAA)!= 0) {
                           if (Restricted_is_OK){
                                   Keep_searching = false;
                           }
                       }


                }
                 REGIONS=REGIONS_TMP;
                 MUTATIONS[regions[reg]][mutationPlace] = MUTATIONS[regions[reg]][mutationPlace] + 1;
            }

            sequence=REGIONS["FWR1"] +REGIONS["CDR1"] +REGIONS["FWR2"] +REGIONS["CDR2"] +REGIONS["FWR3"]  +REGIONS["CDR3"] +REGIONS["FWR4"];
            protein_sequence=DNAtoprotein(sequence);
            AA_REGIONS= Separating(MUTATIONS, protein_sequence);






    }


}


///4. Effect of sequence changes

////4.1. Mutation changes affinity


////4.2. NGly changes affinity
