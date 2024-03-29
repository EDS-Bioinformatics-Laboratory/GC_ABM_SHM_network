#ifndef SHM_H
#define SHM_H

#include <iostream>
#include <fstream>
#include <map>
#include <algorithm> // for replace function
#include <vector>
#include <string>
#include <random>
#include "bcr.h"
#include <boost/math/distributions/inverse_gamma.hpp>


///0. Fixes
extern int global_seed;
extern string base_Path;
extern std::map<std::string, std::vector<float>> Seq_affinity;
extern std::map<std::string, std::map<std::string,std::map<int, std::map<std::string, std::string>>>> Restricted_db;
extern std::map<std::vector<std::string>, std::vector<float>> CDRs_affinity;
extern std::map<std::vector<std::string>, std::string> FWRs_db;
extern std::map<std::string, std::vector<float>> Seq_affinity_original;
extern std::map<std::string, std::map<std::string,std::map<int, std::map<std::string, std::string>>>> Restricted_db_original;
extern std::map<std::vector<std::string>, std::vector<float>> CDRs_affinity_original;
extern std::map<std::vector<std::string>, std::string> FWRs_db_original;



void confirmBCR(vector<float> &BCReceptor, std::string sequence, std::map<std::string, std::string> AA_REGIONS, int mother_ID, int cell_ID);
void setGermBCR(vector<float> &BCReceptor, vector<float> &Germ_BCR,  int mother_ID);
///1. Load FASTA files

static std::string fasta_path = "/home/rgarcia/Escritorio/NGly_scripts/Control_VJ_GC_model"; ///RRR5k "/home/rgarcia/Escritorio/NGly_scripts/MariaGC_1000_processed_fastas/Fastas_ABM/Fastas_r"
static std::vector<std::string> fastas;
static std::vector<std::string> fastas_remaining;
void findFastas(std::string path);
bool checkDNA(std::string DNA);
void findAndReplaceAll(std::string & data, std::string toSearch, std::string replaceStr);
void getSequence(std::map<std::string, std::string> & REGIONS, std::map<std::string, std::vector<int>> & MUTATIONS, std::string &sequence,  std::string &germline_name, int &ORF);
void add_seq_to_Restricted (std::map<std::string, std::string> & AA_REGIONS, std::string sequence_name);
std::string mapToString(std::map<std::string, std::vector<int>>  &m);
std::string mapToString2(std::map<std::string, std::string>  &m);
///2. Affinities
///
///
std::string DNAtoprotein (std::string newSequence);
std::map<std::string, std::string> Separating(std::map<std::string, vector<int>> MUTATIONS, std::string protein_sequence);
std::vector <int> Blacklisting(std::map<std::string, vector<int>> MUTATIONS);
std::vector<int> NGly_positions (std::string sequence, std::vector<int> blacklist, std::vector<int> & FULL_NGly_sites);
bool isThereAnStop (std::string sequence);

static double lambda = 0.4; ///1.5 for acpa


static double alpha_FWR1 = 0.17;
static double alpha_FWR2 = 0.17;
static double alpha_FWR3 = 0.29;
static double alpha_FWR4 = 0.05;
static double alpha_CDR1 = 0.10;
static double alpha_CDR2 = 0.03;
static double alpha_CDR3 = 0.19;


static double alpha_FWR = alpha_FWR1 + alpha_FWR2 + alpha_FWR3 +alpha_FWR4;


static double beta_FWR1 = 0.73;
static double beta_FWR2 = 0.67;
static double beta_FWR3 = 0.73;
static double beta_FWR4  = 0.60;
static double gamma_CDR1 = 0.79;
static double gamma_CDR2 = 0.76;
static double gamma_CDR3 = 0.75;


static double delta = 0.5; //0.18


static std::map<std::string, double> RATES_OF_MUTATIONS = {
    { "lambda",lambda },
    { "alpha_FWR", alpha_FWR },
    { "alpha_FWR1", alpha_FWR1 },
    { "alpha_FWR2", alpha_FWR2 },
    { "alpha_FWR3", alpha_FWR3 },
    { "alpha_FWR4", alpha_FWR4 },
    { "alpha_CDR1", alpha_CDR1 },
    { "alpha_CDR2", alpha_CDR2 },
    { "alpha_CDR3", alpha_CDR3 },
    { "beta_FWR1", beta_FWR1 },
    { "beta_FWR2", beta_FWR2 },
    { "beta_FWR3", beta_FWR3 },
    { "beta_FWR4", beta_FWR4 },
    { "gamma_CDR1", gamma_CDR1 },
    { "gamma_CDR2", gamma_CDR2 },
    { "gamma_CDR3", gamma_CDR3 },
    { "delta", delta }
};
std::vector<int> whichCellsToMutate(int new_n_cells, double lambda);


void randomMutWhere(int num_of_mutations, std::map<std::string, std::string> & REGIONS, std::map<std::string, std::vector<int>> & MUTATIONS, std::map<std::string, std::vector<int>> & SILENT_MUTATIONS, std::string & sequence,std::string & protein_sequence, std::map<std::string, double> RATES_OF_MUTATIONS, std::string & sequence_name);

///
///
///
#endif

