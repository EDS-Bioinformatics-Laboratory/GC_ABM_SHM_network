#ifndef BCR_H
#define BCR_H
#include <vector>
#include "parameters.h"
#include <string>
#include <map>
void initialize_Seeds (parameters &p);

///RRR
extern int global_seed;
/// RRR

class BCR
{
    
    //This class includes the implementation of B-cell receptor
public:
    
    //Shape-space method
    BCR(parameters &p);
    vector<float> BCReceptor;
    vector<float> Germ_BCReceptor;
    double pMut; // Probability of mutation
    int    nMutFromGermline; // Number of mutations from the beginning founder cell
    long double mutateBCR(parameters &p, int Delta_Mut_number, int Delta_NGly_number, std::map<std::string, std::string> AA_REGIONS,std::map<std::string, std::string> AA_GERM_REGIONS, std::map<std::string, std::string> AA_MOM_REGIONS, std::string sequence_name);
    double getMyAffinity4Ag(parameters& p);
    string print_BCR();

};

#endif // BCR_H
