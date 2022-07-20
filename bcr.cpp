#include "bcr.h"
#include "random.h"
#include <cmath>
#include <sstream>
#include <fstream>///R-SHM
#include "cell.h" ///R-SHM
#include "SHM.h" ///R-SHM
#include <iostream>///R-SHM

using namespace std;
vector <vector <int>> Shape_Space_Seeder;


std::default_random_engine generator(global_seed +11213);///R-SHM

// Constructor of BCR for the shape-space method
BCR::BCR(parameters& p)
{   //BCR constructor, sets the number of mutations from germline to 0, takes the size of BCR from parameter file (default: 4) and sets the BCR receptor using a random seed from a pool of randomly produced seeds (default: 100).
    nMutFromGermline = 0.;
    BCReceptor.resize(p.par[BCR_Length], 0);

}


//#Recheck danial: This can get improved
long double BCR::mutateBCR(parameters& p, int Delta_Mut_number, int Delta_NGly_number, std::map<std::string, std::string> AA_REGIONS, std::map<std::string, std::string> AA_GERM_REGIONS, std::map<std::string, std::string> AA_MOM_REGIONS,  std::string germline_name)
{
    long double aff1= getMyAffinity4Ag(p);
    if(BCReceptor.size() != 4) cerr << "Error:MutateBCR: BCR size is: "<< BCReceptor.size()<< endl;

    bool Lethal=false;
    ///R-SHM if(random::randomDouble(1) < pMut)

    if(isThereAnStop(AA_REGIONS["FWR1"]+ AA_REGIONS["CDR1"]+ AA_REGIONS["FWR4"]+ AA_REGIONS["CDR2"]+ AA_REGIONS["FWR3"]+ AA_REGIONS["CDR3"]+ AA_REGIONS["FWR4"])) {
        Lethal=true;
    }
    vector<string> regiones_F={"FWR1","FWR2", "FWR3","FWR4"};

    for (unsigned int zi = 0; zi < 4; zi++) {

        for (unsigned int i=0; i < AA_REGIONS[regiones_F[zi]].size(); i++) {
            if(Restricted_db[germline_name][regiones_F[zi]][i][string(1,AA_REGIONS[regiones_F[zi]][i])].compare("Lethal") == 0){
                Lethal=true;
            };
        }
    }

    if(FWRs_db[{"FWR1", AA_REGIONS["FWR1"]}]=="L" || FWRs_db[{"FWR2", AA_REGIONS["FWR2"]}]=="L" || FWRs_db[{"FWR3", AA_REGIONS["FWR3"]}]=="L" || FWRs_db[{"FWR4", AA_REGIONS["FWR4"]}]=="L" ){
        Lethal=true;
         //cout << "Lethal FWR" <<endl;
    };

    int diferencias=0;
    vector<string> regiones_C={"CDR1","CDR2", "CDR3"};
    for (unsigned int zi = 0; zi < 3; zi++) {

        for (unsigned int i=0; i < AA_REGIONS[regiones_C[zi]].size(); i++) {
            if(AA_REGIONS[regiones_C[zi]][i] != AA_GERM_REGIONS[regiones_C[zi]][i]){
                diferencias++;
            };
        }
    }

    int diferencias_mom=0;
    for (unsigned int zi = 0; zi < 3; zi++) {

        for (unsigned int i=0; i < AA_MOM_REGIONS[regiones_C[zi]].size(); i++) {
            if(AA_MOM_REGIONS[regiones_C[zi]][i] != AA_GERM_REGIONS[regiones_C[zi]][i]){
                diferencias_mom++;
            };
        }
    }

    int diferencias_daughter_mom=0;
    for (unsigned int zi = 0; zi < 3; zi++) {

        for (unsigned int i=0; i < AA_MOM_REGIONS[regiones_C[zi]].size(); i++) {
            if(AA_MOM_REGIONS[regiones_C[zi]][i] != AA_REGIONS[regiones_C[zi]][i]){
                diferencias_daughter_mom++;
//                if (Delta_Mut_number == 3) {
//                   cout << "Something wrong with Deltas2: " << AA_MOM_REGIONS[regiones_C[zi]][i]<< " vs " << AA_REGIONS[regiones_C[zi]][i]<< endl;
//                   cout << AA_MOM_REGIONS[regiones_C[zi]] <<  " vs " << AA_REGIONS[regiones_C[zi]]<< endl;
//                }
            };
        }
    }


    if (diferencias_daughter_mom != Delta_Mut_number) {
        cout << "Something wrong with Deltas: " << Delta_Mut_number<< " vs " << diferencias_daughter_mom<< endl;
        cout <<AA_MOM_REGIONS[regiones_C[0]]<< " "  <<AA_MOM_REGIONS[regiones_C[1]]<< " " <<AA_MOM_REGIONS[regiones_C[2]]<< endl;
        cout <<AA_REGIONS[regiones_C[0]]<< " "  <<AA_REGIONS[regiones_C[1]]<< " " <<AA_REGIONS[regiones_C[2]]<< endl;
    }


    vector<float> tmpBCR= Germ_BCReceptor;
    if (!Lethal) {
        double step=0;
        double max_step=1;
        short int randomL; //random location
        short int randomC; //random change
        short int randomNo; //random NO change
        double aff;

        std::random_device rd;

        std::normal_distribution<float> distribution1(1,0.1);
        std::normal_distribution<float> distribution0(1,0.1);





        if (Delta_Mut_number > 0) {
             int num_steps;
             if(diferencias == (diferencias_mom + Delta_Mut_number)) {
                num_steps = Delta_Mut_number;
             } else {
//                 cout << diferencias <<  " vs " << diferencias_mom<< " + " << Delta_Mut_number << " + " << AA_REGIONS["CDR1"] << " + " << AA_REGIONS["CDR2"]<< " + " << AA_REGIONS["CDR3"]<< endl;
                num_steps = diferencias;
                BCReceptor=Germ_BCReceptor;
                nMutFromGermline = 0;
                //if (diferencias_mom != diferencias) {
                //    cout <<"Special case!"<<endl;
                //    cout << "diferencias_mom: " << diferencias_mom<< endl;
                //    cout << "diferencias_germ: " << diferencias<< endl;
                //    cout << "Delta_Mut_number: " << Delta_Mut_number<< endl;
                //    cout << "Delta_Mut_number: " << diferencias_daughter_mom<< endl;
                //    cout << AA_REGIONS["CDR1"] << AA_REGIONS["CDR2"] << AA_REGIONS["CDR3"] << endl;
                //    cout << AA_MOM_REGIONS["CDR1"] << AA_MOM_REGIONS["CDR2"] << AA_MOM_REGIONS["CDR3"] << endl;
                //    cout << AA_GERM_REGIONS["CDR1"] << AA_GERM_REGIONS["CDR2"] << AA_GERM_REGIONS["CDR3"] << endl;
                //}
             }

             for( int a = 0; a < num_steps; a = a + 1 ) {
                 //cout <<"activado"<<endl;
                 aff = double(int(getMyAffinity4Ag(p) * 10))/ 10;


                 //+0.5 anadido por mi     //#temporary Danial: To make the probability of mutation independent of the affinity and keeping the average affinity of output cells still high, it is possible to change the mutation step size
                 nMutFromGermline += 1; //Increase numer of mutations from the beginning founder
                 randomL = random::randomInteger(4);
                 randomC = random::randomInteger(2);
                 randomNo = random::randomInteger(20);
                 step=0;

         ////hacer que maximo sea 9, minimo 0, y repetir jugada si tal

                 while ((isgreaterequal(BCReceptor[randomL], 9.0) && randomC==1)|| (islessequal(BCReceptor[randomL], 0.0) && randomC==0))
                 {
                     randomL = random::randomInteger(4);
                     randomC = random::randomInteger(2);
                     randomNo = random::randomInteger(20);
                 }
                 if (randomC==1)
                 {
                     while((islessequal(step, 0.0))) {
                           step=distribution1(generator);
                           if (randomNo >= 23) {
                               step=distribution0(generator);}

                     }

                     if(isgreaterequal(BCReceptor[randomL]+step, 9.0)) {
                         step=9.0-BCReceptor[randomL];
                     }

                     BCReceptor[randomL] += step;

                    }
                 if (randomC==0)
                 {
                     while((islessequal(step, 0.0))) {
                           step=distribution1(generator);
                           if (randomNo >= 23) {
                               step=distribution0(generator);}

                     }
                     if(islessequal(BCReceptor[randomL]-step, 0.0)) {
                         step=BCReceptor[randomL]-0.0;
                     }

                          BCReceptor[randomL] -= step;

                 }

             }
        }




        ////////////////////////// NGLY
        ///
//        std::normal_distribution<float> distribution1n(1,0.15);
//        std::normal_distribution<float> distribution0n(1,0.15);
//        for( int b = 2000; b < Delta_NGly_number; b = b + 1 ) {
//            cout<<"NGly raro";
//        randomL = random::randomInteger(4);
//        randomC = random::randomInteger(2);
//        randomNo = random::randomInteger(20);


//////hacer que maximo sea 9, minimo 0, y repetir jugada si tal

//        while ((isgreaterequal(BCReceptor[randomL], 9.0) && randomC==1)|| (islessequal(BCReceptor[randomL], 0.0) && randomC==0))
//     {
//            randomL = random::randomInteger(4);
//            randomC = random::randomInteger(2);
//            randomNo = random::randomInteger(20);
//        }
//        if (randomC==1)
//        {
//            while((islessequal(step, 0.0))) {
//                  step=distribution1n(generator);
//                  if (randomNo >= 10) {
//                      step=distribution0n(generator);}

//            }

//            if(isgreaterequal(BCReceptor[randomL]+step, 9.0)) {
//                step=9.0-BCReceptor[randomL];
//            }
//            BCReceptor[randomL] += step;}
//        if (randomC==0)
//        {
//            while((islessequal(step, 0.0))) {
//                 step=distribution1n(generator);
//                  if (randomNo >= 10) {
//                      step=distribution0n(generator);}

//            }
//            if(islessequal(BCReceptor[randomL]-step, 0.0)) {
//                step=BCReceptor[randomL]-0.0;
//            }
//            BCReceptor[randomL] -= step;}
//        }
    } else {
        BCReceptor[0]=9.0;
        BCReceptor[1]=9.0;
        BCReceptor[2]=9.0;
        BCReceptor[3]=9.0;
    }




    long double aff2= getMyAffinity4Ag(p);
    return double(aff2-aff1); // Returns the change in affinity due to the mutation
}
//This function calculates affinity using shape-space concept by computing the distance between Ag and bcr
double BCR::getMyAffinity4Ag (parameters &p)
{
    double dist = 0.;
    for(unsigned int i = 0; i < BCReceptor.size(); i++)
    {
       dist += fabs(double (3. - BCReceptor[i])); // Philippe: now the target antigen is 3333, modify later
    }
    double dist2=dist*dist;
    return exp(-1. * dist2 / (2.8 * 2.8)); //Danial: #Recheck
}

string BCR::print_BCR(){
    stringstream res;
    int L = (int) BCReceptor.size();
    for(int i = 0; i < L; ++i)
    {
        res<<BCReceptor[i];
    }
    return res.str();
}




