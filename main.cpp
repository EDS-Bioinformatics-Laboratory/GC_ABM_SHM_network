#include "mafalda.h"
#include "random.h"
#include "GC3D.h"
#include <iostream>
#include <vector>
#include <string>
#include <getopt.h>

//#include <sys/stat.h>
using namespace std;

string base_Path=string("/home/rgarcia/Escritorio/NGly_network");
string outputFolder=base_Path + "/Output"; /// I know it should be inputted as command line, don't hate me
int global_seed=12531238;

#include <time.h>
#ifdef _WIN32
#include <windows.h>
#endif
#ifdef __linux__
#include <sys/stat.h>
#endif
#ifdef __APPLE__
#include <sys/param.h>
//#include <boost/filesystem.hpp>
#endif


int main(int argc, char** argv){

    cout << "Starting Mafalda " << endl;
    cout << "------------------------------------- How to use: --------------------------------------------------------\n";
    cout << "Mafalda mafaldaParameterFile.par \n";
    cout << "Mafalda hyphasmaparameterfile.par -h\n";
    cout << "  ... additional options that can be used:\n";
    cout << "  -s seedNumber \n";
    cout << "  -o outputFolder    or -o auto    to create and output in a folder with the parameter file name\n";
    cout << "----------------------------------------------------------------------------------------------------------\n";
    srand(global_seed); ////1233244567887 1233244567885
//    srand(time(NULL));
// Parsing arguments given from command line
// See https://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Option-Example.html#Getopt-Long-Option-Example

    string parfname = string();

// Using hyohasma parameter file
    int takeHyphasmaFile = false;
// Using specific seed for random number generator
    int requested_seed = -1;
    
//Numerb of simulations
    short Nof_simulations = 1;
    
// will be the remaining number of arguments during parsing
    int c = 0;
    while (c != -1){
        // Definition of the list of possible arguments: either "a_word" or "-X argument".
        // Other arguments (not from the predefined list), like parameter files, will also be retrieved at the end.
        static struct option long_options[] = {
            {"hypster",     no_argument,       0, 'h'},   // These options don’t set a flag.
            {"seed",    required_argument, 0, 's'},
            {"outputFolder",    required_argument, 0, 'o'},
            {0, 0, 0, 0}
        };
        
// Parses the arguments one by one
// The index of the identified argument will be put inside.
        int option_index = 0;
// Don't forget to reput the list of -x allowed inside the third argument of getopt
        c = getopt_long (argc, argv, "hs:o:", long_options, &option_index);  // "as:" means, expects -a without argument or -s with argument
        switch (c)
        {
        case 0:
            {if (long_options[option_index].flag != 0)  // i.e. if there was a flag associated. Nothing to do
                cout << long_options[option_index].name << endl;
            if (optarg)
                cout << " with arg " << optarg << endl;
                break;}
        case 'h':
            {cout << "   -h detected -> will read the parameter file as hyphasma file" << endl;
            takeHyphasmaFile = true;
                break;}
        case 's':
            { requested_seed = atoi(optarg); //turn string into int
            cout << "   -s detecetd -> using seed: " << requested_seed << endl;
//printf ("option -s with value `%s'\n", optarg);
                break;}
        case 'o':
            {outputFolder = optarg; //turn string into int
                
            cout << "   -o detected -> Using output folder: " << outputFolder << endl;
//printf ("option -s with value `%s'\n", optarg);
                break;}
                
            case 'i':
            {   Nof_simulations = atoi(optarg); //
                
                cout << "   -i detected -> Number of simulations: " << Nof_simulations << endl;
                //printf ("option -s with value `%s'\n", optarg);
                break;}
                
        case '?':
            { // getopt_long should have printed an error message.
                break;}
        default:
            {// no arguments within the list. Other arguments might be given (see next loop)
            //cout << "no argument given" << endl;
                break;}
        }
    }
    // Additional arguments, that are not in the 'official' list options (the parameter file for instance).
    int nArguments = argc - optind;
    if (optind < argc)
    {
        parfname = string(argv[optind]);
        if(argc > optind+1) cout << "Detected additional parameters: ";
        for(int i = optind+1; i < argc; ++i) // the first one is the parameter file, no need to print here
        {
            cout << "\t" << argv[i];
        }
        cout << endl;
    }
    if(nArguments != 1)
    {
        cerr << "You should give one (and only one) parameter file when you run Mafalda" << endl;
    }

    if((parfname.size() > 0) && (!parfname.substr(parfname.size()-4,4).compare(string(".par")))){
        parfname = parfname.substr(0, parfname.size()-4);
        cout << "   ... The parameter file contained .par at the end => cutted it into " << parfname << endl;
    }

// If the output path is not decleared or set as automatic it sets the output folder name as the parameter file name
    if((outputFolder.size() == 0) || (!outputFolder.compare(string("auto")))){
        outputFolder = parfname;
    }
    

// Print all parameters that is used to analyze file
// Create simulation instance
    if (Nof_simulations>1)
    {
        
    }
    else {
        parameters currentParameterSet;//(takeHyphasmaFile,parfname);
        currentParameterSet.convert_parameters();

        simulation Sim(currentParameterSet);
        Sim.Simulation_ID=Nof_simulations;
        ///initGC3D(argc, argv);  ///////RRRR ver que pasa
        Sim.simulate(*Sim.currentLattice, currentParameterSet);
        currentParameterSet.writeparameters(outputFolder+"/params.txt");
    }
    cout << "Starting OpenGL " << endl;
    
    //#check
    ofstream analyze;
    analyze.open(outputFolder+"/ana_ini.out");
    
// Initialize the visualization, later should put all the visualization together.
    analyze.close();
    return 0;

}
