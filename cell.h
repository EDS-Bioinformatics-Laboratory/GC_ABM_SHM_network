#ifndef CELL_H
#define CELL_H
#include <vector>
#include <string>
#include "vector3d.h"
#include "parameters.h"
#include "network.h"
#include "bcr.h"
#include <utility>
#include <sstream>
//The number Pi
 #define PI 3.141592654

using namespace std;

///R-SHM
#include "SHM.h"
///R-SHM
///
//Declearations
class lattice;
class T_cell;
class output;
class simulation;

//Cell types
enum celltype {empty,FDCell, Stromalcell, TFHC, Centroblast, Centrocyte, Plasmacell,Memorycell, border, cell_type_counter };

string cellToString(celltype v);

//Cell cycles
enum CellCycleState {cycle_G1, cycle_S, cycle_G2, cycle_M, cycle_Divide, cycle_G0, cycle_Ncellstates};

//Cell states
enum Cellstate {founder,unselected,contact_FDC,FDC_selected,contact_TC,TC_selected,recycled,apoptosis,TC_free,TC_connected,Plasma_Out,Plasma_in_GC,Memory_Out,Memory_in_GC,cell_state_counter};

class cell
{
public:
    cell(); //cell with new ID
    cell(cell* copied_cell);    //Constructs a new cell and copies fields from another cell to it
    virtual ~cell(){}   //Deconstructor
    int ID;     // ID of current cell
    int MID;    // ID of the mother cell
    Cellstate cell_state; //Status of the cell
    celltype cell_type; //Cell type
    vector3D position; // Current position
    vector3D polarity; // Current direction of movement
    double persistence_time; // Time left for next turn
    double speed;   //Speed of cell
    stringstream event; //records info of cell
    bool can_move; //A switch to turn moving on/off
    CellCycleState cyclestate; //#Recheck @danial: this is here only for redo function
    
    void getNewPersistentTime(parameters& p); //Update time left to calculate next polarity
    void getRandomPolarity(parameters& p, lattice& l); //Get a random polarity
    virtual void getNewPolarity(parameters& p, lattice& l); //Get a new polarity (not random)
    string print_neighbours(lattice &l); //Prints neighbours of the cell on current lattice l
    void move(parameters &p,  lattice& l , vector<vector3D> &redo_list); //move
    string printcell(); //Prints cell info
    double set_speed(parameters &p); // To change cell speed
};

class B_cell: public cell {
public:

    B_cell(parameters& p); // Get new B_cell with random BCR.
    B_cell(parameters& p, B_cell* Mom_cell); // Copy everything from mother cellto daughter.
    ~B_cell(){}
    BCR myBCR;
    int total_number_of_divisions; //#Recheck @danial: not neccessary
    int nFDCcontacts; // number of FDC contacts.
    double MyAffinity ;
    double pMHC_dependent_number_of_divisions;
    double cycle_state_time; //Time that has passed in the current cycle state
    double time_of_cycle_state_switch;   //Total time of current cycle state that cell has to pass to go to the next cycle
    double Bc_Tc_interaction_clock;
    double Recycling_delay; // Time it takes for B cell to recycle after getting selected by T cells
    double BC_FDC_interaction_clock; // Time since a B_cell becomes CC_free (in sec).
    double TC_selected_clock; // Time since a B_cell becomes selected by a Tcell
    double clock; //#Recheck @danial:change name of this later, Time since LAST interaction (in sec) with FDC. ((For refractory interaction time).
    //    double FDCinteractiontime; //Time since CC became in contact Ag. Swich on when using network
    double retained_Ag; // Ag internalized from interaction with FDC.
    bool IamHighAg; //Recycled cell will become output
    bool   Selected_by_FDC; //Indicates if B-cell rescued by FDC
    bool   Selected_by_TC;
    T_cell* interactingTC; // ID of T cell with who Bcell is interacting.
    CellCycleState cyclestate;
    int nDivisions2do; //Number of divisions left for cell to do
    double delta_Affinity;
    double TCsignalDuration; // Acumulated signal from currently interacting TC (in sec).(As imput to ODE)
    double fdc_interaction_time_history;
    pair<double ,double> Tc_interaction_history;
    bool isResponsive2CXCL12;
    bool isResponsive2CXCL13;
    
    void setMyAffinity(parameters &p);
    void transmit_CCdelay2cycle(parameters &p);
    void ContinueCellCycle(parameters &p);
    void clockreset();//reset clocks if CCrecycled but not differentiated to output
    void timeleft2recycle(parameters & p); //When finished dividing CBs calculate a remining time to differentiate to CC_free.
    void set_Retained_Ag(parameters& p); //#Recheck Update the nFDCcontacts and retained Ag of the B_cell before differentiating (If Ag from previous round should be deleted).
    long double mutate(parameters& p); //Mutation happens in BCR
    void proliferate(parameters &p,lattice &l,double time,vector<B_cell*> &ListB_cell, output &currentoutput, simulation &curent_sim);
    void Resensitize2Chemokines(parameters& p, lattice& l);
    void getNewPolarity(parameters& p, lattice& l);
    bool IsCycling() ;

    string printBcell(); //#Recheck @danial: add all fields to it

    ///R-SHM
    double Affinity_Seq;
    std::map<string, vector<int>> MUTATIONS;
    std::map<string, vector<int>> MOM_MUTATIONS;
    std::map<string, vector<int>> SILENT_MUTATIONS;
    std::map<string, vector<int>> MOM_SILENT_MUTATIONS;
    std::map<std::string, std::string> REGIONS;
    std::map<std::string, std::string> MOM_REGIONS;
    std::map<std::string, std::string> AA_REGIONS;
    std::map<std::string, std::string> AA_MOM_REGIONS;
    std::map<std::string, std::string> AA_GERM_REGIONS;
    std::string sequence;
    std::string sequence_origen;
    std::string protein_sequence;
    int ORF;
    vector<int> All_NGly_germline_sites;
    vector<int> NGly_germline_sites;
    vector<int> NGly_mother_sites;
    vector<int> All_NGly_mutated_sites;
    vector<int> NGly_mutated_sites;
    vector<int> blacklist_nums;
    std::string germline_name;

    ///R-SHM
    ///Elena-network
    int nRecyclings;
    network Bcell_network; // Every B-cell has a copy of the network structure
    double BCL6;
    double IRF4;
    double BLIMP1;
    double BCL6_0; ///RRR
    double IRF4_0; ///RRR
    double BLIMP1_0; ///RRR
    void setBcellTFs();
    void calcNetwork(double integraction_dt, double bcr, double cd40);
    bool TC_signal_start;
    ///Elena-network
};

void redo_move(vector <vector3D> &redo_list,lattice &l);

class T_cell: public cell {
public:
    T_cell(parameters& p);
    ~T_cell(){}
    Cellstate cell_state;
    vector<B_cell*> interactingCC;
    int nIncontactCCs;
    void liberateCC_TC(B_cell *bc);
    string printTcell();
    void getNewPolarity(parameters& p, lattice& l);

};

class FDC: public cell {
public:
     FDC() : cell() {volume = 0; AgperDendrite = 0;}
      ~FDC(){}
    vector<vector3D> occupiedPositions;
    int volume;
    double AgperDendrite;
    bool can_move;
};

class Stromal_cell: public cell {
public:
    Stromal_cell() : cell() {}
    bool can_move;
     ~Stromal_cell(){}
};

class Memory_cell: public cell {
public:
    Memory_cell(parameters & p);
    Memory_cell(parameters & p , B_cell* Bcell);
     ~Memory_cell(){}
    BCR myBCR;
    double MyAffinity ;
    bool can_move;
    bool   Selected_by_FDC; //Indicates if B-cell rescued by FDC
    bool   Selected_by_TC;
    int total_number_of_divisions;
    double  birth_time;
    double fdc_interaction_time_history;
    int nFDCcontacts; // number of FDC contacts.
    pair<double ,double> Tc_interaction_history;
    double retained_Ag;
    double delta_Affinity;
    bool isResponsive2CXCL12;
    bool isResponsive2CXCL13;
    void getNewPolarity(parameters& p, lattice& l);
    ///RRR
    std::string sequence;
    std::string protein_sequence;
    int ORF;
    std::map<string, vector<int>> MUTATIONS;
    std::map<string, vector<int>> SILENT_MUTATIONS;
    std::map<std::string, std::string> REGIONS;
    std::map<std::string, std::string> AA_REGIONS;
    std::map<std::string, std::string> AA_MOM_REGIONS;
    std::map<std::string, std::string> AA_GERM_REGIONS;
    vector<int> NGly_mutated_sites;
    vector<int> All_NGly_mutated_sites;
    std::string germline_name;
    vector<int> blacklist_nums;

    ///RRR
    ///Elena-network
    double BCL6;
    double IRF4;
    double BLIMP1;
    double BCL6_0; ///RRR
    double IRF4_0; ///RRR
    double BLIMP1_0; ///RRR
    ///Elena-network
};

class Plasma_cell: public cell {
public:
    Plasma_cell(parameters & p);
    Plasma_cell(parameters & p , B_cell* Bcell);
     ~Plasma_cell(){}
    BCR myBCR;
    bool   Selected_by_FDC; //Indicates if B-cell rescued by FDC
    bool   Selected_by_TC;
    double MyAffinity ;
    bool can_move;
    int total_number_of_divisions;
    double  birth_time;
    double fdc_interaction_time_history;
    int nFDCcontacts; // number of FDC contacts.
    pair<double ,double> Tc_interaction_history;
    double retained_Ag;
    double delta_Affinity;
    bool isResponsive2CXCL12;
    bool isResponsive2CXCL13;
    void getNewPolarity(parameters& p, lattice& l);
    ///RRR
    std::string sequence;
    std::string protein_sequence;
    int ORF;
    std::map<string, vector<int>> MUTATIONS;
    std::map<string, vector<int>> SILENT_MUTATIONS;
    std::map<std::string, std::string> REGIONS;
    std::map<std::string, std::string> AA_REGIONS;
    std::map<std::string, std::string> AA_MOM_REGIONS;
    std::map<std::string, std::string> AA_GERM_REGIONS;
    vector<int> NGly_mutated_sites;
    vector<int> All_NGly_mutated_sites;
    std::string germline_name;
    vector<int> blacklist_nums;

    ///RRR

    ///Elena-network
    double BCL6;
    double IRF4;
    double BLIMP1;

    double BCL6_0; ///RRR
    double IRF4_0; ///RRR
    double BLIMP1_0; ///RRR
    ///Elena-network
};


#endif // CELL_H


