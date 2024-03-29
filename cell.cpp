//#Recheck @danial: Use of header files has not been checked yet.
#include "cell.h"
#include "lattice.h"
#include "random.h"
#include <sstream>
#include "vector3d.h"
#include "vector"
#include "bcr.h"
#include "math.h"
#include "output.h"
#include "mafalda.h"

///RRR
#include <algorithm>
///RRR

using namespace std;

int getNewId() {
  static int cpt = -1;
  cpt++;
  return cpt;
}

///RRR

std::string cellToString(celltype v)
{
    switch (v)
    {
        case celltype::empty:   return "empty";
        case FDCell:   return "FDCell";
        case Stromalcell: return "Stromalcell";
        case TFHC:   return "TFHC";
        case Centroblast:   return "Centroblast";
        case Centrocyte: return "Centrocyte";
        case Plasmacell:   return "Plasmacell";
        case Memorycell:   return "Memorycell";
        case border: return "border";
        case cell_type_counter: return "cell_type_counter";

        default:      return "Unknown";
    }
}
///RRR


// This function creates a new cell with a new Id
cell::cell() : position(0, 0, 0), polarity(0., 0., 0.) {
  ID = getNewId(); // Get new ID for cell
  MID = -1;  // Mother ID
  cell_state = cell_state_counter; //Status of the cell being set to counter (this means not having an state yet)
  cell_type= cell_type_counter;    //Type of the cell being set to counter (this means not having a type yet)
  persistence_time = 0; // Time left for next turn
  speed = 0.0;  //Speed of cell
  can_move = false; //A switch to turn moving on/off
}

// Create new cell with new ID and copy another cell info to it including position and polarity
cell::cell(cell* copied_cell) : position(0, 0, 0), polarity(0., 0., 0.) {
  // Cell default
  ID = getNewId(); // Get new ID for cell
  MID = copied_cell->ID;    // Mother ID
  cell_state=copied_cell->cell_state; //Status of the cell
  cell_type=copied_cell->cell_type; //Type of the cell
  persistence_time = copied_cell->persistence_time; // Time left for next turn
  speed = copied_cell->speed;   //Speed
  can_move = copied_cell->can_move;
  position = copied_cell->position; //Current position
  polarity = copied_cell->polarity; //Current direction
    
}

//#Recheck danial: this can be done in a projection matrix.
void cell::getRandomPolarity(parameters& p, lattice& l) {
  double n[3];
  n[0] = polarity.X;
  n[1] = polarity.Y;
  n[2] = polarity.Z;

  // Sample theta from  distribution *** Define type in parameters*** and Phi
  // from random distribution [0, 360]
  double theta = l.thetas.get_distribution_value();
  double phi = random::randomDouble(2.0 * PI);
  double tmp[3], ttmp[3];
  // find the phi of the old polarity:
  double nphi = 0.;
  if ((n[0] == 0.) && (n[1] == 0.)) {
    nphi = 0.;
  } else {
    nphi = acos(n[0] / sqrt(n[0] * n[0] + n[1] * n[1]));
  }
  if (n[1] < 0) {
    nphi = 2.0 * PI - nphi;
  }

  // turn n onto the x-z-plane by rotation of -nphi around the z-axis:
  nphi *= -1.0;
  tmp[0] = cos(nphi) * n[0] - sin(nphi) * n[1];
  tmp[1] = sin(nphi) * n[0] + cos(nphi) * n[1];
  tmp[2] = n[2];

  // turn the vector in the z-plane by theta around the y-axis
  ttmp[0] = cos(theta) * tmp[0] + sin(theta) * tmp[2];
  ttmp[1] = tmp[1];
  ttmp[2] = -1.0 * sin(theta) * tmp[0] + cos(theta) * tmp[2];
  // turn back by nphi around z-axis
  nphi *= -1.0;
  tmp[0] = cos(nphi) * ttmp[0] - sin(nphi) * ttmp[1];
  tmp[1] = sin(nphi) * ttmp[0] + cos(nphi) * ttmp[1];
  tmp[2] = ttmp[2];
  // turn the new vector tmp by phi around the old vector n
  ttmp[0] = (cos(phi) + n[0] * n[0] * (1.0 - cos(phi))) * tmp[0] +
            (n[0] * n[1] * (1.0 - cos(phi)) - n[2] * sin(phi)) * tmp[1] +
            (n[0] * n[2] * (1.0 - cos(phi)) + n[1] * sin(phi)) * tmp[2];
  ttmp[1] = (n[0] * n[1] * (1.0 - cos(phi)) + n[2] * sin(phi)) * tmp[0] +
            (cos(phi) + n[1] * n[1] * (1.0 - cos(phi))) * tmp[1] +
            (n[1] * n[2] * (1.0 - cos(phi)) - n[0] * sin(phi)) * tmp[2];
  ttmp[2] = (n[0] * n[2] * (1.0 - cos(phi)) - n[1] * sin(phi)) * tmp[0] +
            (n[1] * n[2] * (1.0 - cos(phi)) + n[0] * sin(phi)) * tmp[1] +
            (cos(phi) + n[2] * n[2] * (1.0 - cos(phi))) * tmp[2];

  double newnorm =
      sqrt(ttmp[0] * ttmp[0] + ttmp[1] * ttmp[1] + ttmp[2] * ttmp[2]);

  for (short a = 0; a < 3; a++) {
    n[a] = ttmp[a] / newnorm;
  }
  polarity.X = n[0];
  polarity.Y = n[1];
  polarity.Z = n[2];
  if (isnan(polarity.X) || isnan(polarity.Y) || isnan(polarity.Z)) {
    cout << "Error in random polarity function" << endl;
    exit(1);
  }
}


// This function calculates a new polarity for cells, default is random but for different class of cells it is affected by chemokines or other factors
void cell::getNewPolarity(parameters& p, lattice& l) {
  
    getRandomPolarity(p, l);  // Set polarity to a random vector.
}
//This function moves the cell on lattice. Since the moving mechanism is the same for different types of cells, there is only one function.
void cell::move(parameters& p, lattice& l, vector<vector3D>& redo_list) {
  
    if (can_move) {
    if (random::randomDouble(1.0) < persistence_time) {
      getNewPolarity(p, l);     // Update polarity
      getNewPersistentTime(p);  // Update persistence time
    }
    if (random::randomDouble(1) < speed) {
      double maxprojection = -99;
      long takethisindex = -1;
      vector3D diff(-1, -1, -1);
      vector<vector3D> neighbours;
      neighbours.reserve(6);
      neighbours = l.getNeighbour_nn(position);
      for (unsigned int i = 0; i < neighbours.size(); i++) {
        if ((l.insideBorders(neighbours[i])) &&
            (l.celltypeat(neighbours[i]) == celltype::empty))  ///RRR
        {
          diff.X = (neighbours[i].X - position.X);
          diff.Y = (neighbours[i].Y - position.Y);
          diff.Z = (neighbours[i].Z - position.Z);
          double scals = getScalarproduct(diff, polarity);
          if (scals > maxprojection) {
            maxprojection = scals;
            takethisindex = i;
          }
        }
      }
      if ((takethisindex >= 0) && (maxprojection >= 0)) {
        // do the movement only if the scalarproduct is positive
        // and an empty neighbour was found:
        l.removecellat(position);
        position = neighbours[takethisindex];
        l.putcellat(this);
      } else {
        redo_list.push_back(position);
      }
    }
  }
}


string cell::printcell() {
  stringstream res;
    //#Recheck @danial: add all fields
    res <<"ID: "<<ID<<" Celltype: "<<cell_type<<" Cellstate: "<<cell_state<<" Position XYZ: "<<position.print()<<" Polarity XYZ: "<<polarity.print()<<endl;
  return res.str();
}

//////////////////////////////B-cell/////////////////////////


// Create a new B_cell with random BCR
B_cell::B_cell(parameters& p)
: cell() /*cell with new ID MID = -1 */, myBCR(p) {
   //ID, MID, can_move, speed, persistence time, cell_type and cell_state are set in the cell() constructor
    total_number_of_divisions = 0;
    nFDCcontacts = 0;
    setMyAffinity(p);
    pMHC_dependent_number_of_divisions = 0.0;
    cycle_state_time = 0.;
    time_of_cycle_state_switch = 0.;
    Bc_Tc_interaction_clock = 0.;
    Recycling_delay=0.0;
    BC_FDC_interaction_clock = 0.;  // Time since a B_cell became CC_free (in sec).
    TC_selected_clock=0.0;
    clock = 0.;         // Time since LAST interaction (in sec) with FDC. ((For
    retained_Ag = 0.;    // Ag internalized from interaction with FDC.
    IamHighAg = false;         // Recycled cell will become output
    Selected_by_FDC = false;   // Indicates if rescued by FDC
    Selected_by_TC = false;
    interactingTC = NULL;
    cyclestate = cycle_Ncellstates;
    nDivisions2do = 0.;
    delta_Affinity = 0.0;
    TCsignalDuration = 0.;  // Acumulated signal from currently interacting TC (in
    fdc_interaction_time_history = 0.0;
    Tc_interaction_history.first = 0.0;
    Tc_interaction_history.second = 0.0;
    isResponsive2CXCL12=false;
    isResponsive2CXCL13=false;


    ///R-SHM
    ///
    MUTATIONS = {
        { "FWR1", vector<int>(0) },
        { "FWR2", vector<int>(0) },
        { "FWR3", vector<int>(0) },
        { "FWR4", vector<int>(0) },
        { "CDR1", vector<int>(0) },
        { "CDR2", vector<int>(0) },
        { "CDR3", vector<int>(0) },
    };
    MOM_MUTATIONS = {
        { "FWR1", vector<int>(0) },
        { "FWR2", vector<int>(0) },
        { "FWR3", vector<int>(0) },
        { "FWR4", vector<int>(0) },
        { "CDR1", vector<int>(0) },
        { "CDR2", vector<int>(0) },
        { "CDR3", vector<int>(0) },
    };
    SILENT_MUTATIONS = {
        { "FWR1", vector<int>(0) },
        { "FWR2", vector<int>(0) },
        { "FWR3", vector<int>(0) },
        { "FWR4", vector<int>(0) },
        { "CDR1", vector<int>(0) },
        { "CDR2", vector<int>(0) },
        { "CDR3", vector<int>(0) },
    };
    MOM_SILENT_MUTATIONS = {
        { "FWR1", vector<int>(0) },
        { "FWR2", vector<int>(0) },
        { "FWR3", vector<int>(0) },
        { "FWR4", vector<int>(0) },
        { "CDR1", vector<int>(0) },
        { "CDR2", vector<int>(0) },
        { "CDR3", vector<int>(0) },
    };
    MOM_REGIONS = {
            { "FWR1", "" },
            { "FWR2", "" },
            { "FWR3", "" },
            { "FWR4", "" },
            { "CDR1", "" },
            { "CDR2", "" },
            { "CDR3", "" },
    };
    REGIONS = {
            { "FWR1", "" },
            { "FWR2", "" },
            { "FWR3", "" },
            { "FWR4", "" },
            { "CDR1", "" },
            { "CDR2", "" },
            { "CDR3", "" },
    };
    sequence="";
    sequence_origen="";
    germline_name="tofill";
    ORF=453;
    getSequence(REGIONS, MUTATIONS, sequence_origen, germline_name, ORF);
    SILENT_MUTATIONS=MUTATIONS;
    sequence=sequence_origen;
    protein_sequence = DNAtoprotein(sequence);
    AA_REGIONS=Separating(MUTATIONS, protein_sequence);
    AA_MOM_REGIONS=AA_REGIONS;
    AA_GERM_REGIONS=AA_REGIONS;
    confirmBCR(myBCR.BCReceptor,protein_sequence, AA_REGIONS, cell::MID,cell::ID);
    setGermBCR(myBCR.BCReceptor,myBCR.Germ_BCReceptor,cell::MID);
    add_seq_to_Restricted(AA_REGIONS, germline_name);
    All_NGly_germline_sites;
    All_NGly_mutated_sites;
    //std::cout << protein_sequence << std::endl;
    NGly_germline_sites= NGly_positions(protein_sequence, Blacklisting(MUTATIONS),  All_NGly_germline_sites);
    NGly_mother_sites= NGly_germline_sites;
    NGly_mutated_sites= NGly_germline_sites;
    All_NGly_mutated_sites= All_NGly_germline_sites;
    blacklist_nums = Blacklisting(MUTATIONS);

    ///R-SHM
    

    ///Elena-network
    Bcell_network.setBaseParameters(); //Elena: network: Set parameters in network using parameters defined inside network class (not from file).
    Bcell_network.initialise();//Elena: network: Set initial TF levels in init vector in network
    setBcellTFs(); //Elena: network: Puts TF levels from init vector (in network) inside Bcell TFs (field)
    BCL6_0=BCL6; ///RRR
    IRF4_0=IRF4; ///RRR
    BLIMP1_0=BLIMP1; ///RRR
    nRecyclings = 0;
    ///Elena-network
}

// Copy from mother to daughter cell

B_cell::B_cell(parameters& p, B_cell* Mom):cell(), myBCR(p) {
    //#Recheck
    MID = Mom->ID;
    Mom->ID = getNewId(); //This is changed becaue every division creates two new cells
    myBCR = Mom->myBCR;
    cell_state = Mom->cell_state;
    cell_type = Mom->cell_type;
    persistence_time = Mom->persistence_time;  // Danial: cahnge this
    speed = Mom->speed;
    can_move = Mom->can_move;
    nFDCcontacts = 0.;
    setMyAffinity(p);
    pMHC_dependent_number_of_divisions = Mom->pMHC_dependent_number_of_divisions; //#Recheck @danial: don't need this field
    total_number_of_divisions = Mom->total_number_of_divisions;  // #Recheck
    cyclestate = Mom->cyclestate; //#Recheck
    interactingTC = NULL;
    cycle_state_time = 0.;
    time_of_cycle_state_switch = 0.;
    Bc_Tc_interaction_clock = 0.;
    Recycling_delay =
    0.;  // Time spent inside CBgoingLZ (time for moving to Light Zone)
    BC_FDC_interaction_clock = 0.;  // Time since a B_cell became CC_free (in sec).
    clock = Mom->clock;
    TC_selected_clock=0.0;
    retained_Ag = 0.;             // Ag internalized from interaction with FDC.
    Selected_by_FDC = false;         // Indicates if rescued by FDC
    Selected_by_TC = false;
    delta_Affinity = 0.0;
    nDivisions2do = Mom->nDivisions2do;
    TCsignalDuration = 0.;  // Acumulated signal from currently interacting TC (in
    // sec).(As imput to ODE)
    IamHighAg = false;
    isResponsive2CXCL12 = Mom->isResponsive2CXCL12;
    isResponsive2CXCL13 = Mom->isResponsive2CXCL13;
    Tc_interaction_history.first = 0.0;
    Tc_interaction_history.second = 0.0;
    fdc_interaction_time_history = 0.0;

    ///R-SHM
    MOM_MUTATIONS =Mom->MUTATIONS;
    MUTATIONS =Mom->MUTATIONS;
    MOM_SILENT_MUTATIONS =Mom->SILENT_MUTATIONS;
    SILENT_MUTATIONS =Mom->SILENT_MUTATIONS;
    REGIONS =Mom->REGIONS;
    MOM_REGIONS =Mom->REGIONS;
    sequence=Mom->sequence;
    ORF=Mom->ORF;
    AA_GERM_REGIONS =Mom->AA_GERM_REGIONS;
    germline_name=Mom->germline_name;
    protein_sequence=Mom->protein_sequence;
    AA_REGIONS=Separating(MUTATIONS, protein_sequence);
    AA_MOM_REGIONS=Mom->AA_REGIONS;
    AA_GERM_REGIONS =Mom->AA_GERM_REGIONS;
    Affinity_Seq=Mom->Affinity_Seq;
    //cout << "HIJA" << endl;
    sequence_origen =Mom->sequence;
    NGly_germline_sites =Mom->NGly_germline_sites;
    All_NGly_germline_sites=Mom->All_NGly_germline_sites;
    NGly_mother_sites =Mom->NGly_mutated_sites;
    NGly_mutated_sites =Mom->NGly_mutated_sites;
    All_NGly_mutated_sites =Mom-> All_NGly_mutated_sites;
    blacklist_nums=Mom-> blacklist_nums;

    ///R-SHM
    ///
    ///Elena-network
    Bcell_network.setBaseParameters(); //Elena: network: Set parameters in network using parameters defined inside network class (not from file).
    Bcell_network.initialise();//Elena: network: Set initial TF levels in init vector in network
    setBcellTFs(); //Elena: network: Puts TF levels from init vector (in network) inside Bcell TFs (field)
    BCL6_0=BCL6; ///RRR
    IRF4_0=IRF4; ///RRR
    BLIMP1_0=BLIMP1; ///RRR
    nRecyclings = Mom->nRecyclings;
    ///Elena-network
}

///Elena-network
void B_cell::calcNetwork(double integraction_dt, double bcr, double cd40) //Elena: network: Initialize network parameters, simulate TF levels for next time step.
{
if(BCL6 ==NAN || BLIMP1==NAN || IRF4==NAN)
    cerr<<"Error: B_cell::calcNetwork: NAN TF levels"<<endl;
if(BCL6 < 0 || BLIMP1 < 0 || IRF4 < 0 )
    cerr<<"Error: B_cell::calcNetwork: Negative TF levels!"<<endl;
Bcell_network.setDinamicParameters(BLIMP1,BCL6,IRF4,bcr,cd40);//Add current TF levels and signal strength (on/off) as parameters to the network.
Bcell_network.initialise(); //Elena: Puts current Bcell TF levels in init vector (in network)
Bcell_network.simulate(integraction_dt); // Integrates next TFs levels using the dt given through argument
setBcellTFs(); //Put resulting TF values inside Bcell.
}

void B_cell::setBcellTFs() //Elena: network: Refresh cell TFs with network output.
{
BCL6 =  Bcell_network.val.at(0);
IRF4 =  Bcell_network.val.at(1);
BLIMP1 =  Bcell_network.val.at(2);;

if(BCL6 ==NAN || BLIMP1==NAN || IRF4==NAN)
    cerr<<"Error: B_cell::setBcellTFs: NAN TF levels"<<endl;
if(BCL6 < 0 || BLIMP1 < 0 || IRF4 < 0 )
    cerr<<"Error: network::setBcellTFs: Negative TF levels!"<<endl;
if(BCL6 > 1000 || BLIMP1 > 1000 || IRF4 > 1000)
    cerr<<"Error: B_cell::setBcellTFs: TF levels way above steady state !"<<endl;
}
///Elena-network

void B_cell::ContinueCellCycle(parameters& p) {
  switch (cyclestate) {
    case cycle_G1: {
      cyclestate = cycle_S;
      time_of_cycle_state_switch =
          random::cell_cycle_time(p.par[c_S], cycle_S);
      break;
    };
    case cycle_S: {
      cyclestate = cycle_G2;
      time_of_cycle_state_switch =
          random::cell_cycle_time(p.par[c_G2], cycle_G2);
      break;
    }
    case cycle_G2: {
      cyclestate = cycle_M;
      time_of_cycle_state_switch =
          random::cell_cycle_time(p.par[c_M], cycle_M);
      break;
    }
    case cycle_M: {
      cyclestate = cycle_Divide;
      break;
    }
    case cycle_Divide: {
      if (nDivisions2do <= 0) {
        cyclestate = cycle_G0;
        time_of_cycle_state_switch = 0;
      }
      break;
    }
    case cycle_G0: {
      break;
    }
    default: {
      cerr << "cell::ContinueCellCycle Error: Cell ID " << ID << " type "
           << cell_type << " in cell cycle state " << cyclestate << endl;
      break;
    }
  }
}



void B_cell::setMyAffinity(parameters& p) {
  MyAffinity = myBCR.getMyAffinity4Ag(p);

  if (isThereAnStop(protein_sequence) || FWRs_db[{"FWR1", AA_REGIONS["FWR1"]}]=="L"  || FWRs_db[{"FWR2", AA_REGIONS["FWR2"]}]=="L"  || FWRs_db[{"FWR3", AA_REGIONS["FWR3"]}]=="L"  || FWRs_db[{"FWR4", AA_REGIONS["FWR4"]}]=="L"){
        MyAffinity = 0;
        cell_state=apoptosis;
  }
}

long double B_cell::mutate(parameters& p) {
    double delta_Affinity = 0;
    ///int number_of_mutations=whichCellsToMutate(1, (RATES_OF_MUTATIONS["lambda"] + double ((0. - RATES_OF_MUTATIONS["lambda"]) * pow(MyAffinity,1))))[0];
    int number_of_mutations=whichCellsToMutate(1, RATES_OF_MUTATIONS["lambda"])[0];

    if (myBCR.BCReceptor[0] == 9 && myBCR.BCReceptor[1] == 9 && myBCR.BCReceptor[2] == 9 && myBCR.BCReceptor[3] == 9 ) {
        cout <<"A cell with 9,9,9,9 is dividing!" << endl;
    }
    if(number_of_mutations!=0) {
        long double aff1= myBCR.getMyAffinity4Ag(p);
        std::string mother_prot= protein_sequence;
        AA_MOM_REGIONS=AA_REGIONS;
        MOM_REGIONS=REGIONS;
        MOM_MUTATIONS = MUTATIONS;
        MOM_SILENT_MUTATIONS = SILENT_MUTATIONS;

        randomMutWhere(number_of_mutations, REGIONS, MUTATIONS, SILENT_MUTATIONS, sequence, protein_sequence, RATES_OF_MUTATIONS, germline_name);

        protein_sequence = DNAtoprotein(sequence);
        double testing;

        AA_REGIONS=Separating(MUTATIONS, protein_sequence);
        blacklist_nums=Blacklisting(MUTATIONS);

        NGly_mutated_sites= NGly_positions(protein_sequence, blacklist_nums, All_NGly_mutated_sites);

        if (Seq_affinity.find(protein_sequence) == Seq_affinity.end() && (CDRs_affinity.find({AA_REGIONS["CDR1"],AA_REGIONS["CDR2"],AA_REGIONS["CDR3"]}) == CDRs_affinity.end())) {


            int Delta_Num_Mut=0;

            for( int a = 0; a < mother_prot.length(); a = a + 1) {
                if (mother_prot[a] != protein_sequence[a] && std::find(blacklist_nums.begin(), blacklist_nums.end(), a) == blacklist_nums.end()){
                    Delta_Num_Mut=Delta_Num_Mut+1;
                    //cout <<"diff en AA " <<  a <<endl;
                }
            }


            int Delta_Num_NGly=0;
            std::vector<int> diff;
            //not need to sort since it already sorted
            std::set_difference(NGly_mutated_sites.begin(), NGly_mutated_sites.end(), NGly_mother_sites.begin(), NGly_mother_sites.end(),
                std::inserter(diff, diff.begin()));
            Delta_Num_NGly=diff.size();

            delta_Affinity=myBCR.mutateBCR(p,Delta_Num_Mut,Delta_Num_NGly, AA_REGIONS, AA_GERM_REGIONS,AA_MOM_REGIONS, germline_name);
            testing=delta_Affinity;

            confirmBCR(myBCR.BCReceptor,protein_sequence, AA_REGIONS, MID, ID);
        } else {

            confirmBCR(myBCR.BCReceptor,protein_sequence, AA_REGIONS, MID, ID);
            testing=myBCR.getMyAffinity4Ag(p)-aff1;
        }

        long double aff2= myBCR.getMyAffinity4Ag(p);
        delta_Affinity= aff2 - aff1;
        if (testing != delta_Affinity) {
            cout << "Diferencias entre " << delta_Affinity << "y" << testing << endl;
            cout << "Numero nuevas mutaciones " << number_of_mutations << endl;
            cout << AA_REGIONS["FWR1"] << " " << AA_REGIONS["CDR1"] << " " << AA_REGIONS["FWR2"] << " " << AA_REGIONS["CDR2"] << " " << AA_REGIONS["FWR3"] << " " << AA_REGIONS["CDR3"] << " " << AA_REGIONS["FWR4"] << " "<< endl;
            cout << AA_MOM_REGIONS["FWR1"] << " " << AA_MOM_REGIONS["CDR1"] << " " << AA_MOM_REGIONS["FWR2"] << " " << AA_MOM_REGIONS["CDR2"] << " " << AA_MOM_REGIONS["FWR3"] << " " << AA_MOM_REGIONS["CDR3"] << " " << AA_MOM_REGIONS["FWR4"] << " "<< endl;
            for (auto i = blacklist_nums.begin(); i != blacklist_nums.end(); ++i) {
                std::cout << *i << ' ';
            }
            std::cout << endl;
            if (Seq_affinity.find(protein_sequence)== Seq_affinity.end()) {
                cout <<"NO ESTÁ"<<endl;
            }
            cout << "Seq_db " << Seq_affinity[protein_sequence][0] << "," << Seq_affinity[protein_sequence][1] << ","<< Seq_affinity[protein_sequence][2] << ","<< Seq_affinity[protein_sequence][3] << endl;
            cout << "CDR_db " << CDRs_affinity[{AA_REGIONS["CDR1"],AA_REGIONS["CDR2"],AA_REGIONS["CDR3"]}][0] << "," << CDRs_affinity[{AA_REGIONS["CDR1"],AA_REGIONS["CDR2"],AA_REGIONS["CDR3"]}][1] << ","<< CDRs_affinity[{AA_REGIONS["CDR1"],AA_REGIONS["CDR2"],AA_REGIONS["CDR3"]}][2] << ","<< CDRs_affinity[{AA_REGIONS["CDR1"],AA_REGIONS["CDR2"],AA_REGIONS["CDR3"]}][3] << endl;

        }
    }


  return delta_Affinity;  // Mutates BCR
}

//#Recheck @danial: change by a distribution
// Set the rate of differentiation CB2CC and CC2CB depending on the cell type
void B_cell::timeleft2recycle(parameters& p) {
  switch (cell_type) {
    case Centrocyte: {
      if (Recycling_delay > 0) {
        cerr << " Error: B_cell::" << cell_state
             << " timeleft2recycle(): " << Recycling_delay
             << " CC cell did not finish its remaining time To Differentiate"
             << endl;
      }

      double shift = 0.;
      double sigmoi;
      double delaytmp = 6.;
      double width = delaytmp * 0.1;
      double delaylength = 3.0 * delaytmp;
      while (delaylength <= 0. || delaylength >= 2. * delaytmp) {
        sigmoi = random::randomDouble(1);
        if ((sigmoi == 1.) || (sigmoi == 0.)) {
          shift = 0.;
        } else {
          shift = width * log((1. - sigmoi) / sigmoi);
        }
        delaylength = (delaytmp + shift);
      }
      Recycling_delay = delaylength;
      break;
    }
    case Plasmacell:
    case Memorycell: {
      Recycling_delay = 0;
      break;
    }
    default: {
      cerr << "Error: B_cell::timeleft2recycle cell that is not a B-cell "
              "trying to differentiate!!"
           << endl;
    }
  }
}

void B_cell::clockreset() {
  BC_FDC_interaction_clock = 0.;
  clock = 0.;
  TCsignalDuration = 0.;
  Recycling_delay = 0.;
}

// Update the nFDCcontacts of the B_cell before differentiating (in case Ag from
// the previous round should be deleted)
void B_cell::set_Retained_Ag(parameters& p) {
//  if (p.par[DeleteAgInFreshCC])  // Parameter that determines if Ag from the
//                                 // previous round should be deleted
//  {
//    retained_Ag = 0;
//    nFDCcontacts = 0;
//    Selected_by_FDC = false;
//    Selected_by_TC=false;
//
//  } else {
//    // Selected_by_FDC = true; //ELENA: DO THEN SELECTEDFDC field =true?
//  }
}

//#Recheck @danial: check fields
string B_cell::printBcell() {
  stringstream res;

  res << "ID= " << ID << " "
      << "Cell type: " << cell_type << " "
      << "cell_state: " << cell_state << " "
      << "Cycle_state: " << cyclestate << " "
      << "cycle_time= " << cycle_state_time << " "
      << "time_to_switch= " << time_of_cycle_state_switch << " "
      << "MID= " << MID << " "
      << "nDivs2do= " << nDivisions2do << " "
      << "Position: " << position.print() << " "
      << "Polarity: " << polarity.print() << " "
      << "tp:" << persistence_time << " "
      << "CXCL12= " << isResponsive2CXCL12 << " "
      << "CXCL13= " << isResponsive2CXCL13 << " "
      << "can_move: " << can_move << " "
      << " "
      << "BCR: " << myBCR.print_BCR() << " "
      << "Affinity: " << MyAffinity << " "
      << "Mutations= " << myBCR.nMutFromGermline << " "
      << "Selected_by_TC: " << Selected_by_TC << " "
      << "pMHC_divisions= " << pMHC_dependent_number_of_divisions << " "
      << "Bc_Tc_interaction_clock= " << Bc_Tc_interaction_clock << " "
      << "BC_FDC_interaction_clock: " << BC_FDC_interaction_clock << " "
      << "clock: " << clock << " "
      << "nFDCcontacts= " << nFDCcontacts << " "
      << "retained_Ag= " << retained_Ag << " "
      << " FDC_Selected: " << Selected_by_FDC << " "
      << "TC_signal_Duration= " << TCsignalDuration << " "
      << "Interacting_TC: " << interactingTC << " "
      << "Individual_delay= " << Recycling_delay << " "
      << "High_Ag: " << IamHighAg << " "
      << "________________________________" << endl;
  return res.str();
}

T_cell::T_cell(parameters& p) : cell() {
  cell_state = TC_free;
  cell_type= TFHC;    //Type of the cell being set to counter (this means not having a type yet)
  nIncontactCCs = 0;
  speed = p.par[Tcell_speed];
  persistence_time = p.par[Tcell_tp]; // Time left for next turn
  can_move = true; //A switch to turn moving on/off
}

// Free T and B cells that are interacting by their ID.
void T_cell::liberateCC_TC(B_cell* bc) {
  nIncontactCCs -= 1;
  int ID = bc->ID;
  interactingCC.erase(std::remove_if(interactingCC.begin(), interactingCC.end(),
                                [&ID](const B_cell* x) { return x->ID == ID; }),
                      interactingCC.end());
  if (nIncontactCCs <= 0) {
    cell_state = TC_free;
    nIncontactCCs = 0;
    can_move = true;
  }
  //    free CC
  bc->interactingTC = NULL;
}

string T_cell::printTcell() {
  stringstream res;
  res << "My T-cell ID: " << ID
      << "; Number interacting CCs: " << interactingCC.size() << "; " << endl;
  for (unsigned int i = 0; i < interactingCC.size(); i = i + 1) {
    res << "interacting CC IDs: " << interactingCC.at(i)->ID
        << ", Affinity: " << interactingCC.at(i)->MyAffinity
        << "; Retained Ag: " << interactingCC.at(i)->retained_Ag << ";" << endl;
  }
  return res.str();
}

// Create a new Plasma cell with random BCR
Plasma_cell::Plasma_cell(parameters& p)
    : cell(),
      myBCR(p)
       {
  retained_Ag = 0.;
  can_move=true;
  total_number_of_divisions=0;
  birth_time=0;
  fdc_interaction_time_history=0;
  Tc_interaction_history.first=0;
  Tc_interaction_history.second=0;
  nFDCcontacts=0; // number of FDC contacts.
  Selected_by_FDC=0; //Indicates if B-cell rescued by FDC
  Selected_by_TC=0;
  delta_Affinity=0;
}

// Copy fields from Bcell into Plasma cell.
Plasma_cell::Plasma_cell(parameters& p, B_cell* Bcell) : cell(Bcell), myBCR(p) {
  MID = Bcell->ID;
  cell_type = Plasmacell;
  cell_state=Plasma_in_GC;
  MyAffinity = Bcell->MyAffinity;
  myBCR=Bcell->myBCR; //#Recheck @danial: enough?!
  myBCR.BCReceptor = Bcell->myBCR.BCReceptor;
  myBCR.nMutFromGermline = Bcell->myBCR.nMutFromGermline;
  can_move = true;
  position = Bcell->position; //#Recheck @danial: recheck assignments
  polarity = Bcell->polarity;
  speed=p.par[Plasmacell_speed];
  isResponsive2CXCL12 = false;
  isResponsive2CXCL13 = false;
  //#Recheck @danial: Are these useful?
  fdc_interaction_time_history = Bcell->fdc_interaction_time_history;
  Tc_interaction_history.first = Bcell->Tc_interaction_history.first;
  Tc_interaction_history.second = Bcell->Tc_interaction_history.second;
  nFDCcontacts = Bcell->nFDCcontacts ; // number of FDC contacts.
  retained_Ag = Bcell->retained_Ag;
  total_number_of_divisions = Bcell->total_number_of_divisions;
  delta_Affinity = Bcell->delta_Affinity;
  Selected_by_FDC= Bcell->Selected_by_FDC; //Indicates if B-cell rescued by FDC
  Selected_by_TC=Bcell->Selected_by_TC;
  ///RRR
  sequence=Bcell->sequence;
  REGIONS=Bcell->REGIONS;
  protein_sequence=Bcell->protein_sequence;
  MUTATIONS = Bcell->MUTATIONS;
  ORF=Bcell->ORF;
  SILENT_MUTATIONS = Bcell->SILENT_MUTATIONS;
  AA_REGIONS=Separating(MUTATIONS, protein_sequence);
  AA_MOM_REGIONS=Bcell->AA_REGIONS;
  NGly_mutated_sites =Bcell->NGly_mutated_sites;
  All_NGly_mutated_sites = Bcell->All_NGly_mutated_sites;
  germline_name=Bcell->germline_name;
  blacklist_nums=Bcell->blacklist_nums;
  AA_GERM_REGIONS =Bcell->AA_GERM_REGIONS;
  ///RRR

    ///Elena-network
    BCL6 = Bcell->BCL6;
    IRF4 = Bcell->IRF4;
    BLIMP1 = Bcell->BLIMP1;
    BCL6_0=BCL6; ///RRR
    IRF4_0=IRF4; ///RRR
    BLIMP1_0=BLIMP1; ///RRR
    ///Elena-network

}

// Creates a new M_cell with random BCR
Memory_cell::Memory_cell(parameters& p)
    : cell(),
      myBCR(p) {
  // Calculate cell's affinity (NOT NAN!!!!!)
  MyAffinity = myBCR.getMyAffinity4Ag(p);
  retained_Ag = 0.;
  can_move=true;
  total_number_of_divisions=0;
  birth_time=0;
  fdc_interaction_time_history=0;
  Tc_interaction_history.first=0;
  Tc_interaction_history.second=0;
  nFDCcontacts=0; // number of FDC contacts.
  Selected_by_FDC=0; //Indicates if B-cell rescued by FDC
  Selected_by_TC=0;
  delta_Affinity=0;
}

// Copy some fields from mother to Memory Bcell.
Memory_cell::Memory_cell(parameters& p, B_cell* Bcell) : cell(Bcell), myBCR(p) {

  MID = Bcell->ID;
  cell_type = Memorycell;
  cell_state=Plasma_in_GC;
  MyAffinity = Bcell->MyAffinity;
  myBCR=Bcell->myBCR; //#Recheck @danial: enough?!
  myBCR.BCReceptor = Bcell->myBCR.BCReceptor;
  myBCR.nMutFromGermline = Bcell->myBCR.nMutFromGermline;
  can_move = true;
  position = Bcell->position; //#Recheck @danial: recheck assignments
  polarity = Bcell->polarity;
  speed=p.par[Plasmacell_speed];
  isResponsive2CXCL12 = false;
  isResponsive2CXCL13 = false;
  //#Recheck @danial: Are these useful?
  fdc_interaction_time_history = Bcell->fdc_interaction_time_history;
  Tc_interaction_history.first = Bcell->Tc_interaction_history.first;
  Tc_interaction_history.second = Bcell->Tc_interaction_history.second;
  nFDCcontacts = Bcell->nFDCcontacts ; // number of FDC contacts.
  retained_Ag = Bcell->retained_Ag;
  total_number_of_divisions = Bcell->total_number_of_divisions;
  delta_Affinity = Bcell->delta_Affinity;
  Selected_by_FDC= Bcell->Selected_by_FDC; //Indicates if B-cell rescued by FDC
  Selected_by_TC=Bcell->Selected_by_TC;
  ///RRR
  sequence=Bcell->sequence;
  REGIONS=Bcell->REGIONS;
  protein_sequence=Bcell->protein_sequence;
  MUTATIONS = Bcell->MUTATIONS;
  ORF=Bcell->ORF;
  SILENT_MUTATIONS = Bcell->SILENT_MUTATIONS;
  AA_REGIONS=Separating(MUTATIONS, protein_sequence);
  AA_MOM_REGIONS=Bcell->AA_REGIONS;
  NGly_mutated_sites =Bcell->NGly_mutated_sites;
  All_NGly_mutated_sites = Bcell->All_NGly_mutated_sites;
  germline_name=Bcell->germline_name;
  blacklist_nums=Bcell->blacklist_nums;
  AA_GERM_REGIONS =Bcell->AA_GERM_REGIONS;
  ///RRR



    ///Elena-network
    BCL6 = Bcell->BCL6;
    IRF4 = Bcell->IRF4;
    BLIMP1 = Bcell->BLIMP1;
    BCL6_0=BCL6; ///RRR
    IRF4_0=IRF4; ///RRR
    BLIMP1_0=BLIMP1; ///RRR
    ///Elena-network
}

//#Recheck danial: this function can be written in a much more efficient way.
void redo_move(vector<vector3D>& redo_list, lattice& l) {
  
  int counter = int(redo_list.size()) - 1;
    
//    cout<<"time="<<time<<" counter="<<counter<<endl;
//    if (counter>0)
//    {
//        cell* ctmp = l.cellat(redo_list[counter]);
//        cout<<"Cell_id="<<ctmp->ID<<" type="<<ctmp->cell_type<<endl;
//
//    }
    
  while (counter >= 0) {
    // check the condition of last cell in the redo list
    cell* c1 = l.cellat(redo_list[counter]);
    bool exchange = false;
    if (c1 != NULL) {
      if (c1->cell_type == celltype::empty) {
        cout << "empty node has been selected for swaping" << endl;
      }
      else if (c1->cell_type == TFHC){
          T_cell* tmp_c = (T_cell*)l.cellat(redo_list[counter]);
          if (tmp_c->cell_state == TC_free)
              exchange = true;
          tmp_c=NULL;
          delete tmp_c;

      }
      else if (c1->cell_type == Centroblast)  {
        B_cell* tmp_c = (B_cell*)l.cellat(redo_list[counter]);
         if (not(tmp_c->cyclestate == cycle_M))
             exchange = true;
          tmp_c=NULL;
          delete tmp_c;
      }
      else if (c1->cell_type == Centrocyte){
          B_cell* tmp_c = (B_cell*)l.cellat(redo_list[counter]);
          if (not(tmp_c->cell_state == contact_FDC || tmp_c->cell_state == contact_TC))
              exchange = true;
          tmp_c=NULL;
          delete tmp_c;
      }
      else if (c1->cell_type == Plasmacell || c1->cell_type == Memorycell) {
        exchange = true;
      }
      else if (c1->cell_type == border) {
        cout << "border has been selected for swaping" << endl;
      }
    }
    // check a posible neighbour for swap
    bool exchange2 = false;
    if (exchange) {
      vector3D swaping_neighbour = l.get_nn_directed2(c1);
      if (swaping_neighbour.X != -1)  // to check if it finds a destination
      {
        if ((l.celltypeat(swaping_neighbour) == celltype::empty)&&(l.insideBorders(swaping_neighbour))) {
//            cout<<"Empty destination in neighbourhood, why swap?"<<" Id= "<<c1->ID<<" type="<<c1->cell_type<<endl;
        } else {
          cell* c2 = l.cellat(swaping_neighbour);
          switch (c2->cell_type) {
            case FDCell: {
              break;
            }
            case Stromalcell: {
              break;
            }
            case TFHC: {
                T_cell* tmp_c = (T_cell*)l.cellat(swaping_neighbour);
              if (tmp_c->cell_state == TC_free) {
                exchange2 = true;
              tmp_c=NULL;
              delete tmp_c;
              }
              break;
            }
            case Centroblast: {
                
                B_cell* tmp_c = (B_cell*)l.cellat(swaping_neighbour);
                if (not(tmp_c->cyclestate == cycle_M))
                    exchange = true;
                tmp_c=NULL;
                delete tmp_c;
              break;
            }
            case Centrocyte: {
                B_cell* tmp_c = (B_cell*)l.cellat(swaping_neighbour);

              if (not(tmp_c->cell_state == contact_FDC ))
                  if (not(tmp_c->cell_state == contact_TC))
                      exchange2 = true;
              tmp_c=NULL;
              delete tmp_c;
              break;
            }
            case Plasmacell:
              exchange2 = true;
              break;
            case Memorycell:
              exchange2 = true;
              break;
            case border: {
              cout << "Swaping neighbour is border." << endl;
              break;
            }
            case cell_type_counter:
              break;
            default:
              break;
          }

          // find swaping neighbour in redo list
          int index = -1;
          for (int j = 0; j <= counter; j++) {
            if (redo_list.at(j).X == c2->position.X)
              if (redo_list.at(j).Y == c2->position.Y)
                if (redo_list.at(j).Z == c2->position.Z)
                    index = j;
          }

          // swap cells
          if ((exchange) && (exchange2) && (index > -1)) {
            if (getScalarproduct(c1->polarity, c2->polarity) < 1e-6) {
              vector3D tmp_pos = c1->position;
              l.removecellat(c1->position);
              l.removecellat(c2->position);
              c1->position = c2->position;
              l.putcellat(c1);
              c2->position = tmp_pos;
              l.putcellat(c2);
              //                            cout<<"succeful swap"<<endl;
            }
            
          }
          if (index > -1) {
              redo_list.erase(redo_list.begin() + index);
            counter--;
          }
        }
      }
    }
    redo_list.pop_back();
    counter--;
  }
}
//#Recheck danial:improvement
void B_cell::proliferate(parameters& p, lattice& l, double time,
                          vector<B_cell*>& ListB_cell, output& currentoutput,
                          simulation& currentsimulation) {
  if (nDivisions2do > 0) {
    vector3D freeSpace2Divide = l.get_position_mitosis(position);
    if (freeSpace2Divide.X != -1 && freeSpace2Divide.Y != -1 &&
        freeSpace2Divide.Z != -1) {
      // Danial: Here we clean the cell that was dividing and use it as one of
      // daughter girls , hence, we need to change its ID and then record every thing
      // and then clean it. Recorded data are for mother cell, so we store ID
      // of mother in MID of daughter_Bcell to later store it in current Bcell "this"
      // daughter cells as MID.
      // Danial: The ID of Mother cell changes inside the constructore of B cell
        
      B_cell* daughter_Bcell = new B_cell(p, this);  // Create a copy of Bcell
        
      //#event  writing data of divided cell
      currentoutput.close_event(this, currentsimulation.sim_output, time);
      currentoutput.write_event(this, currentsimulation.sim_output);

        
      // Cleaning divided cell to use as a daughter cell
      MID = daughter_Bcell->MID ; //MID of "this" Bcell is equal to MID of other daughter cell which is called daughter_Bcell
      // clear stringstream #event
      event.str(string());
      event << ID << "," << time << "," << MID << ",";
      daughter_Bcell->event << daughter_Bcell->ID << "," << time << "," << daughter_Bcell->MID << ",";
      nFDCcontacts = 0.;
      interactingTC = NULL;
      cycle_state_time = 0.;
      Bc_Tc_interaction_clock = 0.;
      Recycling_delay = 0.;  // Time spent inside CBgoingLZ (time for moving to Light Zone)
      BC_FDC_interaction_clock = 0.;  // Time since a B_cell became CC_free (in sec).
      TC_selected_clock=0.0;
      Selected_by_FDC = false;  
      Selected_by_TC = false;   
      delta_Affinity = 0.0;
      TCsignalDuration = 0.;  // Acumulated signal from currently interacting TC (in
        // sec).(As imput to ODE)
      Tc_interaction_history.first = 0.0;
      Tc_interaction_history.second = 0.0;
      fdc_interaction_time_history = 0.0;
      //#Check @danial what does this mean then? total number of divisions?
      // should it be same for both daughter cells then?! I think yes.
      total_number_of_divisions += 1;
      daughter_Bcell->total_number_of_divisions += 1;
      //#Recheck @danial
      daughter_Bcell->position = freeSpace2Divide;  // Put daughter cell in free position.
      daughter_Bcell->polarity = polarity; //#Recheck @danial: shouldn't it be random?
      daughter_Bcell->nDivisions2do -= 1;
      nDivisions2do -= 1;
      if (daughter_Bcell->nDivisions2do <= 0) {
        daughter_Bcell->cyclestate = cycle_G0;
      } else {
        daughter_Bcell->cyclestate = cycle_G1;
        daughter_Bcell->time_of_cycle_state_switch =
            random::cell_cycle_time(p.par[c_G1], cycle_G1);
      }
      if (nDivisions2do <= 0) {
        cyclestate = cycle_G0;
      } else {
        cyclestate = cycle_G1;
        time_of_cycle_state_switch =
            random::cell_cycle_time(p.par[c_G1], cycle_G1);
      }
      bool asymmetric_division = false;
      if (random::randomDouble(1) < p.par[pDivideAgAssymetric]) {
        asymmetric_division = true;
      }
      if (time >= p.par[StartMutation]) {
        if (not(asymmetric_division) || not(retained_Ag > 0)) {
          daughter_Bcell->setMyAffinity(p);
          setMyAffinity(p);


          daughter_Bcell->delta_Affinity = daughter_Bcell->mutate(p);
                     // std::cout<< "1" << std::endl;
                     // std::cout<< protein_sequence << std::endl;


          delta_Affinity = mutate(p);

            ///R-SHM
          //  std::cout<< protein_sequence << std::endl;
          //std::cout<< "2" << std::endl;
            /// R-SHM
          daughter_Bcell->setMyAffinity(p);
          setMyAffinity(p);
        }
      }


    ///Elena-network
    double all_bcl6 = BCL6;
    double all_irf4 = IRF4;
    double all_blimp1 = BLIMP1;
    ///Elena-network

      if (asymmetric_division) {
        daughter_Bcell->IamHighAg = false;
        double all_ag = retained_Ag;
        double pitmp = 1.0;
        double shift = 0.;
        double gaussf;
        double width = 0.;
        if (pitmp < 1. - 4.0 * 0.04) {
          width = 0.04 * pitmp;
        } else {
          width = 0.04 * ((1.0 - pitmp)) * pitmp;  // linear switch
        }
        double tmp = 3.0;
        while (tmp < 0. || tmp > 1.) {
          gaussf = random::randomDouble(1);
          if ((gaussf == 1.) || (gaussf == 0.)) {
            shift = 0;
          } else {
            shift = width * log((1. - gaussf) / gaussf);
          }
          tmp = (pitmp + shift);
        }
        pitmp = tmp;
        retained_Ag = pitmp * all_ag;
        daughter_Bcell->retained_Ag = all_ag - retained_Ag;
        ///Elena-network
        //Elena: network: Divide TF levels asymmetrically among daughter cells
        //pitmp = rand() % 2; ////////////////////// DIVIDE TF FULLY ASYMMETRICALLY
        BCL6 = pitmp*all_bcl6;
        daughter_Bcell->BCL6 = all_bcl6 - BCL6;
        IRF4 = pitmp*all_irf4;
        daughter_Bcell->IRF4 = all_irf4 - IRF4;
        BLIMP1 = pitmp*all_blimp1;
        daughter_Bcell->BLIMP1 = all_blimp1 - BLIMP1;
        ///Elena-network
        ///R TEMP
//        BCL6 = all_bcl6/2;
//        daughter_Bcell->BCL6 = all_bcl6 - BCL6;
//        IRF4 = all_irf4/2;
//        daughter_Bcell->IRF4 = all_irf4 - IRF4;
//        BLIMP1 = all_blimp1/2;
//        daughter_Bcell->BLIMP1 = all_blimp1 - BLIMP1;
        /// R TEMP
      } else {
        double all_ag = retained_Ag;
        daughter_Bcell->retained_Ag = double(all_ag / 2);
        retained_Ag = double(all_ag / 2);
        IamHighAg = false;
        daughter_Bcell->IamHighAg = false;
        ///Elena-network
        //Elena: network: Divide TF levels symmetrically among daughter cells
        BCL6 = all_bcl6/2;
        daughter_Bcell->BCL6 = all_bcl6 - BCL6;
        IRF4 = all_irf4/2;
        daughter_Bcell->IRF4 = all_irf4 - IRF4;
        BLIMP1 = all_blimp1/2;
        daughter_Bcell->BLIMP1 = all_blimp1 - BLIMP1;
        ///Elena-network
      }
      l.putcellat(daughter_Bcell);                   // Update lattice
      ListB_cell.push_back(::move(daughter_Bcell));  // Update BcellList
    } else {
        //cout<<"No free space found for division."<<endl;
        
        }
  } else {
    cyclestate = cycle_G0;
  }
}

//#Recheck improve
void B_cell::transmit_CCdelay2cycle(parameters& p) {
  double waited_time = TC_selected_clock;
  double dtphase;
  if (cyclestate == cycle_G1) {
    dtphase = p.par[c_G1];
  } else if (cyclestate == cycle_G2) {
    dtphase = p.par[c_G2];
  } else if (cyclestate == cycle_S) {
    dtphase = p.par[c_S];
  } else {
    dtphase = 0;
  }
  while ((cyclestate != cycle_M) && (waited_time >= dtphase)) {
    waited_time -= dtphase;
    if (cyclestate == cycle_G1) {
      cyclestate = cycle_S;
      dtphase = p.par[c_S];
    } else if (cyclestate == cycle_S) {
      cyclestate = cycle_G2;
      dtphase = p.par[c_G2];
    } else if (cyclestate == cycle_G2) {
      cyclestate = cycle_M;
      dtphase = p.par[c_M];
    }
  }
  cycle_state_time = waited_time;
  if (cyclestate == cycle_G1) {
    time_of_cycle_state_switch =
        random::cell_cycle_time(p.par[c_G1], cycle_G1);
  } else if (cyclestate == cycle_G2) {
    time_of_cycle_state_switch =
        random::cell_cycle_time(p.par[c_G2], cycle_G2);
  } else if (cyclestate == cycle_S) {
    time_of_cycle_state_switch = random::cell_cycle_time(p.par[c_S], cycle_S);
  } else if (cyclestate == cycle_M) {
    time_of_cycle_state_switch = random::cell_cycle_time(p.par[c_M], cycle_M);
  } else {
    time_of_cycle_state_switch = 0;
  }
  if (cycle_state_time > time_of_cycle_state_switch) {
    cycle_state_time = time_of_cycle_state_switch;
  }
}

string cell::print_neighbours(lattice& l) {
  stringstream res;
  vector<vector3D> neighbours_6 = l.getNeighbour_nn(position);
  vector<vector3D> neighbours_12 = l.getNeighbour_diag(position);
  res << "Neighbour_6 info of Cell ID (" << ID << ") with position "
      << position.print() << ":" << endl
      << "1_" << neighbours_6[0].print()
      << " type=" << l.celltypeat(neighbours_6[0]) << endl
      << "2_" << neighbours_6[1].print()
      << " type=" << l.celltypeat(neighbours_6[1]) << endl
      << "3_" << neighbours_6[2].print()
      << " type=" << l.celltypeat(neighbours_6[2]) << endl
      << "4_" << neighbours_6[3].print()
      << " type=" << l.celltypeat(neighbours_6[3]) << endl
      << "5_" << neighbours_6[4].print()
      << " type=" << l.celltypeat(neighbours_6[4]) << endl
      << "6_" << neighbours_6[5].print()
      << " type=" << l.celltypeat(neighbours_6[5]) << endl;
  //    res<<endl<<"Neighbour_12 info:"<<endl
  //    <<"1_"<<neighbours_6[0].print()<<endl
  //    <<"2_"<<neighbours_6[1].print()<<endl
  //    <<"3_"<<neighbours_6[2].print()<<endl
  //    <<"4_"<<neighbours_6[3].print()<<endl
  //    <<"5_"<<neighbours_6[4].print()<<endl
  //    <<"6_"<<neighbours_6[5].print()<<endl
  //    <<"7_"<<neighbours_6[6].print()<<endl
  //    <<"8_"<<neighbours_6[7].print()<<endl
  //    <<"9_"<<neighbours_6[8].print()<<endl
  //    <<"10_"<<neighbours_6[9].print()<<endl
  //    <<"11_"<<neighbours_6[10].print()<<endl
  //    <<"12_"<<neighbours_6[11].print()<<endl;
  return res.str();
}

//#Check
// Set time left for next polarity calculation.
void cell::getNewPersistentTime(parameters& p) {
  //    double value_tp = -1;
  //    double std_tp = -1;
  //    typeParameter typePersistence = N_types_parameters;
  //    switch(cell_type)
  //        {
  //            case Centroblast:
  //                {
  //                    value_tp = p.par[Bcell_tp];
  //                    std_tp = p.par[Bcell_tp_stddev];
  //                    typePersistence = (typeParameter) (int)
  //                    p.par[type_tp_Bcells];
  //                    break;
  //                }
  //            case Centrocyte:
  //                {
  //                    value_tp = p.par[Bcell_tp];
  //                    std_tp = p.par[Bcell_tp_stddev];
  //                    typePersistence = (typeParameter) (int)
  //                    p.par[type_tp_Bcells];
  //                    break;
  //                }
  //            case TFHC:
  //                {
  //                    value_tp = p.par[Tcell_tp];
  //                    std_tp = p.par[Tcell_tp_stddev];
  //                    typePersistence = (typeParameter) (int)
  //                    p.par[type_tp_Tcells];
  //                    break;
  //                }
  //            case Plasmacell:
  //                {
  //                    value_tp = p.par[Plasmacell_tp];
  //                    std_tp = p.par[Plasmacell_tp_stddev];
  //                    typePersistence = (typeParameter) (int)
  //                    p.par[type_tp_Plasmacells];
  //                    break;
  //                }
  //            case Memorycell:
  //                {
  //                    value_tp = p.par[Memorycell_tp];
  //                    std_tp = p.par[Memorycell_tp_stddev];
  //                    typePersistence = (typeParameter) (int)
  //                    p.par[type_tp_Memorycells];
  //                    break;
  //                }
  //            default:
  //                {
  //                    cerr << "ERR: function getNewPersistenettime was called
  //                    from a non-movable cell" << endl;
  //                }
  //        }
  //    if (value_tp == NAN) cerr<<"Error with getNewpersistentTime: parameters
  //    not found"<<endl;
  //    if (value_tp > 60.0 * p.par[dt])
  //        {
  //            tp = double (60.0 * p.par[dt] / value_tp);
  //        }
  //    else
  //        {
  //            tp = 1.0;   // i.e. change polarity in every time step!
  //        }
}

double cell::set_speed(parameters& p) {
  // This function can be used to change the speed of cell
  double p_move = 0;
  return p_move;
}

void T_cell::getNewPolarity(parameters& p, lattice& l) {
    
    getRandomPolarity(p, l);  // Set polarity to a random vector.
    
    if (cell_type == TFHC) {
        double north[3];
        north[0] = 0.;
        north[1] = 0.;
        north[2] = -1.;
        polarity.X = (1.0 - 0.1) * polarity.X + 0.1 * north[0];
        polarity.Y = (1.0 - 0.1) * polarity.Y + 0.1 * north[1];
        polarity.Z = (1.0 - 0.1) * polarity.Z + 0.1 * north[2];
        polarity.getNormalizedVector();
        return;
    }

}

void Plasma_cell::getNewPolarity(parameters& p, lattice& l) {
    
    getRandomPolarity(p, l);  // Set polarity to a random vector.
    if (cell_type == Plasmacell) {
        if (polarity.Z < 0.) {
            polarity.Z *= -1.;
        }
    }
}

void Memory_cell::getNewPolarity(parameters& p, lattice& l) {

    getRandomPolarity(p, l);  // Set polarity to a random vector.
    if (cell_type == Memorycell) {
        if (polarity.Z < 0.) {
            polarity.Z *= -1.;
        }
    }
}

//#Recheck @danial: Improvment
void B_cell::getNewPolarity(parameters &p, lattice &l)
{
        getRandomPolarity(p, l);  // Set polarity to a random vector.
    
        if (isResponsive2CXCL12 == true) {
            if (l.chemoat(CXCL12, position) < p.par[CXCL12crit]) {
                vector3D CXCL12diff = {
                    double(l.chemoat(CXCL12, position.X + 1, position.Y, position.Z) -
                           l.chemoat(CXCL12, position.X - 1, position.Y, position.Z)),
                    double(l.chemoat(CXCL12, position.X, position.Y + 1, position.Z) -
                           l.chemoat(CXCL12, position.X, position.Y - 1, position.Z)),
                    double(l.chemoat(CXCL12, position.X, position.Y, position.Z + 1) -
                           l.chemoat(CXCL12, position.X, position.Y, position.Z - 1))};
                
                // Danial: normalize gradient vector
                double gradient = CXCL12diff.getNorm();
                CXCL12diff = CXCL12diff.getNormalizedVector();
                
                double alpha =
                double(p.par[chemmax] /
                       (1 + exp(p.par[chemosteep] * (p.par[chemohalf] - gradient))));
                
                polarity.X = polarity.X + (alpha * CXCL12diff.X);
                polarity.Y = polarity.Y + (alpha * CXCL12diff.Y);
                polarity.Z = polarity.Z + (alpha * CXCL12diff.Z);
                polarity.getNormalizedVector();
            } else {
                isResponsive2CXCL12 = false;
            }
        }
        if (isResponsive2CXCL13 == true) {
            if (l.chemoat(CXCL13, position) < p.par[CXCL13crit]) {
                vector3D CXCL13diff = {
                    double(l.chemoat(CXCL13, position.X + 1, position.Y, position.Z) -
                           l.chemoat(CXCL13, position.X - 1, position.Y, position.Z)),
                    double(l.chemoat(CXCL13, position.X, position.Y + 1, position.Z) -
                           l.chemoat(CXCL13, position.X, position.Y - 1, position.Z)),
                    double(l.chemoat(CXCL13, position.X, position.Y, position.Z + 1) -
                           l.chemoat(CXCL13, position.X, position.Y, position.Z - 1))};
                
                // Danial: normalize gradient vector
                double gradient = CXCL13diff.getNorm();
                CXCL13diff = CXCL13diff.getNormalizedVector();
                
                double alpha = double(
                                      10 / (1 + exp(p.par[chemosteep] * (p.par[chemohalf] - gradient))));
                
                polarity.X = polarity.X + (CXCL13diff.X * alpha);
                polarity.Y = polarity.Y + (CXCL13diff.Y * alpha);
                polarity.Z = polarity.Z + (CXCL13diff.Z * alpha);
                polarity.getNormalizedVector();
                
            } else {
                isResponsive2CXCL13 = false;
            }
        }
    
}

bool B_cell::IsCycling() { return (cyclestate != cycle_G0); }

void  B_cell::Resensitize2Chemokines(parameters& p, lattice& l)
{
    switch(cell_type){
        case Centrocyte:
        {
            if(l.chemoat(CXCL13, position)< p.par[CXCL13recrit])
            {
                isResponsive2CXCL13 = true;
            }
            break;
        }
        case Centroblast:
        {
            if(l.chemoat(CXCL12, position) < p.par[CXCL12recrit])
            {
                isResponsive2CXCL12 = true;
            }
            break;
        }
        default: {
            isResponsive2CXCL12 = false;
            isResponsive2CXCL13 = false;
            break;}
    }
}
