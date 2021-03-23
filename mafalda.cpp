
#include "mafalda.h"
#include "random.h"
#include "lattice.h"
#include "cell.h"
#include "GC3D.h"
#include <iostream>
#include <vector>
#include <list>
#include <cmath>
#include "bcr.h"
#include "output.h"
#include <algorithm> ///RRR
bool pause1 = false;
using namespace std;
std::mt19937 g(global_seed +19131); ///RRR

vector3D lattice::get_fdc_position(int j)
{
    return fdcposes[j];
}

void  lattice::set_initial_fdc_position()
{
    fdcposes.reserve(200);
//
//        FILE *fdcpos = fopen("/Users/danial/Desktop/fdcpos.csv", "r");
//        int x,y,z;
//        for (int i=0;i<200;i++){
//            fscanf(fdcpos, "%d,%d,%d", &x, &y,&z);
//            fdcposes.push_back(vector3D(x,y,z));
//        }
//        fclose(fdcpos);

    fdcposes.push_back(vector3D(35,38,23));
    fdcposes.push_back(vector3D(41,32,25));
    fdcposes.push_back(vector3D(32,32,15));
    fdcposes.push_back(vector3D(37,35,23));
    fdcposes.push_back(vector3D(42,31,23));
    fdcposes.push_back(vector3D(33,32,18));
    fdcposes.push_back(vector3D(35,40,24));
    fdcposes.push_back(vector3D(33,32,14));
    fdcposes.push_back(vector3D(33,24,19));
    fdcposes.push_back(vector3D(42,31,27));
    fdcposes.push_back(vector3D(31,23,18));
    fdcposes.push_back(vector3D(40,32,19));
    fdcposes.push_back(vector3D(38,40,25));
    fdcposes.push_back(vector3D(29,24,18));
    fdcposes.push_back(vector3D(36,40,19));
    fdcposes.push_back(vector3D(28,36,18));
    fdcposes.push_back(vector3D(35,26,27));
    fdcposes.push_back(vector3D(34,25,15));
    fdcposes.push_back(vector3D(47,30,26));
    fdcposes.push_back(vector3D(33,41,19));
    fdcposes.push_back(vector3D(29,26,15));
    fdcposes.push_back(vector3D(31,28,29));
    fdcposes.push_back(vector3D(25,38,21));
    fdcposes.push_back(vector3D(29,41,23));
    fdcposes.push_back(vector3D(35,46,30));
    fdcposes.push_back(vector3D(42,40,22));
    fdcposes.push_back(vector3D(36,46,29));
    fdcposes.push_back(vector3D(46,25,21));
    fdcposes.push_back(vector3D(16,27,22));
    fdcposes.push_back(vector3D(19,32,27));
    fdcposes.push_back(vector3D(15,27,24));
    fdcposes.push_back(vector3D(16,29,26));
    fdcposes.push_back(vector3D(48,28,21));
    fdcposes.push_back(vector3D(17,30,22));
    fdcposes.push_back(vector3D(34,28,31));
    fdcposes.push_back(vector3D(20,24,28));
    fdcposes.push_back(vector3D(21,28,29));
    fdcposes.push_back(vector3D(29,30,9));
    fdcposes.push_back(vector3D(29,29,9));
    fdcposes.push_back(vector3D(42,21,24));
    fdcposes.push_back(vector3D(48,30,28));
    fdcposes.push_back(vector3D(18,23,24));
    fdcposes.push_back(vector3D(35,43,31));
    fdcposes.push_back(vector3D(40,17,24));
    fdcposes.push_back(vector3D(33,32,15));
    fdcposes.push_back(vector3D(45,29,14));
    fdcposes.push_back(vector3D(23,40,17));
    fdcposes.push_back(vector3D(26,38,14));
    fdcposes.push_back(vector3D(46,22,28));
    fdcposes.push_back(vector3D(49,24,27));
    fdcposes.push_back(vector3D(18,37,27));
    fdcposes.push_back(vector3D(20,20,30));
    fdcposes.push_back(vector3D(48,31,16));
    fdcposes.push_back(vector3D(48,34,17));
    fdcposes.push_back(vector3D(21,21,21));
    fdcposes.push_back(vector3D(21,35,30));
    fdcposes.push_back(vector3D(47,42,20));
    fdcposes.push_back(vector3D(30,18,27));
    fdcposes.push_back(vector3D(48,42,22));
    fdcposes.push_back(vector3D(19,36,21));
    fdcposes.push_back(vector3D(48,22,24));
    fdcposes.push_back(vector3D(40,32,32));
    fdcposes.push_back(vector3D(33,36,32));
    fdcposes.push_back(vector3D(36,21,16));
    fdcposes.push_back(vector3D(37,32,11));
    fdcposes.push_back(vector3D(37,29,9));
    fdcposes.push_back(vector3D(20,25,15));
    fdcposes.push_back(vector3D(32,32,14));
    fdcposes.push_back(vector3D(35,39,24));
    fdcposes.push_back(vector3D(19,18,30));
    fdcposes.push_back(vector3D(30,53,29));
    fdcposes.push_back(vector3D(31,14,26));
    fdcposes.push_back(vector3D(35,39,23));
    fdcposes.push_back(vector3D(28,20,32));
    fdcposes.push_back(vector3D(32,50,30));
    fdcposes.push_back(vector3D(16,23,30));
    fdcposes.push_back(vector3D(29,52,29));
    fdcposes.push_back(vector3D(39,46,32));
    fdcposes.push_back(vector3D(36,38,23));
    fdcposes.push_back(vector3D(41,31,25));
    fdcposes.push_back(vector3D(42,31,24));
    fdcposes.push_back(vector3D(35,37,23));
    fdcposes.push_back(vector3D(41,23,12));
    fdcposes.push_back(vector3D(25,17,19));
    fdcposes.push_back(vector3D(42,32,25));
    fdcposes.push_back(vector3D(41,31,23));
    fdcposes.push_back(vector3D(37,36,23));
    fdcposes.push_back(vector3D(34,11,30));
    fdcposes.push_back(vector3D(35,38,22));
    fdcposes.push_back(vector3D(18,22,17));
    fdcposes.push_back(vector3D(35,38,24));
    fdcposes.push_back(vector3D(47,19,21));
    fdcposes.push_back(vector3D(37,36,10));
    fdcposes.push_back(vector3D(41,32,24));
    fdcposes.push_back(vector3D(36,12,29));
    fdcposes.push_back(vector3D(48,46,21));
    fdcposes.push_back(vector3D(33,32,17));
    fdcposes.push_back(vector3D(46,37,13));
    fdcposes.push_back(vector3D(25,33,10));
    fdcposes.push_back(vector3D(36,35,23));
    fdcposes.push_back(vector3D(32,32,16));
    fdcposes.push_back(vector3D(35,40,23));
    fdcposes.push_back(vector3D(14,35,29));
    fdcposes.push_back(vector3D(41,32,26));
    fdcposes.push_back(vector3D(34,38,23));
    fdcposes.push_back(vector3D(42,31,26));
    fdcposes.push_back(vector3D(42,32,23));
    fdcposes.push_back(vector3D(37,35,22));
    fdcposes.push_back(vector3D(40,32,25));
    fdcposes.push_back(vector3D(36,40,24));
    fdcposes.push_back(vector3D(37,34,23));
    fdcposes.push_back(vector3D(38,35,23));
    fdcposes.push_back(vector3D(31,24,18));
    fdcposes.push_back(vector3D(37,35,24));
    fdcposes.push_back(vector3D(25,46,18));
    fdcposes.push_back(vector3D(30,24,18));
    fdcposes.push_back(vector3D(41,31,24));
    fdcposes.push_back(vector3D(42,31,22));
    fdcposes.push_back(vector3D(43,34,10));
    fdcposes.push_back(vector3D(42,30,23));
    fdcposes.push_back(vector3D(42,32,24));
    fdcposes.push_back(vector3D(42,31,25));
    fdcposes.push_back(vector3D(25,48,30));
    fdcposes.push_back(vector3D(33,32,19));
    fdcposes.push_back(vector3D(32,24,19));
    fdcposes.push_back(vector3D(33,33,18));
    fdcposes.push_back(vector3D(34,32,18));
    fdcposes.push_back(vector3D(32,31,15));
    fdcposes.push_back(vector3D(33,31,18));
    fdcposes.push_back(vector3D(51,45,22));
    fdcposes.push_back(vector3D(32,32,18));
    fdcposes.push_back(vector3D(41,33,25));
    fdcposes.push_back(vector3D(53,25,22));
    fdcposes.push_back(vector3D(52,25,18));
    fdcposes.push_back(vector3D(36,37,24));
    fdcposes.push_back(vector3D(33,24,20));
    fdcposes.push_back(vector3D(37,34,22));
    fdcposes.push_back(vector3D(53,30,17));
    fdcposes.push_back(vector3D(38,34,23));
    fdcposes.push_back(vector3D(29,25,18));
    fdcposes.push_back(vector3D(35,36,23));
    fdcposes.push_back(vector3D(33,31,17));
    fdcposes.push_back(vector3D(36,40,23));
    fdcposes.push_back(vector3D(34,38,24));
    fdcposes.push_back(vector3D(33,33,15));
    fdcposes.push_back(vector3D(33,23,19));
    fdcposes.push_back(vector3D(41,33,24));
    fdcposes.push_back(vector3D(40,31,25));
    fdcposes.push_back(vector3D(34,53,31));
    fdcposes.push_back(vector3D(37,36,24));
    fdcposes.push_back(vector3D(54,37,25));
    fdcposes.push_back(vector3D(37,37,23));
    fdcposes.push_back(vector3D(34,32,14));
    fdcposes.push_back(vector3D(34,39,24));
    fdcposes.push_back(vector3D(40,33,19));
    fdcposes.push_back(vector3D(36,35,22));
    fdcposes.push_back(vector3D(40,31,19));
    fdcposes.push_back(vector3D(37,36,22));
    fdcposes.push_back(vector3D(35,41,24));
    fdcposes.push_back(vector3D(26,12,27));
    fdcposes.push_back(vector3D(33,33,14));
    fdcposes.push_back(vector3D(35,37,24));
    fdcposes.push_back(vector3D(32,24,18));
    fdcposes.push_back(vector3D(35,37,22));
    fdcposes.push_back(vector3D(19,46,27));
    fdcposes.push_back(vector3D(38,10,32));
    fdcposes.push_back(vector3D(36,38,22));
    fdcposes.push_back(vector3D(35,39,22));
    fdcposes.push_back(vector3D(33,31,15));
    fdcposes.push_back(vector3D(34,39,23));
    fdcposes.push_back(vector3D(38,40,24));
    fdcposes.push_back(vector3D(40,32,24));
    fdcposes.push_back(vector3D(31,23,19));
    fdcposes.push_back(vector3D(32,32,17));
    fdcposes.push_back(vector3D(38,39,25));
    fdcposes.push_back(vector3D(42,32,27));
    fdcposes.push_back(vector3D(39,32,19));
    fdcposes.push_back(vector3D(31,32,15));
    fdcposes.push_back(vector3D(42,32,26));
    fdcposes.push_back(vector3D(33,25,19));
    fdcposes.push_back(vector3D(33,31,14));
    fdcposes.push_back(vector3D(32,23,18));
    fdcposes.push_back(vector3D(36,36,23));
    fdcposes.push_back(vector3D(37,40,25));
    fdcposes.push_back(vector3D(33,24,18));
    fdcposes.push_back(vector3D(41,32,23));
    fdcposes.push_back(vector3D(41,31,26));
    fdcposes.push_back(vector3D(40,32,20));
    fdcposes.push_back(vector3D(34,40,24));
    fdcposes.push_back(vector3D(36,37,23));
    fdcposes.push_back(vector3D(36,39,23));
    fdcposes.push_back(vector3D(36,39,24));
    fdcposes.push_back(vector3D(43,31,23));
    fdcposes.push_back(vector3D(32,33,15));
    fdcposes.push_back(vector3D(33,32,16));
    fdcposes.push_back(vector3D(30,23,18));
    fdcposes.push_back(vector3D(36,38,24));
    fdcposes.push_back(vector3D(41,31,27));
    fdcposes.push_back(vector3D(10,30,24));
    fdcposes.push_back(vector3D(35,40,25));

        cout<<"Setting initial FDC positions"<<endl;
        cout<<"Number of FDC positions:"<< fdcposes.size()<<endl;
}


simulation::simulation(parameters &p)
{
    currentOutput = new output();
    ///Elena-network
    EventOutput = new events(outputFolder,"/historyOut.txt", "/historyDead.txt");
    ///Elena-network
    currentOutput->Output_ID=Simulation_ID;
    //#Check @danial: set local address for these.
    currentLattice = new lattice(p, base_Path + "/cxcl12_3d_5micron.sig",base_Path + "/cxcl13_3d_5micron.sig");
    ListB_cell.reserve(50000);
    ListT_cell.reserve(50000);
    ListP_cell.reserve(50000);
    ListM_cell.reserve(50000);
}
//Destructor of Simulation class
simulation::~simulation()
    {
        delete currentLattice;
        int NSC = (int)ListSC.size();
        for(int i = 0; i < NSC; ++i)
            {
                delete ListSC[i];
            }
        int NFDC = (int)ListFDC.size();
        for(int i = 0; i < NFDC; ++i)
            {
                delete ListFDC[i];
            }
        int NBC = (int)ListB_cell.size();
        for(int i = 0; i < NBC; ++i)
            {
                delete ListB_cell[i];
            }
        int NTC = (int)ListT_cell.size();
        for(int i = 0; i < NTC; ++i)
            {
                delete ListT_cell[i];
            }
    }

// Initialize cells
void simulation::InitialCells(lattice& l, parameters& p)
{
// Initialize Stromal cells
/* Creates stromal cells in the dark zone randomley, these cells are immobile and are placed on the lattice so are not transparent */
    for(unsigned int j = 0; j < p.par[InitialNumberSC] ; j= j+1)
        {
            Stromal_cell* newSC = new Stromal_cell();
            newSC->position = l.getFreePosition( p.par[zoneRatioGC],1);
            ListSC.push_back(newSC);
            newSC->cell_type=Stromalcell;
            newSC->can_move=false;
            l.putcellat(newSC);
        }
    cerr << p.par[InitialNumberSC] << " SC generated " << endl;

// Initialized FDCs
/* Creates follicular dendritic cells network, FDCs are immobile. Each Soma has some dendrities, each FDC is placed randomley in light zone and dendrities are placed in 6 directions (2 in each plan). If the dendrit goes beyond the simulation space the ag amount will be distributed over the other dendrities. FDCs are transparent which means they are not placed on the lattice but the amount of Ag is placed, so to work with FDCs you need to work with the amount of Ag in each lattice node*/
    l.set_initial_fdc_position();
    for(unsigned int j = 0; j < p.par[InitialNumberFDC] ; j= j+1)
        {
            FDC* newFDC = new FDC();
            newFDC->position = l.get_fdc_position(j);
            newFDC->can_move=false;
            newFDC->occupiedPositions.reserve(6*p.par[DendriteLength]);
            newFDC->volume=1;
            newFDC->cell_type=FDCell;
            newFDC->occupiedPositions.push_back(newFDC->position);
            int ttmp[3],tmp[3];
            ttmp[0]=newFDC->position.X;
            ttmp[1]=newFDC->position.Y;
            ttmp[2]=newFDC->position.Z;
            tmp[0]=ttmp[0];
            tmp[1]=ttmp[1];
            tmp[2]=ttmp[2];

            for (int i = 1; i <= (p.par[DendriteLength]); i++)
                {
                    for (int j= 0; j<3; j++)
                        {
                            // positive direction
                            tmp[j]= ttmp[j]+i;
                            if (l.insideBorders(vector3D(tmp[0],tmp[1],tmp[2])))
                                {
                                    newFDC->volume +=1;
                                    newFDC->occupiedPositions.push_back(vector3D(tmp[0],tmp[1],tmp[2]));
                                }
                            // negative direction
                            tmp[j]= ttmp[j]-i;
                            if (l.insideBorders(vector3D(tmp[0],tmp[1],tmp[2])))
                                {
                                    newFDC->volume +=1;
                                    newFDC->occupiedPositions.push_back(vector3D(tmp[0],tmp[1],tmp[2]));
                                }
                            tmp[j] = ttmp[j];
                        }
                }
            // Amount of Ag per dendrite
            newFDC->AgperDendrite = double(p.par[AgAmountperFDC])/double(newFDC->volume);
            for (int i = 0; i < newFDC->volume; i++)
                {
                    l.putAgFDCat(newFDC->occupiedPositions.at(i), newFDC, newFDC->AgperDendrite);
                    l.AddTotalAmountAginLattice(newFDC->AgperDendrite);
                }
            ListFDC.push_back(newFDC);
        }
    cerr << p.par[InitialNumberFDC] << " FDC generated" << endl;

//Initialize 100 affinity seeds for incoming CBs


//Initialize the Seeder B cells
/* If we have initial B cells, they are placed randomley in dark zone. B cells start the cycle in G1 phase with an affinity from initial seeds pool*/
    for(unsigned int j = 0; j < p.par[InitialNumberCB] ; j= j+1)
        {
            B_cell* newB_cell = new B_cell(p);
            newB_cell->event << newB_cell->ID << "," << 0 << "," << newB_cell->MID << ","; ///RRRR
            newB_cell->cell_state = founder; //Status of the cell being set to counter (this means not having an state yet)
            newB_cell->cell_type= Centroblast;    //Type of the cell being set to counter (this means not having a type yet)
            newB_cell->persistence_time = p.par[Bcell_tp]; // Time left for next turn
            newB_cell->speed=p.par[Bcell_speed];
            newB_cell->can_move=true;                     // A switch to turn moving on/off
            newB_cell->setMyAffinity(p);
            newB_cell->time_of_cycle_state_switch=random::cell_cycle_time(p.par[c_G1],cycle_G1);
            newB_cell->cyclestate=cycle_G1;
            newB_cell->position = l.getFreePosition( p.par[zoneRatioGC],1);
            newB_cell->polarity= l.get_random_direction();
            if(not (l.insideBorders(newB_cell->position))) cerr<<"Cell at border position: "<<newB_cell->printcell() <<endl;
            l.putcellat(newB_cell);
            newB_cell->nDivisions2do = p.par[nDiv];
            newB_cell->getNewPersistentTime(p); //#Recheck, @danial: neccessary for the moment?!
            newB_cell->isResponsive2CXCL12=true;
            newB_cell->isResponsive2CXCL13=false;

            ///Elena-network
            EventOutput->recordEvent(newB_cell, event_born,0); //Elena:  events: record event born with ID of cell to track cell history
            ///Elena-network

            ListB_cell.push_back(newB_cell);
        }
    cerr << p.par[InitialNumberCB] << " CB generated" << endl;

// Initialize T follicular helper cells
/* These T cells are mobile when they are not interacting with B cells, but they don't divide*/
    for(unsigned int j = 0; j < p.par[InitialNumberTC] ; j= j+1)
        {
            T_cell* newT_cell = new T_cell(p);
            newT_cell->cell_state = TC_free;
            newT_cell->cell_type=TFHC;
            newT_cell->can_move=true;
            newT_cell->position = l.getFreePosition( 0,p.par[zoneRatioGC]);
            newT_cell->polarity= l.get_random_direction();
            newT_cell->getNewPersistentTime(p);
            newT_cell->persistence_time = p.par[Tcell_tp];  // Danial: for now it is constant, maybe add more options later  #Recheck
            newT_cell->speed=p.par[Tcell_speed];
            l.putcellat(newT_cell);
            ListT_cell.push_back(newT_cell);
        }
    cerr << p.par[InitialNumberTC] << " TC generateds" << endl;
}

void simulation::simulate(lattice& l, parameters& p)
{
//    reallocate_memory();

//Create simulation output folder
    currentOutput->createFolder(outputFolder);
    
//Cell initiation
    InitialCells(l, p );

    int Total_time_steps = 504 / p.par[dt]; //#Temporary504
//Time loop

    cout<<Total_time_steps<<endl;
    for(int counter = 0 ; counter <= Total_time_steps; counter++) //danial: #new_debugging
        {

//Time:
            double t= double(counter)*p.par[dt]; // in hours
            //#temporary
            double recording_time_period=1.0; //in hour
            int recording_time_steps = double (recording_time_period / p.par[dt]);
            
            if (not(pause1))
            {

            // Record Output every fix ouhr
            if(fmod(counter,recording_time_steps) < p.par[dt] )
            {

                currentOutput->record_output_time_step(t, *this, p);
                cout<<"t="<<t<<endl;
            }

//Redo the movment for thoes cells which can not move due to cell trafficking
            vector <vector3D> redo_list;
            redo_list.reserve(6000);
// Vector that stores the dead cells to remove later
            going_to_delet.reserve(10000);

            std::random_device rd; ///RRR

// Shuffle the cell lists
            if (ListB_cell.size()>0)
            std::shuffle ( ListB_cell.begin(), ListB_cell.end(), g); ///RRR
            if (ListT_cell.size()>0)
            std::shuffle ( ListT_cell.begin(), ListT_cell.end(), g); ///RRR
        //  if (ListP_cell.size()>0)
        //  std::random_shuffle (ListP_cell.begin(),ListP_cell.end() );
        //  if (ListM_cell.size()>0)
        //  std::random_shuffle ( ListM_cell.begin(), ListM_cell.end() );

// Calculations for Output cells at each time step
            Calc_Out( t, p,l, redo_list);

// Calculations for T cells at each time step
            Calc_TC(p,l,redo_list);

// Calculations for B cells at each time step
            Calc_BC(t,p,l,redo_list,going_to_delet);
                
// Transfer newly differentiated Plasma cells from B cell list to avoid interference in output files.
            transfer_plasma_from_Bcell_list( t,p,l, redo_list);
//Influx of B cells to GC as an option
            BCinflux(t,p,l);

// Display simulation
            //Visualise(t,p); ///R

                
// Redo move
            bool allow_exchange= true; //Danial: added for the moment, later on put in parameters file
            if (allow_exchange)
                {
                    if (redo_list.size()>0)
                        {
                            redo_move(redo_list,l);
                        }
                }


//Remove dead cells
                clean_dead_cells(l);
 
            }
            
            else {
                counter--;
                t= double(counter)*p.par[dt]; // in hours
                // Visualise(t,p);///RRR

            }
        }
            
    cerr << "Simulation finished" << endl;
    
        for(unsigned int i = 0; i < ListB_cell.size(); i++)
        {
            B_cell* Bcell = ListB_cell.at(i);
                currentOutput->close_event(Bcell, sim_output, 504); //#Temporary 504
                currentOutput->write_event(Bcell, sim_output);
        }
    
    cerr << "writen Output files" << endl;


//    writeEvents
    cerr << "writen Events" << endl;
    currentOutput->write_event_2file(sim_output);

    currentOutput->Plasma_output(504,*this, p);  //Temporary504
    currentOutput->Memory_output(504,*this, p);  //Temporary504

    currentOutput->~output();

    fstream database;
    database.open("bcinflow09/Sequence_db.csv",  fstream::out);
    fstream myStream;
    for(auto& kv : Seq_affinity) {
      database << kv.first  <<  ',' << kv.second[0]  <<  ','  << kv.second[1]  <<  ','  << kv.second[2]  <<  ','  << kv.second[3]  <<  ',' << "\n";
    }
    database.close();


    fstream database2;
    database2.open("bcinflow09/Restricted_db.csv",  fstream::out);
    fstream myStream2;
    for(auto& kv : Restricted_db) {
      for(auto& kv2 : kv.second) {
        for(auto& kv3 : kv2.second) {
            for(auto& kv4 : kv3.second) {
              database2 << kv.first <<  ',' << kv2.first <<  ',' << kv3.first <<  ',' << kv4.first <<  ','  << kv4.second <<  ',' << "\n";
            }      }
      }

    }
    database2.close();

    fstream databaseCDR;
    databaseCDR.open("bcinflow09/CDRs_db.csv",  fstream::out);
    fstream myStreamCDR;
    for(auto& kv : CDRs_affinity) {
      databaseCDR << kv.first[0]  <<  ',' << kv.first[1]  <<  ',' << kv.first[2]  <<  ',' << kv.second[0]  <<  ','  << kv.second[1]  <<  ','  << kv.second[2]  <<  ','  << kv.second[3]  <<  ',' << "\n";
    }
    databaseCDR.close();

    fstream databaseFWR;
    databaseFWR.open("bcinflow09/FWRs_db.csv",  fstream::out);
    fstream myStreamFWR;
    for(auto& kv : FWRs_db) {
      databaseFWR << kv.first[0]  <<  ',' << kv.first[1]  <<  ','  << kv.second << "\n";
    }
    databaseFWR.close();

 //#Temporary
    //#Recheck @danial: delete all dynamically allocated memories here
}

// Influx of B cells into the GC
/* B cells influx to GC by a probability that can change with time, they find a random position in the whole GC to enter*/
void simulation::BCinflux(double time,parameters &p, lattice &l)
{
    double pBCinflux = double ((p.par[rateCBinflow]* p.par[dt])) /double((1.0 + exp((time - p.par[timeStopCBinflow])/p.par[smoothnessStopCBinflow])));
    if (random::randomDouble(1.) < pBCinflux)
        {
            B_cell* Bcell = new B_cell(p);
            // Initialize
            Bcell->cell_state=founder;
            Bcell->cell_type = Centroblast;
            Bcell->persistence_time=p.par[Bcell_tp];
            Bcell->speed= p.par[Bcell_speed];
            Bcell->can_move=true;
            Bcell->setMyAffinity(p);
            Bcell->cyclestate=cycle_G1;
            Bcell->time_of_cycle_state_switch=random::cell_cycle_time(p.par[c_G1],cycle_G1);
            Bcell->position = l.getFreePosition(0,1); // Take free position
            if(not(l.insideBorders(Bcell->position)))
            {cerr<<"Influx of B-cell at border position: "<<Bcell->printcell() <<endl;
                exit(1);}
            Bcell->polarity=l.get_random_direction();
            
            Bcell->nDivisions2do =p.par[nDivinflow]; //#Temporary
            
            Bcell->getNewPersistentTime(p);
            Bcell->myBCR.pMut= p.par[pmutAfterStartMut];
            Bcell->isResponsive2CXCL12=true;
            Bcell->isResponsive2CXCL13=false;
            //Danial: #check #event_record
            Bcell->event<<Bcell->ID<<","<<time<<","<<Bcell->MID<<",";

            ///Elena-network
            EventOutput->recordEvent(Bcell, event_born, time); //Elena:  events: record event born with ID of cell to track cell history
            ///Elena-network

            l.putcellat(Bcell);
            ListB_cell.push_back(Bcell);
        }
}

// Calculation of T cells
void simulation::Calc_TC(parameters &p, lattice &l, vector<vector3D> &redo_list)
{
    int N_T_cell = int (ListT_cell.size());
    for(int i = 0; i<N_T_cell; i++)
        {
            T_cell* Tcell = ListT_cell.at(i);
            switch (Tcell->cell_state)
                {
                    case TC_free:
                        {
                            Tcell->can_move=true;
                            Tcell->move(p,l,redo_list);
                            break;
                        }
                    case TC_connected:
                        {
                            Tcell->can_move=false;
//Tcells in contact to CCs do not move!
                            sort(Tcell->interactingCC.begin(),Tcell->interactingCC.end(), [](const B_cell* x, const B_cell* y){ return (x->retained_Ag > y->retained_Ag);}); //sorts descendingly based on retained Ag
                            //#check @Danil: Miachel is using the Number of times that cells picked up Ag to help cells not retained Ag
                            Tcell->polarity.X= double( Tcell->interactingCC[0]->position.X - Tcell->position.X);
                            Tcell->polarity.Y=double(Tcell->interactingCC[0]->position.Y - Tcell->position.Y);
                            Tcell->polarity.Z=double(Tcell->interactingCC[0]->position.Z - Tcell->position.Z);
                            break;
                        }
                }
        }
}

// Visualization function
void simulation::Visualise(double t, parameters &p)
{
    
      if (t<=1)
        {
//            glm::vec3 cameraPosition(10.0f, 20.0f, 10.0f+ cameraDistance);
//            gluLookAt(100.0,100.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
            

//            nextToDisplay(&ListB_cell, &ListT_cell, &ListFDC, &ListP_cell, NULL, &ListSC, currentLattice, t);
//
//            display();
//             glutMainLoopEvent();
        }
    else if (t>5 && t<85.)
        {
//            gluPerspective(60.0, 1.0, 1, 500); // Note : deph test works only if the first plane is > 0
          
//            gluLookAt(10, 10, 10,  0, 0, 0, 0, 1, 0);

            if(fmod(t+1e-9,1)< p.par[dt])
                {
//                    nextToDisplay(&ListB_cell, &ListT_cell, &ListFDC, &ListP_cell, NULL, &ListSC, currentLattice, t);
//                    display();
//                     glutMainLoopEvent();
                }
        }
    else if (t>85. && t< 87.)
        {
//            nextToDisplay(&ListB_cell, &ListT_cell, &ListFDC, &ListP_cell, NULL, &ListSC, currentLattice, t);
//            display();
//             glutMainLoopEvent();
        }
    else if (t>=87)
        {
            if(fmod(t+1e-9,1)< p.par[dt])
                {
//                    nextToDisplay(&ListB_cell, &ListT_cell, &ListFDC, &ListP_cell, NULL, &ListSC, currentLattice, t);
//                    display();
//                     glutMainLoopEvent();
                }
        }
}

// Calculation of B cells
void simulation::Calc_BC(double t,parameters &p,lattice &l, vector<vector3D>&redo_list, vector<int> &going_to_delet)
{
    for(unsigned int i = 0; i < ListB_cell.size(); i++)
        {
            B_cell* Bcell = ListB_cell.at(i);
            Bcell->clock += 1;
            ///Elena-network
            //Elena: network: Calculate TF levels of Bcell every dt.
            if((Bcell->cell_state != contact_FDC)&&(Bcell->cell_state != contact_TC)) 
            {
                double networkdt = p.par[dt];
                Bcell->calcNetwork(networkdt , 0, 0);
            }
            ///Elena-network

//Centrocytes_____________________________________________________________________________________________________
            if (Bcell->cell_type==Centrocyte)
            {
                switch (Bcell->cell_state)
                {
                    case unselected:
                    {
                        Bcell->clock +=1;
                        Bcell->Resensitize2Chemokines(p, l);
                        Bcell->BC_FDC_interaction_clock += p.par[dt];
                        
                        if (Bcell->BC_FDC_interaction_clock > p.par[collectionFDCperiod])
                        {
                            if (Bcell->retained_Ag> 100)
                            {
                                Bcell->retained_Ag=100;
                            }
                            if (Bcell->retained_Ag <= 0)
                            {
                                Bcell->cell_state = apoptosis;
                                Bcell->isResponsive2CXCL12=false;
                                Bcell->isResponsive2CXCL13=true;
                                ///Elena-network
                                EventOutput->recordEvent(Bcell, event_die, t); //Elena:  events: record event die with ID of cell to track cell history
                                ///Elena-network
                            }
                            else
                            {
                                //CCs selected by FDCs here
                                Bcell->cell_state = FDC_selected;
                                //#check @danial , where do we reset this?!

                                ///Elena-network
                                EventOutput->recordEvent(Bcell, event_FDC_selected, t); //Elena:  events: record event die with ID of cell to track cell history
                                ///Elena-network

                                Bcell->Selected_by_FDC=true;
                                Bcell->can_move=true;
                            }
                        }
                        else if (Bcell->clock>p.par[testDelay])
                        {
                            vector3D fdc_position(-1,-1,-1);
                            vector<vector3D> neighbours = l.getNeighbour_nn(Bcell->position);
                            for (unsigned int j=0; j<neighbours.size();j++)
                            {
                                if (l.insideBorders(neighbours[j]))
                                {
                                    if (l.getAgat(neighbours[j])>0.)
                                    {
                                        fdc_position=neighbours[j];
                                    }
                                }
                            }
                            if (fdc_position.X != -1)
                            {
                                bool suppress_next_interaction = false;
                                Bcell->setMyAffinity(p);
                                double binding_probability= Bcell->MyAffinity;
                                //#recheck @danial: later distinguish successful and unsuccessful FDC contacts.
                                Bcell->nFDCcontacts += 1;

                                if(random::randomDouble(1)<binding_probability)
                                {
                                    short success;
                                    if ((l.getAgat(fdc_position)>= 1.)&& (random::randomDouble(1) < double(l.getAgat(fdc_position)/ p.par[agSaturation])))
                                    {
                                        l.removeAgAt(fdc_position); //Danial: #Check  this, still remove antigen happens after binding before testing Psel !!!!!! why? #Important
                                        success=1;
                                    }
                                    else
                                    {
                                        success= 0;
                                    }
                                    if (success == 1)
                                    {
                                        Bcell->cell_state=contact_FDC;
                                        Bcell->can_move=false;

                                        ///Elena-network
                                        EventOutput->recordEvent(Bcell, event_catch_Ag, t);//Elena: events: record history output at Contact FDC cell sate
                                        ///Elena-network
                                    }
                                    else
                                    {
                                        suppress_next_interaction = true;
                                    }
                                }
                                else
                                {
                                    suppress_next_interaction = true;
                                }
                                if (suppress_next_interaction)
                                {
                                    Bcell->clock = 0;
                                }
                            }
                        }
                        
                        if (Bcell->cell_state==unselected)
                        {
                            Bcell->move(p, l, redo_list); //cehck that in contact ones don't move in move function
                        }
                        break;
                    }
                    case contact_FDC:
                    {
                        Bcell->BC_FDC_interaction_clock += p.par[dt];
                        Bcell->fdc_interaction_time_history += p.par[dt];

                        ///Elena-network
                        double interaction_time =  p.par[dt];
                        ///Elena-network

                        //Psel is to selecct only a ratio of fdc-contact cells to pick up Ag.
                        if (random::randomDouble(1) < p.par[pSel])
                        {
                            Bcell->retained_Ag += 1.0;
                            Bcell->cell_state=unselected;
                            Bcell->can_move=true;
                            Bcell->clock=0;

                            ///Elena-network
                            Bcell->calcNetwork(interaction_time, 1, 0); //bcr0 > 0 means there is interaction with FDC. Intensity of interaction is =1.
                            ///Elena-network

                            ///Elena-network
                            EventOutput->recordEvent(Bcell, event_unselected, t);//Elena: events: record history output at become unselected cell state.
                            ///Elena-network
                        }
                        Bcell->isResponsive2CXCL13=false;
                        ///Elena-network
                        //Elena: network: If NO internalized ag then calculate network with bcr = 0
                        Bcell->calcNetwork(interaction_time, 0, 0); //bcr0 = 1 means there is interaction. Intensity of interaction is to be decided.
                        ///Elena-network
                        break;
                    }
                        
                    case FDC_selected:
                    {
                        Bcell->Resensitize2Chemokines(p, l);
                        vector3D tc_position (-1,-1,-1);
                        vector <vector3D> tmp_neighbours = l.getNeighbour_nn(Bcell->position);
                        vector <vector3D> TC_neighbours;
                        for (unsigned int k=0 ;k<tmp_neighbours.size();k++)
                        {
                            if (l.insideBorders(tmp_neighbours[k]))
                            {
                                if (l.celltypeat(tmp_neighbours[k]) == TFHC)
                                {
                                    TC_neighbours.push_back(tmp_neighbours[k]);
                                }
                            }
                        }
                        if (TC_neighbours.size()>0)
                        {
                            short x = random::randomInteger(int(TC_neighbours.size()));
                            tc_position=TC_neighbours[x];
                        }
                        if (tc_position.X == -1)
                        {
                            Bcell->move(p, l, redo_list);
                        }
                        else
                        {
                            //bind CC to TC
                            Bcell->cell_state= contact_TC;
                            Bcell->can_move=false;
                            Bcell->Bc_Tc_interaction_clock=0;
                            T_cell* TC = (T_cell*) l.cellat(tc_position);
                            Bcell->interactingTC = TC;

                            ///Elena-network
                            EventOutput->recordEvent(Bcell, event_start_contact_TC, t);//Elena: events: record history output at Contact TC state.
                            ///Elena-network

                            TC->nIncontactCCs += 1;
                            TC->cell_state= TC_connected;
                            TC->interactingCC.push_back(Bcell);
                        }
                        break;
                    }
                    case contact_TC:
                    {
                        Bcell->isResponsive2CXCL13=false;
                        //duration of contact
                        Bcell->Tc_interaction_history.first += p.par[dt];
                        Bcell->Bc_Tc_interaction_clock += p.par[dt];
                        T_cell* TC = (T_cell*) Bcell->interactingTC;
                        if (TC==NULL)
                        {
                            cout<<"Accessing null Tcell from CC"<<endl;
                        }
                        else if (TC->ID!=Bcell->interactingTC->ID)
                        {
                            cout<<"Accessing wrong TC from CC"<<endl;
                        }
                        else
                        {
                            vector3D CC_neighbour = l.get_nn_directed2(TC);
                            if (CC_neighbour.X!=-1)
                            {
                                cell* cellthere = l.grid.at(CC_neighbour.X).at(CC_neighbour.Y).at(CC_neighbour.Z);
                                if (cellthere!=NULL)
                                {
                                    B_cell* neighbour_CC = (B_cell*) l.cellat(CC_neighbour);
                                    if (neighbour_CC->ID==Bcell->ID)
                                    {
                                        //TC and CC are face2face
                                        Bcell->TCsignalDuration += p.par[dt];
                                        //duration of signal
                                        Bcell->Tc_interaction_history.second += p.par[dt];

                                        ///Elena-network
                                        //Elena: network: Calculate TF levels during TFHC-Bcell interaction.
                                        double interaction_time = p.par[dt];
                                        // Bcell->calcNetwork(interaction_time, 0, 50); // Elena: network: fixed cd40 = 50.
                                        Bcell->calcNetwork(interaction_time, 0, Bcell->MyAffinity*50); //Elena: network: dynamic cd40 proportional to affinity of Bcell [0,1] corrected to [0,50]
                                        ///Elena-network

                                        ///Elena-network
                                        if(Bcell->TC_signal_start)//Elena: events: to record the start of TC_signal
                                        {
                                        EventOutput->recordEvent(Bcell, event_start_signaling_TC, t);//Elena: events: record history output at start of signalling TC state.
                                        Bcell->TC_signal_start = false;
                         		}
                                        ///Elena-network
                                    }
                                }
                                else
                                {
                                    cout<<"time= "<<t<<" There is no CC in the directed position, TC pos="<<TC->position.print()<<" CC pos="<<Bcell->position.print()<<" Directedpos="<<CC_neighbour.print()<<" CCID="<<Bcell->ID<<" TCID="<<TC->ID<<" cell_state="<<Bcell->cell_state<<endl;
                                }
                            }
                            else
                            {
                                //  cout<<"The CC neighbour is negative -1"<<endl;
                            }
                        }
                        if ( Bcell->TCsignalDuration > p.par[tcRescueTime])
                        {
                            //CC_TC Selection
                            if (Bcell->retained_Ag>100)
                            {
                                Bcell->retained_Ag=100;
                            }
                            Bcell->timeleft2recycle(p);
                            double pMHC = Bcell->retained_Ag;
                            double ag_factor = pow(pMHC,p.par[pMHCdepHill]);
                            //Record selected CC mutation frequencies
                            //Number of pmhc dpendent divisions
                            Bcell->pMHC_dependent_number_of_divisions= p.par[pMHCdepMin] + (p.par[pMHCdepMax] - p.par[pMHCdepMin]) * ag_factor / (ag_factor + pow(p.par[pMHCdepK],p.par[pMHCdepHill]));
                            double ndivtmp = 2.0; //total number of divisions hyphasma temporary
                            if(Bcell->pMHC_dependent_number_of_divisions>=0)
                            {
                                ndivtmp=Bcell->pMHC_dependent_number_of_divisions;
                            }
                            Bcell->nDivisions2do= int(ndivtmp);
                            ndivtmp -= double(Bcell->nDivisions2do);
                            if (random::randomDouble(1) < ndivtmp)
                            {
                                ++Bcell->nDivisions2do;
                            }
                            if (Bcell->nDivisions2do>12)
                            {
                                Bcell->nDivisions2do=12;
                            }
                            
                            Bcell->cell_state= TC_selected;
                            
                            //#Check @danial


                            ///Elena-network
                            //Elena: network: events: stop signal
                            EventOutput->recordEvent(Bcell, event_stop_signaling_TC, t);//Elena: events: record history output at stop signalling TC event
                            Bcell->TC_signal_start = true;
                            ///Elena-network

                            Bcell->Selected_by_TC=true;

                            Bcell->can_move=true;
                            Bcell->TC_selected_clock=0.0;

                        }
                        if (not(Bcell->cell_state == TC_selected) && Bcell->Bc_Tc_interaction_clock > p.par[tcTime] )
                        {
                            //apoptotic cells mutation frequency
                            Bcell->cell_state = apoptosis;
                            Bcell->can_move=true;
                            Bcell->isResponsive2CXCL13=true;
                            Bcell->isResponsive2CXCL12=false;

                            ///Elena-network  
                            //Elena: network: events: die
                            EventOutput->recordEvent(Bcell, event_die, t); //Elena: events: Record event die with ID of cell to track cell history
                            ///Elena-network  
                        }
                        
                        if (Bcell->cell_state == TC_selected || Bcell->cell_state == apoptosis )
                        {
                            if (not(TC==NULL))
                            {
                                TC->liberateCC_TC(Bcell); //hyphasma only liberate TC, cc has been changed to apop

                                ///Elena-network 
                                //Elena: network: events: Liberate bc and TC
                               EventOutput->recordEvent(Bcell, event_stop_contact_TC, t); //Elena:  events: record event die with ID of cell to track cell history
                                ///Elena-network 

                            }
                            Bcell->can_move=true;

                        }
                        break;
                    }
                    case TC_selected:
                    {
                        Bcell->isResponsive2CXCL13=false;
                        Bcell->TC_selected_clock += p.par[dt];
                        if (Bcell->TC_selected_clock > Bcell->Recycling_delay)
                        {
                            if (random::randomDouble(1)< p.par[p_dif]) ////recheck (there is another one, are they same?)
                            {
//Here cells recycle
                                Bcell->isResponsive2CXCL12=true;
                                Bcell->isResponsive2CXCL13=false;
                                Bcell->cell_type=Centroblast;
                                
                                //#Rechcek @Danial:We set this here, while we have to decide after division, in practice it is the same since we set iamhigh for daughter B-cell to false and this cell diffs to plasma.
                                if (Bcell->retained_Ag > 0.)
                                {
                                    Bcell->IamHighAg = true;
                                }
                                
                                //recycling cells mutation frequency
                                //#Check @danial Reset after recording everything, that is now after division happens
                                Bcell->cell_state=recycled;

                                ///Elena-network 
				//Elena: network: events: Liberate bc and TC
                                EventOutput->recordEvent(Bcell, event_recycling, t); //Elena:  events: record event die with ID of cell to track cell history
                                Bcell->nRecyclings = Bcell->nRecyclings + 1;
                                ///Elena-network 

                                Bcell->setMyAffinity(p);
                                
                                //#recheck
//                                Bcell->myBCR.pMut= p.par[pmutAfterStartMut] + double ((0. - p.par[pmutAfterStartMut]) * pow(Bcell->MyAffinity,p.par[pmutAffinityExponent]));
                                
                                Bcell->cyclestate=cycle_G1;
                                Bcell->transmit_CCdelay2cycle(p);
                                if(Bcell->nDivisions2do <= 0)
                                {
                                    Bcell->cyclestate = cycle_G0;
                                }
                                Bcell->Recycling_delay=0;
                            }
                            else
                            {
                                Bcell->move(p, l, redo_list);
                            }
                        }
                        else
                        {
                            Bcell->move(p, l, redo_list);
                        }
                        break;
                    }
                        
                    case apoptosis:
                    {
                        if(random::randomDouble(1)<p.par[macrophage])
                        {
                            //#record_event
                            currentOutput->close_event(Bcell, sim_output,t);
                            currentOutput->write_event( Bcell, sim_output);
                            going_to_delet.push_back(Bcell->ID);
                        }
                        else
                        {
                            Bcell->isResponsive2CXCL12=false;
                            Bcell->isResponsive2CXCL13=true;
                            Bcell->Resensitize2Chemokines(p,l);
                            Bcell->move(p, l, redo_list);
                        }
                        break;
                    }
                    default:
                        break;
                }
            }
//___________________________________________________________________________________________________________
            
//Centroblasts_______________________________________________________________________________________________
            if (Bcell->cell_type==Centroblast)
                {
                if(Bcell->cell_state==apoptosis) {
                    currentOutput->close_event(Bcell, sim_output,t);
                    currentOutput->write_event( Bcell, sim_output);
                    going_to_delet.push_back(Bcell->ID);
                } else {
                        // Increase cell cycle time
                        Bcell->cycle_state_time += p.par[dt];
                        // Resensitize
                        Bcell->Resensitize2Chemokines(p, l);
                        // Switch cycle state
                        if ((Bcell->cycle_state_time >= Bcell->time_of_cycle_state_switch) && (Bcell->cyclestate < cycle_Divide))
                            {
                                Bcell->ContinueCellCycle(p);
                                Bcell->cycle_state_time=0.0;
                            }
                        // #Sequential After finishing the cell cycle B cell divides in same time step
                        if (Bcell->cyclestate==cycle_Divide)
                            {
                                 ///Elena-network
                                 EventOutput->recordEvent(Bcell, event_divide, t); //Elena: events: record history output at division state.
                                 ///Elena-network


                                 ///Elena-network
                                 Bcell->proliferate(p,l,t,ListB_cell,*currentOutput,*this);  ///Commented because Elena inserts the daughterBcell thing, uncommented by me

                                 ///RRR B_cell* daughterBcell = Bcell->proliferate(p,l,t,ListB_cell,*currentOutput,*this);
                                 ///RRR if(daughterBcell)
                                 ///RRR {
                                 ///RRR EventOutput->recordEvent(daughterBcell, event_born, t); //Elena: events: record history output at born sate.
                                 EventOutput->recordEvent(Bcell, event_born, t); //Elena: events: record history output at born state.
                                 ///RRR }
                                 ///Elena-network

                            }
                        else if (Bcell->cyclestate==cycle_G0)
                            {
                                // p.par[p_dif] == 0.02

                                if ( random::randomDouble(1) < p.par[p_dif]) //recheck (there is another one, are they same?)
                                    {
                                        ///Elena-network
                                        // Elena: network: 1- Uncomment if output based on Iamhigh rule.
                                        //                            if ( Bcell->retained_Ag > 0. && Bcell->IamHighAg)
                                        //                                                       {
                                        //                                                           Bcell->cell_type=Plasmacell;
                                        //                                                           Bcell->cell_state=Plasma_in_GC;//Elena: Change cell state. Important for output!
                                        //                                                       }

                                        //Elena: network:  3- Uncomennt if Plasma cell depends on BLIMP1high and Memory cell output based on I amHighAg and not BLIMP1 high.
                                        ///Elena-network

                                        ///RRR elimino from here onwards:
                                    if ( Bcell->retained_Ag > 0.)
                                    {
                                        if(Bcell->BLIMP1 >= 8) //Elena: network parameter defined in parameterfile.
                                        {
                                            //cout<<"Plasma"<<endl;
                                            Bcell->cell_type=Plasmacell;
                                            Bcell->cell_state=Plasma_in_GC;//Elena: Change cell state. Important for output!
                                            EventOutput->recordEvent(Bcell, event_become_plasma, t);
                                        }
                                       else
                                        {
                                            if (Bcell->IamHighAg)
                                            {
                                                //cout<<"Memory"<<endl;
                                                Bcell->cell_type=Memorycell;
                                                Bcell->cell_state=Memory_in_GC;//Elena: Change cell state. Important for output!
                                                EventOutput->recordEvent(Bcell, event_become_memory, t);
                                            }
                                        }
                                     ///Elena-network

                                        ///RRR until here

                                        ///Elena -network adapted
//                                       if ( Bcell->retained_Ag > 0 && Bcell->BLIMP1 >= 8 ) {  /// p.par[BLIMP1th] = 8; in this version we dont have these params in par file
//                                           Bcell->cell_type=Plasmacell;
//                                           Bcell->cell_state=Plasma_in_GC;
//                                           EventOutput->recordEvent(Bcell, event_become_plasma, t);
//                                       } else if (Bcell->retained_Ag > 0 && Bcell->IamHighAg ) {
//                                           Bcell->cell_type=Memorycell;
//                                           Bcell->cell_state=Memory_in_GC;
//                                           EventOutput->recordEvent(Bcell, event_become_memory, t);
                                        ///Elena -network adapted
                                       } else {
                                                //CCs created here
                                                Bcell->isResponsive2CXCL12=false;
                                                Bcell->isResponsive2CXCL13=true;
                                                Bcell->cell_type=Centrocyte;
                                                Bcell->cell_state=unselected;
                                                Bcell->Selected_by_FDC=false;
                                                Bcell->Selected_by_TC=false;
                                                //#check @danial: in a scenario that Antigens remain inside the cell, take care of this
                                                Bcell->nFDCcontacts=0;
                                                Bcell->retained_Ag=0;
                                                Bcell->clock =0 ; //not sure
                                                ///Elena-network
                                                EventOutput->recordEvent(Bcell, event_unselected, t);
                                                ///Elena-network

                                            }
                                    }

                                }
                        if ((Bcell->cell_type==Centroblast) && (Bcell->cyclestate!=cycle_M))
                            {
                                Bcell->move(p, l, redo_list);
                            }
                    }
//______________________________________________________________________________________
                }
        }
}

void simulation::clean_dead_cells(lattice &l)
{
    //Danial: This is not sustainable, can be improved by using template functions for erase (since c++2a)
    long tmp_size_1= ListB_cell.size();
    for (unsigned int j=0; j < going_to_delet.size();j++)
    {
        if (ListB_cell.size()==1)
        {
            if (going_to_delet[j]==ListB_cell.at(0)->ID)
            {
                l.removecellat(ListB_cell.at(0)->position);
                delete ListB_cell[0];
                ListB_cell.pop_back();
            }
            if (ListB_cell.size()>0)
            {
                cout<<"Error in deleting dead B cells (1)."<<endl;
                exit(1);
            }
        }
        else if (ListB_cell.size()> 1)
        {
            for (unsigned int ks=0; ks<ListB_cell.size();ks++)
            {
                long tmp_size_2= ListB_cell.size();
                if (ListB_cell.at(ks)!=NULL)
                {
                    if (ListB_cell.at(ks)->ID==going_to_delet[j])
                    {
                        l.removecellat(ListB_cell.at(ks)->position);
                        delete ListB_cell.at(ks);
                        ListB_cell.at(ks)=NULL;
                        tmp_size_2--;
                    }
                }
            }
        }
        // Alternative deleting method
        //  ListB_cell.erase(remove_if(ListB_cell.begin(), ListB_cell.end(),[&ID,&pos,&tmp_ID](const B_cell* x  ) { if (x->ID==ID){pos=x->position;} return (x->ID == ID);}), ListB_cell.end());
    }

    ListB_cell.erase(remove_if(ListB_cell.begin(), ListB_cell.end(),[](const B_cell* x  ) { return (x==NULL);}), ListB_cell.end());
    if (tmp_size_1 != (going_to_delet.size()+ListB_cell.size()))
    {
        cout<<"Error in deleting dead B cells (2)."<<endl;
        exit(1);
    }
    going_to_delet.clear();
}


void simulation::Calc_Out(double t,parameters &p,lattice &l, vector<vector3D>&redo_list)
{


    for(unsigned int i = 0; i < ListP_cell.size(); i++)
    {
        Plasma_cell* Plasma = ListP_cell.at(i);
        if (not(Plasma->cell_state==Plasma_Out))
        {
        if(l.is_at_border(Plasma->position))
            {
                l.removecellat(Plasma->position);
                Plasma->cell_state=Plasma_Out;
            }
            else
            {
                Plasma->move(p, l, redo_list);
            }
        }
    }
    
//
    //_______________________________________________________________________________________________________________
    ///Elena: Memory output: remove memory cells if next to border! As for memory cells
    for(unsigned int i = 0; i < ListM_cell.size(); i++)
    {
        Memory_cell* Memory = ListM_cell.at(i);

        if (not(Memory->cell_state==Plasma_Out))
        {
        if(l.is_at_border(Memory->position))
            {
                l.removecellat(Memory->position);
                Memory->cell_state=Plasma_Out;
            }
            else
            {
                Memory->move(p, l, redo_list);
            }
        }
    }
     ///Elena
}
void simulation::transfer_plasma_from_Bcell_list(double t,parameters &p,lattice &l, vector<vector3D>&redo_list){
    
    /*Danial: This function is to transfer B cells which differentiate to plasma cell from B cell list to plasma cell list
     We could do this in calc_out but that causes the problem of progressing the simulation for some B cells twice in a row*/
    
    //Output cells___________________________________________________________________________________________________
    
    for(unsigned int i = 0; i < ListB_cell.size(); i++)
    {
        B_cell* Bcell = ListB_cell.at(i);
        if (Bcell->cell_type==Plasmacell)
        {
            Plasma_cell* new_Plasma = new Plasma_cell(p, Bcell);
            new_Plasma->birth_time=t;
            new_Plasma->ID=Bcell->ID;
            new_Plasma->MID=Bcell->MID;
            l.removecellat(Bcell->position);
            l.putcellat(new_Plasma);
            ListP_cell.push_back(new_Plasma);
            delete ListB_cell.at(i);
            ListB_cell.at(i)=NULL;
        }

        ///Elena: Memory output: Put memory cells in Memory B cell list for output and exit of GC!
        if (Bcell->cell_type==Memorycell)
               {
                   Memory_cell* new_Memory = new Memory_cell(p, Bcell);
                   new_Memory->birth_time=t;
                   new_Memory->MID=Bcell->MID;
                   l.removecellat(Bcell->position);
                   l.putcellat(new_Memory);
                   ListM_cell.push_back(new_Memory);
                   delete ListB_cell.at(i);
                   ListB_cell.at(i)=NULL;
               }
        ///Elena
    }
    
    ListB_cell.erase(remove_if(ListB_cell.begin(), ListB_cell.end(),[](const B_cell* x  ) { return (x==NULL);}), ListB_cell.end());
    
    for(unsigned int i = 0; i < ListB_cell.size(); i++)
    {
        B_cell* Bcell = ListB_cell.at(i);
        
        if (Bcell==NULL)
        {
            cout<<"Erorr, NULL BC"<<endl;
        }
    }
    
    for(unsigned int i = 0; i < ListP_cell.size(); i++)
    {
        Plasma_cell* Plasma = ListP_cell.at(i);
        
        if (Plasma==NULL)
        {
            cout<<"Erorr, NULL PC"<<endl;
        }
    }
    
    ///Elena: Memory output:
    for(unsigned int i = 0; i < ListM_cell.size(); i++)
    {
        Memory_cell* Memory = ListM_cell.at(i);

        if (Memory==NULL)
        {
            cout<<"Erorr, NULL MC"<<endl;
        }
    }
     ///Elena
}
