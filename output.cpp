#include "output.h"
#include <set>
#include <numeric>
#include "mafalda.h"
#include "cell.h"
#include <experimental/filesystem> ///RRR
#include "SHM.h"
//#include <boost/algorithm/string/join.hpp>
using namespace std;





// ministats
vector<ministat> Bcell_counts;  // number of Bcells
vector<ministat> Plasma_counts;

// Affinities
// This is for recording (average) affinities of B cells versus time, for plasma
// cells there is a different file of records
ministat Bcell_affinity;

output::output(string _parfname) : parfname(_parfname) { initialize_fileds(); }
output::~output() {}

string output::currentDateTime() {
  time_t now = time(0);
  struct tm tstruct;
  char buf[80];
  tstruct = *localtime(&now);
  strftime(buf, sizeof(buf), "%y-%m-%d_%H-%M-%S", &tstruct);
  return buf;
}

void output::createFolder(string folderName) {


#ifdef __WIN32__
  const char *p = tmp.str();
  const WCHAR *pwcsName;
  int nChars = MultiByteToWideChar(CP_ACP, 0, p, -1, NULL, 0);
  pwcsName = new WCHAR[nChars];
  MultiByteToWideChar(CP_ACP, 0, p, -1, (LPWSTR)pwcsName, nChars);
  CreateDirectory(pwcsName, NULL);
  delete[] pwcsName;
#endif

//#Recheck
#ifdef __linux__
  //#Danial: Fully changed
  stringstream tmp;
  tmp <<"/" << "sim_" << Output_ID << "_" << currentDateTime();
        folderName = folderName + tmp.str();
       output_path = folderName;

      if (!(std::experimental::filesystem::create_directory(output_path))) { ///RRR tmp.str() Y otracosa
        std::cerr << "Error creating directory " << std::endl; ///RRR
      }


#endif

#ifdef __APPLE__
    //#Danial: Fully changed
    stringstream tmp;
    tmp <<"/" << "sim_" << Output_ID << "_" << currentDateTime();
          folderName = folderName + tmp.str();
         output_path = folderName;

        if (!(std::experimental::filesystem::create_directory(output_path))) { ///RRR tmp.str() Y otracosa
          std::cerr << "Error creating directory " << std::endl; ///RRR
        }
      char cstr [folderName.size()+1];
      strcpy(cstr,folderName.c_str());
      mkdir(cstr, 0777);
#endif
}

void output::createFolder(string folderName, parameters &p) {
  //#Recheck
  // Output folder
  stringstream tmp;
  tmp << base_Path << "/" << p.parameter_file_name << "sim_" << Output_ID
      << "_" << currentDateTime();
  output_path = tmp.str();
//#Recheck
#ifdef _WIN32
  const char *p = tmp.str();
  const WCHAR *pwcsName;
  int nChars = MultiByteToWideChar(CP_ACP, 0, p, -1, NULL, 0);
  pwcsName = new WCHAR[nChars];
  MultiByteToWideChar(CP_ACP, 0, p, -1, (LPWSTR)pwcsName, nChars);
  CreateDirectory(pwcsName, NULL);
  delete[] pwcsName;
#endif
    

//#Recheck
#ifdef __linux__
  if (!(std::experimental::filesystem::create_directory(tmp.str()))) { ///RRR tmp.str() Y otracosa
    std::cerr << "Error creating directory " << std::endl; ///RRR
  }
#endif

#ifdef __APPLE__
  char cstr[tmp.str().size() + 1];
  strcpy(cstr, tmp.str().c_str());
  mkdir(cstr, 0777);
#endif
}

void output::initialize_fileds() {
  // 0-9 -> CC states (7) + CB (2)
  for (int i = 0; i < 8; i++) {
    Bcell_counts.push_back(ministat());  // number of Bcells
  }
}

void output::clear_fileds() {
  for (int i = 0; i < 8; i++) {
    Bcell_counts.clear();  //.at(i).clear_ministat();   //number of Bcells
  }
  for (int i = 0; i < 8; i++) {
    Bcell_counts.push_back(ministat());  // number of Bcells
  }

  // Affinities
  Bcell_affinity.clear_ministat();
}

// Take fields from simulation into master observer variable to create file.
void output::record_output_time_step(double currentTime, simulation &currentSim,
                                     parameters &p) {
  // This function records data about B cells population and affinity versus
  // time, this does not include Plasma cells.

  // Bcell data

    int CB = 0;
    int CC = 0;
    int apop =0;
    vector<int> CC_IDs;
  for (int i = 0; i < currentSim.ListB_cell.size(); i++) {
    B_cell *Bcell = currentSim.ListB_cell.at(i);
    // Check integrity of the list

    if (Bcell->cell_state > 7) {
      cout << "Error, wrong cell in BC list," << Bcell->cell_state << endl;
      exit(1);
    }
    Bcell_counts[Bcell->cell_state].add(1);
    if (Bcell->cell_state == 7) {
        apop = apop +1;
    } else if ( Bcell->cell_type == 4) {
        CB = CB +1;
    } else if ( Bcell->cell_type == 5) {
        CC = CC +1;
        CC_IDs.push_back(Bcell->ID);
    }

    // Bcell Affinity
    Bcell->setMyAffinity(p);
    Bcell_affinity.add(Bcell->MyAffinity);
  }

  FILE *Bcell_time_data;
  string folder1 = output_path + "/Bcell_time.csv";
  char *s1 = const_cast<char *>(folder1.c_str());
  ;
  Bcell_time_data = fopen(s1, "a");
  static bool tmp = true;
  if (tmp) {
    tmp = false;
    fprintf(Bcell_time_data, "%s",
            "time,founder,unselected,contact_FDC,FDC_selected,contact_TC,TC_"
            "Selected_by_TC,recycled,apoptosis,affinity,std_affinity,CC_IDs,CB,CC,Pr,Apoptotic\n");
  }

  fprintf(Bcell_time_data, "%f,", currentTime);
  for (int i = 0; i < 8; i++) {
    fprintf(Bcell_time_data, "%f,", Bcell_counts[i].sum);
  }

  fprintf(Bcell_time_data, "%.16G,%.16G,", Bcell_affinity.average(),
          Bcell_affinity.stddev());
  fprintf(Bcell_time_data, "%s,",
          "Prueba");
  fprintf(Bcell_time_data, "%d,%d,%s,%d\n", CB,CC,"Pr",apop);
  fclose(Bcell_time_data);  //#Recheck take care of bins in gle file
  clear_fileds();

  std::ostringstream oss;

  if (!CC_IDs.empty())
  {
    // Convert all but the last element to avoid a trailing ","
    std::copy(CC_IDs.begin(), CC_IDs.end()-1,
        std::ostream_iterator<int>(oss, ";"));

    // Now add the last element with no delimiter
    oss << CC_IDs.back();
  }
  std::string oss_s = oss.str();
  oss.str(std::string());

  FILE *CC_time_data;
  string folderCC = output_path + "/CC_time.csv";
  char *s1CC = const_cast<char *>(folderCC.c_str());
  ;
  CC_time_data = fopen(s1CC, "a");
  static bool tmpCC = true;
  if (tmpCC) {
    tmpCC = false;
    fprintf(CC_time_data, "%s",
            "time,CC_IDs\n");
  }

  fprintf(CC_time_data, "%f,", currentTime);

  fprintf(CC_time_data, "%s\n", oss_s.c_str());
  fclose(CC_time_data);  //#Recheck take care of bins in gle file





  // Bcell seq data
  vector<string> fastas_totales;

      for (const auto & entry : std::experimental::filesystem::directory_iterator(fasta_path)) {
          if (!entry.path().extension().compare(".fasta")) {
               fastas_totales.push_back(entry.path());
          }

      }

  string sequence;

  for (int indiv = 0; indiv < fastas_totales.size(); indiv++) {
      sequence= std::experimental::filesystem::path(fastas_totales[indiv]).filename();
      CB = 0;
      CC = 0;
      apop =0;

      for (int i = 0; i < currentSim.ListB_cell.size(); i++) {
        B_cell *Bcell2 = currentSim.ListB_cell.at(i);
        // Check integrity of the list

        if ((sequence.compare(Bcell2->germline_name)) == 0 ) {
            if (Bcell2->cell_state > 7) {
              cout << "Error, wrong cell in BC list," << Bcell2->cell_state << endl;
              exit(1);
            }
            Bcell_counts[Bcell2->cell_state].add(1);
            if (Bcell2->cell_state == 7) {
                apop = apop +1;
            } else if ( Bcell2->cell_type == 4) {
                CB = CB +1;
            } else if ( Bcell2->cell_type == 5) {
                CC = CC +1;
            }

            // Bcell Affinity
            Bcell2->setMyAffinity(p);
            Bcell_affinity.add(Bcell2->MyAffinity);
        }


      }

      FILE *Bcell_seq_time_data;
      string folder2 = output_path + "/Bcell_time_seq.csv";
      char *s2 = const_cast<char *>(folder2.c_str());

      Bcell_seq_time_data = fopen(s2, "a");
      static bool tmp2 = true;
      if (tmp2) {
        tmp2 = false;
        fprintf(Bcell_seq_time_data, "%s",
                "Fasta, time,founder,unselected,contact_FDC,FDC_selected,contact_TC,TC_"
                "Selected_by_TC,recycled,apoptosis,affinity,std_affinity,Prueba,CB,CC,Pr,Apoptotic\n");
      }
      fprintf(Bcell_seq_time_data, "%s,", sequence.c_str());
      fprintf(Bcell_seq_time_data, "%f,", currentTime);
      for (int i = 0; i < 8; i++) {
        fprintf(Bcell_seq_time_data, "%f,", Bcell_counts[i].sum);
      }

      fprintf(Bcell_seq_time_data, "%.16G,%.16G,", Bcell_affinity.average(),
              Bcell_affinity.stddev());
      fprintf(Bcell_seq_time_data, "%s,",
              "Prueba");
      fprintf(Bcell_seq_time_data, "%d,%d,%s,%d\n", CB,CC,"Pr",apop);
      fclose(Bcell_seq_time_data);  //#Recheck take care of bins in gle file
      clear_fileds();
  }
















}

void output::write_event(cell *Cellx, stringstream &sim_output) {
  sim_output << Cellx->event.str() << endl;
}

void output::write_event_2file(stringstream &sim_output) {
  FILE *event_data;
  string folder1 = output_path + "/event_data.csv";
  event_data = fopen(folder1.c_str(), "a");

  static bool tmp = false;

  if (not(tmp)) {
    fprintf(event_data, "%s",
            "ID,Born_time,MID,States,Affinity,N_of_Ags,N_of_divisions,N_of_"
            "Mutations,delta_aff,FDC_interaction_nums,FDC_interaction_time_avg,"
            "TC_interaction_time,TC_signaling_time,FDC_selected,Selected_by_TC, Original_fasta, Sequence, Protein sequence, Status, Replacement_mutations, Silent_Mutations, Relevant_NGly_sites, All_NGly_sites, Death_time, AA_regions, NT_REGIONS, Blacklists, BLIMP1, BCL6, IRF4, BLIMP1_0, BCL6_0, IRF4_0\n");
    tmp = true;
  }

  fprintf(event_data, "%s", sim_output.str().c_str());
  fclose(event_data);  //#Recheck take care of bins in gle file
}
// for B cells
//cell states: 0-founder,
//1-unselected,
//2-contact_FDC,
//3-FDC_selected,
//4-contact_TC,
//5-TC_selected,
//6-recycled,
//7-apoptosis,
//8-TC_free,
//9-TC_connected,
//10-Plasma_Out,
//11-Plasma_in_GC,
//12-cell_state_counter
void output::close_event(B_cell *Cellx, stringstream &sim_output, double time) {
  Cellx->event << Cellx->cell_state << "," << Cellx->MyAffinity << ","
               << Cellx->retained_Ag << "," << Cellx->total_number_of_divisions
               << "," << Cellx->myBCR.nMutFromGermline << ","
               << Cellx->delta_Affinity << ",";
  Cellx->event << Cellx->nFDCcontacts << ",";

  if (Cellx->nFDCcontacts == 0) {
    Cellx->event << Cellx->fdc_interaction_time_history << ",";
  } else {
    Cellx->event << double(Cellx->fdc_interaction_time_history /
                           double(Cellx->nFDCcontacts))
                 << ",";
  }

  Cellx->event << Cellx->Tc_interaction_history.first << ","
               << Cellx->Tc_interaction_history.second << ",";

  std::stringstream ss;
  for(size_t i = 0; i < Cellx->NGly_mutated_sites.size(); ++i)
  {
    if(i != 0)
      ss << "-";
    ss << Cellx->NGly_mutated_sites[i];
  }
  std::string s = ss.str();
  ss.str(std::string());

  std::stringstream sts;
  for(size_t i = 0; i < Cellx->All_NGly_mutated_sites.size(); ++i)
  {
    if(i != 0)
      sts << "-";
    sts << Cellx->All_NGly_mutated_sites [i];
  }
  std::string st = sts.str();
  sts.str(std::string());

  std::stringstream blk;
  for(size_t i = 0; i < Cellx->blacklist_nums.size(); ++i)
  {
    if(i != 0)
      blk << "-";
    blk << Cellx->blacklist_nums [i];
  }
  std::string blks = blk.str();
  blk.str(std::string());

  Cellx->event << Cellx->Selected_by_FDC << "," << Cellx->Selected_by_TC << ","<< Cellx->germline_name << ","  << Cellx->sequence << "," << Cellx->protein_sequence << ","<< cellToString(Cellx->cell_type) << "," << mapToString(Cellx->MUTATIONS) << "," << mapToString(Cellx->SILENT_MUTATIONS) << "," << s  << "," << st  << ","<<  time << ","<<   mapToString2(Cellx->AA_REGIONS).c_str() << ","<<   mapToString2(Cellx->REGIONS).c_str() << ","<< blks<< ","<< Cellx->BLIMP1 << ","<< Cellx->BCL6<< ","<< Cellx->IRF4<< ","<< Cellx->BLIMP1_0 << ","<< Cellx->BCL6_0<< ","<< Cellx->IRF4_0;
}

void output::Plasma_output(double currentTime, simulation &currentSim,
                           parameters &p) {
  FILE *Plasma_cells_data;
    
  string folder1 = output_path + "/event_data.csv";
  Plasma_cells_data = fopen(folder1.c_str(), "a");
 // fprintf(Plasma_cells_data, "%s",
 //          "ID,Born_time,MID,States,Affinity,N_of_Ags,N_of_divisions,N_of_Mutations,"
  //         "delta_aff,FDC_interaction_nums,FDC_interaction_time_total,TC_"
 //          "interaction_time,TC_signaling_time,FDC_selected,Selected_by_TC,Sequence,Protein_Sequence, Status, Mutations, Relevant_NGly_sites, Death_time\n"); ///RRR
  // Plasma data

      std::stringstream ss;
      std::string s;

  for (int j = 0; j < currentSim.ListP_cell.size(); j++) {
    Plasma_cell *Plasma = currentSim.ListP_cell.at(j);



    for(size_t i = 0; i < Plasma->NGly_mutated_sites.size(); ++i)
    {
      if(i != 0)
        ss << "-";
      ss << Plasma->NGly_mutated_sites[i];
    }
    s = ss.str();
    ss.str(std::string());

    std::stringstream sts;
    for(size_t i = 0; i < Plasma->All_NGly_mutated_sites.size(); ++i)
    {
      if(i != 0)
        sts << "-";
      sts << Plasma->All_NGly_mutated_sites [i];
    }
    std::string st = sts.str();
    sts.str(std::string());


    std::stringstream blk;
    for(size_t i = 0; i < Plasma->blacklist_nums.size(); ++i)
    {
      if(i != 0)
        blk << "-";
      blk << Plasma->blacklist_nums [i];
    }
    std::string blks = blk.str();
    blk.str(std::string());


      fprintf(Plasma_cells_data, "%d,%f,%d,%d,%.16G,%f,%d,%d,%.16G,%d,%f,%f,%f,%d,%d,%s,%s,%s,%s,%s,%s,%s,%s,%f,%s,%s,%s,%f,%f,%f,%f,%f,%f\n",
            Plasma->ID, Plasma->birth_time, Plasma->MID,Plasma->cell_state, Plasma->MyAffinity,
            Plasma->retained_Ag, Plasma->total_number_of_divisions,
            Plasma->myBCR.nMutFromGermline, Plasma->delta_Affinity, Plasma->nFDCcontacts,
            Plasma->fdc_interaction_time_history,
            Plasma->Tc_interaction_history.first,
            Plasma->Tc_interaction_history.second, Plasma->Selected_by_FDC,Plasma->Selected_by_TC,Plasma->germline_name.c_str(), Plasma->sequence.c_str(), Plasma->protein_sequence.c_str(), "Plasma", mapToString(Plasma->MUTATIONS).c_str(), mapToString(Plasma->SILENT_MUTATIONS).c_str(), s.c_str(), st.c_str(),currentTime, mapToString2(Plasma->AA_REGIONS).c_str(), mapToString2(Plasma->REGIONS).c_str(), blks.c_str(),Plasma->BLIMP1,Plasma->BCL6,Plasma->IRF4,Plasma->BLIMP1_0,Plasma->BCL6_0,Plasma->IRF4_0);  ///RRR added sequence in fprintf (2 changes)
  }
  fclose(Plasma_cells_data);  //#Recheck take care of bins in gle file

  fstream database;
  database.open(output_path + "/Sequence_db.csv",  fstream::out);
  fstream myStream;
  for(auto& kv : Seq_affinity) {
    database << kv.first  <<  ',' << kv.second[0]  <<  ','  << kv.second[1]  <<  ','  << kv.second[2]  <<  ','  << kv.second[3]  <<  ',' << "\n";
  }
  database.close();


  fstream database2;
  database2.open(output_path + "/Restricted_db.csv",  fstream::out);
  fstream myStream2;
  for(auto& kv : Restricted_db) {
    for(auto& kv2 : kv.second) {
      for(auto& kv3 : kv2.second) {
          for(auto& kv4 : kv3.second) {
            database2 << kv.first <<  ',' << kv2.first <<  ',' << kv3.first <<  ',' << kv4.first <<  ','  << kv4.second <<  ',' << "\n";
            ////if (kv4.second.compare("")==0) {
               ////cout << "FALLOnoenEscribir"<<endl;
            ////}
          }      }
    }

  }
  database2.close();

  fstream databaseCDR;
  databaseCDR.open("bcinflow09/CDRs_db.csv",  fstream::out);
  fstream myStreamCDR;
  for(auto& kv : CDRs_affinity) {
    databaseCDR << kv.first[0]  <<  ',' << kv.first[1]  <<  ',' << kv.first[2]  <<  ',' << kv.second[0]  <<  ','  << kv.second[1]  <<  ','  << kv.second[2]  <<  ','  <<            kv.second[3]  <<  ',' << "\n";
  }
  databaseCDR.close();

  fstream databaseFWR;
  databaseFWR.open("bcinflow09/FWRs_db.csv",  fstream::out);
  fstream myStreamFWR;
  for(auto& kv : FWRs_db) {
    databaseFWR << kv.first[0]  <<  ',' << kv.first[1]  <<  ','  << kv.second << "\n";
  }
  databaseFWR.close();
}


 ///Elena
void output::Memory_output(double currentTime, simulation &currentSim,
                           parameters &p) {
  // Danial: This function only writes down the data of plasma cells at the end
  // of the simulation.
  /*The order is
   1-Time of production(differetiation)
   2-ID
   3-ID of mother B-cell
   4-Total number of divisions
   5-Total amount of Ag
   6-Affinity
   7-Total number of mutations
   */

  FILE *Memory_cells_data;
  string folder2 = output_path + "/event_data.csv";

  //    char *ss1 = const_cast<char*>(folder1.c_str());;

  Memory_cells_data = fopen(folder2.c_str(), "a");



  std::stringstream ss;
  std::string s;
  for (int j = 0; j < currentSim.ListM_cell.size(); j++) {
    Memory_cell *Memory = currentSim.ListM_cell.at(j);

    for(size_t i = 0; i < Memory->NGly_mutated_sites.size(); ++i)
    {
      if(i != 0)
        ss << "-";
      ss << Memory->NGly_mutated_sites[i];
    }
    s = ss.str();
    ss.str(std::string());

    std::stringstream sts;
    for(size_t i = 0; i < Memory->All_NGly_mutated_sites.size(); ++i)
    {
      if(i != 0)
        sts << "-";
      sts << Memory->All_NGly_mutated_sites [i];
    }
    std::string st = sts.str();
    sts.str(std::string());


    std::stringstream blk;
    for(size_t i = 0; i < Memory->blacklist_nums.size(); ++i)
    {
      if(i != 0)
        blk << "-";
      blk << Memory->blacklist_nums [i];
    }
    std::string blks = blk.str();
    blk.str(std::string());


      fprintf(Memory_cells_data, "%d,%f,%d,%d,%.16G,%f,%d,%d,%.16G,%d,%f,%f,%f,%d,%d,%s,%s,%s,%s,%s,%s,%s,%s,%f,%s,%s,%s,%f,%f,%f,%f,%f,%f\n",
            Memory->ID, Memory->birth_time, Memory->MID,Memory->cell_state, Memory->MyAffinity,
            Memory->retained_Ag, Memory->total_number_of_divisions,
            Memory->myBCR.nMutFromGermline, Memory->delta_Affinity, Memory->nFDCcontacts,
            Memory->fdc_interaction_time_history,
            Memory->Tc_interaction_history.first,
            Memory->Tc_interaction_history.second, Memory->Selected_by_FDC,Memory->Selected_by_TC,Memory->germline_name.c_str(), Memory->sequence.c_str(), Memory->protein_sequence.c_str(), "Memory", mapToString(Memory->MUTATIONS).c_str(), mapToString(Memory->SILENT_MUTATIONS).c_str(), s.c_str(), st.c_str(),currentTime, mapToString2(Memory->AA_REGIONS).c_str(),mapToString2(Memory->REGIONS).c_str(),  blks.c_str(),Memory->BLIMP1,Memory->BCL6,Memory->IRF4,Memory->BLIMP1_0,Memory->BCL6_0,Memory->IRF4_0);  ///RRR added sequence in fprintf (2 changes)
  }





  fclose(Memory_cells_data);  //#Recheck take care of bins in gle file

  fstream database;
  database.open(output_path + "/Sequence_db.csv",  fstream::out);
  fstream myStream;
  for(auto& kv : Seq_affinity) {
    database << kv.first  <<  ',' << kv.second[0]  <<  ','  << kv.second[1]  <<  ','  << kv.second[2]  <<  ','  << kv.second[3]  <<  ',' << "\n";
  }
  database.close();


  fstream database2;
  database2.open(output_path + "/Restricted_db.csv",  fstream::out);
  fstream myStream2;
  for(auto& kv : Restricted_db) {
    for(auto& kv2 : kv.second) {
      for(auto& kv3 : kv2.second) {
          for(auto& kv4 : kv3.second) {
            database2 << kv.first <<  ',' << kv2.first <<  ',' << kv3.first <<  ',' << kv4.first <<  ','  << kv4.second <<  ',' << "\n";
            ////if (kv4.second.compare("")==0) {
               ////cout << "FALLOnoenEscribir"<<endl;
            ////}
          }      }
    }

  }
  database2.close();

  fstream databaseCDR;
  databaseCDR.open("bcinflow09/CDRs_db.csv",  fstream::out);
  fstream myStreamCDR;
  for(auto& kv : CDRs_affinity) {
    databaseCDR << kv.first[0]  <<  ',' << kv.first[1]  <<  ',' << kv.first[2]  <<  ',' << kv.second[0]  <<  ','  << kv.second[1]  <<  ','  << kv.second[2]  <<  ','  <<            kv.second[3]  <<  ',' << "\n";
  }
  databaseCDR.close();

  fstream databaseFWR;
  databaseFWR.open("bcinflow09/FWRs_db.csv",  fstream::out);
  fstream myStreamFWR;
  for(auto& kv : FWRs_db) {
    databaseFWR << kv.first[0]  <<  ',' << kv.first[1]  <<  ','  << kv.second << "\n";
  }
  databaseFWR.close();

}

 ///Elena
