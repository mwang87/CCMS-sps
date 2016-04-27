//Module Includes
#include "SpectralLibrary.h"
#include "utils.h"
#include "projectionutils.h"

// Header Include
#include "ExecSpectralLibrarySLGFCreation.h"
#include "ExecFdrPeptide.h"

// System Include
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

using namespace specnets;
using namespace std;


namespace specnets
{

  ExecSpectralLibrarySLGFCreation::ExecSpectralLibrarySLGFCreation(void){
    m_name = "ExecSpectralLibrarySLGFCreation";
    m_type = "ExecSpectralLibrarySLGFCreation";
  }


  ExecSpectralLibrarySLGFCreation::ExecSpectralLibrarySLGFCreation(const ParameterList & inputParams){
    m_name = "ExecSpectralLibrarySLGFCreation";
    m_type = "ExecSpectralLibrarySLGFCreation";
    inputParams.print(cout);
    
    //Grabbing the location of the model file
    m_model_file_name = inputParams.getValue("model_file");    
    
    //Annotated MGF files most likely comes from some other spectral library
    //And the annotations are embedded in the file, so don't bother with other
    //annotation file
    string input_annotated_mgf_file = inputParams.getValue("annotatedmgf");
    
    //Characters we are excluding
    aminoacidexclusions = inputParams.getValue("aminoacidexclusion");
    
    //File for amino acid masses
    m_aminoacids_masses_file = inputParams.getValue("aminoacidmassfile");
    
    //Specfiying amino acid masses file
    m_aminoacids_masses_file = inputParams.getValue("aminoacidmassfile", "");
    
    //Output file for decoy library
    m_output_decoy_mgf_file_name = inputParams.getValue("outputdecoylibrary", "");
    
    m_output_SLGF_lib = inputParams.getValue("outputSLGFlibrary", "");
    
    //Splitting inputs
    stringSplit(input_annotated_mgf_file, m_mgf_file_names);
    
    m_input_histogram_prefix = inputParams.getValue("histogram_input_prefix", "");
    
    
    //Starting and Ending indicies
    m_SLGF_creation_start_idx = inputParams.getValueInt("SLGF_START_INDEX", -1);
    m_SLGF_creation_end_idx = inputParams.getValueInt("SLGF_END_INDEX", -1);
   
  }


  ExecSpectralLibrarySLGFCreation::~ExecSpectralLibrarySLGFCreation(void){
  }


  ExecBase * ExecSpectralLibrarySLGFCreation::clone(const ParameterList & inputParams) const{
    return new ExecSpectralLibrarySLGFCreation(inputParams);
  }

  bool ExecSpectralLibrarySLGFCreation::invoke(void){
    vector<int> charge_filter;
    
    if(m_SLGF_creation_start_idx != -1 && m_SLGF_creation_end_idx != -1){
        //Filtering library indicies
        SpectralLibraryGFSearch m_temp;
        for(int i = m_SLGF_creation_start_idx; i < m_SLGF_creation_end_idx; i++){
            m_temp.push_back(m_library[i]);
        }
        m_library.resize(0);
        for(int i = 0; i < m_temp.size(); i++){
            m_library.push_back(m_temp[i]);
        }
    }
        
    
    
    m_library.createlibrary(-1, 
                          -1, 
                          model, 
                          ionsToExtract, 
                          allIons, 
                          aminoacidexclusions, 
                          charge_filter,
                          false);
    
    
                          
    DEBUG_MSG("Done loading Library");     

    m_library.load_distributions(m_input_histogram_prefix);
    m_library.generate_distributions(model, ionsToExtract, allIons);
    
    if(m_output_SLGF_lib.length() != 0){
        m_library.SaveSpecSet_mgf(m_output_SLGF_lib.c_str());
    }
    
    return true;
  }

  bool ExecSpectralLibrarySLGFCreation::loadInputData(void){      
    
    load_aminoacid_masses();
      
    //Loading models and setting up ions
    model.LoadModel(m_model_file_name.c_str());
    allIons = "all";
    
    /*
    ionsToExtract.push_back("a");
    ionsToExtract.push_back("b");
    ionsToExtract.push_back("b-iso");
    ionsToExtract.push_back("b-NH3");
    ionsToExtract.push_back("b-H2O");
    ionsToExtract.push_back("b-H2O-NH3");
    ionsToExtract.push_back("b-H2O-H2O");
    ionsToExtract.push_back("b++");
    ionsToExtract.push_back("b++-H2O");
    ionsToExtract.push_back("b++-H2O-H2O");
    ionsToExtract.push_back("b++-NH3");
    ionsToExtract.push_back("b++-iso");
    ionsToExtract.push_back("y");
    ionsToExtract.push_back("y-iso");
    ionsToExtract.push_back("y-NH3");
    ionsToExtract.push_back("y-H2O");
    ionsToExtract.push_back("y-H2O-NH3");
    ionsToExtract.push_back("y-H2O-H2O");
    ionsToExtract.push_back("y++");
    ionsToExtract.push_back("y++-H2O");
    ionsToExtract.push_back("y++-H2O-H2O");
    ionsToExtract.push_back("y++-NH3");
    ionsToExtract.push_back("y++-iso");
    ionsToExtract.push_back("P++");
    ionsToExtract.push_back("P++-H2O");
    ionsToExtract.push_back("P++-NH3");
    ionsToExtract.push_back("P++-H2O-H2O"); */
    
    //High Abundance Ions
    /*ionsToExtract.push_back("a");
    ionsToExtract.push_back("b");
    ionsToExtract.push_back("b-iso");
    ionsToExtract.push_back("b-NH3");
    ionsToExtract.push_back("b-H2O");
    ionsToExtract.push_back("b++");
    ionsToExtract.push_back("b++-NH3");
    ionsToExtract.push_back("b++-H2O");
    ionsToExtract.push_back("y");
    ionsToExtract.push_back("y-iso");
    ionsToExtract.push_back("y-NH3");
    ionsToExtract.push_back("y-H2O");
    ionsToExtract.push_back("y++");
    ionsToExtract.push_back("y++-NH3");
    ionsToExtract.push_back("y++-H2O");*/
    
    //Low Abundance Ions
    ionsToExtract.push_back("a");
    ionsToExtract.push_back("b");
    ionsToExtract.push_back("b++");
    ionsToExtract.push_back("y");
    ionsToExtract.push_back("y++");
    
    //Loading annotated MGF
    for(int file_idx = 0; file_idx < m_mgf_file_names.size(); file_idx++){
        cout<<"Loading MGF\n"<<endl;
        string mgf_file_name = m_mgf_file_names[file_idx];
        m_library.LoadSpecSet_additionalmgf(mgf_file_name.c_str());
        cout<<"LOADED\t"<<mgf_file_name<<"\t as an annotated mgf file"<<endl;
    }
    
    return true;
  }


  bool ExecSpectralLibrarySLGFCreation::saveOutputData(void){
        return true;
  }


  bool ExecSpectralLibrarySLGFCreation::saveInputData(std::vector<std::string> & filenames){
    return true;
  }

  bool ExecSpectralLibrarySLGFCreation::loadOutputData(void){
    return true;
  }

  std::vector<ExecBase *> const & ExecSpectralLibrarySLGFCreation::split(int numSplit){
    return m_subModules;
  }

  bool ExecSpectralLibrarySLGFCreation::merge(void){
    return true;
  }


  bool ExecSpectralLibrarySLGFCreation::validateParams(std::string & error){
    return true;
  }
  
  void ExecSpectralLibrarySLGFCreation::load_aminoacid_masses(){
      if(m_aminoacids_masses_file.length() != 0){
          AAJumps jumps(1);
          jumps.loadJumps(m_aminoacids_masses_file.c_str(), true);  //load globally
      }
  }
  
}
