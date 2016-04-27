//Module Includes
#include "SpectralLibrary.h"
#include "utils.h"
#include "projectionutils.h"

// Header Include
#include "ExecSpectralLibrarySLGFTraining.h"
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

  ExecSpectralLibrarySLGFTraining::ExecSpectralLibrarySLGFTraining(void){
    m_name = "ExecSpectralLibrarySLGFTraining";
    m_type = "ExecSpectralLibrarySLGFTraining";
  }


  ExecSpectralLibrarySLGFTraining::ExecSpectralLibrarySLGFTraining(const ParameterList & inputParams){
    m_name = "ExecSpectralLibrarySLGFTraining";
    m_type = "ExecSpectralLibrarySLGFTraining";
    inputParams.print(cout);
    
    //Reading out p_value_filter
    m_pvalue_filter = inputParams.getValueFloat("pvalue_cutoff", -1);
    m_fdr_filter = inputParams.getValueFloat("fdr_cutoff", -1);
    
    //Characters we are excluding
    aminoacidexclusions = inputParams.getValue("aminoacidexclusion");
    
    //File for amino acid masses
    m_aminoacids_masses_file = inputParams.getValue("aminoacidmassfile");
    
    //Grabbing the location of the model file
    m_model_file_name = inputParams.getValue("model_file");    
    
    //stringSplit(inputParams.getValue("trainannotatedmgf"), m_training_mgf_file_names);
    stringSplit(inputParams.getValue("trainrawnames"), m_training_spectrum_file_names);
    stringSplit(inputParams.getValue("trainpsmnames"), m_training_psm_file_names);
    
    m_output_file_path = inputParams.getValue("histogram_output_prefix");
    
    m_output_mgf_file_name = inputParams.getValue("output_mgf", "");
    
    m_refined_training_spectral_library = inputParams.getValue("SPECTRAL_LIBRARY", "");
  }


  ExecSpectralLibrarySLGFTraining::~ExecSpectralLibrarySLGFTraining(void){
  }


  ExecBase * ExecSpectralLibrarySLGFTraining::clone(const ParameterList & inputParams) const{
    return new ExecSpectralLibrarySLGFTraining(inputParams);
  }

  bool ExecSpectralLibrarySLGFTraining::invoke(void){
    

    vector<int> charge_filter;
    charge_filter.push_back(1);
    charge_filter.push_back(2);
    charge_filter.push_back(3);
    charge_filter.push_back(4);
    charge_filter.push_back(5);
    
    
    DEBUG_MSG("FDR\t"<<m_fdr_filter);
    
    m_training_library.createlibrary(-1, 
                        m_pvalue_filter, 
                        m_fdr_filter,
                        model, 
                        ionsToExtract, 
                        allIons, 
                        aminoacidexclusions, 
                        charge_filter,
                        false);
    
    if(m_refined_training_spectral_library.length() == 0)
        m_distilled_spectral_lib = m_training_library;
    
    m_distilled_spectral_lib.createlibrary(-1, 
                        -1, 
                        -1,
                        model, 
                        ionsToExtract, 
                        allIons, 
                        aminoacidexclusions, 
                        charge_filter,
                        true);
                        
    if(m_output_mgf_file_name.length() > 0){
        m_distilled_spectral_lib.SaveSpecSet_mgf(m_output_mgf_file_name.c_str());
    }
    
    //Setting per peak tolerance for MS2
    m_distilled_spectral_lib.setPeakTolerance(0.5, false);
    m_training_library.setPeakTolerance(0.5, false);
    
    DEBUG_MSG("DONE CREATING TRAINING LIBRARY");
    
    m_distilled_spectral_lib.train_distribution(m_training_library, model, ionsToExtract, allIons);
    //m_distilled_spectral_lib.save_distributions(m_output_file_path);
    
    return true;
  }

  bool ExecSpectralLibrarySLGFTraining::loadInputData(void){      
    
    
    
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
    
    //Old Set
    /*ionsToExtract.push_back("b");
    ionsToExtract.push_back("b-iso");
    ionsToExtract.push_back("b++");
    ionsToExtract.push_back("y");
    ionsToExtract.push_back("y-iso");
    ionsToExtract.push_back("y++");*/
    
    ionsToExtract.push_back("a");
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
    ionsToExtract.push_back("y++-H2O");
    
    
    //Loading the PSMs
    PeptideSpectrumMatchSetSpectralLibraryLoader m_training_psms;
    
    for(int file_idx = 0; file_idx < m_training_psm_file_names.size(); file_idx++){
        string annotation_file_name = m_training_psm_file_names[file_idx];
        PeptideSpectrumMatchSetSpectralLibraryLoader psm_set;
        psm_set.loadMSGFDBResultsFile(annotation_file_name.c_str());
        
        for(int i = 0; i < psm_set.size(); i++){
            m_training_psms.push_back(psm_set[i]);
        }
    }
    
    //Grouping PSMs per file
    map<string, PeptideSpectrumMatchSet> file_psm_map;
    for(int psm_idx = 0; psm_idx < m_training_psms.size(); psm_idx++){
        file_psm_map[m_training_psms[psm_idx]->m_spectrumFile].push_back(m_training_psms[psm_idx]);
    }
    
    int real_total_spectra = 0;
    //Looping through mzxml/annotation files to load
    for(int file_idx = 0; file_idx < m_training_spectrum_file_names.size(); file_idx++){
        //Determining file type by file extension, only supporting plkbin and mzxml
        string spectra_file_name = m_training_spectrum_file_names[file_idx];
        int dotpos = spectra_file_name.find_last_of('.');
        int slashpos = spectra_file_name.find_last_of('/');
        string extension = spectra_file_name.substr(dotpos+1);
        string filename_trunc = spectra_file_name;
        if(slashpos != string::npos)
            filename_trunc = spectra_file_name.substr(slashpos+1);
        
        if(strcmp(extension.c_str(), "mzXML") == 0 || strcmp(extension.c_str(), "mzxml") == 0){
            SpecSet tempspecs;
            cout<<"LOADING\t"<<spectra_file_name<<"\t as an mzxml file...\t";
            //LoadMzxml(spectra_file_name.c_str(), tempspecs, NULL, 0);
            
            for(int i = 0; i < tempspecs.size(); i++){
                tempspecs[i].fileName = filename_trunc;
            }
            
            real_total_spectra += tempspecs.size();
            
            //PeptideSpectrumMatchSetSpectralLibraryLoader psm_set;
            //psm_set.loadSpecnetsResultsFile(annotation_file_name.c_str());
            
            //int original_specs_size = m_training_library.size();
            //m_training_library.appendSpecSet(tempspecs);
            
            
            for(int psm_idx = 0; psm_idx < file_psm_map[filename_trunc].size(); psm_idx++){
                int scannum = file_psm_map[filename_trunc][psm_idx]->m_scanNum;
                //DEBUG_MSG("BEFORE\t"<<m_training_library.size());
                m_training_library.insert(m_training_library.end(), tempspecs.begin() + scannum - 1, tempspecs.begin() + scannum);
                
                file_psm_map[filename_trunc][psm_idx]->m_spectrum = &m_training_library[m_training_library.size()-1];
                m_training_library[m_training_library.size()-1].psmList.push_back(file_psm_map[filename_trunc][psm_idx]);
                //DEBUG_MSG("AFTER\t"<<m_training_library.size());
                
            }
            
            //for(int psm_idx = 0; psm_idx < psm_set.size(); psm_idx++){
            //    psm_set[psm_idx]->m_spectrum = &(m_training_library[psm_set[psm_idx]->m_scanNum - 1 + original_specs_size]);
            //    psm_set[psm_idx]->m_spectrum->psmList.push_back(psm_set[psm_idx]);
            //}
            
            
            cout<<"DONE with "<<m_training_library.size()<< " out of a total of "<<real_total_spectra<<" spectra"<<endl;
            
            
        }
        else if(strcmp(extension.c_str(), "pklbin") == 0){
            //cout<<"LOADED\t"<<spectra_file_name<<"\t as a pklbin file"<<endl;
            //m_training_library.LoadSpecSet_pklbin_with_annotation(spectra_file_name.c_str(), annotation_file_name.c_str());
            
        }
        else if(strcmp(extension.c_str(), "mgf") == 0){
            //cout<<"LOADED\t"<<spectra_file_name<<"\t as a mgf file"<<endl;
            //m_training_library.LoadSpecSet_mgf_with_annotation(spectra_file_name.c_str(), annotation_file_name.c_str());
        }
        else{
            cerr<<"CANNOTLOAD\t"<<spectra_file_name<<endl;
        }
    }
    
    //Loading Spectral Library
    if(m_refined_training_spectral_library.length() > 0){
        DEBUG_MSG("Loading Spectral Library "<<m_refined_training_spectral_library);
        m_distilled_spectral_lib.LoadSpecSet_mgf(m_refined_training_spectral_library.c_str());
    }
    
    /*
    map<string, int> spectrum_start_idx;
    for(int file_idx = 0; file_idx < m_training_psm_file_names.size(); file_idx++){
        string annotation_file_name = m_training_psm_file_names[file_idx];
        PeptideSpectrumMatchSetSpectralLibraryLoader psm_set;
        psm_set.loadMSGFDBResultsFile(annotation_file_name.c_str());
        
        for(int psm_idx = 0; psm_idx < psm_set.size(); psm_idx++){
            //DEBUG_MSG(psm_set[psm_idx]->m_spectrumFile);
            int base_idx = -1;
            if(spectrum_start_idx.find(psm_set[psm_idx]->m_spectrumFile) != spectrum_start_idx.end()){
                //We know the actual start
                base_idx = spectrum_start_idx[psm_set[psm_idx]->m_spectrumFile];
            }
            else{
                for(int spec_idx = 0; spec_idx < m_training_library.size(); spec_idx++){
                    if(m_training_library[spec_idx].fileName == psm_set[psm_idx]->m_spectrumFile){
                        base_idx = spec_idx;
                        spectrum_start_idx[psm_set[psm_idx]->m_spectrumFile] = base_idx;
                        DEBUG_MSG("FOUND\t"<<psm_set[psm_idx]->m_spectrumFile);
                        break;
                    }
                }
            }
            
            if(base_idx == -1) continue;
            
            psm_set[psm_idx]->m_spectrum = &(m_training_library[psm_set[psm_idx]->m_scanNum - 1 + base_idx]);
            psm_set[psm_idx]->m_spectrum->psmList.push_back(psm_set[psm_idx]);
        }
    }*/
    
    return true;
  }


  bool ExecSpectralLibrarySLGFTraining::saveOutputData(void){
        return true;
  }


  bool ExecSpectralLibrarySLGFTraining::saveInputData(std::vector<std::string> & filenames){
    return true;
  }

  bool ExecSpectralLibrarySLGFTraining::loadOutputData(void){
    return true;
  }

  std::vector<ExecBase *> const & ExecSpectralLibrarySLGFTraining::split(int numSplit){
    return m_subModules;
  }

  bool ExecSpectralLibrarySLGFTraining::merge(void){
    return true;
  }


  bool ExecSpectralLibrarySLGFTraining::validateParams(std::string & error){
    return true;
  }
  
  void ExecSpectralLibrarySLGFTraining::load_aminoacid_masses(){
      if(m_aminoacids_masses_file.length() != 0){
          AAJumps jumps(1);
          jumps.loadJumps(m_aminoacids_masses_file.c_str(), true);  //load globally
      }
  }
  
}
