#ifndef SAMPLE_H
#define SAMPLE_H

#include <exception>
#include <iostream>

#include "TColor.h"
#include "TString.h"


// Encapsulates sample information. Each sample is
// accessed by its unique id as defined in the table
// on the TWiki.
//
// Author: Matthias Schroeder
//         matthias.schroeder@AT@desy.de
//         November 2013
class Sample {
public:
  // Compact sample name without spaces, i.e. suitable
  // for file names etc.
  static TString toTString(unsigned int id);

  // Ntuple file name
  static std::vector<TString> fileNameFullSample(unsigned int id, std::vector<double> &xSecVec);

  // Ntuple file name of a small subset of the full sample
  // for quick tests
  static TString fileNameSubSample(unsigned int id);

  // Full-blown sample name including TLatex commands for
  // plot legends etc.
  static TString label(unsigned int id);

  // A color associated with each sample, e.g. for plotting
  static int color(unsigned int id);
  
//  static vector<double> getxSec(){ return xSecVec; }

private:
  static TString path_;

  static void checkId(unsigned int id);

//  static vector<double> xSecVec;
};


TString Sample::path_ = "root://cmseos:1094//eos/uscms/store/user/lpcsusyhad/PHYS14_720_Dec23_2014/";
//TString Sample::path_ = "root://cmsxrootd-site.fnal.gov//store/user/lpcsusyhad/PHYS14_720_Dec23_2014/";

std::vector<TString> Sample::fileNameFullSample(unsigned int id, std::vector<double> &xSecVec) {
  xSecVec.clear();
  checkId(id);

  std::vector<TString> nameVec;
  TString name("");
//  if(      id ==  1  ) name += "HTMHT-Run2012.root";
//  else if( id ==  2  ) name += "MuHad.root";

  if( id == 11  ){
     name = path_+"benwu/PU20bx25_WJetsToLNu_HT-400to600_madgraph-tauola/WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola/PHYS14_PU20bx25_PHYS14_25_V1-FLAT/141224_220815/0000/*.root"; nameVec.push_back(name); xSecVec.push_back(55.61);
     name = path_+"benwu/PU20bx25_WJetsToLNu_HT-600toInf_madgraph-tauola/WJetsToLNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola/PHYS14_PU20bx25_PHYS14_25_V1-FLAT/141229_230015/0000/stopFlatNtuples_*.root"; nameVec.push_back(name); xSecVec.push_back(18.81);
  } else if( id == 12  ){
     name = path_+"lhx/PU20bx25_TTJets_MSDecaysCKM_madgraph-tauola/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/PHYS14_PU20bx25_PHYS14_25_V1-FLAT/141224_052628/0000/stopFlatNtuples_*.root"; nameVec.push_back(name); xSecVec.push_back(818.8);
  } else if( id == 13 ){
     name = path_+"pastika/ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola/PHYS14_PU20bx25_PHYS14_25_V1-FLAT/141227_223010/0000/stopFlatNtuples_*.root"; nameVec.push_back(name); xSecVec.push_back(4.113);
  } else if( id == 14 ){
     name = path_+"malik/PU20bx25_QCD_HT_250To500_13TeV-madgraph_V1_ext1-v2/QCD_HT_250To500_13TeV-madgraph/PHYS14_PU20bx25_PHYS14_25_V1-FLAT/141225_194319/0000/stopFlatNtuples_*.root"; nameVec.push_back(name); xSecVec.push_back(670500);
     name = path_+"malik/PU20bx25_QCD_HT_500To1000_13TeV-madgraph_V1_ext1-v1/QCD_HT-500To1000_13TeV-madgraph/PHYS14_PU20bx25_PHYS14_25_V1-FLAT/141225_201220/0000/stopFlatNtuples_*.root"; nameVec.push_back(name); xSecVec.push_back(26740);
     name = path_+"malik/PU20bx25_QCD_HT_1000ToInf_13TeV-madgraph_V1_ext1-v1/QCD_HT_1000ToInf_13TeV-madgraph/PHYS14_PU20bx25_PHYS14_25_V1-FLAT/141225_194928/0000/stopFlatNtuples_*.root"; nameVec.push_back(name); xSecVec.push_back(769.7);
//     name = path_+"malik/PU20bx25_QCD_HT_250To500_13TeV-madgraph*/QCD_HT_250To500_13TeV-madgraph/PHYS14_PU20bx25_PHYS14_25_V1-FLAT/*/0000/stopFlatNtuples_*.root"; nameVec.push_back(name); xSecVec.push_back(670500);
//     name = path_+"malik/PU20bx25_QCD_HT_500To1000_13TeV-madgraph*/QCD_HT_500To1000_13TeV-madgraph/PHYS14_PU20bx25_PHYS14_25_V1-FLAT/*/0000/stopFlatNtuples_*.root"; nameVec.push_back(name); xSecVec.push_back(26740);
//     name = path_+"malik/PU20bx25_QCD_HT_1000ToInf_13TeV-madgraph*/QCD_HT_1000ToInf_13TeV-madgraph/PHYS14_PU20bx25_PHYS14_25_V1-FLAT/*/0000/stopFlatNtuples_*.root"; nameVec.push_back(name); xSecVec.push_back(769.7);
  } else if( id == 21  ){
     name = path_+"lhx/PU20bx25_T1tttt_mGl-1500_mLSP-100-madgraph-tauola/SMS-T1tttt_2J_mGl-1500_mLSP-100_Tune4C_13TeV-madgraph-tauola/PHYS14_PU20bx25_PHYS14_25_V1-FLAT/141225_024211/0000/stopFlatNtuples_1.root"; nameVec.push_back(name); xSecVec.push_back(0.0141903);
  } else if( id == 22  ){
     name = path_+"lhx/PU20bx25_T1tttt_mGl-1200_mLSP-800-madgraph-tauola/SMS-T1tttt_2J_mGl-1200_mLSP-800_Tune4C_13TeV-madgraph-tauola/PHYS14_PU20bx25_PHYS14_25_V1-FLAT/141225_023727/0000/stopFlatNtuples_1.root"; nameVec.push_back(name); xSecVec.push_back(0.0856418);
  }

  return nameVec;
}


TString Sample::fileNameSubSample(unsigned int id) {
  checkId(id);

  TString name("");
  if(      id ==  1  ) name += "HTMHT-Run2012.root";
  else if( id ==  2  ) name += "MuHad.root";
  else if( id == 11  ) name += "Summer12-WJetsHT400toInf_1.root";
  else if( id == 12  ) name += "Summer12-TTJets-SemiLep_1.root";
  else if( id == 13  ) name += "Summer12-ZJetsHT400toInf.root";
  else if( id == 14  ) name += "Summer12-QCDHT250toInf.root";
  else if( id == 21  ) name += "Summer12-SUSY_LM6.root";
  else if( id == 22  ) name += "Summer12-SUSY_LM9.root";

  return path_+name;
}


// Return the label for a given sample
TString Sample::label(unsigned int id) {
  checkId(id);

  TString label("");
  if(      id ==  1  ) label += "Data";
  else if( id ==  2  ) label += "Data";
  else if( id == 11  ) label += "W(l#nu)+Jets";
  else if( id == 12  ) label += "t#bar{t}+Jets";
  else if( id == 13  ) label += "Z(#nu#bar{#nu})+Jets";
  else if( id == 14  ) label += "QCD";
  else if( id == 21  ) label += "LM6";
  else if( id == 22  ) label += "LM9";

  return label;
}


TString Sample::toTString(unsigned int id) {
  checkId(id);

  TString str("");
  if(      id ==  1  ) str += "Data";
  else if( id ==  2  ) str += "Data";
  else if( id == 11  ) str += "WJets";
  else if( id == 12  ) str += "TTJets";
  else if( id == 13  ) str += "ZJets";
  else if( id == 14  ) str += "QCD";
  else if( id == 21  ) str += "LM6";
  else if( id == 22  ) str += "LM9";

  return str;
}


int Sample::color(unsigned int id) {
  checkId(id);

  int color = kBlack;
  if(      id ==  1  ) color = kBlack;
  else if( id ==  1  ) color = kBlack;
  else if( id == 11  ) color = kGreen+1;
  else if( id == 12  ) color = kRed;
  else if( id == 13  ) color = kYellow+2;
  else if( id == 14  ) color = kRed+3;
  else if( id == 21  ) color = kBlue;
  else if( id == 22  ) color = kBlue+3;

  return color;
}
  

void Sample::checkId(unsigned int id) {
  if( id != 1 && id != 2 && !(id >= 11 && id <= 16) && !(id >= 21 && id <= 22) ) {
    std::cerr << "\n\nERROR invalid sample id " << id << std::endl;
    throw std::exception();
  } 
}

#endif
