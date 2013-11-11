//---------------------------
// GEN Diffractive Z Analyzer
//---------------------------

#define _USE_MATH_DEFINES
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "TH1F.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TLorentzVector.h"
#include <math.h>

class SDDYAnalyzer: public edm::EDAnalyzer {
  public:
    /// Constructor
    SDDYAnalyzer(const edm::ParameterSet& pset);

    /// Destructor
    virtual ~SDDYAnalyzer();

    void analyze(const edm::Event & event, const edm::EventSetup& eventSetup);
    void fillHistos(int index);

    virtual void beginJob(); 
    virtual void endJob() ;

  private:
    // Input from cfg file
    edm::InputTag genParticlesTag_;
    double Ebeam_;
    int particle1Id_;
    int particle2Id_;
    int index;
    bool debug;
    bool cmsAccept;

    // Histograms
    std::vector<TH1F*> hVectorPartEta;
    std::vector<TH1F*> hVectorPartPt;
    std::vector<TH1F*> hVectorPartPhi;
    std::vector<TH1F*> hVectorPartVx;
    std::vector<TH1F*> hVectorPartVy;
    std::vector<TH1F*> hVectorPartVz;
    std::vector<TH1F*> hVectorPartPx;
    std::vector<TH1F*> hVectorPartPy;
    std::vector<TH1F*> hVectorPartPz;
    std::vector<TH1F*> hVectorPartPf;
    std::vector<TH1F*> hVectorPDGId;
    std::vector<TH1F*> hVectorPartEnergy;
    std::vector<TH1F*> hVectorPartEnergyVsEta;
    std::vector<TH1F*> hVectorPartPDGId;

    std::vector<TH1F*> hVectorHFPlusE;
    std::vector<TH1F*> hVectorHFMinusE;
    std::vector<TH1F*> hVectorCastorE;

    std::vector<TH1F*> hVectorLepton1Eta;
    std::vector<TH1F*> hVectorLepton1Phi;
    std::vector<TH1F*> hVectorLepton1Pt;
    std::vector<TH1F*> hVectorLepton1Energy;
    std::vector<TH1F*> hVectorLepton1px;
    std::vector<TH1F*> hVectorLepton1py;
    std::vector<TH1F*> hVectorLepton1pz;
    std::vector<TH1F*> hVectorLepton1pf;
    std::vector<TH1F*> hVectorLepton1vx;
    std::vector<TH1F*> hVectorLepton1vy;
    std::vector<TH1F*> hVectorLepton1vz;

    std::vector<TH1F*> hVectorLepton2Eta;
    std::vector<TH1F*> hVectorLepton2Phi;
    std::vector<TH1F*> hVectorLepton2Pt;
    std::vector<TH1F*> hVectorLepton2Energy;
    std::vector<TH1F*> hVectorLepton2px;
    std::vector<TH1F*> hVectorLepton2py;
    std::vector<TH1F*> hVectorLepton2pz;
    std::vector<TH1F*> hVectorLepton2pf;
    std::vector<TH1F*> hVectorLepton2vx;
    std::vector<TH1F*> hVectorLepton2vy;
    std::vector<TH1F*> hVectorLepton2vz;

    std::vector<TH1F*> hVectorDileptonEta;
    std::vector<TH1F*> hVectorDileptonPhi;
    std::vector<TH1F*> hVectorDileptonPt;
    std::vector<TH1F*> hVectorDileptonM;

    std::vector<TH1F*> hVectorLeptonDeltaEta;
    std::vector<TH1F*> hVectorLeptonDeltaPhi;
    std::vector<TH1F*> hVectorLeptonDeltaPt;
    std::vector<TH1F*> hVectorLeptonVertexD;

    std::vector<TH1F*> hVectortprotonplus;
    std::vector<TH1F*> hVectortprotonminus;

    std::vector<TH1F*> hVectorXiZplus;
    std::vector<TH1F*> hVectorXiZminus;
    std::vector<TH1F*> hVectorXidiffplus;
    std::vector<TH1F*> hVectorXidiffminus;
    std::vector<TH1F*> hVectorXiprotonplus;
    std::vector<TH1F*> hVectorXiprotonminus;
    std::vector<TH1F*> hVectorXiDiffproton;
    std::vector<TH1F*> hVectorXiAll;
    std::vector<TH1F*> hVectorXiAllplus;
    std::vector<TH1F*> hVectorXiAllminus;
    std::vector<TH1F*> hVectorXiDiffAll;

    std::vector<TH1F*> hVectorGapSizeAll;
    std::vector<TH1F*> hVectorGapSizeZ;


    TH1F* hPartEta;
    TH1F* hPartPt;
    TH1F* hPartPhi;
    TH1F* hPartVx;
    TH1F* hPartVy;
    TH1F* hPartVz;
    TH1F* hPartPx;
    TH1F* hPartPy;
    TH1F* hPartPz;
    TH1F* hPartPf;
    TH1F* hPartPDGId;
    TH1F* hPartEnergy;
    TH1F* hPartEnergyVsEta;

    TH1F* hHFPlusE;
    TH1F* hHFMinusE;
    TH1F* hCastorE;

    TH1F* hLepton1Eta;
    TH1F* hLepton1Phi;
    TH1F* hLepton1Pt;
    TH1F* hLepton1Energy;
    TH1F* hLepton1px;
    TH1F* hLepton1py;
    TH1F* hLepton1pz;
    TH1F* hLepton1pf;
    TH1F* hLepton1vx;
    TH1F* hLepton1vy;
    TH1F* hLepton1vz;

    TH1F* hLepton2Eta;
    TH1F* hLepton2Phi;
    TH1F* hLepton2Pt;
    TH1F* hLepton2Energy;
    TH1F* hLepton2px;
    TH1F* hLepton2py;
    TH1F* hLepton2pz;
    TH1F* hLepton2pf;
    TH1F* hLepton2vx;
    TH1F* hLepton2vy;
    TH1F* hLepton2vz;

    TH1F* hDileptonEta;
    TH1F* hDileptonPhi;
    TH1F* hDileptonPt;
    TH1F* hDileptonM;

    TH1F* hLeptonDeltaEta;
    TH1F* hLeptonDeltaPhi;
    TH1F* hLeptonDeltaPt;
    TH1F* hLeptonVertexD;

    TH1F* htprotonplus;
    TH1F* htprotonminus;

    TH1F* hXiZplus;
    TH1F* hXiZminus;
    TH1F* hXidiffplus;
    TH1F* hXidiffminus;
    TH1F* hXiprotonplus;
    TH1F* hXiprotonminus;
    TH1F* hXiDiffproton;
    TH1F* hXiAll;
    TH1F* hXiAllplus;
    TH1F* hXiAllminus;
    TH1F* hXiDiffAll;

    TH1F* hVectorGapSizeAll;
    TH1F* hVectorGapSizeZ;

    // Histogram variables
    double sumHFMinusGEN, sumCastorGEN, sumHFPlusGEN;
    double deltaeta, deltaphi, deltapt, vertex_d;
    double genEPlusPz, genEMinusPz;
    double xi, xi_minus, xi_plus, xiZ_plus, xiZ_minus, xiZ, xi_diff;
    double xi_diff_plus, xi_diff_minus;
    double l1eta, l1phi, l1pt, l1energy;
    double l2eta, l2phi, l2pt, l2energy;
    double l1px, l1py, l1pz, l1pf;
    double l2px, l2py, l2pz, l2pf;
    double l1vx,  l1vy, l1vz;
    double l2vx, l2vy, l2vz;
    double xiProtonPlus;
    double xiProtonMinus;
    double xiDiffProton;
    double proton_eta_plus, proton_phi_plus, proton_pt_plus, proton_energy_plus;
    double proton_px_plus, proton_py_plus, proton_pz_plus, proton_pf_plus;
    double proton_eta_minus, proton_phi_minus, proton_pt_minus, proton_energy_minus;
    double proton_px_minus, proton_py_minus, proton_pz_minus, proton_pf_minus;
    double dibosonEta, dibosonPhi, dibosonPt, dibosonM;
    double diboson_px, diboson_py, diboson_pz, diboson_pf, diboson_energy;
    double t_plus, t_minus;
    double pz_cut;
    double pf_gen;
    double gapsizeall, gapsizeZ;

    bool lepton1, lepton2, dilepton, leptonAcceptance;
    bool single_gap, double_gap;
    bool zboson, zboson_diff;
    bool HF_CASTOR_gap;
    bool protonplus;
    bool protonminus;
    bool ptcut;

    std::vector <std::string> Group1;

    std::vector<double> energy_genAll;
    std::vector<double> px_genAll;
    std::vector<double> py_genAll;
    std::vector<double> pz_genAll;
    std::vector<double> etaAll;
    std::vector<double> phiAll;
    std::vector<double> ptAll;
    std::vector<double> vxAll;
    std::vector<double> vyAll;
    std::vector<double> vzAll;
    std::vector<double> pf_genAll;
    std::vector<double> pt_genAll;
    std::vector<double> pdgIdAll;

};

/// Constructor
SDDYAnalyzer::SDDYAnalyzer(const edm::ParameterSet& pset)
{
  genParticlesTag_ = pset.getParameter<edm::InputTag>("GenParticleTag");
  particle1Id_ = pset.getParameter<int>("Particle1Id");
  particle2Id_ = pset.getParameter<int>("Particle2Id");
  Ebeam_ = pset.getParameter<double>("EBeam");
  debug = pset.getUntrackedParameter<bool>("debug",false);
  cmsAccept = pset.getUntrackedParameter<bool>("cmsAccept",false);
  if(debug){
    std::cout << ">>> First particle Id: " << particle1Id_ << std::endl;
    std::cout << ">>> Second particle Id: " << particle2Id_ << std::endl;
  }
}

/// Destructor
SDDYAnalyzer::~SDDYAnalyzer(){
}

void SDDYAnalyzer::beginJob(){

  edm::Service<TFileService> fs;
  TH1::SetDefaultSumw2(true);

  std::string step0 = "CMS_particles";
  std::string step1 = "CMS_dilepton";
  std::string step2 = "CMS_dileptonPt";
  std::string step3 = "CMS_BosonZ";
  std::string step4 = "CMS_BosonZGapCastor";
  std::string step5 = "CMS_BosonZGapCastorAndHF";

  Group1.push_back(step0);
  Group1.push_back(step1);
  Group1.push_back(step2);
  Group1.push_back(step3);
  Group1.push_back(step4);
  Group1.push_back(step5);

  for (std::vector<std::string>::size_type j=0; j<Group1.size(); j++){

    char hPartEtaN[300];
    sprintf(hPartEtaN,"hPartEta_%s",Group1.at(j).c_str());
    hPartEta = fs->make<TH1F>(hPartEtaN,";#eta; #N Events",500,-15.,15.);
    hVectorPartEta.push_back(hPartEta);

    char hPartPhiN[300];
    sprintf(hPartPhiN,"hPartPhi_%s",Group1.at(j).c_str());
    hPartPhi = fs->make<TH1F>(hPartPhiN,";#phi; #N Events",100,-1.2*M_PI,1.2*M_PI);
    hVectorPartPhi.push_back(hPartPhi);

    char hPartPtN[300];
    sprintf(hPartPtN,"hPartPt_%s",Group1.at(j).c_str());
    hPartPt = fs->make<TH1F>(hPartPtN,";pT [GeV]; #N Events",350,0.,350.);
    hVectorPartPt.push_back(hPartPt);

    char hPartVxN[300];
    sprintf(hPartVxN,"hPartVx_%s",Group1.at(j).c_str());
    hPartVx = fs->make<TH1F>(hPartVxN,";Vertex X position [mm]; #N Events",100,0.,20.);
    hVectorPartVx.push_back(hPartVx);

    char hPartVyN[300];
    sprintf(hPartVyN,"hPartVy_%s",Group1.at(j).c_str());
    hPartVy = fs->make<TH1F>(hPartVyN,";Vertex Y position [mm]; #N Events",100,0.,20.);
    hVectorPartVy.push_back(hPartVy);

    char hPartVzN[300];
    sprintf(hPartVzN,"hPartVz_%s",Group1.at(j).c_str());
    hPartVz = fs->make<TH1F>(hPartVzN,";Vertex Z position [mm]; #N Events",100,0.,20.);
    hVectorPartVz.push_back(hPartVz);

    char hPartPxN[300];
    sprintf(hPartPxN,"hPartPx_%s",Group1.at(j).c_str());
    hPartPx = fs->make<TH1F>(hPartPxN,";#p_{x} [GeV]; #N Events",350,0.,350.);
    hVectorPartPx.push_back(hPartPx);

    char hPartPyN[300];
    sprintf(hPartPyN,"hPartPy_%s",Group1.at(j).c_str());
    hPartPy = fs->make<TH1F>(hPartPyN,";#p_{y} [GeV]; #N Events",350,0.,350.);
    hVectorPartPy.push_back(hPartPy);

    char hPartPzN[300];
    sprintf(hPartPzN,"hPartPz_%s",Group1.at(j).c_str());
    hPartPz = fs->make<TH1F>(hPartPzN,";#p_{z} [GeV]; #N Events",350,0.,350.);
    hVectorPartPz.push_back(hPartPz);

    char hPartPfN[300];
    sprintf(hPartPfN,"hPartPf_%s",Group1.at(j).c_str());
    hPartPf = fs->make<TH1F>(hPartPfN,";#p_{f} [GeV]; #N Events",800,0.,800.);
    hVectorPartPf.push_back(hPartPf);

    char hPartEnergyN[300];
    sprintf(hPartEnergyN,"hPartEnergy_%s",Group1.at(j).c_str());
    hPartEnergy = fs->make<TH1F>(hPartEnergyN,";#p_{f} [GeV]; #N Events",800,0.,800.);
    hVectorPartEnergy.push_back(hPartEnergy);

    char hPartEnergyVsEtaN[300];
    sprintf(hPartEnergyVsEtaN,"hPartEnergyVsEta_%s",Group1.at(j).c_str());
    hPartEnergyVsEta = fs->make<TH1F>(hPartEnergyVsEtaN,";#eta; Energy [GeV]",500,-15.,15.);
    hVectorPartEnergyVsEta.push_back(hPartEnergyVsEta);

    char hPartPDGIdN[300];
    sprintf(hPartPDGIdN,"hPartPDGId_%s",Group1.at(j).c_str());
    hPartPDGId = fs->make<TH1F>(hPartPDGIdN,";pdg Id; #N Events",5000,-2500,2500);
    hVectorPartPDGId.push_back(hPartPDGId);

    char hHFPlusEN[300];
    sprintf(hHFPlusEN,"hHFPlusE_%s",Group1.at(j).c_str());
    hHFPlusE = fs->make<TH1F>(hHFPlusEN,"HF Plus GEN Energy; #HF^{+} #sum Energy [GeV]; #N Events",2000,0.,1000.);
    hVectorHFPlusE.push_back(hHFPlusE);

    char hHFMinusEN[300];
    sprintf(hHFMinusEN,"hHFMinusE_%s",Group1.at(j).c_str());
    hHFMinusE = fs->make<TH1F>(hHFMinusEN,"HF Minus GEN Energy; #HF^{-} #sum Energy [GeV]; #N Events",2000,0.,1000.);
    hVectorHFMinusE.push_back(hHFMinusE);

    char hCastorEN[300];
    sprintf(hCastorEN,"hCastorE_%s",Group1.at(j).c_str());
    hCastorE = fs->make<TH1F>(hCastorEN,"CASTOR GEN Energy; #CASTOR^{-} #sum Energy [GeV]; #N Events",2000,0.,1000.);
    hVectorCastorE.push_back(hCastorE);

    char hLepton1EtaN[300];
    sprintf(hLepton1EtaN,"hLepton1Eta_%s",Group1.at(j).c_str());
    hLepton1Eta = fs->make<TH1F>(hLepton1EtaN,";#eta; #N Events",500,-15.,15.);
    hVectorLepton1Eta.push_back(hLepton1Eta);

    char hLepton1PhiN[300];
    sprintf(hLepton1PhiN,"hLepton1Phi_%s",Group1.at(j).c_str());
    hLepton1Phi = fs->make<TH1F>(hLepton1PhiN,";#phi; #N Events",100,-1.2*M_PI,1.2*M_PI);
    hVectorLepton1Phi.push_back(hLepton1Phi);

    char hLepton1PtN[300];
    sprintf(hLepton1PtN,"hLepton1Pt_%s",Group1.at(j).c_str());
    hLepton1Pt = fs->make<TH1F>(hLepton1PtN,";pT [GeV]; #N Events",350,0.,350.);
    hVectorLepton1Pt.push_back(hLepton1Pt);

    char hLepton1EnergyN[300];
    sprintf(hLepton1EnergyN,"hLepton1Energy_%s",Group1.at(j).c_str());
    hLepton1Energy = fs->make<TH1F>(hLepton1EnergyN,";Energy [GeV]; #N Events",1200,0,600);
    hVectorLepton1Energy.push_back(hLepton1Energy);

    char hLepton1pxN[300];
    sprintf(hLepton1pxN,"hLepton1px_%s",Group1.at(j).c_str());
    hLepton1px = fs->make<TH1F>(hLepton1pxN,";#p_{x} [GeV]; #N Events",350,0.,350.);
    hVectorLepton1px.push_back(hLepton1px);

    char hLepton1pyN[300];
    sprintf(hLepton1pyN,"hLepton1py_%s",Group1.at(j).c_str());
    hLepton1py = fs->make<TH1F>(hLepton1pyN,";#p_{y} [GeV]; #N Events",350,0.,350.);
    hVectorLepton1py.push_back(hLepton1py);

    char hLepton1pzN[300];
    sprintf(hLepton1pzN,"hLepton1pz_%s",Group1.at(j).c_str());
    hLepton1pz = fs->make<TH1F>(hLepton1pzN,";#p_{z} [GeV]; #N Events",350,0.,350.);
    hVectorLepton1pz.push_back(hLepton1pz);

    char hLepton1pfN[300];
    sprintf(hLepton1pfN,"hLepton1pf_%s",Group1.at(j).c_str());
    hLepton1pf = fs->make<TH1F>(hLepton1pfN,";#p_{f} [GeV]; #N Events",800,0.,800.);
    hVectorLepton1pf.push_back(hLepton1pf);

    char hLepton1vxN[300];
    sprintf(hLepton1vxN,"hLepton1vx_%s",Group1.at(j).c_str());
    hLepton1vx = fs->make<TH1F>(hLepton1vxN,";Vertex X position [mm]; #N Events",100,0.,20.);
    hVectorLepton1vx.push_back(hLepton1vx);

    char hLepton1vyN[300];
    sprintf(hLepton1vyN,"hLepton1vy_%s",Group1.at(j).c_str());
    hLepton1vy = fs->make<TH1F>(hLepton1vyN,";Vertex Y position [mm]; #N Events",100,0.,20.);
    hVectorLepton1vy.push_back(hLepton1vy);

    char hLepton1vzN[300];
    sprintf(hLepton1vzN,"hLepton1vz_%s",Group1.at(j).c_str());
    hLepton1vz = fs->make<TH1F>(hLepton1vzN,";Vertex Z position [mm]; #N Events",100,0.,20.);
    hVectorLepton1vz.push_back(hLepton1vz);

    char hLepton2EtaN[300];
    sprintf(hLepton2EtaN,"hLepton2Eta_%s",Group1.at(j).c_str());
    hLepton2Eta = fs->make<TH1F>(hLepton2EtaN,";#eta; #N Events",500,-15.,15.);
    hVectorLepton2Eta.push_back(hLepton2Eta);

    char hLepton2PhiN[300];
    sprintf(hLepton2PhiN,"hLepton2Phi_%s",Group1.at(j).c_str());
    hLepton2Phi = fs->make<TH1F>(hLepton2PhiN,";#phi; #N Events",100,-1.2*M_PI,1.2*M_PI);
    hVectorLepton2Phi.push_back(hLepton2Phi);

    char hLepton2PtN[300];
    sprintf(hLepton2PtN,"hLepton2Pt_%s",Group1.at(j).c_str());
    hLepton2Pt = fs->make<TH1F>(hLepton2PtN,";pT [GeV]; #N Events",350,0.,350.);
    hVectorLepton2Pt.push_back(hLepton2Pt);

    char hLepton2EnergyN[300];
    sprintf(hLepton2EnergyN,"hLepton2Energy_%s",Group1.at(j).c_str());
    hLepton2Energy = fs->make<TH1F>(hLepton2EnergyN,";Energy [GeV]; #N Events",1200,0,600);
    hVectorLepton2Energy.push_back(hLepton2Energy);

    char hLepton2pxN[300];
    sprintf(hLepton2pxN,"hLepton2px_%s",Group1.at(j).c_str());
    hLepton2px = fs->make<TH1F>(hLepton2pxN,";#p_{x} [GeV]; #N Events",350,0.,350.);
    hVectorLepton2px.push_back(hLepton2px);

    char hLepton2pyN[300];
    sprintf(hLepton2pyN,"hLepton2py_%s",Group1.at(j).c_str());
    hLepton2py = fs->make<TH1F>(hLepton2pyN,";#p_{y} [GeV]; #N Events",350,0.,350.);
    hVectorLepton2py.push_back(hLepton2py);

    char hLepton2pzN[300];
    sprintf(hLepton2pzN,"hLepton2pz_%s",Group1.at(j).c_str());
    hLepton2pz = fs->make<TH1F>(hLepton2pzN,";#p_{z} [GeV]; #N Events",350,0.,350.);
    hVectorLepton2pz.push_back(hLepton2pz);

    char hLepton2pfN[300];
    sprintf(hLepton2pfN,"hLepton2pf_%s",Group1.at(j).c_str());
    hLepton2pf = fs->make<TH1F>(hLepton2pfN,";#p_{f} [GeV]; #N Events",800,0.,800.);
    hVectorLepton2pf.push_back(hLepton2pf);

    char hLepton2vxN[300];
    sprintf(hLepton2vxN,"hLepton2vx_%s",Group1.at(j).c_str());
    hLepton2vx = fs->make<TH1F>(hLepton2vxN,";Vertex X position [mm]; #N Events",100,0.,20.);
    hVectorLepton2vx.push_back(hLepton2vx);

    char hLepton2vyN[300];
    sprintf(hLepton2vyN,"hLepton2vy_%s",Group1.at(j).c_str());
    hLepton2vy = fs->make<TH1F>(hLepton2vyN,";Vertex Y position [mm]; #N Events",100,0.,20.);
    hVectorLepton2vy.push_back(hLepton2vy);

    char hLepton2vzN[300];
    sprintf(hLepton2vzN,"hLepton2vz_%s",Group1.at(j).c_str());
    hLepton2vz = fs->make<TH1F>(hLepton2vzN,";Vertex Z position [mm]; #N Events",100,0.,20.);
    hVectorLepton2vz.push_back(hLepton2vz);

    char hDileptonEtaN[300];
    sprintf(hDileptonEtaN,"hDileptonEta_%s",Group1.at(j).c_str());
    hDileptonEta = fs->make<TH1F>(hDileptonEtaN,";#eta; #N Events",500,-15.,15.);
    hVectorDileptonEta.push_back(hDileptonEta);

    char hDileptonPhiN[300];
    sprintf(hDileptonPhiN,"hDileptonPhi_%s",Group1.at(j).c_str());
    hDileptonPhi = fs->make<TH1F>(hDileptonPhiN,";#phi; #N Events",100,-1.2*M_PI,1.2*M_PI);
    hVectorDileptonPhi.push_back(hDileptonPhi);

    char hDileptonPtN[300];
    sprintf(hDileptonPtN,"hDileptonPt_%s",Group1.at(j).c_str());
    hDileptonPt = fs->make<TH1F>(hDileptonPtN,";pT [GeV]; #N Events",350,0.,350.);
    hVectorDileptonPt.push_back(hDileptonPt);

    char hDileptonMN[300];
    sprintf(hDileptonMN,"hDileptonM_%s",Group1.at(j).c_str());
    hDileptonM = fs->make<TH1F>(hDileptonMN,";#M_{ll}; #N Events",2000,0,500.);
    hVectorDileptonM.push_back(hDileptonM);

    char hLeptonDeltaEtaN[300];
    sprintf(hLeptonDeltaEtaN,"hLeptonDeltaEta_%s",Group1.at(j).c_str());
    hLeptonDeltaEta = fs->make<TH1F>(hLeptonDeltaEtaN,";#Delta#eta Lepton; #N Events",100,0.,20.);
    hVectorLeptonDeltaEta.push_back(hLeptonDeltaEta);

    char hLeptonDeltaPhiN[300];
    sprintf(hLeptonDeltaPhiN,"hLeptonDeltaPhi_%s",Group1.at(j).c_str());
    hLeptonDeltaPhi = fs->make<TH1F>(hLeptonDeltaPhiN,";#Delta#phi Lepton; #N Events",100,0.,20.);
    hVectorLeptonDeltaPhi.push_back(hLeptonDeltaPhi);

    char hLeptonDeltaPtN[300];
    sprintf(hLeptonDeltaPtN,"hLeptonDeltaPt_%s",Group1.at(j).c_str());
    hLeptonDeltaPt = fs->make<TH1F>(hLeptonDeltaPtN,";#Delta pT Lepton [GeV]; #N Events",500,0.,200.);
    hVectorLeptonDeltaPt.push_back(hLeptonDeltaPt);

    char hLeptonVertexDN[300];
    sprintf(hLeptonVertexDN,"hLeptonVertexD_%s",Group1.at(j).c_str());
    hLeptonVertexD = fs->make<TH1F>(hLeptonVertexDN,";#Delta Vertex Lepton position [mm]; #N Events",100,0.,20.);
    hVectorLeptonVertexD.push_back(hLeptonVertexD);

    char htprotonplusN[300];
    sprintf(htprotonplusN,"htprotonplus_%s",Group1.at(j).c_str());
    htprotonplus = fs->make<TH1F>(htprotonplusN,";|t| [GeV], positive scattering p; #N Events",1000,0.,10.);
    hVectortprotonplus.push_back(htprotonplus);

    char htprotonminusN[300];
    sprintf(htprotonminusN,"htprotonminus_%s",Group1.at(j).c_str());
    htprotonminus = fs->make<TH1F>(htprotonminusN,";|t| [GeV], negative scattering p; #N Events",1000,0.,10.);
    hVectortprotonminus.push_back(htprotonminus);

    char hXiZplusN[300];
    sprintf(hXiZplusN,"hXiZplus_%s",Group1.at(j).c_str());
    hXiZplus = fs->make<TH1F>(hXiZplusN,";#xi^{+}; #N Events",50000,0,50000);
    hVectorXiZplus.push_back(hXiZplus);

    char hXiZminusN[300];
    sprintf(hXiZminusN,"hXiZminus_%s",Group1.at(j).c_str());
    hXiZminus = fs->make<TH1F>(hXiZminusN,";#xi^{-}; #N Events",50000,0,50000);
    hVectorXiZminus.push_back(hXiZminus);

    char hXidiffplusN[300];
    sprintf(hXidiffplusN,"hXidiffplus_%s",Group1.at(j).c_str());
    hXidiffplus = fs->make<TH1F>(hXidiffplusN,";|#xi_{particles}^{+} - #xi_{Z}^{+}|; #N Events",50000,0,50000);
    hVectorXidiffplus.push_back(hXidiffplus);

    char hXidiffminusN[300];
    sprintf(hXidiffminusN,"hXidiffminus_%s",Group1.at(j).c_str());
    hXidiffminus = fs->make<TH1F>(hXidiffminusN,";|#xi_{particles}^{-} - #xi_{Z}^{-}|; #N Events",50000,0,50000);
    hVectorXidiffminus.push_back(hXidiffminus);

    char hXiprotonplusN[300];
    sprintf(hXiprotonplusN,"hXiprotonplus_%s",Group1.at(j).c_str());
    hXiprotonplus = fs->make<TH1F>(hXiprotonplusN,"; #xi^{+}; #N Events",1000,0,1.);
    hVectorXiprotonplus.push_back(hXiprotonplus);

    char hXiprotonminusN[300];
    sprintf(hXiprotonminusN,"hXiprotonminus_%s",Group1.at(j).c_str());
    hXiprotonminus = fs->make<TH1F>(hXiprotonminusN,"; #xi^{-}; #N Events",1000,0,1.);
    hVectorXiprotonminus.push_back(hXiprotonminus);

    char hXiDiffprotonN[300];
    sprintf(hXiDiffprotonN,"hXiDiffproton_%s",Group1.at(j).c_str());
    hXiDiffproton = fs->make<TH1F>(hXiDiffprotonN,"; |#xi^{+} - #xi^{-}|; #N Events",1000,0,1.);
    hVectorXiDiffproton.push_back(hXiDiffproton);

    char hXiAllplusN[300];
    sprintf(hXiAllplusN,"hXiAllplus_%s",Group1.at(j).c_str());
    hXiAllplus = fs->make<TH1F>(hXiAllplusN,"; #xi_{all particles}^{+}; #N Events",1000,0,1.);
    hVectorXiAllplus.push_back(hXiAllplus);

    char hXiAllminusN[300];
    sprintf(hXiAllminusN,"hXiAllminus_%s",Group1.at(j).c_str());
    hXiAllminus = fs->make<TH1F>(hXiAllminusN,"; #xi_{all particles}^{-}; #N Events",1000,0,1.);
    hVectorXiAllminus.push_back(hXiAllminus);

    char hXiAllN[300];
    sprintf(hXiAllN,"hXiAll_%s",Group1.at(j).c_str());
    hXiAll = fs->make<TH1F>(hXiAllN,";#xi_{all particles}^{+} + #xi_{all particles}^{-}; #N Events",1000,0,1.);
    hVectorXiAll.push_back(hXiAll);

    char hXiDiffAllN[300];
    sprintf(hXiDiffAllN,"hXiDiffAll_%s",Group1.at(j).c_str());
    hXiDiffAll = fs->make<TH1F>(hXiDiffAllN,";|#xi_{all particles}^{+} - #xi_{all particles}^{-}|; #N Events",1000,0,1.);
    hVectorXiDiffAll.push_back(hXiDiffAll);

    char hGapSizeAllN[300];
    sprintf(hGapSizeAllN,"hGapSizeAll_%s",Group1.at(j).c_str());
    hGapSizeAll = fs->make<TH1F>(hGapSizeAllN,";#Delta#eta all particles; #N Events",1000,0,10);
    hVectorGapSizeAll.push_back(hGapSizeAll);

    char hGapSizeZN[300];
    sprintf(hGapSizeZN,"hGapSizeZ_%s",Group1.at(j).c_str());
    hGapSizeZ = fs->make<TH1F>(hGapSizeZN,";#Delta#eta, 60 < M_{ll} < 110 [GeV]; #N Events",1000,0,10);
    hVectorGapSizeZ.push_back(hGapSizeZ);

  }


}

void SDDYAnalyzer::endJob(){
}

void SDDYAnalyzer::fillHistos(int index){

  for (std::vector<std::string>::size_type j=0; j<etaAll.size(); j++){
    hVectorPartEta.at(index)->Fill(etaAll.at(j));
    hVectorPartPhi.at(index)->Fill(phiAll.at(j));
    hVectorPartPt.at(index)->Fill(ptAll.at(j));
    hVectorPartVx.at(index)->Fill(vxAll.at(j));
    hVectorPartVy.at(index)->Fill(vyAll.at(j));
    hVectorPartVz.at(index)->Fill(vzAll.at(j));
    hVectorPartPx.at(index)->Fill(px_genAll.at(j));
    hVectorPartPy.at(index)->Fill(py_genAll.at(j));
    hVectorPartPz.at(index)->Fill(pz_genAll.at(j));
    hVectorPartPf.at(index)->Fill(pf_genAll.at(j));
    hVectorPartEnergy.at(index)->Fill(energy_genAll.at(j));
    hVectorPartEnergyVsEta.at(index)->Fill(etaAll.at(j),energy_genAll.at(j));
    hVectorPartPDGId.at(index)->Fill(pdgIdAll.at(j));
  }

  hVectorHFPlusE.at(index)->Fill(sumHFPlusGEN);
  hVectorHFMinusE.at(index)->Fill(sumHFMinusGEN);
  hVectorCastorE.at(index)->Fill(sumCastorGEN);

  hVectorLepton1Eta.at(index)->Fill(l1eta);
  hVectorLepton1Phi.at(index)->Fill(l1phi);
  hVectorLepton1Pt.at(index)->Fill(l1pt);
  hVectorLepton1Energy.at(index)->Fill(l1energy);
  hVectorLepton1px.at(index)->Fill(l1px);
  hVectorLepton1py.at(index)->Fill(l1py);
  hVectorLepton1pz.at(index)->Fill(l1pz);
  hVectorLepton1pf.at(index)->Fill(l1pf);
  hVectorLepton1vx.at(index)->Fill(l1vx);
  hVectorLepton1vy.at(index)->Fill(l1vy);
  hVectorLepton1vz.at(index)->Fill(l1vz);

  hVectorLepton2Eta.at(index)->Fill(l2eta);
  hVectorLepton2Phi.at(index)->Fill(l2phi);
  hVectorLepton2Pt.at(index)->Fill(l2pt);
  hVectorLepton2Energy.at(index)->Fill(l2energy);
  hVectorLepton2px.at(index)->Fill(l2px);
  hVectorLepton2py.at(index)->Fill(l2py);
  hVectorLepton2pz.at(index)->Fill(l2pz);
  hVectorLepton2pf.at(index)->Fill(l2pf);
  hVectorLepton2vx.at(index)->Fill(l2vx);
  hVectorLepton2vy.at(index)->Fill(l2vy);
  hVectorLepton2vz.at(index)->Fill(l2vz);

  hVectorDileptonEta.at(index)->Fill(dibosonEta);
  hVectorDileptonPhi.at(index)->Fill(dibosonPhi);
  hVectorDileptonPt.at(index)->Fill(dibosonPt);
  hVectorDileptonM.at(index)->Fill(dibosonM);

  hVectorLeptonDeltaEta.at(index)->Fill(deltaeta);
  hVectorLeptonDeltaPhi.at(index)->Fill(deltaphi);
  hVectorLeptonDeltaPt.at(index)->Fill(deltapt);
  hVectorLeptonVertexD.at(index)->Fill(vertex_d);

  hVectortprotonplus.at(index)->Fill(t_plus);
  hVectortprotonminus.at(index)->Fill(t_minus);

  hVectorXiZplus.at(index)->Fill(xiZ_plus);
  hVectorXiZminus.at(index)->Fill(xiZ_minus);
  hVectorXidiffplus.at(index)->Fill(xi_diff_plus);
  hVectorXidiffminus.at(index)->Fill(xi_diff_minus);
  hVectorXiprotonplus.at(index)->Fill(xiProtonPlus);
  hVectorXiprotonminus.at(index)->Fill(xiProtonMinus);
  hVectorXiDiffproton.at(index)->Fill(xiDiffProton);
  hVectorXiAll.at(index)->Fill(xi);
  hVectorXiDiffAll.at(index)->Fill(xi_diff);
  hVectorXiAllplus.at(index)->Fill(xi_plus);
  hVectorXiAllminus.at(index)->Fill(xi_minus);
  hVectorGapSizeZ.at(index)->Fill(gapsizeZ);
  hVectorGapSizeAll.at(index)->Fill(gapsizeall);

}

//void SDDYAnalyzer::fillHistos

void SDDYAnalyzer::analyze(const edm::Event & ev, const edm::EventSetup&){

  bool debug_proton = true;

  sumHFMinusGEN = 0.; sumCastorGEN = 0.; sumHFPlusGEN = 0.;
  deltaeta = 0.; deltaphi = 0.; deltapt = 0.; vertex_d = 0.;
  genEPlusPz = 0.; genEMinusPz = 0.;
  xi = -999.; xi_minus = -999.; xi_plus = -999.; xiZ_plus = -999.; xiZ_minus = -999.; xiZ = -999.; xi_diff = -999.;
  xi_diff_plus = -999.; xi_diff_minus = -999.;
  l1eta = 0.; l1phi = 0.; l1pt = 0.; l1energy = 0.;
  l2eta = 0.; l2phi = 0.; l2pt = 0.; l2energy = 0.;
  l1px = 0.; l1py = 0.; l1pz = 0.; l1pf = 0.;
  l2px = 0.; l2py = 0.; l2pz = 0.; l2pf = 0.;
  l1vx = 0.; l1vy = 0.; l1vz = 0.;
  l2vx = 0.; l2vy = 0.; l2vz = 0.;
  xiDiffProton = -999.;
  xiProtonPlus = -999.;
  xiProtonMinus = -999.;
  proton_eta_plus = 0.; proton_phi_plus = 0.; proton_pt_plus = 0.; proton_energy_plus = 0.;
  proton_px_plus = 0.; proton_py_plus = 0.; proton_pz_plus = 0.; proton_pf_plus = 0.;
  proton_eta_minus = 0.; proton_phi_minus = 0.; proton_pt_minus = 0.; proton_energy_minus = 0.;
  proton_px_minus = 0.; proton_py_minus = 0.; proton_pz_minus = 0.; proton_pf_minus = 0.;
  dibosonEta = 0.; dibosonPhi = 0.; dibosonPt = 0.; dibosonM = 0.;
  diboson_px = 0.; diboson_py = 0.; diboson_pz = 0.; diboson_pf = 0., diboson_energy = 0.;
  t_plus = 0.; t_minus = -999.;
  pz_cut = 0.70*Ebeam_;
  gapsizeall = 0.; gapsizeZ = 0.;
  double s = (2*Ebeam_)*(2*Ebeam_);

  energy_genAll.clear();
  px_genAll.clear();
  py_genAll.clear();
  pz_genAll.clear();
  etaAll.clear();
  phiAll.clear();
  ptAll.clear();
  vxAll.clear();
  vyAll.clear();
  vzAll.clear();
  pf_genAll.clear();
  pt_genAll.clear();
  pdgIdAll.clear();

  lepton1 = false; lepton2 = false; dilepton = false, leptonAcceptance = false;
  single_gap = false; double_gap = false;
  zboson = false; zboson_diff = false;
  HF_CASTOR_gap = false;
  protonplus = false;
  protonminus = false;
  ptcut = false;

  // Generator Particles
  edm::Handle<reco::GenParticleCollection> genParticles;
  ev.getByLabel(genParticlesTag_, genParticles);

  reco::GenParticleCollection::const_iterator particle1 = genParticles->end();
  reco::GenParticleCollection::const_iterator particle2 = genParticles->end();

  reco::GenParticleCollection::const_iterator proton1 = genParticles->end();
  reco::GenParticleCollection::const_iterator proton2 = genParticles->end();

  double pz1max = 0.;
  double pz2min = 0.;
  double M_x = 0.;

  for(reco::GenParticleCollection::const_iterator genpart = genParticles->begin(); genpart != genParticles->end(); ++genpart){

    pf_gen = 0.;
    double pz = genpart->pz();

    if(genpart->status() == 1){ // only final state particles.

      // Identifying Protons
      if(genpart->pdgId() == 2212){

	// Proton Scattered Positive Side
	if(genpart->pz() > pz_cut){

	  if(pz > pz1max){
	    pz1max=pz;
	    proton1 = genpart;
	    proton_energy_plus = proton1->energy();
	    proton_px_plus = proton1->px();
	    proton_py_plus = proton1->py();
	    proton_pz_plus = proton1->pz();
	    proton_pt_plus = proton1->pt();
	    proton_eta_plus = proton1->eta();
	    proton_phi_plus = proton1->phi();
	    proton_pf_plus = sqrt(proton_px_plus*proton_px_plus+proton_py_plus*proton_py_plus+proton_pz_plus*proton_pz_plus);
	    protonplus=true;

	    xiProtonPlus = 1 - (proton_pz_plus/Ebeam_);
	    math::XYZTLorentzVector pi1(0.,0.,Ebeam_,Ebeam_);
	    math::XYZTLorentzVector pf1(proton_px_plus,proton_py_plus,proton_pz_plus,proton_energy_plus); // 4 momentum of p3
	    math::XYZTLorentzVector vec_t_plus = (pf1 - pi1);
	    t_plus=fabs(vec_t_plus.mag2());

	    if (debug_proton){
	      std::cout << "\n< Proton Plus >" << std::endl;
	      std::cout << "Quadrimomentum("<< proton_px_plus << "," << proton_py_plus << "," << proton_pz_plus << "," << proton_energy_plus << ") [GeV]"<< std::endl;
	      std::cout << "Eta: "<< proton_eta_plus << ", phi: " << proton_phi_plus << ", pT: " << proton_pt_plus << " [GeV]"<< std::endl;
	      std::cout << "t p-plus: " << t_plus << " [GeV]" << std::endl;
	      //std::cout << "t p-plus: " << proton1->momentum().perp2() << " [GeV]" << std::endl;
	      std::cout << "xi plus: " << xiProtonPlus << std::endl;
	      std::cout << "Vertex Rho: " << proton1->vertex().Rho() << std::endl;
	      std::cout << "Vertex Z = " << proton1->vertex().Z() << std::endl;
	      std::cout << "Vertex P("<< proton1->vx() << "," << proton1->vy() << "," <<  proton1->vz() << ") mm" << std::endl;
	    }
	  }
	}

	// Proton Scattered Negative Side

	if(genpart->pz() < -pz_cut){
	  if(pz < pz2min){
	    pz2min=pz;
	    proton2 = genpart;
	    proton_energy_minus = proton2->energy();
	    proton_px_minus = proton2->px();
	    proton_py_minus = proton2->py();
	    proton_pz_minus = proton2->pz();
	    proton_pt_minus = proton2->pt();
	    proton_eta_minus = proton2->eta();
	    proton_phi_minus = proton2->phi();
	    proton_pf_minus = sqrt(proton_px_minus*proton_px_minus+proton_py_minus*proton_py_minus+proton_pz_minus*proton_pz_minus);
	    protonminus=true;

	    xiProtonMinus = 1 + (proton_pz_minus/Ebeam_);
	    math::XYZTLorentzVector pi2(0.,0.,-Ebeam_,Ebeam_);
	    math::XYZTLorentzVector pf2(proton_px_minus,proton_py_minus,proton_pz_minus,proton_energy_minus); // 4 momentum of p3
	    math::XYZTLorentzVector vec_t_minus = (pf2 - pi2);
	    t_minus=fabs(vec_t_minus.mag2());

	    if (debug_proton){
	      std::cout << "\n< Proton Minus >" << std::endl;
	      std::cout << "Quadrimomentum("<< proton_px_minus << "," << proton_py_minus << "," << proton_pz_minus << "," << proton_energy_minus << ") [GeV]"<< std::endl;
	      std::cout << "Eta: "<< proton_eta_minus << ", phi: " << proton_phi_minus << ", pT: " << proton_pt_minus << "[GeV]"<< std::endl;
	      std::cout << "t p-minus: " << t_minus << " [GeV]" << std::endl;
	      //std::cout << "t p-minus: " << proton2->momentum().perp2() << " [GeV]" << std::endl;
	      std::cout << "xi minus: " << xiProtonMinus<< std::endl;
	      std::cout << "Vertex Rho: " << proton2->vertex().Rho() << std::endl;
	      std::cout << "Vertex Z = " << proton2->vertex().Z() << std::endl;
	      std::cout << "Vertex P("<< proton2->vx() << "," << proton2->vy() << "," <<  proton2->vz() << ") mm" << std::endl;
	    }
	  }
	}

	xiDiffProton = fabs(xiProtonMinus - xiProtonPlus);

      }

      if(!cmsAccept || (cmsAccept && (genpart->eta() <=5.2 || genpart->eta() >= -6.2))){

	// identifying leptons
	if((particle1 == genParticles->end())&&(abs(genpart->pdgId()) == abs(particle1Id_))) {particle1 = genpart;lepton1=true;continue;}
	if((particle2 == genParticles->end())&&(abs(genpart->pdgId()) == abs(particle2Id_))) {particle2 = genpart;lepton2=true;continue;}


	// Other particles, not proton. Detector noise approximation cut ~3 GeV.
	if (fabs(genpart->pz()) > 3. && genpart->pdgId() != 2212){
	  if (genpart->eta() <= -5.2 && genpart->eta() >= -6.2) sumCastorGEN += genpart->energy();
	  if (genpart->eta() <= -3 && genpart->eta() >= -5.2) sumHFMinusGEN += genpart->energy();
	  if (genpart->eta() >= 3 && genpart->eta() <= 5.2) sumHFPlusGEN += genpart->energy();

	  pf_gen = sqrt(genpart->px()*genpart->px()+genpart->py()*genpart->py()+genpart->pz()*genpart->pz());

	  energy_genAll.push_back(genpart->energy());
	  px_genAll.push_back(genpart->px());
	  py_genAll.push_back(genpart->py());
	  pz_genAll.push_back(genpart->pz());
	  etaAll.push_back(genpart->eta());
	  phiAll.push_back(genpart->phi());
	  ptAll.push_back(genpart->pt());
	  vxAll.push_back(genpart->vx());
	  vyAll.push_back(genpart->vy());
	  vzAll.push_back(genpart->vz());
	  pf_genAll.push_back(pf_gen);
	  pdgIdAll.push_back(genpart->pdgId());

	  genEPlusPz += (genpart->energy() + genpart->pz());
	  genEMinusPz += (genpart->energy() - genpart->pz());
	  M_x += genpart->mass();

	}
      }
    }

  }

  if (lepton1 && lepton2) dilepton = true;
  if (protonplus || protonminus) single_gap = true;
  if (protonplus && protonminus) double_gap = true;

  // Xi all particles
  xi_minus = genEMinusPz/(2*Ebeam_); 
  xi_plus = genEPlusPz/(2*Ebeam_);
  xi = xi_minus+xi_plus;
  xi_diff = fabs(xi_minus-xi_plus);

  // GapSize all particles
  gapsizeall = log(s/M_x);

  if (debug){
    std::cout << "\n< Xi Particles >" << std::endl;
    std::cout << "Xi Plus: " << xi_plus << std::endl;
    std::cout << "Xi Minus: "<< xi_minus << std::endl;
    std::cout << "Xi Diff: " << xi_diff << std::endl;
  }


  // Fill All CMS Particles and Proton
  SDDYAnalyzer::fillHistos(0);

  if(dilepton){

    // distance vertex dilepton
    vertex_d = sqrt( (particle2->vx()-particle1->vx())*(particle2->vx()-particle1->vx()) + (particle2->vy()-particle1->vy())*(particle2->vy()-particle1->vy()) + (particle2->vz()-particle1->vz())*(particle2->vz()-particle1->vz()));

    l1eta = particle1->eta(); l1phi = particle1->phi(), l1pt = particle1->pt(), l1energy = particle1->energy();
    l2eta = particle2->eta(); l2phi = particle2->phi(), l2pt = particle2->pt(), l2energy = particle2->energy();

    l1px = particle1->px(); l1py = particle1->py(); l1pz = particle1->pz();
    l2px = particle2->px(); l2py = particle2->py(); l2pz = particle2->pz();

    l1pf = sqrt(l1px*l1px+l1py*l1py+l1pz*l1pz);
    l2pf = sqrt(l2px*l2px+l2py*l2py+l2pz*l2pz);

    l1vx = particle1->vx();  l1vy = particle1->vy(); l1vz = particle1->vz();
    l2vx = particle2->vx();  l2vy = particle2->vy(); l2vz = particle2->vz();

    deltaeta = fabs(particle1->eta()-particle2->eta()); deltapt = fabs(particle1->pt()-particle2->pt());
    deltaphi = fabs(particle1->phi() - particle2->phi());
    if(deltaphi > M_PI){deltaphi = 2.0*M_PI - deltaphi;}

    math::XYZTLorentzVector myboson(particle1->px() + particle2->px(),
	particle1->py() + particle2->py(),
	particle1->pz() + particle2->pz(),
	particle1->energy() + particle2->energy());

    dibosonEta = myboson.eta();
    dibosonPhi = myboson.phi();
    dibosonPt = myboson.pt();
    dibosonM = myboson.M();

    diboson_px = myboson.px();
    diboson_py = myboson.py();
    diboson_pz = myboson.pz();
    diboson_energy = myboson.energy();
    diboson_pf = sqrt(diboson_px*diboson_px+diboson_py*diboson_py+diboson_pz*diboson_pz);

    // GapSize all particles
    gapsizeZ = log(s/myboson.M());

    // Spectrum selection
    if(myboson.M() >= 60 && myboson.M() <= 110) {
      zboson = true;
    }
    // xi Z positive scattering and xi all particles without Z
    if (myboson.eta() > 0) {
      xiZ_plus = (particle1->energy() + particle2->energy() + particle1->pz() + particle2->pz())/2*Ebeam_;
      xi_diff_plus = fabs(xi_plus - xiZ_plus);
    }
    // xi Z negative scattering and xi all particles without Z
    if (myboson.eta() < 0 ){
      xiZ_minus = (particle1->energy() + particle2->energy() - particle1->pz() - particle2->pz())/2*Ebeam_;
      xi_diff_minus = fabs(xi_minus - xiZ_minus);
    }

    // pT lepton cut (CMS approximation)
    if(particle1->pt() > 10. && particle2->pt() > 10.) ptcut = true;

    // Eta Acceptance CMS
    if((particle1->eta() > -2.5 && particle1->eta()< 2.5) && (particle2->eta() > -2.5 && particle2->eta()< 2.5) ) leptonAcceptance = true;

    if (debug){
      std::cout << "< Dilepton >" << std::endl;
      std::cout << "Lepton1("<< l1px << "," << l1py << "," << l1pz << "," << l1energy << ") [GeV]"<< std::endl;
      std::cout << "Eta: "<< l1eta << ", phi: " << l1phi << ", pT: " << l1pt << "[GeV], pf: " << l1pf << "[GeV]"<<std::endl;
      std::cout << "Vertex P("<< l1vx << "," << l1vy << "," <<  l1vz << ") mm" << std::endl;
      std::cout << "\nLepton2("<< l2px << "," << l2py << "," << l2pz << "," << l2energy << ") [GeV]"<< std::endl;
      std::cout << "Eta: "<< l2eta << ", phi: " << l2phi << ", pT: " << l2pt << "[GeV], pf: " << l2pf << "[GeV]"<< std::endl;
      std::cout << "Vertex P("<< l2vz << "," << l2vy << "," <<  l2vz << ") mm" << std::endl;
      std::cout << "Boson Z("<< diboson_px << "," << diboson_py << "," << diboson_pz << "," << diboson_energy << ") [GeV]" << std::endl;
      std::cout << "Eta: "<< dibosonEta << ", phi: " << dibosonPhi << ", pT: " << dibosonPt << "[GeV], pf: " << diboson_pf << "[GeV]"<<std::endl;
      std::cout << "Xi Z, plus: " << xiZ_plus << std::endl;
      std::cout << "Xi Z, Minus: "<< xiZ_minus << std::endl;
      std::cout << "Xi_all - Xi_z (plus): " << xi_diff_plus << std::endl;
      std::cout << "Xi_all - Xi_z (minus): " << xi_diff_minus << std::endl;
    }

  }

  // Negative very restricted GAP?
  if(sumHFMinusGEN == 0. && sumCastorGEN == 0.) {
    HF_CASTOR_gap = true;
  }

  // Fill Histograms

  // Dilepton CMS
  if (leptonAcceptance) SDDYAnalyzer::fillHistos(1);

  // Dilepton CMS and pTCut
  if (leptonAcceptance && ptcut) SDDYAnalyzer::fillHistos(2);

  // CMS Boson Z full selection
  if (leptonAcceptance && ptcut && zboson) SDDYAnalyzer::fillHistos(3);

  // CMS Boson Z less restricted gap cut
  if (leptonAcceptance && ptcut && zboson && sumCastorGEN==0) SDDYAnalyzer::fillHistos(4);

  // CMS Boson Z very restricted gap cut
  if (leptonAcceptance && ptcut && zboson && HF_CASTOR_gap) SDDYAnalyzer::fillHistos(5);

}	

DEFINE_FWK_MODULE(SDDYAnalyzer);
