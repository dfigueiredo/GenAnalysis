//---------------------------
// GEN Diffractive Z Analyzer
//---------------------------

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

    std::vector<TH1F*> hVectorDiLeptonVertexD;
    std::vector<TH1F*> hVectorHFPlusE;
    std::vector<TH1F*> hVectorHFMinusE;
    std::vector<TH1F*> hVectorCastorE;
    std::vector<TH1F*> hVectorDileptonEta;
    std::vector<TH1F*> hVectorDileptonPhi;
    std::vector<TH1F*> hVectorDileptonPt;
    std::vector<TH1F*> hVectorDiLeptonM;

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

    TH1F* hDiLeptonVertexD;
    TH1F* hHFPlusE;
    TH1F* hHFMinusE;
    TH1F* hCastorE;
    TH1F* hDileptonEta;
    TH1F* hDileptonPhi;
    TH1F* hDileptonPt;
    TH1F* hDiLeptonM;

    TH1F* hEntries;

    // Counters
    int nevents;
    int zboson_counter;
    int zbosondiff_counter;
    int index;

    // Histogram variables
    double sumHFMinusGEN;
    double sumCastorGEN;
    double sumHFPlusGEN;
    double pf_gen;
    double deltaeta;
    double deltaphi;
    double deltapt;
    double vertex_d;
    double xi;
    double xi_minus;
    double xi_plus;
    double l1eta;
    double l1phi;
    double l1pt;
    double l1energy;
    double l1vx;
    double l1vy;
    double l1vz;
    double l2eta;
    double l2phi;
    double l2pt;
    double l2energy;
    double l2vx;
    double l2vy;
    double l2vz;
    double l1px;
    double l1py;
    double l1pz;
    double l1pf;
    double l2px;
    double l2py;
    double l2pz;
    double l2pf;
    double xiZ;
    double xiProtonPlus;
    double xiProtonMinus;
    double proton_px_plus;
    double proton_py_plus;
    double proton_pz_plus;
    double proton_pf_plus;
    double proton_px_minus;
    double proton_py_minus;
    double proton_pz_minus;
    double proton_pf_minus;
    double proton_energy_plus;
    double proton_energy_minus;
    double pz_cut;
    double genEPlusPz;
    double genEMinusPz;
    double t_plus;
    double t_minus;

    bool lepton1;
    bool lepton2;
    bool single_gap;
    bool double_gap;
    bool zboson;
    bool zboson_diff;
    bool HF_CASTOR_gap;
    bool dilepton;
    bool protonplus;
    bool protonminus;
    bool leptonAccept;
    bool ptcut;

    std::vector <std::string> Group1;
    std::vector <std::string> Group2;

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

  // Counters Start
  nevents = 0;
  zboson_counter = 0;
  zbosondiff_counter = 0;

  edm::Service<TFileService> fs;
  TH1::SetDefaultSumw2(true);

  std::string step0 = "without_cuts_AllParticles";
  std::string step1 = "CMS_AllParticles";


  Group1.push_back(step0);
  Group1.push_back(step1);

  hEntries = fs->make<TH1F>("hEntries","hEntries",10,0.,10.);

  for (std::vector<std::string>::size_type j=0; j<Group1.size(); j++){

    char hPartEtaN[300];
    sprintf(hPartEtaN,"hPartEta_%s",Group1.at(j).c_str());
    hPartEta = fs->make<TH1F>(hPartEtaN,"Title...",20,0,500.);
    hVectorPartEta.push_back(hPartEta);

    char hPartPhiN[300];
    sprintf(hPartPhiN,"hPartPhi_%s",Group1.at(j).c_str());
    hPartPhi = fs->make<TH1F>(hPartPhiN,"Title...",20,0,500.);
    hVectorPartPhi.push_back(hPartPhi);

    char hPartPtN[300];
    sprintf(hPartPtN,"hPartPt_%s",Group1.at(j).c_str());
    hPartPt = fs->make<TH1F>(hPartPtN,"Title...",20,0,500.);
    hVectorPartPt.push_back(hPartPt);

    char hPartVxN[300];
    sprintf(hPartVxN,"hPartVx_%s",Group1.at(j).c_str());
    hPartVx = fs->make<TH1F>(hPartVxN,"Title...",20,0,500.);
    hVectorPartVx.push_back(hPartVx);

    char hPartVyN[300];
    sprintf(hPartVyN,"hPartVy_%s",Group1.at(j).c_str());
    hPartVy = fs->make<TH1F>(hPartVyN,"Title...",20,0,500.);
    hVectorPartVy.push_back(hPartVy);

    char hPartVzN[300];
    sprintf(hPartVzN,"hPartVz_%s",Group1.at(j).c_str());
    hPartVz = fs->make<TH1F>(hPartVzN,"Title...",20,0,500.);
    hVectorPartVz.push_back(hPartVz);

    char hPartPxN[300];
    sprintf(hPartPxN,"hPartPx_%s",Group1.at(j).c_str());
    hPartPx = fs->make<TH1F>(hPartPxN,"Title...",20,0,500.);
    hVectorPartPx.push_back(hPartPx);

    char hPartPyN[300];
    sprintf(hPartPyN,"hPartPy_%s",Group1.at(j).c_str());
    hPartPy = fs->make<TH1F>(hPartPyN,"Title...",20,0,500.);
    hVectorPartPy.push_back(hPartPy);

    char hPartPzN[300];
    sprintf(hPartPzN,"hPartPz_%s",Group1.at(j).c_str());
    hPartPz = fs->make<TH1F>(hPartPzN,"Title...",20,0,500.);
    hVectorPartPz.push_back(hPartPz);

    char hPartPfN[300];
    sprintf(hPartPfN,"hPartPf_%s",Group1.at(j).c_str());
    hPartPf = fs->make<TH1F>(hPartPfN,"Title...",20,0,500.);
    hVectorPartPf.push_back(hPartPf);

    char hPartEnergyN[300];
    sprintf(hPartEnergyN,"hPartEnergy_%s",Group1.at(j).c_str());
    hPartEnergy = fs->make<TH1F>(hPartEnergyN,"Title...",20,0,500.);
    hVectorPartEnergy.push_back(hPartEnergy);

    char hPartEnergyVsEtaN[300];
    sprintf(hPartEnergyVsEtaN,"hPartEnergyVsEta_%s",Group1.at(j).c_str());
    hPartEnergyVsEta = fs->make<TH1F>(hPartEnergyVsEtaN,"Title...",20,0,500.);
    hVectorPartEnergyVsEta.push_back(hPartEnergyVsEta);

    char hPartPDGIdN[300];
    sprintf(hPartPDGIdN,"hPartPDGId_%s",Group1.at(j).c_str());
    hPartPDGId = fs->make<TH1F>(hPartPDGIdN,"Title...",20,0,500.);
    hVectorPartPDGId.push_back(hPartPDGId);

  }

  /*
     for (std::vector<std::string>::size_type j=0; j<Group2.size(); j++){

     char hDiLeptonVertexDN[300];
     sprintf(hDiLeptonVertexDN,"hDiLeptonVertexD_%s",Group2.at(j).c_str());
     hDiLeptonVertexD = fs->make<TH1F>(hDiLeptonVertexDN,"Title...",20,0,500.);
     hVectorDiLeptonVertexD.push_back(hDiLeptonVertexD);

     char hHFPlusEN[300];
     sprintf(hHFPlusEN,"hHFPlusE_%s",Group2.at(j).c_str());
     hHFPlusE = fs->make<TH1F>(hHFPlusE,"Title...",20,0,500.);
     hVectorHFPlusE.push_back(hHFPlusE);

     char hHFMinusEN[300];
     sprintf(hHFMinusEN,"hHFMinusE_%s",Group2.at(j).c_str());
     hHFMinusE = fs->make<TH1F>(hHFMinusE,"Title...",20,0,500.);
     hVectorHFMinusE.push_back(hHFMinusE);

     char hCastorEN[300];
     sprintf(hCastorEN,"hCastorE_%s",Group2.at(j).c_str());
     hCastorE = fs->make<TH1F>(hCastorE,"Title...",20,0,500.);
     hVectorCastorE.push_back(hCastorE);

     char hDileptonEtaN[300];
     sprintf(hDileptonEtaN,"hDiLeptonEta_%s",Group2.at(j).c_str());
     hDiLeptonEta = fs->make<TH1F>(hDiLeptonEta,"Title...",20,0,500.);
     hVectorDiLeptonEta.push_back(hDiLeptonEta);

     char hDileptonPhiN[300];
     sprintf(hDileptonPhiN,"hDiLeptonPhi_%s",Group2.at(j).c_str());
     hDiLeptonPhi = fs->make<TH1F>(hDiLeptonPhi,"Title...",20,0,500.);
     hVectorDiLeptonPhi.push_back(hDiLeptonPhi);

     char hDileptonPtN[300];
     sprintf(hDileptonPtN,"hDiLeptonPt_%s",Group2.at(j).c_str());
     hDiLeptonPt = fs->make<TH1F>(hDiLeptonPt,"Title...",20,0,500.);
     hVectorDiLeptonPt.push_back(hDiLeptonPt);

     char hDileptonMN[300];
     sprintf(hDileptonMN,"hDiLeptonM_%s",Group2.at(j).c_str());
     hDiLeptonM = fs->make<TH1F>(hDiLeptonM,"Title...",20,0,500.);
     hVectorDiLeptonM.push_back(hDiLeptonM);

     }
   */

}

void SDDYAnalyzer::endJob(){

  //Counters	
  /*
  hEntries->SetBinContent(2,nevents);
  hEntries->SetBinContent(4,zboson_counter);
  hEntries->SetBinContent(6,zbosondiff_counter);
  */

}

void SDDYAnalyzer::fillHistos(int index){

  if (debug) std::cout << "Filling histograms..." << std::endl;

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

}

//void SDDYAnalyzer::fillHistos

void SDDYAnalyzer::analyze(const edm::Event & ev, const edm::EventSetup&){

  sumHFMinusGEN = 0.;
  sumCastorGEN = 0.;
  sumHFPlusGEN = 0.;
  proton_energy_plus = 0.;
  proton_energy_minus = 0.;
  deltaeta = 0.;
  deltaphi = 0.;
  deltapt = 0.;
  vertex_d = 0.;
  xi = 0.;
  xi_minus = 0.;
  xi_plus = 0.;
  l1eta = 0.;
  l1phi = 0.;
  l1pt = 0.;
  l1energy = 0.;
  l1vx = 0.;
  l1vy = 0.;
  l1vz = 0.;
  l2eta = 0.;
  l2phi = 0.;
  l2pt = 0.;
  l2energy = 0.;
  l2vx = 0.;
  l2vy = 0.;
  l2vz = 0.;
  l1px = 0.;
  l1py = 0.;
  l1pz = 0.;
  l1pf = 0.;
  l2px = 0.;
  l2py = 0.;
  l2pz = 0.;
  l2pf = 0.;
  xiZ = 0.;
  xiProtonPlus = -999;
  xiProtonMinus = -999;
  proton_px_plus = 0.;
  proton_py_plus = 0.;
  proton_pz_plus = 0.;
  proton_pf_plus = 0.;
  proton_px_minus = 0.;
  proton_py_minus = 0.;
  proton_pz_minus = 0.;
  proton_pf_minus = 0.;
  t_plus = 0.;
  t_minus = 0.;
  pz_cut = 0.7*Ebeam_;
  genEPlusPz = 0.;
  genEMinusPz = 0.;
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

  lepton1 = false;
  lepton2 = false;
  single_gap = false;
  double_gap = false;
  zboson = false;
  zboson_diff = false;
  HF_CASTOR_gap = false;
  dilepton = false;
  protonplus = false;
  protonminus = false;
  ptcut = false;

  nevents++;

  // Generator Particles
  edm::Handle<reco::GenParticleCollection> genParticles;
  ev.getByLabel(genParticlesTag_, genParticles);

  reco::GenParticleCollection::const_iterator particle1 = genParticles->end();
  reco::GenParticleCollection::const_iterator particle2 = genParticles->end();

  reco::GenParticleCollection::const_iterator proton1 = genParticles->end();
  reco::GenParticleCollection::const_iterator proton2 = genParticles->end();

  int counter_proton1 = 0;
  int counter_proton2 = 0;
  int counter_total_proton = 0;

  for(reco::GenParticleCollection::const_iterator genpart = genParticles->begin(); genpart != genParticles->end(); ++genpart){

    pf_gen = 0.;

    if (debug) std::cout << ">>>>>>> pid,status,px,py,px,e= "  << genpart->pdgId() << " , " << genpart->status() << " , " << genpart->px() << " , " << genpart->py() << " , " << genpart->pz() << " , " << genpart->energy() << std::endl;

    if(genpart->status() != 1) continue; // only final state particles.
    if(cmsAccept && genpart->eta() >=5.2 && genpart->eta() <= -6.2) continue; // CMS Acceptance

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

    // identifying leptons
    if((particle1 == genParticles->end())&&(abs(genpart->pdgId()) == abs(particle1Id_))) {particle1 = genpart;lepton1=true;continue;}
    if((particle2 == genParticles->end())&&(abs(genpart->pdgId()) == abs(particle2Id_))) {particle2 = genpart;lepton2=true;continue;}

    // identifying protons
    if((genpart->pdgId() == 2212)){
        counter_total_proton++;
      if(genpart->pz() > 0.){
        counter_proton1++;
	proton1 = genpart;
	proton_energy_plus = genpart->energy();
	proton_px_plus = genpart->px(); 
	proton_py_plus = genpart->py();
	proton_pz_plus = genpart->pz();
	proton_pf_plus = sqrt(genpart->px()*genpart->px()+genpart->py()*genpart->py()+genpart->pz()*genpart->pz());
	protonplus=true;
      }
    } else if((genpart->pdgId() == 2212)){ 
      if(genpart->pz() < 0.){
        counter_proton2++;
	proton2 = genpart;
	proton_energy_minus = genpart->energy();
	proton_px_minus = genpart->px();
	proton_py_minus = genpart->py();
	proton_pz_minus = genpart->pz();
	proton_pf_minus = sqrt(genpart->px()*genpart->px()+genpart->py()*genpart->py()+genpart->pz()*genpart->pz());
	protonminus = true;
      }
    }

    if (fabs(genpart->pz()) < pz_cut){
      genEPlusPz += (genpart->energy() + genpart->pz());
      genEMinusPz += (genpart->energy() - genpart->pz());         
    }

    if (genpart->eta() <= -5.2 && genpart->eta() >= -6.2) sumCastorGEN += genpart->energy();
    if (genpart->eta() <= -3 && genpart->eta() >= -5.2) sumHFMinusGEN += genpart->energy();
    if (genpart->eta() >= 3 && genpart->eta() <= 5.2) sumHFPlusGEN += genpart->energy();

  }

  SDDYAnalyzer::fillHistos(0);

  if(protonplus){
    if(debug) std::cout << "Proton 1: " << proton1->pt() << " " << proton1->eta() << " " << proton1->phi() << std::endl;
    if(proton_pz_plus > 0.){
      xiProtonPlus = ( 1 - (proton_pz_plus/Ebeam_) );
      math::XYZTLorentzVector p2(0,0,Ebeam_,Ebeam_);
      math::XYZTLorentzVector p3(proton_px_plus,proton_py_plus,proton_pz_plus,proton_energy_plus); // 4 momentum of p3
      math::XYZTLorentzVector vec_t = (p3 - p2);
      t_plus=vec_t.M2();
    }
  }        

  if(protonminus){
    if(debug) std::cout << "Proton 2: " << proton2->pt() << " " << proton2->eta() << " " << proton2->phi() << std::endl;        
    if(proton_pz_minus < 0.){
      xiProtonMinus = (proton_pz_minus < 0.) ? ( 1 + (proton_pz_minus/Ebeam_) ) : -1.;
      math::XYZTLorentzVector p1(0,0,-Ebeam_,Ebeam_);
      math::XYZTLorentzVector pm(proton_px_minus,proton_py_minus,proton_pz_minus,proton_energy_minus); // 4 momentum of p3
      math::XYZTLorentzVector vec_t = (pm - p1);
      t_minus=vec_t.M2();
    }
  }

  if (lepton1 && lepton2) dilepton = true;
  if (protonplus || protonminus) single_gap = true;
  if (protonplus && protonminus) double_gap = true;

  if(dilepton){
    math::XYZTLorentzVector myboson(particle1->px() + particle2->px(),
	particle1->py() + particle2->py(),
	particle1->pz() + particle2->pz(),
	particle1->energy() + particle2->energy());
    if(particle1->pt() > 10. && particle2->pt() > 10.) ptcut = true;
    if((particle1->eta() > -2.5 && particle1->eta()< 2.5) && (particle2->eta() > -2.5 && particle2->eta()< 2.5) ) leptonAccept = true;
    if(myboson.M() >= 60 && myboson.M() <= 110) {
      zboson_counter++;
      zboson = true;
    }
    if(sumHFMinusGEN == 0. && sumCastorGEN == 0.) {
      zbosondiff_counter++;
      HF_CASTOR_gap = true;
    }

  } 

  std::cout << "Counter Total Proton: " << counter_total_proton << std::endl;
  std::cout << "Counter Proton 1: " << counter_proton1 << std::endl;
  std::cout << "Counter Proton 2: " << counter_proton2 << std::endl;

}	

DEFINE_FWK_MODULE(SDDYAnalyzer);
