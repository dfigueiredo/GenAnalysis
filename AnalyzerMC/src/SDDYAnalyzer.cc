////////// Header section /////////////////////////////////////////////
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "TH1F.h"

class SDDYAnalyzer: public edm::EDAnalyzer {
  public:
    /// Constructor
    SDDYAnalyzer(const edm::ParameterSet& pset);

    /// Destructor
    virtual ~SDDYAnalyzer();

    // Operations

    void analyze(const edm::Event & event, const edm::EventSetup& eventSetup);
    void fillHistos(int index);

    //virtual void beginJob(const edm::EventSetup& eventSetup) ;
    virtual void beginJob(); 
    virtual void endJob() ;

  private:
    // Input from cfg file
    edm::InputTag genParticlesTag_;
    double Ebeam_;
    int particle1Id_;
    int particle2Id_;

    // Histograms
    std::vector<TH1F*> hVectorPartEta;
    std::vector<TH1F*> hVectorPartPt;
    std::vector<TH1F*> hVectorPartPhi;
    std::vector<TH1F*> hVectorCMSEta;
    std::vector<TH1F*> hVectorCMSPt;
    std::vector<TH1F*> hVectorCMSPhi;
    std::vector<TH1F*> hVectorCMSHFMinusE;
    std::vector<TH1F*> hVectorCMSHFPlusE;
    std::vector<TH1F*> hVectorCMSCastorE;
    std::vector<TH1F*> hVectorCMSEnergyvsEta;
    std::vector<TH1F*> hVectorCMSBosonPt;
    std::vector<TH1F*> hVectorCMSBosonEta;
    std::vector<TH1F*> hVectorCMSBosonPhi;
    std::vector<TH1F*> hVectorCMSBosonM;
    
    
    TH1F* hPartEta;
    TH1F* hPartPt;
    TH1F* hPartPhi;
    TH1F* hCMSEta;
    TH1F* hCMSPt;
    TH1F* hCMSPhi;
    TH1F* hCMSHFMinusE;
    TH1F* hCMSHFPlusE;
    TH1F* hCMSCastorE;
    TH1F* hCMSEnergyvsEta;
    TH1F* hCMSBosonPt;
    TH1F* hCMSBosonEta;
    TH1F* hCMSBosonPhi;
    TH1F* hCMSBosonM;  
    TH1F* hEntriesCut;

    // Counters
    int nevents;
    int zboson;
    int zbosondiff;
    int index;
    
    // Histogram variables
    double sumHFMinusGEN;
    double sumCastorGEN;
    double sumHFPlusGEN;
    double px_gen;
    double py_gen;
    double pz_gen;
    double energy_gen;
    std::vector <std::string> Folders;

};

////////// Source code ////////////////////////////////////////////////
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

/// Constructor
SDDYAnalyzer::SDDYAnalyzer(const edm::ParameterSet& pset)
{
  genParticlesTag_ = pset.getParameter<edm::InputTag>("GenParticleTag");
  particle1Id_ = pset.getParameter<int>("Particle1Id");
  particle2Id_ = pset.getParameter<int>("Particle2Id");
  Ebeam_ = pset.getParameter<double>("EBeam");

  debug = pset.getUntrackedParameter<bool>("debug",false);
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
  zboson = 0;
  zbosondiff = 0;
	
  edm::Service<TFileService> fs;
  TH1::SetDefaultSumw2(true);

  std::string step0 = "without_cuts";
  std::string step1 = "CMS";
  std::string step2 = "CMS_step1";
  std::string step3 = "CMS_step2";
  std::string step4 = "CMS_step3";

  Folders.push_back(step0);
  Folders.push_back(step1);
  Folders.push_back(step2);
  Folders.push_back(step3);
  Folders.push_back(step4);

  for (std::vector<std::string>::size_type j=0; j<Folders.size(); j++){

    char hPartEtaN[300];
    sprintf(hPartEtaN,"hPartEta_%s",Folders.at(j).c_str());
    hPartEta = fs->make<TH1F>(hPartEtaN,"Title...",2000,-10.,10.);
    hVectorPartEta.push_back(hPartEta);
    d
  }
/*
  //Booking Histograms
  hPartEta = fs->make<TH1F>("hPartEta","hPartEta",2000,-10.,10.);
  hPartPt = fs->make<TH1F>("hPartPt","hPartPt",100,0.,100.);
  hPartPhi = fs->make<TH1F>("hPartPhi","hPartPhi",50,-3.141592,3.141592);
  hPartCMSEta = fs->make<TH1F>("hPartCMSEta","hPartCMSEta",2000,-10.,10.);
  hPartCMSPt = fs->make<TH1F>("hPartCMSPt","hPartCMSPt",100,0.,100.);
  hPartCMSPhi = fs->make<TH1F>("hPartCMSPhi","hPartCMSPhi",50,-3.141592,3.141592);
  hEnergyvsEtaCMS = fs->make<TH1F>("hEnergyvsEtaCMS","hEnergyvsEtaCMS",100,-15.0,15.0);
  hHFMinusEGEN = fs->make<TH1F>("hHFMinusEGEN","hHFMinusGen",1000,0,1000);
  hHFPlusEGEN = fs->make<TH1F>("hHFMinusEGEN","hHFMinusGen",1000,0,1000);
  hCastorEGEN = fs->make<TH1F>("hCastorEGEN","hCastorEGen",1000,0,1000);
  hBosonPt = fs->make<TH1F>("hBosonPt","hBosonPt",100,0.,50.);
  hBosonEta = fs->make<TH1F>("hBosonEta","hBosonEta",100,-5.,5.);
  hBosonPhi = fs->make<TH1F>("hBosonPhi","hBosonPhi",50,-3.141592,3.141592);
  hBosonM = fs->make<TH1F>("hBosonM","hBosonM",100,40.,150.);
  hBosonPtDiff = fs->make<TH1F>("hBosonPtDiff","hBosonPtDiff",100,0.,50.);
  hBosonEtaDiff = fs->make<TH1F>("hBosonEtaDiff","hBosonEtaDiff",100,-5.,5.);
  hBosonPhiDiff = fs->make<TH1F>("hBosonPhiDiff","hBosonPhiDiff",50,-3.141592,3.141592);
  hBosonMDiff = fs->make<TH1F>("hBosonMDiff","hBosonMDiff",100,40.,150.);
  hEntriesCut = fs->make<TH1F>("hEntriesCut","hEntriesCut",10,0.,10.);
  hEnergyvsEtaCMS->Scale(1/(float)nevents);
  */

}

void SDDYAnalyzer::endJob(){

  //Counters	
  hEntriesCut->SetBinContent(2,nevents);
  hEntriesCut->SetBinContent(4,zboson);
  hEntriesCut->SetBinContent(6,zbosondiff);
  
}

void SDDYAnalyzer::fillHistos(int index){
	
	hVectorPartEta.at(index)->Fill(10);
	
	
}

void SDDYAnalyzer::analyze(const edm::Event & ev, const edm::EventSetup&){

  //switches
  bool lepton1 = false;
  bool lepton2 = false;
  debug printout = false;
  
 sumHFMinusGEN = 0.;
 sumCastorGEN = 0.;
 sumHFPlusGEN = 0.;
 px_gen = 0.;
 py_gen = 0.;
 pz_gen = 0.;
 energy_gen = 0;;
  nevents++;

  // Generator Particles
  edm::Handle<reco::GenParticleCollection> genParticles;
  ev.getByLabel(genParticlesTag_, genParticles);

  reco::GenParticleCollection::const_iterator particle1 = genParticles->end();
  reco::GenParticleCollection::const_iterator particle2 = genParticles->end();
  
  reco::GenParticleCollection::const_iterator proton1 = genParticles->end();
  reco::GenParticleCollection::const_iterator proton2 = genParticles->end();



  for(reco::GenParticleCollection::const_iterator genpart = genParticles->begin(); genpart != genParticles->end(); ++genpart){

    if (printout) std::cout << ">>>>>>> pid,status,px,py,px,e= "  << genpart->pdgId() << " , " << genpart->status() << " , " << genpart->px() << " , " << genpart->py() << " , " << genpart->pz() << " , " << genpart->energy() << std::endl;

    if(genpart->status() != 1) continue; // check if genparticle survives.

    energy_gen = genpart->energy();
    px_gen = genpart->px();
    py_gen = genpart->py();
    pz_gen = genpart->pz();
    
    fillHistos(0);
    fillHistos(1);
    fillHistos(2);
    fillHistos(3);
    fillHistos(4);

    hPartPt->Fill(genpart->pt());
    hPartEta->Fill(genpart->eta());
    hPartPhi->Fill(genpart->phi());

    if (genpart->eta() > 5.2 || genpart->eta() < -6.2) continue;

    hPartCMSPt->Fill(genpart->pt());
    hPartCMSEta->Fill(genpart->eta());
    hPartCMSPhi->Fill(genpart->phi());
    hEnergyvsEtaCMS->Fill(genpart->eta(),genpart->energy());

    if (genpart->eta() <= -5.2 && genpart->eta() >= -6.2) sumCastorGEN += genpart->energy();
    if (genpart->eta() <= -3 && genpart->eta() >= -5.2) sumHFMinusGEN += genpart->energy();
    if (genpart->eta() >= 3 && genpart->eta() <= 5.2) sumHFPlusGEN += genpart->energy();

    if((particle1 == genParticles->end())&&(abs(genpart->pdgId()) == abs(particle1Id_))) {particle1 = genpart;lepton1=true;continue;}
    if((particle2 == genParticles->end())&&(abs(genpart->pdgId()) == abs(particle2Id_))) {particle2 = genpart;lepton2=true;continue;}

  }

  if (lepton1 && lepton2){
    if (debug) std::cout << "Di lepton selected" << std::endl;
    if (particle1->pt() > 10. && particle2->pt() > 10.) {
      if (debug) std::cout << "Di lepton pT cut" << std::endl;
      if((particle1->eta() > -2.5 && particle1->eta()< 2.5) && (particle2->eta() > -2.5 && particle2->eta()< 2.5) ){
	if (debug) std::cout << "Di lepton CMS acceptance cut" << std::endl;

	math::XYZTLorentzVector myboson(particle1->px() + particle2->px(),
	    particle1->py() + particle2->py(),
	    particle1->pz() + particle2->pz(),
	    particle1->energy() + particle2->energy());

	if (myboson.M() >= 60 && myboson.M() <= 110){
          zboson++;
	  hBosonPt->Fill(myboson.pt());
	  hBosonEta->Fill(myboson.eta());
	  hBosonPhi->Fill(myboson.phi());
	  hBosonM->Fill(myboson.M());
	  hHFMinusEGEN->Fill(sumHFMinusGEN);
	  hHFPlusEGEN->Fill(sumHFPlusGEN);
	  hCastorEGEN->Fill(sumCastorGEN);
	}

	if ((myboson.M() >= 60 && myboson.M() <= 110) && sumHFMinusGEN == 0. && sumCastorGEN == 0.){
          zbosondiff++;
	  hBosonPtDiff->Fill(myboson.pt());
	  hBosonEtaDiff->Fill(myboson.eta());
	  hBosonPhiDiff->Fill(myboson.phi());
	  hBosonMDiff->Fill(myboson.M());
	}

      }
    }
  }
}	

DEFINE_FWK_MODULE(SDDYAnalyzer);
