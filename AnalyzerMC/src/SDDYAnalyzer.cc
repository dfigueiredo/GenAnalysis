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
    TH1F* hPartEta;
    TH1F* hPartPt;
    TH1F* hPartPhi;
    TH1F* hPartCMSEta;
    TH1F* hPartCMSPt;
    TH1F* hPartCMSPhi;
    TH1F* hHFMinusEGEN;
    TH1F* hHFPlusEGEN;
    TH1F* hCastorEGEN;
    TH1F* hEnergyvsEtaCMS;
    TH1F* hBosonPt;
    TH1F* hBosonEta;
    TH1F* hBosonPhi;
    TH1F* hBosonM;
    TH1F* hBosonPtDiff;
    TH1F* hBosonEtaDiff;
    TH1F* hBosonPhiDiff;
    TH1F* hBosonMDiff;
    TH1F* hEntriesCut;

    int nevents;
    int zboson;
    int zbosondiff;
    bool debug;

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
  edm::Service<TFileService> fs;
  //TH1::SetDefaultSumw2(true);


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

  nevents = 0;
  zboson = 0;
  zbosondiff = 0;

}

void SDDYAnalyzer::endJob(){
  hEnergyvsEtaCMS->Scale(1/(float)nevents);	
  hEntriesCut->SetBinContent(2,nevents);
  hEntriesCut->SetBinContent(4,zboson);
  hEntriesCut->SetBinContent(6,zbosondiff);
}

void SDDYAnalyzer::analyze(const edm::Event & ev, const edm::EventSetup&){

  bool debug = false;
  bool lepton1 = false;
  bool lepton2 = false;
  nevents++;

  // Generator Particles
  edm::Handle<reco::GenParticleCollection> genParticles;
  ev.getByLabel(genParticlesTag_, genParticles);

  reco::GenParticleCollection::const_iterator particle1 = genParticles->end();
  reco::GenParticleCollection::const_iterator particle2 = genParticles->end();

  double sumHFMinusGEN = 0.;
  double sumCastorGEN = 0.;
  double sumHFPlusGEN = 0.;

  for(reco::GenParticleCollection::const_iterator genpart = genParticles->begin(); genpart != genParticles->end(); ++genpart){

    if (debug) std::cout << ">>>>>>> pid,status,px,py,px,e= "  << genpart->pdgId() << " , " << genpart->status() << " , " << genpart->px() << " , " << genpart->py() << " , " << genpart->pz() << " , " << genpart->energy() << std::endl;

    if(genpart->status() != 1) continue; // check if genparticle survives.

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
