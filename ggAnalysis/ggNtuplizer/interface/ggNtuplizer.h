#ifndef ggNtuplizer_h
#define ggNtuplizer_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Particle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/CiCPhotonID.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoEgamma/PhotonIdentification/interface/GEDPhoIDTools.h"
#include "TTree.h"
#include "TH1F.h"

#include <memory>
#include <fstream>
#include <map>

using namespace std;
using namespace edm;
using namespace reco;
using namespace pat;

class ggNtuplizer : public edm::EDAnalyzer {

   public:

  explicit ggNtuplizer(const edm::ParameterSet&);
  ~ggNtuplizer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
 private:
  
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  void clearVectors();
  
  void getHandles(const edm::Event & event,
		  edm::Handle<std::vector<reco::GenParticle> > & genParticles,
		  edm::Handle<VertexCollection>                & recVtxs,
		  edm::Handle<VertexCollection>                & recVtxsBS,
		  edm::Handle<double>                          & rhoHandle, 
		  edm::Handle<edm::View<pat::Electron> >       & electronHandle,
		  edm::Handle<edm::View<pat::Photon> >         & photonHandle,
		  edm::Handle<EcalRecHitCollection>            & EBReducedRecHits,
		  edm::Handle<EcalRecHitCollection>            & EEReducedRecHits,
		  edm::Handle<EcalRecHitCollection>            & ESRecHits,
		  edm::Handle<reco::PhotonCollection>          & recoPhotonHandle,
		  edm::Handle<TrackCollection>                 & tracksHandle,
		  edm::Handle<GsfElectronCollection>           & gsfElectronHandle,
		  edm::Handle<PFCandidateCollection>           & pfAllCandidates
		  );

  Bool_t   doGenParticles_;
  Bool_t   runOnParticleGun_;  

  InputTag vtxLabel_;  
  InputTag vtxBSLabel_;
  InputTag rhoLabel_;
  InputTag generatorLabel_;
  InputTag puCollection_;
  InputTag genParticlesCollection_;
  InputTag electronCollection_;
  InputTag photonCollection_;
  InputTag ebReducedRecHitCollection_;
  InputTag eeReducedRecHitCollection_;
  InputTag esReducedRecHitCollection_;
  InputTag recophotonCollection_;
  InputTag tracklabel_;
  InputTag gsfElectronlabel_;
  InputTag pfAllParticles_;

  TTree   *tree_;
  TH1F    *hEvents_; 
  TH1F    *hPU_;
  TH1F    *hPUTrue_;

  EcalClusterLazyTools *lazyTool;
  CiCPhotonID          *cicPhotonId_;

  Int_t          run_;
  Long64_t       event_;
  Int_t          lumis_;
  Bool_t         isData_;
  vector<float>  pdf_;
  Float_t        pthat_;
  Float_t        processID_;
  Int_t          nVtx_;
  Float_t        rho_;
  // genParticle
  Int_t          nMC_;
  vector<int>    mcPID;
  vector<float>  mcVtx_x;
  vector<float>  mcVtx_y;
  vector<float>  mcVtx_z;
  vector<float>  mcPt;
  vector<float>  mcMass;
  vector<float>  mcEta;
  vector<float>  mcPhi;
  vector<float>  mcE;
  vector<float>  mcEt;
  vector<int>    mcGMomPID;
  vector<int>    mcMomPID;
  vector<float>  mcMomPt;
  vector<float>  mcMomMass;
  vector<float>  mcMomEta;
  vector<float>  mcMomPhi;
  vector<int>    mcIndex;
  vector<int>    mcDecayType;
  vector<int>    mcParentage;
  vector<int>    mcStatus;
  // PU
  Int_t          nPUInfo_;  
  vector<int>    nPU_;
  vector<int>    puBX_;
  vector<float>  puTrue_;
  // Electron
  Int_t          nEle_;
  vector<int>    eleCharge_;
  vector<int>    eleChargeConsistent_;
  vector<float>  eleEn_;
  vector<float>  eleSCEn_;
  vector<float>  eleESEn_;
  vector<float>  eleD0_;
  vector<float>  eleDz_;
  vector<float>  elePt_;
  vector<float>  eleEta_;
  vector<float>  elePhi_;
  vector<float>  eleSCEta_;
  vector<float>  eleSCPhi_;
  vector<float>  eleSCRawEn_;
  vector<float>  eleSCEtaWidth_;
  vector<float>  eleSCPhiWidth_;
  vector<float>  eleHoverE_;
  vector<float>  eleEoverP_;
  vector<float>  eleEoverPInv_;
  vector<float>  eleBrem_;
  vector<float>  eledEtaAtVtx_;
  vector<float>  eledPhiAtVtx_;
  vector<float>  eleSigmaIEtaIEta_;
  vector<float>  eleSigmaIEtaIPhi_;
  vector<float>  eleSigmaIPhiIPhi_;
  vector<int>    eleConvVeto_;
  vector<int>    eleMissHits_;
  vector<float>  eleESEffSigmaRR_;
  vector<float>  elePFChIso_;
  vector<float>  elePFPhoIso_;
  vector<float>  elePFNeuIso_;
  // Photon
  Int_t          nPho_;      
  vector<float>  phoE_;
  vector<float>  phoEt_;
  vector<float>  phoEta_;
  vector<float>  phoPhi_;
  vector<float>  phoSCE_;
  vector<float>  phoSCRawE_;
  vector<float>  phoESEn_;
  vector<float>  phoSCEta_;
  vector<float>  phoSCPhi_;
  vector<float>  phoSCEtaWidth_;
  vector<float>  phoSCPhiWidth_;
  vector<float>  phoSCBrem_;
  vector<int>    phohasPixelSeed_;
  vector<int>    phoEleVeto_;
  vector<float>  phoR9_;
  vector<float>  phoHoverE_;
  vector<float>  phoSigmaIEtaIEta_;
  vector<float>  phoSigmaIEtaIPhi_;
  vector<float>  phoSigmaIPhiIPhi_;
  vector<float>  phoE1x3_;
  vector<float>  phoE2x2_;
  vector<float>  phoE2x5Max_;
  vector<float>  phoE5x5_;
  vector<float>  phoESEffSigmaRR_;
  vector<float>  phoPFChIso_;
  vector<float>  phoPFPhoIso_;
  vector<float>  phoPFNeuIso_;
  vector<float>  phoPFChWorstIso_;
  vector<float> pfPhoIsodR03_;
  vector<float> pfChgIsodR03_;
  vector<float> pfNeuIsodR03_;
  vector<float> pfPhoIsoFrix1_;
  vector<float> pfChgIsoFrix1_;
  vector<float> pfNeuIsoFrix1_;
  vector<float> pfPhoIsoFrix2_;
  vector<float> pfChgIsoFrix2_;
  vector<float> pfNeuIsoFrix2_;
  vector<float> pfPhoIsoFrix3_;
  vector<float> pfChgIsoFrix3_;
  vector<float> pfNeuIsoFrix3_;
  vector<float> pfPhoIsoFrix4_;
  vector<float> pfChgIsoFrix4_;
  vector<float> pfNeuIsoFrix4_;
  vector<float> pfPhoIsoFrix5_;
  vector<float> pfChgIsoFrix5_;
  vector<float> pfNeuIsoFrix5_; 
  vector<float> pfPhoIsoFrix6_;
  vector<float> pfChgIsoFrix6_;
  vector<float> pfNeuIsoFrix6_;
  vector<float> pfPhoIsoFrix7_;
  vector<float> pfChgIsoFrix7_;
  vector<float> pfNeuIsoFrix7_;
  vector<float> pfPhoIsoFrix8_;
  vector<float> pfChgIsoFrix8_;
  vector<float> pfNeuIsoFrix8_; 
  vector<int> phoLoose_;
  vector<int> phoMedium_;
  vector<int> phoTight_;
 // Physics objects handles
  Handle<std::vector<reco::GenParticle> > genParticlesHandle_;
  Handle<VertexCollection>                recVtxs_;
  Handle<VertexCollection>                recVtxsBS_;
  Handle<double>                          rhoHandle_;
  Handle<View<pat::Electron> >            electronHandle_;
  Handle<View<pat::Photon> >              photonHandle_;
  Handle<EcalRecHitCollection>            EBReducedRecHits_;
  Handle<EcalRecHitCollection>            EEReducedRecHits_;
  Handle<EcalRecHitCollection>            ESReducedRecHits_;
  Handle<reco::PhotonCollection>          recoPhotonHandle_;
  Handle<TrackCollection>                 tracksHandle_;
  Handle<GsfElectronCollection>           gsfElectronHandle_;
  Handle<PFCandidateCollection>           pfAllCandidates_; 

};

#endif
