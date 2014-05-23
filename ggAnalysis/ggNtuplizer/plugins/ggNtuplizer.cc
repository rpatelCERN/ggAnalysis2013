#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"
#include "ggAnalysis/ggNtuplizer/interface/GenParticleParentage.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

ggNtuplizer::ggNtuplizer(const edm::ParameterSet& ps) {

  doGenParticles_            = ps.getParameter<bool>("doGenParticles");
  runOnParticleGun_          = ps.getParameter<bool>("runOnParticleGun");

  vtxLabel_                  = ps.getParameter<InputTag>("VtxLabel");
  vtxBSLabel_                = ps.getParameter<InputTag>("VtxBSLabel");
  rhoLabel_                  = ps.getParameter<InputTag>("rhoLabel");
  generatorLabel_            = ps.getParameter<InputTag>("generatorLabel");
  puCollection_              = ps.getParameter<InputTag>("pileupCollection");
  genParticlesCollection_    = ps.getParameter<InputTag>("genParticleSrc");
  electronCollection_        = ps.getParameter<InputTag>("electronSrc");
  photonCollection_          = ps.getParameter<InputTag>("photonSrc");
  ebReducedRecHitCollection_ = ps.getParameter<InputTag>("ebReducedRecHitCollection");
  eeReducedRecHitCollection_ = ps.getParameter<InputTag>("eeReducedRecHitCollection");
  esReducedRecHitCollection_ = ps.getParameter<InputTag>("esReducedRecHitCollection");
  recophotonCollection_      = ps.getParameter<InputTag>("recoPhotonSrc");
  tracklabel_                = ps.getParameter<InputTag>("TrackLabel");
  gsfElectronlabel_          = ps.getParameter<InputTag>("gsfElectronLabel");
  pfAllParticles_            = ps.getParameter<InputTag>("PFAllCandidates");

  cicPhotonId_ = new CiCPhotonID(ps);

  Service<TFileService> fs;
  hEvents_ = fs->make<TH1F>("hEvents",    "total processed and skimmed events",   2,  0,   2);
  hPU_     = fs->make<TH1F>("hPU",        "number of pileup",                   200,  0, 200);
  hPUTrue_ = fs->make<TH1F>("hPUTrue",    "number of true pilepu"             , 1000, 0, 200);
  tree_    = fs->make<TTree>("EventTree", "Event data");

  tree_->Branch("run", &run_, "run/I");
  tree_->Branch("event", &event_, "event/L");
  tree_->Branch("lumis", &lumis_, "lumis/I");
  tree_->Branch("isData", &isData_, "isData/O");
  tree_->Branch("nVtx", &nVtx_, "nVtx/I");
  tree_->Branch("rho", &rho_, "rho/F");
  if (doGenParticles_) {
    tree_->Branch("pdf", &pdf_);
    tree_->Branch("pthat", &pthat_, "pthat/F");
    tree_->Branch("processID", &processID_, "processID/F");
    // genParticle
    tree_->Branch("nMC", &nMC_, "nMC/I");
    tree_->Branch("mcPID", &mcPID);
    tree_->Branch("mcVtx_x", &mcVtx_x);
    tree_->Branch("mcVtx_y", &mcVtx_y);
    tree_->Branch("mcVtx_z", &mcVtx_z);
    tree_->Branch("mcPt", &mcPt);
    tree_->Branch("mcMass", &mcMass);
    tree_->Branch("mcEta", &mcEta);
    tree_->Branch("mcPhi", &mcPhi);
    tree_->Branch("mcE", &mcE);
    tree_->Branch("mcEt", &mcEt);
    tree_->Branch("mcGMomPID", &mcGMomPID);
    tree_->Branch("mcMomPID", &mcMomPID);
    tree_->Branch("mcMomPt", &mcMomPt);
    tree_->Branch("mcMomMass", &mcMomMass);
    tree_->Branch("mcMomEta", &mcMomEta);
    tree_->Branch("mcMomPhi", &mcMomPhi);
    tree_->Branch("mcIndex", &mcIndex);
    tree_->Branch("mcDecayType", &mcDecayType); //-999:non W or Z, 1:hardronic, 2:e, 3:mu, 4:tau
    tree_->Branch("mcParentage", &mcParentage); // 16*lepton + 8*boson + 4*non-prompt + 2*qcd + exotics
    tree_->Branch("mcStatus", &mcStatus); // status of the particle
    // Gen MET
    //tree_->Branch("genMET", &genMET_, "genMET/F");
    //tree_->Branch("genMETPhi", &genMETPhi_, "genMETPhi/F");
    // PU Info
    tree_->Branch("nPUInfo", &nPUInfo_, "nPUInfo/I");
    tree_->Branch("nPU", &nPU_);
    tree_->Branch("puBX", &puBX_);
    tree_->Branch("puTrue", &puTrue_);
  }
  // Electron
  tree_->Branch("nEle",                &nEle_, "nEle/I");
  tree_->Branch("eleCharge",           &eleCharge_);
  tree_->Branch("eleChargeConsistent", &eleChargeConsistent_);
  tree_->Branch("eleEn",               &eleEn_);
  tree_->Branch("eleSCEn",             &eleSCEn_);
  tree_->Branch("eleESEn",             &eleESEn_);
  tree_->Branch("eleD0",               &eleD0_);
  tree_->Branch("eleDz",               &eleDz_);
  tree_->Branch("elePt",               &elePt_);
  tree_->Branch("eleEta",              &eleEta_);
  tree_->Branch("elePhi",              &elePhi_);
  tree_->Branch("eleSCEta",            &eleSCEta_);
  tree_->Branch("eleSCPhi",            &eleSCPhi_);
  tree_->Branch("eleSCRawEn",          &eleSCRawEn_);
  tree_->Branch("eleSCEtaWidth",       &eleSCEtaWidth_);
  tree_->Branch("eleSCPhiWidth",       &eleSCPhiWidth_);
  tree_->Branch("eleHoverE",           &eleHoverE_);
  tree_->Branch("eleEoverP",           &eleEoverP_);
  tree_->Branch("eleEoverPInv",        &eleEoverPInv_);
  tree_->Branch("eleBrem",             &eleBrem_);
  tree_->Branch("eledEtaAtVtx",        &eledEtaAtVtx_);
  tree_->Branch("eledPhiAtVtx",        &eledPhiAtVtx_);
  tree_->Branch("eleSigmaIEtaIEta",    &eleSigmaIEtaIEta_);
  tree_->Branch("eleSigmaIEtaIPhi",    &eleSigmaIEtaIPhi_);
  tree_->Branch("eleSigmaIPhiIPhi",    &eleSigmaIPhiIPhi_);
  tree_->Branch("eleConvVeto",         &eleConvVeto_);
  tree_->Branch("eleMissHits",         &eleMissHits_);
  tree_->Branch("eleESEffSigmaRR",     &eleESEffSigmaRR_);
  tree_->Branch("elePFChIso",          &elePFChIso_);
  tree_->Branch("elePFPhoIso",         &elePFPhoIso_);
  tree_->Branch("elePFNeuIso",         &elePFNeuIso_);
  // Photon
  tree_->Branch("nPho",             &nPho_, "nPho/I");
  tree_->Branch("phoE",             &phoE_);
  tree_->Branch("phoEt",            &phoEt_);
  tree_->Branch("phoEta",           &phoEta_);
  tree_->Branch("phoPhi",           &phoPhi_);
  tree_->Branch("phoSCE",           &phoSCE_);
  tree_->Branch("phoSCRawE",        &phoSCRawE_);
  tree_->Branch("phoESEn",          &phoESEn_);
  tree_->Branch("phoSCEta",         &phoSCEta_);
  tree_->Branch("phoSCPhi",         &phoSCPhi_);
  tree_->Branch("phoSCEtaWidth",    &phoSCEtaWidth_);
  tree_->Branch("phoSCPhiWidth",    &phoSCPhiWidth_);
  tree_->Branch("phoSCBrem",        &phoSCBrem_);
  tree_->Branch("phohasPixelSeed",  &phohasPixelSeed_);
  tree_->Branch("phoEleVeto",       &phoEleVeto_);
  tree_->Branch("phoR9",            &phoR9_);
  tree_->Branch("phoHoverE",        &phoHoverE_);
  tree_->Branch("phoSigmaIEtaIEta", &phoSigmaIEtaIEta_);
  tree_->Branch("phoSigmaIEtaIPhi", &phoSigmaIEtaIPhi_);
  tree_->Branch("phoSigmaIPhiIPhi", &phoSigmaIPhiIPhi_);
  tree_->Branch("phoE1x3",          &phoE1x3_);
  tree_->Branch("phoE2x2",          &phoE2x2_);
  tree_->Branch("phoE2x5Max",       &phoE2x5Max_);
  tree_->Branch("phoE5x5",          &phoE5x5_);
  tree_->Branch("phoESEffSigmaRR",  &phoESEffSigmaRR_);
  tree_->Branch("phoPFChIso",       &phoPFChIso_);
  tree_->Branch("phoPFPhoIso",      &phoPFPhoIso_);
  tree_->Branch("phoPFNeuIso",      &phoPFNeuIso_);
  tree_->Branch("phoPFChWorstIso",  &phoPFChWorstIso_);
  tree_->Branch("pfPhoIsodR03", &pfPhoIsodR03_);
  tree_->Branch("pfChgIsodR03", &pfChgIsodR03_); 
  tree_->Branch("pfNeuIsodR03", &pfNeuIsodR03_);

  tree_->Branch("pfPhoIsoFrix1", &pfPhoIsoFrix1_);
  tree_->Branch("pfPhoIsoFrix2", &pfPhoIsoFrix2_);
  tree_->Branch("pfPhoIsoFrix3", &pfPhoIsoFrix3_);
  tree_->Branch("pfPhoIsoFrix4", &pfPhoIsoFrix4_);
  tree_->Branch("pfPhoIsoFrix5", &pfPhoIsoFrix5_);
  tree_->Branch("pfPhoIsoFrix6", &pfPhoIsoFrix6_);
  tree_->Branch("pfPhoIsoFrix7", &pfPhoIsoFrix7_);
  tree_->Branch("pfPhoIsoFrix8", &pfPhoIsoFrix8_);

  tree_->Branch("pfChgIsoFrix1", &pfChgIsoFrix1_);
  tree_->Branch("pfChgIsoFrix2", &pfChgIsoFrix2_);
  tree_->Branch("pfChgIsoFrix3", &pfChgIsoFrix3_);
  tree_->Branch("pfChgIsoFrix4", &pfChgIsoFrix4_);
  tree_->Branch("pfChgIsoFrix5", &pfChgIsoFrix5_);
  tree_->Branch("pfChgIsoFrix6", &pfChgIsoFrix6_);
  tree_->Branch("pfChgIsoFrix7", &pfChgIsoFrix7_);
  tree_->Branch("pfChgIsoFrix8", &pfChgIsoFrix8_);

  tree_->Branch("pfNeuIsoFrix1", &pfNeuIsoFrix1_);
  tree_->Branch("pfNeuIsoFrix2", &pfNeuIsoFrix2_);
  tree_->Branch("pfNeuIsoFrix3", &pfNeuIsoFrix3_);
  tree_->Branch("pfNeuIsoFrix4", &pfNeuIsoFrix4_);
  tree_->Branch("pfNeuIsoFrix5", &pfNeuIsoFrix5_);
  tree_->Branch("pfNeuIsoFrix6", &pfNeuIsoFrix6_);
  tree_->Branch("pfNeuIsoFrix7", &pfNeuIsoFrix7_);
  tree_->Branch("pfNeuIsoFrix8", &pfNeuIsoFrix8_);

  tree_->Branch("phoLoose",&phoLoose_);
  tree_->Branch("phoMedium",&phoMedium_);
  tree_->Branch("phoTight",&phoTight_);
}

ggNtuplizer::~ggNtuplizer() {
  delete cicPhotonId_;
}

void ggNtuplizer::getHandles(const edm::Event & event,
			     edm::Handle<std::vector<reco::GenParticle> > & genParticlesHandle,
			     edm::Handle<VertexCollection>                & recVtxs,
			     edm::Handle<VertexCollection>                & recVtxsBS,
			     edm::Handle<double>                          & rhoHandle,
			     edm::Handle<edm::View<pat::Electron> >       & electronHandle,
			     edm::Handle<edm::View<pat::Photon> >         & photonHandle,
			     edm::Handle<EcalRecHitCollection>            & EBReducedRecHits,
			     edm::Handle<EcalRecHitCollection>            & EEReducedRecHits,
			     edm::Handle<EcalRecHitCollection>            & ESReducedRecHits,
			     edm::Handle<reco::PhotonCollection>          & recoPhotonHandle,
			     edm::Handle<TrackCollection>                 & tracksHandle,
			     edm::Handle<GsfElectronCollection>           & gsfElectronHandle,
			     edm::Handle<PFCandidateCollection>           & pfAllCandidates) { 

  if (doGenParticles_) event.getByLabel(genParticlesCollection_, genParticlesHandle);
  event.getByLabel(vtxLabel_,                  recVtxs);
  event.getByLabel(vtxBSLabel_,                recVtxsBS);
  event.getByLabel(rhoLabel_,                  rhoHandle);
  event.getByLabel(electronCollection_,        electronHandle);
  event.getByLabel(photonCollection_,          photonHandle);
  event.getByLabel(ebReducedRecHitCollection_, EBReducedRecHits);
  event.getByLabel(eeReducedRecHitCollection_, EEReducedRecHits);
  event.getByLabel(esReducedRecHitCollection_, ESReducedRecHits);
  event.getByLabel(recophotonCollection_,      recoPhotonHandle);
  event.getByLabel(tracklabel_,                tracksHandle);
  event.getByLabel(gsfElectronlabel_,          gsfElectronHandle);
  event.getByLabel(pfAllParticles_,            pfAllCandidates);
}

void ggNtuplizer::analyze(const edm::Event& e, const edm::EventSetup& es) {

  this->getHandles(e, genParticlesHandle_, recVtxs_, recVtxsBS_, rhoHandle_, electronHandle_, 
		   photonHandle_, EBReducedRecHits_, EEReducedRecHits_, ESReducedRecHits_, recoPhotonHandle_,
		   tracksHandle_, gsfElectronHandle_, pfAllCandidates_);
  
  clearVectors();
  hEvents_->Fill(0.5);
  GEDPhoIDTools GEDIdTool(e);
  //GEDIdTool.~GEDPhoIDTools();
  lazyTool = new EcalClusterLazyTools(e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);

  run_    = e.id().run();
  event_  = e.id().event();
  lumis_  = e.luminosityBlock();
  isData_ = e.isRealData();
  rho_    = *(rhoHandle_.product());

  cicPhotonId_->configure(recVtxsBS_, tracksHandle_, gsfElectronHandle_, pfAllCandidates_, rho_); 

  // vertex
  nVtx_ = 0;
  math::XYZPoint pv(0, 0, 0);
  if (recVtxs_.isValid()) {
    for (size_t i=0; i<recVtxs_->size(); ++i) {
      if (!((*recVtxs_)[i].isFake())) {
	if (nVtx_ == 0) pv.SetXYZ((*recVtxs_)[i].x(), (*recVtxs_)[i].y(), (*recVtxs_)[i].z());
	nVtx_++;
      }
    }
  }

  // PDF information, pthat and processID
  nPUInfo_ = 0; 
  if (!isData_) {
    
    pthat_ = -99;
    processID_ = -99;

    Handle<GenEventInfoProduct> pdfInfoHandle;
    if (e.getByLabel(generatorLabel_, pdfInfoHandle)) {
      if (pdfInfoHandle->pdf()) {
	pdf_.push_back(pdfInfoHandle->pdf()->id.first);    // PDG ID of incoming parton #1
	pdf_.push_back(pdfInfoHandle->pdf()->id.second);   // PDG ID of incoming parton #2
	pdf_.push_back(pdfInfoHandle->pdf()->x.first);     // x value of parton #1
	pdf_.push_back(pdfInfoHandle->pdf()->x.second);    // x value of parton #2
	pdf_.push_back(pdfInfoHandle->pdf()->xPDF.first);  // PDF weight for parton #1
	pdf_.push_back(pdfInfoHandle->pdf()->xPDF.second); // PDF weight for parton #2
	pdf_.push_back(pdfInfoHandle->pdf()->scalePDF);    // scale of the hard interaction
      }
    }

    Handle<GenEventInfoProduct> genEventScale;
    if (e.getByLabel(generatorLabel_, genEventScale)) {
      if (genEventScale->hasBinningValues()) pthat_ = genEventScale->binningValues()[0];
      processID_ = genEventScale->signalProcessID();
    }

    Handle< vector<PileupSummaryInfo> > pileupHandle;
    if (e.getByLabel(puCollection_, pileupHandle)) {
      vector<PileupSummaryInfo>::const_iterator PVI;
      for (PVI = pileupHandle->begin(); PVI != pileupHandle->end(); ++PVI) {
	if (PVI->getBunchCrossing() == 0) {
	  hPU_->Fill(PVI->getPU_NumInteractions());
	  hPUTrue_->Fill(PVI->getTrueNumInteractions());
	}
	nPU_.push_back(PVI->getPU_NumInteractions());
	puTrue_.push_back(PVI->getTrueNumInteractions());
	puBX_.push_back(PVI->getBunchCrossing());
	nPUInfo_++;
      }
    }
  }
  
  // GenParticle
  if (!isData_ && genParticlesHandle_.isValid() ) {
    
    nMC_ = 0;
    int genIndex = 0;
    
    for (vector<GenParticle>::const_iterator ip = genParticlesHandle_->begin(); ip != genParticlesHandle_->end(); ++ip) {
      genIndex++;
      
      int status = ip->status() - 10*(ip->status()/10);
      bool stableFinalStateParticle = status == 1 && ip->pt() > 5.0;
      
      // keep all the photons with pT > 5.0 and all leptons;
      bool photonOrLepton = 
	(status == 1 && ip->pdgId() == 22 && ip->pt() > 5.0 ) ||
	(status == 1 && ( abs(ip->pdgId()) >= 11 && abs(ip->pdgId()) <= 16 ))  ||
	(status < 10 && abs(ip->pdgId()) == 15 );
      // select also Z, W, H, and top
      bool heavyParticle =  
	(ip->pdgId() == 23 || abs(ip->pdgId()) == 24 || ip->pdgId() == 25 || 
	 abs(ip->pdgId()) == 6 || abs(ip->pdgId()) == 5);
      
      if ( stableFinalStateParticle || heavyParticle || photonOrLepton ) {
	const Candidate *p = (const Candidate*)&(*ip);
	if (!runOnParticleGun_ && !p->mother()) continue;
	
	reco::GenParticleRef partRef = reco::GenParticleRef(genParticlesHandle_,
							    ip-genParticlesHandle_->begin());
	genpartparentage::GenParticleParentage particleHistory(partRef);

	mcPID    .push_back(p->pdgId());
	mcVtx_x  .push_back(p->vx());
	mcVtx_y  .push_back(p->vy());
	mcVtx_z  .push_back(p->vz());
	mcPt     .push_back(p->pt());
	mcMass   .push_back(p->mass());
	mcEta    .push_back(p->eta());
	mcPhi    .push_back(p->phi());
	mcE      .push_back(p->energy());
	mcEt     .push_back(p->et());
	mcParentage.push_back(particleHistory.hasLeptonParent()*16   + 
			      particleHistory.hasBosonParent()*8     + 
			      particleHistory.hasNonPromptParent()*4 +
			      particleHistory.hasQCDParent()*2       +
			      particleHistory.hasExoticParent());
	mcStatus.push_back(p->status());

	int mcDecayType_ = -999;
	// if genParticle is W or Z, check its decay type
	if ( ip->pdgId() == 23 || abs(ip->pdgId()) == 24 ) {
	  for (size_t k=0; k < p->numberOfDaughters(); ++k) {
	    const Candidate *dp = p->daughter(k);
	    if (abs(dp->pdgId())<=6)
	      mcDecayType_ = 1;
	    else if (abs(dp->pdgId())==11 || abs(dp->pdgId())==12)
	      mcDecayType_ = 2;
	    else if (abs(dp->pdgId())==13 || abs(dp->pdgId())==14)
	      mcDecayType_ = 3;
	    else if (abs(dp->pdgId())==15 || abs(dp->pdgId())==16)
	      mcDecayType_ = 4;
	  }
	}
	mcDecayType.push_back(mcDecayType_);
	int mcGMomPID_ = -999;
	int mcMomPID_  = -999;
	float mcMomPt_    = -999.;
	float mcMomMass_  = -999.;
	float mcMomEta_   = -999.;
	float mcMomPhi_   = -999.;
	if ( particleHistory.hasRealParent() ) {
	  reco::GenParticleRef momRef = particleHistory.parent();
	  if ( momRef.isNonnull() && momRef.isAvailable() ) {
	    mcMomPID_  = momRef->pdgId();
	    mcMomPt_   = momRef->pt();
	    mcMomMass_ = momRef->mass();
	    mcMomEta_  = momRef->eta();
	    mcMomPhi_  = momRef->phi();
	    
	    // get Granny
	    genpartparentage::GenParticleParentage motherParticle(momRef);
	    if ( motherParticle.hasRealParent() ) {
	      reco::GenParticleRef granny = motherParticle.parent();
	      mcGMomPID_ = granny->pdgId();
	    }
	  }
	}
	mcGMomPID.push_back(mcGMomPID_);
	mcMomPID.push_back(mcMomPID_);
	mcMomPt.push_back(mcMomPt_);
	mcMomMass.push_back(mcMomMass_);
	mcMomEta.push_back(mcMomEta_);
	mcMomPhi.push_back(mcMomPhi_);
	
	mcIndex.push_back(genIndex-1);
	nMC_++;
      } // save info on particles of interest
    } // loop over gen-level particles
  }
  
  // electrons
  nEle_ = 0;
  if (electronHandle_.isValid()) {
    for (View<pat::Electron>::const_iterator iEle = electronHandle_->begin(); iEle != electronHandle_->end(); ++iEle) {

      eleCharge_          .push_back(iEle->charge());
      eleChargeConsistent_.push_back((Int_t)iEle->isGsfCtfScPixChargeConsistent());
      eleEn_              .push_back(iEle->energy());
      eleD0_              .push_back(iEle->gsfTrack()->dxy(pv));
      eleDz_              .push_back(iEle->gsfTrack()->dz(pv));
      elePt_              .push_back(iEle->pt());
      eleEta_             .push_back(iEle->eta());
      elePhi_             .push_back(iEle->phi());
      eleSCEn_            .push_back(iEle->superCluster()->energy());
      eleESEn_            .push_back(iEle->superCluster()->preshowerEnergy());
      eleSCEta_           .push_back(iEle->superCluster()->eta());
      eleSCPhi_           .push_back(iEle->superCluster()->phi());
      eleSCRawEn_         .push_back(iEle->superCluster()->rawEnergy());
      eleSCEtaWidth_      .push_back(iEle->superCluster()->etaWidth());
      eleSCPhiWidth_      .push_back(iEle->superCluster()->phiWidth());
      eleHoverE_          .push_back(iEle->hcalOverEcalBc());
      eleEoverP_          .push_back(iEle->eSuperClusterOverP());
      eleEoverPInv_       .push_back(fabs(1./iEle->ecalEnergy()-1./iEle->trackMomentumAtVtx().R()));
      eleBrem_            .push_back(iEle->fbrem());
      eledEtaAtVtx_       .push_back(iEle->deltaEtaSuperClusterTrackAtVtx());
      eledPhiAtVtx_       .push_back(iEle->deltaPhiSuperClusterTrackAtVtx());
      eleSigmaIEtaIEta_   .push_back(iEle->sigmaIetaIeta());
      eleSigmaIEtaIPhi_   .push_back(iEle->sigmaIetaIphi());
      eleSigmaIPhiIPhi_   .push_back(iEle->sigmaIphiIphi());
      eleConvVeto_        .push_back((Int_t)iEle->passConversionVeto()); // ConvVtxFit || missHit == 0
      eleMissHits_        .push_back(iEle->gsfTrack()->trackerExpectedHitsInner().numberOfHits());
      eleESEffSigmaRR_    .push_back(lazyTool->eseffsirir(*((*iEle).superCluster())));
      elePFChIso_         .push_back(iEle->chargedHadronIso());
      elePFPhoIso_        .push_back(iEle->photonIso());
      elePFNeuIso_        .push_back(iEle->neutralHadronIso());

      nEle_++;
    }
  } 

  // photons
  nPho_ = 0;
  if (photonHandle_.isValid()) {
    for (View<pat::Photon>::const_iterator iPho = photonHandle_->begin(); iPho != photonHandle_->end(); ++iPho) {

      phoE_             .push_back(iPho->energy());
      phoEt_            .push_back(iPho->et());
      phoEta_           .push_back(iPho->eta());
      phoPhi_           .push_back(iPho->phi());
      phoSCE_           .push_back((*iPho).superCluster()->energy());
      phoSCRawE_        .push_back((*iPho).superCluster()->rawEnergy());
      phoESEn_          .push_back((*iPho).superCluster()->preshowerEnergy());
      phoSCEta_         .push_back((*iPho).superCluster()->eta());
      phoSCPhi_         .push_back((*iPho).superCluster()->phi());
      phoSCEtaWidth_    .push_back((*iPho).superCluster()->etaWidth());
      phoSCPhiWidth_    .push_back((*iPho).superCluster()->phiWidth());
      phoSCBrem_        .push_back((*iPho).superCluster()->phiWidth()/(*iPho).superCluster()->etaWidth());
      phohasPixelSeed_  .push_back((Int_t)iPho->hasPixelSeed());
      phoEleVeto_       .push_back((Int_t)iPho->passElectronVeto());
      phoR9_            .push_back(iPho->r9());
      phoHoverE_        .push_back(iPho->hadTowOverEm());
      phoSigmaIEtaIEta_ .push_back(iPho->see());
      phoSigmaIEtaIPhi_ .push_back(iPho->sep());
      phoSigmaIPhiIPhi_ .push_back(iPho->spp());
      phoE1x3_          .push_back(lazyTool->e1x3(*((*iPho).superCluster()->seed())));
      phoE2x2_          .push_back(lazyTool->e2x2(*((*iPho).superCluster()->seed())));
      phoE2x5Max_       .push_back(lazyTool->e2x5Max(*((*iPho).superCluster()->seed())));
      phoE5x5_          .push_back(lazyTool->e5x5(*((*iPho).superCluster()->seed())));
      phoESEffSigmaRR_  .push_back(lazyTool->eseffsirir(*((*iPho).superCluster())));
      phoPFChIso_       .push_back(iPho->chargedHadronIso());
      phoPFPhoIso_      .push_back(iPho->photonIso());
      phoPFNeuIso_      .push_back(iPho->neutralHadronIso());

      size_t rightRecoPho = -1;
      for (size_t iv = 0; iv < recoPhotonHandle_->size(); ++iv) {
	reco::PhotonRef recophoRef2(recoPhotonHandle_, iv);
	if (deltaR(iPho->eta(), iPho->phi(), recophoRef2->eta(), recophoRef2->phi()) < 0.01) rightRecoPho = iv;
	 reco::PhotonRef aPho(reco::PhotonRef(recoPhotonHandle_, iv));
	//GEDIdTool->setPhotonP4();
      }
      //reco::PhotonRef photonRef (reco::PhotonRef( recoPhotonHandle_,rightRecoPho));
      reco::PhotonRef recophoRef(recoPhotonHandle_, rightRecoPho);
      reco::Vertex pv = recVtxs_->at(0); 
      GEDIdTool.setPhotonP4(recophoRef, pv);
      float pfPhoIsodR03=GEDIdTool.SolidConeIso(0.3, reco::PFCandidate::gamma); 
      float pfChgIsodR03=GEDIdTool.SolidConeIso(0.3, reco::PFCandidate::h);
      float pfNeuIsodR03=GEDIdTool.SolidConeIso(0.3, reco::PFCandidate::h0);
 //     cout<<"pho Iso "<<pfPhoIsodR03<<" chg Iso "<<pfChgIsodR03<<" neu "<<pfNeuIsodR03<<endl;
      pfPhoIsodR03_.push_back(pfPhoIsodR03);
      pfChgIsodR03_.push_back(pfChgIsodR03);	
      pfNeuIsodR03_.push_back(pfNeuIsodR03);
      std::vector<double>IsoRings;
      GEDIdTool.FrixioneIso(0.1, 8, reco::PFCandidate::gamma,IsoRings);
      pfPhoIsoFrix1_.push_back(IsoRings[0]);
      pfPhoIsoFrix2_.push_back(IsoRings[1]);
      pfPhoIsoFrix3_.push_back(IsoRings[2]);
      pfPhoIsoFrix4_.push_back(IsoRings[3]);
      pfPhoIsoFrix5_.push_back(IsoRings[4]);
      pfPhoIsoFrix6_.push_back(IsoRings[5]);
      pfPhoIsoFrix7_.push_back(IsoRings[6]);
      pfPhoIsoFrix8_.push_back(IsoRings[7]);
      IsoRings.resize(0);
      GEDIdTool.FrixioneIso(0.1, 8, reco::PFCandidate::h,IsoRings);	
      pfChgIsoFrix1_.push_back(IsoRings[0]);
      pfChgIsoFrix2_.push_back(IsoRings[1]);
      pfChgIsoFrix3_.push_back(IsoRings[2]);
      pfChgIsoFrix4_.push_back(IsoRings[3]);
      pfChgIsoFrix5_.push_back(IsoRings[4]);
      pfChgIsoFrix6_.push_back(IsoRings[5]);
      pfChgIsoFrix7_.push_back(IsoRings[6]);
      pfChgIsoFrix8_.push_back(IsoRings[7]);
      IsoRings.resize(0);

      GEDIdTool.FrixioneIso(0.1, 8, reco::PFCandidate::h0,IsoRings);    
      pfNeuIsoFrix1_.push_back(IsoRings[0]);
      pfNeuIsoFrix2_.push_back(IsoRings[1]);
      pfNeuIsoFrix3_.push_back(IsoRings[2]);
      pfNeuIsoFrix4_.push_back(IsoRings[3]);
      pfNeuIsoFrix5_.push_back(IsoRings[4]);
      pfNeuIsoFrix6_.push_back(IsoRings[5]);
      pfNeuIsoFrix7_.push_back(IsoRings[6]);
      pfNeuIsoFrix8_.push_back(IsoRings[7]);
//      GEDIdTool.defaultCuts();
      int passloose=(int)GEDIdTool.CutBasedID(GEDPhoIDTools::Loose,true);
      int passMed=(int)GEDIdTool.CutBasedID(GEDPhoIDTools::Medium,true);
      int passTight=(int)GEDIdTool.CutBasedID(GEDPhoIDTools::Tight,true);
      phoLoose_.push_back(passloose);
      phoMedium_.push_back(passMed);
      phoTight_.push_back(passTight);
      //if(pass)cout<<"passes "<<endl;
//     cout<<"Iso Rings "<<IsoRings.size()<<endl;
      //std::vector<float> vtxIsolations03 = cicPhotonId_->pfTkIsoWithVertex(recophoRef, 0.3, 0.02, 0.02, 0.0, 0.2, 0.1, reco::PFCandidate::h);
//std::vector<float> vtxIsolations03 = pfTkIsoWithVertex(*recophoRef, 0.3, 0.02, 0.02, 0.0, 0.2, 0.1, reco::PFCandidate::h);
   std::vector<float> vtxIsolations03 = cicPhotonId_->pfTkIsoWithVertex(recophoRef, 0.3, 0.0, 0.0, 0.0, 0.2, 0.1, reco::PFCandidate::h);
 //std::vector<float> vtxIsolations03 = cicPhotonId_->pfTkIsoWithVertex(recophoRef, 0.3, 0.0, 0.0, 0.0, 0.2, 0.1, reco::PFCandidate::h);
     //std::vector< std::vector<float> > vtxIsolations03 = cicPhotonId_->pfTkIsoWithVertex(recophoRef, 0.3, 0.0, 0.0, 0.0, 0.2, 0.1, reco::PFCandidate::h);
      //cout<<iPho->chargedHadronIso()<<" "<<iPho->chargedHadronIsoWrongVtx()<<" "<<vtxIsolations03[0]<<" "<<*max_element(vtxIsolations03.begin(), vtxIsolations03.end())<<endl;
      phoPFChWorstIso_.push_back(*max_element(vtxIsolations03.begin(), vtxIsolations03.end()));

      nPho_++;
    }
  }
  
  hEvents_->Fill(1.5);
  tree_->Fill();

  delete lazyTool;
}

void ggNtuplizer::beginJob() {
}

void ggNtuplizer::endJob() {
}

void ggNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void ggNtuplizer::clearVectors() {

  mcPID.clear();
  mcVtx_x.clear();
  mcVtx_y.clear();
  mcVtx_z.clear();
  mcPt.clear();
  mcMass.clear();
  mcEta.clear();
  mcPhi.clear();
  mcE.clear();
  mcEt.clear();
  mcGMomPID.clear();
  mcMomPID.clear();
  mcMomPt.clear();
  mcMomMass.clear();
  mcMomEta.clear();
  mcMomPhi.clear();
  mcIndex.clear();
  mcDecayType.clear();
  mcParentage.clear();
  mcStatus.clear();

  pdf_.clear();
  nPU_.clear();
  puBX_.clear();
  puTrue_.clear();

  eleCharge_.clear();
  eleChargeConsistent_.clear();
  eleEn_.clear();
  eleSCEn_.clear();
  eleESEn_.clear();
  eleD0_.clear();
  eleDz_.clear();
  elePt_.clear();
  eleEta_.clear();
  elePhi_.clear();
  eleSCEta_.clear();
  eleSCPhi_.clear();
  eleSCRawEn_.clear();
  eleSCEtaWidth_.clear();
  eleSCPhiWidth_.clear();
  eleHoverE_.clear();
  eleEoverP_.clear();
  eleEoverPInv_.clear();
  eleBrem_.clear();
  eledEtaAtVtx_.clear();
  eledPhiAtVtx_.clear();
  eleSigmaIEtaIEta_.clear();
  eleSigmaIEtaIPhi_.clear();
  eleSigmaIPhiIPhi_.clear();
  eleConvVeto_.clear();
  eleMissHits_.clear();
  eleESEffSigmaRR_.clear();
  elePFChIso_.clear();
  elePFPhoIso_.clear();
  elePFNeuIso_.clear();

  phoE_.clear();
  phoEt_.clear();
  phoEta_.clear();
  phoPhi_.clear();
  phoSCE_.clear();
  phoSCRawE_.clear();
  phoESEn_.clear();
  phoSCEta_.clear();
  phoSCPhi_.clear();
  phoSCEtaWidth_.clear();
  phoSCPhiWidth_.clear();
  phoSCBrem_.clear();
  phohasPixelSeed_.clear();
  phoEleVeto_.clear();
  phoR9_.clear();
  phoHoverE_.clear();
  phoSigmaIEtaIEta_.clear();
  phoSigmaIEtaIPhi_.clear();
  phoSigmaIPhiIPhi_.clear();
  phoE1x3_.clear();
  phoE2x2_.clear();
  phoE2x5Max_.clear();
  phoE5x5_.clear();
  phoESEffSigmaRR_.clear();
  phoPFChIso_.clear();
  phoPFPhoIso_.clear();
  phoPFNeuIso_.clear();
  phoPFChWorstIso_.clear();
  pfPhoIsodR03_.clear();
  pfChgIsodR03_.clear();
  pfNeuIsodR03_.clear();

  pfPhoIsoFrix1_.clear();
  pfPhoIsoFrix2_.clear();
  pfPhoIsoFrix3_.clear();
  pfPhoIsoFrix4_.clear();
  pfPhoIsoFrix5_.clear();
  pfPhoIsoFrix6_.clear();
  pfPhoIsoFrix7_.clear();
  pfPhoIsoFrix8_.clear();

  pfChgIsoFrix1_.clear();
  pfChgIsoFrix2_.clear();
  pfChgIsoFrix3_.clear();
  pfChgIsoFrix4_.clear();
  pfChgIsoFrix5_.clear();
  pfChgIsoFrix6_.clear();
  pfChgIsoFrix7_.clear();
  pfChgIsoFrix8_.clear();

  pfNeuIsoFrix1_.clear();
  pfNeuIsoFrix2_.clear();
  pfNeuIsoFrix3_.clear();
  pfNeuIsoFrix4_.clear();
  pfNeuIsoFrix5_.clear();
  pfNeuIsoFrix6_.clear();
  pfNeuIsoFrix7_.clear();
  pfNeuIsoFrix8_.clear();
 
       phoLoose_.clear();
      phoMedium_.clear();
      phoTight_.clear();

}

DEFINE_FWK_MODULE(ggNtuplizer);
