import FWCore.ParameterSet.Config as cms
from HiggsAnalysis.HiggsTo2photons.hggPhotonIDCuts_cfi import *

ggNtuplizer = cms.EDAnalyzer("ggNtuplizer",
                             hggPhotonIDConfiguration = hggPhotonIDCuts,
                             doGenParticles   = cms.bool(True),
                             runOnParticleGun = cms.bool(False),
                             genParticleSrc   = cms.InputTag("genParticles"),
                             generatorLabel   = cms.InputTag("generator"),
                             pileupCollection = cms.InputTag("addPileupInfo"),
                             VtxLabel         = cms.InputTag("offlinePrimaryVertices"),
                             VtxBSLabel       = cms.InputTag("offlinePrimaryVerticesWithBS"),
                             rhoLabel         = cms.InputTag("fixedGridRhoFastjetAll"),
                             electronSrc      = cms.InputTag("cleanPatElectrons"),
                             photonSrc        = cms.InputTag("cleanPatPhotons"),
                             ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
                             eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
                             esReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsES"),
                             recoPhotonSrc             = cms.InputTag("gedPhotons"),
                             TrackLabel                = cms.InputTag("generalTracks"),
                             gsfElectronLabel          = cms.InputTag("gsfElectrons"),
                             PFAllCandidates           = cms.InputTag("particleFlow")
)
