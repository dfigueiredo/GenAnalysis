import FWCore.ParameterSet.Config as cms

process = cms.Process('Analysis')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
 #      fileNames = cms.untracked.vstring("file:/storage2/eliza/TOTEM/pompyt_jpsiminus_GEN1509.root")
 #        fileNames = cms.untracked.vstring("file:/storage2/eliza/TOTEM/pompyt_jpsiplus_GEN.root")
 #        fileNames = cms.untracked.vstring("file:pompyt_jpsiPlus_GEN_TEST.root")
 #        fileNames = cms.untracked.vstring("file:pompyt_jpsiminus_GEN_13TeV.root")

)

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('PhysicsTools.HepMCCandAlgos.genParticles_cfi')

process.genParticles.abortOnUnknownPDGCode = False

process.Jpsi = cms.EDAnalyzer("SDJpsiAnalyzer",      #"SDDYAnalyzer",
	GenParticleTag = cms.InputTag("genParticles"),
	Particle1Id = cms.int32(13),
	Particle2Id = cms.int32(13),
	debug = cms.untracked.bool(True)
)

process.add_(cms.Service("TFileService",
               fileName = cms.string("histo_pompyt_jpsiplus_GEN_8TeV.root")
		#fileName = cms.string("Pompyt_SDMinusMuMu_histos_2609_TEST.root")
	)
)

process.analysis = cms.Path(process.genParticles*process.Jpsi)
