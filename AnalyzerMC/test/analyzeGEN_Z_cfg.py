#flake8: noqa

'''
    >>----------------------<<
    Diffractive Z GEN Analyzer
    >>----------------------<<
    
    Goal:
    Hadron Level analysis
    
    Usage:
    cmsRun analyzeGEN_Z_cfg.py
    
    Example:
    cmsRun analyzeGEN_Z_cfg.py Run=muon Type=diffractive
    
    Optional arguments:
    Run = muon or electron
    Type = diffractive or nondiffractive
    Authors: D. Figueiredo, R. Arciadiacono and N. Cartiglia
    
'''

import FWCore.ParameterSet.Config as cms
import os, sys
import atexit

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
options.register('Run','muon',VarParsing.multiplicity.singleton, VarParsing.varType.string,"Option to Run: muon or electron")
options.register('Type','diffractive',VarParsing.multiplicity.singleton, VarParsing.varType.string,"diffractive or nondiffractive")
options.parseArguments()

process = cms.Process('Analysis')

class config: pass
config.outputFile = "output.root" # output file name
config.Ebeam = 3500.0 # Energy only one beam. CM = 2*EBeam
config.NumberOfEvents = 5000 # number of events
config.debug = False # track errors
config.drawer = False # draw particle GEN chain

##
## PDG ID
##
## electron: 11
## positron: -11
## muon: 13
## antimuon: -13
##

if options.Run == "electron":
    print("")
    print("########")
    print("Electron")
    print("########")
    print("")
    config.pdgid1 = 11
    config.pdgid2 = -11
    if options.Type == "nondiffractive":
       config.input = "file:/storage1/dmf/TestSamples/DyToEEPU2010/DyToEE.root"
    elif options.Type == "diffractive":
       config.input = "file:/storage1/dmf/Samples/DiffractiveZ/GEN/pompyt_minus_Z_M20_cff_py_GEN_7TeV.root"
    else:
       print("")
       print("")
       raise RuntimeError, "Unknown option Type. EXIT! YOU NEED TO SETUP WITH ONE OF THE CORRECT OPTIONS."
       print("")

elif options.Run == "muon":
    print("")
    print("####")
    print("Muon")
    print("####")
    print("")
    config.pdgid1 = 13
    config.pdgid2 = -13
    if options.Type == "nondiffractive":
        config.input = "file:/storage1/dmf/TestSamples/DyToMuMuPU2010/DyToMuMu.root"
    elif options.Type == "diffractive":
        config.input = "file:/storage1/dmf/Samples/DiffractiveZ/GEN/pompyt_minus_Z_M20_cff_py_GEN_7TeV.root"
    else:
        print("")
        print("")
        raise RuntimeError, "Unknown option Type. EXIT! YOU NEED TO SETUP WITH ONE OF THE CORRECT OPTIONS."
        print("")

else:
    print("")
    print("")
    raise RuntimeError, "Unknown option Run. EXIT! YOU NEED TO SETUP WITH ONE OF THE CORRECT OPTIONS."
    print("")

print("")
print(">>> Input Options:")
print("Type of GEN Analysis: %s" % options.Run)
print("Type of process: %s" % options.Type)
print("PDG Id Lepton 1: %s" % config.pdgid1)
print("PDG Id Lepton 2: %s" % config.pdgid2)
print("Energy beam: %.2f" % config.Ebeam)
print("Debug: %s" % config.debug)
print("Output file: %s" % config.outputFile)
print("")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(config.NumberOfEvents)
)

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(config.input)
)

###process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
###process.load('PhysicsTools.HepMCCandAlgos.genParticles_cfi')
###process.genParticles.abortOnUnknownPDGCode = False

process.SDDY = cms.EDAnalyzer("SDDYAnalyzer",
	GenParticleTag = cms.InputTag("genParticles"),
	Particle1Id = cms.int32(config.pdgid1),
	Particle2Id = cms.int32(config.pdgid2),
    EBeam = cms.double(config.Ebeam),
	debug = cms.untracked.bool(config.debug)
)

process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
                                   src = cms.InputTag("genParticles"),
                                   printP4 = cms.untracked.bool(False),
                                   printPtEtaPhi = cms.untracked.bool(False),
                                   printVertex = cms.untracked.bool(False),
                                   printStatus = cms.untracked.bool(True),
                                   printIndex = cms.untracked.bool(True),
                                   status = cms.untracked.vint32(2,1)
                                   )


process.add_(cms.Service("TFileService",
                         fileName = cms.string(config.outputFile)
                         )
             )

###process.analysis = cms.Path(process.genParticles*process.SDDY)

if (config.drawer):
   process.analysis = cms.Path(process.genParticles*process.SDDY*process.printTree)
else:
   process.analysis = cms.Path(process.SDDY)
