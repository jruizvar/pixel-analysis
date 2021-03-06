import FWCore.ParameterSet.Config as cms

process = cms.Process('FASTSIMWDIGI')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('FastSimulation.Configuration.EventContent_cff')
#process.load('FastSimulation.PileUpProducer.PileUpSimulator_NoPileUp_cff')
process.load('SLHCUpgradeSimulations.Geometry.mixLowLumPU_FastSim14TeV_cff')
#process.load('FastSimulation.Configuration.Geometries_MC_cff')
process.load('FastSimulation.Configuration.Geometries_cff')
process.load('SLHCUpgradeSimulations.Geometry.Phase1_R39F16_cmsSimIdealGeometryXML_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('FastSimulation.Configuration.FamosSequences_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedParameters_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Include the RandomNumberGeneratorService definition
process.load("FastSimulation.Configuration.RandomServiceInitialization_cff")
process.RandomNumberGeneratorService.simSiStripDigis = cms.PSet(
      initialSeed = cms.untracked.uint32(1234567),
      engineName = cms.untracked.string('HepJamesRandom'))
process.RandomNumberGeneratorService.simSiPixelDigis = cms.PSet(
      initialSeed = cms.untracked.uint32(1234567),
      engineName = cms.untracked.string('HepJamesRandom'))

process.load('SLHCUpgradeSimulations.Geometry.Digi_Phase1_R39F16_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.load('SLHCUpgradeSimulations.Geometry.fakeConditions_Phase1_R39F16_cff')
process.load("SLHCUpgradeSimulations.Geometry.recoFromSimDigis_cff")
process.load("SLHCUpgradeSimulations.Geometry.upgradeTracking_phase1_cff")

## for fastsim we need these ################################
process.TrackerGeometricDetESModule.fromDDD=cms.bool(True)
process.TrackerDigiGeometryESModule.fromDDD=cms.bool(True)
process.simSiPixelDigis.ROUList =  ['famosSimHitsTrackerHits']
process.simSiStripDigis.ROUList =  ['famosSimHitsTrackerHits']
process.load("SimGeneral.TrackingAnalysis.trackingParticles_cfi")
process.mergedtruth.simHitCollections.tracker = ['famosSimHitsTrackerHits']
process.mergedtruth.simHitCollections.pixel = []
process.mergedtruth.simHitCollections.muon = []
process.mergedtruth.simHitLabel = 'famosSimHits'
## make occupancies more similar to full simulation
process.famosSimHits.ParticleFilter.etaMax = 3.0
process.famosSimHits.ParticleFilter.pTMin = 0.05
process.famosSimHits.TrackerSimHits.pTmin = 0.05
process.famosSimHits.TrackerSimHits.firstLoop = False
#############################################################
process.Timing =  cms.Service("Timing")

# If you want to turn on/off pile-up, default is no pileup
#process.famosPileUp.PileUpSimulator.averageNumber = 50.00
### if doing inefficiency at <PU>=50
#process.simSiPixelDigis.AddPixelInefficiency = 20

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.8 $'),
    annotation = cms.untracked.string('SLHCUpgradeSimulations/Configuration/python/FourMuPt_1_50_cfi.py nevts:10'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Other statements
process.GlobalTag.globaltag = 'DESIGN42_V11::All'

process.famosSimHits.SimulateCalorimetry = True
process.famosSimHits.SimulateTracking = True
process.GaussVtxSmearingParameters.type = cms.string("Gaussian")
process.famosSimHits.VertexGenerator = process.GaussVtxSmearingParameters
process.famosPileUp.VertexGenerator = process.GaussVtxSmearingParameters

####################
process.load("Configuration.Generator.PythiaUESettings_cfi")
process.generator = cms.EDFilter("Pythia6GeneratorFilter",
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaPylistVerbosity = cms.untracked.int32(0),
    filterEfficiency = cms.untracked.double(1.0),
    comEnergy = cms.double(14000.0),
    PythiaParameters = cms.PSet(
        process.pythiaUESettingsBlock,
        processParameters = cms.vstring('MSEL      = 0     ! User defined processes', 
            'MSUB(81)  = 1     ! qqbar to QQbar', 
            'MSUB(82)  = 1     ! gg to QQbar', 
            'MSTP(7)   = 6     ! flavour = top', 
            'PMAS(6,1) = 175.  ! top quark mass'),
        # This is a vector of ParameterSet names to be read, in this order
        parameterSets = cms.vstring('pythiaUESettings', 
            'processParameters')
    ),
    ExternalDecays = cms.PSet(
        Tauola = cms.untracked.PSet(
             UseTauolaPolarization = cms.bool(True),
             InputCards = cms.PSet
             (
                pjak1 = cms.int32(0),
                pjak2 = cms.int32(0),
                mdtau = cms.int32(0)
             )
        ),
        parameterSets = cms.vstring('Tauola')
    )
)

process.load('Configuration.StandardSequences.Validation_cff')

process.anal = cms.EDAnalyzer("EventContentAnalyzer")

process.load('FastSimulation.CaloRecHitsProducer.CaloRecHits_cff')
from FastSimulation.CaloRecHitsProducer.CaloRecHits_cff import *
from RecoLocalCalo.HcalRecAlgos.hcalRecAlgoESProd_cfi import *
# Calo Towers
from RecoJets.Configuration.CaloTowersRec_cff import *

# You may not want to simulate everything for your study
process.famosSimHits.SimulateCalorimetry = True
process.famosSimHits.SimulateTracking = False

#Load Scales
process.load("L1TriggerConfig.L1ScalesProducers.L1CaloInputScalesConfig_cff")
process.load("L1TriggerConfig.L1ScalesProducers.L1CaloScalesConfig_cff")

process.load("SLHCUpgradeSimulations.L1CaloTrigger.SLHCCaloTrigger_cff")

process.ecalRecHit.doDigis = True
process.hbhereco.doDigis = True
process.horeco.doDigis = True
process.hfreco.doDigis = True


process.load("SLHCUpgradeSimulations.L1CaloTrigger.SLHCCaloTriggerAnalysisCalibrated_cfi")

process.p1 = cms.Path(process.generator+
                      process.famosWithTrackerAndCaloHits+
                      process.simEcalTriggerPrimitiveDigis+
                      process.simHcalTriggerPrimitiveDigis+
                      process.SLHCCaloTrigger
)


#process.p1 = cms.Path(process.ProductionFilterSequence*process.famosWithEverything)
process.source = cms.Source("EmptySource")

# To write out events
process.load("FastSimulation.Configuration.EventContent_cff")
process.o1 = cms.OutputModule("PoolOutputModule",
                              outputCommands = cms.untracked.vstring('drop *_*_*_*',
                                                                     'keep *_L1Calo*_*_*',
                                                                     'keep *_*SLHCL1ExtraParticles*_*_*',
                                                                     'keep *_*l1extraParticles*_*_*',
                                                                     'keep *_genParticles_*_*',
                                                                     'keep recoGen*_*_*_*',
                                                                     'keep *_simEcalTriggerPrimitiveDigis_*_*',
                                                                     'keep *_simHcalTriggerPrimitiveDigis_*_*'
                                                                     ),
    fileName = cms.untracked.string('SLHCOutput.root')
)
process.outpath = cms.EndPath(process.o1)

# Make the job crash in case of missing product
process.options = cms.untracked.PSet( Rethrow = cms.untracked.vstring('ProductNotFound') )
