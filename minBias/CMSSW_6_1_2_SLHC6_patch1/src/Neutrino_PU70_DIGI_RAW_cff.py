# Auto generated configuration file
# using: 
# Revision: 1.14 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step2 -s DIGI,L1,DIGI2RAW --python_filename Neutrino_PU35_DIGI_RAW_cff.py --no_exec --fileout file:Neutrino_PU35_DIGI_RAW.root --magField 38T_PostLS1 --geometry ExtendedPhase2TkBE5D --conditions POSTLS261_V2::All --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1,SLHCUpgradeSimulations/Configuration/HCalCustoms.customise_HcalPhase0,SLHCUpgradeSimulations/Configuration/phase2TkCustomsBE5D.customise,SLHCUpgradeSimulations/Configuration/phase2TkCustomsBE5D.l1EventContent,SLHCUpgradeSimulations/Configuration/customise_mixing.customise_NoCrossing --eventcontent FEVTDEBUG --datatier GEN-SIM-DIGI-RAW --filein /Neutrino_Pt2to20_gun/UpgFall13-POSTLS261_V2-v1/GEN-SIM --pileup_input MinBias --pileup AVE_140_BX_25ns -n 1
import FWCore.ParameterSet.Config as cms

process = cms.Process('DIGI2RAW')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mix_E8TeV_AVE_16_BX_25ns_cfi')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('/Neutrino_Pt2to20_gun/UpgFall13-POSTLS261_V2-v1/GEN-SIM')
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.14 $'),
    annotation = cms.untracked.string('step2 nevts:1'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    fileName = cms.untracked.string('file:Neutrino_PU70_DIGI_RAW.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW')
    )
)

# Additional output definition

# Other statements
from MinBiasSource import MinBias
process.mix.input.nbPileupEvents.averageNumber = cms.double(70.000000)
process.mix.bunchspace = cms.int32(25)
process.mix.minBunch = cms.int32(-12)
process.mix.maxBunch = cms.int32(3)
process.mix.input.fileNames = MinBias
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS261_V2::All', '')

# Path and EndPath definitions
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Schedule definition
process.schedule = cms.Schedule(process.digitisation_step,process.L1simulation_step,process.digi2raw_step,process.endjob_step,process.FEVTDEBUGoutput_step)

# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.phase2TkCustomsBE5D
from SLHCUpgradeSimulations.Configuration.phase2TkCustomsBE5D import customise,l1EventContent 

#call to customisation function customise imported from SLHCUpgradeSimulations.Configuration.phase2TkCustomsBE5D
process = customise(process)

#call to customisation function l1EventContent imported from SLHCUpgradeSimulations.Configuration.phase2TkCustomsBE5D
process = l1EventContent(process)

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.HCalCustoms
from SLHCUpgradeSimulations.Configuration.HCalCustoms import customise_HcalPhase0 

#call to customisation function customise_HcalPhase0 imported from SLHCUpgradeSimulations.Configuration.HCalCustoms
process = customise_HcalPhase0(process)

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1 

#call to customisation function customisePostLS1 imported from SLHCUpgradeSimulations.Configuration.postLS1Customs
process = customisePostLS1(process)

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.customise_mixing
from SLHCUpgradeSimulations.Configuration.customise_mixing import customise_NoCrossing 

#call to customisation function customise_NoCrossing imported from SLHCUpgradeSimulations.Configuration.customise_mixing
process = customise_NoCrossing(process)

# End of customisation functions
