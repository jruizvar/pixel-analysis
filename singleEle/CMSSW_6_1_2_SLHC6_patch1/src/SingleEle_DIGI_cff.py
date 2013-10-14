# Auto generated configuration file
# using: 
# Revision: 1.14 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step1 --filein /store/group/comm_trigger/L1TrackTrigger/BE5D_612_SLHC6_patch1/singleEle/SingleEle_GEN_SIM.root --fileout file:SingleEle_DIGI_RAW.root --pileup_input /store/group/comm_trigger/L1TrackTrigger/BE5D_612_SLHC6_patch1/MinBias/MinBias_BE5D_TuneZ2star_GEN_SIM_99.root --mc --eventcontent FEVTDEBUG --pileup AVE_140_BX_25ns --customise SLHCUpgradeSimulations/Configuration/combinedCustoms.cust_phase2_BE5D+noCrossing --datatier GEN-SIM-DIGI-RAW --conditions POSTLS261_V2::All --step DIGI,L1,DIGI2RAW --magField 38T_PostLS1 --python_filename SingleEle_DIGI_cff.py -n 1 --geometry ExtendedPhase2TkBE5D
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
    input = cms.untracked.int32(3)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('/store/group/comm_trigger/L1TrackTrigger/BE5D_612_SLHC6_patch1/singleEle/SingleEle_GEN_SIM.root')
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.14 $'),
    annotation = cms.untracked.string('step1 nevts:1'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    fileName = cms.untracked.string('file:SingleEle_DIGI_RAW.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW')
    )
)

# Additional output definition

# Other statements
from MinBiasSource import MinBias
process.mix.input.nbPileupEvents.averageNumber = cms.double(140.000000)
process.mix.bunchspace = cms.int32(25)
process.mix.minBunch = cms.int32(-12)
process.mix.maxBunch = cms.int32(3)
process.mix.input.fileNames = MinBias
#process.mix.input.fileNames = cms.untracked.vstring(['/store/group/comm_trigger/L1TrackTrigger/BE5D_612_SLHC6_patch1/MinBias/MinBias_BE5D_TuneZ2star_GEN_SIM_99.root'])
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
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1
from SLHCUpgradeSimulations.Configuration.phase2TkCustomsBE5D import customise,l1EventContent
from SLHCUpgradeSimulations.Configuration.HCalCustoms import customise_HcalPhase0
from SLHCUpgradeSimulations.Configuration.customise_mixing import customise_NoCrossing
process=customisePostLS1(process)
process=customise(process)
process=l1EventContent(process)
process=customise_HcalPhase0(process)
process=customise_NoCrossing(process)

# End of customisation functions
