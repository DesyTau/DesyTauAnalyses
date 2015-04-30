import FWCore.ParameterSet.Config as cms

#isData = False
#skim = 0
#year = 2015
#period = 'PHYS14'
#sampleName = 'MonteCarlo'
# sampleName = 'TTJets', "QCD", "DYJetsToLL_M50"

#if "@SKIM@".lower()=="0":
#    skim = 0

#sampleName = "@SAMPLE_NAME@"

process = cms.Process("TreeProducer")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
#process.load('Configuration.StandardSequences.Geometry_cff')
#process.load('Configuration.Geometry.GeometryIdeal_cff')
#process.load('Configuration.StandardSequences.MagneticField_cff')
#process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

process.load("Configuration.Geometry.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.GlobalTag.globaltag = cms.string('PHYS14_25_V2::All')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
# MINIAOD
# WJetsToLNu
#        '/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/02215B44-2D70-E411-90A3-0025905A60B8.root',
#        '/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0603D444-2D70-E411-AF03-002618943922.root',
#        '/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/08947C88-3570-E411-974E-002618FDA26D.root',
#        '/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0E8B81D9-2D70-E411-94AB-0025905A4964.root'
# 
# DYJetsToLL_M-50
        '/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0432E62A-7A6C-E411-87BB-002590DB92A8.root',
        '/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/06C61714-7E6C-E411-9205-002590DB92A8.root',
        '/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0EAD09A8-7C6C-E411-B903-0025901D493E.root',
        '/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1E4D0DAE-7C6C-E411-B488-002590DB923C.root'
#
# TTJets
#        '/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00C90EFC-3074-E411-A845-002590DB9262.root',
#        '/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00D3EAF1-3174-E411-A5B2-0025904B144E.root',
#        '/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/02EF3EFC-0475-E411-A9DB-002590DB9166.root'
#
# QCD_Pt-50to80  
#        '/store/mc/Phys14DR/QCD_Pt-50to80_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_castor_PHYS14_25_V1-v2/00000/164247BB-6577-E411-A063-00259073E50A.root',
#        '/store/mc/Phys14DR/QCD_Pt-50to80_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_castor_PHYS14_25_V1-v2/00000/1EF05EBB-6577-E411-BD3A-20CF305616EF.root',
#        '/store/mc/Phys14DR/QCD_Pt-50to80_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_castor_PHYS14_25_V1-v2/00000/24BF07BB-6577-E411-89E4-00259073E3AE.root',
#        '/store/mc/Phys14DR/QCD_Pt-50to80_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_castor_PHYS14_25_V1-v2/00000/303020BC-6577-E411-A7A7-00259073E516.root'


    ),
##    dropDescendantsOfDroppedBranches=cms.untracked.bool(False),
##    inputCommands=cms.untracked.vstring(
##        'keep *',
##        'drop patTaus_*_*_*',
##        'drop *PFTau*_*_*_*'
##    )
)
 
#####################################################
  
#--------------------------------------------------------------------------------
# produce PAT-tuple
process.load("PhysicsTools/PatAlgos/patSequences_cff")
# configure pat::Jet production
# (enable L2L3Residual corrections in case running on Data)
#jetCorrections = ( 'L1FastJet', 'L2Relative', 'L3Absolute')
#from PhysicsTools.PatAlgos.tools.jetTools import *
#switchJetCollection(
#    process,
#    jetSource = cms.InputTag('ak4PFJetsCHS'),
#    jetCorrections = ( 'AK4PFchs', jetCorrections, "" ),
#    outputModules = []
#)

#process.patJets.addTagInfos = cms.bool(True)
#process.patJets.tagInfoSources = cms.VInputTag("impactParameterTagInfosAOD","secondaryVertexTagInfosAOD","softMuonTagInfosAOD")
#process.patJets.tagInfoSources = cms.VInputTag("secondaryVertexTagInfosEI")


#--------------------------------------------------------------------------------
# switch to HPS PFTaus (and disable all "cleaning" cuts)
#from PhysicsTools.PatAlgos.tools.tauTools import *
#switchToPFTauHPS(process)

#--------------------------------------------------------------------------------
# select "good" reconstructed vertices
#process.load("TauAnalysis/RecoTools/recoVertexSelection_cff")

#--------------------------------------------------------------------------------
# electron Id (CSA14) ->
# input tags for electron id are hardcoded
# in class MyRootMaker    
#process.load("EgammaAnalysis.ElectronTools.electronIdMVAProducer_CSA14_cfi")
#process.load("RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cff")

# specify correct sources for MVA electron Id.
#process.electronIDValueMapProducer.ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB")
#process.electronIDValueMapProducer.eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE")
#process.electronIDValueMapProducer.esReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsES")

#----------------------------------------------------------------------------------
# produce PU Jet Ids & MVA MET
#----------------------------------------------------------------------------------
process.tauPreSelectionDiTau = cms.EDFilter("PATTauSelector",
    src = cms.InputTag("slimmedTaus"),
    cut = cms.string('pt > 40. && abs(eta) < 2.5 && tauID("decayModeFindingNewDMs") > 0.5')
)
 
 
process.tauPreSelectionTauEle = cms.EDFilter("PATTauSelector",
    src = cms.InputTag("slimmedTaus"),
    cut = cms.string('pt > 15. && abs(eta) < 2.5 && tauID("decayModeFinding") > 0.5')
)
 
 
process.tauPreSelectionTauMu = cms.EDFilter("PATTauSelector",
    src = cms.InputTag("slimmedTaus"),
    cut = cms.string('pt > 15. && abs(eta) < 2.5 && tauID("decayModeFinding") > 0.5')
)

process.muonPreSelectionMuEle = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("slimmedMuons"),
    cut = cms.string('pt > 9. && abs(eta) < 2.5 && isPFMuon && (isGlobalMuon || isTrackerMuon)')
)
 
 
process.muonPreSelectionTauMu = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("slimmedMuons"),
    cut = cms.string('pt > 20. && abs(eta) < 2.5 && isPFMuon && (isGlobalMuon || isTrackerMuon)')
)

process.electronPreSelectionMuEle = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("slimmedElectrons"),
    cut = cms.string('pt > 13. && abs(eta) < 2.5')
)
 
 
process.electronPreSelectionTauEle = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("slimmedElectrons"),
    cut = cms.string('pt > 20. && abs(eta) < 2.5')
)

process.leptonSkimSequence = cms.Sequence(process.electronPreSelectionMuEle + process.electronPreSelectionTauEle + process.muonPreSelectionTauMu + process.muonPreSelectionMuEle + process.tauPreSelectionTauMu + process.tauPreSelectionTauEle + process.tauPreSelectionDiTau)

process.ak4PFJets = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(1),
    doAreaFastjet = cms.bool(True),
    voronoiRfact = cms.double(-0.9),
    maxBadHcalCells = cms.uint32(9999999),
    doAreaDiskApprox = cms.bool(False),
    maxRecoveredEcalCells = cms.uint32(9999999),
    jetType = cms.string('PFJet'),
    minSeed = cms.uint32(14327),
    Ghost_EtaMax = cms.double(5.0),
    doRhoFastjet = cms.bool(False),
    jetAlgorithm = cms.string('AntiKt'),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.4),
    maxBadEcalCells = cms.uint32(9999999),
    useDeterministicSeed = cms.bool(True),
    doPVCorrection = cms.bool(False),
    maxRecoveredHcalCells = cms.uint32(9999999),
    rParam = cms.double(0.4),
    maxProblematicHcalCells = cms.uint32(9999999),
    doOutputJets = cms.bool(True),
    src = cms.InputTag("packedPFCandidates"),
    inputEtMin = cms.double(0.0),
    srcPVs = cms.InputTag(""),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    maxProblematicEcalCells = cms.uint32(9999999),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0)
)

process.calibratedAK4PFJetsForPFMVAMEt = cms.EDProducer("PFJetCorrectionProducer",
    src = cms.InputTag("ak4PFJets"),
    correctors = cms.vstring('ak4PFL1FastL2L3')
)

process.puJetIdForPFMVAMEt = cms.EDProducer("PileupJetIdProducer",
    algos = cms.VPSet(cms.PSet(
        tmvaVariables = cms.vstring('nvtx', 
            'jetPt', 
            'jetEta', 
            'jetPhi', 
            'dZ', 
            'beta', 
            'betaStar', 
            'nCharged', 
            'nNeutrals', 
            'dR2Mean', 
            'ptD', 
            'frac01', 
            'frac02', 
            'frac03', 
            'frac04', 
            'frac05'),
        tmvaMethod = cms.string('JetID'),
        cutBased = cms.bool(False),
        tmvaWeights = cms.string('RecoJets/JetProducers/data/TMVAClassificationCategory_JetID_MET_53X_Dec2012.weights.xml.gz'),
        tmvaSpectators = cms.vstring(),
        label = cms.string('full'),
        version = cms.int32(-1),
        JetIdParams = cms.PSet(
            Pt2030_Tight = cms.vdouble(0.3, 0.4, 0.7, 0.8),
            Pt2030_Loose = cms.vdouble(0.0, 0.0, 0.2, 0.6),
            Pt3050_Medium = cms.vdouble(0.3, 0.2, 0.7, 0.8),
            Pt1020_Tight = cms.vdouble(-0.2, 0.2, 0.2, 0.6),
            Pt2030_Medium = cms.vdouble(0.2, 0.2, 0.5, 0.7),
            Pt010_Tight = cms.vdouble(0.5, 0.6, 0.6, 0.9),
            Pt1020_Loose = cms.vdouble(-0.4, -0.4, -0.4, 0.4),
            Pt010_Medium = cms.vdouble(0.2, 0.4, 0.2, 0.6),
            Pt1020_Medium = cms.vdouble(-0.3, 0.0, 0.0, 0.5),
            Pt010_Loose = cms.vdouble(0.0, 0.0, 0.0, 0.2),
            Pt3050_Loose = cms.vdouble(0.0, 0.0, 0.6, 0.2),
            Pt3050_Tight = cms.vdouble(0.5, 0.4, 0.8, 0.9)
        ),
        impactParTkThreshold = cms.double(0.0)
    )),
    inputIsCorrected = cms.bool(True),
    vertexes = cms.InputTag("offlineSlimmedPrimaryVertices"),
    produceJetIds = cms.bool(True),
    jec = cms.string('AK4PF'),
    residualsFromTxt = cms.bool(False),
    applyJec = cms.bool(True),
    jetids = cms.InputTag("pileupJetIdCalculator"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    jets = cms.InputTag("calibratedAK4PFJetsForPFMVAMEt"),
    runMvas = cms.bool(True)
)

process.mvaMETDiTau = cms.EDProducer("PFMETProducerMVA",
    useType1 = cms.bool(True),
    srcLeptons = cms.VInputTag(cms.InputTag("tauPreSelectionDiTau"), cms.InputTag("tauPreSelectionDiTau")),
    verbosity = cms.int32(0),
    srcCorrJets = cms.InputTag("calibratedAK4PFJetsForPFMVAMEt"),
    srcVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    corrector = cms.string('ak4PFL1Fastjet'),
    loadMVAfromDB = cms.bool(False),
    dZcut = cms.double(0.1),
    inputFileNames = cms.PSet(
        DPhi = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrphi_7_2_X_MINIAOD_BX25PU20_Mar2015.root'),
        CovU2 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru2cov_7_2_X_MINIAOD_BX25PU20_Mar2015.root'),
        U = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrmet_7_2_X_MINIAOD_BX25PU20_Mar2015.root'),
        CovU1 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru1cov_7_2_X_MINIAOD_BX25PU20_Mar2015.root')
    ),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll"),
    srcMVAPileupJetId = cms.InputTag("puJetIdForPFMVAMEt","fullDiscriminant"),
    srcUncorrJets = cms.InputTag("ak4PFJets"),
    minNumLeptons = cms.int32(0),
    globalThreshold = cms.double(-1.0),
    permuteLeptons = cms.bool(True),
    minCorrJetPt = cms.double(-1.0),
    inputRecords = cms.PSet(
        DPhi = cms.string('PhiCor'),
        CovU2 = cms.string('CovU2'),
        U = cms.string('RecoilCor'),
        CovU1 = cms.string('CovU1')
    ),
    srcPFCandidates = cms.InputTag("packedPFCandidates")
)

process.mvaMETTauMu = cms.EDProducer("PFMETProducerMVA",
    useType1 = cms.bool(True),
    srcLeptons = cms.VInputTag(cms.InputTag("tauPreSelectionTauMu"), cms.InputTag("muonPreSelectionTauMu")),
    verbosity = cms.int32(0),
    srcCorrJets = cms.InputTag("calibratedAK4PFJetsForPFMVAMEt"),
    srcVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    corrector = cms.string('ak4PFL1Fastjet'),
    loadMVAfromDB = cms.bool(False),
    dZcut = cms.double(0.1),
    inputFileNames = cms.PSet(
        DPhi = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrphi_7_2_X_MINIAOD_BX25PU20_Mar2015.root'),
        CovU2 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru2cov_7_2_X_MINIAOD_BX25PU20_Mar2015.root'),
        U = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrmet_7_2_X_MINIAOD_BX25PU20_Mar2015.root'),
        CovU1 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru1cov_7_2_X_MINIAOD_BX25PU20_Mar2015.root')
    ),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll"),
    srcMVAPileupJetId = cms.InputTag("puJetIdForPFMVAMEt","fullDiscriminant"),
    srcUncorrJets = cms.InputTag("ak4PFJets"),
    minNumLeptons = cms.int32(0),
    globalThreshold = cms.double(-1.0),
    permuteLeptons = cms.bool(True),
    minCorrJetPt = cms.double(-1.0),
    inputRecords = cms.PSet(
        DPhi = cms.string('PhiCor'),
        CovU2 = cms.string('CovU2'),
        U = cms.string('RecoilCor'),
        CovU1 = cms.string('CovU1')
    ),
    srcPFCandidates = cms.InputTag("packedPFCandidates")
)

process.mvaMETTauEle = cms.EDProducer("PFMETProducerMVA",
    useType1 = cms.bool(True),
    srcLeptons = cms.VInputTag(cms.InputTag("tauPreSelectionTauEle"), cms.InputTag("electronPreSelectionTauEle")),
    verbosity = cms.int32(0),
    srcCorrJets = cms.InputTag("calibratedAK4PFJetsForPFMVAMEt"),
    srcVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    corrector = cms.string('ak4PFL1Fastjet'),
    loadMVAfromDB = cms.bool(False),
    dZcut = cms.double(0.1),
    inputFileNames = cms.PSet(
        DPhi = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrphi_7_2_X_MINIAOD_BX25PU20_Mar2015.root'),
        CovU2 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru2cov_7_2_X_MINIAOD_BX25PU20_Mar2015.root'),
        U = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrmet_7_2_X_MINIAOD_BX25PU20_Mar2015.root'),
        CovU1 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru1cov_7_2_X_MINIAOD_BX25PU20_Mar2015.root')
    ),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll"),
    srcMVAPileupJetId = cms.InputTag("puJetIdForPFMVAMEt","fullDiscriminant"),
    srcUncorrJets = cms.InputTag("ak4PFJets"),
    minNumLeptons = cms.int32(0),
    globalThreshold = cms.double(-1.0),
    permuteLeptons = cms.bool(True),
    minCorrJetPt = cms.double(-1.0),
    inputRecords = cms.PSet(
        DPhi = cms.string('PhiCor'),
        CovU2 = cms.string('CovU2'),
        U = cms.string('RecoilCor'),
        CovU1 = cms.string('CovU1')
    ),
    srcPFCandidates = cms.InputTag("packedPFCandidates")
)

process.mvaMETMuEle = cms.EDProducer("PFMETProducerMVA",
    useType1 = cms.bool(True),
    srcLeptons = cms.VInputTag(cms.InputTag("muonPreSelectionMuEle"), cms.InputTag("electronPreSelectionMuEle")),
    verbosity = cms.int32(0),
    srcCorrJets = cms.InputTag("calibratedAK4PFJetsForPFMVAMEt"),
    srcVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    corrector = cms.string('ak4PFL1Fastjet'),
    loadMVAfromDB = cms.bool(False),
    dZcut = cms.double(0.1),
    inputFileNames = cms.PSet(
        DPhi = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrphi_7_2_X_MINIAOD_BX25PU20_Mar2015.root'),
        CovU2 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru2cov_7_2_X_MINIAOD_BX25PU20_Mar2015.root'),
        U = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrmet_7_2_X_MINIAOD_BX25PU20_Mar2015.root'),
        CovU1 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru1cov_7_2_X_MINIAOD_BX25PU20_Mar2015.root')
    ),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll"),
    srcMVAPileupJetId = cms.InputTag("puJetIdForPFMVAMEt","fullDiscriminant"),
    srcUncorrJets = cms.InputTag("ak4PFJets"),
    minNumLeptons = cms.int32(0),
    globalThreshold = cms.double(-1.0),
    permuteLeptons = cms.bool(True),
    minCorrJetPt = cms.double(-1.0),
    inputRecords = cms.PSet(
        DPhi = cms.string('PhiCor'),
        CovU2 = cms.string('CovU2'),
        U = cms.string('RecoilCor'),
        CovU1 = cms.string('CovU1')
    ),
    srcPFCandidates = cms.InputTag("packedPFCandidates")
)

process.patMvaMetDiTau = process.patMETs.clone(
    metSource = cms.InputTag('mvaMETDiTau'),
    addMuonCorrections = cms.bool(False),
    genMETSource = cms.InputTag('genMetTrue'),
    addGenMET = cms.bool(False)
    )
process.patMvaMetTauMu = process.patMETs.clone(
    metSource = cms.InputTag('mvaMETTauMu'),
    addMuonCorrections = cms.bool(False),
    genMETSource = cms.InputTag('genMetTrue'),
    addGenMET = cms.bool(False)
    )
process.patMvaMetTauEle = process.patMETs.clone(
    metSource = cms.InputTag('mvaMETTauEle'),
    addMuonCorrections = cms.bool(False),
    genMETSource = cms.InputTag('genMetTrue'),
    addGenMET = cms.bool(False)
    )
process.patMvaMetMuEle = process.patMETs.clone(
    metSource = cms.InputTag('mvaMETMuEle'),
    addMuonCorrections = cms.bool(False),
    genMETSource = cms.InputTag('genMetTrue'),
    addGenMET = cms.bool(False)
    )

process.mvaMetSequence  = cms.Sequence(process.leptonSkimSequence* process.ak4PFJets * process.calibratedAK4PFJetsForPFMVAMEt * 
                                       process.puJetIdForPFMVAMEt * process.mvaMETDiTau * process.mvaMETTauMu * 
                                       process.mvaMETTauEle * process.mvaMETMuEle * process.patMvaMetDiTau * 
                                       process.patMvaMetTauMu * process.patMvaMetTauEle * process.patMvaMetMuEle)
'''
process.pileupJetIdFull = cms.EDProducer("PileupJetIdProducer",
    produceJetIds = cms.bool(True),
    runMvas = cms.bool(True),
    inputIsCorrected = cms.bool(False),
    vertexes = cms.InputTag("offlineSlimmedPrimaryVertices"),
    residualsTxt = cms.FileInPath('RecoJets/JetProducers/data/download.url'),
    jec = cms.string('AK5PF'),
    residualsFromTxt = cms.bool(False),
    applyJec = cms.bool(True),
    jetids = cms.InputTag(""),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    jets = cms.InputTag("ak4PFJets"),
    algos = cms.VPSet(cms.PSet(
        tmvaVariables = cms.vstring('nvtx', 
            'dZ', 
            'beta', 
            'betaStar', 
            'nCharged', 
            'nNeutrals', 
            'dR2Mean', 
            'ptD', 
            'frac01', 
            'frac02', 
            'frac03', 
            'frac04', 
            'frac05'),
        tmvaMethod = cms.string('JetIDMVAHighPt'),
        cutBased = cms.bool(False),
        #tmvaWeights = cms.string('CondFormats/JetMETObjects/data/TMVAClassificationCategory_JetID_53X_Dec2012.weights.xml'),
        tmvaWeights = cms.string('RecoJets/JetProducers/data/TMVAClassificationCategory_JetID_MET_53X_Dec2012.weights.xml.gz'),
        tmvaSpectators = cms.vstring('jetPt', 
            'jetEta', 
            'jetPhi'),
        label = cms.string('full53x'),
        version = cms.int32(-1),
        JetIdParams = cms.PSet(
            Pt2030_Tight = cms.vdouble(0.73, 0.05, -0.26, -0.42),
            Pt2030_Loose = cms.vdouble(-0.63, -0.6, -0.55, -0.45),
            Pt3050_Medium = cms.vdouble(0.1, -0.36, -0.54, -0.54),
            Pt1020_MET = cms.vdouble(0.3, -0.2, -0.4, -0.4),
            Pt2030_Medium = cms.vdouble(0.1, -0.36, -0.54, -0.54),
            Pt010_Tight = cms.vdouble(-0.83, -0.81, -0.74, -0.81),
            Pt1020_Tight = cms.vdouble(-0.83, -0.81, -0.74, -0.81),
            Pt3050_MET = cms.vdouble(0.0, 0.0, -0.1, -0.2),
            Pt010_MET = cms.vdouble(0.0, -0.6, -0.4, -0.4),
            Pt1020_Loose = cms.vdouble(-0.95, -0.96, -0.94, -0.95),
            Pt010_Medium = cms.vdouble(-0.83, -0.92, -0.9, -0.92),
            Pt1020_Medium = cms.vdouble(-0.83, -0.92, -0.9, -0.92),
            Pt2030_MET = cms.vdouble(0.0, 0.0, 0.0, 0.0),
            Pt010_Loose = cms.vdouble(-0.95, -0.96, -0.94, -0.95),
            Pt3050_Loose = cms.vdouble(-0.63, -0.6, -0.55, -0.45),
            Pt3050_Tight = cms.vdouble(0.73, 0.05, -0.26, -0.42)
        ),
        impactParTkThreshold = cms.double(1.0)
    ))
)
process.puJetIdSequence = cms.Sequence(process.pileupJetIdFull)
'''
##################################################
# Main
process.makeroottree = cms.EDAnalyzer("NTupleMaker",
# data, year, period, skim
IsData = cms.untracked.bool(False),
Year = cms.untracked.uint32(2015),
Period = cms.untracked.string("PHYS14"),
Skim = cms.untracked.uint32(0),
# switches of collections
GenParticles = cms.untracked.bool(True),
Trigger = cms.untracked.bool(True),
RecPrimVertex = cms.untracked.bool(True),
RecBeamSpot = cms.untracked.bool(True),
RecTrack = cms.untracked.bool(False),
RecPFMet = cms.untracked.bool(True),
RecMvaMet = cms.untracked.bool(True),                                      
RecMuon = cms.untracked.bool(True),
RecPhoton = cms.untracked.bool(False),
RecElectron = cms.untracked.bool(True),
RecTau = cms.untracked.bool(True),
RecJet = cms.untracked.bool(True),
# collections
MuonCollectionTag = cms.InputTag("slimmedMuons"), 
ElectronCollectionTag = cms.InputTag("slimmedElectrons"),
TauCollectionTag = cms.InputTag("slimmedTaus"),
JetCollectionTag = cms.InputTag("slimmedJets"),
MetCollectionTag = cms.InputTag("slimmedMETs"),
MvaMetCollectionsTag = cms.VInputTag("patMvaMetDiTau", "patMvaMetTauMu", "patMvaMetTauEle", "patMvaMetMuEle"),
TrackCollectionTag = cms.InputTag("generalTracks"),
GenParticleCollectionTag = cms.InputTag("prunedGenParticles"),
TriggerObjectCollectionTag = cms.InputTag("selectedPatTrigger"),
BeamSpotCollectionTag =  cms.InputTag("offlineBeamSpot"),
PVCollectionTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
# trigger info
HLTriggerPaths = cms.untracked.vstring(
'HLT_IsoMu24_eta2p1_IterTrk02_v',
'HLT_IsoTkMu24_eta2p1_IterTrk02_v',
'HLT_Ele27_eta2p1_WP85_Gsf_v',
'HLT_Ele32_eta2p1_WP85_Gsf_v',
'HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v',
'HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v',
'HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v',
'HLT_PFMET170_NoiseCleaned_v'
),
TriggerProcess = cms.untracked.string("HLT"),
# tracks
RecTrackPtMin = cms.untracked.double(0.5),
RecTrackEtaMax = cms.untracked.double(2.4),
RecTrackNum = cms.untracked.int32(0),
# muons
RecMuonPtMin = cms.untracked.double(10.),
RecMuonEtaMax = cms.untracked.double(2.4),
RecMuonHLTriggerMatching = cms.untracked.vstring(
'HLT_IsoMu24_eta2p1_IterTrk02_v.*:hltL3crIsoL1sMu20Eta2p1L1f0L2f20QL3f24QL3crIsoRhoFiltered0p15IterTrk02',
'HLT_IsoTkMu24_eta2p1_IterTrk02_v.*:hltL3fL1sMu20L1Eta2p1f0TkFiltered24QL3crIsoRhoFiltered0p15IterTrk02',
'HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v.*:hltL1Mu12EG7L3IsoMuFiltered23',
'HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v.*:hltL1sL1Mu5EG20ORL1Mu5IsoEG18L3IsoFiltered8'
),
RecMuonNum = cms.untracked.int32(0),
# photons
RecPhotonPtMin = cms.untracked.double(20.),
RecPhotonEtaMax = cms.untracked.double(10000.),
RecPhotonHLTriggerMatching = cms.untracked.vstring(),
RecPhotonNum = cms.untracked.int32(0),
# electrons
RecElectronPtMin = cms.untracked.double(10.),
RecElectronEta = cms.untracked.double(2.4),
RecElectronHLTriggerMatching = cms.untracked.vstring(
'HLT_Ele27_eta2p1_WP85_Gsf_v.*:hltEle27WP85GsfTrackIsoFilter',
'HLT_Ele32_eta2p1_WP85_Gsf_v.*:hltEle32WP85GsfTrackIsoFilter',
'HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v.*:hltMu23Ele12GsfTrackIsoLegEle12GsfCaloIdTrackIdIsoMediumWPFilter',
'HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v.*:hltMu8Ele23GsfTrackIsoLegEle23GsfCaloIdTrackIdIsoMediumWPFilter'
),
RecElectronNum = cms.untracked.int32(0),
# taus
RecTauPtMin = cms.untracked.double(20),
RecTauEtaMax = cms.untracked.double(2.3),                                      
RecTauHLTriggerMatching = cms.untracked.vstring(
'HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v.*:hltPFTau50TrackPt30LooseAbsOrRelIso'
),
RecTauFloatDiscriminators = cms.untracked.vstring(
#'againstElectronLoose',
#'againstElectronLooseMVA5',
#'againstElectronMVA5category',
#'againstElectronMVA5raw',
#'againstElectronMedium',
#'againstElectronMediumMVA5',
#'againstElectronTight',
#'againstElectronTightMVA5',
#'againstElectronVLooseMVA5',
#'againstElectronVTightMVA5',
#'againstMuonLoose',
#'againstMuonLoose2',
#'againstMuonLoose3',
#'againstMuonLooseMVA',
#'againstMuonMVAraw',
#'againstMuonMedium',
#'againstMuonMedium2',
#'againstMuonMediumMVA',
#'againstMuonTight',
#'againstMuonTight2',
#'againstMuonTight3',
#'againstMuonTightMVA',
#'byCombinedIsolationDeltaBetaCorrRaw3Hits',
#'byIsolationMVA3newDMwLTraw',
#'byIsolationMVA3newDMwoLTraw',
#'byIsolationMVA3oldDMwLTraw',
#'byIsolationMVA3oldDMwoLTraw',
#'byLooseCombinedIsolationDeltaBetaCorr3Hits',
#'byLooseIsolationMVA3newDMwLT',
#'byLooseIsolationMVA3newDMwoLT',
#'byLooseIsolationMVA3oldDMwLT',
#'byLooseIsolationMVA3oldDMwoLT',
#'byMediumCombinedIsolationDeltaBetaCorr3Hits',
#'byMediumIsolationMVA3newDMwLT',
#'byMediumIsolationMVA3newDMwoLT',
#'byMediumIsolationMVA3oldDMwLT',
#'byMediumIsolationMVA3oldDMwoLT',
#'byTightCombinedIsolationDeltaBetaCorr3Hits',
#'byTightIsolationMVA3newDMwLT',
#'byTightIsolationMVA3newDMwoLT',
#'byTightIsolationMVA3oldDMwLT',
#'byTightIsolationMVA3oldDMwoLT',
#'byVLooseIsolationMVA3newDMwLT',
#'byVLooseIsolationMVA3newDMwoLT',
#'byVLooseIsolationMVA3oldDMwLT',
#'byVLooseIsolationMVA3oldDMwoLT',
#'byVTightIsolationMVA3newDMwLT',
#'byVTightIsolationMVA3newDMwoLT',
#'byVTightIsolationMVA3oldDMwLT',
#'byVTightIsolationMVA3oldDMwoLT',
#'byVVTightIsolationMVA3newDMwLT',
#'byVVTightIsolationMVA3newDMwoLT',
#'byVVTightIsolationMVA3oldDMwLT',
#'byVVTightIsolationMVA3oldDMwoLT',
#'chargedIsoPtSum',
#'decayModeFinding',
#'decayModeFindingNewDMs',
#'neutralIsoPtSum',
#'puCorrPtSum'
),
RecTauBinaryDiscriminators = cms.untracked.vstring(),
RecTauNum = cms.untracked.int32(0),
# jets
RecJetPtMin = cms.untracked.double(30.),
RecJetEtaMax = cms.untracked.double(2.5),
RecJetHLTriggerMatching = cms.untracked.vstring(),
RecJetBtagDiscriminators = cms.untracked.vstring(
'jetBProbabilityBJetTags',
'jetProbabilityBJetTags',
'trackCountingHighPurBJetTags',
'trackCountingHighEffBJetTags',
'simpleSecondaryVertexHighEffBJetTags',
'simpleSecondaryVertexHighPurBJetTags',
'combinedInclusiveSecondaryVertexV2BJetTags',
'pfCombinedSecondaryVertexBJetTags',
'combinedMVABJetTags'
),
RecJetNum = cms.untracked.int32(0),
SampleName = cms.untracked.string("MC") 
)

#process.patJets.addBTagInfo = cms.bool(True)

process.p = cms.Path(
#                     process.particleFlowPtrs +
#                     process.makePatElectrons +
#                     process.makePatMuons + 
#                     process.makePatJets +
#	              process.selectPrimaryVertex *
#                     process.patDefaultSequence * 
#                     process.egmGsfElectronIDSequence *
#                     process.mvaTrigV025nsCSA14 * 
#                     process.mvaNonTrigV025nsCSA14 * 
    #process.ak4PFJets*
    process.mvaMetSequence*
    #process.puJetIdSequence*
    process.makeroottree
    )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("output.root")
                                   )

#processDumpFile = open('MyRootMaker.dump', 'w')
#print >> processDumpFile, process.dumpPython()




