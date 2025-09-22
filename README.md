# PixelVertexingAlgos

### Prepare input file

    cmsrel CMSSW_15_1_0_pre6
    cd CMSSW_15_1_0_pre6/src/
    cmsenv
    git cms-init
    git cms-addpkg HLTrigger/NGTScouting 
    git cms-cherry-pick-pr 48935
    git cms-remote add elenavernazza
    git fetch elenavernazza
    git cms-cherry-pick 6d92324d2ed41e279615a59a24f006284ebc1db9
    scram b -j 0

    ALL_FILES_TTBAR='file:/eos/cms/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2580000/05a0431f-a326-47fa-aaaa-99ffb75200f0.root,file:/eos/cms/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2580000/096b5ea2-83af-4217-862e-b38c1c5ea504.root,file:/eos/cms/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2580000/0c1ceb85-79bf-4ab4-84e8-57170773551c.root,file:/eos/cms/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2580000/0d3c7ca2-352f-45d3-9218-78cb8c6c3611.root,file:/eos/cms/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2580000/0dab023d-1aca-4d6f-b402-bc529280a969.root,file:/eos/cms/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2580000/0f371acb-3d0f-4ed1-a71d-05d036a89766.root,file:/eos/cms/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2580000/10d79a7b-b04e-4971-ba10-099c657c1d71.root,file:/eos/cms/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2580000/187d121a-11e0-440b-a88a-3c1a177e08e4.root,file:/eos/cms/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2580000/19bfd771-c7a2-4728-83de-0715f493c715.root,file:/eos/cms/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2580000/28321f2f-ea4a-4b4c-85d0-0ca7939b0b15.root,file:/eos/cms/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2580000/283c2e1f-b2be-487c-8154-e11cc1152392.root,file:/eos/cms/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2580000/28764c84-3189-49a3-9096-00d08c1bcc73.root,file:/eos/cms/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2580000/291dbd7f-886f-4bac-b6b3-d866d47649eb.root,file:/eos/cms/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2580000/2a9305cc-3ff4-497e-80fa-eefe1a2fd149.root,file:/eos/cms/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2580000/2a9793f1-bd97-450c-98e6-86bb0f3f668e.root,file:/eos/cms/store/relval/CMSSW_15_1_0_pre5/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2580000/2f2128d9-494c-4d58-a23e-e7e14fe2bc73.root'

    cmsDriver.py step2 -s L1P2GT,HLT:NGTScouting,DQM:vertexingMonitorHLT+trackingMonitorHLT,VALIDATION:@hltValidation,NANO:@NGTScouting \
            --conditions auto:phase2_realistic_T33 \
            --datatier DQMIO,NANOAODSIM \
            -n 100 \
            --eventcontent DQMIO,NANOAODSIM \
            --geometry ExtendedRun4D110 \
            --era Phase2C17I13M9 \
            --procModifier alpaka,ngtScouting \
            --filein $ALL_FILES_TTBAR \
            --nThreads 0 \
            --process HLTX \
            --inputCommands='keep *, drop *_hlt*_*_HLT, drop triggerTriggerFilterObjectWithRefs_l1t*_*_HLT'

### Apply clustering algorithm

    python3 ClusteringByDensity.py --file step2_L1P2GT_HLT_DQM_VALIDATION_NANO_inNANOAODSIM.root