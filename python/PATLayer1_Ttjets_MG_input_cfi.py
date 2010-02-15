import FWCore.ParameterSet.Config as cms

readFiles = cms.untracked.vstring(
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_100.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_101.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_102.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_103.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_104.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_105.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_106.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_107.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_108.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_109.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_10.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_110.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_111.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_112.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_113.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_114.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_115.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_116.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_118.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_119.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_11.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_120.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_121.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_122.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_123.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_124.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_125.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_126.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_127.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_128.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_129.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_12.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_130.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_131.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_132.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_133.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_134.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_135.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_136.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_137.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_138.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_139.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_13.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_140.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_141.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_142.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_143.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_144.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_145.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_147.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_148.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_149.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_14.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_150.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_151.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_152.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_153.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_154.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_155.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_156.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_157.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_15.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_160.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_161.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_162.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_163.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_164.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_165.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_166.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_167.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_169.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_170.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_171.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_172.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_173.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_174.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_175.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_176.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_177.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_178.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_179.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_17.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_180.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_181.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_182.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_183.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_184.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_185.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_186.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_187.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_188.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_189.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_18.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_190.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_191.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_192.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_193.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_194.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_195.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_196.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_197.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_198.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_199.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_19.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_1.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_200.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_201.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_202.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_203.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_204.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_205.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_206.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_207.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_208.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_209.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_20.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_210.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_21.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_22.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_23.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_24.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_25.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_26.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_27.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_28.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_29.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_2.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_30.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_31.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_32.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_33.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_34.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_35.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_36.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_37.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_38.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_39.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_3.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_40.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_41.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_42.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_43.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_45.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_46.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_47.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_48.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_49.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_4.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_50.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_51.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_52.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_53.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_54.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_55.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_56.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_57.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_58.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_59.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_60.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_61.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_62.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_63.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_64.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_66.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_67.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_69.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_6.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_70.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_71.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_72.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_73.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_74.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_76.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_77.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_78.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_79.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_7.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_80.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_81.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_83.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_84.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_86.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_87.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_88.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_89.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_8.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_90.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_91.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_92.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_93.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_94.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_95.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_96.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_98.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_99.root',
'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/pvmulder/CMSSW223/Common/PATLayer1/Ttjets-madgraph/PATLayer1_9.root'
)
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",duplicateCheckMode = cms.untracked.string('noDuplicateCheck') ,fileNames = readFiles#,
    #inputCommands = cms.untracked.vstring(
    #    "drop patJets_selectedPatJets__SKIMsemilep"
    #)
)
