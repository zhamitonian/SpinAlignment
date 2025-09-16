############################################
# Analysis for Inclusive gadron production
############################################

# Import and mdst loading
import basf2 as b2
from basf2 import LogLevel
import modularAnalysis as ma
from stdCharged import stdPi, stdK
from stdPhotons import stdPhotons
from variables import variables
import variables.collections as vc
import variables.utils as vu

import sys

from ROOT import Belle2

outputFile = "analysis.root"

#b2.conditions.reset()
#### necessary for the 5Sscan datasets on the top of analysis GT
b2.conditions.prepend_globaltag("data_beam_conditions_proc13prompt")
b2.conditions.prepend_globaltag(ma.getAnalysisGlobaltag())
b2.set_log_level(b2.LogLevel.ERROR)

mymain = b2.Path()

# load input root file
ma.inputMdst(environmentType='default', filename='', path=mymain)

# define several aliases
variables.addAlias('MissingMomCMS', 'missingMomentumOfEventCMS')
variables.addAlias('CosHelicity', 'cosHelicityAngleMomentum')
variables.addAlias('CMS_p', 'useCMSFrame(p)')

#skim
variables.addAlias('SkimHad', 'SoftwareTriggerResult(software_trigger_cut&skim&accept_hadron)')

######################################
# information at detector level
######################################

# from BELLE2-NOTE-PH-2019-025
# FEI: analysis/examples/GraFEI/graFEI_ntupliser_Upsilon.py L95
# good track 
// check Belle2 (internal wiki -> Physics Performance Webhome)
track_cut = 'pt > 0.1 and abs(dr) < 1 and abs(dz) < 3'
ma.fillParticleList('pi+:track', track_cut, path=mymain)

# good cluster & photon
// 添加 fakephoton supperssion  backgroundphoton supperssion变量
cluster_cut = 'E > 0.1 and abs(clusterTiming) < 200'
photon_cut = 'E > 0.1 and abs(clusterTiming) < 200 and thetaInCDCAcceptance'
ma.fillParticleList('gamma:cluster', cluster_cut, path=mymain)
ma.fillParticleList('gamma:photon', photon_cut, path=mymain)

# Apply event based cuts
nTrack_cut = "nCleanedTracks(" + track_cut + ") >= 3"
ma.applyEventCuts(f"[{nTrack_cut}]", path=mymain)

# energy of tracks and clusters at CMS, EnergyCMS = visibleEnergyOfEventCMS @ Belle II
// check 下面两个visible E 计算结果是不是一样的， 并且和Belle II官方的  visibleEnergyOfEventCMS 计算进行比较 
//  检查程序有无考虑 中子质量     
variables.addAlias('EnergyCMS', 'formula( useCMSFrame(totalEnergyOfParticlesInList(pi+:track)) + useCMSFrame(totalEnergyOfParticlesInList(gamma:cluster)) )')
# visible energy: momenta of tracks & energy of photon at CMS
variables.addAlias('EvisCMS', 'formula( useCMSFrame(sumValueInList(pi+:track, p)) + useCMSFrame(sumValueInList(gamma:photon, E)) )')

# momentum balance Pz of tracks & photons at CMS
# BalancePzCMS vs. missingMomentumOfEventCMS_Pz @ Belle II: sometime =, sometime !=
# missingMomentumOfEventCMS_Pz @ buildEventKinematics: using good track & good cluster
# BalancePzCMS using good track & good photon [!= good cluster]
variables.addAlias('BalancePzCMS', 'formula( useCMSFrame(totalPzOfParticlesInList(pi+:track))+ useCMSFrame(totalPzOfParticlesInList(gamma:photon)) )')

# calorimeter energy sum in ECL
# energy of photon & neutral hadron in ECAL, energy of tracked tracks in ECAL
# analysis/variables/src/MetaVariables.cc L2574
+visible in ECL
variables.addAlias('ECLEnergy', 'formula( totalEnergyOfParticlesInList(gamma:photon) + totalECLEnergyOfParticlesInList(pi+:track) )')
variables.addAlias('ECLEnergyWO', 'formula( totalEnergyOfParticlesInList(gamma:cluster) + totalECLEnergyOfParticlesInList(pi+:track) )')

# averafe vz of tracks for a event
ma.fillParticleList('pi+:forvz', 'pt > 0.1 and dr < 1', path=mymain)
variables.addAlias('AverageVz', 'formula( averageValueInList(pi+:forvz, dz) )')

# calculate event kinematics variables 
ma.buildEventKinematics(['pi+:track', 'gamma:cluster'], path=mymain)

# Builds the event shape enabling explicitly ALL the variables.
# Most of them are actually enabled by default, but here we prefer
# to list explicitly all the flags
ma.buildEventShape(['pi+:track', 'gamma:cluster'], path=mymain)

// formula?
variables.addAlias('HeavyJetMass', 'formula( max(forwardHemisphereMass, backwardHemisphereMass) )')
#variables.addAlias('HeavyJetEnergy', 'conditionalVariableSelector(forwardHemisphereMass > backwardHemisphereMass, forwardHemisphereEnergy, backwardHemisphereEnergy)')
variables.addAlias('HeavyJetE', 'conditionalVariableSelector(forwardHemisphereMass > backwardHemisphereMass, forwardHemisphereEnergy, backwardHemisphereEnergy)')
ma.variablesToEventExtraInfo('pi+:track', {'HeavyJetE':'HeavyJetE'}, option=2, path=mymain)
variables.addAlias('HeavyJetEnergy', 'eventExtraInfo(HeavyJetE)')

# creates a list of charged kaon & pion
kaon_cut = track_cut + ' and ' + 'kaonID > 0.6'
pion_cut = track_cut + ' and ' + 'pionID > 0.6'
ma.fillParticleList('K+:spin', kaon_cut, path=mymain)
ma.fillParticleList('pi+:spin', pion_cut, path=mymain)

# reconstructs phi -> K+ K-:
ma.reconstructDecay('phi:spin -> K+:spin K-:spin', '0.98 < M < 1.08', path=mymain)

ma.variablesToEventExtraInfo('pi+:track', {'nParticlesInList(pi+:track)':'NosTrack'}, option=2, path=mymain)
variables.addAlias('NosTrack', 'eventExtraInfo(NosTrack)')

ma.variablesToEventExtraInfo('gamma:cluster', {'nParticlesInList(gamma:cluster)':'NosCluster'}, option=2, path=mymain)
variables.addAlias('NosCluster', 'eventExtraInfo(NosCluster)')

ma.variablesToEventExtraInfo('gamma:photon', {'nParticlesInList(gamma:photon)':'NosPhoton'}, option=2, path=mymain)
variables.addAlias('NosPhoton', 'eventExtraInfo(NosPhoton)')

ma.variablesToEventExtraInfo('phi:spin', {'nParticlesInList(phi:spin)':'NosPhi'}, option=2, path=mymain)
variables.addAlias('NosPhi', 'ifNANgiveX(eventExtraInfo(NosPhi), 0)')

# reconstructs K*(892)-> K+ pi-:
ma.reconstructDecay('K*0:spin -> K+:spin pi-:spin', '0.6 < M < 1.2', path=mymain)

ma.variablesToEventExtraInfo('K*0:spin', {'nParticlesInList(K*0:spin)':'NosKstar'}, option=2, path=mymain)
variables.addAlias('NosKstar', 'ifNANgiveX(eventExtraInfo(NosKstar), 0)')

# event shape variables
event_vars  = ['Ecms', 'NosTrack', 'NosCluster', 'NosPhoton', 'SkimHad']
event_vars += ['NosKstar', 'NosPhi', 'MissingMomCMS']
event_vars += ['EnergyCMS', 'EvisCMS', 'BalancePzCMS']
event_vars += ['ECLEnergy', 'ECLEnergyWO', 'AverageVz']
event_vars += ['HeavyJetMass', 'HeavyJetEnergy']
event_vars += ['sphericity', 'aplanarity', 'foxWolframR2']
event_vars += ['thrust', 'thrustAxisCosTheta']
#event_vars += ['visibleEnergyOfEventCMS', 'missingMomentumOfEventCMS_Pz']

is_hadron  = ['Ecms', 'SkimHad']
is_hadron += ['EvisCMS', 'BalancePzCMS', 'HeavyJetMass']
is_hadron += ['ECLEnergyWO', 'AverageVz']

# Select variables that we want to store to ntuple
lab_kinematics = ['px', 'py', 'pz', 'E']
# Kinematic variables in CMS frame
cms_kinematics = vu.create_aliases(lab_kinematics, "useCMSFrame({variable})", "CMS")

part_vars = cms_kinematics + ['p', 'pt', 'theta', 'phi'] + \
            ['kaonID', 'pionID', 'protonID', 'charge', 'nCDCHits'] 
reso_vars = cms_kinematics + ['M', 'CosHelicity', 'CMS_p'] + \
            vu.create_daughter_aliases(part_vars, 0) + \
            vu.create_daughter_aliases(part_vars, 1)

phis_vars = vu.create_aliases_for_selected(list_of_variables=reso_vars, decay_string='^phi', prefix='phi')
kstar_vars = vu.create_aliases_for_selected(list_of_variables=reso_vars, decay_string='^K*0', prefix='kstar')

# Saving variables to ntuple
ma.variablesToNtuple('phi:spin', phis_vars + is_hadron + ['NosPhi'], filename=outputFile, treename='phi', path=mymain)
ma.variablesToNtuple('K*0:spin', kstar_vars + is_hadron + ['NosKstar'], filename=outputFile, treename='kstar', path=mymain)

# Event-wise information for phi list in the event (common)
mymain.add_module('VariablesToEventBasedTree', fileName=outputFile, particleList='pi+:track', treeName='event', event_variables=event_vars)

# Process the events
b2.process(mymain)
# print out the summary
print(b2.statistics)
