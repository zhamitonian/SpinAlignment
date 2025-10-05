#!/usr/bin/env basf2
# -*- coding: utf-8 -*-

############################################
# b2bii test for Inclusive hadron production
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

from STEERING_TOOLS import BelleAnalysisBase


b2.set_log_level(LogLevel.WARNING)

steering_tools = BelleAnalysisBase(belle_version='belle1', analysis_mode='data')

main = b2.create_path()

if len(sys.argv) > 1:
    steering_tools.parse_arguments()
    steering_tools.setup_IO(main)
else:
    input_mdst = "/group/belle/bdata_b/dstprod/dat/e000069/HadronBJ/0127/continuum/08/HadronBJ-e000069r000823-b20090127_0910.mdst" 
    output_root = "b2bii_test_thrust_cal.root"
    steering_tools.setup_IO(main, input_mdst, output_root)
    
steering_tools.setup_common_aliases()
steering_tools.setup_environment()

variables.addAlias('MissingMomCMS', 'missingMomentumOfEventCMS')
variables.addAlias('CosHelicity', 'cosHelicityAngleMomentum')
variables.addAlias('SkimHad', 'SoftwareTriggerResult(software_trigger_cut&skim&accept_hadron)')

track_cut = 'pt > 0.1 and abs(dr) < 1 and abs(dz) < 3'
ma.fillParticleList('pi+:track', track_cut, path=main)

# good cluster & photon
cluster_cut = 'E > 0.1'
photon_cut = 'E > 0.1  and thetaInCDCAcceptance'
#ma.fillParticleList('gamma:cluster', cluster_cut, path=main)
#ma.fillParticleList('gamma:photon', photon_cut, path=main)
ma.cutAndCopyList("gamma:photon", "gamma:mdst", photon_cut, path = main)
ma.cutAndCopyList("gamma:cluster", "gamma:mdst", cluster_cut, path = main)

# Apply event based cuts
nTrack_cut = "nCleanedTracks(" + track_cut + ") >= 3"
ma.applyEventCuts(f"[{nTrack_cut}]", path=main)

variables.addAlias('EnergyCMS', 'formula( useCMSFrame(totalEnergyOfParticlesInList(pi+:track)) + useCMSFrame(totalEnergyOfParticlesInList(gamma:cluster)) )')
#variables.addAlias('EvisCMS', 'formula( useCMSFrame(sumValueInList(pi+:track, E)) + useCMSFrame(sumValueInList(gamma:photon, E)) )')
variables.addAlias('EvisCMS', 'formula( useCMSFrame(totalEnergyOfParticlesInList(pi+:track)) + useCMSFrame(totalEnergyOfParticlesInList(gamma:photon)) )')

variables.addAlias('BalancePzCMS', 'formula( useCMSFrame(totalPzOfParticlesInList(pi+:track))+ useCMSFrame(totalPzOfParticlesInList(gamma:photon)) )')

variables.addAlias('ECLEnergy', 'formula( totalEnergyOfParticlesInList(gamma:photon) + totalECLEnergyOfParticlesInList(pi+:track) )')
variables.addAlias('ECLEnergyWO', 'formula( totalEnergyOfParticlesInList(gamma:cluster) + totalECLEnergyOfParticlesInList(pi+:track) )')

# averafe vz of tracks for a event
ma.fillParticleList('pi+:forvz', 'pt > 0.1 and dr < 1', path=main)
variables.addAlias('AverageVz', 'formula( averageValueInList(pi+:forvz, dz) )')

# calculate event kinematics variables 
ma.buildEventKinematics(['pi+:track', 'gamma:cluster'], path=main)

#ma.buildEventShape(['pi+:track', 'gamma:cluster'], path=main)
# test thrust calculation with charged tracks only
ma.buildEventShape(['pi+:track'], path=main)

variables.addAlias('HeavyJetMass', 'formula( max(forwardHemisphereMass, backwardHemisphereMass) )')
variables.addAlias('HeavyJetE', 'conditionalVariableSelector(forwardHemisphereMass > backwardHemisphereMass, forwardHemisphereEnergy, backwardHemisphereEnergy)')
ma.variablesToEventExtraInfo('pi+:track', {'HeavyJetE':'HeavyJetE'}, option=2, path=main)
variables.addAlias('HeavyJetEnergy', 'eventExtraInfo(HeavyJetE)')

# creates a list of charged kaon & pion
kaon_cut = track_cut + ' and ' + 'kaonID > 0.6'
pion_cut = track_cut + ' and ' + 'pionID > 0.6'
ma.fillParticleList('K+:spin', kaon_cut, path=main)
ma.fillParticleList('pi+:spin', pion_cut, path=main)

# reconstructs phi -> K+ K-:
ma.reconstructDecay('phi:spin -> K+:spin K-:spin', '0.98 < M < 1.08', path=main)

ma.variablesToEventExtraInfo('pi+:track', {'nParticlesInList(pi+:track)':'NosTrack'}, option=2, path=main)
variables.addAlias('NosTrack', 'eventExtraInfo(NosTrack)')

ma.variablesToEventExtraInfo('gamma:cluster', {'nParticlesInList(gamma:cluster)':'NosCluster'}, option=2, path=main)
variables.addAlias('NosCluster', 'ifNANgiveX(eventExtraInfo(NosCluster),0)')

ma.variablesToEventExtraInfo('gamma:photon', {'nParticlesInList(gamma:photon)':'NosPhoton'}, option=2, path=main)
variables.addAlias('NosPhoton', 'ifNANgiveX(eventExtraInfo(NosPhoton), 0)')

ma.variablesToEventExtraInfo('phi:spin', {'nParticlesInList(phi:spin)':'NosPhi'}, option=2, path=main)
variables.addAlias('NosPhi', 'ifNANgiveX(eventExtraInfo(NosPhi), 0)')

# reconstructs K*(892)-> K+ pi-:
ma.reconstructDecay('K*0:spin -> K+:spin pi-:spin', '0.6 < M < 1.2', path=main)

ma.variablesToEventExtraInfo('K*0:spin', {'nParticlesInList(K*0:spin)':'NosKstar'}, option=2, path=main)
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

"""
is_hadron  = ['Ecms', 'SkimHad']
is_hadron += ['EvisCMS', 'BalancePzCMS', 'HeavyJetMass']
is_hadron += ['ECLEnergyWO', 'AverageVz']

# Select variables that we want to store to ntuple
lab_kinematics = ['px', 'py', 'pz', 'E']
# Kinematic variables in CMS frame
cms_kinematics = vu.create_aliases(lab_kinematics, "useCMSFrame({variable})", "CMS")

part_vars = cms_kinematics + ['p', 'pt', 'theta', 'phi'] + \
            ['kaonID', 'pionID', 'protonID', 'charge', 'nCDCHits'] 
reso_vars = cms_kinematics + ['M', 'CosHelicity', 'p_CMS'] + \
            vu.create_daughter_aliases(part_vars, 0) + \
            vu.create_daughter_aliases(part_vars, 1)

phis_vars = vu.create_aliases_for_selected(list_of_variables=reso_vars, decay_string='^phi', prefix='phi')
kstar_vars = vu.create_aliases_for_selected(list_of_variables=reso_vars, decay_string='^K*0', prefix='kstar')

# Saving variables to ntuple
ma.variablesToNtuple('phi:spin', 
                     phis_vars + is_hadron + ['NosPhi'], 
                     filename=steering_tools.output_file, 
                     treename='phi', 
                     path=main)
ma.variablesToNtuple('K*0:spin', 
                     kstar_vars + is_hadron + ['NosKstar'], 
                     filename=steering_tools.output_file, 
                     treename='kstar', 
                     path=main)
"""

kinematics = ['p', 'theta', 'phi']

main.add_module('VariablesToEventBasedTree', 
                fileName=steering_tools.output_file, 
                particleList='pi+:track', 
                treeName='event_trk', 
                event_variables=event_vars,
                variables = kinematics)
main.add_module('VariablesToEventBasedTree', 
                fileName=steering_tools.output_file, 
                particleList='gamma:cluster', 
                treeName='event_cls', 
                variables = kinematics)
main.add_module('VariablesToEventBasedTree', 
                fileName=steering_tools.output_file, 
                particleList='gamma:photon', 
                treeName='event_pho', 
                variables = kinematics)


b2.process(main)
print(b2.statistics)
