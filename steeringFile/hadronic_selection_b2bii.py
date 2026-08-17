#!/usr/bin/env basf2
# -*- coding: utf-8 -*-

############################################
# b2bii hadronic selection
# version : v1.0.1
# Author  : Zheng Wang
# Date    : 2026.08.17
############################################

# Import and mdst loading
import basf2 as b2
from basf2 import LogLevel
import modularAnalysis as ma
from variables import variables
from b2biiConversion import convertBelleMdstToBelleIIMdst
import os
import sys

main = b2.create_path()
b2.set_log_level(LogLevel.DEBUG)

os.environ["USE_GRAND_REPROCESS_DATA"] = "1"
os.environ["PGUSER"] = "g0db"

b2.conditions.disable_globaltag_replay()
b2.conditions.globaltags=['B2BII','BellePID',
                              'b2bii_beamParameters_with_smearing']
b2.conditions.prepend_globaltag('Legacy_CollisionAxisCMS_Belle')
b2.conditions.prepend_globaltag('analysis_b2bii')

input_mdst = "/group/belle/bdata_b/dstprod/dat/e000069/HadronBJ/0127/continuum/08/HadronBJ-e000069r000823-b20090127_0910.mdst" 
output_root = "hadronic_b2bii.root"
convertBelleMdstToBelleIIMdst(
                input_mdst, 
                enableNisKsFinder=False, 
                enableEvtcls=True, 
                HadronA=True, 
                HadronB=True, 
                applySkim=True,
                path=main
            )

## track list
track_cut = "pt > 0.1 and abs(dr) < 2 and abs(dz) < 4"
ma.fillParticleList("pi+:track", track_cut, path = main)
## photon list
photon_cut = "E > 0.1 and clusterBelleQuality == 0"
ma.cutAndCopyList('gamma:photon', 'gamma:mdst', photon_cut, path = main) # gamma:mdst <=> Mdst_gamma_Manager

## event cut
variables.addAlias("Evis_track", "sumValueInList(pi+:track, useCMSFrame(E))")
variables.addAlias("Evis_gamma", "sumValueInList(gamma:photon, useCMSFrame(E))")
variables.addAlias("Evis", "formula(Evis_track + Evis_gamma)") # to make Evis a basf2 event variable

## saving variables same as in Evtcls_hadron_info_Manager
variables.addAlias("Ntrk", "countInList(pi+:track)")
variables.addAlias("Ncls", "countInList(gamma:photon)")
variables.addAlias("Psum", "sumValueInList(pi+:track, useCMSFrame(p))")
variables.addAlias("Pz", "sumValueInList(pi+:track, useCMSFrame(pz))")

ma.buildEventShape(["pi+:track", "gamma:photon"], path=main)
variables.addAlias("heavyJetMass", "formula( max(forwardHemisphereMass, backwardHemisphereMass))")
variables.addAlias("thrustAxisPhi", "atan(thrustAxisY / thrustAxisX)")

variables.addAlias('p_CMS', 'useCMSFrame(p)')
variables.addAlias('theta_CMS', 'useCMSFrame(cosTheta)')
variables.addAlias('phi_CMS', 'useCMSFrame(phi)')

event_vars = ["Ntrk", "Ncls", "Psum", "Evis", "Pz", 
              "thrust", "thrustAxisCosTheta", "thrustAxisPhi", "foxWolframR2", "heavyJetMass"]
kinematic_vars = ["p", "cosTheta", "phi", "p_CMS", "theta_CMS", "phi_CMS"]

main.add_module('VariablesToEventBasedTree', 
                fileName=output_root, 
                particleList='pi+:track', 
                treeName='track', 
                event_variables=event_vars,
                variables = kinematic_vars)

main.add_module('VariablesToEventBasedTree', 
                fileName=output_root, 
                particleList='gamma:photon', 
                treeName='photon', 
                variables = kinematic_vars)

b2.process(main)
print(b2.statistics)


"""
# version log 
## v1.0.0 

## v1.0.1
- use sumValueInList instead of totalEnergyOfParticlesInList to calculate Evis, making it a basf2 event variable, so that it can be used in event var
- date : 2026-08-17
"""