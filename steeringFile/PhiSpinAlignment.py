#!/usr/bin/env basf2
# -*- coding: utf-8 -*-

############################################
# belle qqbar -> phi + anything
#                    \-> K+ K- phi

# Version : v1.0.0
# Data    : 2026.03.02
# Author  : Zhen Wang

# Equal to belle1 SpinAlignment.cc v2.1.2
# And based on hadronic_selection_b2bii.py v1.0.0
############################################

import basf2 as b2
from basf2 import LogLevel
import modularAnalysis as ma
from variables import variables
from b2biiConversion import convertBelleMdstToBelleIIMdst
import os

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
#variables.addAlias("Evis", "(sumValueInList(pi+:track, useCMSFrame(E)) + sumValueInList(gamma:photon, useCMSFrame(E)))") # two Evis definitions has some difference ???
variables.addAlias('Evis', 'formula( useCMSFrame(totalEnergyOfParticlesInList(pi+:track)) + useCMSFrame(totalEnergyOfParticlesInList(gamma:photon)) )')
#ma.applyEventCuts("Evis > 7.0", path = main) # don't know why could not apply event cut directly here
ma.applyCuts("gamma:photon", "Evis > 7.0", path = main)
ma.applyCuts("pi+:track", "Evis > 7.0", path = main)

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
