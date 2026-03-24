#!/usr/bin/env basf2
# -*- coding: utf-8 -*-

############################################
# b2bii test for hadronic selection
############################################

# Import and mdst loading
import basf2 as b2
from basf2 import LogLevel
import modularAnalysis as ma
from variables import variables
from b2biiConversion import convertBelleMdstToBelleIIMdst
import os

b2.set_log_level(LogLevel.DEBUG)

os.environ["USE_GRAND_REPROCESS_DATA"] = "1"
os.environ["PGUSER"] = "g0db"

b2.conditions.disable_globaltag_replay()
b2.conditions.globaltags=['B2BII','BellePID',
                          'b2bii_beamParameters_with_smearing' ,
                          'Legacy_CollisionAxisCMS_Belle',
                          'analysis_b2bii']

main = b2.create_path()

input_mdst = "/group/belle/bdata_b/dstprod/dat/e000069/HadronBJ/0127/continuum/08/HadronBJ-e000069r000823-b20090127_0910.mdst" 
output_root = "./test.root"
convertBelleMdstToBelleIIMdst(input_mdst, applySkim= True, 
                              generatorLevelReconstruction = False, 
                              path = main ,HadronB =True, 
                              enableEvtcls =True)

## track list
track_cut = "pt > 0.1 and abs(dr) < 2.0 and abs(dz) < 4.0"
ma.fillParticleList("pi+:mdst", track_cut, path = main)
## photon list
photon_cut = "E > 0.1 and clusterBelleQuality == 0"
ma.applyCuts("gamma:mdst", photon_cut, path = main)  # gamma:mdst <=> Mdst_gamma_Manager
#ecl_cut = "17 < theta/acos(-1) * 180 < 150"
#ma.cutAndCopyList("gamma:ecl", "gamma:mdst", ecl_cut, path = main)

## event cut
variables.addAlias("Evis", "(sumValueInList(pi+:mdst, useCMSFrame(E)) + sumValueInList(gamma:mdst, useCMSFrame(E)))")
ma.applyEventCuts("Evis > 7.0", path = main)

## saving variables same as in Evtcls_hadron_info_Manager
variables.addAlias("Ntrk", "countInList(pi+:mdst)")
variables.addAlias("Ncls", "countInList(gamma:mdst)")
variables.addAlias("Psum", "sumValueInList(pi+:mdst, useCMSFrame(p))")
#variables.addAlias("Esum", "sumValueInList(gamma:ecl, useCMSFrame(E))")
variables.addAlias("Pz", "sumValueInList(pi+:mdst, useCMSFrame(pz))")

ma.buildEventShape(["pi+:mdst", "gamma:mdst"], path=main)
variables.addAlias("heavyJetMass", "formula( max(forwardHemisphereMass, backwardHemisphereMass))")
variables.addAlias("thrustAxisPhi", "atan(thrustAxisY / thrustAxisX)")

event_vars = ["Ntrk", "Ncls", "Psum", "Esum", "Evis", "Pz", 
              "thrust", "thrustAxisCosTheta", "thrustAxisPhi", "foxWolframR2", "heavyJetMass"]
kinematic_vars = ["p", "cosTheta", "phi"]

main.add_module('VariablesToEventBasedTree', 
                fileName=output_root, 
                particleList='pi+:mdst', 
                treeName='track', 
                event_variables=event_vars,
                variables = kinematic_vars)

main.add_module('VariablesToEventBasedTree', 
                fileName=output_root, 
                particleList='gamma:mdst', 
                treeName='photon', 
                variables = kinematic_vars)

b2.process(main)
print(b2.statistics)
