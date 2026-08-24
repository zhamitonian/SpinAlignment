#!/usr/bin/env basf2
# -*- coding: utf-8 -*-

############################################
# b2bii hadronic selection
############################################

# Import and mdst loading
import basf2 as b2
from basf2 import LogLevel
import modularAnalysis as ma
from variables import variables
from b2biiConversion import convertBelleMdstToBelleIIMdst
import os
import sys
import ROOT
from array import array

main = b2.create_path()
b2.set_log_level(LogLevel.DEBUG)

os.environ["USE_GRAND_REPROCESS_DATA"] = "1"
os.environ["PGUSER"] = "g0db"

b2.conditions.disable_globaltag_replay()
b2.conditions.globaltags=['B2BII','BellePID',
                              'b2bii_beamParameters_with_smearing']
b2.conditions.prepend_globaltag('Legacy_CollisionAxisCMS_Belle')
b2.conditions.prepend_globaltag('analysis_b2bii')

input_mdst = "/group/belle/bdata_b/mcprod/dat/e000065/evtgen/uds/00/all/0127/continuum/06/evtgen-uds-00-all-e000065r000626-b20090127_0910.mdst"
output_root = "test.root"
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
variables.addAlias('atcKPi', 'atcPIDBelle(3,2)')
track_cut = "cosTheta > -0.511 && cosTheta < 0.842 && atcKPi < 0.4"
ma.fillParticleList("pi+:track", track_cut, path = main)
ma.fillParticleList("pi-:track", track_cut, path = main)


class PionTracksToSingleTree(b2.Module):
    def __init__(self, file_name: str, tree_name: str = 'track'):
        super().__init__()
        self._file_name = file_name
        self._tree_name = tree_name

        self._root_file = None
        self._tree = None

        self.expNo = array('i', [0])
        self.runNo = array('i', [0])
        self.evtNo = array('i', [0])

        Vfloat = ROOT.std.vector('float')
        self.pip_px = Vfloat()
        self.pip_py = Vfloat()
        self.pip_pz = Vfloat()
        self.pip_p = Vfloat()
        self.pip_cosTheta = Vfloat()
        self.pip_atcKPi = Vfloat()

        self.pim_px = Vfloat()
        self.pim_py = Vfloat()
        self.pim_pz = Vfloat()
        self.pim_p = Vfloat()
        self.pim_cosTheta = Vfloat()
        self.pim_atcKPi = Vfloat()

    @staticmethod
    def _store_obj(type_name: str, store_name: str | None = None):
        from ROOT import Belle2
        if store_name is None:
            try:
                return Belle2.PyStoreObj(type_name)
            except Exception:
                return None
        try:
            return Belle2.PyStoreObj(type_name, store_name)
        except Exception:
            return None

    @staticmethod
    def _iter_particles(list_name: str):
        plist = PionTracksToSingleTree._store_obj('ParticleList', list_name)
        if plist is None or not plist.isValid():
            return []
        pl = plist.obj()
        try:
            return list(pl)
        except TypeError:
            pass
        try:
            return [pl.getParticle(i) for i in range(pl.getListSize())]
        except Exception:
            return []

    def initialize(self):
        self._root_file = ROOT.TFile.Open(self._file_name, 'RECREATE')
        self._tree = ROOT.TTree(self._tree_name, self._tree_name)

        self._tree.Branch('expNo', self.expNo, 'expNo/I')
        self._tree.Branch('runNo', self.runNo, 'runNo/I')
        self._tree.Branch('evtNo', self.evtNo, 'evtNo/I')

        self._tree.Branch('pip_px', self.pip_px)
        self._tree.Branch('pip_py', self.pip_py)
        self._tree.Branch('pip_pz', self.pip_pz)
        self._tree.Branch('pip_p', self.pip_p)
        self._tree.Branch('pip_cosTheta', self.pip_cosTheta)
        self._tree.Branch('pip_atcKPi', self.pip_atcKPi)

        self._tree.Branch('pim_px', self.pim_px)
        self._tree.Branch('pim_py', self.pim_py)
        self._tree.Branch('pim_pz', self.pim_pz)
        self._tree.Branch('pim_p', self.pim_p)
        self._tree.Branch('pim_cosTheta', self.pim_cosTheta)
        self._tree.Branch('pim_atcKPi', self.pim_atcKPi)

    def event(self):
        self.pip_px.clear()
        self.pip_py.clear()
        self.pip_pz.clear()
        self.pip_p.clear()
        self.pip_cosTheta.clear()
        self.pip_atcKPi.clear()

        self.pim_px.clear()
        self.pim_py.clear()
        self.pim_pz.clear()
        self.pim_p.clear()
        self.pim_cosTheta.clear()
        self.pim_atcKPi.clear()

        evt_meta = self._store_obj('EventMetaData')
        if evt_meta is not None and evt_meta.isValid():
            emd = evt_meta.obj()
            self.expNo[0] = int(emd.getExperiment())
            self.runNo[0] = int(emd.getRun())
            self.evtNo[0] = int(emd.getEvent())
        else:
            self.expNo[0] = 0
            self.runNo[0] = 0
            self.evtNo[0] = 0

        for p in self._iter_particles('pi+:track'):
            self.pip_px.push_back(float(variables.evaluate('px', p)))
            self.pip_py.push_back(float(variables.evaluate('py', p)))
            self.pip_pz.push_back(float(variables.evaluate('pz', p)))
            self.pip_p.push_back(float(variables.evaluate('p', p)))
            self.pip_cosTheta.push_back(float(variables.evaluate('cosTheta', p)))
            self.pip_atcKPi.push_back(float(variables.evaluate('atcKPi', p)))

        for p in self._iter_particles('pi-:track'):
            self.pim_px.push_back(float(variables.evaluate('px', p)))
            self.pim_py.push_back(float(variables.evaluate('py', p)))
            self.pim_pz.push_back(float(variables.evaluate('pz', p)))
            self.pim_p.push_back(float(variables.evaluate('p', p)))
            self.pim_cosTheta.push_back(float(variables.evaluate('cosTheta', p)))
            self.pim_atcKPi.push_back(float(variables.evaluate('atcKPi', p)))

        self._tree.Fill()

    def terminate(self):
        if self._root_file is None:
            return
        self._root_file.cd()
        if self._tree is not None:
            self._tree.Write()
        self._root_file.Close()

main.add_module(PionTracksToSingleTree(output_root, tree_name='track'))

b2.process(main)
print(b2.statistics)
