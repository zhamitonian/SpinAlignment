#!/bin/sh

MODULE=LambdaAna
#MODULE=userana
#HBKFILE=/ghi/fs01/belle/bdata2/users/guanyh/LambdaAna/test.hbk
HBKFILE=MCtest.hbk
# OUTFILE=test.dat
#INFILE=bfsa01:/bdata/dstprod/dat/e000007/HadronB/0416/on_resonance/01/e000007r000100-b20020416_1604.mdst 
#INFILE=/group/belle/bdata_b/dstprod/dat/e000069/HadronBJ/0127/5S_scan/00/HadronBJ-e000069r000012-b20090127_0910.mdst
# INFILE=bfsa07:/bdata/dstprod/dat/e000031/HadronBJ/0429/on_resonance/01/HadronBJ-e000031r000152-b20040429_2253.mdst
# INFILE=bfsa02:/bdata/dstprod/dat/e000051/HadronBJ/0307/on_resonance/01/HadronBJ-e000051r000199-b20070307_1108.mdst
#INFILE=/group/belle/bdata_b/dstprod/dat/e000069/HadronBJ/0127/continuum/08/HadronBJ-e000069r000823-b20090127_0910.mdst
#INFILE=/group/belle/bdata_b/mcprod/dat/e000065/evtgen/charm/00/all/0127/continuum/06/evtgen-charm-00-all-e000065r000626-b20090127_0910.mdst
INFILE=/group/belle/bdata_b/mcprod/dat/e000065/evtgen/uds/00/all/0127/continuum/06/evtgen-uds-00-all-e000065r000662-b20090127_0910.mdst
#INFILE=/group/belle/bdata_b/dstprod/dat/e000061/HadronBJ/1107/on_resonance/01/HadronBJ-e000061r000100-b20081107_1418.mdst

#NEVENT=10000000
NEVENT=10000

basf <<EOF
path create main
path create analysis
path add_module main fix_mdst
path add_condition main >:0:analysis
path add_condition main <=:0:KILL
path add_module analysis ${MODULE}
module put_parameter LambdaAna output_filename\RootMC_tree.root
#module put_parameter LambdaAna output_filename_mcTree\RootMC_tree_mcTree.root
module put_parameter LambdaAna isMCSample\1
module put_parameter LambdaAna rmMCTree\0
module put_parameter LambdaAna rmTree\0
module put_parameter LambdaAna rmMixTree\0

initialize
histogram define ${HBKFILE}
process_event ${INFILE} ${NEVENT}
#process_url http://bweb3/montecarlo.php?ex=69&rs=883&re=902&ty=evtgen-uds&dt=continuum&bl=caseB&st=0
output close
terminate
EOF
