#!/bin/bash
source /sw/belle/local/etc/bashrc_general
export BASF_USER_IF=basfsh.so
export BASF_USER_INIT=user_init.so

basf << EOF

module register fix_mdst HadronB_dataMC_comp
path create main
path create Analysis

path add_module main fix_mdst HadronB_dataMC_comp

path add_condition main >:0:Analysis
path add_condition main =<:0:KILL
path add_module Analysis HadronB_dataMC_comp

module put_parameter HadronB_dataMC_comp output_filename\./HadronBorJ_e35r588_641.root
module put_parameter HadronB_dataMC_comp isMCSample\0

initialize
process_url http://bweb3/mdst.php?ex=35&rs=588&re=1000&skm=HadronBorJ&dt=continuum&bl=caseB
#process_url http://bweb3/montecarlo.php?ex=7&rs=1490&re=1690&ty=evtgen-uds&dt=continuum&bl=caseB&st=20 
#process_event /group/belle/bdata_b/dstprod/dat/e000069/HadronBJ/0127/continuum/08/HadronBJ-e000069r000823-b20090127_0910.mdst 

terminate
EOF
exit
