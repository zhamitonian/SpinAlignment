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

module put_parameter HadronB_dataMC_comp output_filename\./evtgen-uds_e61r728_733_st2.root
module put_parameter HadronB_dataMC_comp isMCSample\1
module put_parameter HadronB_dataMC_comp rmMCTree\1
module put_parameter HadronB_dataMC_comp rmTree\0
module put_parameter HadronB_dataMC_comp rmMixTree\1

initialize
process_url http://bweb3/montecarlo.php?ex=61&rs=731&re=731&ty=evtgen-uds&dt=continuum&bl=caseB&st=2

terminate
EOF
exit
