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

module put_parameter HadronB_dataMC_comp output_filename\./test.root
module put_parameter HadronB_dataMC_comp isMCSample\1
module put_parameter HadronB_dataMC_comp rmMCTree\1
module put_parameter HadronB_dataMC_comp rmTree\0
module put_parameter HadronB_dataMC_comp rmMixTree\1

initialize
process_url http://bweb3/montecarlo.php?ex=65&rs=1008&re=1009&ty=evtgen-charm&dt=1S_scan&bl=caseB&st=0

terminate
EOF
exit
