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

module put_parameter HadronB_dataMC_comp output_filename\./HadronBorJ_e43r592_604.root
module put_parameter HadronB_dataMC_comp isMCSample\0

initialize
process_url http://bweb3/mdst.php?ex=43&rs=592&re=604&skm=HadronBorJ&dt=continuum&bl=caseB

terminate
EOF
exit
