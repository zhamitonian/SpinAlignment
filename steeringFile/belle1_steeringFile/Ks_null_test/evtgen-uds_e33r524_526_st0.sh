#!/bin/bash
source /sw/belle/local/etc/bashrc_general
export BASF_USER_IF=basfsh.so
export BASF_USER_INIT=user_init.so

basf << EOF

module register fix_mdst KsSpinAlignment_NullTest
path create main
path create Analysis

path add_module main fix_mdst KsSpinAlignment_NullTest

path add_condition main >:0:Analysis
path add_condition main =<:0:KILL
path add_module Analysis KsSpinAlignment_NullTest

module put_parameter KsSpinAlignment_NullTest output_filename\./test.root
module put_parameter KsSpinAlignment_NullTest isMCSample\1
module put_parameter KsSpinAlignment_NullTest rmMCTree\1
module put_parameter KsSpinAlignment_NullTest rmTree\0
module put_parameter KsSpinAlignment_NullTest rmMixTree\1

initialize
process_url http://bweb3/montecarlo.php?ex=33&rs=524&re=526&ty=evtgen-uds&dt=continuum&bl=caseB&st=0 

terminate
EOF
exit
