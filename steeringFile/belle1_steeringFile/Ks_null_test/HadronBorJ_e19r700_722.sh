#!/bin/bash
source /sw/belle/local/etc/bashrc_general
export BASF_USER_IF=basfsh.so
export BASF_USER_INIT=user_init.so

basf << EOF

module register fix_mdst KsSpinAlignment_NullTest
path create main
path create Analysis

path add_module main fix_mdst

path add_condition main >:0:Analysis
path add_condition main =<:0:KILL
path add_module Analysis KsSpinAlignment_NullTest

module put_parameter KsSpinAlignment_NullTest output_filename\./HadronBorJ_e19r700_722.root
module put_parameter KsSpinAlignment_NullTest isMCSample\0

initialize
process_url http://bweb3/mdst.php?ex=19&rs=700&re=700&skm=HadronBorJ&dt=continuum&bl=caseB

terminate
EOF
exit
