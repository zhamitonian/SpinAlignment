#!/bin/bash
source /sw/belle/local/etc/bashrc_general
export BASF_USER_IF=basfsh.so
export BASF_USER_INIT=user_init.so

basf << EOF

module register fix_mdst SpinAlignment
path create main
path create Analysis

path add_module main fix_mdst

path add_condition main >:0:Analysis
path add_condition main =<:0:KILL
path add_module Analysis SpinAlignment

module put_parameter SpinAlignment output_filename\./test.root
module put_parameter SpinAlignment isMCSample\0

initialize
process_url http://bweb3/mdst.php?ex=73&rs=448&re=449&skm=HadronBorJ&dt=continuum&bl=caseB

terminate
EOF
exit
