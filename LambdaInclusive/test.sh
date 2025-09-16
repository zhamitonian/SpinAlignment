#!/bin/sh 

source /sw/belle/local/etc/bashrc_general
export BASF_USER_IF=basfsh.so
export BASF_USER_INIT=user_init.so
export BASF_NPROCESS=0
export BELLE_MESSAGE_LEVEL=DDEBUG

basf << EOF >& test.log

module register fix_mdst 
path create main
path add_module main fix_mdst 
path add_condition main <:0:KILL

initialize

process_event /group/belle/bdata_b/mcprod/dat/e000043/evtgen/bsbs/00/all/0127/5S_onresonance/10/evtgen-bsbs-00-all-e000043r001013-b20090127_0910.mdst 0

terminate
EOF

