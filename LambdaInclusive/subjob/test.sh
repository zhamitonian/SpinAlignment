#! /bin/bash
source /sw/belle/local/etc/bashrc_general
export BASF_USER_IF=basfsh.so
export BASF_USER_INIT=user_init.so
#export BASF_NPROCESS=0
#export BELLE_MESSAGE_LEVEL=DDEBUG
#export BELLE_MESSAGE_LEVEL=INFO
#export BELLE_MSG_MAX_SHLVL=1

basf << EOF >& /home/belle/guanyh/belle_ana/LambdaAna/subjob/log/continuum/genericMC/exp17_rs486_re544_evtgen-charm_21.log


module register fix_mdst LambdaAna
path create main
path create Analysis

path add_module main fix_mdst

path add_condition main >:0:Analysis
path add_condition main =<:0:KILL
path add_module Analysis LambdaAna

#data type: MC=1 for simulation,MC=0 for experiment.
module put_parameter LambdaAna  output_filename\/group/belle/users/guanyh/LambdaAna/hbk/continuum/genericMC_thrustSmaler0.8/exp17_rs486_re544_evtgen-charm_21_tree.root
module put_parameter LambdaAna isMCSample\1
module put_parameter LambdaAna rmMCTree\1
module put_parameter LambdaAna rmTree\0
module put_parameter LambdaAna rmMixTree\1
##mother particle: Y1S/Y2S/Y3S/Y4S/Y4Scon(40)/Y5S
#module put_parameter  LambdaAna YnsNo\40

initialize
histogram define /group/belle/users/guanyh/LambdaAna/hbk/continuum/genericMC_thrustSmaler0.8/exp17_rs486_re544_evtgen-charm_21.hbk

#MC or experiment data for procession
process_url http://bweb3/montecarlo.php?ex=17&rs=486&re=544&ty=evtgen-charm&dt=continuum&bl=100&st=21

terminate
EOF
exit
###end of script
