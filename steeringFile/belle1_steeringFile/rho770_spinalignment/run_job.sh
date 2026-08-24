#! /bin/bash
##!/bin/sh 
source /sw/belle/local/etc/bashrc_general
export BASF_USER_IF=basfsh.so
export BASF_USER_INIT=user_init.so
#export BASF_NPROCESS=0
#export BELLE_MESSAGE_LEVEL=DDEBUG
#export BELLE_MESSAGE_LEVEL=INFO
#####################
# Belle Environment #
# #####################
#export BELLE_LEVEL=b20090127_0910
#export BELLE_DEBUG=opt
#export BELLE_MSG_MAX_SHLVL=1
### use i386 version rather than RHEL4 64bit version for exp5x MC prod
#export USE_I386_EXP5X=
### avoid glibc detected error
#export MALLOC_CHECK_=0
##
#if [ -f /sw/belle/local/etc/bashrc_general ]; then
#    . /sw/belle/local/etc/bashrc_general
#fi

#####################
# Belle Environment #
#####################
#setenv BELLE_LEVEL b20090127_0910
#setenv BELLE_DEBUG opt
#setenv BELLE_MSG_MAX_SHLVL 1
## use i386 version rather than RHEL4 64bit version for exp5x MC prod
#setenv USE_I386_EXP5X
## avoid glibc detected error
#setenv MALLOC_CHECK_ 0

#
# *** by default, ROOT version 5.22 is used. If you want to use
# *** different version, please uncomment and modify below
#
# setenv ROOTSYS /sw/belle/cern/root_v5.32.02

#if ( -f /sw/belle/local/etc/cshrc_general ) then
#    source /sw/belle/local/etc/cshrc_general
#endif

basf << EOF #>& ./exp71_rs2249_re2348_evtgen-uds_0.log

module register fix_mdst SpinAlignment
path create main
#path create Analysis
path add_module main fix_mdst SpinAlignment

path add_condition main >:0:EXIT
path add_condition main =<:0:KILL

#data type: MC=1 for simulation,MC=0 for experiment.

module put_parameter SpinAlignment  output_filename\./test.root
#module put_parameter SpinAlignment  output_filename\/gpfs/home/belle2/wangz/Work/SpinAlignment/offline/draw/thrust_test/rootFiles/no_kl.root

#module put_parameter SpinAlignment isMCSample\1
#module put_parameter SpinAlignment rmMCTree\0
##mother particle: Y1S/Y2S/Y3S/Y4S/Y4Scon(40)/Y5S
#module put_parameter  SpinAlignment YnsNo\40

initialize
#histogram define exp71_rs2249_re2348_evtgen-uds_0.hbk

#MC or experiment data for procession
#process_event /group/belle/bdata_b/dstprod/dat/e000069/HadronBJ/0127/continuum/08/HadronBJ-e000069r000823-b20090127_0910.mdst 1000
process_event /group/belle/bdata_b/mcprod/dat/e000065/evtgen/uds/00/all/0127/continuum/06/evtgen-uds-00-all-e000065r000626-b20090127_0910.mdst 1000
#process_url http://bweb3/montecarlo.php?ex=7&rs=1490&re=1690&ty=evtgen-uds&dt=continuum&bl=caseB&st=20 


terminate
EOF
###end of script
