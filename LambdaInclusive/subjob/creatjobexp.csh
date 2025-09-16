#! /bin/tcsh -f

# Arguments
if( $7 == '' ) then
  echo "Usage : $0 [expNo] [runSTART] [runEND] [skimType] [dataType] [blType] [YnsNo]"
  exit
endif

set expNo = $1
set runSTART = $2
set runEND = $3
set skimType = $4
set dataType = $5
set blType = $6
set YnsNo = $7

echo '#! /bin/bash'
echo 'source /sw/belle/local/etc/bashrc_general'
echo 'export BASF_USER_IF=basfsh.so'
echo 'export BASF_USER_INIT=user_init.so'
echo '#export BASF_NPROCESS=0'
echo '#export BELLE_MESSAGE_LEVEL=DDEBUG'
echo '#export BELLE_MESSAGE_LEVEL=INFO'
echo '#export BELLE_MSG_MAX_SHLVL=1'
echo ''
#echo 'set data_type=$1'
#echo 'set expno=$2'
#echo ''
set fileDir = /group/belle/users/guanyh/LambdaAna
#set fileDir = /home/belle/guanyh/belle_ana/LambdaAna 
set fileName = exp${expNo}_rs${runSTART}_re${runEND}_${skimType}
set RootName = exp${expNo}_rs${runSTART}_re${runEND}_${skimType}_tree.root
echo 'basf << EOF >& /home/belle/guanyh/belle_ana/LambdaAna/subjob/log/'${dataType}'/data/'${fileName}'.log'
echo ''
#echo 'nprocess set 2'
echo ''
echo 'module register fix_mdst LambdaAna'
echo 'path create main'
echo 'path create Analysis'
echo ''
echo 'path add_module main fix_mdst'
echo ''
echo 'path add_condition main >:0:Analysis'
echo 'path add_condition main =<:0:KILL'
echo 'path add_module Analysis LambdaAna'
echo ''
echo '#data type: MC=1 for simulation,MC=0 for experiment.'
echo '#module put_parameter ISRetaomega  mc\\0'
echo 'module put_parameter LambdaAna output_filename\\'${fileDir}'/hbk/'${dataType}'/data_thrustSmaler0.8/'${RootName}''
echo 'module put_parameter LambdaAna isMCSample\\0'
echo '##mother particle: Y1S/Y2S/Y3S/Y4S/Y4Scon(40)/Y5S'
echo '#module put_parameter  ISRetaomega  YnsNo\\'${YnsNo}
echo ''
echo 'initialize'
#echo 'histogram define ./hbk/'${dataType}'/data/'${fileName}'.hbk'
#echo 'histogram define /ghi/fs01/belle/bdata2/users/guanyh/LambdaAna/hbk/'${dataType}'/data/'${fileName}'.hbk'
echo 'histogram define '${fileDir}'/hbk/'${dataType}'/data_thrustSmaler0.8/'${fileName}'.hbk'
echo ''

set http = "http://bweb3/mdst.php"
set cmnd = "ex="${expNo}"&rs="${runSTART}"&re="${runEND}"&skm="${skimType}"&dt="${dataType}"&bl="${blType}
echo '#MC or experiment data for procession'
echo "process_url "$http"?"$cmnd

echo ''
echo 'terminate'
echo 'EOF'
echo 'exit'
echo '###end of script'

exit
