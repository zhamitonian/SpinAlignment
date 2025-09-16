#! /bin/tcsh -f 
# Arguments
# 1: absolute path to evtgsim mdst files
# 2: file name

if( $3 == '' ) then
  echo "Usage : ./autorunMC.csh [PathToEvtgsimMdstFiles][log and hbk's filename][YnsNo]"
  exit
endif
set evtgsimmdstpath=$1

if ( $2 == '') then 
  set filename=$1
else
  set filename=$2
endif
set YnsNo=$3

foreach file (`\ls -1 $evtgsimmdstpath/evtgen_exp*.mdst `)
   set expname=`echo $file |awk -F '/' '{print $2}'| sed 's/\.mdst//g'`
   echo $expname 
   bsub -q s ./runmc.csh $expname $filename $YnsNo
end
