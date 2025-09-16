#! /bin/bash -f

#argument: $0 [expNo] [runStart] [runEND] [eventType] [dataType] [runsPerJob] [YnsNo] [stNo] 
let stNo=0
# streamNo region:[0,5]
while [ $stNo -le 0 ]
do
  ## Y(1S)
  ./autorunExpdata.csh 65 1008 1232 evtgen-charm 1S_scan  60 1 $stNo
  ./autorunExpdata.csh 65 1008 1232 evtgen-uds   1S_scan  60 1 $stNo
  
  ## Y(2S)
  ./autorunExpdata.csh 67 1002 1123 evtgen-charm 2S_scan  40 2 $stNo
  ./autorunExpdata.csh 67 1002 1123 evtgen-uds   2S_scan  40 2 $stNo
  ./autorunExpdata.csh 71 303  696  evtgen-charm 2S_scan  70 2 $stNo
  ./autorunExpdata.csh 71 303  696  evtgen-uds   2S_scan  70 2 $stNo
  
  ## Y(3S)
  ./autorunExpdata.csh 49 1001 1227 evtgen-charm 3S_scan 120 3 $stNo
  ./autorunExpdata.csh 49 1001 1227 evtgen-uds   3S_scan 120 3 $stNo

  let stNo+=1
done
