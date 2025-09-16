#! /bin/bash -f

#argument: $0 [expNo] [runStart] [runEND] [eventType] [dataType] [runsPerJob] [YnsNo] [stNo] 
let stNo=0
# streamNo:[0,5]
while [ $stNo -le 0 ]
do
  ### for 5S_onresonance
  ./autorunGenericMC.csh 43 1013 1034 evtgen-charm   5S_onresonance 100 5 $stNo
  ./autorunGenericMC.csh 43 1013 1034 evtgen-uds     5S_onresonance 100 5 $stNo
  ./autorunGenericMC.csh 43 1013 1034 evtgen-bsbs    5S_onresonance 100 5 $stNo
  ./autorunGenericMC.csh 43 1013 1034 evtgen-nonbsbs 5S_onresonance 100 5 $stNo

  ./autorunGenericMC.csh 53 1    272  evtgen-charm   5S_onresonance 50  5 $stNo
  ./autorunGenericMC.csh 53 1    272  evtgen-uds     5S_onresonance 50  5 $stNo
  ./autorunGenericMC.csh 53 1    272  evtgen-bsbs    5S_onresonance 50  5 $stNo
  ./autorunGenericMC.csh 53 1    272  evtgen-nonbsbs 5S_onresonance 50  5 $stNo

  ./autorunGenericMC.csh 67 98   696  evtgen-charm   5S_onresonance 100 5 $stNo
  ./autorunGenericMC.csh 67 98   696  evtgen-uds     5S_onresonance 100 5 $stNo
  ./autorunGenericMC.csh 67 98   696  evtgen-bsbs    5S_onresonance 100 5 $stNo
  ./autorunGenericMC.csh 67 98   696  evtgen-nonbsbs 5S_onresonance 100 5 $stNo

  ./autorunGenericMC.csh 69 1    1309 evtgen-charm   5S_onresonance 103 5 $stNo
  ./autorunGenericMC.csh 69 1    1309 evtgen-uds     5S_onresonance 103 5 $stNo
  ./autorunGenericMC.csh 69 1    1309 evtgen-bsbs    5S_onresonance 103 5 $stNo
  ./autorunGenericMC.csh 69 1    1309 evtgen-nonbsbs 5S_onresonance 103 5 $stNo

  #./autorunGenericMC.csh 71 27   221  evtgen-charm   5S_onresonance 300 5 $stNo
  #./autorunGenericMC.csh 71 27   221  evtgen-uds     5S_onresonance 300 5 $stNo
  #./autorunGenericMC.csh 71 27   221  evtgen-bsbs    5S_onresonance 300 5 $stNo
  #./autorunGenericMC.csh 71 27   221  evtgen-nonbsbs 5S_onresonance 300 5 $stNo

  #./autorunGenericMC.csh 71 2001 2244 evtgen-charm   5S_onresonance 300 5 $stNo
  #./autorunGenericMC.csh 71 2001 2244 evtgen-uds     5S_onresonance 300 5 $stNo
  #./autorunGenericMC.csh 71 2001 2244 evtgen-bsbs    5S_onresonance 300 5 $stNo
  #./autorunGenericMC.csh 71 2001 2244 evtgen-nonbsbs 5S_onresonance 300 5 $stNo

  #### for 5S_continuum
  #./autorunGenericMC.csh 43 1002 1012 evtgen-charm 5S_scan 100 5 $stNo
  #./autorunGenericMC.csh 43 1002 1012 evtgen-uds   5S_scan 100 5 $stNo
  #./autorunGenericMC.csh 61 1210 1373 evtgen-charm 5S_scan 200 5 $stNo
  #./autorunGenericMC.csh 61 1210 1373 evtgen-uds   5S_scan 200 5 $stNo
  #./autorunGenericMC.csh 73 22   916  evtgen-charm 5S_scan 500 5 $stNo
  #./autorunGenericMC.csh 73 22   916  evtgen-uds   5S_scan 500 5 $stNo

  let stNo+=1
done
