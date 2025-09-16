#! /bin/bash -f

#argument: $0 [expNo] [runStart] [runEND] [eventType] [dataType] [runsPerJob] [YnsNo] [stNo] 
let stNo=1
# streamNo:[0,5]
while [ $stNo -le 1 ]
do
  ./autorunGenericMC_bookMCTree.csh 7  1 2865 evtgen-charged on_resonance 800 4 1$stNo
  ./autorunGenericMC_bookMCTree.csh 7  1 2865 evtgen-mixed   on_resonance 800 4 1$stNo
  ./autorunGenericMC_bookMCTree.csh 7  1 2865 evtgen-charm   on_resonance 800 4 1$stNo
  ./autorunGenericMC_bookMCTree.csh 7  1 2865 evtgen-uds     on_resonance 800 4 1$stNo

  ./autorunGenericMC_bookMCTree.csh 9  1 1220 evtgen-charged on_resonance 650 4 1$stNo
  ./autorunGenericMC_bookMCTree.csh 9  1 1220 evtgen-mixed   on_resonance 650 4 1$stNo
  ./autorunGenericMC_bookMCTree.csh 9  1 1220 evtgen-charm   on_resonance 650 4 1$stNo
  ./autorunGenericMC_bookMCTree.csh 9  1 1220 evtgen-uds     on_resonance 650 4 1$stNo

  ./autorunGenericMC_bookMCTree.csh 11 1 1287 evtgen-charged on_resonance 500 4 1$stNo
  ./autorunGenericMC_bookMCTree.csh 11 1 1287 evtgen-mixed   on_resonance 500 4 1$stNo
  ./autorunGenericMC_bookMCTree.csh 11 1 1287 evtgen-charm   on_resonance 500 4 1$stNo
  ./autorunGenericMC_bookMCTree.csh 11 1 1287 evtgen-uds     on_resonance 500 4 1$stNo

  ./autorunGenericMC_bookMCTree.csh 13 1 1627 evtgen-charged on_resonance 500 4 1$stNo
  ./autorunGenericMC_bookMCTree.csh 13 1 1627 evtgen-mixed   on_resonance 500 4 1$stNo
  ./autorunGenericMC_bookMCTree.csh 13 1 1627 evtgen-charm   on_resonance 500 4 1$stNo
  ./autorunGenericMC_bookMCTree.csh 13 1 1627 evtgen-uds     on_resonance 500 4 1$stNo

  ./autorunGenericMC_bookMCTree.csh 15 1 1437 evtgen-charged on_resonance 500 4 1$stNo
  ./autorunGenericMC_bookMCTree.csh 15 1 1437 evtgen-mixed   on_resonance 500 4 1$stNo
  ./autorunGenericMC_bookMCTree.csh 15 1 1437 evtgen-charm   on_resonance 500 4 1$stNo
  ./autorunGenericMC_bookMCTree.csh 15 1 1437 evtgen-uds     on_resonance 500 4 1$stNo

  ./autorunGenericMC_bookMCTree.csh 17 1 937  evtgen-charged on_resonance 500 4 1$stNo
  ./autorunGenericMC_bookMCTree.csh 17 1 937  evtgen-mixed   on_resonance 500 4 1$stNo
  ./autorunGenericMC_bookMCTree.csh 17 1 937  evtgen-charm   on_resonance 500 4 1$stNo
  ./autorunGenericMC_bookMCTree.csh 17 1 937  evtgen-uds     on_resonance 500 4 1$stNo

  ./autorunGenericMC_bookMCTree.csh 19 1 1643 evtgen-charged on_resonance 500 4 1$stNo
  ./autorunGenericMC_bookMCTree.csh 19 1 1643 evtgen-mixed   on_resonance 500 4 1$stNo
  ./autorunGenericMC_bookMCTree.csh 19 1 1643 evtgen-charm   on_resonance 500 4 1$stNo
  ./autorunGenericMC_bookMCTree.csh 19 1 1643 evtgen-uds     on_resonance 500 4 1$stNo

  ./autorunGenericMC_bookMCTree.csh 21 1 324  evtgen-charged on_resonance 400 4 1$stNo
  ./autorunGenericMC_bookMCTree.csh 21 1 324  evtgen-mixed   on_resonance 400 4 1$stNo
  ./autorunGenericMC_bookMCTree.csh 21 1 324  evtgen-charm   on_resonance 400 4 1$stNo
  ./autorunGenericMC_bookMCTree.csh 21 1 324  evtgen-uds     on_resonance 400 4 1$stNo

  ./autorunGenericMC_bookMCTree.csh 23 1 607  evtgen-charged on_resonance 700 4 1$stNo
  ./autorunGenericMC_bookMCTree.csh 23 1 607  evtgen-mixed   on_resonance 700 4 1$stNo
  ./autorunGenericMC_bookMCTree.csh 23 1 607  evtgen-charm   on_resonance 700 4 1$stNo
  ./autorunGenericMC_bookMCTree.csh 23 1 607  evtgen-uds     on_resonance 700 4 1$stNo

  ./autorunGenericMC_bookMCTree.csh 25 1 2122 evtgen-charged on_resonance 500 4 1$stNo
  ./autorunGenericMC_bookMCTree.csh 25 1 2122 evtgen-mixed   on_resonance 500 4 1$stNo
  ./autorunGenericMC_bookMCTree.csh 25 1 2122 evtgen-charm   on_resonance 500 4 1$stNo
  ./autorunGenericMC_bookMCTree.csh 25 1 2122 evtgen-uds     on_resonance 500 4 1$stNo

  ./autorunGenericMC_bookMCTree.csh 27 1 1632 evtgen-charged on_resonance 500 4 1$stNo
  ./autorunGenericMC_bookMCTree.csh 27 1 1632 evtgen-mixed   on_resonance 500 4 1$stNo
  ./autorunGenericMC_bookMCTree.csh 27 1 1632 evtgen-charm   on_resonance 500 4 1$stNo
  ./autorunGenericMC_bookMCTree.csh 27 1 1632 evtgen-uds     on_resonance 500 4 1$stNo

  ./autorunGenericMC_bookMCTree.csh 31 1 1715 evtgen-charged on_resonance 600 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 31 1 1715 evtgen-mixed   on_resonance 600 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 31 1 1715 evtgen-charm   on_resonance 600 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 31 1 1715 evtgen-uds     on_resonance 600 4 $stNo

  ./autorunGenericMC_bookMCTree.csh 33 1 870  evtgen-charged on_resonance 450 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 33 1 870  evtgen-mixed   on_resonance 450 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 33 1 870  evtgen-charm   on_resonance 450 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 33 1 870  evtgen-uds     on_resonance 450 4 $stNo

  ./autorunGenericMC_bookMCTree.csh 35 1 687  evtgen-charged on_resonance 700 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 35 1 687  evtgen-mixed   on_resonance 700 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 35 1 687  evtgen-charm   on_resonance 700 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 35 1 687  evtgen-uds     on_resonance 700 4 $stNo

  ./autorunGenericMC_bookMCTree.csh 37 1 1913 evtgen-charged on_resonance 500 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 37 1 1913 evtgen-mixed   on_resonance 500 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 37 1 1913 evtgen-charm   on_resonance 500 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 37 1 1913 evtgen-uds     on_resonance 500 4 $stNo

  ./autorunGenericMC_bookMCTree.csh 39 1 1357 evtgen-charged on_resonance 500 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 39 1 1357 evtgen-mixed   on_resonance 500 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 39 1 1357 evtgen-charm   on_resonance 500 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 39 1 1357 evtgen-uds     on_resonance 500 4 $stNo

  ./autorunGenericMC_bookMCTree.csh 41 1 1261 evtgen-charged on_resonance 650 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 41 1 1261 evtgen-mixed   on_resonance 650 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 41 1 1261 evtgen-charm   on_resonance 650 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 41 1 1261 evtgen-uds     on_resonance 650 4 $stNo

  ./autorunGenericMC_bookMCTree.csh 43 1 1149 evtgen-charged on_resonance 600 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 43 1 1149 evtgen-mixed   on_resonance 600 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 43 1 1149 evtgen-charm   on_resonance 600 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 43 1 1149 evtgen-uds     on_resonance 600 4 $stNo

  ./autorunGenericMC_bookMCTree.csh 45 1 450  evtgen-charged on_resonance 500 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 45 1 450  evtgen-mixed   on_resonance 500 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 45 1 450  evtgen-charm   on_resonance 500 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 45 1 450  evtgen-uds     on_resonance 500 4 $stNo

  ./autorunGenericMC_bookMCTree.csh 47 1 881  evtgen-charged on_resonance 450 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 47 1 881  evtgen-mixed   on_resonance 450 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 47 1 881  evtgen-charm   on_resonance 450 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 47 1 881  evtgen-uds     on_resonance 450 4 $stNo

  ./autorunGenericMC_bookMCTree.csh 49 1 922  evtgen-charged on_resonance 500 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 49 1 922  evtgen-mixed   on_resonance 500 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 49 1 922  evtgen-charm   on_resonance 500 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 49 1 922  evtgen-uds     on_resonance 500 4 $stNo

  ./autorunGenericMC_bookMCTree.csh 51 1 1777 evtgen-charged on_resonance 600 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 51 1 1777 evtgen-mixed   on_resonance 600 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 51 1 1777 evtgen-charm   on_resonance 600 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 51 1 1777 evtgen-uds     on_resonance 600 4 $stNo

  ./autorunGenericMC_bookMCTree.csh 55 1 1749 evtgen-charged on_resonance 500 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 55 1 1749 evtgen-mixed   on_resonance 500 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 55 1 1749 evtgen-charm   on_resonance 500 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 55 1 1749 evtgen-uds     on_resonance 500 4 $stNo

  ./autorunGenericMC_bookMCTree.csh 61 1 1207 evtgen-charged on_resonance 650 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 61 1 1207 evtgen-mixed   on_resonance 650 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 61 1 1207 evtgen-charm   on_resonance 650 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 61 1 1207 evtgen-uds     on_resonance 650 4 $stNo

  ./autorunGenericMC_bookMCTree.csh 63 1 783  evtgen-charged on_resonance 400 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 63 1 783  evtgen-mixed   on_resonance 400 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 63 1 783  evtgen-charm   on_resonance 400 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 63 1 783  evtgen-uds     on_resonance 400 4 $stNo

  ./autorunGenericMC_bookMCTree.csh 65 1 801  evtgen-charged on_resonance 401 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 65 1 801  evtgen-mixed   on_resonance 401 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 65 1 801  evtgen-charm   on_resonance 401 4 $stNo
  ./autorunGenericMC_bookMCTree.csh 65 1 801  evtgen-uds     on_resonance 401 4 $stNo

  let stNo+=1
done
