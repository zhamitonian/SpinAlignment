#! /bin/bash -f

#argument: $0 [expNo] [runStart] [runEND] [eventType] [dataType] [runsPerJob] [YnsNo] [stNo] 
let stNo=0
# streamNo:[0,5]
while [ $stNo -le 0 ]
do
  ./autorunGenericMC_bookMCTree.csh 7  1490 2293 evtgen-charm  continuum 5000 40 2$stNo
  ./autorunGenericMC_bookMCTree.csh 7  1490 1536 evtgen-uds    continuum 500 40 2$stNo
  ./autorunGenericMC_bookMCTree.csh 7  2163 2293 evtgen-uds    continuum 500 40 2$stNo
  ./autorunGenericMC_bookMCTree.csh 11 165  1367 evtgen-charm  continuum 1000 40 2$stNo
  ./autorunGenericMC_bookMCTree.csh 11 165  1367 evtgen-uds    continuum 1000 40 2$stNo
  ./autorunGenericMC_bookMCTree.csh 13 454  1587 evtgen-charm  continuum 1000 40 2$stNo
  ./autorunGenericMC_bookMCTree.csh 13 454  1587 evtgen-uds    continuum 1000 40 2$stNo
  ./autorunGenericMC_bookMCTree.csh 15 714  1314 evtgen-charm  continuum 700  40 2$stNo
  ./autorunGenericMC_bookMCTree.csh 15 714  1314 evtgen-uds    continuum 700  40 2$stNo
  ./autorunGenericMC_bookMCTree.csh 17 486  544  evtgen-charm  continuum 100  40 2$stNo
  ./autorunGenericMC_bookMCTree.csh 17 486  544  evtgen-uds    continuum 100  40 2$stNo
  ./autorunGenericMC_bookMCTree.csh 19 1    1709 evtgen-charm  continuum 1000 40 2$stNo
  ./autorunGenericMC_bookMCTree.csh 19 1    1709 evtgen-uds    continuum 1000 40 2$stNo
  ./autorunGenericMC_bookMCTree.csh 23 244  356  evtgen-charm  continuum 200  40 2$stNo
  ./autorunGenericMC_bookMCTree.csh 23 244  356  evtgen-uds    continuum 200  40 2$stNo
  ./autorunGenericMC_bookMCTree.csh 25 1468 1599 evtgen-charm  continuum 200  40 2$stNo
  ./autorunGenericMC_bookMCTree.csh 25 1468 1599 evtgen-uds    continuum 200  40 2$stNo
  ./autorunGenericMC_bookMCTree.csh 27 290  1251 evtgen-charm  continuum 500 40 2$stNo
  ./autorunGenericMC_bookMCTree.csh 27 290  1251 evtgen-uds    continuum 500 40 2$stNo

  ./autorunGenericMC_bookMCTree.csh 31 937  1064 evtgen-charm continuum 200  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 31 937  1064 evtgen-uds   continuum 200  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 33 505  624  evtgen-charm continuum 200  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 33 505  624  evtgen-uds   continuum 200  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 35 588  641  evtgen-charm continuum 100  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 35 588  641  evtgen-uds   continuum 100  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 37 659  1664 evtgen-charm continuum 500  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 37 659  1664 evtgen-uds   continuum 500  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 39 636  1203 evtgen-charm continuum 300  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 39 636  1203 evtgen-uds   continuum 300  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 41 703  1157 evtgen-charm continuum 250  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 41 703  1157 evtgen-uds   continuum 250  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 43 559  972  evtgen-charm continuum 300  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 43 559  972  evtgen-uds   continuum 300  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 45 383  421  evtgen-charm continuum 100  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 45 383  421  evtgen-uds   continuum 100  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 47 550  622  evtgen-charm continuum 100  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 47 550  622  evtgen-uds   continuum 100  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 49 553  706  evtgen-charm continuum 200  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 49 553  706  evtgen-uds   continuum 200  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 51 1312 1805 evtgen-charm continuum 250  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 51 1312 1805 evtgen-uds   continuum 250  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 55 793  1677 evtgen-charm continuum 500  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 55 793  1677 evtgen-uds   continuum 500  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 61 668  739  evtgen-charm continuum 100  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 61 668  739  evtgen-uds   continuum 100  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 63 618  679  evtgen-charm continuum 100  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 63 618  679  evtgen-uds   continuum 100  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 65 626  687  evtgen-charm continuum 100  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 65 626  687  evtgen-uds   continuum 100  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 67 698  742  evtgen-charm continuum 100  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 67 698  742  evtgen-uds   continuum 100  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 69 823  887 evtgen-charm continuum 100  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 69 1311  1397 evtgen-charm continuum 100  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 69 823  887 evtgen-uds   continuum 100  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 69 1311  1397 evtgen-uds   continuum 100  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 71 2249 2292 evtgen-charm continuum 100  40  $stNo
  ./autorunGenericMC_bookMCTree.csh 71 2249 2292 evtgen-uds   continuum 100  40  $stNo

  let stNo+=1
done
