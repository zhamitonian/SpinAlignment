#argument: $0 [expNo] [runStart] [runEND] [skimType] [dataType] [runsPerJob] [YnsNo] 
## Y(1S)
./autorunExpdata.csh 65 1008 1232 HadronBJ 1S_scan  60 1

## Y(2S)
./autorunExpdata.csh 67 1002 1123 HadronBJ 2S_scan  40 2
./autorunExpdata.csh 71 303  696  HadronBJ 2S_scan  70 2

## Y(3S)
#./autorunExpdata.csh 49 1001 1227 HadronBJ 3S_scan 120 3
