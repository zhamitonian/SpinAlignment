mkdir fit_result/temp_rootFile -p
mkdir fit_result/splot_rootFile -p
find ./fit_result/ -maxdepth 1 -name "*.root" -exec mv {} fit_result/temp_rootFile/ \;
find ./fit_result/temp_rootFile/ -name "*splot*" -exec mv {} fit_result/splot_rootFile/ \;