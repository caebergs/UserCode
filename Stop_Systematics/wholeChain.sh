#!/bin/bash

outputFromVJetsEstimation="/user/caebergs/VJetEstimation/nominal_STbar017/stdout_err_w.txt" ;
#nomDir="/user/caebergs/Leg3_WnJets_DataMC/miguelf/" ;
nomDir="/user/caebergs/Leg3/Leg3Output_5Dreweigh/" ;
systDir="/user/caebergs/VJetEstimation/testChain/samplesSyst/" ;
controlRegionsDir="/user/caebergs/Leg3/Leg3NewTEff/" ;
extrapolatedRV="false" ;

nbOfUseCase=4 ;
timestamp=$(date +%Y%m%d_%H%M%S) ;

mkdir ${timestamp}_workingDir ;
cd ${timestamp}_workingDir ;
echo "cd ${timestamp}_workingDir"
echo -e "\n --> Extracting content from the V+jets estimation ....\n"

# 1 : muon
# 2 : electron
# 3 : both
for channel in 1 2 3 ; do
    sh /user/caebergs/VJetEstimation/extractCorrMatrix.sh ${outputFromVJetsEstimation} ${channel} > corrMatr_${channel}.txt ;
    sh /user/caebergs/VJetEstimation/extractFitResult.sh ${outputFromVJetsEstimation} ${channel} > fitResult_${channel}.txt ;
done ;

echo -e "\n --> Neural Network cut efficiency and systematics ....\n" ;
#Uncertainties on RTT
cat ../makeSystUncertSummaryPlot.C \
    | sed -e "s#^\([ ]*std::string path = \"\)\(.*\)\(\"[ ]*;\)#\1${systDir}\3#" \
    | sed -e "s#^\([ ]*std::string nompath = \"\)\(.*\)\(\"[ ]*;\)#\1${nomDir}\3#" \
    > makeSystUncertSummaryPlot.C ;
#Uncertainties on RV 
cat ../makeSystUncertSummaryPlot_VJets.C \
    | sed -e "s#^\([ ]*std::string path = \"\)\(.*\)\(\"[ ]*;\)#\1${systDir}\3#" \
    | sed -e "s#^\([ ]*std::string nompath = \"\)\(.*\)\(\"[ ]*;\)#\1${nomDir}\3#" \
    > makeSystUncertSummaryPlot_VJets.C ;
for channel in 1 2 3 ; do
# Actual computations of the uncertainties on RTT and RV
    for useCase in 0 1 2 3 ; do
	cat <<EOD | root -l -b > syst_RTT_channel${channel}_${useCase}.txt 2>&1
.L makeSystUncertSummaryPlot.C++
makeSystUncertSummaryPlot(${useCase}, ${channel})
.q
EOD
	cat <<EOD | root -l -b > syst_RV_channel${channel}_${useCase}.txt 2>&1
.L makeSystUncertSummaryPlot_VJets.C++
makeSystUncertSummaryPlot_VJets(${useCase}, ${channel})
.q
EOD
    done ;
done ;

echo " ---> Trend systematics ...." ;
cp ../makeSystUncertTrendPlot_Summary_VJets.C makeSystUncertTrendPlot_Summary_VJets.C ;
for channel in 1 2 3 ; do
    cat makeSystUncertTrendPlot_Summary_VJets.C \
	| sed -e "s#^\([ ]*std::string path = \"\)\(.*\)\(\"[ ]*;\)#\1${systDir}\3#" \
	| sed -e "s#^\([ ]*std::string nompath = \"\)\(.*\)\(\"[ ]*;\)#\1${nomDir}\3#" \
	| sed -e "s/channel[ ]*=[ ]*[0-9]*[ ]*;/channel = ${channel} ;/" \
	> makeSystUncertTrendPlot_Summary_VJets_${channel}.C ;
done ;
#	| sed -e "s/_channelXXX\"/_channel${channel}\"/" \
cat ../makeSystUncertTrendPlot_VJets.C \
    | sed -e "s#^\([ ]*std::string path = \"\)\(.*\)\(\"[ ]*;\)#\1${systDir}\3#" \
    | sed -e "s#^\([ ]*std::string nompath = \"\)\(.*\)\(\"[ ]*;\)#\1${nomDir}\3#" \
    > makeSystUncertTrendPlot_VJets.C ;
cat ../makeSystUncertTrendPlot_JES_VJets.C \
    | sed -e "s#^\([ ]*std::string path = \"\)\(.*\)\(\"[ ]*;\)#\1${systDir}\3#" \
    | sed -e "s#^\([ ]*std::string nompath = \"\)\(.*\)\(\"[ ]*;\)#\1${nomDir}\3#" \
    > makeSystUncertTrendPlot_JES_VJets.C ;

for useCase in 0 1 2 3 ; do
    declare -i syst ;
    for syst in -1 0 1 2 ; do
	if [ "${syst}" = "-1" ] ; then
	    for channel in 1 2 3 ; do
		declare -i index ;
		index=0 ;
#		echo "Production channel : ${channel}  tmp(useCase ${useCase}, syst ${syst})" ;
		cat <<EOD | root -l -b > syst_RV_Trend_syst${syst}_channel${channel}_${useCase}.txt 2>&1
.L makeSystUncertTrendPlot_JES_VJets.C++
makeSystUncertTrendPlot_JES_VJets(${useCase}, ${channel})
.q
EOD
	    done ;
	else
	    for channel in 1 2 3 ; do
#		echo "Production channel : ${channel}  tmp(useCase ${useCase}, syst ${syst})"
		cat <<EOD | root -l -b > syst_RV_Trend_syst${syst}_channel${channel}_${useCase}.txt 2>&1 ;
.L makeSystUncertTrendPlot_VJets.C++
makeSystUncertTrendPlot_VJets(${useCase}, ${syst}, ${channel})
.q
EOD
	    done ;
	fi ;
	for channel in 1 2 3 ; do
	    echo "Resultats channel : ${channel}  tmp(useCase ${useCase}, syst ${syst})" ;
#echo "channel : $channel"
	    declare -i index ;
	    index=0 ;
#cat syst_RV_Trend_syst${syst}_channel${channel}_${useCase}.txt
#		| grep -c "Extrapolated uncert. "
	    if [ -f syst_RV_Trend_syst${syst}_channel${channel}_${useCase}.txt ] ; then 
		for relErr in $(cat syst_RV_Trend_syst${syst}_channel${channel}_${useCase}.txt \
		    | grep "Extrapolated uncert" \
		    | sed -e "s/^.*Extrapolated uncert\. (.*) at [0-9\.]* = \([0-9eE\.+\-]*\)[ ]*$/\1/" ) ; do 
		  echo "${index} : ${relErr}" ;
		  index=${index}+1 ;
		  indexPattern=${index};
		  find5D=$(cat syst_RV_Trend_syst${syst}_channel${channel}_${useCase}.txt \
		      | grep "Extrapolated uncert\. (5D Rew\. (" ) ;
		  if [ "${find5D}" != "" ] ; then 
		      indexPattern="[01]" ;
		  fi ;
		  cat makeSystUncertTrendPlot_Summary_VJets_${channel}.C \
		      | sed -e "s/\(rel_syst_uncert\[${useCase}\]\[$((${syst}+1))\]\[${indexPattern}\]\)=.*; /\1 = ${relErr} ; \//" \
		      > tmp.C ;
		  mv tmp.C makeSystUncertTrendPlot_Summary_VJets_${channel}.C ;
		done ;
	    fi ;
	done ;
    done ;
done ;


echo "Rhaha"
declare -a RTTrelP ;
declare -a RTTrelM ;
declare -a RTTnom ;
declare -a RVrelP ;
declare -a RVrelM ;
declare -a RVnom ;
declare -a RVnomTrend ;
for channel in 1 2 3 ; do
    for useCase in 0 1 2 3 ; do
	RTTrelP[${channel}*${nbOfUseCase}+${useCase}]=$(cat syst_RTT_channel${channel}_${useCase}.txt | grep "Rel\. Comb\. (+) (stat\.+syst\.) uncert : [0-9]*\.[0-9]*" | sed -e "s/Rel\. Comb\. (+) (stat\.+syst\.) uncert : \([0-9\.]*\)/\1/")
	RTTrelM[${channel}*${nbOfUseCase}+${useCase}]=$(cat syst_RTT_channel${channel}_${useCase}.txt | grep "Rel\. Comb\. (-) (stat\.+syst\.) uncert : [0-9]*\.[0-9]*" | sed -e "s/Rel\. Comb\. (-) (stat\.+syst\.) uncert : \([0-9\.]*\)/\1/")
	RTTnom[${channel}*${nbOfUseCase}+${useCase}]=$(cat syst_RTT_channel${channel}_${useCase}.txt | grep "\\$R_{TT}\\$ & \$[0-9]*\.[0-9]*\\\pm[0-9]*\.[0-9]*^{+[0-9]*\.[0-9]*}_{-[0-9]*\.[0-9]*}\\$" | head -n 1 | sed -e "s/\\\$R_{TT}\\\$ \& \\\$\([0-9\.]*\)\\\\pm[0-9\.]*^{+[0-9\.]*}_{-[0-9\.]*}\\\$/\1/")
	RVrelP[${channel}*${nbOfUseCase}+${useCase}]=$(cat syst_RV_channel${channel}_${useCase}.txt | grep "Rel\. Comb\. (+) (stat\.+syst\.) uncert : [0-9]*\.[0-9]*" | sed -e "s/Rel\. Comb\. (+) (stat\.+syst\.) uncert : \([0-9\.]*\)/\1/")
	RVrelM[${channel}*${nbOfUseCase}+${useCase}]=$(cat syst_RV_channel${channel}_${useCase}.txt | grep "Rel\. Comb\. (-) (stat\.+syst\.) uncert : [0-9]*\.[0-9]*" | sed -e "s/Rel\. Comb\. (-) (stat\.+syst\.) uncert : \([0-9\.]*\)/\1/")
	RVnom[${channel}*${nbOfUseCase}+${useCase}]=$(cat syst_RV_channel${channel}_${useCase}.txt | grep "\\$R_{V}\\$ & \$[0-9]*\.[0-9]*\\\pm[0-9]*\.[0-9]*^{+[0-9]*\.[0-9]*}_{-[0-9]*\.[0-9]*}\\$" | head -n 1 | sed -e "s/\\\$R_{V}\\\$ \& \\\$\([0-9\.]*\)\\\\pm[0-9\.]*^{+[0-9\.]*}_{-[0-9\.]*}\\\$/\1/")
	echo "Channel ${channel} , UseCase ${useCase}" ;
	echo "  RTT : ${RTTnom[${channel}*${nbOfUseCase}+${useCase}]} +rel ${RTTrelP[${channel}*${nbOfUseCase}+${useCase}]} -rel ${RTTrelM[${channel}*${nbOfUseCase}+${useCase}]}" ;
	echo "  RV : ${RVnom[${channel}*${nbOfUseCase}+${useCase}]} +rel ${RVrelP[${channel}*${nbOfUseCase}+${useCase}]} -rel ${RVrelM[${channel}*${nbOfUseCase}+${useCase}]}" ;
    done ;
done ;
echo "Ici"
for channel in 1 2 3 ; do
    var1=${RVnom[$((${channel}*${nbOfUseCase}+0))]} ;
    var2=${RVnom[$((${channel}*${nbOfUseCase}+1))]} ;
    var3=${RVnom[$((${channel}*${nbOfUseCase}+2))]} ;
    var4=${RVnom[$((${channel}*${nbOfUseCase}+3))]} ;
    cat makeSystUncertTrendPlot_Summary_VJets_${channel}.C \
	| sed -e "s/\(  const double RV\[NbOfUseCase\] = {\)[0-9eE\.+\-]*,[0-9eE\.+\-]*,[0-9eE\.+\-]*,[0-9eE\.+\-]*\(};\)/\1${var1},${var2},${var3},${var4}\2/" \
	> tmp.C
    mv tmp.C makeSystUncertTrendPlot_Summary_VJets_${channel}.C ;
    grep -H "const double RV" makeSystUncertTrendPlot_Summary_VJets_${channel}.C ;
done ;
echo "La"
for channel in 1 2 3 ; do
    for useCase in 0 1 2 3 ; do
	cat <<EOD | root -l -b > syst_RV_Trend_Summary_channel${channel}_${useCase}.txt 2>&1
.L makeSystUncertTrendPlot_Summary_VJets_${channel}.C++g
makeSystUncertTrendPlot_Summary_VJets(${useCase})
.q
EOD
	relP=$(cat syst_RV_Trend_Summary_channel${channel}_${useCase}.txt | grep "Rel\. Comb\. (+) syst\. uncert : [0-9eE\.+\-]*" | sed -e "s/Rel\. Comb\. (+) syst\. uncert : \([0-9eE\.+\-]*\)/\1/")
	relM=$(cat syst_RV_Trend_Summary_channel${channel}_${useCase}.txt | grep "Rel\. Comb\. (-) syst\. uncert : [0-9eE\.+\-]*" | sed -e "s/Rel\. Comb\. (-) syst\. uncert : \([0-9eE\.+\-]*\)/\1/")
	nom=$(cat syst_RV_Trend_Summary_channel${channel}_${useCase}.txt | grep "R_{V} = [0-9eE\.+\-]*^{+[0-9eE\.+\-]*}_{-[0-9eE\.+\-]*}" | head -n 1 | sed -e "s/R_{V} = \([0-9eE\.+\-]*\)^{+[0-9eE\.+\-]*}_{-[0-9eE\.+\-]*}$/\1/")
	echo "  RV combined (ch ${channel} , case ${useCase}) : ${nom} +rel ${relP} -rel ${relM}" ;
	if [ "${extrapolatedRV}" = "true" ] ; then
	    echo "Replacing    ${RVrelP[${channel}*${nbOfUseCase}+${useCase}]}   by   ${relP}" ;
	    RVrelP[${channel}*${nbOfUseCase}+${useCase}]=${relP} ;
	    echo "Replacing    ${RVrelM[${channel}*${nbOfUseCase}+${useCase}]}   by   ${relM}" ;
	    RVrelM[${channel}*${nbOfUseCase}+${useCase}]=${relM} ;
	    echo "Replacing    ${RVnom[${channel}*${nbOfUseCase}+${useCase}]}   by   ${nom}" ;
	    RVnom[${channel}*${nbOfUseCase}+${useCase}]=${nom} ;
	fi ;
    done
done
echo "Ici aussi"
#exit 1 ;


#Extrapolation of uncertainties
#cat ../makeSystUncertTrendPlot_Summary_VJets.C \
#| sed -e "s///" \
#> makeSystUncertTrendPlot_Summary_VJets.C
#
#root -l makeSystUncertTrendPlot_Summary_VJets.C++
#

echo "Individuel"
# 0 : muon
# 1 : electron
# 2 : 
cat ../TotalEstimatedNumbers_Errors_IndivChannel.C \
    | sed -e "s#\(^[ ]*std::string path = \"\)\(.*\)\(\"[ ]*;\)#\1${nomDir}\3#" \
    > tmp.txt
mv tmp.txt TotalEstimatedNumbers_Errors_IndivChannel.C ;
for useCase in 0 1 2 3 ; do
    cat TotalEstimatedNumbers_Errors_IndivChannel.C \
	| sed -e  "s/\( RTT_syst_err_up[${useCase}]  = (channel==0 ? \)[0-9\.]* : [0-9\.]*);/ \1${RTT_relP[1*${nbOfUseCase}+${useCase}]} : ${RTT_relP[2*${nbOfUseCase}+${useCase}]});/" \
	| sed -e  "s/\( RTT_syst_err_low[${useCase}] = (channel==0 ? \)[0-9\.]* : [0-9\.]*);/ \1${RTT_relP[1*${nbOfUseCase}+${useCase}]} : ${RTT_relP[2*${nbOfUseCase}+${useCase}]});/" \
	| sed -e "s/\( double RV_syst_err_up\[4\]  = \){[0-9\.]*\*y_vjets,[0-9\.]*\*y_vjets,[0-9\.]*\*y_vjets,[0-9\.]*\*y_vjets}; /\1{${RVrelP[3*${nbOfUseCase}+0]} *y_vjets,${RVrelP[3*${nbOfUseCase}+1]} *y_vjets,${RVrelP[3*${nbOfUseCase}+2]} *y_vjets,${RVrelP[3*${nbOfUseCase}+3]} *y_vjets}; /" \
	| sed -e "s/\( double RV_syst_err_low\[4\] = \){[0-9\.]*\*y_vjets,[0-9\.]*\*y_vjets,[0-9\.]*\*y_vjets,[0-9\.]*\*y_vjets}; /\1{${RVrelM[3*${nbOfUseCase}+0]}*y_vjets,${RVrelM[3*${nbOfUseCase}+1]}*y_vjets,${RVrelM[3*${nbOfUseCase}+2]}*y_vjets,${RVrelM[3*${nbOfUseCase}+3]}*y_vjets}; /" \
	> tmp.txt ;
    mv tmp.txt TotalEstimatedNumbers_Errors_IndivChannel.C ;
done
echo "Groupe"
cat ../TotalEstimatedNumbers_Errors_new.C \
    | sed -e "s#\(^[ ]*std::string path = \"\)\(.*\)\(\"[ ]*;\)#\1${nomDir}\3#" \
    | sed -e "s/\( double RTT_syst_err_up\[4\]  = \){[0-9\.]*,[0-9\.]*,[0-9\.]*,[0-9\.]*}; /\1{${RTTrelP[3*${nbOfUseCase}+0]}*y_ttlike,${RTTrelP[3*${nbOfUseCase}+1]}*y_ttlike,${RTTrelP[3*${nbOfUseCase}+2]}*y_ttlike,${RTTrelP[3*${nbOfUseCase}+3]}*y_ttlike}; /" \
    | sed -e "s/\( double RTT_syst_err_low\[4\] = \){[0-9\.]*,[0-9\.]*,[0-9\.]*,[0-9\.]*}; /\1{${RTTrelM[3*${nbOfUseCase}+0]}*y_ttlike,${RTTrelM[3*${nbOfUseCase}+1]}*y_ttlike,${RTTrelM[3*${nbOfUseCase}+2]}*y_ttlike,${RTTrelM[3*${nbOfUseCase}+3]}*y_ttlike}; /" \
    | sed -e "s/\( double RV_syst_err_up\[4\]  = \){[0-9\.]*\*y_vjets,[0-9\.]*\*y_vjets,[0-9\.]*\*y_vjets,[0-9\.]*\*y_vjets}; /\1{${RVrelP[3*${nbOfUseCase}+0]} *y_vjets,${RVrelP[3*${nbOfUseCase}+1]} *y_vjets,${RVrelP[3*${nbOfUseCase}+2]} *y_vjets,${RVrelP[3*${nbOfUseCase}+3]} *y_vjets}; /" \
    | sed -e "s/\( double RV_syst_err_low\[4\] = \){[0-9\.]*\*y_vjets,[0-9\.]*\*y_vjets,[0-9\.]*\*y_vjets,[0-9\.]*\*y_vjets}; /\1{${RVrelM[3*${nbOfUseCase}+0]}*y_vjets,${RVrelM[3*${nbOfUseCase}+1]}*y_vjets,${RVrelM[3*${nbOfUseCase}+2]}*y_vjets,${RVrelM[3*${nbOfUseCase}+3]}*y_vjets}; /" \
    > tmp.txt ;
mv tmp.txt TotalEstimatedNumbers_Errors_new.C

cat ../ControlRegions.C \
| sed -e "s#\(^[ ]*std::string path = \"\)\(.*\)\(\"[ ]*;\)#\1${controlRegionsDir}\3#" \
> tmp.txt
mv tmp.txt ControlRegions.C ;


for channel in 1 2 3 ; do
#    cp TotalEstimatedNumbers_Errors_new.C TotalEstimatedNumbers_Errors_new_${channel}.C ;
    declare -i idx=0 ;
    idx=0;
    echo -e -n "channel : $channel\n" ;
    ( cat fitResult_${channel}.txt \
	| while read ligne ; do #in $(cat fitResult_${channel}.txt) ; do
#echo "a"
	    varName=$( echo "${ligne}" | sed -e    "s/^[ ]*\([^ ]*\)[ ]*\([^ ]*\)[ ]*\([^ ]*\)[ ]*(\([^ ]*\),\([^ ]*\))[ ]*<none>$/\1/" ) ;
	    echo -n "varName=$varName   "
	    estimation=$( echo "${ligne}" | sed -e "s/^[ ]*\([^ ]*\)[ ]*\([^ ]*\)[ ]*\([^ ]*\)[ ]*(\([^ ]*\),\([^ ]*\))[ ]*<none>$/\3/" ) ;
	    echo -n "estimation=$estimation   "
	    errPlus=$( echo "${ligne}" | sed -e    "s/^[ ]*\([^ ]*\)[ ]*\([^ ]*\)[ ]*\([^ ]*\)[ ]*(\([^ ]*\),\([^ ]*\))[ ]*<none>$/\4/" ) ;
	    echo -n "errPlus=$errPlus   "
	    errMinus=$( echo "${ligne}" | sed -e   "s/^[ ]*\([^ ]*\)[ ]*\([^ ]*\)[ ]*\([^ ]*\)[ ]*(\([^ ]*\),\([^ ]*\))[ ]*<none>$/\5/" ) ;
	    echo -n "errMinux=$errMinus   "
	    errP=$( echo "${errPlus}" | sed -e "s/e/\*10\^/" | sed -e "s/\+//g" ) ;
	    echo -n "errP=$errP   "
	    errM=$( echo "${errMinus}" | sed -e "s/e/\*10\^/" | sed -e "s/\+//g" ) ;
	    echo -n "errM=$errM   "
	    err=$( echo "scale=90; if ( ${errP} + ( -1 * ${errM} ) > 0 ) ${errP} else ( -1 * ${errM} ) " | tee test.txt| bc | sed -e 's/\\//' | tr -d '\n' | sed -e "s/[0]*$//" ) ;
#cat test.txt
	    echo -e -n "err=$err\n"
#systUncert : eXbq sur VJets method ????
cat ControlRegions.C \
| sed -e "s/^  nbOfEvents\[[ ]*${idx}[ ]*\] =[ ]*\([0-9\.]*\);/  nbOfEvents[${idx}] = ${estimation} ; \/\/ ${varName} \/\/ /" \
| sed -e "s/^  statUncert\[[ ]*${idx}[ ]*\] =[ ]*\([0-9\.]*\);/  statUncert[${idx}] = ${err} ; \/\/ ${varName} \/\/ /" \
> tmp.txt
mv tmp.txt ControlRegions_${channel}.C

	    if [ "${channel}" = "3" ]; then 
		cat TotalEstimatedNumbers_Errors_new.C \
		    | sed -e "s/^  nbOfEvents\[[ ]*${idx}[ ]*\] =[ ]*\([0-9]*\);/  nbOfEvents[${idx}] = ${estimation} ; \/\/ ${varName} \/\/ /" \
		    | sed -e "s/^  statUncert\[[ ]*${idx}[ ]*\] =[ ]*\([0-9]*\);/  statUncert[${idx}] = ${err} ; \/\/ ${varName} \/\/ /" \
		    > tmp.txt
		mv tmp.txt TotalEstimatedNumbers_Errors_new.C #TotalEstimatedNumbers_Errors_new_${channel}.C ;
	    elif [ "${channel}" = "1" ]; then
		cat  TotalEstimatedNumbers_Errors_IndivChannel.C \
		    | sed -e "s/^  nbOfEvents\[[ ]*${idx}[ ]*\] = (channel==0 ?[ ]*\([0-9]*\)[ ]*:[ ]*\([0-9]*\)[ ]*);/  nbOfEvents[${idx}] = (channel==0 ? ${estimation} : \2 ) ; \/\/ ${varName} \/\/ /" \
		    | sed -e "s/^  statUncert\[[ ]*${idx}[ ]*\] = (channel==0 ?[ ]*\([0-9]*\)[ ]*:[ ]*\([0-9]*\)[ ]*);/  statUncert[${idx}] = (channel==0 ? ${err} : \2 ); \/\/ ${varName} \/\/ /" \
		    > tmp.txt
		mv tmp.txt TotalEstimatedNumbers_Errors_IndivChannel.C
	    elif [ "${channel}" = "1" ]; then
		cat  TotalEstimatedNumbers_Errors_IndivChannel.C \
		    | sed -e "s/^  nbOfEvents\[[ ]*${idx}[ ]*\] = (channel==0 ?[ ]*\([0-9]*\)[ ]*:[ ]*\([0-9]*\)[ ]*);/  nbOfEvents[${idx}] = (channel==0 ? \1 : ${estimation} ) ; \/\/ ${varName} \/\/ /" \
		    | sed -e "s/^  statUncert\[[ ]*${idx}[ ]*\] = (channel==0 ?[ ]*\([0-9]*\)[ ]*:[ ]*\([0-9]*\)[ ]*);/  statUncert[${idx}] = (channel==0 ? \1 : ${err} ); \/\/ ${varName} \/\/ /" \
		    > tmp.txt
		mv tmp.txt TotalEstimatedNumbers_Errors_IndivChannel.C
	    fi
	    idx=${idx}+1 ;
	done
    ) ;
    echo "" ;
done ;

#  double RTT_syst_err_up[4]  = {0.000742,0.000565,0.000727,0.001017}; // Up
#  double RTT_syst_err_low[4] = {0.000892,0.000581,0.001316,0.001055}; // Down
#  double RV_syst_err_up[4]  = {0.0482035*y_vjets,0.0226508*y_vjets,0.0233062*y_vjets,0.0149049*y_vjets}; // Up
#  double RV_syst_err_low[4] = {0.2622430*y_vjets,0.3797950*y_vjets,0.3956500*y_vjets,0.4722110*y_vjets}; // Down


#| sed -e "s#  systUncert\[[ ]*${idx}[ ]*\] =[ ]*\([0-9]*\);##"

#exit 1;
echo -e "\n --> Computing the total numbers ....\n" ;

for useCase in 0 1 2 3 ; do
    for channel in 1 2 3 ; do
        cat <<EOD | root -l -b > ControlRegions_channel_${channel}_${useCase}_${useCase}.txt 2>&1
.L ControlRegions_${channel}.C++
Int_t bin = (${useCase}==0?18:19)*0
bin
Bool_t UseWNJets = kFALSE
UseWNJets
ControlRegions("corrMatr_${channel}.txt", ${useCase}, bin, UseWNJets, ${channel})
.q
EOD
	if [ "${channel}" == "3" ] ; then
	    cat <<EOD | root -l -b > TotNumbers_channel_${channel}_${useCase}.txt 2>&1
.L TotalEstimatedNumbers_Errors_new.C++
Int_t bin = (${useCase}==0?18:19)*0
bin
Bool_t UseWNJets = kFALSE
UseWNJets
TotalEstimatedNumbers_Errors_new("corrMatr_${channel}.txt", ${useCase}, bin, UseWNJets)
.q
EOD
	else
	    cat <<EOD | root -l -b > TotNumbers_channel_${channel}_${useCase}.txt 2>&1
.L TotalEstimatedNumbers_Errors_IndivChannel.C++
Int_t bin = (${useCase}==0?18:19)*0+1
bin
Bool_t UseWNJets = kFALSE
UseWNJets
TotalEstimatedNumbers_Errors_IndivChannel("corrMatr_${channel}.txt", ${useCase}, bin, UseWNJets, ${channel}-1)
.q
EOD
	fi;
	echo "UseCase : ${useCase} , channel : ${channel} : $(cat TotNumbers_channel_${channel}_${useCase}.txt | grep 'Ntotal = ')" ;
    done
done

