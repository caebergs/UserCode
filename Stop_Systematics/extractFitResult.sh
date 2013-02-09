#!/bin/bash
fichier="stdout_err_w.txt"
if [ "$1" != "" ]; then fichier="$1" ; fi
channel="3"
if [ "$2" != "" ]; then channel="$2" ; fi
declare -i lineStart=0 ;
declare -i lineStop=0 ;
declare -i lineFile=0 ;
linesMatrixFollows=$(cat ${fichier} | grep -n "    Floating Parameter  InitialValue    FinalValue (+HiError,-LoError)    GblCorr." | cut -d":" -f1 | tr "\n" "\t") ;
#echo "${linesMatrixFollows}"
linesMatrixFollowsPruned=$(echo "${linesMatrixFollows}" | sed -e "s/^.*\t\([0-9]*\)\t\([0-9]*\)\t\([0-9]*\)\t$/\\${channel}/") ; #select the channel'th latest
#echo "${linesMatrixFollowsPruned}"
for lineFile in ${linesMatrixFollowsPruned} ; do
    lineStart=${lineFile}+2 ;
#    echo "lineFile = ${lineFile}" ;
    lineStop=${lineStart} ;
    tmpLine=$(cat ${fichier} | sed -n -e "${lineStop},+0 p") ;
    while [ "${tmpLine}" != "" ]; do
	lineStop=${lineStop}+1 ;
	tmpLine=$(cat ${fichier} | sed -n -e "${lineStop},+0 p") ;
    done
    declare -i length=${lineStop}-${lineStart}-1 ;
    cat ${fichier} | sed -n -e "${lineStart},+${length} p" ;
#    matrixSize=$(cat ${fichier} | sed -n -e "${lineFile},+0 p" |  sed -e "s/^\([0-9]*\)x\([0-9]*\) matrix is as follows$/\1/") ;
#    echo "matrixSize = ${matrixSize}" ;
#    numLines=$(((4+matrixSize)*(((matrixSize-1)/5)+1))) ;
#    echo "numLines = ${numLines}" ;
#    start=3 ;
#    for (( start=5 ; start < $((${matrixSize}+5))  ; start=${start}+1 )) do
#     ligne=$(cat ${fichier} | sed -n -e "${lineFile},+${numLines} p" | cut -d"|" -f2- | sed -n -e "${start}~$((${matrixSize}+4)) p" | tr -d "\n") ;
#echo "${ligne}" ;
#    done
done

#for a in 
#cat stdout_err_w.txt | grep " matrix is as follows" | sed -e "s/^\([0-9]*\)x\([0-9]*\) matrix is as follows$/\1/"
