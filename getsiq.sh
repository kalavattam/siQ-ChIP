#!/bin/bash

echo "We assume your annotations are sorted as 'sort -k1,1 -k2,2n' otherwise this will fail"
#MAYBE: Write a function that checks that a BED is properly sorted
echo "Make sure your annotations are linked 'ln -s SORTEDannotations.bed ./Annotations.bed' HERE"
annot=./Annotations.bed  #TODO Make this into a keyword parameter
det_exp="./EXPlayout"  #TODO Make this into a keyword parameter
widths=30  # Bin size for signal processing  #TODO Make this into a keyword parameter
egs=12157105  # Effective genome size for model organism  #TODO Make this into a keyword parameter


if [[ ! -f "${annot}" ]]; then
    echo "Error: Missing ${annot}" >&2
    # exit 1
    touch "${annot}"
fi

echo "Proceeding..."
#Sanity check for the experiment details/layout file
if ! test -f "${det_exp}" ; then
    echo "Error: Missing experiment layout file" >&2
    exit 1
fi

z=$(awk '/getTracks/ { print NR }' "${det_exp}")
x=$(awk '/getResponse/ { print NR }' "${det_exp}")
c=$(awk '/getFracts/ { print NR }' "${det_exp}")
v=$(awk '/END/ { print NR }' "${det_exp}")
if [[ ${z} -lt ${x} ]] && [[ ${x} -lt ${c} ]] && [[ ${c} -lt ${v} ]]; then  # Everything is good
    if [[ 2 -gt 1 ]]; then  # For debugging  #ON
        ##############
        ## SIQ PART ##
        ##############
        #loop all by names construct like: "IPfile-INPUTfile-PARAMSfile-outputNAME"
        #you can list these in a single line or not, up to you. but the tic and quote 
        # symbols must be correct! ... tic = ` and quote = "
        z=$(awk '/getTracks/ { print NR + 1 }' "${det_exp}")
        x=$(awk '/getResponse/ { print NR - 1 }' "${det_exp}")
        
        if [[ ${z} -le ${x} ]]; then
            cmd="sed -n ${z},${x}p ${det_exp}"  # echo "${cmd}"
            pare=$(eval "${cmd}" | sed -e 's/ /,/g' | sed -e 's/,$//g')  #DONE: Change delimiter from dash/hyphen to comma or semicolon; file names are far less likely to contain commas or semicolons than dashes/hyphens
            echo "${pare}"
                for w in ${pare} ; do
                    # w="IP1.bed,input1.bed,params1.in,exp1siq"  #TEST
                    ipfi=$(echo "${w}" | sed -e 's/,/ /g'| awk '{ print $1 }')
                    inputfi=$(echo "${w}" | sed -e 's/,/ /g'| awk '{ print $2 }')
                    paramfi=$(echo "${w}" | sed -e 's/,/ /g'| awk '{ print $3 }')
                    namefi=$(echo "${w}" | sed -e 's/,/ /g'| awk '{ print $4 }')
                    echo "bash ./runCrunch.sh ${ipfi} ${inputfi} ${paramfi} ${namefi} ${widths} ${egs}"
                    bash ./runCrunch.sh "${ipfi}" "${inputfi}" "${paramfi}" "${namefi}" "${widths}" "${egs}"
                done
        fi #bounce if empty
    fi
fi

#pairs is a list of files and names, - is used to separate the names. in an entry ONE-TWO-THREE, the files ONE and TWO are compared and the data is written to a file named THREE
if [[ 0 -gt 1 ]]; then # for debugg -toggle on if here --- 0 for off 2 for on  # 2
    #pairs=`echo "SIQLnormK27dmso.bed-SIQLnormK27cbp.bed-Ldms2cbpK27NS SIQLnormK27dmso.bed-SIQLnormK27a485.bed-Ldms2a48K27NS SIQLnormK18dmso.bed-SIQLnormK18cbp.bed-Ldms2cbpK18NS SIQLnormK18dmso.bed-SIQLnormK18a485.bed-Ldms2a48K18NS"`
    #--->dothis SIQLnormK27dmso.bed SIQLnormK27cbp.bed SIQLnormK27a485.bead

    z=$(awk '/getResponse/ { print NR + 1 }' "${det_exp}")
    x=$(awk '/getFracts/ { print NR - 1 }' "${det_exp}")

    if [[ ${z} -le ${x} ]] ; then
        cmd="sed -n ${z},${x}p ${det_exp}"
        pairs=$(eval "${cmd}" | sed -e 's/ /,/g' | sed -e 's/,$//g')
        echo ${pairs}
        for w in ${pairs} ; do
            hifi=$(echo "${w}" | sed -e 's/,/ /g'| awk '{ print $1 }')
            lofi=$(echo "${w}" | sed -e 's/,/ /g'| awk '{ print $2 }')
            filename=$(echo "${w}" | sed -e 's/,/ /g'| awk '{ print $3 }')
            ./WGfrechet.sh "${filename}" "${hifi}" "${lofi}"
            nAn=$(wc -l "${annot}" | awk '{ print $1 }')
            
            if [[ ${nAn} -gt 0 ]]; then
                #something to discuss: +/-50 on center but not what was detected?
                awk '{ print $2, $12 - 50, $13 + 50, $7 }' ${filename} | sort -k1,1 -k2,2n > ${filename}-preanno
                gfortran -O3 -fbounds-check -o readHist.exe binReads.f90
                ./readHist.exe "${filename}-preanno" "${annot}"
                mv matches.coords "${filename}-anno"
                awk '{ print $2, $12 - 50, $13 + 50, $5 }' ${filename} | sort -k1,1 -k2,2n > ${filename}-frechNS-preanno
                gfortran -O3 -fbounds-check -o readHist.exe binReads.f90
                ./readHist.exe "${filename}-frechNS-preanno" "${annot}" 
                #NaN can happen if a peak 'disapear' in the experiment track so we give a 1 to it as this indicates lost info
                sed -i 's/NaN/1/g' matches.coords
                mv matches.coords "${filename}-frechNS-anno"
            fi #only non-empty Annotations
        done
    fi #bounce if empty
fi

##########
##########
if [[ 0 -gt 1 ]]; then #for debugg  # 2
    #compute fractional composition, the user has to adapt scripts to get siq and mass heatmaps
    z=$(grep -n getFracts ${det_exp} | sed -e 's/:/ /g' | awk '{ print $1 + 1 }')
    x=$(grep -n END ${det_exp} | sed -e 's/:/ /g' | awk '{ print $1 - 1 }')
    if [[ $z -le $x ]]; then
        nAn=$(wc -l Annotations.bed | awk '{ print $1 }')
        if [[ ${nAn} -gt 0 ]]; then
            cmd=$(echo "sed -n ${z},${x}p ${det_exp}")
            spare=$(eval "${cmd}" | sed -e 's/ /,/g' | sed -e 's/,$//g')
            
            #define a name file here like all others, loop over sets.
            #the code below will work, but update to catch out last file as output!!!!
            for w in ${spare} ; do 
                files=$(echo "${w}" | sed -e 's/,/ /g')
                declare -a list
                list=( $(echo ${files}) )
                nw=$(echo ${#list[@]})
                listed=$(
                    for (( i=0; i < nw; i++)); do
                        if [[ ${i} -lt $(( nw - 1 )) ]]; then
                            echo "${list[i]}"
                        else
                            echo "> ${list[i]}"
                        fi
                    done
                )
                #echo $listed
                cmd="./makeFracts.sh ${listed}"
                eval "${cmd} "
                sed -E -i 's/^[0-9][0-9]_|^[0-9]_//g' ${list[$((nw - 1))]}
            done   ##########How to deal with the naming situation?
        fi #bounce if Annotations.bed is empty
    fi #bounce if empty
fi #full make fractions

# else
# echo "Your experiment details file is out of order.";
# fi
# else 
# echo "You have no experiment details file."
# fi
#
# ## Shoulder plots
# #gfortran -O3 -fbounds-check -o shoulder.exe shoulders.f90
# else
# echo "I don't see your Annotations.bed file here.
#  Please link to an empty file if you don't have any annotations. 
#  For example, run 'touch Annotations.bed' and try again. Your data will be built without annotations."
# fi
