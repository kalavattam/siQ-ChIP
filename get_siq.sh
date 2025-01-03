#!/bin/bash

#  get_siq.sh
#  BD

#  Run script in interactive/test mode (true) or command-line mode (false)
interactive=false

#  Exit on errors, unset variables, or pipe failures if not in "interactive
#+ mode"
if ! ${interactive}; then set -euo pipefail; fi

#  Set the path to the directory containing scripts
if ${interactive}; then
    ## WARNING: If interactive=true, change path as needed ##
    dir_scr="${HOME}/repos/siQ-ChIP"
else
    dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

#TODO: Add 'dry_run' code throughout script


#  Define functions ===========================================================
function check_set_flag() {
    local name="${1}"
    local value="${2}"

    if [[ -z "${value}" ]]; then
        echo "Error: Flag --${name} is not set." >&2
        return 1
    fi

    case "${value}" in
        true|false) return 0 ;;
        *)
            echo \
                "Error: Flag --${name} must be Boolean 'true' or 'false'," \
                "but got '${value}'." >&2
            return 1
            ;;
    esac
}


function check_supplied_arg() {
    local name="${1}"
    local asgn="${2}"

    if [[ -z "${asgn}" ]]; then
        echo "Error: --${name} is required." >&2
        return 1
    fi
}


function check_exists_file() {
    local name="${1}"
    local item="${2}"

    if [[ ! -f "${item}" ]]; then
        echo \
            "Error: File associated with --${name} does not exist:" \
            "'${item}'." >&2
        return 1
    fi
}


function check_nonempty_file() {
    local name="${1}"
    local item="${2}"

    if [[ ! -s "${item}" ]]; then
        echo \
            "Error: File associated with --${name} is empty: '${item}'." >&2
        return 1
    fi
}


function check_int_nonneg() {
    local name="${1}"
    local value="${2}"

    #  Return error message if value is empty
    if [[ -z "${value}" ]]; then
        echo \
            "Error: --${name} is empty but must be assigned a non-negative" \
            "integer." >&2
        return 1
    fi

    #  Validate that value is a non-negative integer
    if ! [[ "${value}" =~ ^[0-9]+$ ]]; then    
        echo \
            "Error: --${name} was assigned '${value}' but must be a" \
            "non-negative integer." >&2
        return 1
    fi
}


function check_exists_dir() {
    local name="${1}"
    local item="${2}"

    if ! [[ -d "${item}" ]]; then
        echo \
            "Error: Directory associated with --${name} does not exist:" \
            "'${item}'." >&2
        return 1
    fi
}


#  Function to compile Fortran scripts if necessary
function compile_fortran() {
    local fil_src="${1}"  # Path to Fortran source file
    local fil_bin="${2}"  # Path to output binary file

    if [[ ! -f "${fil_bin}" || "${fil_src}" -nt "${fil_bin}" ]]; then
        if ${dry_run}; then
            echo "Would compile: ${fil_src} -> ${fil_bin}"
        else
            gfortran -O3 -fbounds-check -o "${fil_bin}" "${fil_src}"
        fi
    fi
}


#  Initialize argument variables, parse and check arguments, etc. =============
#  Initialize argument variables, assigning default values where applicable
verbose=false
dry_run=false
exp_lay="EXPlayout"
bed_ano=""
siz_bin=30
siz_gen=12157105
dir_out=""
raw=false

#  Assign variable for help message
show_help=$(cat << EOM
Usage:
  get_siq.sh
    [--verbose] [--dry_run] --exp_lay <str> [--bed_ano <str>] --siz_bin <int>
    --siz_gen <int> --dir_out <str> [--raw]

Description:
  #TODO

Arguments:
   -h, --help     Print this help message and exit.
   -v, --verbose  Run script in 'verbose mode' (optional).
  -dr, --dry_run  Run script in 'dry-run mode' (optional).
  -el, --exp_lay  Experiment layout file specifying BED IP, BED input, and TXT
                  parameter files for generating siQ-ChIP tracks, and any track
                  comparisons (default: ${exp_lay}).
  -ba, --bed_ano  BED file of annotations for analyzing fragment distributions
                  from paired-end alignments (optional).
  -bs, --siz_bin  Bin size in bp for computing coverage (default: ${siz_bin}).
  -sg, --siz_gen  Effective genome size for model organism (default:
                  ${siz_gen}).
  -do, --dir_out  The directory to write various outfiles.
   -r, --raw      Write intermediate non-normalized compressed bedGraph files.

Dependencies:
  - AWK
  - Bash or Zsh
  - bc (Bash Calculator)
  - BEDTools
  - GFortran (GNU Fortran)
  - Gnuplot

Notes:
  - BED infiles (IP, input, annotations, etc.) must be coordinate-sorted, e.g.,
    via 'sort -k1,1 -k2,2n'.
  - If any multi-mapping alignments are included in IP and input BED files,
    compute effective genome size with 'faCount' (e.g., 12157105 for S.
    cerevisiae R64 genome).
  - If multi-mapping alignments are excluded, compute with the 'khmer' script
    'unique-kmers.py' [e.g., 11624332 for S. cerevisae R64 genome (50-mers)].
  - ...

Example:
  \`\`\`
  dir_bas="\${HOME}/repos"  ## WARNING: Change as needed ##
  dir_rep="\${dir_bas}/protocol_chipseq_signal_norm"
  dir_dat="\${dir_rep}/data"
  dir_cvg="\${dir_dat}/processed/compute_coverage/bowtie2_global_flag-2_mapq-1"
  dir_exp="\${dir_cvg}/siqchip"
  dir_doc="\${dir_exp}/docs"

  pth_exp="\${dir_doc}/EXPlayout"
  pth_ano="\${dir_doc}/Annotations.bed"
  siz_bin=30
  siz_gen=12157105

  bash "\${HOME}/path/to/repos/siQ-ChIP/get_siq.sh"
      --verbose
      --exp_lay "\${pth_exp}"
      --bed_ano "\${pth_ano}"
      --siz_bin \${siz_bin}
      --siz_gen \${siz_gen}
      --dir_out "\${dir_exp}"
  \`\`\`
EOM
)

#  Parse arguments
if [[ -z "${1:-}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    echo "${show_help}"
    exit 0
fi

while [[ "$#" -gt 0 ]]; do
    case "${1}" in
         -v|--verbose) verbose=true;   shift 1 ;;
        -dr|--dry_run) dry_run=true;   shift 1 ;;
        -el|--exp_lay) exp_lay="${2}"; shift 2 ;;
        -ba|--bed_ano) bed_ano="${2}"; shift 2 ;;
        -bs|--siz_bin) siz_bin="${2}"; shift 2 ;;
        -sg|--siz_gen) siz_gen="${2}"; shift 2 ;;
        -do|--dir_out) dir_out="${2}"; shift 2 ;;
        -io|--raw) raw=true;   shift 1 ;;
        *)
            echo "## Unknown parameter passed: ${1} ##" >&2
            echo "" >&2
            echo "${show_help}" >&2
            exit 1
            ;;
    esac
done

#  Check arguments
check_set_flag      "verbose" ${verbose}

check_set_flag      "dry_run" ${dry_run}

check_supplied_arg  "exp_lay" "${exp_lay}"
check_exists_file   "exp_lay" "${exp_lay}"
check_nonempty_file "exp_lay" "${exp_lay}"

if [[ -n "${bed_ano}" ]]; then
    if \
           check_exists_file   "bed_ano" "${bed_ano}" \
        && check_nonempty_file "bed_ano" "${bed_ano}"
    then
        has_ano=true
    else
        exit 1
    fi
else
    echo "##########"
    echo "## Note ##"
    echo "##########"
    echo ""
    echo \
        "A BED file of annotations for analyzing fragment distributions was" \
        "not provided; will not run computations for response or fractional" \
        "composition."
    echo ""
    echo ""

    has_ano=false
fi

check_supplied_arg "siz_bin" "${siz_bin}"
check_int_nonneg   "siz_bin" "${siz_bin}"

check_supplied_arg "siz_gen" "${siz_gen}"
check_int_nonneg   "siz_gen" "${siz_gen}"

check_supplied_arg "dir_out" "${dir_out}"
check_exists_dir   "dir_out" "${dir_out}"

check_set_flag     "raw"     ${raw}


#  Do the main work ===========================================================
if ${verbose}; then
    echo "##################################################"
    echo "## Argument variable assignments for get_siq.sh ##"
    echo "##################################################"
    echo ""
    echo "verbose=${verbose}"
    echo "dry_run=${dry_run}"
    echo "exp_lay=${exp_lay}"
    echo "bed_ano=${bed_ano:-#N/A}"
    echo "siz_bin=${siz_bin}"
    echo "siz_gen=${siz_gen}"
    echo "dir_out=${dir_out}"
    echo ""
    echo ""
fi

if ${verbose}; then
    echo "###################"
    echo "## Proceeding... ##"
    echo "###################"
    echo ""
    echo ""
fi

z=$(awk '/getTracks/ { print NR }' "${exp_lay}")
x=$(awk '/getResponse/ { print NR }' "${exp_lay}")
c=$(awk '/getFracts/ { print NR }' "${exp_lay}")
v=$(awk '/END/ { print NR }' "${exp_lay}")
if [[ ${z} -lt ${x} && ${x} -lt ${c} && ${c} -lt ${v} ]]; then  # Everything is good
    if [[ 2 -gt 1 ]]; then  # For debugging  #0  #ON
        ##############
        ## SIQ PART ##
        ##############
        #loop all by names construct like: "IPfile-INPUTfile-PARAMSfile-outputNAME"
        #you can list these in a single line or not, up to you. but the tic and quote 
        # symbols must be correct! ... tic = ` and quote = "
        z=$(awk '/getTracks/ { print NR + 1 }' "${exp_lay}")
        x=$(awk '/getResponse/ { print NR - 1 }' "${exp_lay}")
        
        if [[ ${z} -le ${x} ]]; then  # Skip if empty
            cmd="sed -n ${z},${x}p ${exp_lay}"  # echo "${cmd}"
            grp_smp=$(eval "${cmd}" | sed -e 's/ /,/g' | sed -e 's/,$//g')
            
            if ${verbose}; then
                echo "#######################################"
                echo "## Processing experiment layout file ##"
                echo "#######################################"
                echo ""
                echo "${grp_smp}"
                echo ""
                echo ""
            fi
            
            for grp in ${grp_smp}; do
                fil_ip=$(echo "${grp}" | awk -F ',' '{ print $1 }')
                fil_in=$(echo "${grp}" | awk -F ',' '{ print $2 }')
                fil_prm=$(echo "${grp}" | awk -F ',' '{ print $3 }')
                str_stm=$(echo "${grp}" | awk -F ',' '{ print $4 }')
                
                #TODO: Parallelize call to run_crunch.sh with SLURM and GNU Parallel
                # shellcheck disable=SC2046
                bash "${dir_scr}/run_crunch.sh" \
                    $(if ${verbose}; then echo "--verbose"; fi) \
                    --fil_ip  "${fil_ip}" \
                    --fil_in  "${fil_in}" \
                    --fil_prm "${fil_prm}" \
                    --str_stm "${str_stm}" \
                    --siz_bin "${siz_bin}" \
                    --siz_gen "${siz_gen}" \
                    --dir_out "${dir_out}" \
                    $(if ${raw}; then echo "--raw"; fi)
            done
        fi
    fi
fi

#pairs is a list of files and names, - is used to separate the names. in an entry ONE-TWO-THREE, the files ONE and TWO are compared and the data is written to a file named THREE
if ${has_ano} && [[ 0 -gt 1 ]]; then # for debugg -toggle on if here --- 0 for off 2 for on  # 2  #OFF
    #pairs=`echo "SIQLnormK27dmso.bed-SIQLnormK27cbp.bed-Ldms2cbpK27NS SIQLnormK27dmso.bed-SIQLnormK27a485.bed-Ldms2a48K27NS SIQLnormK18dmso.bed-SIQLnormK18cbp.bed-Ldms2cbpK18NS SIQLnormK18dmso.bed-SIQLnormK18a485.bed-Ldms2a48K18NS"`
    #--->dothis SIQLnormK27dmso.bed SIQLnormK27cbp.bed SIQLnormK27a485.bead

    compile_fortran "${dir_scr}/binReads.f90" "${dir_scr}/bin_fragments"

    z=$(awk '/getResponse/ { print NR + 1 }' "${exp_lay}")
    x=$(awk '/getFracts/ { print NR - 1 }' "${exp_lay}")

    if [[ ${z} -le ${x} ]] ; then
        cmd="sed -n ${z},${x}p ${exp_lay}"
        pairs=$(eval "${cmd}" | sed -e 's/ /,/g' | sed -e 's/,$//g')
        echo "${pairs}"
        for w in ${pairs} ; do
            hifi=$(echo "${w}" | sed -e 's/,/ /g' | awk '{ print $1 }')
            lofi=$(echo "${w}" | sed -e 's/,/ /g' | awk '{ print $2 }')
            filename=$(echo "${w}" | sed -e 's/,/ /g' | awk '{ print $3 }')
            bash "${dir_scr}/WGfrechet.sh" "${filename}" "${hifi}" "${lofi}"
            
            nAn=$(awk 'END { print NR }' "${bed_ano}")
            
            if [[ ${nAn} -gt 0 ]]; then
                #something to discuss: +/-50 on center but not what was detected?
                #TODO: Add optional parameters for +/- posint (bp)
                awk '{ print $2, $12 - 50, $13 + 50, $7 }' "${filename}" \
                    | sort -k1,1 -k2,2n \
                        > "${filename}-preanno"
                
                "${dir_scr}/bin_fragments" "${filename}-preanno" "${bed_ano}"
                mv "matches.coords" "${filename}-anno"

                awk '{ print $2, $12 - 50, $13 + 50, $5 }' "${filename}" \
                    | sort -k1,1 -k2,2n \
                        > "${filename}-frechNS-preanno"

                "${dir_scr}/bin_fragments" "${filename}-frechNS-preanno" "${bed_ano}" 
                #NaN can happen if a peak 'disapear' in the experiment track so we give a 1 to it as this indicates lost info
                sed -i 's/NaN/1/g' "matches.coords"
                mv "matches.coords" "${filename}-frechNS-anno"
            fi #only non-empty Annotations
        done
    fi #bounce if empty
fi

##########
##########
if ${has_ano} && [[ 0 -gt 1 ]]; then #for debugg  # 2  #OFF
    #compute fractional composition, the user has to adapt scripts to get siq and mass heatmaps
    z=$(grep -n getFracts ${exp_lay} | sed -e 's/:/ /g' | awk '{ print $1 + 1 }')
    x=$(grep -n END ${exp_lay} | sed -e 's/:/ /g' | awk '{ print $1 - 1 }')
    if [[ ${z} -le ${x} ]]; then
        # nAn=$(wc -l "${bed_ano}" | awk '{ print $1 }')
        nAn=$(awk 'END { print NR }' "${bed_ano}")
        if [[ ${nAn} -gt 0 ]]; then
            cmd="sed -n ${z},${x}p ${exp_lay}"
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
                sed -E -i 's/^[0-9][0-9]_|^[0-9]_//g' ${list[$(( nw - 1 ))]}
            done   ##########How to deal with the naming situation?
        fi #bounce if annotations BED file is empty
    fi #bounce if empty
fi #full make fractions

# ## Shoulder plots
# #gfortran -O3 -fbounds-check -o shoulder.exe shoulders.f90
# else
# echo "I don't see your annotations BED file here.
#  Please link to an empty file if you don't have any bed_anoations. 
#  For example, run 'touch Annotations.bed' and try again. Your data will be built without bed_anoations."
# fi
