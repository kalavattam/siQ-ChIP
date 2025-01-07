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

#TODO: Implement argument for env_nam (only used when running jobs with SLURM)


#  Define functions ===========================================================
function check_flag() {
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


#  Set paths, values, and parameters for interactive mode
function set_interactive() {
    #  Set hardcoded paths, values, etc.
    ## WARNING: If interactive=true, change values as needed ##
    dir_bas="${HOME}/repos"
    dir_rep="${dir_bas}/protocol_chipseq_signal_norm"
    dir_pro="${dir_rep}/data/processed"
    dir_cvg="${dir_pro}/compute_coverage/bowtie2_global_flag-2_mapq-1"
    dir_exp="${dir_cvg}/siqchip"
    dir_doc="${dir_exp}/docs"

    #  Set hardcoded argument assignments
    verbose=true
    dry_run=false
    exp_lay="${dir_doc}/layout_exp.txt"
    bed_ano=""
    siz_bin=1
    siz_gen=12157105
    dir_out="${dir_exp}"
    raw=true
    submit="serial"
    nam_job="run_crunch"
    max_job=12
    time="0:30:00"
}


#  Initialize argument variables, parse and check arguments, etc. =============
#  Initialize argument variables, assigning default values where applicable
verbose=false
dry_run=false
exp_lay=""
bed_ano=""
siz_bin=30
siz_gen=12157105
dir_out=""
raw=false
submit="serial"
nam_job="run_crunch"
max_job=12
time="0:30:00"

#  Assign variable for help message
show_help=$(cat << EOM
Usage:
  get_siq.sh
    [--verbose] [--dry_run] --exp_lay <str> [--bed_ano <str>] --siz_bin <int>
    --siz_gen <int> --dir_out <str> [--raw] --nam_job <str> --submit <str>
    --max_job <int> --time <str>

Description:
  #TODO

Arguments:
   -h, --help     Print this help message and exit.
   -v, --verbose  Run script in 'verbose mode' (optional).
  -dr, --dry_run  Run script in 'dry-run mode' (optional).
  -el, --exp_lay  Experiment layout file specifying BED IP, BED input, and TXT
                  parameter files for generating siQ-ChIP tracks, and any track
                  comparisons.
  -ba, --bed_ano  BED file of annotations for analyzing fragment distributions
                  from paired-end alignments (optional).
  -bs, --siz_bin  Bin size in bp for computing coverage (default: ${siz_bin}).
  -sg, --siz_gen  Effective genome size for model organism (default:
                  ${siz_gen}).
  -do, --dir_out  The directory to write various outfiles.
   -r, --raw      Write intermediate non-normalized (raw) compressed bedGraph
                  files (optional).
  -nj, --nam_job  Names of jobs (default: "${nam_job}").
   -s, --submit   Specifies the method for submitting and executing jobs;
                  options: "slurm", "gnu", "serial" (default: "${submit}").
                  Details:
                    + "slurm": Submit jobs via a SLURM scheduler.
                    + "gnu": Use GNU Parallel for job execution.
                    + "serial": Run jobs sequentially in a single process.
  -mj, --max_job  The maximum number of jobs to run at one time (required if
                  '--submit slurm' or '--submit gnu' is specified, and ignored
                  if '--submit serial' is specified; default: ${max_job}).
  -tm, --time     The length of time, in 'mm:ss', 'h:mm:ss', or 'hh:mm:ss'
                  format, for the SLURM job (required if '--submit slurm' is
                  specified, ignored if not; default: "${time}").

Dependencies:
  - AWK
  - Bash or Zsh
  - bc (Bash Calculator)
  - BEDTools
  - Conda (if '--submit slurm')
  - GFortran (GNU Fortran)
  - GNU Parallel (if '--submit gnu')
  - Gnuplot
  - SLURM (if '--submit slurm')

Notes:
  - BED infiles (IP, input, annotations, etc.) must be coordinate-sorted, e.g.,
    via 'sort -k1,1 -k2,2n'.
  - Notes on effective genome size:
    + If any multi-mapping alignments are included in IP and input BED files,
      compute effective genome size with 'faCount' (e.g., 12157105 for S.
      cerevisiae R64 genome).
    + If multi-mapping alignments are excluded, compute with the 'khmer' script
      'unique-kmers.py' [e.g., 11624332 for S. cerevisae R64 genome (50-mers)].
  - Notes on calling 'run_crunch.sh' with SLURM ('sbatch'; i.e.,
    '--submit slurm'):
    + Provide the script directory path ("\${dir_scr}") as a first positional
      argument.
    + Provide the name of the Conda environment ("\${env_nam:-env_siqchip}")
      as a second positional argument.
    + These positional arguments are not required for serial execution or when
      using GNU Parallel.
  - ...

Example:
  \`\`\`
  dir_bas="\${HOME}/repos"
  dir_rep="\${dir_bas}/protocol_chipseq_signal_norm"
  dir_dat="\${dir_rep}/data"
  dir_cvg="\${dir_dat}/processed/compute_coverage/bowtie2_global_flag-2_mapq-1"
  dir_exp="\${dir_cvg}/siqchip"
  dir_doc="\${dir_exp}/docs"

  pth_exp="\${dir_doc}/EXPlayout"
  pth_ano="\${dir_doc}/Annotations.bed"
  siz_bin=10
  siz_gen=12157105
  raw=true
  submit="slurm"

  bash "\${HOME}/path/to/repos/siQ-ChIP/get_siq.sh"
      --verbose
      --exp_lay "\${pth_exp}"
      --bed_ano "\${pth_ano}"
      --siz_bin \${siz_bin}
      --siz_gen \${siz_gen}
      --dir_out "\${dir_exp}"
      \$(if \${raw}; then echo "--raw"; fi)
      --submit "\${submit}"
  \`\`\`
EOM
)

#  Parse arguments
if [[ -z "${1:-}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    echo "${show_help}"
    if ! ${interactive}; then exit 0; fi
fi

if ${interactive}; then
    set_interactive
else
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
             -v|--verbose) verbose=true;    shift 1 ;;
            -dr|--dry_run) dry_run=true;    shift 1 ;;
            -el|--exp_lay) exp_lay="${2}";  shift 2 ;;
            -ba|--bed_ano) bed_ano="${2}";  shift 2 ;;
            -bs|--siz_bin) siz_bin="${2}";  shift 2 ;;
            -sg|--siz_gen) siz_gen="${2}";  shift 2 ;;
            -do|--dir_out) dir_out="${2}";  shift 2 ;;
             -r|--raw)     raw=true;        shift 1 ;;
            -nj|--nam_job) nam_job="${2}";  shift 2 ;;
             -s|--submit)  submit="${2,,}"; shift 2 ;;
            -mj|--max_job) max_job="${2}";  shift 2 ;;
            -tm|--time)    time="${2}";     shift 2 ;;
            *)
                echo "## Unknown parameter passed: ${1} ##" >&2
                echo "" >&2
                echo "${show_help}" >&2
                if ! ${interactive}; then exit 1; fi
                ;;
        esac
    done
fi

#  Check arguments
check_flag      "verbose" ${verbose}

check_flag      "dry_run" ${dry_run}

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
    echo "#########################"
    echo "## Note on annotations ##"
    echo "#########################"
    echo ""
    echo \
        "A BED file of annotations for analyzing fragment distributions was" \
        "not provided; will not run computations for response and fractional" \
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

check_flag     "raw"     ${raw}

check_supplied_arg "nam_job" "${nam_job}"

check_supplied_arg "submit"  "${submit}"
case "${submit}" in
    slurm|gnu)
        check_supplied_arg "max_job" "${max_job}"
        check_int_nonneg   "max_job" "${max_job}"

        if [[ "${submit}" == "slurm" ]]; then
            check_supplied_arg "time" "${time}"
            if [[ ! "${time}" =~ ^([0-9]{1,2}:)?[0-5][0-9]:[0-5][0-9]$ ]]; then
                echo \
                    "Error: --time was assigned '${time}', an invalid" \
                    "format. Expected format is 'mm:ss', 'h:mm:ss', or" \
                    "'hh:mm:ss'." >&2
                if ! ${interactive}; then exit 1; fi
            fi
        elif [[ "${submit}" == "gnu" ]]; then
            unset time
        fi
        ;;
    serial)
        unset max_job time
        ;;
    *)
        echo \
            "Error: --submit was assigned '${submit}', but it must be one of" \
            "'slurm', 'gnu', or 'serial'." >&2
        if ! ${interactive}; then exit 1; fi
        ;;
esac


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
    echo "raw=${raw}"
    echo "nam_job=${nam_job}"
    echo "submit=${submit}"
    echo "max_job=${max_job:-#N/A}"
    echo "time=${time:-#N/A}"
    echo ""
    echo ""
fi


#  Compute normalized, raw (optional), and siQ-scaled coverage ----------------
#  Collect line numbers for layout markers
ln_trk=$(awk '/getTracks/ { print NR }' "${exp_lay}")    # echo "${ln_trk}"
ln_rsp=$(awk '/getResponse/ { print NR }' "${exp_lay}")  # echo "${ln_rsp}"
ln_frc=$(awk '/getFracts/ { print NR }' "${exp_lay}")    # echo "${ln_frc}"
ln_end=$(awk '/END/ { print NR }' "${exp_lay}")          # echo "${ln_end}"

#  Validate line numbers for layout markers
for ln_nam in ln_trk ln_rsp ln_frc ln_end; do
    ln_val="${!ln_nam}"
    if [[ -z "${ln_val}" ]]; then
        echo \
            "Error: Layout marker for '${ln_nam}' not found in ${exp_lay}." >&2
        if ! ${interactive}; then exit 1; fi
    fi
    
    if [[ ! ${ln_val} =~ ^[0-9]+$ ]]; then
        echo \
            "Error: Invalid line number for '${ln_nam}' layout marker in" \
            "${exp_lay}." >&2
        if ! ${interactive}; then exit 1; fi
    fi
done

if [[
       ${ln_trk} -lt ${ln_rsp} \
    && ${ln_rsp} -lt ${ln_frc} \
    && ${ln_frc} -lt ${ln_end}
]]; then
    if [[ 2 -gt 1 ]]; then  # Debugging: 0 is "off", 2 is "on"
        ln_trk=$(awk '/getTracks/ { print NR + 1 }' "${exp_lay}")    # echo "${ln_trk}"
        ln_rsp=$(awk '/getResponse/ { print NR - 1 }' "${exp_lay}")  # echo "${ln_rsp}"
        
        #TODO: Modularize processing of layout file, array assignment, etc.
        if [[ ${ln_trk} -le ${ln_rsp} ]]; then
            #  Construct a call to 'sed'
            cmd="sed -n ${ln_trk},${ln_rsp}p ${exp_lay}"  # echo "${cmd}"
            
            #  Extract the processed lines into a variable (already
            #+ comma-separated)
            grp_smp=$(eval "${cmd}" | sed -e 's/ /,/g' -e 's/,$//g')

            #  Check contents of processed experiment layout file
            if ${verbose}; then
                echo "#######################################"
                echo "## Processing experiment layout file ##"
                echo "#######################################"
                echo ""
                echo "${grp_smp}"
                echo ""
                echo ""
            fi

            #  Initialize arrays for grp_smp records and fields
            unset arr_ip  && typeset -a arr_ip
            unset arr_in  && typeset -a arr_in
            unset arr_prm && typeset -a arr_prm
            unset arr_stm && typeset -a arr_stm

            #  Read each record from grp_smp, splitting fields into arrays
            while IFS=, read -r fld_ip fld_in fld_prm fld_stm; do
                if [[
                       -z "${fld_ip}"  || -z "${fld_in}" \
                    || -z "${fld_prm}" || -z "${fld_stm}"
                ]]; then
                    echo \
                        "Error: Malformed line in experiment layout:" \
                        "${fld_ip},${fld_in},${fld_prm},${fld_stm}" >&2
                    if ! ${interactive}; then exit 1; fi
                fi

                arr_ip+=( "${fld_ip}" )
                arr_in+=( "${fld_in}" )
                arr_prm+=( "${fld_prm}" )
                arr_stm+=( "${fld_stm}" )
            done <<< "${grp_smp}"

            #  Validate arrays are populated
            if [[
                  ${#arr_ip[@]}  -eq 0 || ${#arr_in[@]}  -eq 0 \
               || ${#arr_prm[@]} -eq 0 || ${#arr_stm[@]} -eq 0
            ]]; then
                echo \
                    "Error: One or more of arrays 'arr_ip', 'arr_in'," \
                    "'arr_prm', or 'arr_stm' are empty." >&2
                if ! ${interactive}; then exit 1; fi
            fi

            #  Check contents of arrays
            if ${verbose}; then
                echo "########################"
                echo "## Contents of arrays ##"
                echo "########################"
                echo ""
                echo "\${fld_ip[*]}=( ${arr_ip[*]} )"
                echo "\${fld_in[*]}=( ${arr_in[*]} )"
                echo "\${fld_prm[*]}=( ${arr_prm[*]} )"
                echo "\${fld_stm[*]}=( ${arr_stm[*]} )"
                echo ""
                
                for idx in "${!arr_ip[@]}"; do
                    echo \
                        "## Array contents and variable assignments for" \
                        "index ${idx} ##"
                    echo ""
                    echo "fld_ip=${arr_ip[idx]}"
                    echo "fld_in=${arr_in[idx]}"
                    echo "fld_prm=${arr_prm[idx]}"
                    echo "fld_stm=${arr_stm[idx]}"
                    echo ""
                    echo ""
                done
            fi

            case "${submit}" in
                serial)
                    if ${verbose} || ${dry_run}; then
                        echo "#####################################"
                        echo "## Serial call(s) to run_crunch.sh ##"
                        echo "#####################################"
                        echo ""
                    fi
                    
                    for idx in "${!arr_ip[@]}"; do
                        if ${verbose} || ${dry_run}; then
                            echo "bash ${dir_scr}/run_crunch.sh \\"
                            echo "    $(if ${verbose}; then echo "--verbose"; fi) \\"
                            echo "    $(if ${dry_run}; then echo "--dry_run"; fi) \\"
                            echo "    --fil_ip  ${arr_ip[idx]} \\"
                            echo "    --fil_in  ${arr_in[idx]} \\"
                            echo "    --fil_prm ${arr_prm[idx]} \\"
                            echo "    --str_stm ${arr_stm[idx]} \\"
                            echo "    --siz_bin ${siz_bin} \\"
                            echo "    --siz_gen ${siz_gen} \\"
                            echo "    --dir_out ${dir_out} \\"
                            echo "    $(if ${raw}; then echo "--raw"; fi)"
                            echo "         > ${dir_out}/logs/${nam_job}.${arr_stm[idx]}.stdout.txt \\"
                            echo "        2> ${dir_out}/logs/${nam_job}.${arr_stm[idx]}.stderr.txt"
                            echo ""
                        fi

                        if ! ${dry_run}; then
                            #TODO: Test new use of arrays here
                            # shellcheck disable=SC2046
                            bash "${dir_scr}/run_crunch.sh" \
                                $(if ${verbose}; then echo "--verbose"; fi) \
                                $(if ${dry_run}; then echo "--dry_run"; fi) \
                                --fil_ip  "${arr_ip[idx]}" \
                                --fil_in  "${arr_in[idx]}" \
                                --fil_prm "${arr_prm[idx]}" \
                                --str_stm "${arr_stm[idx]}" \
                                --siz_bin "${siz_bin}" \
                                --siz_gen "${siz_gen}" \
                                --dir_out "${dir_out}" \
                                $(if ${raw}; then echo "--raw"; fi) \
                                     > "${dir_out}/logs/${nam_job}.${arr_stm[idx]}.stdout.txt" \
                                    2> "${dir_out}/logs/${nam_job}.${arr_stm[idx]}.stderr.txt"
                        fi
                    done

                    if ${verbose} || ${dry_run}; then echo ""; fi
                    ;;

                slurm)
                    if ${verbose} || ${dry_run}; then
                        echo "##########################"
                        echo "## Call to SLURM sbatch ##"
                        echo "##########################"
                        echo ""
                        echo "sbatch \\"
                        echo "    --job-name=${nam_job} \\"
                        echo "    --nodes=1 \\"
                        echo "    --cpus-per-task=1 \\"
                        echo "    --time=${time} \\"
                        echo "    --output=${dir_out}/logs/${nam_job}.%A-%a.stdout.txt \\"
                        echo "    --error=${dir_out}/logs/${nam_job}.%A-%a.stderr.txt \\"
                        echo "    --array=1-${#arr_ip[@]}%${max_job:-12} \\"
                        echo "    ${dir_scr}/run_crunch.sh \\"
                        echo "        ${dir_scr} \\"
                        echo "        ${env_nam:-env_siqchip} \\"
                        echo "        $(if ${verbose}; then echo "--verbose"; fi) \\"
                        echo "        $(if ${dry_run}; then echo "--dry_run"; fi) \\"
                        echo "        --fil_ip $(echo "${arr_ip[*]}" | tr ' ' ',') \\"
                        echo "        --fil_in $(echo "${arr_in[*]}" | tr ' ' ',') \\"
                        echo "        --fil_prm $(echo "${arr_prm[*]}" | tr ' ' ',') \\"
                        echo "        --str_stm $(echo "${arr_stm[*]}" | tr ' ' ',') \\"
                        echo "        --siz_bin ${siz_bin} \\"
                        echo "        --siz_gen ${siz_gen} \\"
                        echo "        --dir_out ${dir_out} \\"
                        echo "        $(if ${raw}; then echo "--raw"; fi)"
                        echo ""
                        echo ""
                    fi
                    
                    if ! ${dry_run}; then
                        # shellcheck disable=SC2046,SC2086
                        sbatch \
                            --job-name=${nam_job} \
                            --nodes=1 \
                            --cpus-per-task=1 \
                            --time=${time} \
                            --output=${dir_out}/logs/${nam_job}.%A-%a.stdout.txt \
                            --error=${dir_out}/logs/${nam_job}.%A-%a.stderr.txt \
                            --array=1-${#arr_ip[@]}%${max_job:-12} \
                            ${dir_scr}/run_crunch.sh \
                                ${dir_scr} \
                                ${env_nam:-env_siqchip} \
                                $(if ${verbose}; then echo "--verbose"; fi) \
                                $(if ${dry_run}; then echo "--dry_run"; fi) \
                                --fil_ip $(echo "${arr_ip[*]}" | tr ' ' ',') \
                                --fil_in $(echo "${arr_in[*]}" | tr ' ' ',') \
                                --fil_prm $(echo "${arr_prm[*]}" | tr ' ' ',') \
                                --str_stm $(echo "${arr_stm[*]}" | tr ' ' ',') \
                                --siz_bin ${siz_bin} \
                                --siz_gen ${siz_gen} \
                                --dir_out ${dir_out} \
                                $(if ${raw}; then echo "--raw"; fi)
                    fi
                    ;;

                gnu)
                    #  Construct base command for GNU Parallel
                    cmd="bash ${dir_scr}/run_crunch.sh"
                    
                    #  Dynamically add optional flags to command
                    if ${verbose}; then cmd+=" --verbose"; fi
                    if ${dry_run}; then cmd+=" --dry_run"; fi
                    if ${raw}; then cmd+=" --raw"; fi

                    #  Add argument placeholders
                    cmd+=" --fil_ip {1} --fil_in {2} --fil_prm {3}"
                    cmd+=" --str_stm {4} --siz_bin {5} --siz_gen {6}"
                    cmd+=" --dir_out {7}"

                    if ${dry_run} || ${verbose}; then
                        {
                            echo "#############################"
                            echo "## Call(s) to GNU Parallel ##"
                            echo "#############################"
                            echo ""
                            parallel --colsep ' ' --jobs "${max_job:-12}" --dry-run \
                                "${cmd}" \
                                ::: "${arr_ip[@]}" \
                                :::+ "${arr_in[@]}" \
                                :::+ "${arr_prm[@]}" \
                                :::+ "${arr_stm[@]}" \
                                ::: "${siz_bin}" \
                                ::: "${siz_gen}" \
                                ::: "${dir_out}"
                            echo ""
                        } \
                             > >(tee -a "${dir_out}/logs/${nam_job}.stdout.txt") \
                            2> >(tee -a "${dir_out}/logs/${nam_job}.stderr.txt")
                    fi
                    #TODO: Troubleshoot errors saying std{out,err} files don't exist
                    #      Prob b/c forgot to 'mkdir logs' in $dir_exp...

                    if ! ${dry_run}; then
                        parallel --colsep ' ' --jobs "${max_job:-12}" \
                            "${cmd}" \
                            ::: "${arr_ip[@]}" \
                            :::+ "${arr_in[@]}" \
                            :::+ "${arr_prm[@]}" \
                            :::+ "${arr_stm[@]}" \
                            ::: "${siz_bin}" \
                            ::: "${siz_gen}" \
                            ::: "${dir_out}"
                    fi
                    #TODO: Need to handle writing of std{out,err} files; otherwise, works!
                    ;;
            esac
            # unset ln_trk ln_rsp ln_frc ln_end cmd grp_smp
        fi
    fi
fi

#pairs is a list of files and names, - is used to separate the names. in an entry ONE-TWO-THREE, the files ONE and TWO are compared and the data is written to a file named THREE
if ${has_ano} && [[ 0 -gt 1 ]]; then # for debugg -toggle on if here --- 0 for off 2 for on  # 2  #OFF
    #pairs=`echo "SIQLnormK27dmso.bed-SIQLnormK27cbp.bed-Ldms2cbpK27NS SIQLnormK27dmso.bed-SIQLnormK27a485.bed-Ldms2a48K27NS SIQLnormK18dmso.bed-SIQLnormK18cbp.bed-Ldms2cbpK18NS SIQLnormK18dmso.bed-SIQLnormK18a485.bed-Ldms2a48K18NS"`
    #--->dothis SIQLnormK27dmso.bed SIQLnormK27cbp.bed SIQLnormK27a485.bead

    compile_fortran "${dir_scr}/binReads.f90" "${dir_scr}/bin_fragments"

    #TODO: Give variables clearer names
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
    #TODO: Give variables clearer names
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
