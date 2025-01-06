#!/bin/bash

#  run_crunch.sh
#  BD

#  Run script in interactive/test mode (true) or command-line mode (false)
interactive=false

#  Exit on errors, unset variables, or pipe failures if not in "interactive
#+ mode"
if ! ${interactive}; then set -euo pipefail; fi

#  Set path to directory containing scripts
if ${interactive}; then
    ## WARNING: If interactive=true, change path as needed ##
    dir_scr="${HOME}/repos/siQ-ChIP"
    env_nam=""  # Unused in interactive mode
elif [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    dir_scr="${1:-}"
    env_nam="${2:-}"
    if [[ -z "${dir_scr}" ]]; then
        echo \
            "Error: Variable 'dir_scr' must be provided as the first" \
            "argument (positional) when using SLURM." >&2
        exit 1
    fi

    if [[ -z "${env_nam}" ]]; then
        echo \
            "Error: Variable 'env_nam' must be provided as the second" \
            "argument (positional) when using SLURM." >&2
        exit 1
    fi

    shift 2
else
    dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    env_nam=""  # Unused in non-SLURM mode
fi

#  Validate path to directory containing scripts
if [[ ! -d "${dir_scr}" ]]; then
    echo "Error: Script directory is invalid: '${dir_scr}'." >&2
    exit 1
fi

#TODO: Add 'dry_run' code throughout script
#MAYBE: Have conda not be a dependency of siQ-ChIP, although it only is when
#       SLURM is used...
#INPROGRESS: Change f90 script names to snake case
#TODO: Starting at line 650, if verbose=true, then print that files are being
#      written
#TODO: Some kind of check for max_job !> no. threads when '--para gnu'


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
        if ${dry_run:-false}; then
            echo "Dry run: Would compile '${fil_src}' -> '${fil_bin}'."
        else
            if ! \
                gfortran -O3 -fbounds-check -o "${fil_bin}" "${fil_src}"
            then
                echo "Error: Failed to compile '${fil_src}'." >&2
                return 1
            fi
        fi
    fi
}


#  Function to write a compressed bedGraph file
function write_bdg_gz() {
    local fil_in="${1}"     # Input data file (e.g., 'IP.data' or 'IN.data')
    local fil_out="${2}"    # Output compressed bedGraph file: '*.bdg.gz'
    local fct_nrm="${3:-}"  # Norm. factor: Empty or 1 for no norm.

    #  Use awk to process the file
    # shellcheck disable=SC2002
    cat "${fil_in}" \
        | awk \
            -v OFS="\t" -v nl="${fct_nrm}" '
            {
                if (nl != "") { print $1, $2, $3, $4 / nl }
                else { print $1, $2, $3, $4 }
            }
        '  \
        | gzip \
            > "${fil_out}"
}


#  Function to print error for writing compressed bedGraph coverage files
function print_err_cvg() {
    local type="${1}"
    echo \
        "Error with no exit: Encountered issue writing *.bdg.gz file of" \
        "${type} coverage." >&2
}


#  Helper function to write compressed bedGraph file, handle errors, and
#+ optionally clean up data infile
function write_check_bdg() {
    local fil_src=""
    local fil_out=""
    local typ_cvg=""
    local n_lines=""
    local rmv_src=false

    #  Parse keyword parameters
    while [[ $# -gt 0 ]]; do
        case "$1" in
             -s|--fil_src) fil_src="${2}"; shift 2 ;;
             -o|--fil_out) fil_out="${2}"; shift 2 ;;
            -tc|--typ_cvg) typ_cvg="${2}"; shift 2 ;;
            -nl|--n_lines) n_lines="${2}"; shift 2 ;;
            -rs|--rmv_src) rmv_src=true;   shift 1 ;;
            *)
                echo "## Unknown parameter passed: ${1} ##" >&2
                return 1
                ;;
        esac
    done

    #  Validate required parameters
    check_supplied_arg "fil_src" "${fil_src}"
    check_exists_file  "fil_src" "${fil_src}"

    check_supplied_arg "fil_out" "${fil_out}"

    check_supplied_arg "typ_cvg" "${typ_cvg}"

    #  If supplied, validate denominator for normalization 
    if [[ -n "${n_lines}" ]]; then check_int_nonneg "n_lines" "${n_lines}"; fi

    #  Write compressed bedGraph file
    # shellcheck disable=SC2086
    if ! write_bdg_gz "${fil_src}" "${fil_out}" ${n_lines}; then
        print_err_cvg "${typ_cvg}"
        return 1
    fi

    #  Optionally remove source file if compressed bedGraph file was
    #+ successfully written
    if ${rmv_src} && [[ -f "${fil_out}" ]]; then
        rm "${fil_src}" || {
            echo \
                "Warning: --rmv_src was specified but could not remove" \
                "source file '${fil_src}'." >&2
        }
    fi

    return 0
}


function check_unity() {
    local fil_in="${1}"
    local rng_gt="${2:-0.98}"
    local rng_lt="${3:-1.02}"
    local reader

    #  Determine whether to read file with zcat or cat
    if [[ "${fil_in}" == *.gz ]]; then
        reader="zcat"
    else
        reader="cat"
    fi

    #  Process the file and check its sum
    ${reader} "${fil_in}" \
        | awk -v rng_gt="${rng_gt}" -v rng_lt="${rng_lt}" '{
            sum += $4
        } END {
            if (sum >= rng_gt && sum <= rng_lt) {
                print "File sums to approximately unity:", sum
            } else {
                print "File does not sum to unity:", sum
            }
        }'
}


#  Initialize argument variables, parse and check arguments, etc. =============
#  Initialize argument variables, assigning default values where applicable
verbose=false
dry_run=false
fil_ip=""
fil_in=""
fil_prm=""
str_stm=""
siz_bin=30
siz_gen=12157105
dir_out=""
raw=false

#  Assign variable for help message
show_help=$(cat << EOM
Usage:
  run_crunch.sh
    [dir_scr] [env_nam] [--verbose] [--dry_run] --fil_ip <str> --fil_in <str>
    --fil_prm <int> --str_stm <str> --siz_bin <int> --siz_gen <int>
    --dir_out <str> [--raw]

Description:
  This script automates the calculation of coverage signal tracks from sorted
  BED files of fragments from paired-end alignments, returning normalized
  (proportional) coverage, metadata values, siQ-scaled coverage, and optionally
  intermediate non-normalized (raw) coverage. The different kinds of coverage
  are written to compressed bedGraph files; metadata values are written to a
  TXT file.

Positional arguments:
  1, dir_scr      Directory containing scripts. When running with SLURM, must
                  be supplied as first positional argument. Omit argument for
                  serial or GNU Parallel execution.
  2, env_nam      Name of Conda environment to activate. When running with
                  SLURM, must be supplied as second positional argument. Omit
                  argument for serial or GNU Parallel execution.

Keyword arguments:
   -h, --help     Print this help message and exit.
   -v, --verbose  Run script in 'verbose mode' (optional).
  -dr, --dry_run  Run script in 'dry-run mode' (optional).
  -fp, --fil_ip   Path to the IP (main) BED file of aligned fragments
                  (chromosome, start, end, length).
  -fn, --fil_in   Path to the input (control) BED file of aligned fragments
                  (chromosome, start, end, length).
  -fr, --fil_prm  Path to the parameter file for siQ-ChIP measurements
                  corresponding to the IP and input BED files.
  -fs, --str_stm  Custom string, or 'stem', for siQ-ChIP output files.
  -bs, --siz_bin  Bin size in bp for computing coverage (default: ${siz_bin}).
  -sg, --siz_gen  Effective genome size of model organism (default:
                  ${siz_gen}).
  -do, --dir_out  The directory to write various outfiles.
   -r, --raw      Write intermediate non-normalized compressed bedGraph files.

Dependencies:
  - AWK
  - Gzip
  - ...

Notes:
  - If using SLURM (e.g., with 'sbatch'), the path to the script directory must
    be provided as the first positional argument, and a valid Conda environment
    must be provided as the second positional argument.
  - IP and input BED files must be coordinate-sorted, e.g., via
    'sort -k1,1 -k2,2n'.
  - #TODO: Notes about 'fil_prm'
  - Notes on effective genome size:
    + If any multi-mapping alignments are included in IP and input BED files,
      compute effective genome size with 'faCount' (e.g., 12157105 for S.
      cerevisiae R64 genome).
    + If multi-mapping alignments are excluded, compute with the 'khmer' script
      'unique-kmers.py' [e.g., 11624332 for S. cerevisae R64 genome (50-mers)].
  - ...
  - If '--raw' is specified, intermediate non-normalized compressed bedGraph
    files are written in addition to compressed bedGraph files of normalized
    and siQ-scaled coverage.

Example:
  \`\`\`
  echo "Processing experiment layout file:"
  echo "\${grp_smp}"

  #  SLURM usage  #TODO: Complete this
  sbatch run_crunch.sh \
    "\${HOME}/path/to/script/directory" \
    "env_siqchip" \
    --fil_ip <fil_ip> \
    --fil_in <fil_in> \
    --fil_prm <fil_prm> \
    --str_stm <str_stm> \
    --siz_bin <siz_bin> \
    --siz_gen <siz_gen> \
    --dir_out <dir_out>

  #  Serial usage
  for grp in \${grp_smp}; do
    fil_ip=\$(echo "\${grp}" | awk -F ',' '{ print \$1 }')
    fil_in=\$(echo "\${grp}" | awk -F ',' '{ print \$2 }')
    fil_prm=\$(echo "\${grp}" | awk -F ',' '{ print \$3 }')
    str_stm=\$(echo "\${grp}" | awk -F ',' '{ print \$4 }')
    bash "\${dir_scr}/run_crunch.sh"
      --fil_ip  "\${fil_ip}"
      --fil_in  "\${fil_in}"
      --fil_prm "\${fil_prm}"
      --str_stm "\${str_stm}"
      --siz_bin "\${siz_bin}"
      --siz_gen "\${siz_gen}"
      --dir_out "\${dir_out}"
      --raw
  done
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
        -fp|--fil_ip)  fil_ip="${2}";  shift 2 ;;
        -fn|--fil_in)  fil_in="${2}";  shift 2 ;;
        -fr|--fil_prm) fil_prm="${2}"; shift 2 ;;
        -fs|--str_stm) str_stm="${2}"; shift 2 ;;
        -bs|--siz_bin) siz_bin="${2}"; shift 2 ;;
        -sg|--siz_gen) siz_gen="${2}"; shift 2 ;;
        -do|--dir_out) dir_out="${2}"; shift 2 ;;
         -r|--raw)     raw=true;       shift 1 ;;
        *)
            echo "## Unknown argument passed: ${1} ##" >&2
            echo "" >&2
            echo "${show_help}" >&2
            exit 1
            ;;
    esac
done

#  Perform additional parsing and validation for SLURM array job submission
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    if [[ -z "${dir_scr}" ]]; then
        echo \
            "Error: Variable 'dir_scr' is empty; 'dir_scr' should be" \
            "assigned the path to the script directory as the first" \
            "positional argument when running with SLURM." >&2
        exit 1
    fi

    #  Validate Conda and supplied environment, then activate environment
    if ! command -v conda &> /dev/null; then
        echo \
            "Error: Conda is not installed or not accessible in the current" \
            "environment." >&2
        exit 1
    fi

    if [[ -n "${env_nam}" ]]; then
        if ! conda env list | grep -q "^${env_nam}\b"; then
            echo \
                "Error: Conda environment '${env_nam}' does not exist." >&2
            exit 1
        fi
    else
        echo \
            "Error: Variable 'env_nam' is empty; please provide a valid" \
            "Conda environment name." >&2
        exit 1
    fi

    if [[ "${CONDA_DEFAULT_ENV}" != "${env_nam}" ]]; then
        eval "$(conda shell.bash hook)"
        if ! conda activate "${env_nam}"; then
            echo \
                "Error: Failed to activate environment '${env_nam}'." >&2
                exit 1
        fi
    fi

    #  Use SLURM_ARRAY_TASK_ID to determine current task
    IFS=',' read -r -a arr_fil_ip <<< "${fil_ip}"
    IFS=',' read -r -a arr_fil_in <<< "${fil_in}"
    IFS=',' read -r -a arr_fil_prm <<< "${fil_prm}"
    IFS=',' read -r -a arr_str_stm <<< "${str_stm}"

    #  Validate arrays
    for arr_nam in arr_fil_ip arr_fil_in arr_fil_prm arr_str_stm; do
        arr_len=$(eval "echo \${#${arr_nam}[@]}")

        if [[ ${arr_len} -eq 0 ]]; then
            echo \
                "Error: Array '${arr_nam}' is empty. Ensure input arguments" \
                "are correct." >&2
            exit 1
        fi
    done

    #  Debug array contents
    if ${verbose}; then
        echo "####################"
        echo "## Array contents ##"
        echo "####################"
        echo ""
        for arr_nam in arr_fil_ip arr_fil_in arr_fil_prm arr_str_stm; do
            echo "Array '${arr_nam}':"
            # Use eval to expand the array and iterate over its contents
            eval "for elem in \${${arr_nam}[@]}; do echo \"  \${elem}\"; done"
            echo ""
        done
        echo ""
        echo ""
    fi

    #  Ensure task index is valid
    idx=$(( SLURM_ARRAY_TASK_ID - 1 ))
    if [[ "${idx}" -ge "${#arr_fil_ip[@]}" ]]; then
        echo \
            "Error: SLURM_ARRAY_TASK_ID (${SLURM_ARRAY_TASK_ID}) exceeds" \
            "length of arrays." >&2
        exit 1
    fi

    #  Assign variable values for given task
    fil_ip="${arr_fil_ip[idx]}"
    fil_in="${arr_fil_in[idx]}"
    fil_prm="${arr_fil_prm[idx]}"
    str_stm="${arr_str_stm[idx]}"

    #  Debug SLURM_ARRAY_TASK_ID-dependent variable assignments 
    if ${verbose}; then
        echo "###########################################################"
        echo "## Variable assignments based on \${SLURM_ARRAY_TASK_ID} ##"
        echo "###########################################################"
        echo ""
        echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
        echo "fil_ip=${fil_ip}"
        echo "fil_in=${fil_in}"
        echo "fil_prm=${fil_prm}"
        echo "str_stm=${str_stm}"
        echo ""
        echo ""
    fi
fi

#  Check keyword arguments and flags
check_set_flag      "verbose" ${verbose}

check_set_flag      "dry_run" ${dry_run}

check_supplied_arg  "fil_ip" "${fil_ip}"
check_exists_file   "fil_ip" "${fil_ip}"
check_nonempty_file "fil_ip" "${fil_ip}"

check_supplied_arg  "fil_in" "${fil_in}"
check_exists_file   "fil_in" "${fil_in}"
check_nonempty_file "fil_in" "${fil_in}"

check_supplied_arg  "fil_prm" "${fil_prm}"
check_exists_file   "fil_prm" "${fil_prm}"
check_nonempty_file "fil_prm" "${fil_prm}"

check_supplied_arg  "str_stm" "${str_stm}"

check_supplied_arg  "siz_bin" "${siz_bin}"
check_int_nonneg    "siz_bin" "${siz_bin}"

check_supplied_arg  "siz_gen" "${siz_gen}"
check_int_nonneg    "siz_gen" "${siz_gen}"

check_supplied_arg  "dir_out" "${dir_out}"
check_exists_dir    "dir_out" "${dir_out}"

check_set_flag      "raw"     ${raw}


#  Do the main work ===========================================================
if ${verbose}; then
    echo "#####################################################"
    echo "## Argument variable assignments for run_crunch.sh ##"
    echo "#####################################################"
    echo ""
    echo "dir_scr=${dir_scr}"
    echo "env_nam=${env_nam:-#N/A}"
    echo "verbose=${verbose}"
    echo "dry_run=${dry_run}"
    echo "fil_ip=${fil_ip}"
    echo "fil_in=${fil_in}"
    echo "fil_prm=${fil_prm}"
    echo "str_stm=${str_stm}"
    echo "siz_bin=${siz_bin}"
    echo "siz_gen=${siz_gen}"
    echo "dir_out=${dir_out}"
    echo "raw=${raw}"
    echo ""
    echo ""
fi

#  If necessary, compile required Fortran binaries 
compile_fortran "${dir_scr}/tracks.f90" "${dir_scr}/tracks"
compile_fortran "${dir_scr}/getalpha.f90" "${dir_scr}/get_alpha"
compile_fortran "${dir_scr}/merge_tracks.f90" "${dir_scr}/merge_tracks"

#  Extract first chromosome of model organism from IP BED file, skipping
#+ headers or empty lines if present
chr=$(awk 'NR > 1 && $1 !~ /^#/ && NF >= 3 { print $1; exit }' "${fil_ip}")

#  Calculate the number of lines (regions) in the IP and input files
n_ip=$(awk 'END { print NR }' "${fil_ip}")  # No. lines in IP BED file
n_in=$(awk 'END { print NR }' "${fil_in}")  # No. lines in input BED file

#  Compute the average region length from the input (control) file (4th
#+ column values; note: variable 'len_avg' is currently unused)
len_avg=$(awk '{ sum += $4 } END { print sum / NR }' "${fil_in}")

#  Compute depth factor for input normalization
fct_dep=$(
    echo \
        "((${n_in} * ${siz_bin}) / ${siz_gen}) /" \
        "(1 - (${siz_bin} / ${siz_gen}))" \
    | bc -l
)

#  Process IP and input files to generate binned coverage files
"${dir_scr}/tracks" \
    --fil_ip="${fil_ip}" \
    --fil_in="${fil_in}" \
    --bin_ip="${siz_bin}" \
    --bin_in="${siz_bin}" \
    --dat_ip="${dir_out}/IP_${str_stm}.data" \
    --dat_in="${dir_out}/in_${str_stm}.data" \
    --avg_in="${dir_out}/in_avg_${str_stm}.data" \
    --chr="${chr}"

#  (Add some whitespace to output)
echo ""
echo ""

#  Compute scaling factor alpha
#TODO: Implement keyword arguments
alf=$("${dir_scr}/get_alpha" "${fil_prm}" "${n_in}" "${n_ip}")

#  Save metadata such as file line counts and alpha to a stem-specific TXT file
if ${verbose}; then
    echo -e \
        "file_ip\tfile_in\tchr_1st\tn_lin_ip\tn_lin_in\tlen_avg_in\tfct_dep\talpha" \
            | tee -a "${dir_out}/meta_${str_stm}.txt"
    echo -e \
        "${fil_ip}\t${fil_in}\t${chr}\t${n_ip}\t${n_in}\t${len_avg}\t${fct_dep}\t${alf}" \
            | tee -a "${dir_out}/meta_${str_stm}.txt"
    echo ""
    echo ""
else
    echo -e \
        "file_ip\tfile_in\tchr_1st\tn_lin_ip\tn_lin_in\tlen_avg_in\tfct_dep\talpha" \
            >> "${dir_out}/meta_${str_stm}.txt"
    echo -e \
        "${fil_ip}\t${fil_in}\t${chr}\t${n_ip}\t${n_in}\t${len_avg}\t${fct_dep}\t${alf}" \
            >> "${dir_out}/meta_${str_stm}.txt"
fi

#  Combine IP and input data using alpha and depth factors
"${dir_scr}/merge_tracks" \
    --fil_ip="${dir_out}/IP_${str_stm}.data" \
    --fil_in="${dir_out}/in_${str_stm}.data" \
    --fil_siq="${dir_out}/merged_siq_${str_stm}.data" \
    --factr="${alf}" \
    --dep="${fct_dep}" \
    --chr="${chr}"

#  (Add some whitespace to output)
echo ""
echo ""

#  Optionally write non-normalized intermediate IP coverage
if ${raw}; then
    write_check_bdg \
        --fil_src "${dir_out}/IP_${str_stm}.data" \
        --fil_out "${dir_out}/raw_IP_${str_stm}.bdg.gz" \
        --typ_cvg "IP intermediate"
fi

#  Optionally write non-normalized intermeditae input coverage
if ${raw}; then
    write_check_bdg \
        --fil_src "${dir_out}/in_${str_stm}.data" \
        --fil_out "${dir_out}/raw_in_${str_stm}.bdg.gz" \
        --typ_cvg "input intermediate"
fi

#  Normalize and write IP coverage
write_check_bdg \
    --fil_src "${dir_out}/IP_${str_stm}.data" \
    --fil_out "${dir_out}/norm_IP_${str_stm}.bdg.gz" \
    --typ_cvg "IP normalized" \
    --n_lines "${n_ip}" \
    --rmv_src

#  Normalize and write input coverage
write_check_bdg \
    --fil_src "${dir_out}/in_${str_stm}.data" \
    --fil_out "${dir_out}/norm_in_${str_stm}.bdg.gz" \
    --typ_cvg "input normalized" \
    --n_lines "${n_in}" \
    --rmv_src

#  Write compressed bedGraph of siQ-scaled coverage
write_check_bdg \
    --fil_src "${dir_out}/merged_siq_${str_stm}.data" \
    --fil_out "${dir_out}/siq_${str_stm}.bdg.gz" \
    --typ_cvg "siQ-scaled" \
    --rmv_src

#  Check that IP and input normalized coverage sum to unity
echo "#########################"
echo "## Check coverage sums ##"
echo "#########################"
echo ""

if ${raw}; then
    echo "## Raw IP ${str_stm} bedGraph file ##"
    check_unity "${dir_out}/raw_IP_${str_stm}.bdg.gz"
    echo ""

    echo "## Raw input ${str_stm} bedGraph file ##"
    check_unity "${dir_out}/raw_in_${str_stm}.bdg.gz"
    echo ""
fi

echo "## Normalized IP ${str_stm} bedGraph file ##"
check_unity "${dir_out}/norm_IP_${str_stm}.bdg.gz"
echo ""

echo "## Normalized input ${str_stm} bedGraph file ##"
check_unity "${dir_out}/norm_in_${str_stm}.bdg.gz"
echo ""

echo "## siQ-scaled ${str_stm} bedGraph file ##"
check_unity "${dir_out}/siq_${str_stm}.bdg.gz"
echo ""
echo ""
