#!/bin/bash

#  run_crunch.sh
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


#  Function to write a compressed bedGraph file
function write_bdg_gz() {
    local fil_in="${1}"   # Input data file (e.g., 'IP.data' or 'IN.data')
    local fil_out="${2}"  # Output compressed bedGraph file: '*.bdg.gz'
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
        "Error: Encountered issue writing *.bdg.gz file of ${type}" \
        "coverage." >&2
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
    check_exists_file "fil_src" "${fil_src}"

    check_supplied_arg "fil_out" "${fil_out}"

    check_supplied_arg "typ_cvg" "${typ_cvg}"

    #  If supplied, validate normalization number 
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
            echo "Warning: Could not remove source file '${fil_src}'." >&2
        }
    fi

    return 0
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
    [--verbose] [--dry_run] --fil_ip <str> --fil_in <str> --fil_prm <int>
    --str_stm <str> --siz_bin <int> --siz_gen <int> --dir_out <str> [--raw]

Description:
  This script automates the calculation of coverage signal tracks from sorted
  BED files of fragments from paired-end alignments, returning normalized
  (proportional) coverage, metadata values, siQ-scaled coverage, and optionally
  intermediate non-normalized coverage. The different kinds of coverage are
  written to compressed bedGraph files; the metadata values are written to a
  TXT file.

Arguments:
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
  -sg, --siz_gen  Effective genome size for model organism (default:
                  ${siz_gen}).
  -do, --dir_out  The directory to write various outfiles.
   -r, --raw      Write intermediate non-normalized compressed bedGraph files.

Dependencies:
  - AWK
  - Gzip
  - ...

Notes:
  - IP and input BED files must be coordinate-sorted, e.g., via
    'sort -k1,1 -k2,2n'.
  - #TODO: Note about 'fil_prm'
  - If any multi-mapping alignments are included in IP and input BED files,
    compute effective genome size with 'faCount' (e.g., 12157105 for S.
    cerevisiae R64 genome).
  - If multi-mapping alignments are excluded, compute with the 'khmer' script
    'unique-kmers.py' [e.g., 11624332 for S. cerevisae R64 genome (50-mers)].
  - ...

Example:
  \`\`\`
  echo "Processing experiment layout file:"
  echo "\${grp_smp}"

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

#  Check arguments
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
#TODO: Change f90 script name to snake case
compile_fortran "${dir_scr}/getalpha.f90" "${dir_scr}/get_alpha"
#TODO: Change f90 script name to snake case
compile_fortran "${dir_scr}/mergetracks.f90" "${dir_scr}/merge_tracks"

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
#TODO: Update tracks.f90 to allow specification of dir_out
#TODO: Update tracks.f90 to handle keyword arguments
"${dir_scr}/tracks" \
    "${fil_ip}" "${fil_in}" \
    "${siz_bin}" "${siz_bin}" \
    "${chr}"

#  (Add some whitespace to output)
echo ""
echo ""

#  Compute scaling factor alpha
alf=$("${dir_scr}/get_alpha" "${fil_prm}" "${n_in}" "${n_ip}")

#  Save metadata such as file line counts and alpha to a stem-specific TXT file
if ${verbose}; then
    echo -e \
        "file_ip\tfile_in\tchr_1st\tn_lin_ip\tn_lin_in\tlen_avg_in\tfct_dep\talpha" \
            | tee -a "alpha_etc_${str_stm}.txt"
    echo -e \
        "${fil_ip}\t${fil_in}\t${chr}\t${n_ip}\t${n_in}\t${len_avg}\t${fct_dep}\t${alf}" \
            | tee -a "alpha_etc_${str_stm}.txt"
    echo ""
    echo ""
else
    echo -e \
        "file_ip\tfile_in\tchr_1st\tn_lin_ip\tn_lin_in\tlen_avg_in\tfct_dep\talpha" \
            >> "alpha_etc_${str_stm}.txt"
    echo -e \
        "${fil_ip}\t${fil_in}\t${chr}\t${n_ip}\t${n_in}\t${len_avg}\t${fct_dep}\t${alf}" \
            >> "alpha_etc_${str_stm}.txt"
fi

#  Combine IP and input data using alpha and depth factors
#TODO: Update mergetracks.f90 to allow specification of dir_out
#TODO: Update mergetracks.f90 to handle keyword arguments
"${dir_scr}/merge_tracks" \
    "${dir_scr}/IP.data" \
    "${dir_scr}/IN.data" \
    "${alf}" \
    "${fct_dep}" \
    "${chr}"

#  (Add some whitespace to output)
echo ""
echo ""

#  Write compressed bedGraph of siQ-scaled coverage
write_check_bdg \
    --fil_src "mergedSIQ.data" \
    --fil_out "${dir_out}/siq_${str_stm}.bdg.gz" \
    --typ_cvg "siQ-scaled" \
    --rmv_src

#  Optionally write non-normalized intermediate IP coverage
if ${raw}; then
    write_check_bdg \
        --fil_src "IP.data" \
        --fil_out "${dir_out}/raw_IP_${str_stm}.bdg.gz" \
        --typ_cvg "IP intermediate"
fi

#  Normalize and write IP coverage
write_check_bdg \
    --fil_src "IP.data" \
    --fil_out "${dir_out}/norm_IP_${str_stm}.bdg.gz" \
    --typ_cvg "IP normalized" \
    --n_lines "${n_ip}" \
    --rmv_src

#  Optionally write non-normalized intermeditae input coverage
if ${raw}; then
    write_check_bdg \
        --fil_src "IN.data" \
        --fil_out "${dir_out}/raw_in_${str_stm}.bdg.gz" \
        --typ_cvg "input intermediate"
fi

#  Normalize and write input coverage
write_check_bdg \
    --fil_src "IN.data" \
    --fil_out "${dir_out}/norm_in_${str_stm}.bdg.gz" \
    --typ_cvg "input normalized" \
    --n_lines "${n_in}" \
    --rmv_src
