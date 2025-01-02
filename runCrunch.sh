#!/bin/bash

#  Count the number of arguments passed to the script
narg=$(echo $@ | wc | awk '{ print $2 }')  #TODO: Parse keyword parameters

#TODO: Make the writing of intermediate coverage files optional, specified by some flag
#  Check if the required number of arguments (4) is provided
if [[ ${narg} -lt 6 ]]; then  #TODO: Parse keyword parameters
    #  Display usage instructions if fewer than 4 arguments are provided
    echo "Run this script as:
    ./runCrunch.sh IPFILE INPUTFILE PARAMFILE TAG WIDTHS EGS

    IPFILE: Path to the main (IP) BED file of aligned fragments (chr, start,
            end, length).
    INPUTFILE: Path to the control (input) BED file of aligned fragments (chr,
               start, end, length).
    PARAMFILE: Path to the parameter file for siQ-ChIP measurements.
    TAG: Custom name for the output siQ-ChIP file.
    WIDTHS: Bin widths for processing IP and INPUT files.
    EGS: Effective genome size to compute depth factor for input normalization.
      - If any multi-mapping alignments are included in IP and input BED files,
        compute with 'faCount' (e.g., 12157105 for S. cerevisiae)
      - If multi-mapping alignments are excluded, compute with the 'khmer'
        script 'unique-kmers.py' [e.g., 11624332 for S. cerevisae (50-mers)]

    Note: The path and name of IPFILE and INPUTFILE must each be less than 65
          characters.

    Currently, you are missing one or more of these arguments.
    "
else 
    #  Assign arguments to variables for clarity and later usage
    ipfile=${1}         # Path to IP file
    infile=${2}         # Path to input (control) file
    paramfile=${3}      # Path to parameter file
    tag=${4}            # Custom name for output file
    widths=${5:-30}     # Set bin width for processing
    egs=${6:-12157105}  # Effective genome size for model organism (e.g.,
                        # '12157105' for S. cerevisiae including multi-mapping
                        # alignments)

    #NOTE: Since we parse positional parameters, all parameters must be
    #      explicitly specifed by users
    #TODO: Refactor driver script to parse keyword parameters

    #  Determine first chromosome of model organism, extracting it from sorted
    #+ IP BED file
    rchr=$(awk 'NR == 2 { print $1 }' "${ipfile}")

    #  Calculate the number of lines (regions) in the IP and input files
    nip=$(wc -l "${ipfile}" | awk '{ print $1 }')  # Number of lines in IP file
    nin=$(wc -l "${infile}" | awk '{ print $1 }')  # Number of lines in input (control) file

    #  Compute the average region length from the input (control) file (4th
    #+ column values; note: variable 'legs' is currently unused)
    legs=$(awk '{ sum += $4 } END { print sum / NR }' "${infile}")

    #  Compute depth factor for input normalization
    dep=$(  #NOTE: Works with both Bash and Zsh now
        echo "((${nin} * ${widths}) / ${egs}) / (1 - (${widths} / ${egs}))" | bc -l
    )
    echo "${dep}"  # Print the computed depth factor for reference

    #  Compile the `tracks.f90` program and execute it; `tracks.f90` processes
    #+ the IP and input files, generating binned coverage files
    gfortran -O3 -fbounds-check tracks.f90                             
    ./a.out "${ipfile}" "${infile}" "${widths}" "${widths}" "${rchr}"  # Execute with input parameters

    #  Compile and run `getalpha.f90` to compute the alpha scaling factor
    gfortran -O3 -fbounds-check getalpha.f90
    a=$(./a.out "${paramfile}" "${nin}" "${nip}")               # Capture alpha value output from 'getalpha.f90'
    echo "I called alpha already: ${nin} ${nip} ${a}"           # Print debug message for alpha calculation
    echo "alpha"$'\t'"line_IP"$'\t'"line_in" > "alpha_${tag}.txt"
    echo "${a}"$'\t'"${nip}"$'\t'"${nin}" >> "alpha_${tag}.txt"  # Save alpha and file line counts to a tag-specific file

    #  Print the depth factor again for reference
    echo "${dep} is dep"

    #  Compile and run 'mergetracks.f90' to combine IP and input data using
    #+ alpha and depth factors
    gfortran -O3 -fbounds-check mergetracks.f90        # Compile with optimizations and bounds-checking
    ./a.out IP.data IN.data "${a}" "${dep}" "${rchr}"  # Pass the processed data files, scaling factor, and model first chromosome
    
    #  Give merged outfile tag-specific name with '.bdg.gz' extension
    # mv mergedSIQ.data "siq_${tag}.bed"  #INPROGRESS: Make it into compressed bedGraph, and give it appropriate name  #DEKHO
    # shellcheck disable=SC2002
    if ! \
        cat mergedSIQ.data \
            | awk -v OFS="\t" '{ print $1, $2, $3, $4 }' \
            | gzip \
                > "siq_${tag}.bdg.gz"
    then
        echo \
            "Error: Encountered issue writing bedGraph file of siQ-scaled" \
            "coverage" >&2
    else
        echo "You created the following file: siq_${tag}.bed"   # Confirm output creation
    fi

    #  Compute normalized coverage for IP and input files
    nl_ip=$(wc -l < "${ipfile}")  # Count lines in IP BED file
    nl_in=$(wc -l < "${infile}")  # Count lines in input BED file
    
    # shellcheck disable=SC2002
    {  #TODO: Modularize bedGraph writing
        #  Write compressed bedGraph of non-normalized intermediate IP coverage
        cat IP.data \
            | awk -v OFS="\t" '{ print $1, $2, $3, $4 }' \
            | gzip \
                > "inter_IP_${tag}.bdg.gz"

        #  Write compressed bedGraph of normalized IP coverage
        cat IP.data \
            | awk \
                -v OFS="\t" -v nl="${nl_ip}" \
                '{ print $1, $2, $3, $4 / nl }' \
            | gzip \
                > "norm_IP_${tag}.bdg.gz"

        #  Write compressed bedGraph of non-normalized intermeditae input
        #+ coverage
        cat IN.data \
            | awk -v OFS="\t" '{ print $1, $2, $3, $4 }' \
            | gzip \
                > "inter_in_${tag}.bdg.gz"

        #  Write compressed bedGraph of normalized input coverage
        cat IN.data \
            | awk \
                -v OFS="\t" -v nl="${nl_in}" \
                '{ print $1, $2, $3, $4 / nl }' \
            | gzip \
                > "norm_in_${tag}.bdg.gz"
    }
fi
