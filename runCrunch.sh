#  Set the bin width for processing (used for normalization and calculations later).
widths=30

#  Count the number of arguments passed to the script.
narg=`echo $@ | wc | awk '{ print $2 }'`

#  Check if the required number of arguments (4) is provided.
if [ ${narg} -lt 4 ]; then
    # Display usage instructions if fewer than 4 arguments are provided.
    echo 'Run this script as:
    ./runCrunch.sh IPFILE INPUTFILE PARAMfile TAG

    IPFILE: Path to the main (IP) data file.
    INPUTFILE: Path to the control (input) data file.
    PARAMfile: Path to the parameter file for siQ-ChIP measurements.
    TAG: A custom name for the output siQ-ChIP file.

    Note: The path and name of IPFILE and INPUTFILE must each be less than 65 characters.

    Currently, you are missing one or more of these arguments. <----
    '
else 
    #  Assign arguments to variables for clarity and later usage
    ipfile=${1}     # Path to the IP file
    infile=${2}     # Path to the input (control) file
    paramfile=${3}  # Path to the parameter file
    tag=${4}        # Custom name for the output file

    #  Calculate the number of lines (regions) in the IP and input files
    nip=`wc -l ${ipfile} | awk '{ print $1 }'`  # Number of lines in IP file
    nin=`wc -l ${infile} | awk '{ print $1 }'`  # Number of lines in input (control) file

    #  Compute the average region length from the input (control) file (4th
    #+ column values; currently unused)
    legs=`awk '{ sum += $4 } END { print sum/NR }' ${infile}`

    #  Compute the depth factor for input normalization
    dep=`echo $nin*$widths/3200000000./$begin:math:text$1-$widths/3200000000$end:math:text$ | bc -l` 
    echo ${dep}  # Print the computed depth factor for reference

    #  Compile the `tracks.f90` program and execute it; `tracks.f90` processes
    #+ the IP and input files, generating binned coverage files
    gfortran -O3 -fbounds-check tracks.f90            # Compile with optimizations and bounds checking
    ./a.out ${ipfile} ${infile} ${widths} ${widths}   # Execute with input parameters

    #  Compile and run `getalpha.f90` to compute the alpha scaling factor
    gfortran -O3 -fbounds-check getalpha.f90
    a=`./a.out ${paramfile} ${nin} ${nip}`            # Capture alpha value output from 'getalpha.f90'
    echo "i called alpha already" ${nin} ${nip} ${a}  # Print debug message for alpha calculation
    echo ${a} ${nin} ${nip} > $tag.alpha              # Save alpha and file line counts to a tag-specific file

    #  Print the depth factor again for reference
    echo "${dep} is dep"

    #  Compile and run 'mergetracks.f90' to combine the IP and input data using the alpha and depth factors
    gfortran -O3 -fbounds-check mergetracks.f90
    ./a.out IP.data IN.data ${a} ${dep}  # Pass the processed data files and scaling factors
    
    #  Rename the merged output file to a tag-specific name with a '.bed' extension
    mv mergedSIQ.data ${tag}.bed
    echo "you created the file: " ${tag}.bed  # Confirm output creation

    #  Compute normalized coverage for IP and input files 
    nl=`wc -l ${ipfile} | awk '{ print $1 }'`  # Count lines in the IP file
    awk -v var=${nl} '{ print $1,$2,$3,$4 }' IP.data > IntermedIP_${tag}.bed  # Non-normalized intermediate IP coverage
    awk -v var=${nl} '{ print $1,$2,$3,$4/var }' IP.data > NormCovIP_${tag}.bed   # Normalize IP coverage

    nl=`wc -l ${infile} | awk '{ print $1 }'`  # Count lines in the input file
    awk -v var=${nl} '{ print $1,$2,$3,$4 }' IN.data > IntermedIN_${tag}.bed  # Non-normalized intermediate input coverage
    awk -v var=${nl} '{ print $1,$2,$3,$4/var }' IN.data > NormCovIN_${tag}.bed   # Normalize input coverage and save
fi
