program tracks

   ! Rough documentation:
   !   'tracks' processes IP and input BED files and generates binned output
   !   data.

   ! Declare variables
   real(8), dimension(100000) :: bin
   real*8 :: sum, itops
   ! bin: Array for binning data
   ! sum: Cumulative sum
   ! itops: Max bin value
   integer nleng(2), ilft, irght, leng, leftend, rightend, Reason, lpos
   ! nleng(2): Array for bin sizes of the two input files
   ! ilft, irght: Start and end positions of the genomic region being processed
   ! leng: Length of the region (irght - ilft)
   ! leftend, rightend: Boundaries of the current binning range
   ! Reason: I/O status code after reading input
   ! lpos: Tracks current genomic position in the bin array
   character(len = 6) :: unk, rchr                                                                ! Chromosome identifiers
   character(len = 5) :: unk2, rchr2                                                              ! Alternate chromosome identifiers (unused)
   character(len = 44) :: siq                                                                     ! Placeholder for HMM annotations (unused)
   character(len = 256) :: arg, path(2)                                                            ! Infile paths
   logical :: file_exists                                                                         ! File existence flag
   integer :: inot, iread
   character(len = 256) :: fil_ip, fil_in, dat_ip, dat_in, avg_in                                 ! File paths
   integer :: bin_ip, bin_in                                                                      ! Bin sizes

   ! Assignments by Brad...
   inot = 0                                                                                       ! Integer variable initialized to 0 (unused)
   itops = 0                                                                                      ! Tracks maximum value in bin array
   iread = 0                                                                                      ! Integer variable initialized to 0 (unused)

   ! Call argument parser
   call parse_args(fil_ip, fil_in, bin_ip, bin_in, dat_ip, dat_in, avg_in, rchr)

   ! Debug or verify arguments (optional)
   write(*, *) "#################################"
   write(*, *) "## Parsed arguments for tracks ##"
   write(*, *) "#################################"
   write(*, *) ""
   write(*, *) "fil_ip=", fil_ip
   write(*, *) "fil_in=", fil_in
   write(*, *) "bin_ip=", bin_ip
   write(*, *) "bin_in=", bin_in
   write(*, *) "dat_ip=", dat_ip
   write(*, *) "dat_in=", dat_in
   write(*, *) "avg_in=", avg_in
   write(*, *) "   chr=", rchr
   write(*, *) ""
   write(*, *) ""


   ! Check that IP and input BED files exist
   inquire(file = fil_ip, exist = file_exists)
   if (.not. file_exists) then
      write(6, *) "Error: IP BED file not found: ", fil_ip
      stop 1
   endif

   inquire(file = fil_in, exist = file_exists)
   if (.not. file_exists) then
      write(6, *) "Error: Input BED file not found: ", fil_in
      stop 1
   endif

   ! Assign file paths to the path array
   path(1) = fil_ip
   path(2) = fil_in

   ! Assign bin sizes to nleng
   nleng(1) = bin_ip
   nleng(2) = bin_in

   write(*, *) "path(1) = ", path(1)
   write(*, *) "path(2) = ", path(2)
   write(*, *) ""
   write(*, *) "nleng(1) = ", nleng(1)
   write(*, *) "nleng(2) = ", nleng(2)
   write(*, *) ""
   write(*, *) ""

   ! Check that output directories exist
   call check_exists_dir(dat_ip)
   call check_exists_dir(dat_in)
   call check_exists_dir(avg_in)

   ! Loop through the two input files
   do k = 1, 2
      open(12, file = path(k))                                                                    ! Open file for reading

      bin = 0d0                                                                                   ! Initialize bin array to zero
      leftend = 0                                                                                 ! Reset left boundary of current region
      rightend = 0                                                                                ! Reset right boundary of current region

      ! Open output files based on input file type
      if (k .eq. 1) open(88, file = dat_ip)                                                       ! File for treatment data
      if (k .eq. 2) open(88, file = dat_in)                                                       ! File for control data
      ! if (k .eq. 2) open(89, file = avg_in)                                                     ! File for averaged control data
   12345 continue  ! Label: Loop back after processing single input line    
      read(12, *, IOSTAT = Reason) unk, ilft, irght, leng                                         ! Read input line: chrom, start, end, length
      if (Reason .gt. 0) then
         write(*, *) 'There was an error in input file ', path(k)                                 ! Handle I/O errors
         stop                                                                                     ! Stop execution if there's an error
      elseif (Reason.eq.0) then
   23456 continue  ! Label: Begin to process successfully read input line
         ! Brad: "Only work the given chr by using this line"
         ! if(unk /= "chr1")go to 54321

         ! Process regions only if chromosome identifier is valid
         if (SCAN(unk(5:6), "_") .eq. 0) then                                                     ! Check that 5th and 6th characters of 'unk' do not contain underscore
            ! Handle new, non-overlapping regions
            if (ilft .gt. rightend .and. unk == rchr) then                                        ! If new region starts after current region
               ! Write out current bin data
               if (itops.gt.0) then                                                               ! Execute if current bin has data
                  lpos = leftend                                                                  ! Initialize position to left boundary
                  ave = 0d0                                                                       ! Reset total average accumulator
                  count = 0d0                                                                     ! Reset range counter

                  do i = 1, rightend - leftend, nleng(k)                                          ! Iterate over bins with step size of bin length
                     sum = 0                                                                      ! Reset sum for current range
                     iend = min(rightend - leftend, i + nleng(k) - 1)                             ! Define end of range
                     idub = 0                                                                     ! Reset bin counter for this range

                     do j = i, iend                                                               ! Loop through indices in range
                        lpos = lpos + 1                                                           ! Update genomic position
                        idub = idub + 1                                                           ! Increment number of bins in range
                        sum = sum + bin(j)                                                        ! Add bin value to sum
                     enddo

                     write(88, *) unk, lpos - idub, lpos - 1, dble(sum)                           ! Output range data: chrom, start, end, sum
                     ave = ave + dble(sum) / dble(idub)                                           ! Add average value for range to total average
                     count = count + 1d0                                                          ! Increment range counter
                  enddo

                  ! if (k .eq. 2) then
                  !    write(89, *) unk, leftend, rightend, ave / count                           ! Output average for control data
                  ! endif
               endif

               bin = 0                                                                            ! Reset bin array
               itops = 0                                                                          ! Reset maximum bin value
               
               ! Update boundaries for new region
               leftend = ilft
               rightend = ilft + leng

               do i = 1, rightend - leftend
                  bin(i) = bin(i) + 1d0 / dble(leng)                                              ! Initialize bins for new region with normalized values
                  itops = max(itops, bin(i))                                                      ! Update maximum bin value
               enddo

               go to 12345  ! Return to main processing loop
            elseif (ilft .le. rightend .and. unk == rchr) then                                    ! If region overlaps with current one
               ! Update bins if regions overlap
               iindicate = 0
               rightend = max(rightend, ilft + leng)                                              ! Extend right boundary to include new region
               ! Brad: "leftend = 1 in position so ilft - leftend + 1 is start location"
               ! if(1 + ilft - leftend .lt. 1)then
               !    write(*, *) unk, ilft, leftend, "ji", rchr
               !    stop
               ! endif
               if ((1 + ilft - leftend) .lt. 100000 .and. (rightend - leftend) .le. 100000) then  ! Ensure bin indices stay within array bounds (1 to 100000)
                  do i = min(100000, 1 + ilft - leftend), min(100000, 1 + ilft - leftend + leng)  ! Update bins safely within allowed range for current genomic region  ! rightend - leftend
                     bin(i) = bin(i) + 1d0 / dble(leng)                                           ! Update overlapping bins with normalized values
                     itops = max(itops, bin(i))                                                   ! Update maximum bin value
                  enddo   
                  go to 12345                                                                     ! Continue processing next region
               else  ! Brad: "The bounds of bin will be exceeded so reset bins now; can happen on inputs"
                  ! Handle bin overflow if region is too large  ! Brad: "Final update and dump it out"
                  do i = min(100000, 1 + ilft - leftend), min(100000, 1 + ilft - leftend + leng)  ! Handle bins that exceed array bounds
                     bin(i) = bin(i) + 1d0 / dble(leng)                                           ! Add normalized values to bins
                     itops = max(itops, bin(i))                                                   ! Update maximum bin value
                  enddo

                  ! Write out current bins
                  if (itops .gt. 1 * 0) then                                                      ! Check that there are data in current bins (not likely needed)  ! Can probably just change this to 0?
                     lpos = leftend                                                               ! Start genomic position counter at left boundary
                     ave = 0d0                                                                    ! Initialize total average accumulator for region
                     count = 0d0                                                                  ! Initialize counter for number of processed ranges

                     do i = 1, ilft, nleng(k)                                                     ! Loop through all bins in region in steps of bin size 'nleng(k)'
                        sum = 0                                                                   ! Reset cumulative sum for current range
                        idub = 0                                                                  ! Reset counter for number of bins in range
                        iend = min(ilft - leftend, i + nleng(k) - 1)                              ! Determine end index for range, ensuring 'iend' doesn't exceed region bounds ('ilft - leftend') or cause array overflow

                        do j = min(i, 100000), min(iend, 100000)                                  ! Loop through bins in current range
                           lpos = lpos + 1                                                        ! Increment genomic position for each bin processed
                           idub = idub + 1                                                        ! Increment bin count for range
                           sum = sum + bin(j)                                                     ! Add value of current bin to cumulative sum
                        enddo
                        ! write(88, *) unk, lpos - nleng(k) / 2, dble(sum) / dble(nleng(k))

                        if (idub .gt. 0) then                                                     ! Process range only if bins are in it
                           ! Output chromosome, range, and sum for bins in range
                           write(88, *) unk, lpos - idub, lpos - 1, dble(sum)
                           ave = ave + dble(sum) / dble(idub)                                     ! Add average of range to total average
                           count = count + 1d0                                                    ! Increment processed range counter
                        endif
                     enddo
                     iindicate = 1

                     ! if (k .eq. 2) then                                                         ! Special case: If processing control file (k = 2)
                     !    ! Output average signal for entire region to separate file
                     !    write(89, *) unk, leftend, rightend, ave / count                        ! Write average
                     ! endif
                  endif

                  ! Check that there are no updates to the bins
                  if (iindicate .eq. 0) then
                     bin = 0                                                                      ! Reset the bin array to zero
                     itops = 0                                                                    ! Reset the maximum bin value tracker
                     
                     ! Reset the genomic region boundaries to zero  ! Brad: "Set new bins"
                     leftend = 0
                     rightend = 0
                  else                                                                            ! Otherwise, process and shift bin data  ! Brad: "ilft should hold the last fragment read and"
                     ! Brad: "----------------- WORKING HERE TO SHIFT BIN correctly"
                     ! Variable 'ilft' holds the left boundary of the last fragment
                     
                     itops = 0                                                                    ! Clear the current max bin value  ! Brad: "Clear this"
                     j = 0                                                                        ! Initialize counter 'j' for shifting bins

                     ! Loop through the portion of the bin array that needs to be shifted
                     do i = 1 + ilft - leftend, 100000
                        j = j + 1                                                                 ! Increment the new bin index
                        bin(j) = bin(i)                                                           ! Shift values from the old bin array ('bin(i)') to the new array ('bin(j)')
                        itops = max(itops,bin(j))                                                 ! Update the maximum bin value tracker
                     enddo                                                                        ! End shifting of bins  ! Brad: "Anything outside should have been zero anyway"

                     ! Clear the remainder of the bin array, beyond the shifted portion
                     ! Brad: "We need to implement better bookkeeping right here, imho"
                     do i = j + 1, 100000
                        bin(i) = 0d0   ! Reset these bins to zero
                     enddo

                     ! Update the left and right boundaries for the genomic region
                     leftend = ilft
                     rightend = irght                                                             ! rightend  ! Brad: "same one"
                  endif                                                                           ! End of 'iindicate' check

               ! Update the chromosome being processed to the current region's chromosome
               rchr = unk

               ! Jump back to label 12345 to process the next input line
               go to 12345                                                                        ! Brad: "23456"
               endif                                                                              ! End of the block that protects bin boundaries  ! Brad: "Protect bounds of bin"
            endif                                                                                 ! End of condition for overlapping/non-overlapping regions

            ! If the chromosome identifier changes, finalize and reset for the new chromosome
            if (unk /= rchr) then
               ! Print debug information about the current state
               write(*, *) unk, rchr, k, itops

               ! If the current bins contain data, write it out and reset
               if (itops .gt. 1 * 0) then
                  lpos = leftend                                                                  ! Start genomic position counter at the left boundary
                  ave = 0d0                                                                       ! Initialize total average accumulator
                  count = 0d0                                                                     ! Initialize counter for processed ranges

                  ! Loop through all bins in steps of bin size
                  do i = 1, rightend - leftend, nleng(k)
                     sum = 0                                                                      ! Reset cumulative sum for the current range
                     idub = 0                                                                     ! Reset bin count for the current range
                     
                     ! Process each bin within the range
                     do j = i, min(rightend - leftend, i + nleng(k) - 1)                          ! Limit to valid bins  ! Brad: "Double test"
                        lpos = lpos + 1                                                           ! Increment genomic position
                        idub = idub + 1                                                           ! Increment bin count
                        sum = sum + bin(j)                                                        ! Add bin value to cumulative sum
                     enddo                                                                        ! End inner loop over bins

                     ! Write output for the range
                     write(88, *) rchr, lpos - idub, lpos - 1, dble(sum)  ! /dble(idub)
                     ! write(88, *) unk, lpos - nleng(k) / 2, dble(sum) / dble(nleng(k))

                     ! Update the total average and range counter
                     ave = ave + dble(sum) / dble(idub)
                     count = count + 1d0
                  enddo                                                                           ! End loop over all bins in the region

                  ! ! If processing control (k = 2), write average for the region
                  ! if (k .eq. 2) then
                  !    write(89, *) rchr, leftend, rightend, ave / count
                  ! endif
               endif

               ! Reset the bins and associated variables
               bin = 0
               itops = 0
               ! Brad: "Set new bins"
               leftend = 0
               rightend = 0

               ! Update the current chromosome identifier
               rchr = unk

               ! Jump to label 23456 to process the next region
               go to 23456
            endif                                                                                 ! End of chromosome switch check
         else                                                                                     ! If the region doesn't match the current chromosome
            ! Jump back to label 12345 to process the next input line
            go to 12345
         endif                                                                                    ! End of valid line read  ! Brad: "The _ screen"
      endif
      54321 continue                                                                              ! Label: Finalize processing of the current chromosome  ! Brad: "Jumpout landing"
      close(12)                                                                                   ! Close the current input file
      close(88)                                                                                   ! Close the main output file
   enddo                                                                                          ! End of the outer loop over the two input files  ! Brad: "k files loop"

contains

   subroutine parse_args(fil_ip, fil_in, bin_ip, bin_in, dat_ip, dat_in, avg_in, rchr)
      ! Subroutine to parse command-line arguments

      ! Arguments
      character(len = 256), intent(out) :: fil_ip, fil_in, dat_ip, dat_in, avg_in
      integer, intent(out) :: bin_ip, bin_in
      character(len = 6), intent(out) :: rchr

      ! Local variables
      character(len = 256) :: arg, tmp_arg
      integer :: i
      logical :: has_bin_ip, has_bin_in

      ! Initialize defaults
      fil_ip = ""
      fil_in = ""
      bin_ip = 30
      bin_in = 30
      dat_ip = "IP.data"
      dat_in = "in.data"
      avg_in = "in_avg.data"
      rchr = "chr1"
      has_bin_ip = .false.
      has_bin_in = .false.

      ! Check if no arguments are passed or help is requested
      if (iargc() == 0) then
         call print_help()
         stop
      endif

      ! Loop over all command-line arguments
      do i = 1, iargc()
         call getarg(i, arg)

         ! Match known arguments
         if (arg == "-h" .or. arg == "--help") then
            call print_help()
            stop
         elseif (arg(1:4) == "-fp=" .or. arg(1:9) == "--fil_ip=") then
            fil_ip = trim(adjustl(arg(index(arg, "=") + 1:)))
         elseif (arg(1:4) == "-fn=" .or. arg(1:9) == "--fil_in=") then
            fil_in = trim(adjustl(arg(index(arg, "=") + 1:)))
         elseif (arg(1:4) == "-bp=" .or. arg(1:9) == "--bin_ip=") then
            tmp_arg = trim(adjustl(arg(index(arg, "=") + 1:)))
            read(tmp_arg, *) bin_ip
            has_bin_ip = .true.
         elseif (arg(1:4) == "-bn=" .or. arg(1:9) == "--bin_in=") then
            tmp_arg = trim(adjustl(arg(index(arg, "=") + 1:)))
            read(tmp_arg, *) bin_in
            has_bin_in = .true.
         elseif (arg(1:4) == "-dp=" .or. arg(1:9) == "--dat_ip=") then
            dat_ip = trim(adjustl(arg(index(arg, "=") + 1:)))
         elseif (arg(1:4) == "-dn=" .or. arg(1:9) == "--dat_in=") then
            dat_in = trim(adjustl(arg(index(arg, "=") + 1:)))
         elseif (arg(1:4) == "-an=" .or. arg(1:9) == "--avg_in=") then
            avg_in = trim(adjustl(arg(index(arg, "=") + 1:)))
         elseif (arg(1:3) == "-c=" .or. arg(1:6) == "--chr=") then
            rchr = trim(adjustl(arg(index(arg, "=") + 1:)))
         else
            write(6, *) "Error: Unrecognized argument: ", arg
            stop 1
         endif
      enddo

      ! Check required arguments
      if (fil_ip == "") then
         write(6, *) "Error: Missing required argument -fp/--fil_ip."
         stop 1
      endif

      if (fil_in == "") then
         write(6, *) "Error: Missing required argument -fn/--fil_in."
         stop 1
      endif

      if (.not. has_bin_ip) then
         write(6, *) "Error: Missing required argument -bp/--bin_ip."
         stop 1
      endif

      if (bin_ip < 1) then
         write(6, *) "Error: -bp/--bin_ip must be a positive integer."
         stop 1
      endif

      if (.not. has_bin_in) then
         write(6, *) "Error: Missing required argument -bn/--bin_in."
         stop 1
      endif

      if (bin_in < 1) then
         write(6, *) "Error: -bn/--bin_in must be a positive integer."
         stop 1
      endif

      if (dat_ip == "") then
         write(6, *) "Error: Missing required argument -dp/--dat_ip."
         stop 1
      endif

      if (dat_in == "") then
         write(6, *) "Error: Missing required argument -dn/--dat_in."
         stop 1
      endif

      if (avg_in == "") then
         write(6, *) "Error: Missing required argument -an/--avg_in."
         stop 1
      endif

      if (rchr == "") then
         write(6, *) "Error: Missing required argument -c/--chr."
         stop 1
      endif
   end subroutine parse_args

   subroutine print_help()
      write(*, *) "Usage:"
      write(*, *) "  tracks"
      write(*, *) "    --fil_ip=<str> --fil_in=<str> [--bin_ip=<int>]"
      write(*, *) "    [--bin_in=<int>] [--dat_ip=<str>] [--dat_in=<str>]"
      write(*, *) "    [--avg_in=<int>] [--chr=<str>]"
      write(*, *) ""
      write(*, *) "Required arguments:"
      write(*, *) "  -fp, --fil_ip  Path to the IP BED file."
      write(*, *) "  -fn, --fil_in  Path to the input BED file."
      write(*, *) ""
      write(*, *) "Optional arguments:"
      write(*, *) "   -h, --help    Show this help message and exit."
      write(*, *) "  -bp, --bin_ip  Bin size for IP file processing (default: 30)."
      write(*, *) "  -bn, --bin_in  Bin size for input file processing (default: 30)."
      write(*, *) "  -dp, --dat_ip  Outfile for IP data (default: IP.data)."
      write(*, *) "  -dn, --dat_in  Outfile for input data (default: in.data)."
      write(*, *) "  -an, --avg_in  Outfile for averaged input data (default:"
      write(*, *) "                 in_avg.data)."
      write(*, *) "   -c, --chr     First chromosome in model organism"
      write(*, *) "                 (default: chr1)."
      write(*, *) ""
   end subroutine print_help

   subroutine check_exists_dir(file_path)
      character(len = *), intent(in) :: file_path
      character(len = 256) :: dir_path
      logical :: dir_exists
      integer :: last_slash

      ! Find the last slash to extract the directory
      last_slash = index(file_path, "/", back = .true.)
      if (last_slash > 0) then
         dir_path = file_path(:last_slash)
      else
         dir_path = "./"  ! Default to current directory if no slash
      endif

      ! Check if the directory exists
      inquire(file = dir_path, exist = dir_exists)
      if (.not. dir_exists) then
         write(6, *) "Error: Directory does not exist: ", trim(dir_path)
         stop 1
      endif
   end subroutine check_exists_dir

end program