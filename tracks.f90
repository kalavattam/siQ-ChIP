program tracks

! Rough documentation:
! This script takes five positional parameters; the first four are required:
! 1. IP BED file
! 2. input BED file
! 3. IP bin width
! 4. input bin width
! 5. (optional) first chromosome of model organism (default: "I")

! Declare variables
real*8 bin(100000), sum, itops
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
character(len = 62) :: arg, path(2)                                                            ! Infile paths
logical :: file_exists                                                                         ! File existence flag
inot = 0                                                                                       ! Integer variable initialized to 0 (unused)
itops = 0                                                                                      ! Tracks maximum value in bin array
iread = 0                                                                                      ! Integer variable initialized to 0 (unused)
rchr = "I"                                                                                     ! Default behavior: Start processing from chromosome "I"

! Parse the command line arguments to get file paths, bin sizes, and optional rchr
do i = 1, iargc()                                                                              ! Loop through all CL arguments
   if (i .lt. 3) then                                                                          ! First two arguments are infile paths
      call getarg(i, arg)                                                                      ! Get the argument string
      path(i) = arg                                                                            ! Store file path
   elseif (i .le. 4) then                                                                      ! Next two arguments are bin sizes
      call getarg(i, arg)                                                                      ! Get argument string
      read(arg, *) nleng(i - 2)                                                                ! Read and store bin size for each input file
   elseif (i .eq. 5) then                                                                      ! Fifth argument is optional rchr
      call getarg(i, rchr)                                                                     ! Optionally, assign rchr from command-line argument (i.e., if not "I")
   endif
enddo

! Loop through the two input files
do k = 1, 2
   inquire(file = path(k), EXIST = file_exists)                                                ! Check that file exists
   if (file_exists) then
      open(12, file = path(k))                                                                 ! Open file for reading
   else
      write(*, *) 'Your first file or path is incorrect'
      stop                                                                                     ! Terminate program if file is missing
   endif

   bin = 0d0                                                                                   ! Initialize bin array to zero
   leftend = 0                                                                                 ! Reset left boundary of current region
   rightend = 0                                                                                ! Reset right boundary of current region

   ! Open output files based on input file type
   if (k.eq.1) open(88, file = 'IP.data')                                                      ! File for treatment data
   if (k.eq.2) open(88, file = 'IN.data')                                                      ! File for control data
   if (k.eq.2) open(89, file = 'INave.data')                                                   ! File for averaged control data
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

               if (k .eq. 2) then
                  write(89, *) unk, leftend, rightend, ave / count                             ! Output average for control data
               endif
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
            !    write(*,*) unk, ilft, leftend, "ji", rchr
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
               if (itops.gt.1*0) then                                                          ! Check that there are data in current bins (not likely needed)
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

                  if (k .eq. 2) then                                                           ! Special case: If processing control file (k = 2)
                     ! Output average signal for entire region to separate file
                     write(89, *) unk, leftend, rightend, ave / count                          ! Write average
                  endif
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

               ! If processing control (k = 2), write average for the region
               if (k .eq. 2) then
                  write(89, *) rchr, leftend, rightend, ave / count
               endif
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

end program