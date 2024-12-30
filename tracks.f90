! Declare variables
real*8 bin(100000), sum, itops  ! Array for binning data, 'bin'; scalars for cumulative sum, 'sum', and max bin value, 'itops'
! Declare integer variables:
! `nleng(2)`: Array for bin sizes of the two input files.
! `ilft`, `irght`: Start and end positions of the genomic region being processed.
! `leng`: Length of the region (`irght - ilft`).
! `leftend`, `rightend`: Boundaries of the current binning range.
! `Reason`: I/O status code after reading input.
! `lpos`: Tracks current genomic position in the `bin` array.
integer nleng(2), ilft, irght, leng, leftend, rightend, Reason, lpos
character(len=6) :: unk, rchr      ! Chromosome identifiers
character(len=5) :: unk2, rchr2    ! Alternate chromosome identifiers (unused)
character(len=44) :: siq           ! Placeholder for HMM annotations (unused)
character(len=62) :: arg, path(2)  ! Infile paths
logical :: file_exists             ! File existence flag
inot = 0       ! Integer variable initialized to 0 (unused)
itops = 0      ! Tracks maximum value in bin array
iread = 0      ! Integer variable initialized to 0 (unused)
rchr = "chr1"  ! Start processing from chromosome "chr1"

! Read command-line arguments to get file paths and bin sizes
do i = 1, iargc()              ! Loop through all CL arguments
   if (i.lt.3) then            ! First two arguments are infile paths
      call getarg(i, arg)      ! Get the argument string
      path(i) = arg            ! Store file path
   else                        ! Third argument onwards specifies bin sizes
      call getarg(i, arg)      ! Get argument string
      read(arg, *) nleng(i-2)  ! Read and store bin size for each input file
   endif
enddo

! Loop through the two input files
do k = 1, 2
   inquire(file=path(k), EXIST=file_exists)  ! Check that file exists
   if (file_exists) then
      open(12, file=path(k))                 ! Open file for reading
   else
      write(*, *) 'Your first file or path is incorrect'
      stop                                   ! Terminate program if file is missing
   endif

   bin = 0d0     ! Initialize bin array to zero
   leftend = 0   ! Reset left boundary of current region
   rightend = 0  ! Reset right boundary of current region

   ! Open output files based on input file type
   if (k.eq.1) open(88, file='IP.data')    ! File for treatment data
   if (k.eq.2) open(88, file='IN.data')    ! File for control data
   if (k.eq.2) open(89, file='INave.data') ! File for averaged control data

12345 continue  ! Label: Loop back after processing single input line
   read(12, *, IOSTAT=Reason) unk, ilft, irght, leng  ! Read input line: chrom, start, end, length
   if (Reason.gt.0) then
      write(*, *) 'There was an error in input file ', path(k)  ! Handle I/O errors
      stop                                                      ! Stop execution if there's an error
   elseif (Reason.eq.0) then
23456 continue  ! Label: Begin to process successfully read input line
      ! Skip regions not on current chromosome
      if (unk /= rchr) go to 54321  ! If chromosome changes, process current chromosome and switch

      ! Process regions only if chromosome identifier is valid
      if (SCAN(unk(5:6), "_").eq.0) then  ! Check that 5th and 6th characters of 'unk' don't contain underscore
         ! Handle new, non-overlapping regions
         if (ilft.gt.rightend.and.unk == rchr) then           ! If new region starts after current region
            ! Write out current bin data
            if (itops.gt.0) then                              ! Execute if current bin has data
               lpos = leftend                                 ! Initialize position counter to left boundary
               ave = 0d0                                      ! Reset total average accumulator
               count = 0d0                                    ! Reset range counter
               do i = 1, rightend-leftend, nleng(k)           ! Iterate over bins with step size of bin length
                  sum = 0                                     ! Reset sum for current range
                  iend = min(rightend-leftend, i+nleng(k)-1)  ! Define end of the range
                  idub = 0                                    ! Reset bin counter for this range
                  do j = i, iend                              ! Loop through indices in range
                     lpos = lpos + 1                          ! Update genomic position
                     idub = idub + 1                          ! Increment number of bins in range
                     sum = sum + bin(j)                       ! Add bin value to sum
                  enddo
                  write(88, *) unk, lpos-idub, lpos-1, dble(sum)  ! Output range data: chromosome, start, end, sum
                  ave = ave + dble(sum) / dble(idub)  ! Add average value for range to total average
                  count = count + 1d0                 ! Increment range counter
               enddo
               if (k.eq.2) then
                  write(89, *) unk, leftend, rightend, ave / count  ! Output average for control data
               endif
            endif
            bin = 0    ! Reset bin array
            itops = 0  ! Reset maximum bin value

            ! Update boundaries for new region
            leftend = ilft
            rightend = ilft + leng
            do i = 1, rightend-leftend
               bin(i) = bin(i) + 1d0 / dble(leng)  ! Initialize bins for new region with normalized values
               itops = max(itops, bin(i))          ! Update maximum bin value
            enddo
            go to 12345  ! Return to main processing loop
         elseif (ilft.le.rightend.and.unk == rchr) then  ! If region overlaps with current one
            ! Update bins if regions overlap
            rightend = max(rightend, ilft + leng)        ! Extend right boundary to include new region
            if (1+ilft-leftend.lt.100000.and.rightend-leftend.le.100000) then        ! Ensure bin indices stay within array bounds (1 to 100000)
               do i = min(100000, 1+ilft-leftend), min(100000, 1+ilft-leftend+leng)  ! Update bins safely within allowed range for current genomic region
                  bin(i) = bin(i) + 1d0 / dble(leng)  ! Update overlapping bins with normalized values
                  itops = max(itops, bin(i))          ! Update maximum bin value
               enddo
               go to 12345                            ! Continue processing next region
            else
               ! Handle bin overflow if region is too large
               do i = min(100000, 1+ilft-leftend), min(100000, 1+ilft-leftend+leng)  ! Handle bins that exceed array bounds
                  bin(i) = bin(i) + 1d0 / dble(leng)  ! Add normalized values to bins
                  itops = max(itops, bin(i))          ! Update maximum bin value
               enddo

               ! Write out current bins
               if (itops.gt.0) then  ! Check that there are data in current bins
                  lpos = leftend     ! Start genomic position counter at left boundary
                  ave = 0d0          ! Initialize total average accumulator for region
                  count = 0d0        ! Initialize counter for number of processed ranges
                  do i = 1, ilft, nleng(k)  ! Loop through all bins in region in steps of bin size 'nleng(k)'
                     sum = 0         ! Reset cumulative sum for current range
                     idub = 0        ! Reset counter for number of bins in range
                     iend = min(ilft-leftend, i+nleng(k)-1)  ! Determine end index for range, ensuring 'iend' doesn't exceed region bounds ('ilft-leftend') or cause array overflow
                     
                     do j = min(i, 100000), min(iend, 100000)  ! Loop through bins in current range
                        lpos = lpos + 1                        ! Increment genomic position for each bin processed
                        idub = idub + 1                        ! Increment bin count for range
                        sum = sum + bin(j)                     ! Add value of current bin to cumulative sum
                     enddo

                     if (idub.gt.0) then  ! Process the range only if bins are in it
                        ! Output chromosome, range, and sum for bins in range
                        write(88, *) unk, lpos-idub, lpos-1, dble(sum)
                        ave = ave + dble(sum) / dble(idub)  ! Add range's average to total average
                        count = count + 1d0                 ! Increment processed range counter
                     endif
                  enddo

                  if (k.eq.2) then  ! Special case: If processing control file (k = 2)
                     ! Output average signal for entire region to separate file
                     write(89, *) unk, leftend, rightend, ave / count  ! Write average
                  endif
               endif

               bin = 0           ! Reset bin array to prepare for next region
               itops = 0         ! Reset maximum bin value for next region
               leftend = ilft    ! Update left boundary for next genomic region
               rightend = irght  ! Update right boundary for next genomic region
            endif
            go to 12345  ! Jump back to process next line of input
         endif
      endif
54321 continue   ! Label: Mark end of processing for current chromosome
      close(12)  ! Close infile for current dataset
      close(88)  ! Close main outfile
   enddo         ! End loop over the two input files
end program