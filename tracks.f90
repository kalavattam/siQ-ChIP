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
inot = 0    ! Integer variable initialized to 0 (unused)
itops = 0   ! Tracks maximum value in bin array
iread = 0   ! Integer variable initialized to 0 (unused)
rchr = "I"  ! Start processing from chromosome "I"
      
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

do k=1,2 !files loop
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
   read(12,*,IOSTAT=Reason) unk, ilft, irght, leng  ! Read input line: chrom, start, end, length
   if (Reason.gt.0) then
      write(*, *) 'There was an error in input file ', path(k)  ! Handle I/O errors
      stop                                                      ! Stop execution if there's an error
   elseif (Reason.eq.0) then
23456 continue  ! Label: Begin to process successfully read input line
      ! only work the given chr by using this line
      ! if(unk /= "chr1")go to 54321

      if(SCAN(unk(5:6), "_").eq.0)then

         if(ilft.gt.rightend.and.unk == rchr)then
            !write out the current bin of data
            if(itops.gt.1*0)then !nonempty
               lpos=leftend
               ave=0d0
               count=0d0
               do i=1,rightend-leftend,nleng(k)
                  sum=0
                  iend=min(rightend-leftend,i+nleng(k)-1)!double test
                  idub=0
                  do j=i,iend
                     lpos=lpos+1
                     idub=idub+1
                     sum=sum+bin(j)
                  enddo
      !                  write(88,*) unk, lpos-nleng(k)/2, dble(sum)/dble(nleng(k))
                  write(88,*) unk, lpos-idub,lpos-1, dble(sum)!/dble(idub)
                  ave=ave+dble(sum)/dble(idub)
                  count=count+1d0
               enddo
               if(k.eq.2)then
                  write(89,*) unk, leftend, rightend, ave/count
               endif
            endif
            bin=0
            itops=0
            !set new bins
            leftend=ilft
            rightend=ilft+leng!same as irght sometimes
            do i=1,rightend-leftend
               bin(i)=bin(i)+1d0/dble(leng)
               itops=max(itops,bin(i))
            enddo
            go to 12345
         elseif(ilft.le.rightend.and.unk == rchr)then
            iindicate=0
            rightend=max(rightend,ilft+leng)!new right end
      !            leftend = 1 in position so ilft-leftend+1 is start location
      !            if(1+ilft-leftend.lt.1)then
      !               write(*,*) unk, ilft, leftend, "ji", rchr
      !               stop
      !            endif
            if(1+ilft-leftend.lt.100000.and.rightend-leftend.le.100000)then
            do i=min(100000,1+ilft-leftend),min(100000,1+ilft-leftend+leng)!rightend-leftend
               bin(i)=bin(i)+1d0/dble(leng)
               itops=max(itops,bin(i))
            enddo   
            go to 12345
            else!the bounds of bin will be exceeded so reset bins now-can happen on inputs
               !final update and dump it out
               do i=min(100000,1+ilft-leftend),min(100000,1+ilft-leftend+leng)!rightend-leftend
                  bin(i)=bin(i)+1d0/dble(leng)
                  itops=max(itops,bin(i))
               enddo

            if(itops.gt.1*0)then!not likely needed
               lpos=leftend
               ave=0d0
               count=0d0
               do i=1,ilft,nleng(k)!rightend-leftend,nleng(k)
                  sum=0
                  idub=0
                  iend=min(ilft-leftend,i+nleng(k)-1)!double test
                  do j=min(i,100000),min(iend,100000)!double test
                     lpos=lpos+1
                     idub=idub+1
                     sum=sum+bin(j)
                  enddo
      !                  write(88,*) unk, lpos-nleng(k)/2, dble(sum)/dble(nleng(k))
                  if(idub.gt.1)then
                  write(88,*) unk, lpos-idub,lpos-1, dble(sum)!/dble(idub)
                  ave=ave+dble(sum)/dble(idub)
                  count=count+1d0
                  endif
               enddo
               iindicate=1
               if(k.eq.2)then
                  write(89,*) unk, leftend, rightend, ave/count
               endif

            endif
            if(iindicate.eq.0)then
               bin=0
            itops=0
            !set new bins
            leftend=0
            rightend=0
            else!ilft should hold the last fragment read and 
      ! ----------------- WORKING HERE TO SHIFT BIN correctly
               itops=0!clear this
               j=0
               do i=1+ilft-leftend,100000
                  j=j+1
                  bin(j)=bin(i)
                  itops=max(itops,bin(j))                  
               enddo!anything outside should have been zero anyway
               !we need to implement better bookkeeping right here, imho
               do i=j+1,100000
                  bin(i)=0d0
               enddo
               leftend=ilft
               rightend=irght!rightend!same one
            endif
            rchr=unk
            go to 12345!23456
            endif!protect bounds of bin
         endif
         if(unk /= rchr)then
            write(*,*) unk, rchr, k, itops
            if(itops.gt.1*0)then
               lpos=leftend
               ave=0d0
               count=0d0
               do i=1,rightend-leftend,nleng(k)
                  sum=0
                  idub=0
                  do j=i,min(rightend-leftend,i+nleng(k)-1)!double test
                     lpos=lpos+1
                     idub=idub+1
                     sum=sum+bin(j)
                  enddo
                  write(88,*) rchr, lpos-idub,lpos-1, dble(sum)!/dble(idub)
      !                  write(88,*) unk, lpos-nleng(k)/2, dble(sum)/dble(nleng(k))
                  ave=ave+dble(sum)/dble(idub)
                  count=count+1d0
               enddo
               if(k.eq.2)then
                  write(89,*) rchr, leftend, rightend, ave/count
               endif
            endif
            bin=0
            itops=0
            !set new bins
            leftend=0
            rightend=0
            rchr=unk
            go to 23456
         endif
      else
         go to 12345
      endif!the _ screen
   endif
   54321 continue!jumpout landing
   close(12)
   close(88)
enddo!k files loop

    end program