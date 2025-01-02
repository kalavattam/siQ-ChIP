
! Rough documentation:
! This script takes five positional parameters; the first four are required:
! 1. IP BED file
! 2. input BED file
! 3. IP bin width
! 4. input bin width
! 5. (optional) first chromosome of model organism (default: "I")

      integer bin(100000)
      real * 8 iphits, inhits, totalcount, count, dep, factr
      integer Reason, Reason2, ilft, irght, ilft2, irght2
      character(len = 5) :: unk,rchr !should have been filtered to chrxx by now
      character(len = 5) :: unk2,rchr2
      character(len = 44) :: siq  ! For HMM annotations
      character(len = 62) :: arg
      character(len = 62) :: path(3)
      logical :: file_exists
      totalcount = 0d0
      count = 0d0
      inot = 0
      iread = 0
      rchr = "I"  ! Default value for rchr is "I"

      ! Handle command-line arguments
      do i = 1, 2  ! First two arguments: Paths to input files
         call getarg(i, arg)
         path(i) = arg
      enddo

      inquire(file = path(1), EXIST = file_exists)
      if (file_exists .eqv. .true.) then
         open(12, file = path(1))      
      else
         write(*, *) 'Your first file or path is incorrect.'
         stop
      endif

      inquire(file = path(2), EXIST = file_exists)
      if (file_exists .eqv. .true.) then
         open(13, file = path(2))      
      else
         write(*, *) 'Your second file or path is incorrect.'
         stop
      endif

      ! Third argument: factr (alpha)
      call getarg(3, arg)
      read(arg, *) factr  ! This is alpha

      ! Fourth argument: dep (expected depth)
      call getarg(4, arg)
      read(arg, *) dep    ! This is expected depth input

      ! Optional fifth argument: rchr
      if (iargc() .eq. 5) then
         call getarg(5, rchr)  ! Optionally, assign value of rchr from command line (i.e., if not "I")
      endif

      open(88, file = 'mergedSIQ.data')

      ! Match intervals as much as possible
12121 continue
      read(12, *, IOSTAT = Reason) unk, ilft, irght, iphits
      if (Reason .gt. 0) then
         write(*,*) 'There was an error in input file ', path(1)
         stop
      elseif (Reason .eq. 0) then 
12344    continue
         read(13, *, IOSTAT = Reason2) unk2, ilft2, irght2, inhits
         if (Reason2 .gt. 0) then
            write(*,*) 'There was an error in input file ', path(2)
            stop
         elseif (Reason2 .eq. 0) then
12345       continue
            if (unk == unk2 .and. ilft .le. irght2 .and. ilft2 .le. irght) then  ! Intersect
!               if (inhits.lt.1)inhits=oin!no boom
!               if (inhits.gt.0d0) then!dont smooth input
!               dep=5.29491465d0! a485
!               dep=5.915557010625d0! cbp
!               dep=5.507106525d0! dmso
!               dep=4.725d0 !k9ac low
!               dep=inhits!use input
               totalcount = totalcount + 1d0
               if (inhits .gt. dep) then  ! 7d0 is 100m-depth estimate
               write(88, *) unk, ilft, irght, factr*iphits/inhits
!!               write(88, *) unk, ilft, irght, factr*iphits/1d0
               else  ! These lines are for smoothing input
!               write(88, *) unk, ilft, irght, factr * (iphits/7d0)
!               write(88, *) unk, ilft, irght, factr * (iphits/inhits)  ! This is to block fakeinputs
               write(88, *) unk, ilft, irght, factr*iphits/dep  ! This has fake input
               count = count + 1d0
               endif
               ! Get a new IP line and return to check
               read(12, *, IOSTAT = Reason) unk, ilft, irght, iphits
               if (Reason .eq. 0)go to 12345  ! Check start
            elseif (unk2 == unk .and. ilft .gt. irght2) then
               oin=inhits!save
               go to 12344!hasnt matched but need a new comp line
            elseif (unk2 == unk .and. ilft2 .gt. irght) then  ! Get new line
               read(12,*,IOSTAT=Reason) unk, ilft, irght, iphits
               if (Reason .eq. 0) go to 12345  ! Check start
            elseif (unk2 /= unk) then
               if (ilft .gt. irght2) then  ! Because 13- went to new chr
                  go to 12121
               endif
               oin = inhits  ! Save
               go to 12344   ! Get new comp
            endif
         endif  ! Close 13 read
      endif  ! Close 12 read
!      write(*,*) "i stop now at"
!      write(*,*) unk, ilft, irght, iphits
!      write(*,*) unk2, ilft2, irght2, inhits      
      close(88)
      close(12)
      close(13)
      write(*, *)"i counted it as this ", totalcount, count
      end program
