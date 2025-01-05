program tracks

   ! Rough documentation:
   !   'mergetracks' processes IP and input BED files to generate 'merged',
   !   siQ-scaled data.

   integer bin(100000)
   real * 8 iphits, inhits, totalcount, count, dep, factr
   integer Reason, Reason2, ilft, irght, ilft2, irght2
   character(len = 5) :: unk, rchr ! Should have been filtered to chrxx by now
   character(len = 5) :: unk2, rchr2
   character(len = 44) :: siq  ! For HMM annotations
   character(len = 62) :: arg
   character(len = 256) :: fil_ip, fil_in, fil_siq
   logical :: file_exists
   totalcount = 0d0
   count = 0d0
   inot = 0
   iread = 0
   ! rchr = "I"  ! Default value for rchr is "I"

   ! Call argument parser
   call parse_args(fil_ip, fil_in, fil_siq, factr, dep, rchr)

   ! Debug or verify arguments (optional)
   write(*, *) "######################"
   write(*, *) "## Parsed arguments ##"
   write(*, *) "######################"
   write(*, *) ""
   write(*, *) " fil_ip=", fil_ip
   write(*, *) " fil_in=", fil_in
   write(*, *) "fil_siq=", fil_siq
   write(*, *) "  factr=", factr
   write(*, *) "    dep=", dep
   write(*, *) "    chr=", rchr

   inquire(file = fil_ip, EXIST = file_exists)
   if(file_exists .eqv. .true.)then
      open(12, file = fil_ip)      
   else
      write(*, *) 'Your first file or path is incorrect.'
      stop
   endif
   
   inquire(file = fil_in, EXIST = file_exists)
   if(file_exists .eqv. .true.)then
      open(13, file = fil_in)      
   else
      write(*,*) 'Your second file or path is incorrect.'
      stop
   endif

   call check_exists_dir(fil_siq)
   open(88, file = fil_siq)

   ! Match intervals as much as possible
12121 continue
   read(12, *, IOSTAT = Reason) unk, ilft, irght, iphits
   if (Reason .gt. 0) then
      write(*, *) 'There was an error in input file ', file_ip
      stop
   elseif (Reason .eq. 0) then 
12344    continue
      read(13, *, IOSTAT = Reason2) unk2, ilft2, irght2, inhits
      if (Reason2 .gt. 0) then
         write(*, *) 'There was an error in input file ', fil_in
         stop
      elseif (Reason2 .eq. 0) then
12345       continue
         if (unk == unk2 .and. ilft .le. irght2 .and. ilft2 .le. irght) then  ! Intersect
            ! if (inhits .lt. 1) inhits = oin  !no boom
            ! if (inhits .gt. 0d0) then  !dont smooth input
            ! dep = 5.29491465d0  ! a485
            ! dep = 5.915557010625d0  ! cbp
            ! dep = 5.507106525d0  ! dmso
            ! dep = 4.725d0  ! k9ac low
            ! dep = inhits  !use input
            totalcount = totalcount + 1d0

            if (inhits .gt. dep) then  ! 7d0 is 100m-depth estimate
               write(88, *) unk, ilft, irght, factr * iphits / inhits
               ! write(88, *) unk, ilft, irght, factr * iphits / 1d0
            else  ! These lines are for smoothing input
               ! write(88, *) unk, ilft, irght, factr * (iphits / 7d0)
               ! write(88, *) unk, ilft, irght, factr * (iphits / inhits)  ! This is to block fakeinputs
               write(88, *) unk, ilft, irght, factr * iphits / dep  ! This has fake input
               count = count + 1d0
            endif

            ! Get a new IP line and return to check
            read(12, *, IOSTAT = Reason) unk, ilft, irght, iphits
            if (Reason .eq. 0) go to 12345  ! Check start
         elseif (unk2 == unk .and. ilft .gt. irght2) then
            oin = inhits!save
            go to 12344!hasnt matched but need a new comp line
         elseif (unk2 == unk .and. ilft2 .gt. irght) then  ! Get new line
            read(12, *, IOSTAT = Reason) unk, ilft, irght, iphits
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
   ! write(*, *) "i stop now at"
   ! write(*, *) unk, ilft, irght, iphits
   ! write(*, *) unk2, ilft2, irght2, inhits      
   close(88)
   close(12)
   close(13)
   write(*, *) "i counted it as this ", totalcount, count
   ! write(*, '(A)') "totalcount" // CHAR(9) // "count"
   ! write(*, '(F10.3, A, F10.3)') totalcount, CHAR(9), count

contains

   subroutine parse_args(fil_ip, fil_in, fil_siq, factr, dep, rchr)
      ! Subroutine to parse command-line arguments
   
      ! Arguments
      character(len = 256), intent(out) :: fil_ip, fil_in, fil_siq
      real(8), intent(out) :: factr, dep
      character(len = 5), intent(out) :: rchr
   
      ! Local variables
      character(len = 256) :: arg, tmp_arg
      integer :: i
      logical :: has_factr, has_dep
   
      ! Initialize defaults
      fil_ip = ""
      fil_in = ""
      fil_siq = ""
      factr = 1.0d0
      dep = 1.0d0
      rchr = "chr1"
      has_factr = .false.
      has_dep = .false.
   
      ! Check if no arguments are passed or help is requested
      if (iargc() == 0) then
         call print_help()
         stop
      endif
   
      ! Loop over all command-line arguments
      do i = 1, iargc()
         call getarg(i, arg)
   
         ! Match known arguments
         ! Match known arguments
         if (arg == "-h" .or. arg == "--help") then
            call print_help()
            stop
         elseif (arg(1:4) == "-fp=" .or. arg(1:9) == "--fil_ip=") then
            fil_ip = trim(adjustl(arg(index(arg, "=") + 1:)))
         elseif (arg(1:4) == "-fn=" .or. arg(1:9) == "--fil_in=") then
            fil_in = trim(adjustl(arg(index(arg, "=") + 1:)))
         elseif (arg(1:4) == "-fs=" .or. arg(1:10) == "--fil_siq=") then
            fil_siq = trim(adjustl(arg(index(arg, "=") + 1:)))
         elseif (arg(1:4) == "-fa=" .or. arg(1:9) == "--factr=") then
            tmp_arg = trim(adjustl(arg(index(arg, "=") + 1:)))
            read(tmp_arg, *) factr
            has_factr = .true.
         elseif (arg(1:3) == "-d=" .or. arg(1:6) == "--dep=") then
            tmp_arg = trim(adjustl(arg(index(arg, "=") + 1:)))
            read(tmp_arg, *) dep
            has_dep = .true.
         elseif (arg(1:3) == "-c=" .or. arg(1:6) == "--chr=") then
            rchr = trim(adjustl(arg(index(arg, "=") + 1:)))
         else
            write(*, *) "Error: Unrecognized argument: ", arg
            stop
         endif
      enddo
   
      ! Check required arguments
      if (fil_ip == "") then
         write(*, *) "Error: Missing required argument -fp/--fil_ip."
         stop
      endif
   
      if (fil_in == "") then
         write(*, *) "Error: Missing required argument -fn/--fil_in."
         stop
      endif
   
      if (fil_siq == "") then
         write(*, *) "Error: Missing required argument -fs/--fil_siq."
         stop
      endif
   
      if (.not. has_factr) then
         write(*, *) "Error: Missing required argument -fa/--factr."
         stop
      endif
   
      if (factr <= 0.0d0) then
         write(*, *) "Error: -fa/--factr must be a positive value greater than zero."
         stop
      endif
   
      if (.not. has_dep) then
         write(*, *) "Error: Missing required argument -d/--dep."
         stop
      endif
   
      if (dep <= 0.0d0) then
         write(*, *) "Error: -d/--dep must be a positive value greater than zero."
         stop
      endif
   
      if (rchr == "") then
         write(*, *) "Error: Missing required argument -c/--chr."
         stop
      endif
   end subroutine parse_args
   
   subroutine print_help()
      write(*, *) "Usage:"
      write(*, *) "  merge_tracks"
      write(*, *) "    --fil_ip=<str> --fil_in=<str> --fil_siq=<str>"
      write(*, *) "    --factr=<flt> --dep=<flt> [--chr=<str>]"
      write(*, *) ""
      write(*, *) "Required arguments:"
      write(*, *) "  -fp, --fil_ip   Path to the IP BED file."
      write(*, *) "  -fn, --fil_in   Path to the input BED file."
      write(*, *) "  -fs, --fil_siq  Outfile for merged, siQ-scaled data."
      write(*, *) "  -fa, --factr    Scaling factor applied to merged data."
      write(*, *) "   -d, --dep      Expected depth."
      write(*, *) ""
      write(*, *) "Optional arguments:"
      write(*, *) "   -h, --help     Show this help message and exit."
      write(*, *) "   -c, --chr      First chromosome in model organism"
      write(*, *) "                  (default: chr1)."
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
         write(*, *) "Error: Directory does not exist: ", trim(dir_path)
         stop
      endif
   end subroutine check_exists_dir

end program