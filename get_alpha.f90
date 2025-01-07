program get_alpha

   ! Rough documentation:
   !   'mergetracks' processes IP and input BED files to generate 'merged'
   !   (IP ÷ input) siQ-scaled data.

   character(len = 256) :: fil_prm
   real(8) :: dep_ip, dep_in
   real(8), dimension(24) :: prm
   real(8) :: empir
   logical :: file_exists
   integer :: i, iostat

   ! Call argument parser
   call parse_args(fil_prm, dep_ip, dep_in)

   ! ! Debug or verify parsed arguments (optional)
   ! write(*, *) "####################################"
   ! write(*, *) "## Parsed arguments for get_alpha ##"
   ! write(*, *) "####################################"
   ! write(*, *) ""
   ! write(*, *) "fil_prm: ", trim(fil_prm)
   ! write(*, *) "dep_ip:  ", dep_ip
   ! write(*, *) "dep_in:  ", dep_in
   ! write(*, *) ""
   ! write(*, *) ""

   ! Check that file exists
   inquire(file = fil_prm, exist = file_exists)
   if (.not. file_exists) then
      write(6, *) "Error: The provided file path does not exist: ", trim(fil_prm)
      stop 1
   endif

   ! Open file and read parameters
   open(unit = 33, file = fil_prm, status = "old", iostat = iostat)
   if (iostat /= 0) then
      write(6, *) "Error: Could not open file: ", trim(fil_prm)
      stop 1
   endif

   ! Read parameter array
   do i = 1, 6
      read(33, *, iostat = iostat) prm(i)
      if (iostat /= 0) then
         write(6, *) "Error: Failed to read parameter at index ", i
         stop
      endif
   end do
   close(33)

   ! Compute empirical scaling factor
   empir = compute_empir(prm, dep_ip, dep_in)

   ! Output the result
   write(*, *) empir

contains

   ! Subroutine to parse command-line arguments
   subroutine parse_args(fil_prm, dep_ip, dep_in)
      character(len = 256), intent(out) :: fil_prm
      real(8), intent(out) :: dep_ip, dep_in
      character(len = 256) :: arg
      character(len = 256) :: tmp_arg
      integer :: i

      ! Initialize defaults
      fil_prm = ""
      dep_ip = -1.0d0
      dep_in = -1.0d0

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
         elseif (arg(1:4) == "-fp=" .or. arg(1:10) == "--fil_prm=") then
            fil_prm = trim(adjustl(arg(index(arg, "=") + 1:)))
         elseif (arg(1:4) == "-di=" .or. arg(1:9) == "--dep_ip=") then
            tmp_arg = trim(adjustl(arg(index(arg, "=") + 1:)))
            read(tmp_arg, *) dep_ip
         elseif (arg(1:4) == "-dn=" .or. arg(1:9) == "--dep_in=") then
            tmp_arg = trim(adjustl(arg(index(arg, "=") + 1:)))
            read(tmp_arg, *) dep_in
         else
            write(6, *) "Error: Unrecognized argument: ", arg
            stop 1
         endif
      end do

      ! Check required arguments
      if (fil_prm == "") then
         write(6, *) "Error: Missing required argument -fp/--fil_prm."
         stop 1
      endif

      if (dep_ip <= 0.0d0) then
         write(6, *) "Error: Invalid or missing value for -di/--dep_ip."
         stop 1
      endif

      if (dep_in <= 0.0d0) then
         write(6, *) "Error: Invalid or missing value for -dn/--dep_in."
         stop 1
      endif
   end subroutine parse_args

   ! Subroutine to print help message
   subroutine print_help()
      write(*, *) "Usage:"
      write(*, *) "  get_alpha"
      write(*, *) "    --fil_prm=<str> --dep_ip=<flt> --dep_in=<flt>"
      write(*, *) ""
      write(*, *) "Required arguments:"
      write(*, *) "  -fp, --fil_prm   Path to the input parameter file."
      write(*, *) "  -di, --dep_ip    Depth for IP."
      write(*, *) "  -dn, --dep_in    Depth for input."
      write(*, *) ""
      write(*, *) "Optional arguments:"
      write(*, *) "   -h, --help        Show this help message and exit."
   end subroutine print_help

   ! Function to compute empirical scaling factor
   real(8) function compute_empir(prm, dep_ip, dep_in)
      real(8), intent(in) :: prm(6)
      real(8), intent(in) :: dep_ip, dep_in
      compute_empir = (prm(1) / (prm(2) - prm(1))) &
         * (prm(4) / prm(3)) &
         * (dep_ip / dep_in) &
         * (prm(6) / prm(5))
   end function compute_empir

end program get_alpha