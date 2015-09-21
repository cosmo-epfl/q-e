!
! Copyright (C) 2005 Andrea Ferretti
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
PROGRAM sumpdos
  !
  IMPLICIT NONE
  !
  ! AUTHOR: Andrea Ferretti
  !
  ! this program reads and sum pdos from different
  ! files (which are related to different atoms)
  !
  ! file names are read from stdin
  ! USAGE: sumpdos <file1> ... <fileN>
  !
  INTEGER, EXTERNAL   :: i_argc             ! function giving no of arguments
  INTEGER             :: ngrid              ! dimension of the energy grid
  INTEGER             :: nfile              ! number of files to sum
  INTEGER             :: nspin              ! number of spin_component


  CHARACTER(256), ALLOCATABLE    :: file(:) ! names of the files to sum
  CHARACTER(256)      :: filein
  CHARACTER(10)       :: cdum, str1, str2

  LOGICAL             :: exist
  REAL                :: efermi = 0.0d0       ! translate the input grid
  REAL, ALLOCATABLE   :: pdos(:,:,:)
  REAL, ALLOCATABLE   :: egrid(:)
  REAL, ALLOCATABLE   :: mysum(:,:)

  INTEGER :: ios, ierr, iarg, ie, isp, ifile, i


!**************************************************************
! User should supply input values here
!
efermi = 0.0d0

!**************************************************************

!
! get the number of arguments (i.e. the number of files)
!
   nfile = i_argc ()
   IF ( nfile == 0 ) THEN
      WRITE(0,"( 'No file to sum' )")
      STOP
   ENDIF

   CALL get_arg ( 1, str1 )
   !
   SELECT CASE ( trim(str1) )
   CASE ( "-h" )
      !
      ! write the manual
      !
      WRITE(0,"(/,'USAGE: sumpdos [-h] [-f <filein>] [<file1> ... <fileN>]', /, &
                &'  Sum the pdos from the file specified in input and write the sum ', /, &
                &'  to stdout', /, &
                &'     -h           : write this manual',/, &
                &'     -f <filein>  : takes the list of pdos files from <filein> ', /, &
                &'                    (one per line) instead of command line',/, &
                &'     <fileM>      : the M-th pdos file', &
                & / )")
      STOP
      !
   CASE ( "-f" )
      !
      ! read file names from file
      !
      CALL get_arg ( 2, filein )
      IF ( len_trim(filein) == 0 ) CALL errore('sumpdos','provide filein name',2)

      INQUIRE( FILE=trim(filein), EXIST=exist )
      IF (.not. exist) CALL errore('sumpdos','file '//trim(filein)//' does not exist',3)
      OPEN( 10, FILE=trim(filein), IOSTAT=ios )
      IF (ios/=0) CALL errore('sumpdos','opening '//trim(filein),abs(ios))

      !
      ! get the number of non-empty lines in the file
      ! (which is assumed to be the number of files to sum)
      !
      ios = 0
      nfile = 0
      !
      DO WHILE ( ios == 0 )
         nfile = nfile + 1
         READ(10, *, IOSTAT=ios ) cdum
         IF ( ios ==0 .and. len_trim(cdum)==0 ) nfile = nfile -1
      ENDDO
      nfile = nfile -1

      !
      IF (nfile ==0 ) CALL errore('sumpdos','no file to sum in '//trim(filein),4)
      !
      ALLOCATE( file(nfile), STAT=ierr )
      IF (ierr/=0) CALL errore('sumpdos','allocating FILE',abs(ierr))
      !
      REWIND(10)

      DO i = 1, nfile
         file(i) = ' '
         DO WHILE( len_trim(file(i)) == 0 )
            READ(10,*, IOSTAT=ios) file(i)
            IF (ios /=0 ) CALL errore('sumpdos','reading from '//trim(filein),i)
         ENDDO
      ENDDO

   CASE DEFAULT

      !
      ! get the names of the files
      ! here we use GETARG
      !
      ALLOCATE( file(nfile), STAT=ierr )
      IF (ierr/=0) CALL errore('sumpdos','allocating FILE',abs(ierr))
      DO iarg = 1, nfile
         CALL get_arg ( iarg, file(iarg) )
      ENDDO

   END SELECT

!
! open the first file and get data about spin
! and grid dimensions
!
   INQUIRE( FILE=trim(file(1)), EXIST=exist )
   IF (.not. exist) CALL errore('sumpdos','file '//trim(file(1))//' does not exist',3)
   !
   WRITE(0,"('Reading dimensions from file: ',a)") trim(file(1))
   !
   OPEN(10, FILE=trim(file(1)), IOSTAT=ios)
      IF (ios/=0) CALL errore("sumpdos", "error opening "//trim(file(1)), 1)
      !
      ! try to understand if we have 1 or 2 spin
      !
      READ(10,*, IOSTAT=ios) cdum, cdum, cdum, str1, str2
      IF (ios/=0) CALL errore("sumpdos", "reading first line of "//trim(file(1)), 1)
      !
      IF ( trim(str1) == 'ldos(E)' ) THEN
          nspin = 1
      ELSEIF ( trim(str1) == 'ldosup(E)' .and.  trim(str2) == 'ldosdw(E)' ) THEN
          nspin = 2
      ELSE
          CALL errore("sumpdos", "wrong fmf in the first line of "//trim(file(1)), 1)
      ENDIF
      !
      ! determine the dimension fo the energy mesh
      ! no further control will be done on the consistency of the energy
      ! grid of each file
      !
      ie = 0
      DO WHILE ( .true.  )
         READ( 10, *, IOSTAT=ios )
         IF ( ios /= 0 ) exit
         ie = ie + 1
      ENDDO
      ngrid = ie

   CLOSE(10)

!
! allocations
!
   ALLOCATE( pdos( ngrid, nspin, nfile), STAT=ierr )
      IF (ierr/=0) CALL errore("sumpdos", "allocating pdos", ierr)
   ALLOCATE( mysum( ngrid, nspin), STAT=ierr )
      IF (ierr/=0) CALL errore("sumpdos", "allocating mysum", ierr)
   ALLOCATE( egrid( ngrid) )
      IF (ierr/=0) CALL errore("sumpdos", "allocating egrid", ierr)


!
! get data
!
   WRITE(0,"('Reading the following ',i5,' files: ')") nfile
   !
   DO ifile = 1, nfile
      !
      INQUIRE( FILE=trim(file(ifile)), EXIST=exist )
      IF (.not. exist) &
         CALL errore('sumpdos','file '//trim(file(ifile))//' does not exist',ifile)
      !
      WRITE(0,"(2x,'Reading file: ',a)") trim(file(ifile))
      OPEN(10, FILE=trim(file(ifile)), IOSTAT=ios)
         IF (ios/=0) CALL errore("sumpdos", "error opening "//trim(file(ifile)), ios )
         !
         READ(10,*, IOSTAT=ios)
         IF (ios/=0) &
            CALL errore("sumpdos", "reading first line in "//trim(file(ifile)), ios )
         !
         ! egrid is overwritten every time
         !
         DO ie = 1, ngrid
            READ(10, *, IOSTAT=ios ) egrid(ie), pdos(ie, 1:nspin, ifile)
            IF (ios/=0) &
            CALL errore("sumpdos", "reading first line in "//trim(file(ifile)), ie )
         ENDDO
     CLOSE(10)
   ENDDO

!
! perform the sum and write
!
   IF ( nspin == 1 ) THEN
      WRITE(6,"('# E (eV)  pdos(E) ')")
   ELSEIF ( nspin == 2) THEN
      WRITE(6,"('# E (eV)  pdos_UP(E)    pdos_DW(E) ')")
   ELSE
      CALL errore("sunpdos", "really sure NSPIN /= 1 or 2 ???", 3 )
   ENDIF

   mysum = 0.0d0
   DO ie=1,ngrid
      DO isp=1,nspin
         mysum(ie,isp) = sum( pdos(ie,isp,:) )
      ENDDO
      WRITE(6,"(3f15.9)") egrid(ie) - efermi, mysum(ie,1:nspin)
   ENDDO

!
! clean
!
   DEALLOCATE( file, STAT=ierr )
      IF (ierr/=0) CALL errore("sumpdos", "deallocating file", ierr)
   DEALLOCATE( pdos, STAT=ierr )
      IF (ierr/=0) CALL errore("sumpdos", "deallocating pdos", ierr)
   DEALLOCATE( mysum, STAT=ierr )
      IF (ierr/=0) CALL errore("sumpdos", "deallocating mysum", ierr)
   DEALLOCATE( egrid, STAT=ierr )
      IF (ierr/=0) CALL errore("sumpdos", "deallocating egrid", ierr)

CONTAINS

!*************************************************
SUBROUTINE errore(routine, msg, ierr)
   !*************************************************
   IMPLICIT NONE
   CHARACTER(*),    INTENT(in) :: routine, msg
   INTEGER,         INTENT(in) :: ierr

   !
   WRITE( UNIT = 0, FMT = '(/,1X,78("*"))')
   WRITE( UNIT = 0, &
          FMT = '(5X,"from ",A," : error #",I10)' ) routine, ierr
   WRITE( UNIT = 0, FMT = '(5X,A)' ) msg
   WRITE( UNIT = 0, FMT = '(1X,78("*"),/)' )
   !
   STOP
   RETURN
END SUBROUTINE errore

END PROGRAM sumpdos

