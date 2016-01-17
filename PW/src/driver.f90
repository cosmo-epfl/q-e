   SUBROUTINE driver()
    USE basis,            ONLY : startingconfig, starting_pot, starting_wfc
    USE input_parameters, ONLY : restart_mode
    USE io_global,        ONLY : ionode, ionode_id
    USE io_files,         ONLY : srvaddress, seqopn, delete_if_present
    USE mp_global,        ONLY : mp_startup, mp_bcast, mp_global_end, intra_image_comm
    USE control_flags,    ONLY : conv_elec, conv_ions, restart, lbfgs, lmd
    USE extrapolation,    ONLY : update_file, update_pot
    USE ions_base,              ONLY : tau
    USE cell_base,              ONLY : alat, at, omega, bg
    USE force_mod,              ONLY : force
    USE ener,                   ONLY : etot
    USE f90sockets,               ONLY: open_socket, readbuffer, writebuffer
    USE cellmd,                 ONLY : lmovecell, at_old, omega_old    
    USE fft_base,               ONLY : dfftp
    USE fft_base,               ONLY : dffts
    USE gvect,                  ONLY : gcutm
    USE gvecs,                  ONLY : gcutms
    USE grid_subroutines,       ONLY : realspace_grid_init
  USE read_input,        ONLY : read_input_file
      USE wvfct,           ONLY : ecutwfc
    IMPLICIT NONE
     
    INTEGER, PARAMETER :: MSGLEN=12
    LOGICAL :: isinit=.false., hasdata=.false., firststep=.true.
    CHARACTER*12 :: header
    CHARACTER*1024 :: parbuffer, host
    INTEGER :: socket, nat, inet, port, replicaid=-1, parbufflen
    REAL*8 :: sigma(3,3)
    REAL *8 :: cellh(3,3), cellih(3,3), vir(3,3), pot, mtxbuffer(9)
    REAL*8, ALLOCATABLE :: combuf(:)
    INTEGER::CALLS=1, str_index, str_index2
    REAL*8:: dist_ang(6),dist_ang_old(6),dist_ang_reset(6)
!Conditions to reset the cell
    LOGICAL:: lang_reset=.false., lang_old=.false., ldist_reset=.false., &
              ldist_old=.false., lvol_old=.false., lvol_reset=.false.
    REAL*8 :: ang_tol=10.d0, dist_tol=0.05d0, vol_tol=0.05d0
    REAL*8 :: omega_reset=0.0d0

!Ecutwfc
    REAL*8:: ecutwfc_orig
 
    ecutwfc_orig = ecutwfc
    lbfgs = .true.
    lmd = .true.
    lvol_old   = .false.
    lvol_reset = .false.
    ldist_old  = .false.
    ldist_reset= .false.
    lang_old   = .false.
    lang_reset = .false.
    ! parses host name, port and socket type, from input (format:srvaddress = [hostname]:[UNIX|port]
    inet=1
    host=srvaddress(1:INDEX(srvaddress,':',back=.true.)-1)//achar(0)
    read(srvaddress(INDEX(srvaddress,':',back=.true.)+1:),*) port

    IF (srvaddress(1:INDEX(srvaddress,':')-1).eq.('UNIX')) THEN
      inet=0
      host=srvaddress(6:INDEX(srvaddress,':',back=.true.)-1)//achar(0)    
    ENDIF

    IF (ionode) write(*,*) " @ DRIVER MODE: Connecting to host:port ", trim(srvaddress(1:INDEX(srvaddress,':',back=.true.)-1)), port

    ! opens socket and starts main loop
    IF (ionode) call open_socket(socket, inet, port, host)            

    driver_loop: DO
      ! since often input configuration can be different from what is passed first from i-PI, it is better to start from random atomic
      starting_pot = "atomic"
      starting_wfc = 'atomic+random'
      startingconfig = 'input'
      restart = .false.
      restart_mode = "from_scratch"

      ! do communication on master node only...
      if (ionode) call readbuffer(socket, header, MSGLEN)
      call mp_bcast(header,ionode_id, intra_image_comm)
      
      if (ionode) write(*,*) " @ DRIVER MODE: Message from server: ", trim(header)

      if (trim(header) == "STATUS") then
         if (ionode) then  ! does not  need init (well, maybe it should, just to check atom numbers and the like... )
            if (hasdata) then
               call writebuffer(socket,"HAVEDATA    ",MSGLEN)               
            else if (isinit) then
               call writebuffer(socket,"READY       ",MSGLEN)
            else 
               call writebuffer(socket,"NEEDINIT    ",MSGLEN)             
            endif
         endif
      else if (trim(header) == "INIT") then
         if (ionode) call readbuffer(socket, nat) ! actually this reads the replica id when the driver is used for multiple independent runs or replicas
         call mp_bcast(nat,ionode_id, intra_image_comm)
         if (nat.ne.replicaid .and. .not. firststep) then
           call delete_if_present('update') ! resets file history 
           call update_file()
         endif
         if (ionode) write(*,*) " @ DRIVER MODE: Receiving replica", nat, "old: ", replicaid
         replicaid = nat

         if (ionode) then 
           call readbuffer(socket, parbufflen) ! length of parameter string
           call readbuffer(socket, parbuffer, parbufflen)
           write(*,*)" @ DRIVER MODE: parbuffer ", parbuffer(1:parbufflen)
           str_index = VERIFY(parbuffer(1:parbufflen),":",.false.)
           str_index2= VERIFY(parbuffer(1:parbufflen),",",.true.)
           write(*,*)" @ DRIVER MODE: message ", trim(parbuffer(str_index:str_index2))
         endif
         call mp_bcast(parbufflen,ionode_id, intra_image_comm)
         call mp_bcast(parbuffer, ionode_id, intra_image_comm)

         isinit=.true.
      else if (trim(header) == "POSDATA") then

         ! receives the positions & the cell data
         ! first reads cell and the number of atoms
         if (ionode) call readbuffer(socket, mtxbuffer, 9)
         cellh = RESHAPE(mtxbuffer, (/3,3/))         
         if (ionode) call readbuffer(socket, mtxbuffer, 9)
         cellih = RESHAPE(mtxbuffer, (/3,3/))
         if (ionode) call readbuffer(socket, nat)
  
         call mp_bcast(cellh,ionode_id, intra_image_comm)  ! must communicate to other nodes
         call mp_bcast(cellih,ionode_id, intra_image_comm)
         call mp_bcast(nat,ionode_id, intra_image_comm)
         if (.not.allocated(combuf)) then
           allocate(combuf(3*nat))
         end if
         if (ionode) call readbuffer(socket, combuf, nat*3)
         call mp_bcast(combuf,ionode_id, intra_image_comm)
                  
         at_old = at
         omega_old = omega
         dist_ang_old=dist_ang
         ! convert the incoming configuration to the internal pwscf format
         cellh=transpose(cellh)                       ! row-major to column-major 
         cellih=transpose(cellih)         
         tau = RESHAPE(combuf, (/ 3 , nat /) )/alat   ! internally positions are in alat 
         at = cellh / alat                            ! and so the cell

         !Check here how much the cell has changed since last and reset
         if(.not.firststep) then
           call dist_latvec2ang(dist_ang,cellh)
           lvol_old   =((omega_old-omega)/omega.gt.vol_tol)
           lvol_reset =((omega_reset-omega)/omega.gt.vol_tol)
           ldist_old  =any(abs(dist_ang(1:3)-dist_ang_old(1:3))  /maxval(dist_ang(1:3)).gt.dist_tol)
           ldist_reset=any(abs(dist_ang(1:3)-dist_ang_reset(1:3))/maxval(dist_ang(1:3)).gt.dist_tol)
           lang_old   =any(abs(dist_ang(4:6)-dist_ang_old(4:6))  .gt.ang_tol)
           lang_reset =any(abs(dist_ang(4:6)-dist_ang_reset(4:6)).gt.ang_tol)
         endif
                  
         if (ionode) write(*,*) " @ DRIVER MODE: Received positions "

         lmovecell=.TRUE. !This is true for most cases... should be eventually read from input files
         CALL recips( at(1,1),at(1,2),at(1,3), bg(1,1),bg(1,2),bg(1,3) )
         CALL volume( alat, at(1,1), at(1,2), at(1,3), omega )
         if (firststep) then ! do run initialisation here, so it is done with the positions sent from i-PI
            at_old = at
            omega_old = omega
            omega_reset=omega
            dist_ang_old=dist_ang
            dist_ang_reset=dist_ang
            CALL setup ()
            if (ionode) write(*,*) " @ DRIVER MODE: Preparing first evaluation "
            lmovecell=.TRUE.
            ! if (ionode) call system("rm -rf pw*") SERIOUSLY????
            call init_run()
            firststep=.false.
         else
             if(ldist_old.or.ldist_reset.or.lang_old.or.lang_reset.or.lvol_old.or.lvol_reset) then
                   ! ... Variable-cell optimization: if cell changes too much 
                   ! ... reinitialize the calculation with G-vectors and plane waves
                   ! ... calculated for the new cell (may differ from the curent
                   ! ... result, using G_vectors and PWs for the starting cell)
                   !
                   if (ionode) write(*,*) " @ DRIVER MODE: reset evaluation "
                     if (ionode) then
                         if(ldist_old) write(*,'(a,3es15.6)') &
     " @ DRIVER MODE: reset by cell length change from previous step : ", & 
     abs(dist_ang(1:3)-dist_ang_old(1:3))  /maxval(dist_ang(1:3)) 
                         if(ldist_reset) write(*,'(a,3es15.6)') &
     " @ DRIVER MODE: reset by cell length change from previous reset: ", &
     abs(dist_ang(1:3)-dist_ang_reset(1:3))/maxval(dist_ang(1:3))
                         if(lang_old) write(*,'(a,3es15.6)') &
     " @ DRIVER MODE: reset by cell angle change from previous step  : ", &
     abs(dist_ang(4:6)-dist_ang_old(4:6))
                         if(lang_reset)  write(*,'(a,3es15.6)') &
     " @ DRIVER MODE: reset by cell angle change from previous reset : ", &
     abs(dist_ang(4:6)-dist_ang_reset(4:6))
                         if(lvol_old) write(*,'(a,3es15.6)') &
     " @ DRIVER MODE: reset by cell vol   change from previous step  : ", &
     abs(dist_ang(4:6)-dist_ang_reset(4:6))
                         if(lvol_reset)  write(*,'(a,3es15.6)') &
     " @ DRIVER MODE: reset by cell vol   change from previous reset : ", &
     abs(dist_ang(4:6)-dist_ang_reset(4:6))
                     endif
                     
                     call delete_if_present('update') ! resets file history
                     call update_file( )
                     at_old = at
                     omega_old = omega
                     omega_reset = omega
                     dist_ang_reset=dist_ang
          
                     CALL clean_pw( .false. )
                     CALL close_files(.true.)
                     call setup()
                     lmovecell=.true.
                     lbfgs=.true.
                     lmd=.true.
                     conv_ions = .false.
                     dfftp%nr1=0; dfftp%nr2=0; dfftp%nr3=0; dffts%nr1=0; dffts%nr2=0; dffts%nr3=0
                     CALL realspace_grid_init (dfftp, at, bg, gcutm )
                     CALL realspace_grid_init (dfftp, at, bg, gcutms )
                     ! if (ionode) call system("rm -rf pw*")
                     CALL init_run()
             else
                     lmovecell=.TRUE.
                     CALL hinit1()
             endif
         end if 

         ! does the actual SCF calculation
         CALL electrons()
         IF ( .NOT. conv_elec ) THEN
           CALL punch( 'all' )
           CALL stop_run( conv_elec )
         ENDIF         
         CALL forces()
         CALL stress(sigma)
         CALLS=CALLS+1
         
         ! converts energy & forces to the format expected by i-pi (so go from Ry to Ha)
         combuf=RESHAPE(force, (/ 3 * nat /) ) * 0.5   ! return force in atomic units    
         pot=etot * 0.5   ! return potential in atomic units
         vir=transpose(sigma)*omega*0.5   ! return virial in atomic units and without the volume scaling         

         ! updates history
         call update_file()
         call update_pot()
         
         hasdata=.true.
      else if (trim(header)=="GETFORCE") then
         ! communicates energy info back to i-pi
         if (ionode) write(*,*) " @ DRIVER MODE: Returning v,forces,stress "
                               
         if (ionode) call writebuffer(socket,"FORCEREADY  ",MSGLEN)         
         if (ionode) call writebuffer(socket,pot)
         if (ionode) call writebuffer(socket,nat)            
         if (ionode) call writebuffer(socket,combuf,3*nat)
         if (ionode) call writebuffer(socket,RESHAPE(vir,(/9/)),9)
         
         ! i-pi can also receive an arbitrary string, that will be printed out to the "extra" 
         ! trajectory file. this is useful if you want to return additional information, e.g.
         ! atomic charges, wannier centres, etc. one must return the number of characters, then
         ! the string. here we just send back zero characters.
         nat = 0
         if (ionode) call writebuffer(socket,nat)
              
         isinit = .false. ! resets init so that we will get replica index again at next step!
         hasdata=.false.
         CALL punch( 'config' )
      endif
    ENDDO driver_loop    
    conv_ions=.true.
contains
!************************************************************************************

  subroutine dist_latvec2ang(dist_ang,latvec)
    !This subroutine will generate the angdeg represenation of the cell from the lattice vectors
    implicit none
     real(8):: dist_ang(6),latvec(3,3),pi,convang
     pi=acos(-1.d0) 
     convang=180.d0/pi
     dist_ang(1)=sqrt(dot_product(latvec(:,1),latvec(:,1)))
     dist_ang(2)=sqrt(dot_product(latvec(:,2),latvec(:,2)))
     dist_ang(3)=sqrt(dot_product(latvec(:,3),latvec(:,3)))
     dist_ang(4)=acos(dot_product(latvec(:,2),latvec(:,3))/(dist_ang(2)*dist_ang(3)))*convang
     dist_ang(5)=acos(dot_product(latvec(:,3),latvec(:,1))/(dist_ang(3)*dist_ang(1)))*convang
     dist_ang(6)=acos(dot_product(latvec(:,1),latvec(:,2))/(dist_ang(1)*dist_ang(2)))*convang
  end subroutine
!************************************************************************************

    
  END SUBROUTINE
