   SUBROUTINE driver()
    USE io_global,        ONLY : stdout, ionode, ionode_id
    USE io_files,         ONLY : srvaddress, iunupdate, seqopn
    USE mp_global,        ONLY : mp_startup, mp_bcast, mp_global_end, intra_image_comm
    USE control_flags,    ONLY : conv_elec, istep, conv_ions, &
                                 gvec_omega_tol, gvec_ang_tol, gvec_dist_tol
    
    USE ions_base,              ONLY : tau, ityp
    USE cell_base,              ONLY : alat, at, omega, bg
    USE cellmd,                 ONLY : omega_old, at_old, lmovecell, calc
    USE force_mod,              ONLY : force
    USE ener,                   ONLY : etot
    USE f90sockets,             ONLY : open_socket, readbuffer, writebuffer
    USE fft_base,               ONLY : dfftp
    USE fft_base,               ONLY : dffts
    ! USE grid_subroutines,       ONLY : realspace_grids_init
    USE gvect,                  ONLY : gcutm
    USE gvecs,                  ONLY : gcutms
    USE symm_base,              ONLY : checkallsym
    USE extrapolation,          ONLY : update_file, update_pot
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: MSGLEN=12
    LOGICAL :: isinit=.false., hasdata=.false., firststep=.true., exst, lgreset
    CHARACTER*12 :: header
    CHARACTER*1024 :: parbuffer, host
    INTEGER :: socket, nat, inet, port, ccmd, i, info, replicaid=-1
    REAL*8 :: sigma(3,3), omega_reset, at_reset(3,3), dist_reset, ang_reset
    REAL *8 :: cellh(3,3), cellih(3,3), vir(3,3), pot, mtxbuffer(9)
    REAL*8, ALLOCATABLE :: combuf(:)!, tauhist(:,:,:)
    REAL*8 :: dist_ang(6), dist_ang_reset(6)

    lmovecell = .true.
    omega_reset = .0
    dist_ang_reset = .0
    omega_old = .0
    at_old = .0
    
    ! parses host name, port and socket type
    inet=1
    host=srvaddress(1:INDEX(srvaddress,':',back=.true.)-1)//achar(0)
    read(srvaddress(INDEX(srvaddress,':',back=.true.)+1:),*) port

    IF (srvaddress(1:INDEX(srvaddress,':')-1).eq.('UNIX')) THEN
      inet=0
      host=srvaddress(6:INDEX(srvaddress,':',back=.true.)-1)//achar(0)    
    ENDIF

    IF (ionode) write(*,*) " @ DRIVER MODE: Connecting to host:port ", trim(srvaddress(6:INDEX(srvaddress,':',back=.true.)-1)), port
    ! opens socket and starts main loop
    IF (ionode) call open_socket(socket, inet, port, host)

    ! Start main socket loop
    driver_loop: DO
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

          if (ionode) call readbuffer(socket, nat) ! actually this reads the replica id         
          call mp_bcast(nat,ionode_id, intra_image_comm)

          ! Here nat is the replicaid
          if (nat.ne.replicaid .and. .not. firststep) then
             ! history = 1 ! resets history -- want to do new-old propagation if replica changed
             ! if (ionode) then !Those files are used for extraplolation
             !    CALL seqopn( iunupdate, 'update', 'FORMATTED', exst )
             !    WRITE( UNIT = iunupdate, FMT = * ) history
             !    WRITE( UNIT = iunupdate, FMT = * ) tauhist
             !    CLOSE( UNIT = iunupdate, STATUS = 'KEEP' )
             ! endif
             ! new way of resetting history (check if something better exists).
            if(firststep) CALL setup()
!            CALL realspace_grids_init (dfftp, dffts,at, bg, gcutm, gcutms )
            ! CALL checkallsym( nat, tau, ityp, dfftp%nr1, dfftp%nr2, dfftp%nr3 )
            ! if(ionode) CALL system("rm -rf dia*")
            
            CALL clean_pw(.FALSE.)
            if(.not.firststep) CALL close_files(.TRUE.)
            firststep = .false.
          endif

          if (ionode) write(*,*) " @ DRIVER MODE: Receiving replica", nat, replicaid
          replicaid = nat

          if (ionode) then 
             call readbuffer(socket, nat) ! length of parameter string -- ignored at present!
             call readbuffer(socket, parbuffer, nat)
          endif
          
         isinit=.true.

      else if (trim(header) == "POSDATA") then

         if(.not.firststep) then
            at_old = at
            omega_old = omega
         end if
         
         ! receives the positions & the cell data
         ! first reads cell and the number of atoms
         if (ionode) call readbuffer(socket, mtxbuffer, 9)
         cellh = RESHAPE(mtxbuffer, (/3,3/))         
         if (ionode) call readbuffer(socket, mtxbuffer, 9)
         cellih = RESHAPE(mtxbuffer, (/3,3/))
         ! Here nat is the number of atoms
         if (ionode) call readbuffer(socket, nat)
  
         call mp_bcast(cellh,ionode_id, intra_image_comm)  ! must communicate to other nodes
         call mp_bcast(cellih,ionode_id, intra_image_comm)
         call mp_bcast(nat,ionode_id, intra_image_comm)
         IF (.not.allocated(combuf)) THEN
           ALLOCATE(combuf(3*nat))
         END IF
         if (ionode) call readbuffer(socket, combuf, nat*3)
         call mp_bcast(combuf,ionode_id, intra_image_comm)

         IF( .not.firststep ) then
            at_old = at
            omega_old = omega
         END IF
                  
         ! convert the incoming configuration to the internal pwscf format
         cellh=transpose(cellh)                       ! row-major to column-major 
         cellih=transpose(cellih)         
         tau = RESHAPE(combuf, (/ 3 , nat /) )/alat   ! internally positions are in alat 
         at = cellh / alat                            ! and so the cell
                  
         IF (ionode) WRITE(*,*) " @ DRIVER MODE: Received positions "
         
         CALL recips( at(1,1),at(1,2),at(1,3), bg(1,1),bg(1,2),bg(1,3) )
         CALL volume( alat, at(1,1), at(1,2), at(1,3), omega )

         CALL cell2dist_ang(dist_ang, cellh)
         lgreset = ((abs(omega_reset - omega)/omega .gt. gvec_omega_tol) .or. &
              & (any(abs(dist_ang(1:3) - dist_ang_reset(1:3)) / maxval(dist_ang(1:3)) .gt. gvec_dist_tol)) .or. &
              & (any(abs(dist_ang(4:6) - dist_ang_reset(4:6)) .gt. gvec_ang_tol)))

         ! refresh the cell stuff
         IF (firststep .or. ( lmovecell .and. lgreset )) THEN
            ! changes needed only if cell moves

            IF (ionode) write(*,*) " @ DRIVER MODE: reinitialize G-vectors "
            
            at_reset = at
            omega_reset = omega
            dist_ang_reset = dist_ang
                     
            IF (firststep) CALL setup()
            
            CALL clean_pw(.FALSE.)
            IF (.not.firststep) CALL close_files(.TRUE.)
            firststep = .false.
            
            CALL init_run()
            
            CALL mp_bcast( at,        ionode_id, intra_image_comm )
            CALL mp_bcast( at_old,    ionode_id, intra_image_comm )
            CALL mp_bcast( omega,     ionode_id, intra_image_comm )
            CALL mp_bcast( omega_old, ionode_id, intra_image_comm )
            CALL mp_bcast( bg,        ionode_id, intra_image_comm )
            !
         ELSE
            CALL update_pot()
            CALL hinit1()
         END IF

                  
         ! Compute everything !
         CALL electrons()
         IF ( .NOT. conv_elec ) THEN
           CALL punch( 'all' )
           CALL stop_run( conv_elec )
         ENDIF         
         CALL forces()
         CALL stress(sigma)
         
         ! converts energy & forces to the format expected by i-pi (so go from Ry to Ha)
         combuf=RESHAPE(force, (/ 3 * nat /) ) * 0.5   ! return force in atomic units    
         pot=etot * 0.5                                ! return potential in atomic units
         vir=transpose(sigma)*omega*0.5                ! return virial in atomic units and without the volume scaling         

         ! updates history
         istep = istep+1
         CALL update_file()
                         
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

  CONTAINS

    SUBROUTINE cell2dist_ang(dist_ang, latvec)
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
    END SUBROUTINE cell2dist_ang
    
  END SUBROUTINE
