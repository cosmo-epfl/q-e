
   SUBROUTINE driver()
    USE io_global,        ONLY : stdout, ionode, ionode_id
    USE io_files,         ONLY : srvaddress, iunupdate, seqopn
    USE mp_global,        ONLY : mp_startup, mp_bcast, mp_global_end, intra_image_comm
    USE control_flags,    ONLY : conv_elec, history, istep, conv_ions
    
    USE ions_base,              ONLY : tau
    USE cell_base,              ONLY : alat, at, omega
    USE force_mod,              ONLY : force
    USE ener,                   ONLY : etot
    USE f90sockets,               ONLY: open_socket, readbuffer, writebuffer
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: MSGLEN=12
    LOGICAL :: isinit=.false., hasdata=.false., firststep=.true., exst
    CHARACTER*12 :: header
    CHARACTER*1024 :: parbuffer, host
    INTEGER :: socket, nat, inet, port, ccmd, i, info, replicaid=-1
    REAL*8 :: sigma(3,3)
    REAL *8 :: cellh(3,3), cellih(3,3), vir(3,3), pot, mtxbuffer(9)
    REAL*8, ALLOCATABLE :: combuf(:), tauhist(:,:,:)

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
         if (nat.ne.replicaid .and. .not. firststep) then
           history = 1 ! resets history -- want to do new-old propagation if replica changed
           if (ionode) then
             CALL seqopn( iunupdate, 'update', 'FORMATTED', exst )
             WRITE( UNIT = iunupdate, FMT = * ) history
             WRITE( UNIT = iunupdate, FMT = * ) tauhist
             CLOSE( UNIT = iunupdate, STATUS = 'KEEP' )
           endif
         endif
         if (ionode) write(*,*) " @ DRIVER MODE: Receiving replica", nat, replicaid, history
         replicaid = nat

         if (ionode) then 
           call readbuffer(socket, nat) ! length of parameter string -- ignored at present!
           call readbuffer(socket, parbuffer, nat)
         endif
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
           allocate(tauhist(3,nat,3) )
           tauhist=0.d0
         end if
         if (ionode) call readbuffer(socket, combuf, nat*3)
         call mp_bcast(combuf,ionode_id, intra_image_comm)
                  
         ! convert the incoming configuration to the internal pwscf format
         cellh=transpose(cellh)                       ! row-major to column-major 
         cellih=transpose(cellih)         
         tau = RESHAPE(combuf, (/ 3 , nat /) )/alat   ! internally positions are in alat 
         at = cellh / alat                            ! and so the cell
                  
         if (ionode) write(*,*) " @ DRIVER MODE: Received positions "

         if (firststep) then ! do run initialisation here, so it is done with the positions sent from i-PI
            if (ionode) write(*,*) " @ DRIVER MODE: Preparing first evaluation "
            call init_run()
            firststep=.false.
         end if 

         CALL hinit1()
         CALL electrons()
         IF ( .NOT. conv_elec ) THEN
           CALL punch( 'all' )
           CALL stop_run( conv_elec )
         ENDIF         
         CALL forces()
         CALL stress(sigma)
         
         ! converts energy & forces to the format expected by i-pi (so go from Ry to Ha)
         combuf=RESHAPE(force, (/ 3 * nat /) ) * 0.5   ! return force in atomic units    
         pot=etot * 0.5   ! return potential in atomic units
         vir=transpose(sigma)*omega*0.5   ! return virial in atomic units and without the volume scaling         

         ! updates history
         istep = istep+1
         history = history+1
         tauhist(:,:,3) = tauhist(:,:,2)
         tauhist(:,:,2) = tauhist(:,:,1)
         tauhist(:,:,1) = tau(:,:)

         if (ionode) then
           CALL seqopn( iunupdate, 'update', 'FORMATTED', exst )
           WRITE( UNIT = iunupdate, FMT = * ) history
           WRITE( UNIT = iunupdate, FMT = * ) tauhist
           CLOSE( UNIT = iunupdate, STATUS = 'KEEP' )
         endif
                
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
1	conv_ions=.true.
    
  END SUBROUTINE
