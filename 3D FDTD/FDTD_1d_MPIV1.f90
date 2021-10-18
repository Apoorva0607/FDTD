! This is a finite-difference time-domain (FDTD) code that solves                                                                               
! Maxwell's equations in 1-D.  A plane wave propagates in free space                                                                            
! along the x-direction.  The unknowns are are Ez and Hy.                                                                                       
! The source is a Gaussian waveform at the left-most position of the                                                                            
! the grid (i = 1).      
Program oneD_FDTD_MPIcode

 implicit none !this ensures the code does not create anynew variable if it is accidently declared
 include 'mpif.h'

 real S,delta,dt,c,mu,eps,nhalf,n0,pi,Caz,Cbz,Day,Dby
 real, allocatable, dimension(:) :: Ez_plot, Ez_plot_idThis
 integer nmax,imax,i,n, ierr, idThis, numP,istart_Ez_store, iend_Ez_store, istart_Ez_update, iend_Ez_update, istart_Hy_store, iend_Hy_store, istart_Hy_update, iend_Hy_update
 integer status(MPI_STATUS_SIZE)
 integer tag1

 real, allocatable, dimension(:):: Ez,Hy

!initialize MPI
 call MPI_INIT(ierr)
!get my rank ( a unique starting with 0) of each processor 
 call MPI_COMM_RANK(MPI_COMM_WORLD,idThis,ierr)
! determine the totoal no. of processors
 call MPI_COMM_SIZE(MPI_COMM_WORLD,numP,ierr)
!pi 
 pi = 3.14159265

! speed of light
                                                                                                                          
 c = 2.99792456E8

 tag1 = 20                                                                                                                                                                                                                                                                    
! permittivity of free space                                                                                                                    
 eps = 8.85418782E-12                                                                                                                           
                                                                                                                                                
! permeability of free space                                                                                                                    
 mu = 4*pi*1E-7                                                                                                                                 
                                                                                                                                                
! space increment                                                                                                                               
 delta = 4.3E-3                                                                                                                                      
                                                                                                                                                
! time step increment (S = 1)                                                                                                                   
S = 0.99                                                                                                                                       
dt = S*delta/c                                                                                                                                    
                                                                                                                                                
! imax is the number of grid cells in the x-direction                                                                                           
imax = 400                                                                                                                                     
                                                                                                                                                
! nmax is the total number of time steps to be run.                                                                                             
nmax = 410

! nhalf is the half-width of the Gaussian Pulse                                                                                                 
nhalf = 20.0                                                                                                                                     
! no is how many time steps the center of the Gaussian is delayed                                                                               
! in time from the start of the simulation                                                                                                      
n0 = nhalf*3

!starting and ending indices for both processors 0 and 1
! Ez and Hy  processor update to account for parellelization
if(idThis == 0) then
  istart_Ez_store = 1
  istart_Ez_update = 2
  iend_Ez_store = imax/2
  iend_Ez_update = imax/2
  istart_Hy_store = 1
  istart_Hy_update = 1
  iend_Hy_store= imax/2
  iend_Hy_update = imax/2 - 1
elseif (idThis == 1) then
  istart_Ez_store = imax/2
  istart_Ez_update = (imax/2) + 1
  iend_Ez_store= imax
  iend_Ez_update = imax - 1
  istart_Hy_store = imax/2
  istart_Hy_update = imax/2
  iend_Hy_store= imax - 1
  iend_Hy_update = imax- 1
endif

! initialize Ez and Hy 

allocate(Ez(istart_Ez_store:iend_Ez_store))
allocate(Hy(istart_Hy_store:iend_Hy_store))
allocate(Ez_plot(1:imax)) ! to hold Ez values across the grid
allocate(Ez_plot_idThis(1:imax)) ! to hold the Ez values on this proc


Ez(:) = 0.0
Hy(:) = 0.0
Ez_plot(:) = 0.0
Ez_plot_idThis(:) = 0.0

!Pre-calculate the updating coefficients                                                                                                       
Caz = 1.0                                                                                                                                      
Cbz = dt/(eps*delta)                                                                                                                              
Day = 1.0                                                                                                                                      
Dby = dt/(mu*delta) 


!**************************************************************************
! Start the time-stepping loop
!**************************************************************************

do n=1,nmax
   
!********************************
! update the Hy fields
!********************************
    do i = istart_Hy_update,iend_Hy_update
        Hy(i) = Day*Hy(i) + Dby*(Ez(i+1) - Ez(i))
    enddo
! to transfer ghost fields from processor 0 to processor 1 for Hy
    if(idThis == 1) then
       call mpi_send(Hy(istart_Hy_update), 1, MPI_REAL, 0, & tag1, MPI_COMM_WORLD, status, ierr)
    endif
! processor 0 recieves Hy value from processor 1
    if(idThis == 0) then
       call mpi_recv(Hy(iend_Hy_store), 1, MPI_REAL, 1, & tag1, MPI_COMM_WORLD, status, ierr)
    endif
    
!********************************    
! Update the Ez fields 
!********************************

    ! implement the source on Ez at i = 1
    Ez(1) = exp(-(REAL(n-n0)/REAL(nhalf))**2)
  
    ! Update the Ez fields 
    do i = istart_Ez_update,iend_Ez_update
        Ez(i)= Caz*Ez(i) + Cbz*(Hy(i) - Hy(i-1))
    enddo
! transfer the ghost fields between processors
! Processor 0 send the Ez componenet after regular update
    if(idThis == 0) then
       call mpi_send(Ez(iend_Ez_update), 1, MPI_REAL, 1, & tag1, MPI_COMM_WORLD, status, ierr)
    endif
! Processor 1 recieve Ez from Proc 0

   if(idThis == 1) then
       call mpi_recv(Ez(istart_Ez_store), 1, MPI_REAL, 0, & tag1, MPI_COMM_WORLD, status, ierr)
    endif
 
!********************************    
! Output
!********************************
   do i = istart_Ez_update, iend_Ez_update
       Ez_plot_idThis(i) = Ez(i)
   enddo

   CALL MPI_REDUCE(Ez_plot_idThis, Ez_plot, imax, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

   if (idThis == 0) then
    
    ! Plot the Ez values at each time step and label the axes.
    open(unit = 20,file ='Ez_plot.dat')
    do i = 1,imax
       write(20,*)Ez_plot(i)
    enddo 
    close(20)
  
   endif


enddo
!**************************************************************************
! End the time-stepping loop
!**************************************************************************
call MPI_FINALIZE(ierr)

stop
end program oneD_FDTD_MPIcode
