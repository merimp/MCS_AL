program MCS_AL
! ##########################################################################
!                MONTE CARLO SIMULATION FOR ATOMIC LIQUIDS
! ##########################################################################
! Author:   Meritxell Malagarriga Pérez
! Date:     February 2023
! Subject:  Distance Learning Programming Project
! ##########################################################################


!----------------------- CALLING MODULES USED -----------------------
use mc_subroutines      ! contains subroutines for both potentials (Lennard-Jones and Stillinger) 
use potential_LJ        ! contains subroutines used when considering the Lennard-Jones potential
use potential_Sti       ! contains subroutines used when considering the Stillinger potential
use potential_LJ_NL     ! contains subroutines used when implementing the Lennard-Jones potential with Neighbour List
use potential_Sti_NL    ! contains subroutines used when implementing the Stillinger potential with Neighbour List
!----------------------- END MODULES USED -----------------------


!----------------------- DECLARATION OF VARIABLES -----------------------
implicit none
real :: t_0, t_f                    ! initial and final times for getting the real time of the simulation
integer*8 :: i, j, iteration        ! counters
integer*8 :: N, eq_sweep, prod_sweep, num_iter, atom_i, atom_j, nint_g, Naccum, nmax                ! simulation parameter integer variables
real*8 :: L, hL, T, rho, DV, mean_alpha, Drmax, rcutoff, rmin, rskin, rn, rmax_g, int_g             ! simulation parameter real variables
real*8, allocatable :: m_coord_ini(:,:), m_coord_trial(:,:), rd_g(:,:), m_dist_sq(:,:), m_position_NL(:,:)
integer, allocatable :: m_neighbours(:,:)
integer*8, allocatable :: Nhist(:)
character( len = 2) :: potential, name
character( len = 30) :: structure_file
logical ::  NL, initial_structure
!----------------------- END DECLARATION OF VARIABLES -----------------------

!----------------------- LIST OF FORMATS FOR READING AND WRITING -----------------------
100 FORMAT (t37, i8)    ! for integer
101 FORMAT (t37, f10.7) ! for reals
102 FORMAT (t37, a)     ! for character
103 FORMAT (t37, L )    ! for logical
!----------------------- END LIST OF FORMATS FOR READING AND WRITING -----------------------

call cpu_time( t_0 )    ! Iinitalizing time count



!----------------------- SET UP OF THE SYSTEM -----------------------
! --- Reading input file: 'parameters.in' WHERE THE PARAMETERS ARE GIVEN
! All quantities in reduced units
open ( unit = 10, file = 'parameters.in', action = 'read')
read (10, *)
read (10, 102) potential            ! Potential used: LJ for Lennard-Jones or ST for Stillinger
read (10, 103) NL                   ! Neigbour list (NL): True for using it or False if not
read (10, 103) initial_structure    ! Initial Structure given: True if given or False if not
if ( initial_structure ) then       ! If an initial structure is given
    read (10, 102) structure_file   ! read the name of the file in the 'parameter.in' file
    read(10,*)                      ! skips the number of atoms of the parameters.in file      
    open (11, file = structure_file, action = 'read')
    read (11, *) N                               ! reads the number of atoms of the .xyz file
    read (11, *)
else
    read (10,*)
    read (10, 100) N                ! Number of atoms in the parameters.in file
endif
read (10, 101) rho                  ! Density 
read (10, 101) T                    ! Temperature
read (10, 101) rcutoff              ! Cutoff distance to consider potential interactions
read (10, 100) eq_sweep             ! Sweeps for equilibrating the system
read (10, 100) prod_sweep           ! Sweeps of production
read (10, 101) Drmax                ! Maximum move for MonteCarlo (MC)
read (10, 101) rmin                 ! Minimum atom-atom distance
read (10, 101) rmax_g               ! Maximum distance for the radial distribution g(r) calculation
read (10, 100) nint_g               ! Number of points considered when calculating g(r)
read (10, 101) rskin                ! Extra distance used by the NL
close (10)
! --- Ends reading the parameter file with the parameters of the system ---


! --- Allocating allocatable variables ------------------------------------
allocate( m_coord_ini(3,N), m_coord_trial(3, N), m_dist_sq(N,N), Nhist(nint_g), rd_g(2,nint_g) )

! --- Setting calculated parameters ---------------------------------------
   ! Box Length of the system
    L = ( real(N)/rho )**(1./3.)        ! Box Length
    hL = L/2.                           ! Half of the Box length

   ! Iterations for the Monte Carlo simulation
    eq_sweep = eq_sweep * N             ! Actual number of equilibration sweeps
    prod_sweep = prod_sweep * N         ! Actual number of production sweeps
    num_iter = eq_sweep + prod_sweep    ! Number of MC iterations 
   ! Radial distribution funcition
    int_g = rmax_g / nint_g             ! Interval length for the radial distribution function g(r) calculation
    do i = 1, nint_g                
        rd_g(1, i) = int_g * (i-1)      ! Interval limits for the radial distribution function g(r) calculation
    enddo                               ! rd_g(2,1:nint_g) are the values for the radial distribution function g(ri) = g(1,nint_g)

! --- Initializing Parameters ---------------------------------------
    Naccum = 0                          ! Number of times the histogram of g(r) is updated
    mean_alpha = 0                      ! Mean success rate of the trial moves (€ [0,1]) in a MC simulation. Proportion of atoms that change their position.

!----------------------- END SET UP OF THE SYSTEM -----------------------




!----------------------- INITIAL STRUCTURE OF THE SYSTEM -----------------------
if ( initial_structure ) then ! FOR A GIVEN INITIAL STRUCTURE (in Angstrom): reads the structure of the 'initial_structure' file
    do i = 1, N
        read (11, *) name, m_coord_ini(:, i)                 ! to read the parameter file: Si X Y Z = name x y z
        do j = 1,3                                            
            m_coord_ini(j, i) = m_coord_ini(j, i) / 2.1     ! to convert the coordinates from Angstrom to reduced units
        enddo
    enddo
    close (11)
else    ! FOR A RANDOM INITIAL STRUCTURE: generates the coordinates using the random number fortran function
    call random_initial_structure(N, hL, L, rmin, m_coord_ini)
endif

! WRITE THE INITIAL STRUCTURE AS A .xyz FILE
open ( unit = 12, file = 'initial_structure.xyz', action = 'write' )
write (12, *) N
write (12, *)
do i = 1, N
    write (12, *) m_coord_ini(:, i)
enddo
close (12)

!----------------------- END INITIAL STRUCTURE OF THE SYSTEM -----------------------




!----------------------- SIMULATION CODE WHEN USING NEIGHBOUR LIST -----------------------
if ( NL ) then
    ! NEIGHBOUR LIST EXTRA PARAMETERS SETUP:
    rn = rcutoff + rskin    ! rn distance with rskin some distance such that the atoms occupying the ‘skin’ between the sphere of radius rcutoff and that 
                            ! of radius rn are those which are likely to move closer than rcutoff to the central atom within the next few Monte Carlo moves
    nmax = 2 * int( rho * 4 * 4.D0 * DATAN(1.D0) * rn**(3.) / 3. )  ! maximum number of expected neighbours per atom
    
    allocate( m_position_NL(3, N), m_neighbours(nmax, N) )
    m_position_NL = 0   ! (3xN) matrix with the positions of each atom when its neighbour list was updated for the last time
    m_neighbours = 0    ! (nmax x N) matrix where each column corresponds to the neighbour list of an atom. Initializing to 0


    ! --------- START OF MONTE CARLO ALGORITHM USING NEIGHBOUR LIST ------------
    do iteration = 1, num_iter
        if ( mod(iteration,N) .eq. 0 ) write (*,*) "Monte Carlo iteration", iteration/N," of", num_iter/N

    ! TRIAL MOVE FOR RANDOM ATOM: atom_i. It is accepted or denied within the subroutine
        call trial_move_MC( m_coord_ini, Drmax, L, N, atom_i, m_coord_trial )

    ! (N x N) MATRIX OF DISTANCES SQUARED BETWEEN ALL ATOMS: position ij is the distance squared between atom_i and atom_j
        if ( iteration .eq. 1 ) then
            call matrix_dist_sq( m_coord_ini, N, hL, L, m_dist_sq )                 ! Creates a NxN matrix of distances squared
        else
            call matrix_dist_sq_update( m_coord_ini, atom_i, N, hL , L, m_dist_sq ) ! Updates the distances squared corresponding to atom_i
        endif

    ! NEIGHBOUR LIST: creates the neighbour list of atom_i if it doesn't exist, it is updated if the position of the atom has changed more than rskin/2
        call NL_sub( m_neighbours, m_position_NL, m_coord_ini, m_coord_trial, m_dist_sq, rn, atom_i, N, hL, L, rskin )

    ! POTENTIAL CALCULATION:
        ! * If Lennard-Jones potential is selected:
        if ( potential .eq. "LJ") then    ! it calls the Lennard-Jones potential (using neighbour list) subroutine
            call DV_LJ_NL( m_coord_ini, m_coord_trial, m_neighbours, rcutoff, hL, L, atom_i, nmax, DV)

        ! * If Stillinger potential is selected:
        elseif ( potential .eq. "ST" ) then
            do j = 1, nmax                          ! atom_j's NL must be updated too since the potential uses it to calculate the three-atom term
                atom_j = m_neighbours(j, atom_i)    ! atom_j are the atoms in the neighbour list of atom_i
                call NL_sub( m_neighbours, m_position_NL, m_coord_ini, m_coord_trial, m_dist_sq, rn, atom_j, N, hL, L, rskin )
            enddo
            ! it calls the Stillinger potential (using neighbour list) subroutine
            call DV_Stillinger_NL( m_coord_ini, m_coord_trial, m_neighbours, hL, L, rcutoff, atom_i, nmax, DV )
        endif  
        
        call decision( DV, T, atom_i, m_coord_ini, m_coord_trial, mean_alpha )  ! it decides whether or not to accept the change in structure (trial move)
                                                                                ! based on the value of the energy change that would result from this change

        if ( ( mod(iteration,N) .eq. 0 ) .and. ( iteration .gt. eq_sweep ) ) &  ! every x iterations and during the production sweeps       
        call Nhistogram( m_dist_sq, N, nint_g, int_g, Naccum, Nhist )           ! update the histogram used afterwards for the radial distribution function
    enddo 
    ! ------------ END OF MONTE CARLO ALGORITHM USING NEIGHBOUR LIST ------------
!----------------------- END SIMULATION CODE WHEN USING NEIGHBOUR LIST -----------------------


!----------------------- SIMULATION CODE WHEN ## NOT ##  USING NEIGHBOUR LIST -----------------------
else !No Neigh list
    ! ------------ START MONTE CARLO ALGORITHM ## NOT ## USING NEIGHBOUR LIST ------------
    do iteration = 1, num_iter
        if ( mod(iteration,N) .eq. 0 ) write (*, *) "Monte Carlo iteration", iteration/N," of", num_iter/N
        
    ! TRIAL MOVE FOR RANDOM ATOM: atom_i. It is accepted or denied within the subroutine
        call trial_move_MC( m_coord_ini, Drmax, L, N, atom_i, m_coord_trial )

    ! (N x N) MATRIX OF DISTANCES SQUARED BETWEEN ALL ATOMS: position ij is the distance between atom_i and atom_j
        if ( iteration .eq. 1 ) then
            call matrix_dist_sq( m_coord_ini, N, hL, L, m_dist_sq )                 ! Creates the matrix
        else
            call matrix_dist_sq_update( m_coord_ini, atom_i, N, hL, L, m_dist_sq )  ! Updates the distances corresponding to atom_i
        endif
        
    ! POTENTIAL CALCULATION:
        ! * If Lennard-Jones potential is selected:
        if ( potential .eq. "LJ" ) then        ! it calls the Lennard-Jones potential (NOT using neighbour list) subroutine
            call DV_LJ( m_coord_ini, m_coord_trial, atom_i, DV, rcutoff, N, hL, L )
        ! * If Stillinger potential is selected:
        elseif ( potential .eq. "ST" ) then    ! it calls the Stillinger potential (NOT using neighbour list) subroutine
            call DV_Stillinger( m_coord_ini, m_coord_trial, hL, L, rcutoff, N, atom_i, DV )
        endif  

        call decision( DV, T, atom_i, m_coord_ini, m_coord_trial, mean_alpha )  ! it decide whether or not to accept the change in structure (trial move)
                                                                                ! based on the value of the energy change that would result from this change

        if ( ( mod(iteration,N) .eq. 0 ) .and. ( iteration .gt. eq_sweep ) ) &  ! every x iterations and during the production sweeps
        call Nhistogram( m_dist_sq, N, nint_g, int_g, Naccum, Nhist )           ! update the histogram used afterwards for the radial distribution function
    enddo
    ! ------------ END MONTE CARLO ALGORITHM ## NOT ## USING NEIGHBOUR LIST ------------
endif

write (*, *) "Mean success rate:", mean_alpha/num_iter  ! Success rate of the trial moves for the MC algorithm
!----------------------- SIMULATION CODE WHEN ## NOT ## USING NEIGHBOUR LIST -----------------------



!----------------------- RADIAL DISTRIBUTION FUNCTION -----------------------
call radial_distribution( N, nint_g, rd_g, rho, Naccum, Nhist ) ! calls the radial distribution subroutine

open ( unit = 19, file = 'gr.dat', action = 'write')        ! writes each point of the radial distribution function into a file named 'gr.dat'
do i = 1, nint_g
    write (19, *) rd_g(:, i)
enddo
close (19)

    ! Plotting the radial distribution function with gnuplot:
open (19, file = 'gr_plot.plt')                             ! creates a file named 'radial_plot.plt'
write (19, *)"set nokey"                                        ! erases the key of the plot
write (19, *)"set grid"                                         ! sets grid into the plot
write (19, *)'set xlabel "r"'                                   ! sets the x label to distance
write (19, *)'set ylabel "g(r)"'                                ! sets the y label to the radial distribution function g(r)
write (19, *)"set title 'Radial Distribution function g(r)'"    ! sets the title of the plot
write (19, *)'plot "gr.dat" u 1:2 w l'                      ! plots the values of the radial distribution function g(r) written before in the 'gr.dat' file
close (19)                                                      ! closes the file which will be used by gnuplot
call system('gnuplot -p radial_plot.plt')                       ! calls gnuplot to plot the radial distribution function
!----------------------- END RADIAL DISTRIBUTION FUNCTION -----------------------



!----------------------- CPU TIME CALCULATION -----------------------
call cpu_time( t_f )
write (*, *) "Total CPU time", t_f-t_0, 's'
!----------------------- END CPU TIME CALCULATION -----------------------

end program MCS_AL
