module mc_subroutines
implicit none
contains

    subroutine random_initial_structure(N, hL, L, rmin, m_coord_ini)
    ! Atom i is placed randomly and random positions are proposed for the following atoms
    ! Before accepting them, a test is made to check what is the shortest distance to a previously positioned atom.
    ! If this shortest distance falls below some threshold rmin (which is typically of the order of 0.9 in reduced units), 
    ! this proposed position is rejected and a new one is generated.
    ! IN : number of atoms N, dimensions of the box (hL, L), minimum distance at which atoms can be placed rmin
    ! OUT : (3xN) matrix with the initial coordinates of the N atoms of the system
        implicit none
        integer*8, intent(in) :: N
        real*8, intent(in) :: hL, L, rmin
        integer :: reject
        integer*8 :: i, j, k 
        real*8 :: r, r_sq
        real*8, intent(inout) :: m_coord_ini(:,:)

        do i = 1, N
            reject = 1
            do while ( reject .eq. 1 )
                do j = 1, 3
                    call random_number( r )     ! call the function that generates random numbers from 0 to 1 (Fortran intrinsic subroutine)
                    m_coord_ini(j, i) = r * L   ! assign a random number inside the interval [0,L] to atom i
                enddo
                reject = 0
                do k = 1, i-1                   ! for all the atoms with lower index than i
                    ! calculate the distance between atoms i and k
                    ! if it is lower than the rmin value reject the coordinates for atom i
                    if ( reject .eq. 1 ) exit  
                    call dist_sq_with_pbc( m_coord_ini(:, i), m_coord_ini(:, k), hL, L, r_sq )  
                    if ( r_sq**(0.5) .lt. rmin ) reject = 1      ! Position of the new atom is rejected as it does not satisfy the threshold condition
                enddo      
            enddo
        enddo
    end subroutine random_initial_structure


    subroutine trial_move_MC( m_coord_ini, Drmax, L, N, atom, m_coord_trial )
    ! Trial move of the Monte Carlo simulation: MOVING THE RANDOM ATOM A DISTANCE WITHIN THE BOX (-L/2, L/2):
    ! atom is moved, with its x,y and z coordinates being changed by respectively 2*Drmax(w_i-0.5) with w_i being a random number with i=x,y,z
    ! Drmax defines the maximum change that can occur for each coordinate
    ! IN : number of atoms N, (3xN) matrix of coordiantes before the change m_coord_ini, maximum change that can occur to each coordinate Drmax, dimension of the box L
    ! OUT : random atom that will change or not position atom, (3xN) matrix of coodinates after changing the coordinates of atom m_coord_trial
        implicit none
        integer*8,intent(in) :: N
        real*8,intent(in) :: m_coord_ini(:,:), Drmax, L
        integer*8 :: m
        real*8 :: w, m_coord_aux
        integer*8,intent(out) :: atom
        real*8,intent(out) :: m_coord_trial(:,:)

        call random_number(w)       ! takes a random number within the interval [0,1)
        atom = FLOOR(N*w) + 1       ! index of the atom that will change position (or not). takes integer from 1 to the number of atoms N
        m_coord_trial = m_coord_ini ! sets trial matrix as the initial coordinates matrix. Dimensions (3xN)

        do m = 1, 3                 ! for the three coordinates of an atom x,y,z:
            call random_number(w)
            m_coord_aux = m_coord_trial(m,atom) + 2 * Drmax * (w - 0.5) ! moves the coordinates of atom with index atom a distance 2*Drmax*(w-0.5)
        ! IMPLEMENTATION OF PERIODIC BOUNDARY CONDITIONS    
            if (m_coord_aux .gt. L) then                    ! if the new position is larger than the dimension of the box
                m_coord_trial(m,atom) = m_coord_aux - L     ! subtract the dimension of the box from the position
            elseif (m_coord_aux .lt. 0) then                ! if the new position is smaller than 0
                m_coord_trial(m,atom) = m_coord_aux + L     ! add L to it so it is placed inside the box
            else                                            ! if the new position is inside the box, take its coordinates
                m_coord_trial(m,atom) = m_coord_aux
            endif
        enddo
    end subroutine trial_move_MC


    subroutine dist_sq_with_pbc(atom_i, atom_j, hL, L, r_sq)
    ! Calculates the distance squared between two points taking in consideration the periodic boundary conditions (PBC)
    ! IN : takes the position of two atoms (atom_i, atom_j), and the dimensions of the box (hl, L)
    ! OUT : gives the distane squared between them r_sq
        implicit none
        real*8, intent(in) :: atom_i(3,1), atom_j(3,1), hL, L 
        integer*8 :: m
        real*8 :: vec_dist(3)
        real*8,intent(out) :: r_sq                              

        r_sq = 0

        vec_dist = atom_i(:,1) - atom_j(:,1)        ! for the two atoms, save the difference between their coordinates in a vector
        do m = 1,3                                  ! for the x, y, z coordinates, do:
            if (vec_dist(m) .gt. hL) then           ! if the distance is greater than half the length of the box:
                vec_dist(m) = vec_dist(m) - L       ! the distance is with respect to the atom in the previous box
            else if (vec_dist(m) .lt.  -hL) then    ! if the distance is lower than minus half the length of the box:
                vec_dist(m) = vec_dist(m) + L       ! the distance is with respect to the atom in the next box
            endif
            r_sq = r_sq + vec_dist(m)*vec_dist(m)   ! adds x^2 + y^2 + z^2
        enddo
    end subroutine dist_sq_with_pbc


    subroutine matrix_dist_sq( m_coor, N, hL, L, m_dist_sq )
    ! Creates a (NxN) matrix where each element a_ij corresponds to the distance squared between atom_i and atom_j with i, j € [1,N]
    ! IN : takes (3xN) matrix with the coordinates of all the N atoms of the system, the dimensions of the box (hL and L)
    ! OUT : gives (NxN) matrix with the distances squared between all atoms
        implicit none
        integer*8, intent(in) :: N
        real*8, intent(in) :: m_coor(:,:), hL, L 
        integer*8 :: i, j
        real*8 :: rij_sq
        real*8, intent(out) :: m_dist_sq(N,N)

        m_dist_sq = 0
        do i = 1, N
            do j = i+1, N
                call dist_sq_with_pbc( m_coor(:,i), m_coor(:,j), hL, L, rij_sq )    ! calculates distance squared between atom_i and atom_j
                m_dist_sq(i,j) = rij_sq             ! assigns the distance squared to the matrix element ij
                m_dist_sq(j,i) = m_dist_sq(i,j)     ! assigns the distance squared to the matrix element ji
            enddo
        enddo
    end subroutine matrix_dist_sq


    subroutine matrix_dist_sq_update( m_coor, atom , N, hL, L, m_dist_sq )
    ! Same as subroutine matrix_dist_sq, but only modifies the distances squared belonging to a given atom named "atom"
    ! (matrix elements: (a_i,atom) and (atom,i) )
    ! IN : (3xN) matrix with the coordinates of all the N atoms of the system, (NxN) matrix with the distances squared between all atoms, dimension of the box (hL, L)
    ! OUT : (NxN) matrix with the distances squared between all atoms with distances squared 
        implicit none
        integer*8, intent(in) :: atom, N
        real*8, intent(in) :: m_coor(:,:), hL, L
        integer*8 :: i
        real*8 :: rij_sq
        real*8, intent(inout) :: m_dist_sq(N,N)

        do i = 1, N
            call dist_sq_with_pbc( m_coor(:,i), m_coor(:,atom), hL, L, rij_sq )
            m_dist_sq(i,atom) = rij_sq
            m_dist_sq(atom,i) = m_dist_sq(i,atom)
        enddo
    end subroutine matrix_dist_sq_update


    subroutine decision( DV, T, atom, m_coord_ini, m_coord_trial, mean_alpha )
    ! Decides whether or not to accept the change in structure (trial move) based on the value of the energy change that would result from the change
    ! If DV .le. 0: change accepted.
    ! If DV > 0: 
    !   1. The Boltzmann factor alpha = exp( -DV (Dr)/kBT ) is calculated
    !   2. Random number w € (0,1] is chosen
    !       2.a. If alpha .qe. w: change accepted   -->     update structure
    !       2.b. If alpha < w: change rejected      -->     keep current structure
    ! IN : potential energy difference DV, temperature T, atom to which the trial move is applied atom, (3xN) initial coordinates matrix, (3xN) trial coordinates matrix
    ! OUT : (3xN) initial coordinates matrix (updated or not with new coordinates), mean alpha which is the proportion of the atoms that have changed position 
        implicit none
        integer*8,intent(in) :: atom
        real*8, intent(in) :: DV, T, m_coord_trial(:,:)
        real*8 :: alpha, w
        real*8, intent(inout) :: m_coord_ini(:,:), mean_alpha

        if (DV .le. 0.) then                                ! If negative or 0 the new structure is accepted
            m_coord_ini(:,atom) = m_coord_trial(:,atom)
            mean_alpha = mean_alpha + 1
        else                                                ! If positive the structure is accepted or not depending on the Boltzmann factor
            alpha = exp( - DV / T)
            call random_number(w)
            if (alpha .ge. w) then                          ! If Boltzmann factor > random number --> strucuture accepted
                m_coord_ini(:,atom) = m_coord_trial(:,atom)
                mean_alpha = mean_alpha + 1
            else
                continue                                    !If Boltzmann factor > random number --> structure rejected
            endif
        endif
    end subroutine decision


    subroutine NL_sub_atom( m_dist_sq, m_neighbours, rn, atom, N)
    ! Creates the neighbour list column of the atom given. 
    !   1. Sets the column with the index of the atom to 0.
    !   2. Checks for all atoms (except for the one cosidered) if the distance between them d=m_dist_sq(i,atom) is lower than rn.
    !       2.a. If d .le. rn : the index of the other atom is added to the neighbour list of the considered one.
    !       2.b. If d > rn : nothing is done
    ! IN : number of atoms N, atom considered, (NxN) matrix of distances squared, distance of cutoff + rskin = rn
    ! INOUT : (nmax x N) matrix where each column is the neighbour list of the atom with the column index
        implicit none
        integer*8, intent(in) :: N, atom
        real*8, intent(in) :: m_dist_sq(:,:), rn
        integer, intent(inout) :: m_neighbours(:,:)
        integer*8 :: i, counter

        m_neighbours(:,atom) = 0
        counter = 0
        do i = 1, N
            if ((m_dist_sq(i,atom) .le. rn*rn) .and. (i .ne. atom)) then
                counter = counter + 1
                m_neighbours(counter,atom) = int(i,4)
            endif
        enddo
    end subroutine NL_sub_atom


    subroutine NL_sub( m_neighbours, m_position_NL, m_coor, m_trial, m_dist_sq, rn, atom_i, N, hL, L, rskin )
    ! Creates the neighbour list of atom_i if it doesn't exist, it is updated if the position of the atom has changed more than rskin/2
    ! IN : number of atoms N, atom to which the neighbour list is referred to atom_i, (3xN) matrix of initial coordinates m_coor, 
    !       (3xN) matrix of trial coordinates m_trial, NxN matrix of distance squared m_dist_sq, dimensions of the box hL and L, 
    !       distance to which consider the atoms of the neighbour list rn, extra distance considered by rn not belonging to cutoff distance rskin,
    ! INOUT : (nmax x N) matrix of the neighbour list of atom_i m_neighbours, (3xN) matrix with the positions of each atom when its 
    !       neighbour list was updated for the last time m_position_NL.
        implicit none
        integer*8, intent(in) :: N, atom_i
        real*8, intent(in) :: m_coor(:,:), m_trial(:,:), m_dist_sq(:,:), hL, L, rn, rskin
        real*8 :: r_sq
        integer, intent(inout) :: m_neighbours(:,:)
        real*8, intent(inout) :: m_position_NL(:,:)

        if ( m_neighbours(1, atom_i) .eq. 0 ) then                              ! if the neighbour list for atom_i has not been created
            m_position_NL(:, atom_i) = m_coor(:, atom_i)                        ! the position of atom_i is first stored
            call NL_sub_atom( m_dist_sq, m_neighbours, rn, atom_i, N )          ! the neighbour list of atom_i is created
        else                                                                    ! if the neighbout list for atom_i has been created before
            call dist_sq_with_pbc( m_trial(:, atom_i), m_position_NL(:, atom_i), hL, L, r_sq )   ! calculate distance between the current position and the position when the neighbour list was updated
            if ( r_sq**(0.5) .gt. rskin/2. ) then                               ! if the atom has moved a distance larger than rskin/2
                m_position_NL(:, atom_i) = m_coor(:, atom_i)                    ! the new position is stored
                call NL_sub_atom( m_dist_sq, m_neighbours, rn, atom_i, N )      ! the neighbour list is updated
            endif
        endif
    end subroutine NL_sub


    subroutine Nhistogram( m_dist_sq, N, nint_g, delta_r, Naccum, Nhist )
    ! Creates a histogram Ni(ri) that records how many times during the simulation two atoms are separated by a distance in the interval (ri, ri + delta_r],
    ! where delta_r is small enough to generate a smooth g(r) function, but large enough that sampling does not require too many steps.
    ! The histogram is accumulated from r0=0=rd_g(1,0), to rmax=rd_g(1,nint_g), which should be lower than L/2 since the value of g(r) becomes meaningless for r > L/2.
    ! The histogram is updated every sweep of Monte Carlo steps (one sweep meaning N MC iterations, with N the number of atoms)
    ! IN : number of points considered when calculating g(r) nint_g, number of atoms N
    ! INOUT : gives a nint_g dimension vector that stores the histogram for each point ri Nhist, number of times the histogram has been updated Naccum
        implicit none
        integer*8, intent(in) :: nint_g, N
        real*8, intent(in) :: delta_r, m_dist_sq(:,:)
        integer*8 :: i, j, k
        integer*8, intent(inout) :: Nhist(:), Naccum
        
        if ( Naccum .eq. 0 ) Nhist = 0  ! if the histogram has not been created before, initialize Nhist
        do i = 1, nint_g                ! for each point ri of g(ri)
            do j = 1, N                 ! check if the distance of all atoms is between the interval (ri, ri + delta_r], if true, add 1 to the histogram of that point ri
                do k = j+1, N           ! (the j,k loop iterates over the upper diagonal part of the distance squared matrix)
                    if ( (sqrt(m_dist_sq(j,k)) .gt. (delta_r*(i - 1))) .and. (sqrt(m_dist_sq(j,k)) .le. (delta_r*i)) ) then
                        Nhist(i) = Nhist(i) + 1
                    endif
                enddo
            enddo
        enddo
        Naccum = Naccum + 1             ! each time the histogram is updated
    end subroutine Nhistogram


    subroutine radial_distribution( N, nint_g, rd_g, rho, Naccum, Nhist )
    ! Ratio of the density of atoms found at a given distance r from a reference atom. It is computed after the MC simulation has finished.
    ! It compares the number fo atom pairs found to have had distances in each interval (ri, ri + delta_r] with that expected for the given density
    ! where (ri, ri + delta_r] =  (rd_g(1,i), rd_g(1,i) + rd_g(1,2)]
    ! IN : number of atoms N, number of points considered when calculating g(r) nint_g, density of the atoms rho,
    !      number of times that the histogram has been updated Naccum
    ! OUT : interval length for the radial distribution function g(r) calculation rd_g,
        implicit none
        integer*8, intent(in) :: nint_g, N, Nhist(:), Naccum
        real*8, intent(in) :: rho
        integer*8 :: i
        real*8 :: pi, val
        real*8, intent(out) :: rd_g(2,nint_g)

        pi = 4.D0 * DATAN(1.D0)
        val = real(Naccum) * real(N)/2. * 4./3. * pi * rho          ! calculates the denominator of the function
        do i = 1, nint_g                                            ! for all ri=nint_g(1,1:nint_g)
            rd_g(2,i) = real(Nhist(i)) / ( val * ((rd_g(1,i) + rd_g(1,2))**3. - (rd_g(1,i))**3.) ) 
        enddo
    end subroutine radial_distribution

end module mc_subroutines