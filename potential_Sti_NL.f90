module potential_Sti_NL
use potential_Sti
implicit none
contains

    subroutine DV_Stillinger_NL( m_coord, m_trial, m_neighbours, hL, L, rcutoff, atom_i, n_max, DV )
    ! Calculates the Stillinger potential energy difference considering atoms atom_i, atom_j and atom_k when moving atom_i from
    !   one position to a trial one, taking into consideration equations 19, 20 and 21 of the tccm_distance_unit_manual_2022.pdf
    !   implementing Neighbour Lists
    ! IN : number of atoms of the system N, atom which position is tested in the interation atom_i, dimensions of the box (hL, L),
    !       cutoff distance for the potential rcutoff, (3xN) matrix of initial coordinates, (3xN) matrix of trial coordinates
    ! OUT : value of the Stillinger potential energy difference DV
        implicit none
        integer*8,intent(in) ::  atom_i, n_max
        integer, intent(in) :: m_neighbours(:,:)
        integer*8 ::  j, atom_j, k, atom_k_of_i, atom_k_of_j
        real*8, intent(in) :: m_coord(:,:), m_trial(:,:), hL, L, rcutoff
        real*8 :: f2val, hijkval, rij2, rik2, rjk2, theta_jik, theta_ijk
        real*8, intent(out) :: DV

        DV = 0
        do j = 1, n_max                         ! from 1 to the end of atom_i neighbour list, do
            atom_j = m_neighbours(j,atom_i)     ! atoms of atoms_i's neighbour list 
            if (atom_j .eq. 0) exit             ! if the index in the neighbour list is 0, exit the loop
            ! ---------------------------- INITIAL COORDINATES ---------------------------- 
            call dist_sq_with_pbc( m_coord(:,atom_i), m_coord(:,atom_j), hL, L, rij2 )  ! calculate the distance squared between atom i and atom j
            if (rij2 .lt. rcutoff*rcutoff) then                                         ! if the distance is lower than the cutoff distance then
                call f2( rij2, rcutoff, f2val )                                         ! calculate f2 function involving atoms i and j
                DV = DV - f2val                                                         ! subtract the value from the energy difference
                do k = j+1, n_max                                                                   ! for all atoms with index greater than j in atom_i's neighbour list
                    atom_k_of_i = m_neighbours(k,atom_i)                                            ! atoms of atoms_i's neighbour list 
                    if (atom_k_of_i .eq. 0) exit                                                    ! if the index in the neighbour list is 0, exit the loop
                    call dist_sq_with_pbc( m_coord(:,atom_i), m_coord(:,atom_k_of_i), hL, L, rik2 ) ! calculate distance squared between atom_i and atom_k 
                    if (rik2 .lt. rcutoff*rcutoff) then                                             ! if the distance is lower than the cutoff distance then
                        call angle_jik( m_coord(:,atom_i),m_coord(:,atom_j), &                      ! calculate angle theta_jik
                        m_coord(:,atom_k_of_i), hL, L, theta_jik ) 
                        call hijk( rij2, rik2, rcutoff, theta_jik, hijkval )                        ! calculate h(rij,rik,theta_jik) function involving atoms i, j, k
                        DV = DV - hijkval                                                           ! subtract the value from the energy difference
                    endif
                enddo

                do k = 1, n_max                             ! from 1 to the end of atom_j neighbour list, do
                    atom_k_of_j = m_neighbours(k,atom_j)    ! atoms of atoms_j's neighbour list 
                    if (atom_k_of_j .eq. 0) exit            ! if the index in the neighbour list is 0, exit the loop
                    if (atom_k_of_j .ne. atom_i) then       ! if the atom in atom_j's neighbour list is different from atom_i then
                        call dist_sq_with_pbc( m_coord(:,atom_j),m_coord(:,atom_k_of_j), hL, L, rjk2 )      ! calculate the distance squared between atom j and atom k
                        if (rjk2 .lt. rcutoff*rcutoff) then                                                 ! if the distance is lower than the cutoff distance then
                            call angle_jik( m_coord(:,atom_j),m_coord(:,atom_i), m_coord(:,atom_k_of_j), hL, L, theta_ijk ) ! calculate the angle with atom_j at the center
                            call hijk( rij2, rjk2, rcutoff, theta_ijk, hijkval )                            ! calculate h(rji, rjk, theta_ijk) function
                            DV = DV - hijkval                                                               ! subtract the value from the energy difference
                        endif
                    endif
                enddo
            endif
            ! ---------------------------- TRIAL COORDINATES ---------------------------- 
            ! the same as before but with the trial coordinates for atom_i
            call dist_sq_with_pbc( m_trial(:,atom_i), m_trial(:,atom_j), hL, L, rij2 ) 
            if (rij2 .lt. rcutoff*rcutoff) then
                call f2( rij2,rcutoff,f2val )
                DV = DV + f2val                             ! here the f2 function is added to the potential since DV = V_trial - V_initial
                do k = j+1, n_max
                    atom_k_of_i = m_neighbours(k,atom_i)    ! atoms of atoms_i's neighbour list 
                    if (atom_k_of_i .eq. 0) exit
                    call dist_sq_with_pbc( m_trial(:,atom_i),m_trial(:,atom_k_of_i), hL, L, rik2 ) 
                    if (rik2 .lt. rcutoff*rcutoff) then
                        call angle_jik( m_trial(:,atom_i), m_trial(:,atom_j),&
                            m_trial(:,atom_k_of_i), hL, L, theta_jik ) 
                        call hijk( rij2, rik2, rcutoff, theta_jik, hijkval )
                        DV = DV + hijkval                   ! h(rij,rik,theta_jik) function is also added
                    endif
                enddo

                do k = 1, n_max
                    atom_k_of_j = m_neighbours(k,atom_j)
                    if (atom_k_of_j .eq.0) exit
                    if (atom_k_of_j .ne. atom_i) then
                        call dist_sq_with_pbc( m_trial(:,atom_j), m_trial(:,atom_k_of_j), hL, L, rjk2 )
                        if (rjk2 .lt. rcutoff*rcutoff) then
                            call angle_jik( m_trial(:,atom_j), m_trial(:,atom_i), m_trial(:,atom_k_of_j), hL, L, theta_ijk )
                            call hijk( rij2, rjk2, rcutoff, theta_ijk, hijkval )
                            DV = DV + hijkval               ! h(rji, rjk, theta_ijk) function is also added
                        endif
                    endif
                enddo
            endif
        enddo
    end subroutine DV_Stillinger_NL

end module potential_Sti_NL