module potential_Sti
use mc_subroutines
implicit none
contains

    subroutine DV_Stillinger( m_coord, m_trial, hL, L, rcutoff, N, atom_i, DV )
    ! Calculates the Stillinger potential energy difference considering atoms atom_i, atom_j and atom_k when moving atom_i from
    !   one position to a trial one, taking into consideration equations 19, 20 and 21 of the tccm_distance_unit_manual_2022.pdf
    ! IN : number of atoms of the system N, atom which position is tested in the interation atom_i, dimensions of the box (hL, L),
    !       cutoff distance for the potential rcutoff, (3xN) matrix of initial coordinates, (3xN) matrix of trial coordinates
    ! OUT : value of the Stillinger potential energy difference DV
        implicit none
        integer*8,intent(in) :: N, atom_i
        real*8,intent(in) :: m_coord(:,:), m_trial(:,:), hL, L, rcutoff
        integer*8 :: j, k
        real*8 :: f2val, hijkval, rij2, rik2, rjk2, theta_jik, theta_ijk
        real*8,intent(out) :: DV

        DV = 0
        do j = 1, N
            
            if (j .ne. atom_i) then     ! for all atoms of the system different from the one considered (atom_i)
            ! ---------------------------- INITIAL COORDINATES ---------------------------- 
                call dist_sq_with_pbc( m_coord(:,atom_i), m_coord(:,j), hL, L, rij2 )   ! calculate distance squared between atom_i and atom_j
                if (rij2 .lt. rcutoff*rcutoff) then                                     ! if the distance is lower than the cutoff distance then
                    call f2( rij2, rcutoff, f2val )                                     ! calculate f2 function involving atoms i and j
                    DV = DV - f2val                                                     ! subtract the value from the energy difference
                    do k = j+1, N                                                       ! for all atoms different from j and different from i
                        if (k .ne. atom_i) then
                            call dist_sq_with_pbc( m_coord(:,atom_i), m_coord(:,k), hL, L, rik2 )   ! calculate distance squared between atom_i and atom_k
                            if (rik2 .lt. rcutoff*rcutoff) then                                     ! if the distance is lower than the cutoff distance then
                                call angle_jik( m_coord(:,atom_i), m_coord(:,j), m_coord(:,k), hL, L, theta_jik ) ! calculate angle theta_jik
                                call hijk( rij2, rik2, rcutoff, theta_jik, hijkval )                ! calculate h(rij,rik,theta_jik) function
                                DV = DV - hijkval                                                   ! subtract the value from the energy difference
                            endif
                        endif
                    enddo

                    do k = 1, N         ! for all atoms of the system different from i and j
                        if ( (k .ne. atom_i) .and. (k .ne. j) ) then    
                            call dist_sq_with_pbc( m_coord(:,j), m_coord(:,k), hL, L, rjk2 )    ! calculate distance squared between atom_k and atom_j
                            if (rjk2 .lt. rcutoff*rcutoff) then                                 ! if the distance is lower than the cutoff distance then
                                call angle_jik( m_coord(:,j), m_coord(:,atom_i), m_coord(:,k), hL, L, theta_ijk ) ! calculate angle theta_ijk
                                call hijk( rij2, rjk2, rcutoff, theta_ijk, hijkval )            ! calculate h(rji, rjk, theta_ijk) function
                                DV = DV - hijkval                                               ! subtract the value from the energy difference
                            endif
                        endif
                    enddo
                endif
            ! ---------------------------- TRIAL COORDINATES ---------------------------- 
                ! the same as before but with the trial coordinates for atom_i
                call dist_sq_with_pbc( m_trial(:,atom_i), m_trial(:,j), hL, L, rij2 ) 
                if (rij2 .lt. rcutoff*rcutoff) then
                    call f2( rij2, rcutoff, f2val )
                    DV = DV + f2val                 ! here the f2 function is added to the potential since DV = V_trial - V_initial
                    do k = j+1, N
                        if (k .ne. atom_i) then
                            call dist_sq_with_pbc( m_trial(:,atom_i), m_trial(:,k), hL, L, rik2 ) 
                            if (rik2 .lt. rcutoff*rcutoff) then
                                call angle_jik( m_trial(:,atom_i), m_trial(:,j), m_trial(:,k), hL, L, theta_jik )
                                call hijk( rij2, rik2, rcutoff, theta_jik, hijkval )
                                DV = DV + hijkval   ! h(rij,rik,theta_jik) function is also added
                            endif
                        endif
                    enddo

                    do k = 1, N
                        if ( (k .ne. atom_i) .and. (k .ne. j) ) then
                            call dist_sq_with_pbc( m_trial(:,j), m_trial(:,k), hL, L, rjk2 )
                            if (rjk2 .lt. rcutoff*rcutoff) then
                                call angle_jik( m_trial(:,j), m_trial(:,atom_i), m_trial(:,k), hL, L, theta_ijk)
                                call hijk( rij2, rjk2, rcutoff, theta_ijk, hijkval)
                                DV = DV + hijkval   ! h(rji, rjk, theta_ijk) function is also added
                            endif
                        endif
                    enddo
                endif
            endif
        enddo
    end subroutine DV_Stillinger

    subroutine f2( rij2, rcutoff, f2val)
    ! f2 function to calculate the Stillinger potential, equation (10) of the tccm_distance_unit_manual_2022.pdf
    ! IN : cutoff distance for considering the potential rcutoff, distance squared between two atoms rij2
    ! OUT : value of the f2 function
        implicit none
        real*8, intent(in) :: rcutoff, rij2
        real*8, parameter :: A=7.049556277, B=0.6022245584
        real*8, intent(out):: f2val

        f2val=A* ( (B * rij2**(-2.)) - 1 ) * EXP( 1 / (sqrt(rij2) - rcutoff) )
    end subroutine f2


    subroutine angle_jik( atom_i, atom_j, atom_k, hL, L, theta )
    ! Angle between three atoms atom_i, atom_j, atom_k being atom_i the central one considering perio
    ! IN : (3x1) vector of coordinates of the three atoms, dimensions of the box (hL and L)
    ! OUT : angle between atoms atom_j - atom_i - atom_k theta
        implicit none
        real*8, intent(in):: atom_i(3,1), atom_j(3,1), atom_k(3,1), hL, L
        integer*8 ::i
        real*8 :: vec_ij(3), vec_ik(3)
        real*8, intent(out) :: theta
    
        vec_ij = atom_j(:,1) - atom_i(:,1)  ! vector from atom_i to atom_j
        vec_ik = atom_k(:,1) - atom_i(:,1)  ! vector from atom_i to atom_k
        do i = 1, 3                             ! for each coordinate x,y,z of the vector from atom_i to atom_j, do:
            if (vec_ij(i) .gt. hL) then         ! if the difference between coordinates is greater than the dimensions of the box
                vec_ij(i) = vec_ij(i) - L       ! subtract the dimension of the box L
            elseif (vec_ij(i) .lt. -hL) then    ! if the difference between coordinates is lower than the dimensions of the box
                vec_ij(i) = vec_ij(i) + L       ! add the dimension of the box L
            endif
            if (vec_ik(i) .gt. hL) then         ! the same for the vector from atom_i to atom_k
                vec_ik(i) = vec_ik(i) - L
            elseif(vec_ik(i) .lt. -hL) then
                vec_ik(i) = vec_ik(i) + L
            endif
        enddo
        ! get the value of theta: theta = arc cos( (vec_ij · vec_ik) / (|vec_ij| · |vec_ik|))
        theta = acos(dot_product(vec_ij,vec_ik) / (sqrt(dot_product(vec_ij,vec_ij)) * sqrt(dot_product(vec_ik,vec_ik)) ))
    end subroutine angle_jik
    

    subroutine hijk( rij2, rik2, rcutoff, theta_jik, hijkval )
    ! h function to calculate the Stillinger potential, equation (12) of the tccm_distance_unit_manual_2022.pdf
    ! IN : distance squared between atom_i and atom_j rij2, distance squared between atom_i and atom_k rik2,
    !       cutoff distance for the potential rcutoff, angle between the three atoms, being atom_i the central atom
    ! OUT : value of the h function
        implicit none
        real*8, intent(in):: rij2, rik2, rcutoff, theta_jik
        real*8, parameter :: lambda=21.0, gamma=1.20
        real*8,intent(out) :: hijkval

        hijkval = lambda * EXP((gamma/(sqrt(rij2)-rcutoff)) + &
        (gamma/(sqrt(rik2)-rcutoff))) * ((COS(theta_jik)+(1./3.))*(COS(theta_jik)+(1./3.)))
    end subroutine hijk

end module potential_Sti