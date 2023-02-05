module potential_LJ
use mc_subroutines
implicit none
contains

    subroutine DV_LJ( m_coord, m_trial, atom, DV, rcutoff, N, hL, L)
    ! Calculates the Lennard-Jones potential energy difference when moving one atom from its initial position to a trial position
    ! Equation (8) of the tccm_distance_unit_manual_2022.pdf
    ! IN : takes the total number of atoms N, the atom that has been changed in the iteration, (3xN) matrix of cordinates before moving the atom 
    !       m_coord, and after moving the atom m_trial, the dimensions of the box (hL and L)
    ! OUT : Lenard Jones potential energy difference DV
        implicit none
        integer*8,intent(in) :: atom, N
        real*8,intent(in) :: m_coord(:,:), m_trial(:,:), rcutoff, hL, L
        integer*8:: i
        real*8 :: rij2_coord, rij2_trial
        real*8,intent(out) :: DV
 
        DV = 0
        do i = 1, N                     ! for all atoms of the system
            if (i .ne. atom) then       ! when the atom is different from the atom being tested

                call dist_sq_with_pbc( m_coord(:,i), m_coord(:,atom), hL, L, rij2_coord )   ! calculate the distance squared between all atoms with 'atom' before moving it
                call dist_sq_with_pbc( m_trial(:,i), m_trial(:,atom), hL, L, rij2_trial )   ! calculate the distance squared between all atoms with 'atom' after moving it

                if (rij2_coord .lt. rcutoff*rcutoff) then       ! if the distance before moving the atom is lower than the cutoff distance
                    if (rij2_trial .lt. rcutoff*rcutoff) then   ! and the distance after moving the atom is lower than the cutoff distance
                        DV = DV + ( ((1/rij2_trial)**(6.) - (1/rij2_trial)**3.) - ((1/rij2_coord)**(6.) - (1/rij2_coord)**3.) ) ! calculate the full DV
                    else                                        ! but if the distance after moving the atom is greater than the cutoff distance
                        DV = DV + ( -((1/rij2_coord)**(6.) - (1/rij2_coord)**3.) )  ! do not consider it for the DV
                    endif
                else                                            ! if the distance before moving the atom is greater than the cutoff distance
                    if (rij2_trial .lt. rcutoff*rcutoff) then   ! and the distance after moving the atom is lower than the cutoff distance
                        DV = DV + ((1/rij2_trial)**(6.) - (1/rij2_trial)**3.)   ! only consider the potential that arises from the trial position
                    else                                        ! but if the distance after moving the atom is also larger than the cutoff distance
                        DV = DV                                 ! the potential energy difference is 0.
                    endif
                endif
            endif
        enddo
        DV = 4 * DV
    end subroutine DV_LJ

end module potential_LJ