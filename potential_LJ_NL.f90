module potential_LJ_NL
use mc_subroutines
implicit none
contains

    subroutine DV_LJ_NL( m_coord, m_trial, m_neighbours, rcutoff, hL, L, atom_i, nmax, DV )
    ! Calculates the Lennard-Jones potential difference between two configurations of the system: m_coord and m_trial using Neighbour Lists
    ! IN : (3xN) matrix of coordinates, (3xN) matrix of trial coordinates, (nmax x N) matrix of neighbour atoms for each atom,
        ! cutoff for the interaction, dimension of the box (hL, L), atom_i that has been changed in the iteration
    ! OUT : Lennard-Jones potential energy difference DV
        implicit none
        integer*8,intent(in) :: atom_i, nmax
        integer,intent(in) :: m_neighbours(:,:)
        real*8,intent(in) :: m_coord(:,:), m_trial(:,:), rcutoff, hL, L
        integer*8 :: j, atom_j
        real*8 :: rij2_coord, rij2_trial
        real*8,intent(out) :: DV

        DV = 0
        do j = 1, nmax
            atom_j = m_neighbours(j,atom_i)     ! atom_j are the neighbour atoms of atom_i (belonging to its neighbour list m_neighbours(:,atom_i))
            if (atom_j .eq. 0) exit             ! if there are no neighbour atoms left, exit the loop

            call dist_sq_with_pbc( m_coord(:,atom_j), m_coord(:,atom_i), hL, L, rij2_coord )    ! calculate distance sq between atom_i and atom_j BEFORE the trial move
            call dist_sq_with_pbc( m_trial(:,atom_j), m_trial(:,atom_i), hL, L, rij2_trial )    ! calculate distance sq between atom_i and atom_j AFTER the trial move

            if (rij2_coord .lt. rcutoff*rcutoff) then           ! if the distance before the trial move is lower than the cutoff distance
                if (rij2_trial .lt. rcutoff*rcutoff) then           ! and the distance after the trial move is also lower than the cutoff distance
                    DV = DV + ( ((1/rij2_trial)**(6.) - (1/rij2_trial)**3.) - ((1/rij2_coord)**(6.) - (1/rij2_coord)**3.) ) ! DV has contributions from before and after the trial move
                else                                                ! but the distance after the trial move is greater than the cutoff distance
                    DV = DV + ( - ((1/rij2_coord)**(6.) - (1/rij2_coord)**3.) ) ! DV only has contribution from before the trial move
                endif
            else                                                ! if the distance before the trial move is greater than the cutoff distance
                if (rij2_trial .lt. rcutoff*rcutoff) then           ! but the distance after the trial move is lower than the cutoff distance
                    DV = DV + ( (1/rij2_trial)**(6.) - (1/rij2_trial)**3. ) ! DV only has contribution from after the trial move
                else                                                ! and the distance after the trial move is also greater than the cutoff distance
                    DV = DV
                endif
            endif
        enddo
        DV = 4 * DV
    end subroutine DV_LJ_NL

end module potential_LJ_NL