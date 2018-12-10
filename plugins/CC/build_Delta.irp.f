subroutine build_Delta(oovv_Delta)  

  BEGIN_DOC
  ! Delta
  END_DOC

  implicit none

  double precision,intent(out) :: oovv_Delta(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)

  integer :: i,j,a,b,Conv_Spin_Index

  do i=1,n_spin_occ
    do j=1,n_spin_occ
      do a=n_spin_occ+1,n_spin_tot
        do b=n_spin_occ+1,n_spin_tot
          oovv_Delta(i,j,a,b) = &
                                fock_matrix_mo(Conv_Spin_Index(a),conv_spin_index(a)) &
                              + fock_matrix_mo(conv_spin_index(b),conv_spin_index(b)) &
                              - fock_matrix_mo(conv_spin_index(i),conv_spin_index(i)) &
                              - fock_matrix_mo(conv_spin_index(j),conv_spin_index(j)) &
        enddo
      enddo
    enddo
  enddo

end subroutine build_Delta
