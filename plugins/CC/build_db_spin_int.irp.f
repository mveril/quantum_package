subroutine build_oooo_db_spin_int(db_int)
  
  BEGIN_DOC
  ! Build OOOO block of antisymmetrized MO integrals over spinorbitals
  END_DOC
  
  implicit none

! Local  variables

  double precision             :: spinint
  integer                      :: i,j,k,l

! Output variables

  double precision,intent(out) :: db_int(n_spin_occ,n_spin_occ,n_spin_occ,n_spin_occ)

! Initialize variables

  db_int(:,:,:,:) = 0d0

  do i=1,n_spin_occ
    do j=1,n_spin_occ
      do k=1,n_spin_occ
        do l=1,n_spin_occ
           db_int(i,j,k,l) = oooo_spinint(i,j,k,l) - oooo_spinint(i,j,l,k)
        enddo
      enddo 
    enddo
  enddo

end subroutine build_oooo_db_spin_int

subroutine build_oovv_db_spin_int(db_int)

  BEGIN_DOC
  ! Build OOVV block of antisymmetrized MO integrals over spinorbitals
  END_DOC

  implicit none

! Local variables 
  double precision             :: spinint
  integer                      :: i,j,a,b

! Output variables 
  double precision,intent(out) :: db_int(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)

! Initialize variables

  db_int(:,:,:,:) = 0d0

  do i=1,n_spin_occ
    do j=1,n_spin_occ
      do a=1,n_spin_virt
        do b=1,n_spin_virt
          db_int(i,j,a,b) = spinint(i,j,n_spin_occ+a,n_spin_occ+b) & 
                          - spinint(i,j,n_spin_occ+b,n_spin_occ+a)
        enddo
      enddo 
    enddo
  enddo

end subroutine build_oovv_db_spin_int

subroutine build_vvvv_db_spin_int(db_int)

  BEGIN_DOC
  ! Build VVVV block of antisymmetrized MO integrals over spinorbitals
  END_DOC

  implicit none

  double precision,intent(out) :: db_int(n_spin_virt,n_spin_virt,n_spin_virt,n_spin_virt)
  double precision :: spinint
  integer :: a,b,c,d

  db_int(:,:,:,:) = 0d0

  do a=1,n_spin_virt
    do b=1,n_spin_virt
      do c=1,n_spin_virt
        do d=1,n_spin_virt
          db_int(a,b,c,d) = spinint(n_spin_occ+a,n_spin_occ+b,n_spin_occ+c,n_spin_occ+d) &
                          - spinint(n_spin_occ+a,n_spin_occ+b,n_spin_occ+d,n_spin_occ+c)
        enddo
      enddo 
    enddo
  enddo

end subroutine build_vvvv_db_spin_int

subroutine build_ovov_db_spin_int(db_int)

  BEGIN_DOC
  ! iajb integrals
  END_DOC

  implicit none
  
  double precision,intent(out) :: db_int(n_spin_occ,n_spin_virt,n_spin_occ,n_spin_virt)
  double precision :: spinint
  integer :: i,a,j,b

  db_int(:,:,:,:) = 0d0
 
  do i=1,n_spin_occ
    do a=1,n_spin_virt
      do j=1,n_spin_occ
        do b=1,n_spin_virt
          db_int(i,a,j,b) = spinint(i,n_spin_occ+a,j,n_spin_occ+b) &
                          - spinint(i,n_spin_occ+a,n_spin_occ+b,j)
        enddo
      enddo 
    enddo
  enddo

end subroutine build_ovov_db_spin_int

subroutine build_vvoo_db_spin_int(anti)

  BEGIN_DOC
  ! abij integrals
  END_DOC

  implicit none
  
  double precision,intent(out) :: db_int(n_spin_virt,n_spin_virt,n_spin_occ,n_spin_occ)
  double precision :: spinint
  integer :: a,b,i,j

  db_int(:,:,:,:) = 0d0

  do a=1,n_spin_virt
    do b=1,n_spin_virt
      do i=1,n_spin_occ
        do j=1,n_spin_occ
          db_int(a,b,i,j) = spinint(n_spin_occ+a,n_spin_occ+b,i,j) &
                          - spinint(n_spin_occ+a,n_spin_occ+b,j,i)
        enddo
      enddo 
    enddo
  enddo

end subroutine build_vvoo_db_spin_int
