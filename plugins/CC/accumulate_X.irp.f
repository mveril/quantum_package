subroutine accumulate_X1(v,oovv_db_spin_int,t2)

  implicit none

  BEGIN_DOC
! Accumulate X1 on vector V
  END_DOC

  double precision,intent(inout) :: V(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)

  double precision,intent(in)    :: oovv_db_spin_int(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)
  double precision,intent(in)    :: t2(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)

  double precision,allocatable   :: X1(:,:,:,:)

  integer :: i,j,a,b,k,l

  allocate(X1(n_spin_occ,n_spin_occ,n_spin_occ,n_spin_occ))

  call build_X1(X1,oovv_db_spin_int,t2)

  do i=1,n_spin_occ
    do j=1,n_spin_occ
      do a=1,n_spin_virt
        do b=1,n_spin_virt
          do k=1,n_spin_occ
            do l=1,n_spin_occ
              v(i,j,a,b) += 0.25d0*X1(k,l,i,j)*t2(k,l,a,b)
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  deallocate(X1)

end subroutine accumulate_X1


subroutine accumulate_X2(v,oovv_db_spin_int,t2)
  
  implicit none
  
  BEGIN_DOC
!   Accumulate X2 on vector V
  END_DOC

  double precision,intent(inout) :: V(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)

  double precision,intent(in)    :: oovv_db_spin_int(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)
  double precision,intent(in)    :: t2(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)

  double precision,allocatable :: X2(:,:)
  integer :: i,j,a,b,c

  allocate(X2(n_spin_virt,n_spin_virt))

  call Build_X2(X2,oovv_db_spin_int,t2)

  do i=1,n_spin_occ
    do j=1,n_spin_occ
      do a=1,n_spin_virt
        do b=1,n_spin_virt
          do c=1,n_spin_virt
            v(i,j,a,b) = v(i,j,a,b)                &
                       - 0.5d0*X2(b,c)*t2(i,j,a,c) &
                       - 0.5d0*X2(a,c)*t2(i,j,c,b)
          enddo
        enddo
      enddo
    enddo
  enddo

  deallocate(X2)

end subroutine Accumulate_X2

subroutine Accumulate_X3(v,oovv_db_spin_int,t2)
  
  implicit none
  
  BEGIN_DOC
!   Accumulate X3 on vector V
  END_DOC
  
  double precision,intent(inout) :: V(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)

  double precision,intent(in)    :: oovv_db_spin_int(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)
  double precision,intent(in)    :: t2(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)

  double precision,allocatable :: X3(:,:)
  integer :: i,j,a,b,k

  allocate(X3(n_spin_virt,n_spin_virt))
  call Build_X3(X3,oovv_db_spin_int,t2)

  do i=1,n_spin_occ
    do j=1,n_spin_occ
      do a=1,n_spin_virt
        do b=1,n_spin_virt
          do k=1,n_spin_occ
            v(i,j,a,b) = v(i,j,a,b)                &
                       - 0.5d0*X3(k,j)*t2(i,k,a,b) &
                       - 0.5d0*X3(k,i)*t2(k,j,a,b)
          enddo
        enddo
      enddo
    enddo
  enddo

  deallocate(X3)

end subroutine Accumulate_X3


subroutine Accumulate_X4(v,oovv_db_spin_int,t2)
  
  implicit none
  
  BEGIN_DOC
!   Accumulate X4 on vector V
  END_DOC
  
  double precision,intent(inout) :: V(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)

  double precision,intent(in)    :: oovv_db_spin_int(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)
  double precision,intent(in)    :: t2(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)

  double precision,allocatable :: X4(:,:,:,:)

  integer :: i,j,a,b,k,c

  allocate(X4(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt))
  call Build_X4(X4,oovv_db_spin_int,t2)

  do i=1,n_spin_occ
    do j=1,n_spin_occ
      do a=1,n_spin_virt
        do b=1,n_spin_virt
          do k=1,n_spin_occ
            do c=1,n_spin_virt
              V(i,j,a,b) = V(i,j,a,b)              &
                         + X4(i,k,a,c)*t2(j,k,b,c) &
                         + X4(i,k,b,c)*t2(k,j,a,c)
            end do
          enddo
        enddo
      enddo
    enddo
  enddo

  deallocate(X4)

end subroutine Accumulate_X4
