subroutine Accumulate_X1(v,ijab_antispinint,t2)

  implicit none

  double precision,intent(inout),dimension(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt) :: V
  double precision,intent(in),dimension(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt) :: ijab_antispinint,t2
  double precision,allocatable :: X1(:,:,:,:)
  integer :: i,j,a,b,k,l

  allocate(X1(n_spin_occ,n_spin_occ,n_spin_occ,n_spin_occ))
  call Build_X1(X1,ijab_antispinint,t2)

  do i=1,n_spin_occ
    do j=1,n_spin_occ
      do a=1,n_spin_virt
        do b=1,n_spin_virt
          do k=1,n_spin_occ
            do l=1,n_spin_occ
              V(i,j,a,b) += 0.25d0*X1(k,l,i,j)*t2(k,l,a,b)
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  deallocate(X1)
end subroutine Accumulate_X1


subroutine Accumulate_X2(v,ijab_antispinint,t2)
  
  implicit none
  
  double precision,intent(inout),dimension(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt) :: V
  double precision,intent(in),dimension(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt) :: ijab_antispinint,t2
  double precision,allocatable :: X2(:,:)
  integer :: i,j,a,b,c

  allocate(X2(n_spin_virt,n_spin_virt))
  call Build_X2(X2,ijab_antispinint,t2)

  do i=1,n_spin_occ
    do j=1,n_spin_occ
      do a=1,n_spin_virt
        do b=1,n_spin_virt
          do c=1,n_spin_virt
            V(i,j,a,b) -= 0.5d0*(X2(b,c)*t2(i,j,a,c)+X2(a,c)*t2(i,j,c,b))
          enddo
        enddo
      enddo
    enddo
  enddo

  deallocate(X2)

end subroutine Accumulate_X2

subroutine Accumulate_X3(v,ijab_antispinint,t2)
  
  implicit none
  
  double precision,intent(inout),dimension(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt) :: V
  double precision,intent(in),dimension(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt) :: ijab_antispinint,t2
  double precision,allocatable :: X3(:,:)
  integer :: i,j,a,b,k

  allocate(X3(n_spin_virt,n_spin_virt))
  call Build_X3(X3,ijab_antispinint,t2)

  do i=1,n_spin_occ
    do j=1,n_spin_occ
      do a=1,n_spin_virt
        do b=1,n_spin_virt
          do k=1,n_spin_occ
            V(i,j,a,b) -= 0.5d0*(X3(k,j)*t2(i,k,a,b)+X3(k,i)*t2(k,j,a,b))
          enddo
        enddo
      enddo
    enddo
  enddo

  deallocate(X3)

end subroutine Accumulate_X3


subroutine Accumulate_X4(v,ijab_antispinint,t2)
  
  implicit none
  
  double precision,intent(inout),dimension(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt) :: V
  double precision,intent(in),dimension(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt) :: ijab_antispinint,t2
  double precision,allocatable :: X4(:,:,:,:)
  integer :: i,j,a,b,k,c

  allocate(X4(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt))
  call Build_X4(X4,ijab_antispinint,t2)

  do i=1,n_spin_occ
    do j=1,n_spin_occ
      do a=1,n_spin_virt
        do b=1,n_spin_virt
          do k=1,n_spin_occ
            do c=1,n_spin_virt
              V(i,j,a,b) += X4(i,k,a,c)*t2(j,k,b,c)+X4(i,k,b,c)*t2(k,j,a,c)
            end do
          enddo
        enddo
      enddo
    enddo
  enddo

  deallocate(X4)

end subroutine Accumulate_X3
