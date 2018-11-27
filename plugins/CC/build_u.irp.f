subroutine build_u(u,ijkl_antispinint,iajb_antispinint,abcd_antispinint,t2)
  
  BEGIN_DOC
  ! u
  END_DOC
  
  implicit none
  
  
  integer :: i,j,a,b
  double precision,intent(out) :: u(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)
  double precision,intent(in) :: ijkl_antispinint(n_spin_occ,n_spin_occ,n_spin_occ,n_spin_occ)
  double precision,intent(in) :: iajb_antispinint(n_spin_occ,n_spin_virt,n_spin_occ,n_spin_virt)
  double precision,intent(in) :: abcd_antispinint(n_spin_virt,n_spin_virt,n_spin_virt,n_spin_virt)
  double precision,intent(in) :: t2(n_spin_occ,n_spin_occ,n_spin_virt,n_spin_virt)
  double precision :: pcd, pkl, pkc 
  integer :: k,l,c,d
  u=0d0
  do i=1,n_spin_occ
    do j=1,n_spin_occ
      do a = 1,n_spin_virt
        do b = 1,n_spin_virt
          pcd=0d0
          pkl=0d0
          pkc=0d0
          do c = 1,n_spin_virt
            do d = 1,n_spin_virt
              !cd loop
              pcd += abcd_antispinint(a,b,c,d)*t2(i,j,c,d)
            end do
          end do
          do k = 1,n_spin_occ
            do l= 1 ,n_spin_occ
            !Kl loop
              pkl += ijkl_antispinint(k,l,i,j)*t2(k,l,a,b)
            end do
            do c = 1,n_spin_virt
              !KC only part
              pkc -=iajb_antispinint(k,b,j,c)*t2(i,k,a,c)
              pkc +=iajb_antispinint(k,a,j,c)*t2(i,k,b,c)
              pkc -=iajb_antispinint(k,a,i,c)*t2(j,k,b,c)
              pkc +=iajb_antispinint(k,b,i,c)*t2(j,k,a,c)
            end do
          end do
          u(i,j,a,b) = 0.5d0*pcd+0.5d0*pkl+pkc
        end do
      end do
    end do
  end do
end subroutine
