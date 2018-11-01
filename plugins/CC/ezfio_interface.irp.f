! DO NOT MODIFY BY HAND
! Created by $QP_ROOT/scripts/ezfio_interface/ei_handler.py
! from file /mnt/beegfs/mveril/git/quantum_package/src/CC/EZFIO.cfg


BEGIN_PROVIDER [ double precision, thresh_cc  ]
  implicit none
  BEGIN_DOC
! Threshold on the convergence of the coupled cluster.
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_cc_thresh_cc(has)
    if (has) then
      call ezfio_get_cc_thresh_cc(thresh_cc)
    else
      print *, 'cc/thresh_cc not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( thresh_cc, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read thresh_cc with MPI'
    endif
  IRP_ENDIF

  call write_time(6)
  if (mpi_master) then
    write(6, *) 'Read  thresh_cc'
  endif

END_PROVIDER

BEGIN_PROVIDER [ integer, n_it_cc_max  ]
  implicit none
  BEGIN_DOC
! Maximum number of coupled cluster iterations
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_cc_n_it_cc_max(has)
    if (has) then
      call ezfio_get_cc_n_it_cc_max(n_it_cc_max)
    else
      print *, 'cc/n_it_cc_max not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( n_it_cc_max, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read n_it_cc_max with MPI'
    endif
  IRP_ENDIF

  call write_time(6)
  if (mpi_master) then
    write(6, *) 'Read  n_it_cc_max'
  endif

END_PROVIDER
