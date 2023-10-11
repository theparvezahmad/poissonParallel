program pois3Dcart
!=======================================================================
! Developer        : Parvez Ahmad <parvez.cfd@gmail.com>
!-----------------------------------------------------------------------
! Purpose                : Solves 3D cartesian Poisson Eqn
! Input                        : None
! Output                : Multiple domain-wise solution files
!-----------------------------------------------------------------------
! This program uses FORTRAN 2003 features
! For best viewing, set tab width = 2
! Uses cartesian topology for 3D MPI Parallelization
! Uses SOR to solve system of equations
!=======================================================================
  use mpi
  implicit none
!-----------------------------------------------------------------------
  logical, parameter, dimension(3)  ::period = [.false., .false., .false.]

  integer, parameter, dimension(3)  ::domMPI = [2, 2, 2], &
                                       nP = [51, 51, 48]

  double precision, parameter, dimension(3) ::pStart = [0.0, 0.0, 0.0], &
    pEnd = [1.0, 1.0, 1.0]

  double precision, parameter ::tol = 1e-4, &
    w = 1.75d0
!-----------------------------------------------------------------------
  type gridType
    double precision, allocatable::p(:)
  end type

  type(gridType) ::x(3)
!-----------------------------------------------------------------------
  character(len=30) :: filename
  logical           :: myperiod, bottom, top, east, west, north, south
!-----------------------------------------------------------------------
  integer             ::liney, planeyz, planexz, planexy, STATUS(mpi_status_size)
  integer             ::i, j, k, cnt, nproc, ierr, id, flag, comm3d, sizDP
  integer, dimension(3)    ::bs, es, bn, en, siz, myDim, myCoord, src, dest
  integer, allocatable, dimension(:, :)::tArr, bb, ee
!-----------------------------------------------------------------------
  double precision, dimension(nP(1), nP(2), nP(3))        ::u, f
  double precision  ::t1, t2, t3, res, totRes(2), lhs, rhs, it1, it2, dx(3)

!==============================MPI======================================
  CALL mpi_init(ierr)
  CALL mpi_comm_rank(mpi_comm_world, id, ierr)
  CALL mpi_comm_size(mpi_comm_world, nproc, ierr)

  call mpi_cart_create(mpi_comm_world, 3, domMPI, period, .false., comm3d, ierr)
  call mpi_cart_get(comm3d, 3, myDim, myperiod, mycoord, ierr)

!do i=1,3
!call mpi_cart_shift(comm3d,i-1,1,src(i),dest(i),ierr)
!enddo

  call mpi_cart_shift(comm3d, 0, 1, src(1), dest(1), ierr)
  call mpi_cart_shift(comm3d, 1, 1, src(2), dest(2), ierr)
  call mpi_cart_shift(comm3d, 2, 1, src(3), dest(3), ierr)

  ALLOCATE (tArr(0:maxval(myDim) - 1, 3))
  allocate (bb(3, 0:nproc - 1), ee(3, 0:nproc - 1))

  east = (mycoord(1) == myDim(1) - 1)
  west = (mycoord(1) == 0)
  north = (mycoord(2) == myDim(2) - 1)
  south = (mycoord(2) == 0)
  top = (mycoord(3) == myDim(3) - 1)
  bottom = (mycoord(3) == 0)

  do k = 1, 3

    it1 = nP(k)/myDim(k)
    it2 = myDim(k) - mod(nP(k), myDim(k))

    tArr(0, 1) = 1
    do i = 0, myDim(k) - 1
      if (i == it2) it1 = it1 + 1
      tArr(i, 2) = tArr(i, 1) + it1 - 1
      tArr(i, 3) = tArr(i, 2) - tArr(i, 1) + 1
      if (i == myDim(k) - 1) exit
      tArr(i + 1, 1) = tArr(i, 2) + 1
    end do

    do i = 0, myDim(k) - 1
      if (i == mycoord(k)) then
        bs(k) = tArr(i, 1)
        es(k) = tArr(i, 2)
        siz(k) = tArr(i, 3)
      end if
    end do

    bn(k) = bs(k)
    en(k) = es(k)

    if ((west .and. k == 1) .or. (south .and. k == 2) .or. (bottom .and. k == 3)) bn(k) = bs(k) + 1
    if ((east .and. k == 1) .or. (north .and. k == 2) .or. (top .and. k == 3)) en(k) = es(k) - 1

  end do

  CALL mpi_type_vector(siz(2), 1, nP(1), mpi_double_precision, liney, ierr)
  CALL mpi_type_commit(liney, ierr)
  CALL mpi_type_extent(mpi_double_precision, sizDP, ierr)
  CALL mpi_type_hvector(siz(3), 1, nP(1)*nP(2)*sizDP, liney, planeyz, ierr)
  CALL mpi_type_vector(siz(3), siz(1), nP(1)*nP(2), mpi_double_precision, planexz, ierr)
  CALL mpi_type_vector(siz(2), siz(1), nP(1), mpi_double_precision, planexy, ierr)

  CALL mpi_type_commit(planeyz, ierr)
  CALL mpi_type_commit(planexz, ierr)
  CALL mpi_type_commit(planexy, ierr)

  CALL mpi_barrier(mpi_comm_world, ierr)
!-----------------------------MPI---------------------------------------

  allocate (x(1)%p(nP(1)))
  allocate (x(2)%p(nP(2)))
  allocate (x(3)%p(nP(3)))

!=========================Grid Setup====================================
  x(1)%p(1) = pStart(1)
  dx(1) = (pEnd(1) - pStart(1))/(nP(1) - 1)
  do i = 1, nP(1) - 1
    x(1)%p(i + 1) = x(1)%p(i) + dx(1)
  end do

  x(2)%p(1) = pStart(2)
  dx(2) = (pEnd(2) - pStart(2))/(nP(2) - 1)
  do i = 1, nP(2) - 1
    x(2)%p(i + 1) = x(2)%p(i) + dx(2)
  end do

  x(3)%p(1) = pStart(3)
  dx(3) = (pEnd(3) - pStart(3))/(nP(3) - 1)
  do i = 1, nP(3) - 1
    x(3)%p(i + 1) = x(3)%p(i) + dx(3)
  end do
!-----------------------------------------------------------------------
  f(bs(1):es(1), bs(2):es(2), bs(3):es(3)) = -5.0d0

  u(bs(1):es(1), bs(2):es(2), bs(3):es(3)) = 0.0d0

  t1 = 1.0/dx(1)**2.0
  t2 = 1.0/dx(2)**2.0
  t3 = 1.0/dx(3)**2.0

  flag = 0
  cnt = 0
  do while (.true.)
    cnt = cnt + 1

    do k = bn(3), en(3)
      do j = bn(2), en(2)
        do i = bn(1), en(1)
          u(i, j, k) = (1 - w)*u(i, j, k) + w*((u(i - 1, j, k) + u(i + 1, j, k))*t1 + &
                                               (u(i, j - 1, k) + u(i, j + 1, k))*t2 + &
                                               (u(i, j, k - 1) + u(i, j, k + 1))*t3 - &
                                               f(i, j, k) &
                                               )/(2*(t1 + t2 + t3))
        end do
      end do
    end do
    !=====================================================================
 call mpi_sendrecv(u(en(1),bn(2),bn(3)),1,planeyz,dest(1),50,u(bn(1)-1,bn(2),bn(3)),1,planeyz,src(1) ,50,mpi_comm_world,STATUS,ierr)
 call mpi_sendrecv(u(bn(1),bn(2),bn(3)),1,planeyz,src(1) ,50,u(en(1)+1,bn(2),bn(3)),1,planeyz,dest(1),50,mpi_comm_world,STATUS,ierr)

 call mpi_sendrecv(u(bn(1),en(2),bn(3)),1,planexz,dest(2),50,u(bn(1),bn(2)-1,bn(3)),1,planexz,src(2) ,50,mpi_comm_world,STATUS,ierr)
 call mpi_sendrecv(u(bn(1),bn(2),bn(3)),1,planexz,src(2) ,50,u(bn(1),en(2)+1,bn(3)),1,planexz,dest(2),50,mpi_comm_world,STATUS,ierr)

 call mpi_sendrecv(u(bn(1),bn(2),en(3)),1,planexy,dest(3),50,u(bn(1),bn(2),bn(3)-1),1,planexy,src(3) ,50,mpi_comm_world,STATUS,ierr)
 call mpi_sendrecv(u(bn(1),bn(2),bn(3)),1,planexy,src(3) ,50,u(bn(1),bn(2),en(3)+1),1,planexy,dest(3),50,mpi_comm_world,STATUS,ierr)

    !=======================Residual & Stopping Criterion=================
    res = 0.0d0
    totRes(2) = totRes(1)

    do k = bn(3), en(3)
      do j = bn(2), en(2)
        do i = bn(1), en(1)
          lhs = (u(i - 1, j, k) + u(i + 1, j, k))*t1 + &
                (u(i, j - 1, k) + u(i, j + 1, k))*t2 + &
                (u(i, j, k - 1) + u(i, j, k + 1))*t3 - &
                2*u(i, j, k)*(t1 + t2 + t3)

          rhs = -5.0d0

          res = res + (lhs - rhs)
        end do
      end do
    end do

    CALL MPI_Reduce(res, totRes(1), 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
    totRes(1) = abs(totRes(1)/product(nP))

    if (id == 0) then
      if (totRes(1) < tol .and. totRes(1) < totRes(2)) then
        print *, 'Poisson converged at', cnt, 'iteration : Residual', totRes(1)
        flag = 1
      end if

      if (cnt > 50 .and. totRes(1) > totRes(2)) then
        print *, 'Poisson diverged at', cnt, 'iteration : Residual', totRes(1)
        stop
      end if
    end if

    call mpi_bcast(flag, 1, mpi_integer, 0, mpi_comm_world, ierr)
    if (flag == 1) exit
    !---------------------------------------------------------------------
  end do

!=====================Write component files=============================
  WRITE (filename, '(a,i2.2,a)') '../output/3D', id + 1, '.dat'
  OPEN (UNIT=15, FILE=filename)

  write (15, 100) 'Title = "Poisson Solution"'
  write (15, 100) 'Variables = "x","y","z","u"'
  write (15, 107) 'Zone k=', es(3) - bs(3) + 1, ',j=', es(2) - bs(2) + 1, ',i=', es(1) - bs(1) + 1, ', DATAPACKING="POINT"'

  do k = bs(3), es(3)
    write (15, *)
    do j = bs(2), es(2)
      write (15, *)
      do i = bs(1), es(1)
        write (15, 101) x(1)%p(i), x(2)%p(j), x(3)%p(k), u(i, j, k)
      end do
    end do
  end do
  close (15)
!-----------------------------------------------------------------------
  call mpi_gather(bs(1), 3, mpi_integer, bb(1, 0), 3, mpi_integer, 0, mpi_comm_world, ierr)
  call mpi_gather(es(1), 3, mpi_integer, ee(1, 0), 3, mpi_integer, 0, mpi_comm_world, ierr)
!=======================Write Domain Info===============================
  if (id == 0) then
    open (unit=15, file="domInfo.txt")
    write (15, 108) nP(1), nP(2), nP(3)
    write (15, 109) nproc
    do i = 0, nproc - 1
      write (15, 110) bb(1, i), ee(1, i), bb(2, i), ee(2, i), bb(3, i), ee(3, i)
    end do
    close (15)
  end if
!-----------------------------------------------------------------------
100 format(a)
101 format(4(F8.4, 1X))
107 format(3(a, i3.3), a)
108 format(3(i4))
109 format(i3)
110 format(6(i4, 1x))
!-----------------------------------------------------------------------
  CALL mpi_type_free(planeyz, ierr)
  CALL mpi_type_free(planexz, ierr)
  CALL mpi_type_free(planexy, ierr)
  CALL mpi_type_free(liney, ierr)
  CALL mpi_finalize(ierr)
!-----------------------------------------------------------------------
  if (id == 0) print *, "Program executed successfully"
end program pois3Dcart
