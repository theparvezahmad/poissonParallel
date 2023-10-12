program poisonSolver_serial
!=======================================================================
! Developer        : Parvez Ahmad <parvez.cfd@gmail.com>
!-----------------------------------------------------------------------
! Purpose                : Solves 3D cartesian Poisson Eqn
! Input                        : None
! Output                : Multiple domain-wise solution files
!-----------------------------------------------------------------------
! This program uses FORTRAN 2003 features
! For best viewing, set tab width = 2
! Uses SOR to solve system of equations
!=======================================================================
  implicit none
!-----------------------------------------------------------------------
  integer, parameter, dimension(3)  ::nP = [51, 51, 48]

  double precision, parameter, dimension(3) ::pStart = [0.0, 0.0, 0.0], &
    pEnd = [1.0, 1.0, 1.0]

  double precision, parameter ::tol = 1e-4, &
    w = 1.75d0
!-----------------------------------------------------------------------
  type gridType
    double precision, allocatable::p(:)
  end type

  type(gridType) ::x(3)
  character(len=30) :: filename
  integer             ::i, j, k, cnt
!-----------------------------------------------------------------------
  double precision, dimension(nP(1), nP(2), nP(3))        ::u, f
  double precision  ::t1, t2, t3, res, totRes(2), lhs, rhs, it1, it2, dx(3)

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
  f = -5.0d0

  u = 0.0d0

  t1 = 1.0/dx(1)**2.0
  t2 = 1.0/dx(2)**2.0
  t3 = 1.0/dx(3)**2.0

  cnt = 0
  do while (.true.)
    cnt = cnt + 1

    do k = 2, nP(3) - 1
      do j = 2, nP(2) - 1
        do i = 2, nP(1) - 1
          u(i, j, k) = (1 - w)*u(i, j, k) + w*((u(i - 1, j, k) + u(i + 1, j, k))*t1 + &
                                               (u(i, j - 1, k) + u(i, j + 1, k))*t2 + &
                                               (u(i, j, k - 1) + u(i, j, k + 1))*t3 - &
                                               f(i, j, k) &
                                               )/(2*(t1 + t2 + t3))
        end do
      end do
    end do
    !=======================Residual & Stopping Criterion=================
    res = 0.0d0
    totRes(2) = totRes(1)

    do k = 2, nP(3) - 1
      do j = 2, nP(2) - 1
        do i = 2, nP(1) - 1
          lhs = (u(i - 1, j, k) + u(i + 1, j, k))*t1 + &
                (u(i, j - 1, k) + u(i, j + 1, k))*t2 + &
                (u(i, j, k - 1) + u(i, j, k + 1))*t3 - &
                2*u(i, j, k)*(t1 + t2 + t3)

          rhs = -5.0d0

          res = res + (lhs - rhs)
        end do
      end do
    end do

    totRes(1) = abs(res/product(nP))

    if (totRes(1) < tol .and. totRes(1) < totRes(2)) then
      print *, 'Poisson converged at', cnt, 'iteration : Residual', totRes(1)
      exit
    end if

    if (cnt > 50 .and. totRes(1) > totRes(2)) then
      print *, 'Poisson diverged at', cnt, 'iteration : Residual', totRes(1)
      stop
    end if
    !---------------------------------------------------------------------
  end do

!=====================Write component files=============================
  WRITE (filename, '(a)') '../output/3D.dat'
  OPEN (UNIT=15, FILE=filename)

  write (15, 100) 'Title = "Poisson Solution"'
  write (15, 100) 'Variables = "x","y","z","u"'
  write (15, 107) 'Zone k=', nP(3), ',j=', nP(2), ',i=', nP(1), ', DATAPACKING="POINT"'

  do k = 1, nP(3)
    write (15, *)
    do j = 1, nP(2)
      write (15, *)
      do i = 1, nP(1)
        write (15, 101) x(1)%p(i), x(2)%p(j), x(3)%p(k), u(i, j, k)
      end do
    end do
  end do
  close (15)
!-----------------------------------------------------------------------
100 format(a)
101 format(4(F8.4, 1X))
107 format(3(a, i3.3), a)
! 108 format(3(i4))
! 109 format(i3)
! 110 format(6(i4, 1x))
!-----------------------------------------------------------------------
  print *, "Program executed successfully"
end program poisonSolver_serial
