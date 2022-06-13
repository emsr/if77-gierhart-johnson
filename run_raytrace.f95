
! gfortran -g -fno-align-commons -o run_raytrace run_raytrace.f95

!
! Main test driver.
!
program run_raytrace

  implicit real*8(A-H,O-Z)

  ! qns surface refrac
  !
  !
  !
  !
  common/rytc/qns,qhc,qha,qhs,qqd

  zero = 0.0
  hrp = 0.0 ! ground plane elev (km presumably)
  h1 = 4.8768 / 1000 ! km
  h2 = 0.3048 * 30000 / 1000 ! km
  hp2 = h2 - hrp
  hp1 = h1 - hrp
  ens = 301 ! Maybe! They may translate from ground to platform?

  ! initialize stuff...
  write(6,*) ''
  qlim = -1.56
  qns = 329.0
  qhc = hp1
  qha = hp2
  qhs = hrp
  zero = rytra(zero)
  ry = trcry(qlim)
  ds0 = qqd
  write(6,*) 'ry = ', ry
  write(6,*) 'qqd = ', qqd

  write(6,*) ''
  zer = 0.0
  qns = 301.0
  qhc = 0.0
  qha = hp2
  qhs = hrp
  zero = rytra(zero)
  ry2 = trcry(zer)
  dlsr = qqd
  write(6,*) 'ry2 = ', ry2
  write(6,*) 'qqd = ', qqd

  write(6,*) ''
  qns = ens
  qhc = zer
  qha = hp1
  qhs = hrp
  ry1 = trcry(zer)
  dlst = qqd
  write(6,*) 'ry1 = ', ry1
  write(6,*) 'qqd = ', qqd

  ! Made up....
  write(6,*) ''
  qns = ens
  qhc = hp1
  qha = hp2
  qhs = hrp
  ray = trcry(0.05d0)
  dist = qqd
  write(6,*) 'ray = ', ray
  write(6,*) 'qqd = ', qqd

  ! Made up point downward....
  write(6,*) ''
  qns = ens
  qhc = hp1
  qha = hp2
  qhs = hrp
  ray = trcry(-0.05d0)
  dist = qqd
  write(6,*) 'ray = ', ray
  write(6,*) 'qqd = ', qqd

end program

!
! Ray tracer.
!
real*8 function rytra(theta_init)

  !implicit real*8(a-h,o-z)
  implicit none
  save

  real*8 trcry
  real*8 refractivity_surf, height_tmit, height_recv, elev_surf, distance
  common/rytc/refractivity_surf, height_tmit, height_recv, elev_surf, distance

  real*8 index_of_refrac, index_of_refrac_tmit, index_of_refrac_recv, refractivity_tmit, refractivity_recv
  real*8 refractivity
  real*8 a, a_numer, a_denom, height, theta, radius, radius_tmit, radius_recv
  real*8 radius_earth, radius_surf, delta_refractivity, ce
  real*8 theta_init, theta_all, ry, zero, x, y, z, w, ate, theta_curr, ca, dn, ct, cx
  real*8 refrac_curr, st, teg, radius_curr, n_curr, tnt
  real*8 theta_final
  integer*4 i, im1, nla, nk, nl
  logical*4 recv_above_layers
  dimension a(25), index_of_refrac(25), refractivity(25), height(25), theta(25), radius(25)

  data height/0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.305, 0.5, 0.7, 1.0, 1.524, 2.0, 3.048, 5.0, &
         7.0, 10.0, 20.0, 30.480, 50.0, 70.0, 90.0, 110.0, 225.0, 350.0, 475.0/

  ! Setting up array of refractivity

  delta_refractivity = -7.32 * exp(0.005577 * refractivity_surf)
  ce = log(refractivity_surf / (refractivity_surf + delta_refractivity))
  radius_earth = 6370.
  zero = 0.0
  radius_surf = radius_earth + elev_surf
  do i = 1, 25
    refractivity(i) = exp(-ce * height(i) ) * refractivity_surf * 1.0e-6
    index_of_refrac(i) = 1.0 + refractivity(i)
    radius(i) = radius_earth + height(i) + elev_surf
  end do
  do i = 2, 25
    im1 = i - 1
    a_numer = log(index_of_refrac(i)) - log(index_of_refrac(im1)) ! could use log1p
    a_denom = log(radius(i)) - log(radius(im1))
    a(i) = a_numer / a_denom
  end do
!----------------------------------------------------------
  write(6,*) 'ENS = ', refractivity_surf
  write(6,*) 'DN = ', delta_refractivity
  write(6,*) 'HC = ', height_tmit
  write(6,*) 'HA = ', height_recv
  write(6,*) 'HS = ', elev_surf
  write(6,*) 'az = ', radius_earth
  write(6,*) 'as = ', radius_surf
  do i = 1, 25
    write(6,*) height(i), refractivity(i), index_of_refrac(i), radius(i), a(i)
  end do
!----------------------------------------------------------
  rytra = zero
  theta_init = zero
  return

  ! Entrance for tracing ray

  entry trcry(theta_init)

      theta_curr = theta_init
      radius_tmit = radius_earth + elev_surf + height_tmit
      refractivity_tmit = 1.0e-6 * refractivity_surf * exp(-ce * height_tmit)
      index_of_refrac_tmit = 1.0 + refractivity_tmit
      radius_recv = radius_earth + elev_surf + height_recv
      refractivity_recv = 1.0e-6 * refractivity_surf * exp(-ce * height_recv)
      index_of_refrac_recv = 1.0 + refractivity_recv
      theta_all = 0.0
      ate = theta_curr

      if (theta_curr .ge. 0.0) goto 41

      if (radius(1) .eq. radius_tmit) then
        teg = 0.0
      else
        x = radius(1) / (2.0 * radius_tmit)
        z = (radius_tmit - radius(1)) / radius(1)
        w = (refractivity(1) - refractivity_tmit) / index_of_refrac_tmit
        teg = -2.0 * asin(sqrt(x * (z - w)))
      end if
      if (theta_curr .lt. teg) &
        theta_curr = teg
      ate = abs(theta_curr)

      if (theta_curr .ge. 0.0) goto 41

      do i = 2, 25
        y = 2.0 * (sin(0.5 * ate))**2
        z = (radius(i) - radius_tmit) / radius_tmit
        w = (refractivity_tmit - refractivity(i)) * cos(ate) / index_of_refrac(i)
        x = y + z - w
        if (x .lt. 0.0) cycle
        ct = sqrt(0.5 * radius_tmit * x / radius(i))
        if (ct .le. 1.0) exit
      end do
      ct = 2.0 * asin(ct)
      theta_all = 2.0 * ct * (-a(i) / (a(i) + 1.0))
      theta(i) = ct
      nk = i + 1
      do i = nk, 25
        radius_curr = radius(i)
        n_curr = index_of_refrac(i)
        if (radius_curr - radius_tmit .gt. 0.0) then
          radius_curr = radius_tmit
          n_curr = index_of_refrac_tmit
        end if
        im1 = i - 1
        x = index_of_refrac(im1) * radius(im1) / (n_curr * radius_curr)
        theta(i) = acos(cos(theta(im1)) * x)
        x = 2.0 * (-a(i)) / (a(i) + 1.0)
        theta_all = theta_all + (theta(i) - theta(im1)) * x
        nla = i
        if (radius_curr .eq. radius_tmit) exit
      end do

   40 continue

      if (nla .lt. 2) then
        nla = 2
      end if
      theta(nla - 1) = ate
      recv_above_layers = .true.
      do i = nla, 25
        im1 = i - 1
        radius_curr = radius(i)
        n_curr = index_of_refrac(i)
        refrac_curr = refractivity(i)
        if (radius_curr - radius_recv .gt. 0.0) then
          n_curr = index_of_refrac_recv
          radius_curr = radius_recv
          refrac_curr = refractivity_recv
        end if
        x = radius_tmit / (2.0 * radius_curr)
        y = 2.0 * (sin(0.5 * ate))**2
        z = (radius_curr - radius_tmit) / radius_tmit
        w = (refractivity_tmit - refrac_curr) * cos(ate) / n_curr
        theta(i) = 2.0 * asin(sqrt(x * (y + z - w)))
        x= -a(i) / (a(i)+1.)
        theta_all = theta_all + (theta(i) - theta(im1)) * x
        theta_final = theta(i)
        if (radius(i) .gt. radius_recv) recv_above_layers = .false.
      end do
      if (recv_above_layers) then
        x = index_of_refrac(25) * radius(25) / radius_recv
        theta_final = acos(cos(theta_final) * x)
      end if
      ca = (theta_final - theta_curr + theta_all)
      distance = radius_surf * ca
      dn = distance * 0.5399568034 ! I think they re-used a dn variable.
      ct = cos(theta_all)
      st = sin(theta_all)
      tnt = tan(theta_final)
      y = index_of_refrac_recv / index_of_refrac_tmit
      x = (ct - st * tnt - y) / (y * tan(theta_curr) - st - ct * tnt)
      x = atan(x)
      cx = theta_curr - x
      rytra = theta_final

      return

   41 nla = 25
      do nl = 2, 25
        if (radius_tmit .le. radius(nl)) then
          nla = nl
          exit
        end if
      end do
      goto 40

end function
