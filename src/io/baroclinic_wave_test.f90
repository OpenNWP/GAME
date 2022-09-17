! This is basically a copy from https://github.com/ClimateGlobalChange/DCMIP2016.

MODULE baroclinic_wave

!=======================================================================
!
!  Date:  July 29, 2015
!
!  Functions for setting up idealized initial conditions for the
!  Ullrich, Melvin, Staniforth and Jablonowski baroclinic instability.
!
!  SUBROUTINE baroclinic_wave_sample(
!    deep,moist,pertt,X,lon,lat,p,z,zcoords,u,v,w,t,phis,ps,rho,q)
!
!  Options:
!     deep    deep atmosphere (1 = yes or 0 = no)
!    moist    include moisture (1 = yes or 0 = no)
!    pertt    type of perturbation (0 = exponential, 1 = stream function)
!        X    Earth scaling factor
!
!  Given a point specified by: 
!      lon    longitude (radians) 
!      lat    latitude (radians) 
!      p/z    pressure (Pa) / height (m)
!  zcoords    1 if z is specified, 0 if p is specified
!
!  the functions will return:
!        p    pressure if z is specified and zcoords = 1 (Pa)
!        u    zonal wind (m s^-1)
!        v    meridional wind (m s^-1)
!        t    temperature (K)
!   thetav    virtual potential temperature (K)
!     phis    surface geopotential (m^2 s^-2)
!       ps    surface pressure (Pa)
!      rho    density (kj m^-3)
!        q    water vapor mixing ratio (kg/kg)
!
!
!  Author: Paul Ullrich
!          University of California, Davis
!          Email: paullrich@ucdavis.edu
!
!=======================================================================

  USE MO_DEFINITIONS, ONLY: wp

  IMPLICIT NONE

!=======================================================================
!    Physical constants
!=======================================================================

  REAL(wp), PARAMETER ::               &
       a     = 6371220.0_wp,           & ! Reference Earth's Radius (m)
       Rd    = 287.0_wp,               & ! Ideal gas const dry air (J kg^-1 K^1)
       g     = 9.80616_wp,             & ! Gravity (m s^2)
       cp    = 1004.5_wp,              & ! Specific heat capacity (J kg^-1 K^1)
       Lvap  = 2.5d6,                 & ! Latent heat of vaporization of water
       Rvap  = 461.5_wp,               & ! Ideal gas constnat for water vapor
       Mvap  = 0.608_wp,               & ! Ratio of molar mass of dry air/water
       pi    = 3.14159265358979_wp,    & ! pi
       p0    = 100000.0_wp,            & ! surface pressure (Pa)
       kappa = 2._wp/7._wp,             & ! Ratio of Rd to cp
       omega = 7.29212d-5,            & ! Reference rotation rate of the Earth (s^-1)
       deg2rad  = pi/180._wp             ! Conversion factor of degrees to radians

!=======================================================================
!    Test case parameters
!=======================================================================
  REAL(wp), PARAMETER ::               &
       T0E        = 310._wp     ,      & ! temperature at equatorial surface (K)
       T0P        = 240._wp     ,      & ! temperature at polar surface (K)
       B          = 2._wp       ,      & ! jet half-width parameter
       K          = 3._wp       ,      & ! jet width parameter
       lapse      = 0.005_wp             ! lapse rate parameter

  REAL(wp), PARAMETER ::               &
       pertu0     = 0.5_wp      ,      & ! SF Perturbation wind velocity (m/s)
       pertr      = 1._wp/6._wp  ,     & ! SF Perturbation radius (Earth radii)
       pertup     = 1.0_wp      ,      & ! Exp. perturbation wind velocity (m/s)
       pertexpr   = 0.1_wp      ,      & ! Exp. perturbation radius (Earth radii)
       pertlon    = pi/9._wp    ,      & ! Perturbation longitude
       pertlat    = 2._wp*pi/9._wp,    & ! Perturbation latitude
       pertz      = 15000._wp   ,      & ! Perturbation height cap
       dxepsilon  = 1.d-5               ! Small value for numerical derivatives
 
  REAL(wp), PARAMETER ::               &
       moistqlat  = 2._wp*pi/9._wp,    & ! Humidity latitudinal width
       moistqp    = 34000._wp,         & ! Humidity vertical pressure width
       moisttr    = 0.1_wp,            & ! Vertical cut-off pressure for humidity
       moistqs    = 1.e-12_wp,         & ! Humidity above cut-off
       moistq0    = 0.018_wp,          & ! Maximum specific humidity
       moistqr    = 0.9_wp,            & ! Maximum saturation ratio
       moisteps   = 0.622_wp,          & ! Ratio of gas constants
       moistT0    = 273.16_wp,         & ! Reference temperature (K)
       moistE0Ast = 610.78_wp            ! Saturation vapor pressure at T0 (Pa) 

CONTAINS

!=======================================================================
!    Generate the baroclinic instability initial conditions
!=======================================================================
  SUBROUTINE baroclinic_wave_test(deep,moist,pertt,X,lon,lat,p,z,zcoords,u,v,t,thetav,phis,ps,rho,q) &
    BIND(c, name = "baroclinic_wave_test")
 
    IMPLICIT NONE

!-----------------------------------------------------------------------
!     input/output params parameters at given location
!-----------------------------------------------------------------------
    INTEGER, INTENT(IN)  :: &
                deep,       & ! Deep (1) or Shallow (0) test case
                moist,      & ! Moist (1) or Dry (0) test case
                pertt         ! Perturbation type

    REAL(wp), INTENT(IN)  :: &
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                X             ! Earth scaling parameter

    REAL(wp), INTENT(INOUT) :: &
                p,            & ! Pressure (Pa)
                z               ! Altitude (m)

    INTEGER, INTENT(IN) :: zcoords     ! 1 if z coordinates are specified
                                       ! 0 if p coordinates are specified

    REAL(wp), INTENT(OUT) :: &
                u,          & ! Zonal wind (m s^-1)
                v,          & ! Meridional wind (m s^-1)
                t,          & ! Temperature (K)
                thetav,     & ! Virtual potential temperature (K)
                phis,       & ! Surface Geopotential (m^2 s^-2)
                ps,         & ! Surface Pressure (Pa)
                rho,        & ! density (kg m^-3)
                q             ! water vapor mixing ratio (kg/kg)

    !------------------------------------------------
    !   Local variables
    !------------------------------------------------
    REAL(wp) :: aref, omegaref
    REAL(wp) :: T0, constH, constC, scaledZ, inttau2, rratio
    REAL(wp) :: inttermU, bigU, rcoslat, omegarcoslat
    REAL(wp) :: eta, qratio, qnum, qden

    !------------------------------------------------
    !   Pressure and temperature
    !------------------------------------------------
    if (zcoords .eq. 1) then
      CALL evaluate_pressure_temperature(deep, X, lon, lat, z, p, t)
    else
      CALL evaluate_z_temperature(deep, X, lon, lat, p, z, t)
    end if

    !------------------------------------------------
    !   Compute test case constants
    !------------------------------------------------
    aref = a / X
    omegaref = omega * X

    T0 = 0.5_wp * (T0E + T0P)

    constH = Rd * T0 / g

    constC = 0.5_wp * (K + 2._wp) * (T0E - T0P) / (T0E * T0P)

    scaledZ = z / (B * constH)

    inttau2 = constC * z * exp(- scaledZ**2)

    ! radius ratio
    if (deep .eq. 0) then
      rratio = 1._wp
    else
      rratio = (z + aref) / aref;
    end if

    !-----------------------------------------------------
    !   Initialize surface pressure
    !-----------------------------------------------------
    ps = p0

    !-----------------------------------------------------
    !   Initialize velocity field
    !-----------------------------------------------------
    inttermU = (rratio * cos(lat))**(K - 1._wp) - (rratio * cos(lat))**(K + 1._wp)
    bigU = g / aref * K * inttau2 * inttermU * t
    if (deep .eq. 0) then
      rcoslat = aref * cos(lat)
    else
      rcoslat = (z + aref) * cos(lat)
    end if

    omegarcoslat = omegaref * rcoslat
    
    u = - omegarcoslat + sqrt(omegarcoslat**2 + rcoslat * bigU)
    v = 0._wp

    !-----------------------------------------------------
    !   Add perturbation to the velocity field
    !-----------------------------------------------------

    ! Exponential type
    if (pertt .eq. 0) then
      u = u + evaluate_exponential(lon, lat, z)

    ! Stream function type
    elseif (pertt .eq. 1) then
      u = u - 1._wp / (2._wp * dxepsilon) *                       &
          ( evaluate_streamfunction(lon, lat + dxepsilon, z)    &
          - evaluate_streamfunction(lon, lat - dxepsilon, z))

      v = v + 1._wp / (2._wp * dxepsilon * cos(lat)) *            &
          ( evaluate_streamfunction(lon + dxepsilon, lat, z)    &
          - evaluate_streamfunction(lon - dxepsilon, lat, z))
    end if

    !-----------------------------------------------------
    !   Initialize surface geopotential
    !-----------------------------------------------------
    phis = 0._wp

    !-----------------------------------------------------
    !   Initialize density
    !-----------------------------------------------------
    rho = p / (Rd * t)

    !-----------------------------------------------------
    !   Initialize specific humidity
    !-----------------------------------------------------
    if (moist .eq. 1) then
      eta = p/p0

      if (eta .gt. moisttr) then
        q = moistq0 * exp(- (lat/moistqlat)**4)          &
                    * exp(- ((eta-1._wp)*p0/moistqp)**2)
      else
        q = moistqs
      end if

      ! Convert virtual temperature to temperature
      t = t / (1._wp + Mvap * q)

    else
      q = 0._wp
    end if

    !-----------------------------------------------------
    !   Initialize virtual potential temperature
    !-----------------------------------------------------
    thetav = t * (1._wp + 0.61_wp * q) * (p0 / p)**(Rd / cp)

  END SUBROUTINE baroclinic_wave_test

!-----------------------------------------------------------------------
!    Calculate pointwise pressure and temperature
!-----------------------------------------------------------------------
  SUBROUTINE evaluate_pressure_temperature(deep, X, lon, lat, z, p, t)

    INTEGER, INTENT(IN)  :: deep ! Deep (1) or Shallow (0) test case

    REAL(wp), INTENT(IN)  :: &
                X,          & ! Earth scaling ratio
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                z             ! Altitude (m)

    REAL(wp), INTENT(OUT) :: &
                p,          & ! Pressure (Pa)
                t             ! Temperature (K)

    REAL(wp) :: aref, omegaref
    REAL(wp) :: T0, constA, constB, constC, constH, scaledZ
    REAL(wp) :: tau1, tau2, inttau1, inttau2
    REAL(wp) :: rratio, inttermT

    !--------------------------------------------
    ! Constants
    !--------------------------------------------
    aref = a / X
    omegaref = omega * X

    T0 = 0.5_wp * (T0E + T0P)
    constA = 1._wp / lapse
    constB = (T0 - T0P) / (T0 * T0P)
    constC = 0.5_wp * (K + 2._wp) * (T0E - T0P) / (T0E * T0P)
    constH = Rd * T0 / g

    scaledZ = z / (B * constH)

    !--------------------------------------------
    !    tau values
    !--------------------------------------------
    tau1 = constA * lapse / T0 * exp(lapse * z / T0) &
         + constB * (1._wp - 2._wp * scaledZ**2) * exp(- scaledZ**2)
    tau2 = constC * (1._wp - 2._wp * scaledZ**2) * exp(- scaledZ**2)

    inttau1 = constA * (exp(lapse * z / T0) - 1._wp) &
            + constB * z * exp(- scaledZ**2)
    inttau2 = constC * z * exp(- scaledZ**2)

    !--------------------------------------------
    !    radius ratio
    !--------------------------------------------
    if (deep .eq. 0) then
      rratio = 1._wp
    else
      rratio = (z + aref) / aref;
    end if

    !--------------------------------------------
    !    interior term on temperature expression
    !--------------------------------------------
    inttermT = (rratio * cos(lat))**K &
             - K / (K + 2._wp) * (rratio * cos(lat))**(K + 2._wp)

    !--------------------------------------------
    !    temperature
    !--------------------------------------------
    t = 1._wp / (rratio**2 * (tau1 - tau2 * inttermT))

    !--------------------------------------------
    !    hydrostatic pressure
    !--------------------------------------------
    p = p0 * exp(- g / Rd * (inttau1 - inttau2 * inttermT))

  END SUBROUTINE evaluate_pressure_temperature

!-----------------------------------------------------------------------
!    Calculate pointwise z and temperature given pressure
!-----------------------------------------------------------------------
  SUBROUTINE evaluate_z_temperature(deep, X, lon, lat, p, z, t)
    
    INTEGER, INTENT(IN)  :: deep ! Deep (1) or Shallow (0) test case

    REAL(wp), INTENT(IN)  :: &
                X,          & ! Earth scaling ratio
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                p             ! Pressure (Pa)

    REAL(wp), INTENT(OUT) :: &
                z,          & ! Altitude (m)
                t             ! Temperature (K)

    INTEGER :: ix

    REAL(wp) :: z0, z1, z2
    REAL(wp) :: p0, p1, p2

    z0 = 0._wp
    z1 = 10000._wp

    CALL evaluate_pressure_temperature(deep, X, lon, lat, z0, p0, t)
    CALL evaluate_pressure_temperature(deep, X, lon, lat, z1, p1, t)

    DO ix = 1, 100
      z2 = z1 - (p1 - p) * (z1 - z0) / (p1 - p0)

      CALL evaluate_pressure_temperature(deep, X, lon, lat, z2, p2, t)

      IF (ABS((p2 - p)/p) .lt. 1.0d-13) THEN
        EXIT
      END IF

      z0 = z1
      p0 = p1

      z1 = z2
      p1 = p2
    END DO

    z = z2

    CALL evaluate_pressure_temperature(deep, X, lon, lat, z, p0, t)

  END SUBROUTINE evaluate_z_temperature

!-----------------------------------------------------------------------
!    Exponential perturbation function
!-----------------------------------------------------------------------
  REAL(wp) FUNCTION evaluate_exponential(lon, lat, z)

    REAL(wp), INTENT(IN)  :: &
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                z             ! Altitude (meters)

    REAL(wp) :: greatcircler, perttaper

    ! Great circle distance
    greatcircler = 1._wp / pertexpr &
      * acos(sin(pertlat) * sin(lat) + cos(pertlat) * cos(lat) * cos(lon - pertlon))

    ! Vertical tapering of stream function
    if (z < pertz) then
      perttaper = 1._wp - 3._wp * z**2 / pertz**2 + 2._wp * z**3 / pertz**3
    else
      perttaper = 0._wp
    end if

    ! Zonal velocity perturbation
    if (greatcircler < 1._wp) then
      evaluate_exponential = pertup * perttaper * exp(- greatcircler**2)
    else
      evaluate_exponential = 0._wp
    end if

  END FUNCTION evaluate_exponential

!-----------------------------------------------------------------------
!    Stream function perturbation function
!-----------------------------------------------------------------------
  REAL(wp) FUNCTION evaluate_streamfunction(lon, lat, z)

    REAL(wp), INTENT(IN)  :: &
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                z             ! Altitude (meters)

    REAL(wp) :: greatcircler, perttaper, cospert

    ! Great circle distance
    greatcircler = 1._wp / pertr &
      * acos(sin(pertlat) * sin(lat) + cos(pertlat) * cos(lat) * cos(lon - pertlon))

    ! Vertical tapering of stream function
    if (z < pertz) then
      perttaper = 1._wp - 3._wp * z**2 / pertz**2 + 2._wp * z**3 / pertz**3
    else
      perttaper = 0._wp
    end if

    ! Horizontal tapering of stream function
    if (greatcircler .lt. 1._wp) then
      cospert = cos(0.5_wp * pi * greatcircler)
    else
      cospert = 0._wp
    end if

    evaluate_streamfunction = &
        (- pertu0 * pertr * perttaper * cospert**4)

  END FUNCTION evaluate_streamfunction

END MODULE baroclinic_wave
