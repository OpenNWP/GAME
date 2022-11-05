! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_definitions

  ! This file contains some definitions.
  
  implicit none
  
  ! setting the floating point precision
  ! single precision
  integer, parameter :: ps = 6                                ! single decimal precision
  integer, parameter :: rs = 37                               ! single exponent precision
  ! double precision
  integer, parameter :: pd = 12                               ! double decimal precision
  integer, parameter :: rd = 37                               ! double exponent precision
  
  integer, parameter :: sp = selected_real_kind(ps,rs)        ! single precision
  integer, parameter :: dp = selected_real_kind(pd,rd)        ! double precision

#ifdef SINGLE_PRECISION
  integer, parameter :: wp = sp                               ! working precision
#endif
#ifndef SINGLE_PRECISION
  integer, parameter :: wp = dp                               ! working precision
#endif

  ! grid quantities
  type t_grid
    
    real(wp), allocatable :: dx(:,:)                          ! horizontal normal distance
    real(wp), allocatable :: dz(:,:)                          ! vertical normal distance
    real(wp), allocatable :: volume(:,:)                      ! volumes of the grid boxes
    real(wp), allocatable :: area_h(:,:)                      ! horizontal areas (areas with horizontal normal)
    real(wp), allocatable :: area_v(:,:)                      ! vertical areas (areas with vertical normal)
    real(wp), allocatable :: z_scalar(:,:)                    ! z-coordinates of the scalar gridpoints
    real(wp), allocatable :: z_vector_h(:,:)                  ! z-coordinates of the horizontal vector points
    real(wp), allocatable :: z_vector_v(:,:)                  ! z-coordinates of the vertical vector points
    real(wp), allocatable :: gravity_potential(:,:)           ! gravity potential
    real(wp), allocatable :: gravity_m_v(:,:)                 ! vertical acceleration due to gravity
    real(wp), allocatable :: slope(:,:)                       ! slope of zeta coordinate surfaces
    real(wp), allocatable :: theta_v_bg(:,:)                  ! background virtual potential temperature
    real(wp), allocatable :: exner_bg(:,:)                    ! background Exner pressure
    real(wp), allocatable :: exner_bg_grad_h(:,:)             ! horizontal gradient of the background Exner pressure
    real(wp), allocatable :: exner_bg_grad_v(:,:)             ! vertical gradient of the background Exner pressure
    real(wp), allocatable :: layer_thickness(:,:)             ! layer thicknesses
    integer,  allocatable :: trsk_indices(:,:)                ! edge indices for computing the TRSK reconstruction
    integer,  allocatable :: trsk_modified_curl_indices(:,:)  ! edge indices for computing the vorticity used in the TRSK reconstruction
    integer,  allocatable :: from_cell(:)                     ! neighbouring cell of an edge in opposite direction of the vector
    integer,  allocatable :: to_cell(:)                       ! neighbouring cell of an edge in opposite direction of the vector
    integer,  allocatable :: adjacent_edges(:,:)              ! neigbouring edges of a cell
    integer,  allocatable :: adjacent_signs(:,:)              ! 1 if neighbouring cell is in the direction of the vecor, -1 otherwise
    integer,  allocatable :: density_to_rhombi_indices(:,:)   ! indices needed for interpolating the densities to the rhombi
    real(wp), allocatable :: lat_c(:)                         ! latitudes of the cell centers
    real(wp), allocatable :: lon_c(:)                         ! longitudes of the cell centers
    real(wp), allocatable :: inner_product_weights(:,:,:)     ! weights needed for computing the inner product
    real(wp), allocatable :: direction(:)                     ! directions of the horizontal vectors
    real(wp), allocatable :: density_to_rhombi_weights(:,:)   ! weights needed for interpolating the densities to the rhombi
    real(wp), allocatable :: trsk_weights(:,:)                ! weights needed for computing the TRSK reconstruction
    real(wp), allocatable :: sfc_albedo(:)                    ! surface albedo
    real(wp), allocatable :: sfc_rho_c(:)                     ! indices needed for interpolating the densities to the rhombi
    real(wp), allocatable :: t_conduc_soil(:)                 ! temperature conductivity of the soil (m**2/s)
    real(wp), allocatable :: roughness_length(:)              ! roughness length at the surface
    integer,  allocatable :: is_land(:)                       ! land-sea mask
    integer,  allocatable :: latlon_interpol_indices(:,:,:)   ! indices for computing the interpolation to the lat-lon grid
    real(wp), allocatable :: latlon_interpol_weights(:,:,:)   ! weights for computing the interpolation to the lat-lon grid
    real(wp), allocatable :: z_soil_interface(:)              ! z-coordinates of the interfaces of the soil layers
    real(wp), allocatable :: z_soil_center(:)                 ! z-coordinates of the centers of the soil layers
    real(wp), allocatable :: t_const_soil(:)                  ! temperature below the lowest soil layer
    real(wp), allocatable :: area_dual_h(:,:)                 ! horizontal dual areas
    real(wp), allocatable :: area_dual_v(:,:)                 ! vertical dual areas
    real(wp), allocatable :: z_vector_dual_h(:,:)             ! z-coordinates of the horizontal dual vectors
    real(wp), allocatable :: z_vector_dual_v(:,:)             ! z-coordinates of the vertical dual vectors
    real(wp), allocatable :: dy(:,:)                          ! tangential grid distances
    real(wp), allocatable :: dz_dual(:,:)                     ! vertical gridpoint distances of the dual grid
    integer,  allocatable :: from_cell_dual(:)                ! neighbouring cell of a dual edge in opposite direction of the vector
    integer,  allocatable :: to_cell_dual(:)                  ! neighbouring cell of a dual edge in opposite direction of the vector
    integer,  allocatable :: vorticity_indices_triangles(:,:) ! indices for computing the triangles at the edges
    integer,  allocatable :: vorticity_signs_triangles(:,:)   ! signs for computing the triangles at the edges
    real(wp), allocatable :: f_vec_h(:)                       ! horizontal Coriolis vector component at the edges
    real(wp), allocatable :: f_vec_v(:)                       ! vertical Coriolis vector component at the edges
    
  end type t_grid
  
  ! state variables (prognostic quantities)
  type t_state
    
    real(wp), allocatable :: rho(:,:,:)            ! mass densities
    real(wp), allocatable :: rhotheta_v(:,:)       ! virtual potential temperature density
    real(wp), allocatable :: theta_v_pert(:,:)     ! virtual potential temperature perturbation
    real(wp), allocatable :: exner_pert(:,:)       ! Exner pressure perturbation
    real(wp), allocatable :: wind_h(:,:)           ! horizontal wind
    real(wp), allocatable :: wind_v(:,:)           ! vertical wind
    real(wp), allocatable :: temperature_soil(:,:) ! temperature in the soil
  
  end type t_state
  
  ! diagnostic quantities
  type t_diag
    
    real(wp), allocatable :: flux_density_h(:,:)                   ! arbitrary horizontal flux density
    real(wp), allocatable :: flux_density_v(:,:)                   ! arbitrary vertical flux density
    real(wp), allocatable :: flux_density_div(:,:)                 ! divergence of the flux density
    real(wp), allocatable :: rel_vort_on_triangles(:,:)            ! relative vorticity on triangles
    real(wp), allocatable :: rel_vort_h(:,:)                       ! horizontal relative vorticity (at edges)
    real(wp), allocatable :: rel_vort_v(:,:)                       ! vertical relative vorticity (at edges)
    real(wp), allocatable :: pot_vort_h(:,:)                       ! horizontal potential vorticity (at edges)
    real(wp), allocatable :: pot_vort_v(:,:)                       ! vertical potential vorticity (at edges)
    real(wp), allocatable :: temperature(:,:)                      ! air temperature
    real(wp), allocatable :: v_squared(:,:)                        ! squared velocity
    real(wp), allocatable :: wind_div(:,:)                         ! divergence of the wind field
    real(wp), allocatable :: curl_of_vorticity_h(:,:)              ! curl of the horizontal vorticity
    real(wp), allocatable :: scalar_placeholder(:,:)               ! placeholder for scalar fields
    real(wp), allocatable :: vector_placeholder_h(:,:)             ! placeholder for horizontal vector fields
    real(wp), allocatable :: vector_placeholder_v(:,:)             ! placeholder for vertical vector fields
    real(wp), allocatable :: n_squared(:,:)                        ! squared Brunt-Väisälä frequency
    real(wp), allocatable :: dv_hdz(:,:)                           ! vertical gradient of horizontal vorticity
    real(wp), allocatable :: scalar_flux_resistance(:)             ! flux resistance at the surface for scalar quantities
    real(wp), allocatable :: power_flux_density_sensible(:)        ! sensible power flux density at the surface
    real(wp), allocatable :: power_flux_density_latent(:)          ! latent power flux density at the surface
    real(wp), allocatable :: roughness_velocity(:)                 ! roughness velocity at the surface
    real(wp), allocatable :: monin_obukhov_length(:)               ! Monin-Obukhov length (m)
    real(wp), allocatable :: temperature_diffusion_heating(:,:)    ! power density resulting from temperature diffusion (W/m**3)
    real(wp), allocatable :: friction_acc_h(:,:)                   ! horizontal friction acceleration
    real(wp), allocatable :: friction_acc_v(:,:)                   ! vertical friction acceleration
    real(wp), allocatable :: heating_diss(:,:)                     ! dissipative heating rate (W/m**3)
    real(wp), allocatable :: molecular_diffusion_coeff(:,:)        ! molecular diffusion coefficient
    real(wp), allocatable :: mass_diffusion_coeff_numerical_h(:,:) ! horizontal effective mass diffusion coefficient
    real(wp), allocatable :: mass_diffusion_coeff_numerical_v(:,:) ! vertical effective mass diffusion coefficient
    real(wp), allocatable :: temp_diffusion_coeff_numerical_h(:,:) ! horizontal effective temperature diffusion coefficient
    real(wp), allocatable :: temp_diffusion_coeff_numerical_v(:,:) ! vertical effective temperature diffusion coefficient
    real(wp), allocatable :: pressure_gradient_decel_factor(:,:)   ! pressure gradient deceleration factor due to condensates
    real(wp), allocatable :: condensates_sediment_heat(:,:)        ! heating rate due to falling condensates
    real(wp), allocatable :: mass_diff_tendency(:,:,:)             ! mass source rate due to mass diffusion
    real(wp), allocatable :: phase_trans_rates(:,:,:)              ! phase transition rates
    real(wp), allocatable :: phase_trans_heating_rate(:,:)         ! phase transition heating rate
    real(wp), allocatable :: viscosity(:,:)                        ! effective viscosity in cell centers
    real(wp), allocatable :: viscosity_rhombi(:,:)                 ! effective viscosity at the edges
    real(wp), allocatable :: viscosity_triangles(:,:)              ! effective viscosity at the triangles
    real(wp), allocatable :: vert_hor_viscosity(:,:)               ! effective visosity for vertical diffusion of horizontal velocity
    real(wp), allocatable :: tke(:,:)                              ! specific turbulent kinetic energy (J/kg)
    real(wp), allocatable :: pgrad_acc_old_h(:,:)                  ! old time step horizontal pressure gradient
    real(wp), allocatable :: pressure_gradient_acc_neg_nl_h(:,:)   ! negative non-linear component of the horizontal pressure gradient
    real(wp), allocatable :: pressure_gradient_acc_neg_nl_v(:,:)   ! vertical non-linear component of the horizontal pressure gradient
    real(wp), allocatable :: pressure_gradient_acc_neg_l_h(:,:)    ! negative linear component of the horizontal pressure gradient
    real(wp), allocatable :: pressure_gradient_acc_neg_l_v(:,:)    ! vertical linear component of the horizontal pressure gradient
    real(wp), allocatable :: pressure_grad_condensates_v(:,:)      ! vertical pressure gradient component created by the condensates
    real(wp), allocatable :: v_squared_grad_h(:,:)                 ! horizontal gradient of squared velocity
    real(wp), allocatable :: v_squared_grad_v(:,:)                 ! vertical gradient of squared velocity
    real(wp), allocatable :: pot_vort_tend_h(:,:)                  ! horizontal component of the vorticity flux term
    real(wp), allocatable :: pot_vort_tend_v(:,:)                  ! vertical component of the vorticity flux term
    real(wp), allocatable :: sfc_sw_in(:)                          ! surface inbound short-wave radiation flux density (W/m**2)
    real(wp), allocatable :: sfc_lw_out(:)                         ! surface outbound long-wave radiation flux density (W/m**2)
    real(wp), allocatable :: radiation_tendency(:,:)               ! radiative flux convergence power density (W/m**3)
    
  end type t_diag

end module mo_definitions













