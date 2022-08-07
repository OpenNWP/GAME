! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module definitions

  ! This file contains some definitions.
                            
  implicit none
  
  private
  
  public :: wp
  
  ! setting the floating pointeger :: precision
  ! single precision
  integer, parameter :: ps = 6
  integer, parameter :: rs = 37
  
  ! real(wp) :: precision
  integer, parameter :: pd = 12
  integer, parameter :: rd = 37
  
  integer, parameter :: sp = selected_real_kind(ps,rs) ! single precission
  integer, parameter :: dp = selected_real_kind(pd,rd) ! real(wp) :: precission
  
  integer, parameter :: wp = dp                        ! working precission
  
  type t_grid
  
    real(wp), allocatable :: normal_distance(:)
    real(wp), allocatable :: volume(:)
    real(wp), allocatable :: area(:)
    real(wp), allocatable :: z_scalar(:)
    real(wp), allocatable :: z_vector(:)
    real(wp), allocatable :: gravity_potential(:)
    real(wp), allocatable :: gravity_m(:)
    real(wp), allocatable :: slope(:)
    real(wp), allocatable :: theta_v_bg(:)
    real(wp), allocatable :: exner_bg(:)
    real(wp), allocatable :: exner_bg_grad(:)
    real(wp), allocatable :: layer_thickness(:)
    integer,  allocatable :: trsk_indices(:)
    integer,  allocatable :: trsk_modified_curl_indices(:)
    integer,  allocatable :: from_index(:)
    integer,  allocatable :: to_index(:)
    integer,  allocatable :: adjacent_vector_indices_h(:)
    integer,  allocatable :: adjacent_signs_h(:)
    integer,  allocatable :: density_to_rhombi_indices(:)
    real(wp), allocatable :: latitude_scalar(:)
    real(wp), allocatable :: longitude_scalar(:)
    real(wp), allocatable :: inner_product_weights(:)
    real(wp), allocatable :: direction(:)
    real(wp), allocatable :: density_to_rhombi_weights(:)
    real(wp), allocatable :: trsk_weights(:)
    real(wp), allocatable :: sfc_albedo(:)
    real(wp), allocatable :: sfc_rho_c(:)
    real(wp), allocatable :: t_conduc_soil(:)
    real(wp), allocatable :: roughness_length(:)
    integer,  allocatable :: is_land(:)
    integer,  allocatable :: latlon_interpol_indices(:)
    real(wp), allocatable :: latlon_interpol_weights(:)
    real(wp), allocatable :: z_soil_interface(:)
    real(wp), allocatable :: z_soil_center(:)
    real(wp), allocatable :: t_const_soil(:)
    real(wp)              :: mean_velocity_area
    real(wp)              :: z_t_const
    real(wp)              :: toa
    integer               :: oro_id
    real(wp)              :: stretching_parameter
    real(wp)              :: radius
    real(wp)              :: eff_hor_res
    integer               :: no_of_oro_layers
  
  end type t_grid
  
  type t_dual_grid
  
    real(wp), allocatable :: area(:)
    real(wp), allocatable :: z_vector(:)
    real(wp), allocatable :: normal_distance(:)
    integer,  allocatable :: from_index(:)
    integer,  allocatable :: to_index(:)
    integer,  allocatable :: vorticity_indices_triangles(:)
    integer,  allocatable :: vorticity_signs_triangles(:)
    real(wp), allocatable :: f_vec(:)
  
  end type t_dual_grid
  
  type t_state
    
    real(wp), allocatable :: rho(:)
    real(wp), allocatable :: rhotheta_v(:)
    real(wp), allocatable :: theta_v_pert(:)
    real(wp), allocatable :: exner_pert(:)
    real(wp), allocatable :: wind(:)
    real(wp), allocatable :: temperature_soil(:)
  
  end type t_state
  
  type t_diagnostics
    
    real(wp), allocatable :: flux_density(:)
    real(wp), allocatable :: flux_density_div(:)
    real(wp), allocatable :: rel_vort_on_triangles(:)
    real(wp), allocatable :: rel_vort(:)
    real(wp), allocatable :: pot_vort(:)
    real(wp), allocatable :: temperature(:)
    real(wp), allocatable :: c_g_p_field(:)
    real(wp), allocatable :: v_squared(:)
    real(wp), allocatable :: wind_div(:)
    real(wp), allocatable :: curl_of_vorticity(:)
    real(wp), allocatable :: scalar_field_placeholder(:)
    real(wp), allocatable :: vector_field_placeholder(:)
    real(wp), allocatable :: u_at_edge(:)
    real(wp), allocatable :: v_at_edge(:)
    real(wp), allocatable :: u_at_cell(:)
    real(wp), allocatable :: v_at_cell(:)
    real(wp), allocatable :: n_squared(:)
    real(wp), allocatable :: dv_hdz(:)
    real(wp), allocatable :: scalar_flux_resistance(:)
    real(wp), allocatable :: power_flux_density_sensible(:)
    real(wp), allocatable :: power_flux_density_latent(:)
    real(wp), allocatable :: roughness_velocity(:)
    real(wp), allocatable :: monin_obukhov_length(:)
  
  end type t_diagnostics
  
  type t_radiation
    
    real(wp), allocatable :: lat_scal(:)
    real(wp), allocatable :: lon_scal(:)
    real(wp), allocatable :: sfc_sw_in(:)
    real(wp), allocatable :: sfc_lw_out(:)
    real(wp), allocatable :: sfc_albedo(:)
    real(wp), allocatable :: temp_sfc(:)
    real(wp), allocatable :: z_scal(:)
    real(wp), allocatable :: z_vect(:)
    real(wp), allocatable :: rho(:)
    real(wp), allocatable :: temp(:)
    real(wp), allocatable :: rad_tend(:)
  
  end type t_radiation
  
  type t_forcings
    
    real(wp), allocatable :: pgrad_acc_old(:)
    real(wp), allocatable :: pressure_gradient_acc_neg_nl(:)
    real(wp), allocatable :: pressure_gradient_acc_neg_l(:)
    real(wp), allocatable :: pressure_grad_condensates_v(:)
    real(wp), allocatable :: v_squared_grad(:)
    real(wp), allocatable :: pot_vort_tend(:)
    real(wp), allocatable :: sfc_sw_in(:)
    real(wp), allocatable :: sfc_lw_out(:)
    real(wp), allocatable :: radiation_tendency(:)
  
  end type t_forcings

end module definitions













