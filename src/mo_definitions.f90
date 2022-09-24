! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_definitions

  ! This file contains some definitions.
  
  implicit none
  
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
    
    real(wp), allocatable :: dx(:,:)
    real(wp), allocatable :: dz(:,:)
    real(wp), allocatable :: volume(:,:)
    real(wp), allocatable :: area_h(:,:)
    real(wp), allocatable :: area_v(:,:)
    real(wp), allocatable :: z_scalar(:,:)
    real(wp), allocatable :: z_vector_h(:,:)
    real(wp), allocatable :: z_vector_v(:,:)
    real(wp), allocatable :: gravity_potential(:,:)
    real(wp), allocatable :: gravity_m_h(:,:)
    real(wp), allocatable :: gravity_m_v(:,:)
    real(wp), allocatable :: slope(:,:)
    real(wp), allocatable :: theta_v_bg(:,:)
    real(wp), allocatable :: exner_bg(:,:)
    real(wp), allocatable :: exner_bg_grad_h(:,:)
    real(wp), allocatable :: exner_bg_grad_v(:,:)
    real(wp), allocatable :: layer_thickness(:,:)
    integer,  allocatable :: trsk_indices(:,:)
    integer,  allocatable :: trsk_modified_curl_indices(:,:)
    integer,  allocatable :: from_cell(:)
    integer,  allocatable :: to_cell(:)
    integer,  allocatable :: adjacent_edges(:,:)
    integer,  allocatable :: adjacent_signs(:,:)
    integer,  allocatable :: density_to_rhombi_indices(:,:)
    real(wp), allocatable :: lat_c(:)
    real(wp), allocatable :: lon_c(:)
    real(wp), allocatable :: inner_product_weights(:,:,:)
    real(wp), allocatable :: direction(:)
    real(wp), allocatable :: density_to_rhombi_weights(:,:)
    real(wp), allocatable :: trsk_weights(:,:)
    real(wp), allocatable :: sfc_albedo(:)
    real(wp), allocatable :: sfc_rho_c(:)
    real(wp), allocatable :: t_conduc_soil(:)
    real(wp), allocatable :: roughness_length(:)
    integer,  allocatable :: is_land(:)
    integer,  allocatable :: latlon_interpol_indices(:,:,:)
    real(wp), allocatable :: latlon_interpol_weights(:,:,:)
    real(wp), allocatable :: z_soil_interface(:)
    real(wp), allocatable :: z_soil_center(:)
    real(wp), allocatable :: t_const_soil(:)
    real(wp), allocatable :: area_dual_h(:,:)
    real(wp), allocatable :: area_dual_v(:,:)
    real(wp), allocatable :: z_vector_dual_h(:,:)
    real(wp), allocatable :: z_vector_dual_v(:,:)
    real(wp), allocatable :: dy(:,:)
    real(wp), allocatable :: dz_dual(:,:)
    integer,  allocatable :: from_cell_dual(:)
    integer,  allocatable :: to_cell_dual(:)
    integer,  allocatable :: vorticity_indices_triangles(:,:)
    integer,  allocatable :: vorticity_signs_triangles(:,:)
    real(wp), allocatable :: f_vec_h(:)
    real(wp), allocatable :: f_vec_v(:)
  
  end type t_grid
  
  type t_state
    
    real(wp), allocatable :: rho(:,:,:)
    real(wp), allocatable :: rhotheta_v(:,:)
    real(wp), allocatable :: theta_v_pert(:,:)
    real(wp), allocatable :: exner_pert(:,:)
    real(wp), allocatable :: wind_h(:,:)
    real(wp), allocatable :: wind_v(:,:)
    real(wp), allocatable :: temperature_soil(:,:)
  
  end type t_state
  
  type t_diag
    
    real(wp), allocatable :: flux_density_h(:,:)
    real(wp), allocatable :: flux_density_v(:,:)
    real(wp), allocatable :: flux_density_div(:,:)
    real(wp), allocatable :: rel_vort_on_triangles(:,:)
    real(wp), allocatable :: rel_vort_h(:,:)
    real(wp), allocatable :: rel_vort_v(:,:)
    real(wp), allocatable :: pot_vort_h(:,:)
    real(wp), allocatable :: pot_vort_v(:,:)
    real(wp), allocatable :: temperature(:,:)
    real(wp), allocatable :: c_g_p_field(:,:)
    real(wp), allocatable :: v_squared(:,:)
    real(wp), allocatable :: wind_div(:,:)
    real(wp), allocatable :: curl_of_vorticity_h(:,:)
    real(wp), allocatable :: curl_of_vorticity_v(:,:)
    real(wp), allocatable :: scalar_placeholder(:,:)
    real(wp), allocatable :: vector_placeholder_h(:,:)
    real(wp), allocatable :: vector_placeholder_v(:,:)
    real(wp), allocatable :: n_squared(:,:)
    real(wp), allocatable :: dv_hdz(:)
    real(wp), allocatable :: scalar_flux_resistance(:)
    real(wp), allocatable :: power_flux_density_sensible(:)
    real(wp), allocatable :: power_flux_density_latent(:)
    real(wp), allocatable :: roughness_velocity(:)
    real(wp), allocatable :: monin_obukhov_length(:)
    real(wp), allocatable :: temperature_diffusion_heating(:,:)
    real(wp), allocatable :: friction_acc_h(:,:)
    real(wp), allocatable :: friction_acc_v(:,:)
    real(wp), allocatable :: heating_diss(:,:)
    real(wp), allocatable :: molecular_diffusion_coeff(:,:)
    real(wp), allocatable :: mass_diffusion_coeff_numerical_h(:,:)
    real(wp), allocatable :: mass_diffusion_coeff_numerical_v(:,:)
    real(wp), allocatable :: temp_diffusion_coeff_numerical_h(:,:)
    real(wp), allocatable :: temp_diffusion_coeff_numerical_v(:,:)
    real(wp), allocatable :: pressure_gradient_decel_factor(:,:)
    real(wp), allocatable :: condensates_sediment_heat(:,:)
    real(wp), allocatable :: mass_diff_tendency(:,:,:)
    real(wp), allocatable :: phase_trans_rates(:,:,:)
    real(wp), allocatable :: phase_trans_heating_rate(:,:)
    real(wp), allocatable :: viscosity(:,:)
    real(wp), allocatable :: viscosity_rhombi(:)
    real(wp), allocatable :: viscosity_triangles(:,:)
    real(wp), allocatable :: vert_hor_viscosity(:)
    real(wp), allocatable :: tke(:,:)
    real(wp), allocatable :: pgrad_acc_old_v(:,:)
    real(wp), allocatable :: pgrad_acc_old_h(:,:)
    real(wp), allocatable :: pressure_gradient_acc_neg_nl_h(:,:)
    real(wp), allocatable :: pressure_gradient_acc_neg_nl_v(:,:)
    real(wp), allocatable :: pressure_gradient_acc_neg_l_h(:,:)
    real(wp), allocatable :: pressure_gradient_acc_neg_l_v(:,:)
    real(wp), allocatable :: pressure_grad_condensates_v(:,:)
    real(wp), allocatable :: v_squared_grad_h(:,:)
    real(wp), allocatable :: v_squared_grad_v(:,:)
    real(wp), allocatable :: pot_vort_tend_h(:,:)
    real(wp), allocatable :: pot_vort_tend_v(:,:)
    real(wp), allocatable :: sfc_sw_in(:)
    real(wp), allocatable :: sfc_lw_out(:)
    real(wp), allocatable :: radiation_tendency(:,:)
  
  end type t_diag
  
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

end module mo_definitions













