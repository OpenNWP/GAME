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

end module definitions













