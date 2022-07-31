! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module vertical_grid

  ! This file contains functions that compute properties of the vertical grid.

  use iso_c_binding
  use definitions, only: wp
  use constants,   only: gravity,surface_temp,tropo_height,lapse_rate,inv_height,t_grad_inv,r_d, &
                         p_0_standard,c_d_p,p_0
  use grid_nml,    only: n_scalars_h,n_scalars,n_vectors_per_layer,n_layers,n_levels, &
                         n_vectors,n_dual_scalars_h,n_dual_scalars,n_vectors_h, &
                         n_dual_vectors,n_dual_vectors_per_layer,toa,n_oro_layers,stretching_parameter
  use geodesy,     only: calculate_vertical_area
  
  implicit none
  
  private
  
  public :: set_z_scalar
  public :: set_gravity_potential
  public :: set_volume
  public :: standard_temp
  public :: standard_pres
  public :: set_z_scalar_dual
  public :: set_background_state
  public :: set_area
  
  contains
  
  subroutine set_z_scalar(z_scalar,oro,max_oro) &
  bind(c,name = "set_z_scalar")

    ! This function sets the z coordinates of the scalar data points.
    
    real(wp), intent(out) :: z_scalar(n_scalars)
    real(wp), intent(in)  :: oro(n_scalars_h)
    real(wp), intent(in)  :: max_oro
    
    ! local variables
    integer(c_int) :: h_index,level_index,layer_index
    real(wp)       :: A,B,sigma_z,z_rel,z_vertical_vector_pre(n_levels)
    
    ! the heights are defined according to z_k = A_k + B_k*oro with A_0 = toa, A_{N_LEVELS} = 0, B_0 = 0, B_{N_LEVELS} = 1
    
    ! loop over all columns
    !$omp parallel do private(h_index,level_index,layer_index,A,B,sigma_z,z_rel,z_vertical_vector_pre)
    do h_index=1,n_scalars_h
    
      ! filling up z_vertical_vector_pre
      do level_index=1,n_levels
        z_rel = 1._wp-(level_index-1._wp)/n_layers ! z/toa
        sigma_z = z_rel**stretching_parameter
        A = sigma_z*toa ! the height without orography
        ! B corrects for orography
        if (level_index>=n_layers-n_oro_layers+1) then
          B = (level_index-(n_layers-n_oro_layers)-1._wp)/n_oro_layers
        else
          B = 0._wp
        endif
        z_vertical_vector_pre(level_index) = A + B*oro(h_index)
      enddo
    
      ! doing a check
      if (h_index==1) then
        if (max_oro>=z_vertical_vector_pre(n_layers-n_oro_layers+1)) then
          write(*,*) "Maximum of orography larger or equal to the height of the lowest flat level."
          call exit(1)
        endif
      endif
    
      ! placing the scalar points in the middle between the preliminary values of the adjacent levels
      do layer_index=0,n_layers-1
        z_scalar(layer_index*n_scalars_h+h_index) = 0.5_wp*( &
        z_vertical_vector_pre(layer_index+1)+z_vertical_vector_pre(layer_index+2))
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine set_z_scalar
  
  subroutine set_gravity_potential(z_scalar,gravity_potential,radius) &
  bind(c,name = "set_gravity_potential")

    ! This subroutine computes the gravity potential.
    real(wp), intent(in)  :: z_scalar(n_scalars)
    real(wp), intent(out) :: gravity_potential(n_scalars)
    real(wp), intent(in)  :: radius
  
    ! local variables
    integer(c_int) :: ji
    
    !$omp parallel do private(ji)
    do ji=1,n_scalars
      gravity_potential(ji) = -gravity*(radius**2/(radius+z_scalar(ji))-radius)
    enddo
    !$omp end parallel do
    
  end subroutine set_gravity_potential
  
  
  subroutine set_volume(volume,z_vector,area,radius) &
  bind(c,name = "set_volume")

    ! This subroutine computes the volumes of the grid boxes.
  
    real(wp), intent(out) :: volume(n_scalars)
    real(wp), intent(in)  :: z_vector(n_vectors)
    real(wp), intent(in)  :: area(n_vectors)
    real(wp), intent(in)  :: radius
  
    ! local variables
    integer(c_int) :: ji,layer_index,h_index
    real(wp) :: radius_1,radius_2, base_area
    
    !$omp parallel do private(ji,layer_index,h_index,radius_1,radius_2,base_area)
    do ji=1,n_scalars
      layer_index = (ji-1)/n_scalars_h
      h_index = ji-layer_index*n_scalars_h
      base_area = area(h_index+(layer_index+1)*n_vectors_per_layer)
      radius_1 = radius+z_vector(h_index+(layer_index+1)*n_vectors_per_layer)
      radius_2 = radius+z_vector(h_index+layer_index*n_vectors_per_layer)
      volume(ji) = base_area/(3._wp*radius_1**2)*(radius_2**3-radius_1**3)
    enddo
    !$omp end parallel do
    
  end subroutine set_volume
  
  function standard_temp(z_height) &
  bind(c,name = "standard_temp")
  
    ! This function returns the temperature in the standard atmosphere.
    
    real(wp), intent(in) :: z_height
    real(wp)             :: standard_temp
    
    ! local variables
    real(wp) :: tropo_temp_standard
    
    tropo_temp_standard = surface_temp-tropo_height*lapse_rate
    
    if (z_height<tropo_height) then
      standard_temp = surface_temp-z_height*lapse_rate
    elseif (z_height<inv_height) then
      standard_temp = tropo_temp_standard
    else
      standard_temp = tropo_temp_standard+t_grad_inv*(z_height-inv_height)
    endif
    
  end function standard_temp

  function standard_pres(z_height) &
  bind(c,name = "standard_pres")
  
    ! This function returns the pressure in the standard atmosphere.
    real(wp), intent(in) :: z_height
    real(wp)             :: standard_pres
 
    ! local variables
    real(wp) :: tropo_temp_standard,pressure_at_inv_standard
    
    tropo_temp_standard = surface_temp-tropo_height*lapse_rate
    
    if (z_height<tropo_height) then
      standard_pres = p_0_standard*(1._wp-lapse_rate*z_height/surface_temp)**(gravity/(r_d*lapse_rate))
    else if (z_height<inv_height) then
      standard_pres = p_0_standard*(1._wp-lapse_rate*tropo_height/surface_temp) &
      **(gravity/(r_d*lapse_rate)) &
      *exp(-gravity*(z_height-tropo_height)/(r_d*tropo_temp_standard))
    else
      pressure_at_inv_standard = p_0_standard*(1._wp-lapse_rate*tropo_height/surface_temp) &
      **(gravity/(r_d*lapse_rate))*exp(-gravity*(inv_height-tropo_height)/(r_d*tropo_temp_standard))
      standard_pres = pressure_at_inv_standard*(1._wp-lapse_rate*(z_height-inv_height)/surface_temp) &
      **(gravity/(r_d*lapse_rate))
    endif
    
  end function standard_pres
  
  subroutine set_z_scalar_dual(z_scalar_dual,z_vector,from_index,to_index,vorticity_indices_triangles) &
  bind(c,name = "set_z_scalar_dual")
  
    ! This subroutine sets the z coordinates of the dual scalar points.
  
    real(wp), intent(out)      :: z_scalar_dual(n_dual_scalars)
    real(wp), intent(in)       :: z_vector(n_vectors)
    integer(c_int), intent(in) :: from_index(n_vectors_h),to_index(n_vectors_h), &
                                  vorticity_indices_triangles(3*n_dual_scalars_h)

    ! local variables
    integer(c_int) :: ji,layer_index,h_index
  
    !$omp parallel do private(ji,layer_index,h_index)
    do ji=1,n_dual_scalars
      layer_index = (ji-1)/n_dual_scalars_h
      h_index = ji-layer_index*n_dual_scalars_h
      z_scalar_dual(ji) &
      = 1._wp/6._wp*( &
      z_vector(layer_index*n_vectors_per_layer+1+from_index(1+vorticity_indices_triangles(3*(h_index-1)+1))) &
      + z_vector(layer_index*n_vectors_per_layer+1+from_index(1+vorticity_indices_triangles(3*(h_index-1)+2))) &
      + z_vector(layer_index*n_vectors_per_layer+1+from_index(1+vorticity_indices_triangles(3*(h_index-1)+3))) &
      + z_vector(layer_index*n_vectors_per_layer+1+to_index(1+vorticity_indices_triangles(3*(h_index-1)+1))) &
      + z_vector(layer_index*n_vectors_per_layer+1+to_index(1+vorticity_indices_triangles(3*(h_index-1)+2))) &
      + z_vector(layer_index*n_vectors_per_layer+1+to_index(1+vorticity_indices_triangles(3*(h_index-1)+3))))
    enddo
    !$omp end parallel do
    
  end subroutine set_z_scalar_dual
  
  subroutine set_background_state(z_scalar,gravity_potential,theta_v_bg,exner_bg) &
  bind(c,name = "set_background_state")

    ! This subroutine sets the hydrostatic background state.
    
    real(wp), intent(in)  :: z_scalar(n_scalars),gravity_potential(n_scalars)
    real(wp), intent(out) :: theta_v_bg(n_scalars),exner_bg(n_scalars)
  
    ! local variables
    integer(c_int) :: h_index,layer_index,scalar_index
    real(wp) :: temperature,pressure,b,c
  
    !$omp parallel do private(h_index,layer_index,scalar_index,temperature,pressure,b,c)
    do h_index=1,n_scalars_h
      ! integrating from bottom to top
      do layer_index=n_layers-1,0,-1
        scalar_index = layer_index*n_scalars_h+h_index
        temperature = standard_temp(z_scalar(scalar_index))
        ! lowest layer
        if (layer_index==n_layers-1) then
          pressure = standard_pres(z_scalar(scalar_index))
          exner_bg(scalar_index) = (pressure/p_0)**(r_d/c_d_p)
          theta_v_bg(scalar_index) = temperature/exner_bg(scalar_index)
        ! other layers
        else
          ! solving a quadratic equation for the Exner pressure
          b = -0.5_wp*exner_bg(scalar_index + n_scalars_h)/standard_temp(z_scalar(scalar_index+n_scalars_h)) &
          *(temperature - standard_temp(z_scalar(scalar_index + n_scalars_h)) &
          + 2._wp/c_d_p*(gravity_potential(scalar_index) - gravity_potential(scalar_index+n_scalars_h)))
          c = exner_bg(scalar_index+n_scalars_h)**2*temperature/standard_temp(z_scalar(scalar_index+n_scalars_h))
          exner_bg(scalar_index) = b + (b**2 + c)**0.5_wp
          theta_v_bg(scalar_index) = temperature/exner_bg(scalar_index)
        endif
      enddo
    enddo
    !$omp end parallel do
  
  end subroutine set_background_state
  
  subroutine set_area(area,z_vector,z_vector_dual,normal_distance_dual,pent_hex_face_unity_sphere,radius) &
  bind(c,name = "set_area")

    ! This function sets the areas of the grid boxes.
    real(wp), intent(out) :: area(n_vectors)
    real(wp), intent(in)  :: z_vector(n_vectors),z_vector_dual(n_dual_vectors), &
                             normal_distance_dual(n_dual_vectors),pent_hex_face_unity_sphere(n_scalars_h), &
                             radius
  
    ! local variables
    integer(c_int) :: ji,layer_index,h_index,dual_vector_index
    real(wp)       :: base_distance,radius_1,radius_2
    
    !$omp parallel do private(ji,layer_index,h_index,dual_vector_index,base_distance,radius_1,radius_2)
    do ji=1,n_vectors
      layer_index = (ji-1)/n_vectors_per_layer
      h_index = ji-layer_index*n_vectors_per_layer
      if (h_index<=n_scalars_h) then
        area(ji) = pent_hex_face_unity_sphere(h_index)*(radius+z_vector(ji))**2
      else
        dual_vector_index = (layer_index+1)*n_dual_vectors_per_layer + h_index - n_scalars_h
        radius_1 = radius+z_vector_dual(dual_vector_index)
        radius_2 = radius+z_vector_dual(dual_vector_index - n_dual_vectors_per_layer)
        base_distance = normal_distance_dual(dual_vector_index)
        area(ji) = calculate_vertical_area(base_distance,radius_1,radius_2)
      endif
    enddo
    !$omp end parallel do
    
  end subroutine set_area
  
end module vertical_grid













