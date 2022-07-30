! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module vertical_grid

  ! This file contains functions that compute properties of the vertical grid.

  use iso_c_binding
  use definitions, only: wp
  use constants,   only: gravity,surface_temp,tropo_height,lapse_rate,inv_height,t_grad_inv,r_d, &
                         p_0_standard,c_d_p,p_0
  use grid_nml,    only: no_of_scalars_h,no_of_scalars,no_of_vectors_per_layer,no_of_layers, &
                         no_of_vectors,no_of_dual_scalars_h,no_of_dual_scalars,no_of_vectors_h,grid_nml_setup
  
  implicit none
  
  private
  
  public :: set_gravity_potential
  public :: set_volume
  public :: standard_temp
  public :: standard_pres
  public :: set_z_scalar_dual
  public :: set_background_state
  
  contains
  
  subroutine set_gravity_potential(z_scalar,gravity_potential,radius) &
  bind(c,name = "set_gravity_potential")

    ! This subroutine computes the gravity potential.
    real(c_double), intent(in)  :: z_scalar(no_of_scalars)
    real(c_double), intent(out) :: gravity_potential(no_of_scalars)
    real(c_double), intent(in)  :: radius
  
    ! local variables
    integer(c_int) :: ji
    
    call grid_nml_setup()
    
    !$omp parallel do private(ji)
    do ji=1,no_of_scalars
      gravity_potential(ji) = -gravity*(radius**2/(radius+z_scalar(ji))-radius)
    enddo
    !$omp end parallel do
    
  end subroutine set_gravity_potential
  
  
  subroutine set_volume(volume,z_vector,area,radius) &
  bind(c,name = "set_volume")

    ! This subroutine computes the volumes of the grid boxes.
  
    real(c_double), intent(out) :: volume(no_of_scalars)
    real(c_double), intent(in)  :: z_vector(no_of_vectors)
    real(c_double), intent(in)  :: area(no_of_vectors)
    real(c_double), intent(in)  :: radius
  
    ! local variables
    integer(c_int) :: ji,layer_index,h_index
    real(c_double) :: radius_1,radius_2, base_area
    
    call grid_nml_setup()
    
    !$omp parallel do private(ji,layer_index,h_index,radius_1,radius_2,base_area)
    do ji=1,no_of_scalars
      layer_index = (ji-1)/no_of_scalars_h
      h_index = ji-layer_index*no_of_scalars_h
      base_area = area(h_index+(layer_index+1)*no_of_vectors_per_layer)
      radius_1 = radius+z_vector(h_index+(layer_index+1)*no_of_vectors_per_layer)
      radius_2 = radius+z_vector(h_index+layer_index*no_of_vectors_per_layer)
      volume(ji) = base_area/(3._wp*radius_1**2)*(radius_2**3-radius_1**3)
    enddo
    !$omp end parallel do
    
  end subroutine set_volume
  
  function standard_temp(z_height) &
  bind(c,name = "standard_temp")
  
    ! This function returns the temperature in the standard atmosphere.
    
    real(c_double), intent(in) :: z_height
    real(c_double)             :: standard_temp
    
    ! local variables
    real(c_double) :: tropo_temp_standard
    
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
    real(c_double), intent(in) :: z_height
    real(c_double)             :: standard_pres
 
    ! local variables
    real(c_double) :: tropo_temp_standard,pressure_at_inv_standard
    
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
  
    real(c_double), intent(out) :: z_scalar_dual(no_of_dual_scalars)
    real(c_double), intent(in)  :: z_vector(no_of_vectors)
    integer(c_int), intent(in)  :: from_index(no_of_vectors_h),to_index(no_of_vectors_h), &
                                   vorticity_indices_triangles(3*no_of_dual_scalars_h)

    ! local variables
    integer(c_int) :: ji,layer_index,h_index
  
    call grid_nml_setup()
  
    !$omp parallel do private(ji,layer_index,h_index)
    do ji=1,no_of_dual_scalars
      layer_index = (ji-1)/no_of_dual_scalars_h
      h_index = ji-layer_index*no_of_dual_scalars_h
      z_scalar_dual(ji) &
      = 1._wp/6._wp*( &
      z_vector(layer_index*no_of_vectors_per_layer+1+from_index(1+vorticity_indices_triangles(3*(h_index-1)+1))) &
      + z_vector(layer_index*no_of_vectors_per_layer+1+from_index(1+vorticity_indices_triangles(3*(h_index-1)+2))) &
      + z_vector(layer_index*no_of_vectors_per_layer+1+from_index(1+vorticity_indices_triangles(3*(h_index-1)+3))) &
      + z_vector(layer_index*no_of_vectors_per_layer+1+to_index(1+vorticity_indices_triangles(3*(h_index-1)+1))) &
      + z_vector(layer_index*no_of_vectors_per_layer+1+to_index(1+vorticity_indices_triangles(3*(h_index-1)+2))) &
      + z_vector(layer_index*no_of_vectors_per_layer+1+to_index(1+vorticity_indices_triangles(3*(h_index-1)+3))))
    enddo
    !$omp end parallel do
    
  end subroutine set_z_scalar_dual
  
  subroutine set_background_state(z_scalar,gravity_potential,theta_v_bg,exner_bg) &
  bind(c,name = "set_background_state")

    ! This subroutine sets the hydrostatic background state.
    
    real(c_double), intent(in)  :: z_scalar(no_of_scalars),gravity_potential(no_of_scalars)
    real(c_double), intent(out) :: theta_v_bg(no_of_scalars),exner_bg(no_of_scalars)
  
    ! local variables
    integer(c_int) :: h_index,layer_index,scalar_index
    real(c_double) :: temperature,pressure,b,c
  
    call grid_nml_setup()
  
    !$omp parallel do private(h_index,layer_index,scalar_index,temperature,pressure,b,c)
    do h_index=1,no_of_scalars_h
      ! integrating from bottom to top
      do layer_index=no_of_layers-1,0,-1
        scalar_index = layer_index*no_of_scalars_h+h_index
        temperature = standard_temp(z_scalar(scalar_index))
        ! lowest layer
        if (layer_index==no_of_layers-1) then
          pressure = standard_pres(z_scalar(scalar_index))
          exner_bg(scalar_index) = (pressure/p_0)**(r_d/c_d_p)
          theta_v_bg(scalar_index) = temperature/exner_bg(scalar_index)
        ! other layers
        else
          ! solving a quadratic equation for the Exner pressure
          b = -0.5_wp*exner_bg(scalar_index + no_of_scalars_h)/standard_temp(z_scalar(scalar_index+no_of_scalars_h)) &
          *(temperature - standard_temp(z_scalar(scalar_index + no_of_scalars_h)) &
          + 2._wp/c_d_p*(gravity_potential(scalar_index) - gravity_potential(scalar_index+no_of_scalars_h)))
          c = exner_bg(scalar_index+no_of_scalars_h)**2*temperature/standard_temp(z_scalar(scalar_index+no_of_scalars_h))
          exner_bg(scalar_index) = b + (b**2 + c)**0.5_wp
          theta_v_bg(scalar_index) = temperature/exner_bg(scalar_index)
        endif
      enddo
    enddo
    !$omp end parallel do
  
  end subroutine set_background_state

end module vertical_grid













