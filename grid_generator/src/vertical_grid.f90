! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module vertical_grid

  ! This file contains functions that compute properties of the vertical grid.

  use iso_c_binding
  use constants, only: gravity,surface_temp,tropo_height,lapse_rate,inv_height,t_grad_inv,r_d, &
                       p_0_standard
  use grid_nml,  only: no_of_scalars_h,no_of_scalars,no_of_vectors_per_layer, &
                       no_of_vectors,grid_nml_setup
  
  implicit none
  
  private
  
  public :: set_gravity_potential
  public :: set_volume
  
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
      volume(ji) = base_area/(3._c_double*radius_1**2)*(radius_2**3-radius_1**3)
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
    
    tropo_temp_standard = surface_temp-tropo_height*lapse_rate;
    
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
      standard_pres = p_0_standard*(1._c_double-lapse_rate*z_height/surface_temp)**(gravity/(r_d*lapse_rate))
    else if (z_height<inv_height) then
      standard_pres = p_0_standard*(1._c_double-lapse_rate*tropo_height/surface_temp) &
      **(gravity/(r_d*lapse_rate)) &
      *exp(-gravity*(z_height-tropo_height)/(r_d*tropo_temp_standard))
    else
      pressure_at_inv_standard = p_0_standard*(1._c_double-lapse_rate*tropo_height/surface_temp) &
      **(gravity/(r_d*lapse_rate))*exp(-gravity*(inv_height-tropo_height)/(r_d*tropo_temp_standard))
      standard_pres = pressure_at_inv_standard*(1._c_double-lapse_rate*(z_height-inv_height)/surface_temp) &
      **(gravity/(r_d*lapse_rate))
    endif
    
  end function standard_pres

end module vertical_grid













