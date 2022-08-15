! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_vorticity_flux

  ! In this module, the vorticity flux term of the Lamb tansformation gets computed.

  use iso_c_binding
  use definitions, only: wp
  use grid_nml,    only: n_vectors_per_layer,n_layers,n_pentagons,n_vectors_h,n_vectors,n_scalars_h,n_scalars
  
  implicit none
  
  contains

  subroutine vorticity_flux(from_index,to_index,out_field,trsk_indices,trsk_modified_curl_indices,trsk_weights, &
                            mass_flux_density,pot_vorticity,inner_product_weights,adjacent_vector_indices_h) &
  bind(c,name = "vorticity_flux")
    
    ! This subroutine computes the vorticity flux term.
    
    integer,  intent(in)  :: from_index(n_vectors_h),to_index(n_vectors_h),trsk_indices(10*n_vectors_h), &
                             trsk_modified_curl_indices(10*n_vectors_h),adjacent_vector_indices_h(6*n_scalars_h)
    real(wp), intent(in)  :: trsk_weights(10*n_vectors_h),mass_flux_density(n_vectors), &
                             pot_vorticity((2*n_layers+1)*n_vectors_h),inner_product_weights(8*n_scalars)
    real(wp), intent(out) :: out_field(n_vectors)
    
    ! local variables
    integer  :: ji,jk,h_index,layer_index,h_index_shifted,n_edges,mass_flux_base_index,pot_vort_base_index
    real(wp) :: vert_weight
    
    !$omp parallel do private(ji,jk,h_index,layer_index,n_edges,vert_weight,h_index_shifted, &
    !$omp mass_flux_base_index,pot_vort_base_index)
    do h_index=1,n_vectors_per_layer
      do layer_index=0,n_layers
      
        ji = layer_index*n_vectors_per_layer + h_index
        
        ! Calculating the horizontal component of the vorticity flux term.
        ! ----------------------------------------------------------------
        if (h_index>=n_scalars_h+1 .and. layer_index<n_layers) then
          out_field(ji) = 0._wp
          h_index_shifted = h_index - n_scalars_h
          mass_flux_base_index = n_scalars_h + layer_index*n_vectors_per_layer
          pot_vort_base_index = n_vectors_h + layer_index*2*n_vectors_h
          
          ! "Standard" component (vertical potential vorticity times horizontal mass flux density).
          ! ---------------------------------------------------------------------------------------
          ! From_index comes before to_index as usual.
          if (from_index(h_index_shifted)<n_pentagons) then
            do jk=1,4
              out_field(ji) = out_field(ji) &
              + trsk_weights(10*(h_index_shifted-1)+jk) &
              *mass_flux_density(mass_flux_base_index + 1+trsk_indices(10*(h_index_shifted-1)+jk)) &
              *pot_vorticity(pot_vort_base_index + 1+trsk_modified_curl_indices(10*(h_index_shifted-1)+jk))
            enddo
          else
            do jk=1,5
              if (jk==3) then
                out_field(ji) = out_field(ji) &
                + trsk_weights(10*(h_index_shifted-1)+jk) &
                *mass_flux_density(mass_flux_base_index + 1+trsk_indices(10*(h_index_shifted-1)+jk)) &
                *0.5 &
                *(pot_vorticity(pot_vort_base_index + 1+trsk_modified_curl_indices(10*(h_index_shifted-1)+jk)) &
                + pot_vorticity(pot_vort_base_index + h_index_shifted))
              else
                out_field(ji) = out_field(ji) &
                + trsk_weights(10*(h_index_shifted-1)+jk) &
                *mass_flux_density(mass_flux_base_index + 1+trsk_indices(10*(h_index_shifted-1)+jk)) &
                *pot_vorticity(pot_vort_base_index + 1+trsk_modified_curl_indices(10*(h_index_shifted-1)+jk))
              endif
            enddo
          endif
          if (to_index(h_index_shifted)<n_pentagons) then
            do jk=6,10
              out_field(ji) = out_field(ji) &
              + trsk_weights(10*(h_index_shifted-1)+jk) &
              *mass_flux_density(mass_flux_base_index + 1+trsk_indices(10*(h_index_shifted-1)+jk)) &
              *pot_vorticity(pot_vort_base_index + 1+trsk_modified_curl_indices(10*(h_index_shifted-1)+jk))
            enddo
          else
            do jk=6,10
              if (jk==8) then
                out_field(ji) = out_field(ji) &
                + trsk_weights(10*(h_index_shifted-1)+jk) &
                *mass_flux_density(mass_flux_base_index + 1+trsk_indices(10*(h_index_shifted-1)+jk)) &
                *0.5 &
                *(pot_vorticity(pot_vort_base_index + 1+trsk_modified_curl_indices(10*(h_index_shifted-1)+jk)) &
                + pot_vorticity(pot_vort_base_index + h_index_shifted))
              else
                out_field(ji) = out_field(ji) &
                + trsk_weights(10*(h_index_shifted-1)+jk) &
                *mass_flux_density(mass_flux_base_index + 1+trsk_indices(10*(h_index_shifted-1)+jk)) &
                *pot_vorticity(pot_vort_base_index + 1+trsk_modified_curl_indices(10*(h_index_shifted-1)+jk))
              endif
            enddo
          endif
            
          ! Horizontal "non-standard" component (horizontal potential vorticity times vertical mass flux density).
          ! ------------------------------------------------------------------------------------------------------
           
          ! effect of layer above
          out_field(ji) = out_field(ji) &
          - 0.5_wp &
          *inner_product_weights(8*(layer_index*n_scalars_h + from_index(h_index_shifted)) + 7) &
          *mass_flux_density(layer_index*n_vectors_per_layer + 1+from_index(h_index_shifted)) &
          *pot_vorticity(h_index_shifted + layer_index*2*n_vectors_h)
          out_field(ji) = out_field(ji) &
          - 0.5_wp &
          *inner_product_weights(8*(layer_index*n_scalars_h + to_index(h_index_shifted)) + 7) &
          *mass_flux_density(layer_index*n_vectors_per_layer + 1+to_index(h_index_shifted)) &
          *pot_vorticity(h_index_shifted + layer_index*2*n_vectors_h)
            ! effect of layer below
          out_field(ji) = out_field(ji) &
          - 0.5_wp &
          *inner_product_weights(8*(layer_index*n_scalars_h + from_index(h_index_shifted)) + 8) &
          *mass_flux_density((layer_index + 1)*n_vectors_per_layer + 1+from_index(h_index_shifted)) &
          *pot_vorticity(h_index_shifted + (layer_index + 1)*2*n_vectors_h)
          out_field(ji) = out_field(ji) &
          - 0.5_wp &
          *inner_product_weights(8*(layer_index*n_scalars_h + to_index(h_index_shifted)) + 8) &
          *mass_flux_density((layer_index + 1)*n_vectors_per_layer + 1+to_index(h_index_shifted)) &
          *pot_vorticity(h_index_shifted + (layer_index + 1)*2*n_vectors_h)
        ! Calculating the vertical component of the vorticity flux term.
        ! --------------------------------------------------------------
        elseif (h_index<=n_scalars_h) then
          out_field(ji) = 0._wp
          
          ! Determining the vertical acceleration due to the vorticity flux term.
          
          ! determining the number of edges
          n_edges = 6
          if (h_index<=n_pentagons) then
            n_edges = 5
          endif
          ! determining the vertical interpolation weight
          vert_weight = 0.5_wp
          if (layer_index==0 .or. layer_index==n_layers) then
            vert_weight = 1
          endif
          if (layer_index >= 1) then
            do jk=1,n_edges
              out_field(ji) = out_field(ji) &
              + vert_weight &
              *inner_product_weights(8*((layer_index - 1)*n_scalars_h + h_index-1)+jk) &
              *mass_flux_density(n_scalars_h + (layer_index-1)*n_vectors_per_layer+1+adjacent_vector_indices_h(6*(h_index-1)+jk)) &
              *pot_vorticity(layer_index*2*n_vectors_h+1+adjacent_vector_indices_h(6*(h_index-1)+jk))
            enddo
          endif
          if (layer_index<=n_layers-1) then
            do jk=1,n_edges
              out_field(ji) = out_field(ji) &
              + vert_weight &
              *inner_product_weights(8*(layer_index*n_scalars_h + h_index-1)+jk) &
              *mass_flux_density(n_scalars_h + layer_index*n_vectors_per_layer+1+adjacent_vector_indices_h(6*(h_index-1)+jk)) &
              *pot_vorticity(layer_index*2*n_vectors_h+1+adjacent_vector_indices_h(6*(h_index-1)+jk))
            enddo
          endif
        endif
      enddo
    enddo
  
  end subroutine vorticity_flux
  
end module mo_vorticity_flux









