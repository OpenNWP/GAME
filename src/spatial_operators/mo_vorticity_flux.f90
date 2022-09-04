! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_vorticity_flux

  ! In this module, the vorticity flux term of the Lamb tansformation gets computed.

  use mo_definitions, only: wp,t_grid,t_diag
  use mo_grid_nml,    only: n_vectors_per_layer,n_layers,n_pentagons,n_vectors_h,n_vectors,n_scalars_h,n_scalars
  
  implicit none
  
  contains

  subroutine vorticity_flux(diag,grid)
    
    ! This subroutine computes the vorticity flux term.
    
    type(t_diag), intent(inout) :: diag ! diagnostic quantities
    type(t_grid), intent(in)    :: grid ! grid quantities
    
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
          diag%pot_vort_tend(ji) = 0._wp
          h_index_shifted = h_index - n_scalars_h
          mass_flux_base_index = n_scalars_h + layer_index*n_vectors_per_layer
          pot_vort_base_index = n_vectors_h + layer_index*2*n_vectors_h
          
          ! "Standard" component (vertical potential vorticity times horizontal mass flux density).
          ! ---------------------------------------------------------------------------------------
          ! from_index comes before to_index as usual.
          if (grid%from_index(h_index_shifted)<n_pentagons) then
            do jk=1,4
              diag%pot_vort_tend(ji) = diag%pot_vort_tend(ji) &
              + grid%trsk_weights(10*(h_index_shifted-1)+jk) &
              *diag%flux_density(mass_flux_base_index + 1+grid%trsk_indices(10*(h_index_shifted-1)+jk)) &
              *diag%pot_vort(pot_vort_base_index + 1+grid%trsk_modified_curl_indices(10*(h_index_shifted-1)+jk))
            enddo
          else
            do jk=1,5
              if (jk==3) then
                diag%pot_vort_tend(ji) = diag%pot_vort_tend(ji) &
                + grid%trsk_weights(10*(h_index_shifted-1)+jk) &
                *diag%flux_density(mass_flux_base_index + 1+grid%trsk_indices(10*(h_index_shifted-1)+jk)) &
                *0.5 &
                *(diag%pot_vort(pot_vort_base_index + 1+grid%trsk_modified_curl_indices(10*(h_index_shifted-1)+jk)) &
                + diag%pot_vort(pot_vort_base_index + h_index_shifted))
              else
                diag%pot_vort_tend(ji) = diag%pot_vort_tend(ji) &
                + grid%trsk_weights(10*(h_index_shifted-1)+jk) &
                *diag%flux_density(mass_flux_base_index + 1+grid%trsk_indices(10*(h_index_shifted-1)+jk)) &
                *diag%pot_vort(pot_vort_base_index + 1+grid%trsk_modified_curl_indices(10*(h_index_shifted-1)+jk))
              endif
            enddo
          endif
          if (grid%to_index(h_index_shifted)<n_pentagons) then
            do jk=6,10
              diag%pot_vort_tend(ji) = diag%pot_vort_tend(ji) &
              + grid%trsk_weights(10*(h_index_shifted-1)+jk) &
              *diag%flux_density(mass_flux_base_index + 1+grid%trsk_indices(10*(h_index_shifted-1)+jk)) &
              *diag%pot_vort(pot_vort_base_index + 1+grid%trsk_modified_curl_indices(10*(h_index_shifted-1)+jk))
            enddo
          else
            do jk=6,10
              if (jk==8) then
                diag%pot_vort_tend(ji) = diag%pot_vort_tend(ji) &
                + grid%trsk_weights(10*(h_index_shifted-1)+jk) &
                *diag%flux_density(mass_flux_base_index + 1+grid%trsk_indices(10*(h_index_shifted-1)+jk)) &
                *0.5 &
                *(diag%pot_vort(pot_vort_base_index + 1+grid%trsk_modified_curl_indices(10*(h_index_shifted-1)+jk)) &
                + diag%pot_vort(pot_vort_base_index + h_index_shifted))
              else
                diag%pot_vort_tend(ji) = diag%pot_vort_tend(ji) &
                + grid%trsk_weights(10*(h_index_shifted-1)+jk) &
                *diag%flux_density(mass_flux_base_index + 1+grid%trsk_indices(10*(h_index_shifted-1)+jk)) &
                *diag%pot_vort(pot_vort_base_index + 1+grid%trsk_modified_curl_indices(10*(h_index_shifted-1)+jk))
              endif
            enddo
          endif
            
          ! Horizontal "non-standard" component (horizontal potential vorticity times vertical mass flux density).
          ! ------------------------------------------------------------------------------------------------------
           
          ! effect of layer above
          diag%pot_vort_tend(ji) = diag%pot_vort_tend(ji) &
          - 0.5_wp &
          *grid%inner_product_weights(8*(layer_index*n_scalars_h + grid%from_index(h_index_shifted)) + 7) &
          *diag%flux_density(layer_index*n_vectors_per_layer + 1+grid%from_index(h_index_shifted)) &
          *diag%pot_vort(h_index_shifted + layer_index*2*n_vectors_h)
          diag%pot_vort_tend(ji) = diag%pot_vort_tend(ji) &
          - 0.5_wp &
          *grid%inner_product_weights(8*(layer_index*n_scalars_h + grid%to_index(h_index_shifted)) + 7) &
          *diag%flux_density(layer_index*n_vectors_per_layer + 1+grid%to_index(h_index_shifted)) &
          *diag%pot_vort(h_index_shifted + layer_index*2*n_vectors_h)
            ! effect of layer below
          diag%pot_vort_tend(ji) = diag%pot_vort_tend(ji) &
          - 0.5_wp &
          *grid%inner_product_weights(8*(layer_index*n_scalars_h + grid%from_index(h_index_shifted)) + 8) &
          *diag%flux_density((layer_index + 1)*n_vectors_per_layer + 1+grid%from_index(h_index_shifted)) &
          *diag%pot_vort(h_index_shifted + (layer_index + 1)*2*n_vectors_h)
          diag%pot_vort_tend(ji) = diag%pot_vort_tend(ji) &
          - 0.5_wp &
          *grid%inner_product_weights(8*(layer_index*n_scalars_h + grid%to_index(h_index_shifted)) + 8) &
          *diag%flux_density((layer_index + 1)*n_vectors_per_layer + 1+grid%to_index(h_index_shifted)) &
          *diag%pot_vort(h_index_shifted + (layer_index + 1)*2*n_vectors_h)
        ! Calculating the vertical component of the vorticity flux term.
        ! --------------------------------------------------------------
        elseif (h_index<=n_scalars_h) then
          diag%pot_vort_tend(ji) = 0._wp
          
          ! Determining the vertical acceleration due to the vorticity flux term.
          
          ! determining the number of edges
          n_edges = 6
          if (h_index<=n_pentagons) then
            n_edges = 5
          endif
          ! determining the vertical interpolation weight
          vert_weight = 0.5_wp
          if (layer_index==0 .or. layer_index==n_layers) then
            vert_weight = 1._wp
          endif
          if (layer_index >= 1) then
            do jk=1,n_edges
              diag%pot_vort_tend(ji) = diag%pot_vort_tend(ji) &
              + vert_weight &
              *grid%inner_product_weights(8*((layer_index - 1)*n_scalars_h + h_index-1)+jk) &
              *diag%flux_density(n_scalars_h+(layer_index-1)*n_vectors_per_layer+ &
                                 1+grid%adjacent_vector_indices_h(6*(h_index-1)+jk)) &
              *diag%pot_vort(layer_index*2*n_vectors_h+1+grid%adjacent_vector_indices_h(6*(h_index-1)+jk))
            enddo
          endif
          if (layer_index<=n_layers-1) then
            do jk=1,n_edges
              diag%pot_vort_tend(ji) = diag%pot_vort_tend(ji) &
              + vert_weight &
              *grid%inner_product_weights(8*(layer_index*n_scalars_h + h_index-1)+jk) &
              *diag%flux_density(n_scalars_h + layer_index*n_vectors_per_layer+ &
                                      1+grid%adjacent_vector_indices_h(6*(h_index-1)+jk)) &
              *diag%pot_vort(layer_index*2*n_vectors_h+1+grid%adjacent_vector_indices_h(6*(h_index-1)+jk))
            enddo
          endif
        endif
      enddo
    enddo
  
  end subroutine vorticity_flux
  
end module mo_vorticity_flux









