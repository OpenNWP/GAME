/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This file contains functions for reading the grid properties as well as setting the time step.
*/

#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include "../game_types.h"
#include "../game_constants.h"
#include "../spatial_operators/spatial_operators.h"
#include "../constituents/constituents.h"
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

int set_grid_properties(Grid *grid, Dualgrid *dualgrid, char grid_file_name[])
{
	/*
	This function reads all the grid properties from the grid netcdf file.
	*/
    int ncid, retval;
    int normal_distance_id, volume_id, area_id, z_scalar_id, z_vector_id, trsk_weights_id, area_dual_id, z_vector_dual_id, f_vec_id, to_index_id, from_index_id,
    to_index_dual_id, from_index_dual_id, adjacent_vector_indices_h_id, trsk_indices_id, trsk_modified_curl_indices_id, adjacent_signs_h_id, direction_id,
    gravity_potential_id, inner_product_weights_id, density_to_rhombi_weights_id, density_to_rhombi_indices_id, normal_distance_dual_id, vorticity_indices_triangles_id,
    vorticity_signs_triangles_id, latitude_scalar_id, longitude_scalar_id, toa_id, radius_id, interpol_indices_id, interpol_weights_id, theta_v_bg_id,
    exner_bg_id, sfc_rho_c_id, sfc_albedo_id, roughness_length_id, is_land_id, t_conductivity_id, no_of_oro_layers_id, stretching_parameter_id;
    if ((retval = nc_open(grid_file_name, NC_NOWRITE, &ncid)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "no_of_oro_layers", &no_of_oro_layers_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "toa", &toa_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "radius", &radius_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "stretching_parameter", &stretching_parameter_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "normal_distance", &normal_distance_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "volume", &volume_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "area", &area_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "z_scalar", &z_scalar_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "theta_v_bg", &theta_v_bg_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "exner_bg", &exner_bg_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "gravity_potential", &gravity_potential_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "z_vector", &z_vector_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "trsk_weights", &trsk_weights_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "area_dual", &area_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "z_vector_dual", &z_vector_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "f_vec", &f_vec_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "to_index", &to_index_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "to_index_dual", &to_index_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "direction", &direction_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "normal_distance_dual", &normal_distance_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "from_index", &from_index_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "from_index_dual", &from_index_dual_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "adjacent_vector_indices_h", &adjacent_vector_indices_h_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "vorticity_indices_triangles", &vorticity_indices_triangles_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "vorticity_signs_triangles", &vorticity_signs_triangles_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "trsk_indices", &trsk_indices_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "trsk_modified_curl_indices", &trsk_modified_curl_indices_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "adjacent_signs_h", &adjacent_signs_h_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "inner_product_weights", &inner_product_weights_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "density_to_rhombi_weights", &density_to_rhombi_weights_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "density_to_rhombi_indices", &density_to_rhombi_indices_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "latitude_scalar", &latitude_scalar_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "longitude_scalar", &longitude_scalar_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "interpol_indices", &interpol_indices_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "interpol_weights", &interpol_weights_id)))
        ERR(retval);
	if ((retval = nc_inq_varid(ncid, "sfc_rho_c", &sfc_rho_c_id)))
	    ERR(retval);
	if ((retval = nc_inq_varid(ncid, "sfc_albedo", &sfc_albedo_id)))
	    ERR(retval);
	if ((retval = nc_inq_varid(ncid, "roughness_length", &roughness_length_id)))
	    ERR(retval);
	if ((retval = nc_inq_varid(ncid, "t_conductivity", &t_conductivity_id)))
	    ERR(retval);
	if ((retval = nc_inq_varid(ncid, "is_land", &is_land_id)))
	    ERR(retval);
    if ((retval = nc_get_var_int(ncid, no_of_oro_layers_id, &(grid -> no_of_oro_layers))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, toa_id, &(grid -> toa))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, radius_id, &(grid -> radius))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, stretching_parameter_id, &(grid -> stretching_parameter))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, normal_distance_id, &(grid -> normal_distance[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, inner_product_weights_id, &(grid -> inner_product_weights[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, volume_id, &(grid -> volume[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, area_id, &(grid -> area[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, z_scalar_id, &(grid -> z_scalar[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, theta_v_bg_id, &(grid -> theta_v_bg[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, exner_bg_id, &(grid -> exner_bg[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, gravity_potential_id, &(grid -> gravity_potential[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, z_vector_id, &(grid -> z_vector[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, trsk_weights_id, &(grid -> trsk_weights[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, area_dual_id, &(dualgrid -> area[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, z_vector_dual_id, &(dualgrid -> z_vector[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, direction_id, &(grid -> direction[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, f_vec_id, &(dualgrid -> f_vec[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, density_to_rhombi_weights_id, &(grid -> density_to_rhombi_weights[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, normal_distance_dual_id, &(dualgrid -> normal_distance[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, latitude_scalar_id, &(grid -> latitude_scalar[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, longitude_scalar_id, &(grid -> longitude_scalar[0]))))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, interpol_weights_id, &(grid -> latlon_interpol_weights[0]))))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, from_index_id, &(grid -> from_index[0]))))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, to_index_id, &(grid -> to_index[0]))))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, from_index_dual_id, &(dualgrid -> from_index[0]))))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, to_index_dual_id, &(dualgrid -> to_index[0]))))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, adjacent_vector_indices_h_id, &(grid -> adjacent_vector_indices_h[0]))))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, vorticity_indices_triangles_id, &(dualgrid -> vorticity_indices_triangles[0]))))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, vorticity_signs_triangles_id, &(dualgrid -> vorticity_signs_triangles[0]))))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, trsk_indices_id, &(grid -> trsk_indices[0]))))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, trsk_modified_curl_indices_id, &(grid -> trsk_modified_curl_indices[0]))))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, adjacent_signs_h_id, &(grid -> adjacent_signs_h[0]))))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, density_to_rhombi_indices_id, &(grid -> density_to_rhombi_indices[0]))))
        ERR(retval);
    if ((retval = nc_get_var_int(ncid, interpol_indices_id, &(grid -> latlon_interpol_indices[0]))))
        ERR(retval);
	if ((retval = nc_get_var_double(ncid, sfc_rho_c_id, &(grid -> sfc_rho_c[0]))))
	    ERR(retval);
	if ((retval = nc_get_var_double(ncid, sfc_albedo_id, &(grid -> sfc_albedo[0]))))
	    ERR(retval);
	if ((retval = nc_get_var_double(ncid, roughness_length_id, &(grid -> roughness_length[0]))))
	    ERR(retval);
	if ((retval = nc_get_var_double(ncid, t_conductivity_id, &(grid -> t_conduc_soil[0]))))
	    ERR(retval);
	if ((retval = nc_get_var_int(ncid, is_land_id, &(grid -> is_land[0]))))
	    ERR(retval);
    if ((retval = nc_close(ncid)))
        ERR(retval);
    #pragma omp parallel for
    for (int i = 0; i < 6*N_SCALS_H; ++i)
    {
        if (grid -> adjacent_vector_indices_h[i] == -1)
        {
        	grid -> adjacent_vector_indices_h[i] = 0;
        }
    }
    
    // calculating the layer thicknesses
    int layer_index, h_index;
    #pragma omp parallel for private(layer_index, h_index)
    for (int i = 0; i < N_SCALARS; ++i)
    {
    	layer_index = i/N_SCALS_H;
    	h_index = i - layer_index*N_SCALS_H;
    	grid -> layer_thickness[i] = grid -> z_vector[h_index + layer_index*N_VECS_PER_LAYER]
    	- grid -> z_vector[h_index + (layer_index + 1)*N_VECS_PER_LAYER];
    }
	
    // determining coordinate slopes
    grad_hor_cov(grid -> z_scalar, grid -> slope, grid -> from_index, grid -> to_index, grid -> normal_distance);
    // computing the gradient of the gravity potential
    grad(grid -> gravity_potential, grid -> gravity_m, grid);
    // computing the gradient of the background Exner pressure
    grad(grid -> exner_bg, grid -> exner_bg_grad, grid);
	
	// fundamental SFC properties
	grid -> z_t_const = -10.0;
	#pragma omp parallel for
	for (int i = 0; i < N_SCALS_H; ++i)
    {
		grid -> t_const_soil[i] = T_0 + 25.0*cos(2.0*grid -> latitude_scalar[i]);
    }
    
    /*
    constructing the soil grid
    --------------------------
    */
    
	double sigma_soil = 0.352;
	
	// the surface is always at zero
	grid -> z_soil_interface[0] = 0;
	for (int i = 1; i < N_SOIL_LAYERS + 1; ++i)
	{
		grid -> z_soil_interface[i] = grid -> z_soil_interface[i - 1] + pow(sigma_soil, N_SOIL_LAYERS - i);
	}
	double rescale_factor = grid -> z_t_const/grid -> z_soil_interface[N_SOIL_LAYERS];
	for (int i = 1; i < N_SOIL_LAYERS + 1; ++i)
	{
		grid -> z_soil_interface[i] = rescale_factor*grid -> z_soil_interface[i];
	}
	for (int i = 0; i < N_SOIL_LAYERS; ++i)
	{
		grid -> z_soil_center[i] = 0.5*(grid -> z_soil_interface[i] + grid -> z_soil_interface[i + 1]);
	}
    
	return 0;
}








