/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

int find_triangle_on_face_index_from_dual_scalar_on_face_index(int, int, int *, int *, int *, int *);
int find_triangle_edge_points_from_dual_scalar_on_face_index(int, int, int, int *, int *, int *, int [][3], int face_edges[][3], int face_edges_reverse[][3]);
extern int find_points_per_edge();
extern int find_scalar_points_per_inner_face();
extern int upscale_scalar_point();
int write_scalar_coordinates(int, int, int, int, int, int, int, double [], double [], double [], double [], double []);
int find_v_vector_indices_for_dual_scalar_z(int [], int [], int [], int, int []);
extern int find_coords_from_triangle_on_face_index();
extern int find_triangle_on_face_index_from_coords();
int find_triangle_indices_from_h_vector_index(int, int, int *, int *, int *, int *, int *, int *, int *, int *, int [][3], int [][3], int [][3]);
int find_triangle_edge_points(int, int, int, int *, int *, int *, int *, int *, int *, int *, int [][3], int [][3], int [][3]);
int build_icosahedron(double [], double [], int [][2], int [][3], int [][3], int [][3]);
int generate_horizontal_generators(double [], double [], double [], double [], double [], double [], double [], int [][3], int [][3], int [][3]);
int coriolis(int [], int [], int [], double [], double [], int [], double [], double [], double [], double [], double [], double [], double [], double [], double [], int [], int [], int [], double [], double [], double);
int calc_cell_area_unity(double [], double [], double [], int [], int []);
int calc_triangle_area_unity(double [], double [], double [], int [][3], int [][3], int [][3]);
int set_vector_h_doubles(int [], int [], double [], double [], double [], double [], double []);
int set_from_to_index(int [], int [], int [][3], int [][3], int [][3], int [][2]);
int set_scalar_h_dual_coords(double [], double [], double [], double [], int [][3], int [][3], int [][3]);
extern int find_adjacent_vector_indices_h();
int set_horizontal_curl_indices(double [], double [], int [], int [], int [], double, int []);
int rhombus_averaging(int [], int [], int [], int [], int [], int [], int [], int [], double [], double [], double [], double [], double [], double [], double [], double [], double [], double);
extern int set_dual_vector_h_atttributes();
int set_from_to_index_dual(int [], int [], int [][3], int [][3]);
extern int calc_vorticity_indices_triangles();
int optimize_to_scvt(double [], double [], double [], double [], int, int [][3], int [][3], int [][3], int [], int [], int []);
extern int read_horizontal_explicit();
int direct_tangential_unity(double [], double [], double [], double [], int [], int [], double [], double);
extern int write_statistics_file();
extern int set_f_vec();
extern int set_sfc_properties();
extern int interpolate_ll();
extern int set_z_vector_and_normal_distance();
extern int calc_z_vector_dual_and_normal_distance_dual();
extern double rad2deg();
extern double deg2rad();
extern int find_min_index();
extern int find_max_index();
extern int find_min_index_exclude();
extern int in_bool_checker();
extern double calculate_vertical_area();
extern double scalar_product_elementary();
extern double scalar_product_elementary_2d();
extern double find_turn_angle();
extern double calculate_distance_cart();
extern double calculate_distance_h();
extern int active_turn();
extern int passive_turn();
extern int calc_local_i();
extern int calc_local_j();
extern int active_turn_x();
extern int normalize_cartesian();
extern int find_geos();
extern int find_global_normal();
extern int find_between_point();
extern int find_geodetic();
extern double find_geodetic_direction();
extern int find_voronoi_center_sphere();
extern double calc_triangle_area();
extern double rel_on_line();
extern int sort_vertex_indices();
extern double calc_spherical_polygon_area();
extern int set_gravity_potential();
extern int set_volume();
extern int calc_inner_product();
extern int set_z_scalar_dual();
extern int grid_nml_setup();
extern int set_background_state();
extern int set_area();
extern int set_z_scalar();
extern int set_area_dual();














