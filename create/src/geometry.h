#ifndef GEOMETRY_H
#define GEOMETRY_H

bool tri_line_intersect(const double vert1[3],
                        const double vert2[3],
                        const double vert3[3],
                        const double normal[3],
                        const double line_start[3],
                        const double line_end[3],
                              double *intersect_point,
                              double &parametric_value,
                              int    &intersect_flag,
                              int    &side_of_intersect);

int tri_hex_intersect(const double tri_vertex1[3],
                      const double tri_vertex2[3],
                      const double tri_vertex3[3], 
                      const double hex_origin[3],
                      const int    start_i,
                      const int    start_j,
                      const int    start_k,
                      const int    end_i,
                      const int    end_j,
                      const int    end_k,
                      const double delta_x,
                      const double delta_y,
                      const double delta_z,
		            int    *hex_index_array);

bool point_tri_move(const double tri_vertex1[3], 
                    const double tri_vertex2[3],
                    const double tri_vertex3[3],
                    const double start_location[3],
                    const double start_direction[3],
                    const double start_distance,
                          double * new_location,
                          double & distance_moved,
                          int & end_edge_index);
bool push_away(const double vert1[3],
               const double vert2[3],
               const double vert3[3],
               const double point[3],
               const double distance_epsilon,
                     double *new_point);
bool heart_of_tri_line_intersect(const double vert1[3],
                                 const double vert2[3],
                                 const double vert3[3],
                                 const double normal[3],
                                 const double line_start[3],
                                 const double line_end[3],
                                 const int    starting_on_edge,
                                       double *intersect_point,
                                       double &parametric_value,
                                       int    &side_of_intersect);
double distance_to_triangle(const double tri_vertex1[3],
                            const double tri_vertex2[3],
                            const double tri_vertex3[3],
                            const double point[3]);
bool line_line_intersect_2d(const double A[2],
                            const double B[2],
                            const double C[2],
                            const double D[2],
                                  double &t);
double calculate_2d_dot_product(const double point_1[2],
                                const double point_2[2]);
double calculate_3d_dot_product(const double point_1[3],
                                const double point_2[3]);
double normalize_vector(const double input_vector[3], double normal_vector[3]);

double Determinant2x2(const double c1[2], double c2[2]);
void GeneralizedProjectPoint(double x[3], double origin[3],
                             double normal[3], double xproj[3]);
double Distance2BetweenPoints(const double x[3], const double y[3]);
double DistanceToLine(double x[3], double p1[3], double p2[3],
                      double &t, double closestPoint[3]);
void ComputeNormalDirection(const double v1[3], const double v2[3],
                            const double v3[3], double n[3]);
int edge_or_vertex_intersect(const double vert1[3],
                             const double vert2[3],
                             const double vert3[3],
                             const double *intersect_point);


#endif
