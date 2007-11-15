// Geometry routines
// written by Boyd Tidwell & Ray Meyers at Elemental Technologies

#include "stdio.h"

#include "math.h"
#include "geometry.h"

#define VTK_TOL 1.e-05
const double mTolerance = 0.0000;
enum {NONE, OUTSIDE, INSIDE};

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

// ---------------------------------------------------------------------
// detect intersection between a line segment and a triangle
// intersection is defined as any line segment pt (including end pts)
//   in common with any triangle pt (interior, edge, vertex)
// vert1, vert2, vert3 = 3 vertices of triangle
// normal = unit vector normal to triangle plane, pointing OUT
// line_start, line_end = end points of directed line segment
// return TRUE if there is an intersection, else FALSE
// if TRUE:
//   intersect_point = pt of intersection
//   parametric_value = intersection pt is this fraction along segment
//   intersection_flag = 0 for INTERIOR, 
//                       1 for EDGE 1 
//                       2 for EDGE 2 
//                       3 for EDGE 3 
//                       4 for VERTEX 1 
//                       5 for VERTEX 2 
//                       6 for VERTEX 3 
//   side_of_intersect = OUTSIDE or INSIDE (enum value)
// ---------------------------------------------------------------------

bool tri_line_intersect(const double vert1[3],
                        const double vert2[3],
                        const double vert3[3],
                        const double normal[3],
                        const double line_start[3],
                        const double line_end[3],
                        double *     intersect_point,
                        double       &parametric_value,
                        int          &intersect_flag,
                        int          &side_of_intersect)
{
  bool line_intersects; 
  int starting_on_edge = 0;

  starting_on_edge = 
             edge_or_vertex_intersect(vert1, vert2, vert3, line_start);

  line_intersects = 
       heart_of_tri_line_intersect(vert1, vert2, vert3, normal, 
                    line_start, line_end, starting_on_edge, 
                    intersect_point, parametric_value, side_of_intersect);

  // check for edge or vertex intersection
  if (line_intersects)
  {
    intersect_flag = 
             edge_or_vertex_intersect(vert1, vert2, vert3, intersect_point);
  }

  return line_intersects;

}  // end of tri_line_intersect


// ---------------------------------------------------------------------
//  push a point away from an edge when it's too close 
//  Returns:
//    FALSE when no "pushing" is required.
//    TRUE when the point has been moved away from the edge.
// ---------------------------------------------------------------------

bool push_away(const double vert1[3],
               const double vert2[3],
               const double vert3[3],
               const double point[3],
               const double distance_epsilon,
                     double *new_point)
{
  double dist_to_edge1, dist_to_edge2, dist_to_edge3;
  double closest_point1[3], closest_point2[3], closest_point3[3];
  double t;  // not used for anything, just a required parameter
  double centroid[3];
  double center_vector[3];
  double normalized_center_vector[3];
  double new_vector[3];
  double edge_vector[3];
  double distance_to_move;
  double dot_product;
  double edge_vector_magnitude;
  double center_vector_magnitude;
  double distance_to_edge;
  double theta;

  dist_to_edge1 = DistanceToLine((double*)point,(double*)vert1,
                                 (double*)vert2,t,closest_point1);
  dist_to_edge2 = DistanceToLine((double*)point,(double*)vert2,
                                 (double*)vert3,t,closest_point2);
  dist_to_edge3 = DistanceToLine((double*)point,(double*)vert3,
                                 (double*)vert1,t,closest_point3);

  if (dist_to_edge1 < distance_epsilon && dist_to_edge1 > 0.0 || 
      dist_to_edge2 < distance_epsilon && dist_to_edge2 > 0.0 || 
      dist_to_edge3 < distance_epsilon && dist_to_edge3 > 0.0)
  {
    // calculate centroid of triangle
    centroid[0] = (vert1[0] + vert2[0] + vert3[0]) / 3;
    centroid[1] = (vert1[1] + vert2[1] + vert3[1]) / 3;
    centroid[2] = (vert1[2] + vert2[2] + vert3[2]) / 3;

    // calc vector to center of triangle
    center_vector[0] = centroid[0] - point[0];
    center_vector[1] = centroid[1] - point[1];
    center_vector[2] = centroid[2] - point[2];

    // normalize the center pointing vector
    normalize_vector(center_vector, normalized_center_vector);

    if (dist_to_edge1 < distance_epsilon)
    {
      edge_vector[0] = vert2[0] - vert1[0]; 
      edge_vector[1] = vert2[1] - vert1[1]; 
      edge_vector[2] = vert2[2] - vert1[2]; 
      distance_to_edge = dist_to_edge1;
    }
    else if (dist_to_edge2 < distance_epsilon)
    {
      edge_vector[0] = vert3[0] - vert2[0]; 
      edge_vector[1] = vert3[1] - vert2[1]; 
      edge_vector[2] = vert3[2] - vert2[2]; 
      distance_to_edge = dist_to_edge2;
    }
    else
    {
      edge_vector[0] = vert1[0] - vert3[0]; 
      edge_vector[1] = vert1[1] - vert3[1]; 
      edge_vector[2] = vert1[2] - vert3[2]; 
      distance_to_edge = dist_to_edge3;
    }

    dot_product = calculate_3d_dot_product(edge_vector, center_vector);

    edge_vector_magnitude = sqrt(edge_vector[0] * edge_vector[0] +
                                 edge_vector[1] * edge_vector[1] +
                                 edge_vector[2] * edge_vector[2]);

    center_vector_magnitude = sqrt(center_vector[0] * center_vector[0] +
                                   center_vector[1] * center_vector[1] +
                                   center_vector[2] * center_vector[2]);

    theta =
        acos(dot_product / (center_vector_magnitude * edge_vector_magnitude));

    distance_to_move = (distance_epsilon - distance_to_edge) / sin(theta); 

    // calculate new location
    new_vector[0] = normalized_center_vector[0] * distance_to_move;
    new_vector[1] = normalized_center_vector[1] * distance_to_move;
    new_vector[2] = normalized_center_vector[2] * distance_to_move;

    new_point[0] = point[0] + new_vector[0];
    new_point[1] = point[1] + new_vector[1];
    new_point[2] = point[2] + new_vector[2];

    return true;
  }
  else
  {
    return false;
  }
}


// ---------------------------------------------------------------------
//  the main part of tri_line_intersect_processing
// ---------------------------------------------------------------------

bool heart_of_tri_line_intersect(const double vert1[3],
                                 const double vert2[3],
                                 const double vert3[3],
                                 const double normal[3],
                                 const double line_start[3],
                                 const double line_end[3],
                                 const int    starting_on_edge,
                                 double *     intersect_point,
                                 double       &parametric_value,
                                 int          &side_of_intersect)
{
  double u;
  double v;
  double t;
  double * intersect_point_ptr = intersect_point;
  int    edge_number_to_check = 0;
  int i;

  double v_dir[3];
  double vector_dot_product;
  bool   line_intersects = false;
  bool   the_lines_intersect = false;

  side_of_intersect = NONE;

  // caluclate vector direction total length
  v_dir[0] = line_end[0] - line_start[0];
  v_dir[1] = line_end[1] - line_start[1];
  v_dir[2] = line_end[2] - line_start[2];

  // find vectors for edges
  double edge1[3];
  double edge2[3];
  double edge3[3];
  for (i = 0; i < 3; i++) {
    edge1[i] = vert2[i] - vert1[i];
    edge2[i] = vert3[i] - vert1[i];
    edge3[i] = vert3[i] - vert2[i];
  }

  // calculate determinate
  double pvec[3];
  pvec[0] = v_dir[1] * edge2[2] - v_dir[2] * edge2[1];
  pvec[1] = v_dir[2] * edge2[0] - v_dir[0] * edge2[2];
  pvec[2] = v_dir[0] * edge2[1] - v_dir[1] * edge2[0];
  double det = edge1[0] * pvec[0] + edge1[1] * pvec[1] + edge1[2] * pvec[2];

  // if determinate is near zero, the ray is in plane of triangle
  if (det >= -mTolerance && det <= mTolerance) 
  {
    double A[2], B[2], C[2], D[2];  // points for 2D lines used in calculations 

    int start_i, end_i;

    // check for 2D intersection for each tri edge
    for (i=0; i < 3; i++)
    {
      edge_number_to_check = 0;

      switch (i)
      {
        case 0:
          start_i = 1;
          end_i = 2;
          break;
    
        case 1:
          start_i = 2;
          end_i = 0;
          break;
    
        case 2:
          start_i = 0;
          end_i = 1;
          break;
      }

      // check edge 1  
      if (vert1[i] == line_start[i] && vert2[i] == line_end[i] &&
          line_start[i] == line_end[i] && starting_on_edge != 1 &&
          starting_on_edge != 4 && starting_on_edge != 5) 
      {
        C[0] = vert1[start_i];
        C[1] = vert1[end_i];
        D[0] = vert2[start_i];
        D[1] = vert2[end_i];
        edge_number_to_check = 1;
      }

check_edge_2:
      // check edge 2  
      if (vert2[i] == line_start[i] && vert3[i] == line_end[i] &&
          line_start[i] == line_end[i] && edge_number_to_check == 0 && 
          starting_on_edge != 2 && starting_on_edge != 5 && 
          starting_on_edge != 6) 
      {
        C[0] = vert2[start_i];
        C[1] = vert2[end_i];
        D[0] = vert3[start_i];
        D[1] = vert3[end_i];
        edge_number_to_check = 2;
      }

check_edge_3:
      // check edge 3  
      if (vert1[i] == line_start[i] && vert3[i] == line_end[i] &&
          line_start[i] == line_end[i] && edge_number_to_check == 0 &&
          starting_on_edge != 3 && starting_on_edge != 4 && 
          starting_on_edge != 6) 
      {
        C[0] = vert1[start_i];
        C[1] = vert1[end_i];
        D[0] = vert3[start_i];
        D[1] = vert3[end_i];
        edge_number_to_check = 3;
      }

      if (edge_number_to_check > 0) 
      {
        A[0] = line_start[start_i];
        A[1] = line_start[end_i];
        B[0] = line_end[start_i];
        B[1] = line_end[end_i];
 
        the_lines_intersect = line_line_intersect_2d(A, B, C, D, t);

        if (the_lines_intersect)
        {    
          // calc point of intersection
          *intersect_point_ptr = t * line_end[0] + ((1-t) * line_start[0]);
          intersect_point_ptr++;
          *intersect_point_ptr = t * line_end[1] + ((1-t) * line_start[1]);
          intersect_point_ptr++;
          *intersect_point_ptr = t * line_end[2] + ((1-t) * line_start[2]);

          parametric_value = t;  // return paramater

          // figure out if intersection is INSIDE or OUTSIDE
          double line_vector[3];
          line_vector[0] = line_end[0] - line_start[0];          
          line_vector[1] = line_end[1] - line_start[1];          
          line_vector[2] = line_end[2] - line_start[2];          

          vector_dot_product =
                       calculate_3d_dot_product(line_vector, normal); 
          if (vector_dot_product > 0)
            side_of_intersect = INSIDE;
          else
            side_of_intersect = OUTSIDE;

           return true;
        }
        else
        {
          // Lines didn't intersect, see if we need to check the 
          // other edges for this iteration
          if (edge_number_to_check == 1)
          {
            edge_number_to_check = 0;
            goto check_edge_2;
          }
          else 
          {
            if (edge_number_to_check == 2)
            {
              edge_number_to_check = 0;
              goto check_edge_3;
            }
          }
        }
      }
    }
    return false;
  }

  // calculate distance from vert1 to ray origin
  double tvec[3];
  for (i = 0; i < 3; i++)
    tvec[i] = line_start[i] - vert1[i];

  // calculate U parameter and test bounds
  double inv_det = 1.0/det;
  u = (tvec[0] * pvec[0] + tvec[1] * pvec[1] + tvec[2] * pvec[2]) * inv_det;
  if (u >= -mTolerance && u <= mTolerance) {
    u = 0.0;
  }
  else {
    if ((u - 1.0) >= -mTolerance && (u - 1.0) <= mTolerance) { 
      u = 1.0;
    }
  }
  if (u < 0.0 || u > 1.0) {
    return false;
  }

    // calculate V parameter and test bounds
  double qvec[3];
  qvec[0] = tvec[1] * edge1[2] - tvec[2] * edge1[1];
  qvec[1] = tvec[2] * edge1[0] - tvec[0] * edge1[2];
  qvec[2] = tvec[0] * edge1[1] - tvec[1] * edge1[0];
  v = (v_dir[0] * qvec[0] + v_dir[1] * qvec[1] + v_dir[2] * qvec[2]) * inv_det;
  if (v >= -mTolerance && v <= mTolerance) 
    v = 0.0;
  else if ((v - 1.0) >= -mTolerance && (v - 1.0) <= mTolerance) 
         v = 1.0;
  if (v < 0.0 || (u + v - 1.0) >= mTolerance) {
    return false;
  }

    // calculate T, ray intersects triangle
  parametric_value = (edge2[0] * qvec[0] + edge2[1] * qvec[1] +
                      edge2[2] * qvec[2]) * inv_det;
  if (parametric_value >= -mTolerance && parametric_value <= mTolerance) 
    parametric_value = 0.0;
  else if ((parametric_value - 1.0) >= -mTolerance 
                         && (parametric_value - 1.0) <= mTolerance) 
         parametric_value = 1.0;

  if (parametric_value <= 1.0 && parametric_value>= 0.0)
  {
    line_intersects = true;

    // calc point of intersection
    *intersect_point_ptr = parametric_value * line_end[0] +
                         ((1 - parametric_value) * line_start[0]);    
    intersect_point_ptr++;
    *intersect_point_ptr = parametric_value * line_end[1] +
                         ((1 - parametric_value) * line_start[1]); 
    intersect_point_ptr++;   
    *intersect_point_ptr = parametric_value * line_end[2] +
                         ((1 - parametric_value) * line_start[2]);    

    double line_vector[3];
    line_vector[0] = line_end[0] - line_start[0];          
    line_vector[1] = line_end[1] - line_start[1];          
    line_vector[2] = line_end[2] - line_start[2];          

    // calculate side of intersection
    vector_dot_product = 
                 calculate_3d_dot_product(line_vector, normal); 
    if (vector_dot_product > 0)
      side_of_intersect = INSIDE;
    else
      side_of_intersect = OUTSIDE;
  }
  else
  {
    line_intersects = false;
  }
  return line_intersects;
}  // end of heart_of_tri_line_intersect


// ---------------------------------------------------------------------
// compute intersection of a triangle with a 3d array of hex cells
// intersection for a triangle/hex pair is defined as
//   any triangle pt (interior, edge, vertex) in common with
//   any hex pt (interior, face, edge, vertex)
// tri_vertex1, tri_vertex2, tri_vertex3 = 3 vertices of triangle
// hex_origin = origin pt for 3d array of hex cells
// start_i, start_j, start_k = indices of first hex in 3d array
//   0,0,0 is a hex with its lower left corner at the origin
// end_i, end_j, end_k = indices of last hex in 3d array
// delta_x, delta_y, delta_z = size of hex cells in each of 3 dimensions
// return # of bins that are intersected
// return hex_index_array = list of hex IDs for each intersected hex
//   hex ID = integer value >= 0
//   hex at (start_i,start_j,start_k) is ID = 0
//   hex at (end_i,end_j,end_k) is ID = NX * NY * NZ
//     where NX = end_i - start_i + 1, similarly for NY, NZ
//   IDs vary fastest in x, next fastest in y, slowest in z
// ---------------------------------------------------------------------

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
		      int *  hex_index_array)
{

  int i, j, k;  // loop control variables
  int * local_hex_index_array = hex_index_array;
  double x_origin, y_origin, z_origin;
  double x_loop_offset, y_loop_offset, z_loop_offset;
  double current_v_start[3];
  double current_v_end[3];
  double current_x_origin, current_y_origin, current_z_origin;

  double parametric_value = -1.0; 
  double intersect_point[3];
  int intersect_side = NONE;
  bool line_intersects = false;
  double tri_normal[3]; 
  int hex_number;
  int number_of_intersections = 0;

  // calculate normal for triangle
  ComputeNormalDirection(tri_vertex1, tri_vertex2, tri_vertex3, tri_normal);

  double min,max;

  min = tri_vertex1[0];
  min = MIN(min,tri_vertex2[0]);
  min = MIN(min,tri_vertex3[0]);
  int x_hex_min = static_cast<int> ((min-hex_origin[0]) / delta_x);
  if (hex_origin[0] + x_hex_min*delta_x == min) x_hex_min--;
  x_hex_min = MAX(x_hex_min,start_i);

  max = tri_vertex1[0];
  max = MAX(max,tri_vertex2[0]);
  max = MAX(max,tri_vertex3[0]);
  int x_hex_max = static_cast<int> ((max-hex_origin[0]) / delta_x);
  x_hex_max = MIN(x_hex_max,end_i);

  min = tri_vertex1[1];
  min = MIN(min,tri_vertex2[1]);
  min = MIN(min,tri_vertex3[1]);
  int y_hex_min = static_cast<int> ((min-hex_origin[1]) / delta_y);
  if (hex_origin[1] + y_hex_min*delta_y == min) y_hex_min--;
  y_hex_min = MAX(y_hex_min,start_j);

  max = tri_vertex1[1];
  max = MAX(max,tri_vertex2[1]);
  max = MAX(max,tri_vertex3[1]);
  int y_hex_max = static_cast<int> ((max-hex_origin[1]) / delta_y);
  y_hex_max = MIN(y_hex_max,end_j);

  min = tri_vertex1[2];
  min = MIN(min,tri_vertex2[2]);
  min = MIN(min,tri_vertex3[2]);
  int z_hex_min = static_cast<int> ((min-hex_origin[2]) / delta_z);
  if (hex_origin[2] + z_hex_min*delta_z == min) z_hex_min--;
  z_hex_min = MAX(z_hex_min,start_k);

  max = tri_vertex1[2];
  max = MAX(max,tri_vertex2[2]);
  max = MAX(max,tri_vertex3[2]);
  int z_hex_max = static_cast<int> ((max-hex_origin[2]) / delta_z);
  z_hex_max = MIN(z_hex_max,end_k);

  x_origin = hex_origin[0];
  y_origin = hex_origin[1];
  z_origin = hex_origin[2];

  int ny = end_j - start_j + 1;
  int nx = end_i - start_i + 1;

  // find all tri/hex edge intersections
  for (k=z_hex_min; k <= z_hex_max; k++)  // z direction
  {
    z_loop_offset = k * delta_z;
    for (j=y_hex_min; j <= y_hex_max; j++)  // y direction
    {
      y_loop_offset = j * delta_y;
      for (i=x_hex_min; i <= x_hex_max; i++)  // x direction
      {
        // check hex edges for intersection
        x_loop_offset = i * delta_x;

	hex_number = (k - start_k)*ny*nx + (j - start_j)*nx + (i - start_i);

        // check edge 1
        current_v_start[0] = x_origin + x_loop_offset;
        current_v_start[1] = y_origin + y_loop_offset;
        current_v_start[2] = z_origin + z_loop_offset;
        current_v_end[0] = x_origin + x_loop_offset + delta_x;
        current_v_end[1] = y_origin + y_loop_offset;
        current_v_end[2] = z_origin + z_loop_offset;
        line_intersects = heart_of_tri_line_intersect(
                       tri_vertex1, tri_vertex2, tri_vertex3,
                       tri_normal, current_v_start, current_v_end, 
                       0, intersect_point, parametric_value, intersect_side);
        if (line_intersects)
        {
          // record hex number in output array
          *local_hex_index_array = hex_number;
          local_hex_index_array++;
          number_of_intersections++; 
          continue;   // no need to check the remaining edges
        }

        // check edge 2
        current_v_start[0] = x_origin  + x_loop_offset + delta_x;
        current_v_end[1] = y_origin + y_loop_offset + delta_y;
        line_intersects = heart_of_tri_line_intersect(
                       tri_vertex1, tri_vertex2, tri_vertex3,
                       tri_normal, current_v_start, current_v_end, 
                       0, intersect_point, parametric_value, intersect_side);
        if (line_intersects)
        {
          // record hex number in output array
          *local_hex_index_array = hex_number;
          local_hex_index_array++;
          number_of_intersections++; 
          continue;   // no need to check the remaining edges
        }

        // check edge 3
        current_v_start[1] = y_origin + y_loop_offset + delta_y;
        current_v_end[0] = x_origin + x_loop_offset;
        line_intersects = heart_of_tri_line_intersect(
                       tri_vertex1, tri_vertex2, tri_vertex3,
                       tri_normal, current_v_start, current_v_end, 
                       0, intersect_point, parametric_value, intersect_side);
        if (line_intersects)
        {
          // record hex number in output array
          *local_hex_index_array = hex_number;
          local_hex_index_array++;
          number_of_intersections++; 
          continue;   // no need to check the remaining edges
        }

        // check edge 4
        current_v_start[0] = x_origin + x_loop_offset;
        current_v_end[1] = y_origin + y_loop_offset;
        line_intersects = heart_of_tri_line_intersect(
                       tri_vertex1, tri_vertex2, tri_vertex3,
                       tri_normal, current_v_start, current_v_end, 
                       0, intersect_point, parametric_value, intersect_side);
        if (line_intersects)
        {
          // record hex number in output array
          *local_hex_index_array = hex_number;
          local_hex_index_array++;
          number_of_intersections++; 
          continue;   // no need to check the remaining edges
        }

        // check edge 5
        current_v_start[1] = y_origin + y_loop_offset;
        current_v_start[2] = z_origin + z_loop_offset + delta_z;
        current_v_end[0] = x_origin + x_loop_offset + delta_x;
        current_v_end[2] = z_origin + z_loop_offset + delta_z;
        line_intersects = heart_of_tri_line_intersect(
                       tri_vertex1, tri_vertex2, tri_vertex3,
                       tri_normal, current_v_start, current_v_end, 
                       0, intersect_point, parametric_value, intersect_side);
        if (line_intersects)
        {
          // record hex number in output array
          *local_hex_index_array = hex_number;
          local_hex_index_array++;
          number_of_intersections++; 
          continue;   // no need to check the remaining edges
        }

        // check edge 6
        current_v_start[0] = x_origin + x_loop_offset + delta_x;
        current_v_end[1] = y_origin + y_loop_offset + delta_y;
        line_intersects = heart_of_tri_line_intersect(
                       tri_vertex1, tri_vertex2, tri_vertex3,
                       tri_normal, current_v_start, current_v_end, 
                       0, intersect_point, parametric_value, intersect_side);
        if (line_intersects)
        {
          // record hex number in output array
          *local_hex_index_array = hex_number;
          local_hex_index_array++;
          number_of_intersections++; 
          continue;   // no need to check the remaining edges
        }

        // check edge 7
        current_v_start[1] = y_origin + y_loop_offset + delta_y;
        current_v_end[0] = x_origin + x_loop_offset;
        line_intersects = heart_of_tri_line_intersect(
                       tri_vertex1, tri_vertex2, tri_vertex3,
                       tri_normal, current_v_start, current_v_end, 
                       0, intersect_point, parametric_value, intersect_side);
        if (line_intersects)
        {
          // record hex number in output array
          *local_hex_index_array = hex_number;
          local_hex_index_array++;
          number_of_intersections++; 
          continue;   // no need to check the remaining edges
        }

        // check edge 8
        current_v_start[0] = x_origin + x_loop_offset;
        current_v_end[1] = y_origin + y_loop_offset;
        line_intersects = heart_of_tri_line_intersect(
                       tri_vertex1, tri_vertex2, tri_vertex3,
                       tri_normal, current_v_start, current_v_end, 
                       0, intersect_point, parametric_value, intersect_side);
        if (line_intersects)
        {
          // record hex number in output array
          *local_hex_index_array = hex_number;
          local_hex_index_array++;
          number_of_intersections++; 
          continue;   // no need to check the remaining edges
        }

        // check edge 9
        current_v_start[1] = y_origin + y_loop_offset;
        current_v_start[2] = z_origin + z_loop_offset;
        // end point is same as for edge 8
        line_intersects = heart_of_tri_line_intersect(
                       tri_vertex1, tri_vertex2, tri_vertex3,
                       tri_normal, current_v_start, current_v_end, 
                       0, intersect_point, parametric_value, intersect_side);
        if (line_intersects)
        {
          // record hex number in output array
          *local_hex_index_array = hex_number;
          local_hex_index_array++;
          number_of_intersections++; 
          continue;   // no need to check the remaining edges
        }

        // check edge 10
        current_v_start[0] = x_origin + x_loop_offset + delta_x;
        current_v_end[0] = x_origin + x_loop_offset + delta_x;
        line_intersects = heart_of_tri_line_intersect(
                       tri_vertex1, tri_vertex2, tri_vertex3,
                       tri_normal, current_v_start, current_v_end, 
                       0, intersect_point, parametric_value, intersect_side);
        if (line_intersects)
        {
          // record hex number in output array
          *local_hex_index_array = hex_number;
          local_hex_index_array++;
          number_of_intersections++; 
          continue;   // no need to check the remaining edges
        }

        // check edge 11
        current_v_start[1] = y_origin + y_loop_offset + delta_y;
        current_v_end[1] = y_origin + y_loop_offset + delta_y;
        line_intersects = heart_of_tri_line_intersect(
                       tri_vertex1, tri_vertex2, tri_vertex3,
                       tri_normal, current_v_start, current_v_end, 
                       0, intersect_point, parametric_value, intersect_side);
        if (line_intersects)
        {
          // record hex number in output array
          *local_hex_index_array = hex_number;
          local_hex_index_array++;
          number_of_intersections++; 
          continue;   // no need to check the remaining edges
        }

        // check edge 12
        current_v_start[0] = x_origin + x_loop_offset;
        current_v_end[0] = x_origin + x_loop_offset;
        line_intersects = heart_of_tri_line_intersect(
                       tri_vertex1, tri_vertex2, tri_vertex3,
                       tri_normal, current_v_start, current_v_end, 
                       0, intersect_point, parametric_value, intersect_side);

        if (line_intersects)
        {
          // record hex number in output array
          *local_hex_index_array = hex_number;
          local_hex_index_array++;
          number_of_intersections++; 
          continue;
        }

        // check for vertices inside hex
        current_x_origin = x_origin + x_loop_offset;
        current_y_origin = y_origin + y_loop_offset;
        current_z_origin = z_origin + z_loop_offset;

        // check vertex 1
        if ((current_x_origin <= tri_vertex1[0] && 
                    tri_vertex1[0] <= current_x_origin + delta_x) &&
            (current_y_origin <= tri_vertex1[1] && 
                    tri_vertex1[1] <= current_y_origin + delta_y) &&
            (current_z_origin <= tri_vertex1[2] && 
                    tri_vertex1[2] <= current_z_origin + delta_z) )
        {
          // record hex number in output array
          *local_hex_index_array = hex_number;
          local_hex_index_array++;
          number_of_intersections++; 
          continue;
        }

        // check vertex 2
        if ((current_x_origin <= tri_vertex2[0] && 
                    tri_vertex2[0] <= current_x_origin + delta_x) &&
            (current_y_origin <= tri_vertex2[1] && 
                    tri_vertex2[1] <= current_y_origin + delta_y) &&
            (current_z_origin <= tri_vertex2[2] && 
                    tri_vertex2[2] <= current_z_origin + delta_z) )
        {
          // record hex number in output array
          *local_hex_index_array = hex_number;
          local_hex_index_array++;
          number_of_intersections++; 
          continue;
        }

        // check vertex 3
        if ((current_x_origin <= tri_vertex3[0] && 
                    tri_vertex3[0] <= current_x_origin + delta_x) &&
            (current_y_origin <= tri_vertex3[1] && 
                    tri_vertex3[1] <= current_y_origin + delta_y) &&
            (current_z_origin <= tri_vertex3[2] && 
                    tri_vertex3[2] <= current_z_origin + delta_z) )
        {
          // record hex number in output array
          *local_hex_index_array = hex_number;
          local_hex_index_array++;
          number_of_intersections++; 
          continue;
        }
      }  // end of x direction loop
    }  // end of y direction loop
  }  // end of z direction loop
  return number_of_intersections;

}  // end of tri_hex_intersect 

// ---------------------------------------------------------------------
// return nearest distance from a point to a triangle
// "nearest distance" is the minimum distance from the point
//   to any point on the triangle (interior, edge, vertex)
// The function distance_to_triangle was adapted from the method
//   vtkTriangle::EvaluatePosition in the Visualization Toolkit (VTK)
//   which is Copyright (c) 1993-2001 Ken Martin, Will Schroeder, Bill Lorensen
//   All rights reserved.
// ---------------------------------------------------------------------

double distance_to_triangle(const double tri_vertex1[3], 
                            const double tri_vertex2[3],
                            const double tri_vertex3[3], 
                            const double point[3]) 
{
  int i, j;
  double n[3], fabsn;
  double rhs[2], c1[2], c2[2];
  double pcoords[3] = {0.0, 0.0, 0.0};;
  double dist2;
  double distance;
  double det;
  double maxComponent;
  int idx=0, indices[2];
  double dist2Point, dist2Line1, dist2Line2;
  double *closest, closestPoint1[3], closestPoint2[3], cp[3];
  double closestPoint[3]; 

  ComputeNormalDirection(tri_vertex1, tri_vertex2, tri_vertex3, n);

  //
  // Project point to plane
  //
  GeneralizedProjectPoint((double*)point,(double*)tri_vertex1,n,cp);

  //
  // Construct matrices.  Since we have over determined system, need to find
  // which 2 out of 3 equations to use to develop equations. (Any 2 should
  // work since we've projected point to plane.)
  //
  for (maxComponent=0.0, i=0; i<3; i++)
    {
    // trying to avoid an expensive call to fabs()
    if (n[i] < 0)
      {
      fabsn = -n[i];
      }
    else
      {
      fabsn = n[i];
      }
    if (fabsn > maxComponent)
      {
      maxComponent = fabsn;
      idx = i;
      }
    }
  for (j=0, i=0; i<3; i++)
    {
    if ( i != idx )
      {
      indices[j++] = i;
      }
    }

  for (i=0; i<2; i++)
    { 
    rhs[i] = cp[indices[i]] - tri_vertex3[indices[i]];
    c1[i] = tri_vertex1[indices[i]] - tri_vertex3[indices[i]];
    c2[i] = tri_vertex2[indices[i]] - tri_vertex3[indices[i]];
    }

  if ( (det = Determinant2x2(c1,c2)) == 0.0 )
    {
    return -1;
    }

  pcoords[0] = Determinant2x2(rhs,c2) / det;
  pcoords[1] = Determinant2x2(c1,rhs) / det;
  pcoords[2] = 1.0 - (pcoords[0] + pcoords[1]);
  //
  // Okay, now find closest point to element
  //
  if ( pcoords[0] >= 0.0 && pcoords[0] <= 1.0 &&
       pcoords[1] >= 0.0 && pcoords[1] <= 1.0 &&
       pcoords[2] >= 0.0 && pcoords[2] <= 1.0 )
    {
    //projection distance
    dist2 = Distance2BetweenPoints(cp,point);
    closestPoint[0] = cp[0];
    closestPoint[1] = cp[1];
    closestPoint[2] = cp[2];
    distance = sqrt(dist2);
    return distance;
    }
  else
    {
    double t;
    if ( pcoords[0] < 0.0 && pcoords[1] < 0.0 )
      {
      dist2Point = Distance2BetweenPoints(point,tri_vertex3);
      dist2Line1 = 
            DistanceToLine((double*)point,(double*)tri_vertex1,
                                 (double*)tri_vertex3,t,closestPoint1);
      dist2Line2 =
            DistanceToLine((double*)point,(double*)tri_vertex3,
                               (double*)tri_vertex2,t,closestPoint2);
      if (dist2Point < dist2Line1)
        {
        dist2 = dist2Point;
        closest = (double*)tri_vertex3;
        }
      else
        {
        dist2 = dist2Line1;
        closest = closestPoint1;
        }
      if (dist2Line2 < dist2)
        {
        dist2 = dist2Line2;
        closest = closestPoint2;
        }
      for (i=0; i<3; i++)
        {
        closestPoint[i] = closest[i];
        }
      }
    else if ( pcoords[1] < 0.0 && pcoords[2] < 0.0 )
      {
      dist2Point = Distance2BetweenPoints(point,tri_vertex1);
      dist2Line1 =
           DistanceToLine((double*)point,(double*)tri_vertex1,
                          (double*)tri_vertex3,t,closestPoint1);
      dist2Line2 =
           DistanceToLine((double*)point,(double*)tri_vertex1,
                          (double*)tri_vertex2,t,closestPoint2);
      if (dist2Point < dist2Line1)
        {
        dist2 = dist2Point;
        closest = (double*)tri_vertex1;
        }
      else
        {
        dist2 = dist2Line1;
        closest = closestPoint1;
        }
      if (dist2Line2 < dist2)
        {
        dist2 = dist2Line2;
        closest = closestPoint2;
        }
      for (i=0; i<3; i++)
        {
        closestPoint[i] = closest[i];
        }
      }
    else if ( pcoords[0] < 0.0 && pcoords[2] < 0.0 )
      {
      dist2Point = Distance2BetweenPoints(point,tri_vertex2);
      dist2Line1 =
            DistanceToLine((double*)point,(double*)tri_vertex2,
                            (double*)tri_vertex3,t,closestPoint1);
        dist2Line2 =
              DistanceToLine((double*)point,(double*)tri_vertex1,
                             (double*)tri_vertex2,t,closestPoint2);
      if (dist2Point < dist2Line1)
        {
        dist2 = dist2Point;
        closest = (double*)tri_vertex2;
        }
      else
        {
        dist2 = dist2Line1;
        closest = closestPoint1;
        }
      if (dist2Line2 < dist2)
        {
        dist2 = dist2Line2;
        closest = closestPoint2;
        }
      for (i=0; i<3; i++)
        {
        closestPoint[i] = closest[i];
        }
      }
    else if ( pcoords[0] < 0.0 )
      {
      dist2 = DistanceToLine((double*)point,(double*)tri_vertex2,
                             (double*)tri_vertex3,t,closestPoint);
      }
    else if ( pcoords[1] < 0.0 )
      {
      dist2 = DistanceToLine((double*)point,(double*)tri_vertex1,
                             (double*)tri_vertex3,t,closestPoint);
      }
    else if ( pcoords[2] < 0.0 )
      {
      dist2 = DistanceToLine((double*)point, (double*)tri_vertex1,
                             (double*)tri_vertex2,t,closestPoint);
      }
    distance = sqrt(dist2);
    return distance;
    }
}


// ---------------------------------------------------------------------
// given a triangle and a point, move the point along the surface a
//   specified direction and distance
// Inputs:
//     tri_vertex1, tri_vertex2, tri_vertex3 = 3 vertices of triangle
//     start_location = location of the point to be moved
//     start_direction = direction to move the point
//     start_distance = distance to move the point over the triangle
// Outputs:
//     new_location = the moved location of the point
//     distance_moved = the actual distance the point was moved.  May differ
//        from start_direction in case where an edge was encountered before
//        moving the start_distance.
//     end_edge_index = the index of the edge or vertex which the new location
//        lies on.  Possible values are:  0 for INTERIOR, 
//                                        1 for EDGE 1 
//                                        2 for EDGE 2 
//                                        3 for EDGE 3 
//                                        4 for VERTEX 1 
//                                        5 for VERTEX 2 
//                                        6 for VERTEX 3 
// Return value:
//     returns TRUE if the point could be successfully moved, FALSE otherwise
// ---------------------------------------------------------------------
bool point_tri_move(const double tri_vertex1[3],
                    const double tri_vertex2[3],
                    const double tri_vertex3[3],
                    const double start_location[3],
                    const double start_direction[3],
                    const double start_distance,
                          double * new_location,
                          double & distance_moved,
                          int & end_edge_index)
{
  double result;
  double new_vector[3];
  double expected_new_location[3];
  bool   line_intersects = false;
  double parametric_value = -1.0; 
  double intersect_point[3];
  int    intersect_side = NONE;

  // check for new location inside triangle
  // calculate expected new_location
  new_vector[0] = (start_direction[0] * start_distance);
  new_vector[1] = (start_direction[1] * start_distance);
  new_vector[2] = (start_direction[2] * start_distance);

  expected_new_location[0] = start_location[0] + new_vector[0];
  expected_new_location[1] = start_location[1] + new_vector[1];
  expected_new_location[2] = start_location[2] + new_vector[2];

  result = distance_to_triangle(tri_vertex1, tri_vertex2, tri_vertex3, 
                                expected_new_location);  
  if (result == 0.0)   // new location is on triangle
  {
    distance_moved = 
        sqrt(Distance2BetweenPoints(start_location, expected_new_location));

    new_location[0] = expected_new_location[0];
    new_location[1] = expected_new_location[1];
    new_location[2] = expected_new_location[2];

    end_edge_index = 
      edge_or_vertex_intersect(tri_vertex1, tri_vertex2, tri_vertex3, new_location);
  }
  else  // expected new location is outside tri, we will hit an edge
  {  
    line_intersects = tri_line_intersect(
                      tri_vertex1, tri_vertex2, tri_vertex3,
                      start_direction, start_location, expected_new_location, 
                      intersect_point, parametric_value,
                      end_edge_index, intersect_side);

    // calculate distnace to new point on edge or vertex
    new_location[0] = intersect_point[0];
    new_location[1] = intersect_point[1];
    new_location[2] = intersect_point[2];
  }

  distance_moved = sqrt(Distance2BetweenPoints(start_location, new_location));

  return line_intersects;
}


bool line_line_intersect_2d(const double A[2], 
                           const double B[2],
                           const double C[2],
                           const double D[2],
                           double & t)
{
  double b_vector[2], c_vector[2], d_vector[2];
  double b_normal[2], d_normal[2];
  double d_normal_dot_b, b_normal_dot_c, d_normal_dot_c;

  // calculate the vector from A to B (B-A)

  b_vector[0] = B[0] - A[0];
  b_vector[1] = B[1] - A[1];

  // calculate normal of vector of b
  b_normal[0] = -b_vector[1];
  b_normal[1] = b_vector[0];

  // calculate the vector from C to D (D-C)
   d_vector[0] = D[0] - C[0];
   d_vector[1] = D[1] - C[1];
 
   // calculate normal of vector of d
   d_normal[0] = -d_vector[1];
   d_normal[1] = d_vector[0];

   // calculate the vector from A to C (C-A)
   c_vector[0] = C[0] - A[0];
   c_vector[1] = C[1] - A[1];

   // calculate dot product of d_normal and b
   d_normal_dot_b = calculate_2d_dot_product(d_normal, b_vector); 

   if (fabs(d_normal_dot_b) >= mTolerance)  // lines are not parallel
   {
     // calculate dot product of d_normal and c
     d_normal_dot_c =  calculate_2d_dot_product(d_normal, c_vector);

     double u;  // parametric value for intersection calculations
     // calculate the parametric parameter t for A + b*t.
     // This is where on line 1 the two lines intersect.
     t = d_normal_dot_c / d_normal_dot_b;

     // calculate dot product of b_normal and c
     b_normal_dot_c =  calculate_2d_dot_product(b_normal, c_vector);

     // calculate the parametric parameter u for C + d*u.
     // This is where on line 1 the two lines intersect.
     u = b_normal_dot_c / d_normal_dot_b;

     // if both t and u are inside of [0.0, 1.0] then the lines intersect
     if ((t > -mTolerance && t <= 1.0 + mTolerance) &&
         (u > -mTolerance && u <= 1.0 + mTolerance))
     {
       // lines intersect, calculate other info to return

       return true;
     }
   }
  return false;
}            


// ---------------------------------------------------------------------
// dot product between two 2-vectors
// ---------------------------------------------------------------------

double calculate_2d_dot_product(const double point_1[2],
                                const double point_2[2])
{
  double result;

  result = (point_1[0] * point_2[0]) + (point_1[1] * point_2[1]);
  return result;
}

// ---------------------------------------------------------------------
// dot product between two 3-vectors
// ---------------------------------------------------------------------

double calculate_3d_dot_product(const double point_1[3],
                                const double point_2[3])
{
  double result;

  result = (point_1[0] * point_2[0]) + (point_1[1] * point_2[1]) +
           (point_1[2] * point_2[2]);
  return result;
}


// ---------------------------------------------------------------------
// normalize a vector
// ---------------------------------------------------------------------
double normalize_vector(const double input_vector[3], double normal_vector[3])
{
  double normal_len;
  normal_len = sqrt(input_vector[0] * input_vector[0] +
                    input_vector[1] * input_vector[1] +
                    input_vector[2] * input_vector[2]);
  if (normal_len != 0)
  {
    normal_vector[0] = input_vector[0] / normal_len;
    normal_vector[1] = input_vector[1] / normal_len;
    normal_vector[2] = input_vector[2] / normal_len;
    return normal_len;
  }
  else
  {
    return 0.0;
  }
}


// ---------------------------------------------------------------------
// determinant of 2x2 matrix - two columns of matrix are input
// ---------------------------------------------------------------------

double Determinant2x2(const double c1[2], double c2[2]) 
{
  return (c1[0]*c2[1] - c2[0]*c1[1]);
}

// ---------------------------------------------------------------------
// Project a point x onto plane defined by origin and normal
// ---------------------------------------------------------------------

void GeneralizedProjectPoint(double x[3], double origin[3],
                             double normal[3], double xproj[3])
{
  double t, xo[3], n2;

  xo[0] = x[0] - origin[0];
  xo[1] = x[1] - origin[1];
  xo[2] = x[2] - origin[2];

  t = calculate_3d_dot_product(normal,xo);
  n2 = calculate_3d_dot_product(normal, normal);

  if (n2 != 0)
  {
    xproj[0] = x[0] - t * normal[0]/n2;
    xproj[1] = x[1] - t * normal[1]/n2;
    xproj[2] = x[2] - t * normal[2]/n2;
  }
  else
  {
    xproj[0] = x[0];
    xproj[1] = x[1];
    xproj[2] = x[2];
  }
}

// ---------------------------------------------------------------------
// Cartesian distance squared between 2 points
// ---------------------------------------------------------------------

double Distance2BetweenPoints(const double x[3], const double y[3])
{
  return ((x[0]-y[0])*(x[0]-y[0]) + (x[1]-y[1])*(x[1]-y[1]) +
          (x[2]-y[2])*(x[2]-y[2]));
}

// ---------------------------------------------------------------------
// Compute distance to finite line. Returns parametric coordinate t
// and point location on line.
// ---------------------------------------------------------------------

double DistanceToLine(double x[3], double p1[3], double p2[3],
                      double &t, double closestPoint[3])
{
  double p21[3], denom, num;
  double *closest;
  double tolerance;
  //
  //   Determine appropriate vectors
  //
  p21[0] = p2[0]- p1[0];
  p21[1] = p2[1]- p1[1];
  p21[2] = p2[2]- p1[2];

  //
  //   Get parametric location
  //
  num = p21[0]*(x[0]-p1[0]) + p21[1]*(x[1]-p1[1]) + p21[2]*(x[2]-p1[2]);
  denom = calculate_3d_dot_product(p21,p21);

  // trying to avoid an expensive fabs
  tolerance = VTK_TOL*num;
  if (tolerance < 0.0)
    {
    tolerance = -tolerance;
    }
  if ( -tolerance < denom && denom < tolerance ) //numerically bad!
    {
    closest = p1; //arbitrary, point is (numerically) far away
    }
  //
  // If parametric coordinate is within 0<=p<=1, then the point is closest to
  // the line.  Otherwise, it's closest to a point at the end of the line.
  //
  else if ( (t=num/denom) < 0.0 )
    {
    closest = p1;
    }
  else if ( t > 1.0 )
    {
    closest = p2;
    }
  else
    {
    closest = p21;
    p21[0] = p1[0] + t*p21[0];
    p21[1] = p1[1] + t*p21[1];
    p21[2] = p1[2] + t*p21[2];
    }

  closestPoint[0] = closest[0];
  closestPoint[1] = closest[1];
  closestPoint[2] = closest[2];
  return sqrt(Distance2BetweenPoints(closest,x));
}

// ---------------------------------------------------------------------
// Compute normal direction to a triangle
// v1,v2,v3 = 3 vertices of triangle
// return n = OUTWARD normal, NOT normalized to unit length
// 3 vertices are ordered so that right-hand rule points outward
// ---------------------------------------------------------------------

void ComputeNormalDirection(const double v1[3], const double v2[3],
                            const double v3[3], double n[3])
{
  double ax, ay, az, bx, by, bz;

  // order is important!!! maintain consistency with triangle vertex order
  ax = v3[0] - v2[0]; ay = v3[1] - v2[1]; az = v3[2] - v2[2];
  bx = v1[0] - v2[0]; by = v1[1] - v2[1]; bz = v1[2] - v2[2];

  n[0] = (ay * bz - az * by);
  n[1] = (az * bx - ax * bz);
  n[2] = (ax * by - ay * bx);
}


// ---------------------------------------------------------------------
// Determine if a given point is on an edge or vertex of a triangle
// ---------------------------------------------------------------------

int edge_or_vertex_intersect(const double vert1[3],
                             const double vert2[3],
                             const double vert3[3],
                             const double *intersect_point)
{
  double dist_to_edge1, dist_to_edge2, dist_to_edge3;
  double closest_point1[3], closest_point2[3], closest_point3[3];
  double t;  // not used for anything, just a required parameter
  int edge_or_vertex = 0;

  // check for edge or vertex intersection
  edge_or_vertex = 0;  // no edge or vertex yet
  dist_to_edge1 = 
         DistanceToLine((double*)intersect_point,(double*)vert1,
                        (double*)vert2,t,closest_point1);
  dist_to_edge2 = 
         DistanceToLine((double*)intersect_point,(double*)vert2,
                        (double*)vert3,t,closest_point2);
  dist_to_edge3 = 
         DistanceToLine((double*)intersect_point,(double*)vert3,
                        (double*)vert1,t,closest_point3);
  if (dist_to_edge1 == 0.0) {
    if (dist_to_edge2 == 0.0) {
      edge_or_vertex = 5;  // vertex 2
    }
    else if (dist_to_edge3 == 0.0) {
       edge_or_vertex = 4;  // vertex 1
    }
    else
    {
      edge_or_vertex = 1;  // edge 1
    }
  }
  else if (dist_to_edge2 == 0.0) {
    if (dist_to_edge3 == 0.0) {
      edge_or_vertex = 6;  // vertex 3
    }
    else
    {
      edge_or_vertex = 2;  // edge 2
    }
  }
    else if (dist_to_edge3 == 0.0) {
    edge_or_vertex = 3;  // edge 3
  }
  return edge_or_vertex;
}
