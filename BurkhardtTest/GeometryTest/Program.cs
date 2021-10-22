using System;

namespace GeometryTest
{
class Program
{
static void Main(string[] args)
//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for GEOMETRY_PRB.
//
//  Discussion:
//
//    GEOMETRY_PRB tests the GEOMETRY library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//   20 January 2018
//
//  Author:
//
//    John Burkardt
//
{
Console.WriteLine("");
Console.WriteLine("GEOMETRY_PRB");
Console.WriteLine("  Test the GEOMETRY library.");

AngleTest.angle_box_2d_test ( );
AngleTest.angle_contains_ray_2d_test ( );
AngleTest.angle_deg_2d_test ( );
AngleTest.angle_half_2d_test ( );
AngleTest.angle_rad_2d_test ( );
AngleTest.angle_rad_3d_test ( );
AngleTest.angle_rad_nd_test ( );
AngleTest.angle_turn_2d_test ( );

AnnulusTest.annulus_sector_centroid_2d_test ( );

BallTest.ball01_sample_2d_test ( );
BallTest.ball01_sample_3d_test ( );
BallTest.ball01_sample_nd_test ( );
BallTest.ball01_volume_test ( ) ;

BasisTest.basis_map_3d_test ( );

BoxTest.box_contains_point_2d_test ( );

BoxTest.box_segment_clip_2d_test ( );
BoxTest.box_ray_int_2d_test ( );
BoxTest.box01_contains_point_2d_test ( );

CircleTest.circle_dia2imp_2d_test ( );

CircleTest.circle_exp_contains_point_2d_test ( );
CircleTest.circle_exp2imp_2d_test ( );

CircleTest.circle_imp_point_dist_2d_test ( );
CircleTest.circle_imp_points_arc_2d_test ( );

CircleTest.circle_llr2imp_2d_test ( );

CircleTest.circle_lune_angle_by_height_2d_test ( );
CircleTest.circle_lune_area_by_angle_2d_test ( );
CircleTest.circle_lune_area_by_height_2d_test ( );
CircleTest.circle_lune_centroid_2d_test ( );
CircleTest.circle_lune_height_by_angle_2d_test ( );

CircleTest.circle_pppr2imp_3d_test ( );

CircleTest.circle_ppr2imp_2d_test ( );

CircleTest.circle_sector_area_2d_test ( );
CircleTest.circle_sector_centroid_2d_test ( );

CircleTest.circle_triangle_area_2d_test ( );

CircleTest.test0155 ( );
CircleTest.test0156 ( );
CircleTest.test016 ( );
CircleTest.test0165 ( );

CircleTest.circles_intersect_area_2d_test ( );
CircleTest.circles_intersect_points_2d_test ( );

CubeTest.test020 ( );

CubeTest.cube01_volume_test ( );

CylinderTest.cylinder_point_dist_3d_test ( );
CylinderTest.cylinder_point_dist_signed_3d_test ( );

CylinderTest.test0202 ( );
CylinderTest.test0203 ( );
CylinderTest.test02035 ( );
CylinderTest.test0204 ( );
DegRadTest.test0205 ( );
DirectionTest.test021 ( );
DirectionTest.test022 ( );

DirectionTest.direction_uniform_nd_test ( );

DiskPointTest.disk_point_dist_3d_test ( );

test0234 ( );
test0235 ( );
test0236 ( );

DualTest.dual_size_3d_test ( );
DualTest.dual_shape_3d_test ( );


test025 ( );

ellipse_area1_test ( );
ellipse_area2_test ( );
ellipse_area3_test ( );
ellipse_point_near_2d_test ( );

test026 ( );

ellipse_points_arc_2d_test ( );

halfplane_contains_point_2d_test ( );

test029 ( );
test030 ( );
test031 ( );
test0315 ( );
test032 ( );
test0321 ( );
test0322 ( );
test0323 ( );
test0325 ( );
test0327 ( );

line_exp_perp_2d_test ( );
line_exp_point_dist_2d ( );

test0336 ( );
test0337 ( );
test034 ( );
test0345 ( );
test0346 ( );

line_imp_point_dist_2d_test ( );

test0351 ( );
test0352 ( );
test038 ( );
test0385 ( );
test03855 ( );
test0386 ( );

lines_exp_int_2d_test ( );

test040 ( );
test041 ( );
test0415 ( );
test0416 ( );
test0418 ( );
test042 ( );
test043 ( );
test044 ( );
test045 ( );
test046 ( );
test047 ( );
test0475 ( );
test0477 ( );
test0478 ( );

parallelogram_contains_point_2d_test ( );
parallelogram_contains_point_2d_test2 ( );
parallelogram_contains_point_3d_test ( );

test0493 ( );
test0495 ( );

plane_exp_normal_3d_test ( );

test051 ( );
test052 ( );
test053 ( );
test054 ( );

plane_imp2normal_3d_test ( );
plane_imp_line_par_int_3d_test ( );

test057 ( );
test058 ( );
test059 ( );
test060 ( );
test061 ( );
test0615 ( );
test0616 ( );
test0617 ( );
test062 ( );
test063 ( );
test064 ( );

points_centroid_2d_test ( );
points_colin_2d_test ( );

test068 ( );
test0685 ( );
test0755 ( );

polygon_angles_2d_test ( );

test076 ( );
test0765 ( );
test078 ( );
test0782 ( );
test0784 ( );

polygon_centroid_3d_test ( );
polygon_contains_point_2d_test ( );
polygon_contains_point_2d_2_test ( );
polygon_contains_point_2d_3_test ( );

test080 ( );
test0803 ( );
test0805 ( );

polygon_solid_angle_3d_test ( );

polyhedron_area_3d_test ( );
polyhedron_centroid_3d_test ( );

test0825 ( );
test083 ( );
test084 ( );

polyline_points_nd_test ( );

test0845 ( );
test0846 ( );

plane_exp_pro3_test ( );

test170 ( );
test171 ( );
test1712 ( );
test1715 ( );

r8_acos_test ( );
r8_asin_test ( );
r8_atan_test ( );

test0243 ( );
test0245 ( );
test173 ( );
test174 ( );
test1745 ( );
test1746 ( );
test1787 ( );
test1893 ( );
test036 ( );
test0365 ( );

segment_point_dist_3d_test ( );
segment_point_near_2d_test ( );
segment_point_near_3d_test ( );
segment_point_near_3d_test2 ( );

test1788 ( );
test1789 ( );
test179 ( );
test180 ( );
test1805 ( );
test0125 ( );
test0126 ( );
test0127 ( );

sphere_dia2imp_3d_test ( );

test182 ( );
test183 ( );
test1835 ( );
test1836 ( );
test187 ( );
test188 ( );
test189 ( );
test1892 ( );
test1895 ( );
test190 ( );
test191 ( );
test192 ( );
test193 ( );

sphere_unit_sample_nd_2_test ( );

test195 ( );
test1955 ( );
test196 ( );
test197 ( );
test198 ( );
test199 ( );
test200 ( );
test201 ( );
test202 ( );
test203 ( );
test2031 ( );
test2032 ( );
test20321 ( );
test20322 ( );

tetrahedron_lattice_layer_point_next_test ( );

test203225 ( );
test20323 ( );
test203232 ( );
test203233 ( );
test203234 ( );
test203235 ( );
test20324 ( );
test20325 ( );

tetrahedron_solid_angles_3d_test ( );

test2033 ( );
test204 ( );
test205 ( );

triangle_angles_2d_test ( );

test20605 ( );
test2061 ( );
test2062 ( );
test209 ( );
test20655 ( );
test2066 ( );
test2094 ( );
test2101 ( );
test21011 ( );
test2067 ( );
test21015 ( );

triangle_contains_line_exp_3d_test ( );
triangle_contains_line_par_3d_test ( );

test207 ( );
test2075 ( );
test208 ( );
test2102 ( );
test2070 ( );
test20701 ( );
test2104 ( );
test2105 ( );
test211 ( );
test2103 ( );
test2071 ( );
test20715 ( );

triangle_point_dist_3d_test ( );
triangle_point_near_2d_test ( );
triangle_quality_2d_test ( );

test212 ( );
test213 ( );

tube_2d_test ( );

vector_directions_nd_test ( );
vector_rotate_2d_test ( );
vector_rotate_3d_test ( );
vector_rotate_base_2d_test ( );
vector_separation_nd_test ( );

voxels_dist_l1_nd_test ( );
voxels_line_3d_test ( );
voxels_region_3d_test ( );
voxels_step_3d_test ( );

wedge01_volume_test ( );

Console.WriteLine("");
Console.WriteLine("GEOMETRY_PRB");
Console.WriteLine("  Normal end of execution.");
Console.WriteLine("");
}
}
}