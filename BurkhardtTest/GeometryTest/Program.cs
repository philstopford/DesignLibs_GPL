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

r8Test.test0234 ( );
dmsradTest.test0235 ( );
DodecahedronTest.test0236 ( );

DualTest.dual_size_3d_test ( );
DualTest.dual_shape_3d_test ( );


EllipseTest.test025 ( );

EllipseTest.ellipse_area1_test ( );
EllipseTest.ellipse_area2_test ( );
EllipseTest.ellipse_area3_test ( );
EllipseTest.ellipse_point_near_2d_test ( );

EllipseTest.test026 ( );

EllipseTest.ellipse_points_arc_2d_test ( );

HalfTest.halfplane_contains_point_2d_test ( );

HalfTest.test029 ( );
HalfTest.test030 ( );
HaversineTest.test031 ( );
HexagonTest.test0315 ( );
HexagonTest.test032 ( );
HexagonTest.test0321 ( );
i4Test.test0322 ( );
i4Test.test0323 ( );
IcosahedronTest.test0325 ( );
LinesTest.test0327 ( );

LinesTest.line_exp_perp_2d_test ( );
LinesTest.line_exp_point_dist_2d ( );

LinesTest.test0336 ( );
LinesTest.test0337 ( );
LinesTest.test034 ( );
LinesTest.test0345 ( );
LinesTest.test0346 ( );

LinesTest.line_imp_point_dist_2d_test ( );

LinesTest.test0351 ( );
LinesTest.test0352 ( );
LinesTest.test038 ( );
LinesTest.test0385 ( );
LinesTest.test03855 ( );
LinesTest.test0386 ( );

LinesTest.lines_exp_int_2d_test ( );

LinesTest.test040 ( );
LinesTest.test041 ( );
LinesTest.test0415 ( );
LinesTest.test0416 ( );
SegmentTest.test0418 ( );
SegmentTest.test042 ( );
SegmentTest.test043 ( );
SegmentTest.test044 ( );
SegmentTest.test045 ( );
LocalMinimumTest.test046 ( );
LocalMinimumTest.test047 ( );
OctahedronTest.test0475 ( );
ParallelogramTest.test0477 ( );
ParallelogramTest.test0478 ( );

ParallelogramTest.parallelogram_contains_point_2d_test ( );
ParallelogramTest.parallelogram_contains_point_2d_test2 ( );
ParallelogramTest.parallelogram_contains_point_3d_test ( );

ParabolaTest.test0493 ( );
ParallelepipedTest.test0495 ( );

PlaneTest.plane_exp_normal_3d_test ( );

PlaneTest.test051 ( );
PlaneTest.test052 ( );
PlaneTest.test053 ( );
PlaneTest.test054 ( );

PlaneTest.plane_imp2normal_3d_test ( );
PlaneTest.plane_imp_line_par_int_3d_test ( );

PlaneTest.test057 ( );
PlaneTest.test058 ( );
PlaneTest.test059 ( );
PlaneTest.test060 ( );
PlaneTest.test061 ( );
PlaneTest.test0615 ( );
PlaneTest.test0616 ( );
PlaneTest.test0617 ( );
PlaneTest.test062 ( );
PlaneTest.test063 ( );
PlaneTest.test064 ( );

PointsTest.points_centroid_2d_test ( );
PointsTest.points_colin_2d_test ( );

SphereTest.test068 ( );
PolarTest.test0685 ( );
PolygonTest.test0755 ( );

PolygonTest.polygon_angles_2d_test ( );

PolygonTest.test076 ( );
PolygonTest.test0765 ( );
PolygonTest.test078 ( );
PolygonTest.test0782 ( );
PolygonTest.test0784 ( );

PolygonTest.polygon_centroid_3d_test ( );
PolygonTest.polygon_contains_point_2d_test ( );
PolygonTest.polygon_contains_point_2d_2_test ( );
PolygonTest.polygon_contains_point_2d_3_test ( );

PolygonTest.test080 ( );
PolygonTest.test0803 ( );
PolygonTest.test0805 ( );

PolygonTest.polygon_solid_angle_3d_test ( );

PolyhedronTest.polyhedron_area_3d_test ( );
PolyhedronTest.polyhedron_centroid_3d_test ( );

PolyhedronTest.test0825 ( );
PolyhedronTest.test083 ( );
PolylineTest.test084 ( );

PolylineTest.polyline_points_nd_test ( );

PolyloopTest.test0845 ( );
PolyloopTest.test0846 ( );

PlaneTest.plane_exp_pro3_test ( );

ProvecTest.test170 ( );
QuadrilateralTest.test171 ( );
QuadrilateralTest.test1712 ( );
QuadrilateralTest.test1715 ( );

r8Test.r8_acos_test ( );
r8Test.r8_asin_test ( );
r8Test.r8_atan_test ( );

r8Test.test0243 ( );
r8Test.test0245 ( );
RadecTest.test173 ( );
RadecTest.test174 ( );
r8Test.test1745 ( );
r8Test.test1746 ( );
DGETest.test1787 ( );
XYTest.test1893 ( );
SegmentTest.test036 ( );
SegmentTest.test0365 ( );

SegmentTest.segment_point_dist_3d_test ( );
SegmentTest.segment_point_near_2d_test ( );
SegmentTest.segment_point_near_3d_test ( );
SegmentTest.segment_point_near_3d_test2 ( );

SimplexTest.test1788 ( );
SimplexTest.test1789 ( );
IcosahedronTest.test179 ( );
SortHeapTest.test180 ( );
SimplexTest.test1805 ( );
SphereTest.test0125 ( );
SphereTest.test0126 ( );
SphereTest.test0127 ( );

SphereTest.sphere_dia2imp_3d_test ( );

SphereTest.test182 ( );
SphereTest.test183 ( );
SphereTest.test1835 ( );
SphereTest.test1836 ( );
SphereTest.test187 ( );
SphereTest.test188 ( );
SphereTest.test189 ( );
SphereTest.test1892 ( );
SphereTest.test1895 ( );
SphereTest.test190 ( );
SphereTest.test191 ( );
SphereTest.test192 ( );
SphereTest.test193 ( );

SphereTest.sphere_unit_sample_nd_2_test ( );

SphereTest.test195 ( );
SphereTest.test1955 ( );
ShapeTest.test196 ( );
ShapeTest.test197 ( );
ShapeTest.test198 ( );
ShapeTest.test199 ( );
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