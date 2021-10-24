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

            AngleTest.angle_box_2d_test();
            AngleTest.angle_contains_ray_2d_test();
            AngleTest.angle_deg_2d_test();
            AngleTest.angle_half_2d_test();
            AngleTest.angle_rad_2d_test();
            AngleTest.angle_rad_3d_test();
            AngleTest.angle_rad_nd_test();
            AngleTest.angle_turn_2d_test();

            AnnulusTest.annulus_sector_centroid_2d_test();

            BallTest.ball01_sample_2d_test();
            BallTest.ball01_sample_3d_test();
            BallTest.ball01_sample_nd_test();
            BallTest.ball01_volume_test();

            BasisTest.basis_map_3d_test();

            BoxTest.box_contains_point_2d_test();

            BoxTest.box_segment_clip_2d_test();
            BoxTest.box_ray_int_2d_test();
            BoxTest.box01_contains_point_2d_test();

            CircleTest.circle_dia2imp_2d_test();

            CircleTest.circle_exp_contains_point_2d_test();
            CircleTest.circle_exp2imp_2d_test();

            CircleTest.circle_imp_point_dist_2d_test();
            CircleTest.circle_imp_points_arc_2d_test();

            CircleTest.circle_llr2imp_2d_test();

            CircleTest.circle_lune_angle_by_height_2d_test();
            CircleTest.circle_lune_area_by_angle_2d_test();
            CircleTest.circle_lune_area_by_height_2d_test();
            CircleTest.circle_lune_centroid_2d_test();
            CircleTest.circle_lune_height_by_angle_2d_test();

            CircleTest.circle_pppr2imp_3d_test();

            CircleTest.circle_ppr2imp_2d_test();

            CircleTest.circle_sector_area_2d_test();
            CircleTest.circle_sector_centroid_2d_test();

            CircleTest.circle_triangle_area_2d_test();

            CircleTest.test0155();
            CircleTest.test0156();
            CircleTest.test016();
            CircleTest.test0165();

            CircleTest.circles_intersect_area_2d_test();
            CircleTest.circles_intersect_points_2d_test();

            CubeTest.test020();

            CubeTest.cube01_volume_test();

            CylinderTest.cylinder_point_dist_3d_test();
            CylinderTest.cylinder_point_dist_signed_3d_test();

            CylinderTest.test0202();
            CylinderTest.test0203();
            CylinderTest.test02035();
            CylinderTest.test0204();
            DegRadTest.test0205();
            DirectionTest.test021();
            DirectionTest.test022();

            DirectionTest.direction_uniform_nd_test();

            DiskPointTest.disk_point_dist_3d_test();

            r8Test.test0234();
            dmsradTest.test0235();
            DodecahedronTest.test0236();

            DualTest.dual_size_3d_test();
            DualTest.dual_shape_3d_test();


            EllipseTest.test025();

            EllipseTest.ellipse_area1_test();
            EllipseTest.ellipse_area2_test();
            EllipseTest.ellipse_area3_test();
            EllipseTest.ellipse_point_near_2d_test();

            EllipseTest.test026();

            EllipseTest.ellipse_points_arc_2d_test();

            HalfTest.halfplane_contains_point_2d_test();

            HalfTest.test029();
            HalfTest.test030();
            HaversineTest.test031();
            HexagonTest.test0315();
            HexagonTest.test032();
            HexagonTest.test0321();
            i4Test.test0322();
            i4Test.test0323();
            IcosahedronTest.test0325();
            LinesTest.test0327();

            LinesTest.line_exp_perp_2d_test();
            LinesTest.line_exp_point_dist_2d();

            LinesTest.test0336();
            LinesTest.test0337();
            LinesTest.test034();
            LinesTest.test0345();
            LinesTest.test0346();

            LinesTest.line_imp_point_dist_2d_test();

            LinesTest.test0351();
            LinesTest.test0352();
            LinesTest.test038();
            LinesTest.test0385();
            LinesTest.test03855();
            LinesTest.test0386();

            LinesTest.lines_exp_int_2d_test();

            LinesTest.test040();
            LinesTest.test041();
            LinesTest.test0415();
            LinesTest.test0416();
            SegmentTest.test0418();
            SegmentTest.test042();
            SegmentTest.test043();
            SegmentTest.test044();
            SegmentTest.test045();
            LocalMinimumTest.test046();
            LocalMinimumTest.test047();
            OctahedronTest.test0475();
            ParallelogramTest.test0477();
            ParallelogramTest.test0478();

            ParallelogramTest.parallelogram_contains_point_2d_test();
            ParallelogramTest.parallelogram_contains_point_2d_test2();
            ParallelogramTest.parallelogram_contains_point_3d_test();

            ParabolaTest.test0493();
            ParallelepipedTest.test0495();

            PlaneTest.plane_exp_normal_3d_test();

            PlaneTest.test051();
            PlaneTest.test052();
            PlaneTest.test053();
            PlaneTest.test054();

            PlaneTest.plane_imp2normal_3d_test();
            PlaneTest.plane_imp_line_par_int_3d_test();

            PlaneTest.test057();
            PlaneTest.test058();
            PlaneTest.test059();
            PlaneTest.test060();
            PlaneTest.test061();
            PlaneTest.test0615();
            PlaneTest.test0616();
            PlaneTest.test0617();
            PlaneTest.test062();
            PlaneTest.test063();
            PlaneTest.test064();

            PointsTest.points_centroid_2d_test();
            PointsTest.points_colin_2d_test();

            SphereTest.test068();
            PolarTest.test0685();
            PolygonTest.test0755();

            PolygonTest.polygon_angles_2d_test();

            PolygonTest.test076();
            PolygonTest.test0765();
            PolygonTest.test078();
            PolygonTest.test0782();
            PolygonTest.test0784();

            PolygonTest.polygon_centroid_3d_test();
            PolygonTest.polygon_contains_point_2d_test();
            PolygonTest.polygon_contains_point_2d_2_test();
            PolygonTest.polygon_contains_point_2d_3_test();

            PolygonTest.test080();
            PolygonTest.test0803();
            PolygonTest.test0805();

            PolygonTest.polygon_solid_angle_3d_test();

            PolyhedronTest.polyhedron_area_3d_test();
            PolyhedronTest.polyhedron_centroid_3d_test();

            PolyhedronTest.test0825();
            PolyhedronTest.test083();
            PolylineTest.test084();

            PolylineTest.polyline_points_nd_test();

            PolyloopTest.test0845();
            PolyloopTest.test0846();

            PlaneTest.plane_exp_pro3_test();

            ProvecTest.test170();
            QuadrilateralTest.test171();
            QuadrilateralTest.test1712();
            QuadrilateralTest.test1715();

            r8Test.r8_acos_test();
            r8Test.r8_asin_test();
            r8Test.r8_atan_test();

            r8Test.test0243();
            r8Test.test0245();
            RadecTest.test173();
            RadecTest.test174();
            r8Test.test1745();
            r8Test.test1746();
            DGETest.test1787();
            XYTest.test1893();
            SegmentTest.test036();
            SegmentTest.test0365();

            SegmentTest.segment_point_dist_3d_test();
            SegmentTest.segment_point_near_2d_test();
            SegmentTest.segment_point_near_3d_test();
            SegmentTest.segment_point_near_3d_test2();

            SimplexTest.test1788();
            SimplexTest.test1789();
            IcosahedronTest.test179();
            SortHeapTest.test180();
            SimplexTest.test1805();
            SphereTest.test0125();
            SphereTest.test0126();
            SphereTest.test0127();

            SphereTest.sphere_dia2imp_3d_test();

            SphereTest.test182();
            SphereTest.test183();
            SphereTest.test1835();
            SphereTest.test1836();
            SphereTest.test187();
            SphereTest.test188();
            SphereTest.test189();
            SphereTest.test1892();
            SphereTest.test1895();
            SphereTest.test190();
            SphereTest.test191();
            SphereTest.test192();
            SphereTest.test193();

            SphereTest.sphere_unit_sample_nd_2_test();

            SphereTest.test195();
            SphereTest.test1955();
            ShapeTest.test196();
            ShapeTest.test197();
            ShapeTest.test198();
            ShapeTest.test199();
            SphereTest.test200();
            SegmentTest.test201();
            EllipseTest.test202();
            TetrahedronTest.test203();
            TetrahedronTest.test2031();
            TetrahedronTest.test2032();
            TetrahedronTest.test20321();
            TetrahedronTest.test20322();

            TetrahedronTest.tetrahedron_lattice_layer_point_next_test();

            TetrahedronTest.test203225();
            TetrahedronTest.test20323();
            TetrahedronTest.test203232();
            TetrahedronTest.test203233();
            TetrahedronTest.test203234();
            TetrahedronTest.test203235();
            TetrahedronTest.test20324();
            TetrahedronTest.test20325();

            TetrahedronTest.tetrahedron_solid_angles_3d_test();

            TetrahedronTest.test2033();
            TransMatTest.test204();
            TransMatTest.test205();

            TriangleTest.triangle_angles_2d_test();

            TriangleTest.test20605();
            TriangleTest.test2061();
            TriangleTest.test2062();
            TriangleTest.test209();
            TriangleTest.test20655();
            TriangleTest.test2066();
            TriangleTest.test2094();
            TriangleTest.test2101();
            TriangleTest.test21011();
            TriangleTest.test2067();
            TriangleTest.test21015();

            TriangleTest.triangle_contains_line_exp_3d_test();
            TriangleTest.triangle_contains_line_par_3d_test();

            TriangleTest.test207();
            TriangleTest.test2075();
            TriangleTest.test208();
            TriangleTest.test2102();
            TriangleTest.test2070();
            TriangleTest.test20701();
            TriangleTest.test2104();
            TriangleTest.test2105();
            TriangleTest.test211();
            TriangleTest.test2103();
            TriangleTest.test2071();
            TriangleTest.test20715();

            TriangleTest.triangle_point_dist_3d_test();
            TriangleTest.triangle_point_near_2d_test();
            TriangleTest.triangle_quality_2d_test();

            TriangleTest.test212();
            TriangleTest.test213();

            TubeTest.tube_2d_test();

            VectorTest.vector_directions_nd_test();
            VectorTest.vector_rotate_2d_test();
            VectorTest.vector_rotate_3d_test();
            VectorTest.vector_rotate_base_2d_test();
            VectorTest.vector_separation_nd_test();

            VoxelTest.voxels_dist_l1_nd_test();
            VoxelTest.voxels_line_3d_test();
            VoxelTest.voxels_region_3d_test();
            VoxelTest.voxels_step_3d_test();

            WedgeTest.wedge01_volume_test();

            Console.WriteLine("");
            Console.WriteLine("GEOMETRY_PRB");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }
    }
}