using System;
using Burkardt.IO;
using Burkardt.Pointset;
using Burkardt.Sampling;
using Burkardt.Stroud;
using Burkardt.Table;
using Burkardt.TriangulationNS;
using Burkardt.Types;

namespace QualityTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for QUALITY_TEST.
            //
            //  Discussion:
            //
            //    QUALITY_TEST tests QUALITY.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            string input_filename;
            int n;
            int dim_num;
            int ns;
            int seed_init;
            double[] z;

            Console.WriteLine("");
            Console.WriteLine("QUALITY_TEST");
            Console.WriteLine("  Test the QUALITY library.");

            test_cvt();

            test_halton();

            test_sphere();

            Console.WriteLine("");
            Console.WriteLine("QUALITY_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test_cvt()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST_CVT carries out tests of a pointset in the unit hypercube.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int flag;
            string input_filename;
            int n;
            int dim_num;
            int ns;
            int seed_init;
            int[] triangle_neighbor = null;
            int[] triangle_node = null;
            int triangle_num = 0;
            double[] z = null;

            Console.WriteLine("");
            Console.WriteLine("TEST_HALTON");
            Console.WriteLine("  Analyze a pointset in the unit hypercube.");
            Console.WriteLine("");
            Console.WriteLine("  We use a built-in sample routine.");

            ns = 100000;
            seed_init = 123456789;
            input_filename = "cvt_02_00100.txt";

            TableHeader h = typeMethods.dtable_header_read(input_filename);
            dim_num = h.m;
            n = h.n;

            Console.WriteLine("");
            Console.WriteLine("  Measures of uniform point dispersion.");
            Console.WriteLine("");
            Console.WriteLine("  The pointset was read from \"" + input_filename + "\"");
            Console.WriteLine("  The sample routine will be SAMPLE_HYPERCUBE_UNIFORM.");
            Console.WriteLine("");
            Console.WriteLine("  Spatial dimension DIM_NUM =      " + dim_num + "");
            Console.WriteLine("  The number of points N =         " + n + "");
            Console.WriteLine("  The number of sample points NS = " + ns + "");
            Console.WriteLine("  The random number SEED_INIT =    " + seed_init + "");
            Console.WriteLine("");

            z = typeMethods.dtable_data_read(input_filename, dim_num, n);

            typeMethods.r8mat_transpose_print_some(dim_num, n, z, 1, 1, 5, 5,
                "  5x5 portion of data read from file:");
            //
            //  For 2D datasets, compute the Delaunay triangulation.
            //
            if (dim_num == 2)
            {
                triangle_node = new int[3 * 2 * n];
                triangle_neighbor = new int[3 * 2 * n];

                flag = Delauney.dtris2(n, 0, ref z, ref triangle_num, ref triangle_node, ref triangle_neighbor);
                Console.WriteLine("");
                Console.WriteLine("  Triangulated data generates " + triangle_num + " triangles.");
            }
            else
            {
                triangle_num = 0;
            }

            if (dim_num == 2)
            {
                test005(n, z, triangle_num, triangle_node);
                test006(n, z, triangle_num, triangle_node);
            }

            test007(dim_num, n, z);
            test01(dim_num, n, z, ns, Burkardt.HyperGeometry.Hypercube.Sample.sample_hypercube_uniform, seed_init);
            test02(dim_num, n, z, ns, Burkardt.HyperGeometry.Hypercube.Sample.sample_hypercube_uniform, seed_init);
            test03(dim_num, n, z, ns, Burkardt.HyperGeometry.Hypercube.Sample.sample_hypercube_uniform, seed_init);
            test04(dim_num, n, z);
            test05(dim_num, n, z, ns, Burkardt.HyperGeometry.Hypercube.Sample.sample_hypercube_uniform, seed_init);
            test06(dim_num, n, z);
            test07(dim_num, n, z, ns, Burkardt.HyperGeometry.Hypercube.Sample.sample_hypercube_uniform, seed_init);
            test08(dim_num, n, z, ns, Burkardt.HyperGeometry.Hypercube.Sample.sample_hypercube_uniform, seed_init);

            if (dim_num == 2)
            {
                test083(n, z, triangle_num, triangle_node);
            }

            test085(dim_num, n, z);
            test09(dim_num, n, z);
            test10(dim_num, n, z, ns, Burkardt.HyperGeometry.Hypercube.Sample.sample_hypercube_uniform, seed_init);
            test11(dim_num, n, z);

        }

        static void test_halton()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST_HALTON carries out tests of a pointset in the unit hypercube.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int flag;
            string input_filename;
            int n;
            int dim_num;
            int ns;
            int seed_init;
            int[] triangle_neighbor = null;
            int[] triangle_node = null;
            int triangle_num = 0;
            double[] z = null;

            Console.WriteLine("");
            Console.WriteLine("TEST_HALTON");
            Console.WriteLine("  Analyze a pointset in the unit hypercube.");
            Console.WriteLine("");
            Console.WriteLine("  We use a built-in sample routine.");

            ns = 100000;
            seed_init = 123456789;
            input_filename = "halton_02_00100.txt";

            TableHeader h = typeMethods.dtable_header_read(input_filename);
            dim_num = h.m;
            n = h.n;

            Console.WriteLine("");
            Console.WriteLine("  Measures of uniform point dispersion.");
            Console.WriteLine("");
            Console.WriteLine("  The pointset was read from \"" + input_filename + "\"");
            Console.WriteLine("  The sample routine will be SAMPLE_HYPERCUBE_UNIFORM.");
            Console.WriteLine("");
            Console.WriteLine("  Spatial dimension DIM_NUM =      " + dim_num + "");
            Console.WriteLine("  The number of points N =         " + n + "");
            Console.WriteLine("  The number of sample points NS = " + ns + "");
            Console.WriteLine("  The random number SEED_INIT =    " + seed_init + "");
            Console.WriteLine("");

            z = typeMethods.dtable_data_read(input_filename, dim_num, n);

            typeMethods.r8mat_transpose_print_some(dim_num, n, z, 1, 1, 5, 5,
                "  5x5 portion of data read from file:");
            //
            //  For 2D datasets, compute the Delaunay triangulation.
            //
            if (dim_num == 2)
            {
                triangle_node = new int[3 * 2 * n];
                triangle_neighbor = new int[3 * 2 * n];

                flag = Delauney.dtris2(n, 0, ref z, ref triangle_num, ref triangle_node, ref triangle_neighbor);
                Console.WriteLine("");
                Console.WriteLine("  Triangulated data generates " + triangle_num + " triangles.");
            }
            else
            {
                triangle_num = 0;
            }

            if (dim_num == 2)
            {
                test005(n, z, triangle_num, triangle_node);
                test005(n, z, triangle_num, triangle_node);
            }

            test007(dim_num, n, z);
            test01(dim_num, n, z, ns, Burkardt.HyperGeometry.Hypercube.Sample.sample_hypercube_uniform, seed_init);
            test02(dim_num, n, z, ns, Burkardt.HyperGeometry.Hypercube.Sample.sample_hypercube_uniform, seed_init);
            test03(dim_num, n, z, ns, Burkardt.HyperGeometry.Hypercube.Sample.sample_hypercube_uniform, seed_init);
            test04(dim_num, n, z);
            test05(dim_num, n, z, ns, Burkardt.HyperGeometry.Hypercube.Sample.sample_hypercube_uniform, seed_init);
            test06(dim_num, n, z);
            test07(dim_num, n, z, ns, Burkardt.HyperGeometry.Hypercube.Sample.sample_hypercube_uniform, seed_init);
            test08(dim_num, n, z, ns, Burkardt.HyperGeometry.Hypercube.Sample.sample_hypercube_uniform, seed_init);

            if (dim_num == 2)
            {
                test083(n, z, triangle_num, triangle_node);
            }

            test085(dim_num, n, z);
            test09(dim_num, n, z);
            test10(dim_num, n, z, ns, Burkardt.HyperGeometry.Hypercube.Sample.sample_hypercube_uniform, seed_init);
            test11(dim_num, n, z);

        }

        static void test_sphere()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST_SPHERE carries out tests of a pointset in the unit sphere.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int flag;
            string input_filename;
            int n;
            int dim_num;
            int ns;
            int seed_init;
            int[] triangle_neighbor = null;
            int[] triangle_node = null;
            int triangle_num = 0;
            double[] z = null;

            Console.WriteLine("");
            Console.WriteLine("TEST_SPHERE");
            Console.WriteLine("  Analyze a pointset in the unit sphere.");
            Console.WriteLine("");
            Console.WriteLine("  We use a built-in sample routine.");

            ns = 100000;
            seed_init = 123456789;
            input_filename = "sphere_02_00100.txt";

            TableHeader h = typeMethods.dtable_header_read(input_filename);
            dim_num = h.m;
            n = h.n;

            Console.WriteLine("");
            Console.WriteLine("  Measures of uniform point dispersion.");
            Console.WriteLine("");
            Console.WriteLine("  The pointset was read from \"" + input_filename + "\"");
            Console.WriteLine("  The sample routine will be SAMPLE_SPHERE_UNIFORM.");
            Console.WriteLine("");
            Console.WriteLine("  The spatial dimension DIM_NUM =     " + dim_num + "");
            Console.WriteLine("  The number of points N =         " + n + "");
            Console.WriteLine("  The number of sample points NS = " + ns + "");
            Console.WriteLine("  The random number SEED_INIT =    " + seed_init + "");
            Console.WriteLine("");

            z = typeMethods.dtable_data_read(input_filename, dim_num, n);

            typeMethods.r8mat_transpose_print_some(dim_num, n, z, 1, 1, 5, 5,
                "  5x5 portion of data read from file:");
            //
            //  For 2D datasets, compute the Delaunay triangulation.
            //
            if (dim_num == 2)
            {
                triangle_node = new int[3 * 2 * n];
                triangle_neighbor = new int[3 * 2 * n];

                flag = Delauney.dtris2(n, 0, ref z, ref triangle_num, ref triangle_node, ref triangle_neighbor);
                Console.WriteLine("");
                Console.WriteLine("  Triangulated data generates " + triangle_num + " triangles.");
            }
            else
            {
                triangle_num = 0;
            }

            if (dim_num == 2)
            {
                test005(n, z, triangle_num, triangle_node);
                test006(n, z, triangle_num, triangle_node);
            }

            test007(dim_num, n, z);
            test01(dim_num, n, z, ns, Burkardt.SphereNS.Sample.sample_sphere_uniform, seed_init);
            test02(dim_num, n, z, ns, Burkardt.SphereNS.Sample.sample_sphere_uniform, seed_init);
            test03(dim_num, n, z, ns, Burkardt.SphereNS.Sample.sample_sphere_uniform, seed_init);
            test04(dim_num, n, z);
            test05(dim_num, n, z, ns, Burkardt.SphereNS.Sample.sample_sphere_uniform, seed_init);
            test06(dim_num, n, z);
            test07(dim_num, n, z, ns, Burkardt.SphereNS.Sample.sample_sphere_uniform, seed_init);
            test08(dim_num, n, z, ns, Burkardt.SphereNS.Sample.sample_sphere_uniform, seed_init);

            if (dim_num == 2)
            {
                test083(n, z, triangle_num, triangle_node);
            }

            test085(dim_num, n, z);
            test09(dim_num, n, z);
            test10(dim_num, n, z, ns, Burkardt.SphereNS.Sample.sample_sphere_uniform, seed_init);
            test11(dim_num, n, z);

        }

        static void test005(int n, double[] z, int nt, int[] triangle)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST005 tests ALPHA_MEASURE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    08 November 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int triangle_order = 3;

            Console.WriteLine("");
            Console.WriteLine("TEST005");
            Console.WriteLine("  ALPHA_MEASURE computes the ALPHA measure of quality.");
            Console.WriteLine("  The minimum angle measure    ALPHA = "
                              + Quality.alpha_measure(n, z, triangle_order, nt, triangle) + "");

            return;
        }

        static void test006(int n, double[] z, int nt, int[] triangle)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST006 tests AREA_MEASURE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    08 November 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int triangle_order = 3;

            Console.WriteLine("");
            Console.WriteLine("TEST006");
            Console.WriteLine("  AREA_MEASURE computes the AREA measure of quality.");
            Console.WriteLine("  The area ratio measure        AREA = "
                              + Quality.area_measure(n, z, triangle_order, nt, triangle) + "");

            return;
        }

        static void test007(int dim_num, int n, double[] z)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST007 tests BETA_MEASURE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("TEST007");
            Console.WriteLine("  BETA_MEASURE computes the BETA measure of quality.");
            Console.WriteLine("  Relative spacing deviation BETA =    "
                              + Quality.beta_measure(dim_num, n, z) + "");

            return;
        }

        static void test01(int dim_num, int n, double[] z, int ns,
                Func<int, int, int, GeometrySampleResult> sample_routine, int seed_init)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 tests CHI_MEASURE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  CHI_MEASURE computes the CHI measure of quality.");
            Console.WriteLine("  The regularity measure         Chi = "
                              + Quality.chi_measure(dim_num, n, z, ns, sample_routine, seed_init) + "");

            return;
        }

        static void test02(int dim_num, int n, double[] z, int ns,
                Func<int, int, int, GeometrySampleResult> sample_routine, int seed_init)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 tests D_MEASURE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("TEST02");
            Console.WriteLine("  D_MEASURE computes the D measure of quality.");
            Console.WriteLine("  2nd moment determinant measure   D = "
                              + Quality.d_measure(dim_num, n, z, ns, sample_routine, seed_init) + "");

            return;
        }

        static void test03(int dim_num, int n, double[] z, int ns,
                Func<int, int, int, GeometrySampleResult> sample_routine, int seed_init)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST03 tests E_MEASURE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("TEST03");
            Console.WriteLine("  E_MEASURE computes the E measure of quality.");
            Console.WriteLine("  The Voronoi energy measure       E = "
                              + Quality.e_measure(dim_num, n, z, ns, sample_routine, seed_init) + "");

            return;
        }

        static void test04(int dim_num, int n, double[] z)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST04 tests GAMMA_MEASURE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("TEST04");
            Console.WriteLine("  GAMMA_MEASURE computes the Gamma measure of quality.");
            Console.WriteLine("  The mesh ratio               Gamma = "
                              + Quality.gamma_measure(dim_num, n, z) + "");

            return;
        }

        static void test05(int dim_num, int n, double[] z, int ns,
                Func<int, int, int, GeometrySampleResult> sample_routine, int seed_init)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST05 tests H_MEASURE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("TEST05");
            Console.WriteLine("  H_MEASURE computes the H measure of quality.");
            Console.WriteLine("  The point distribution norm      H = "
                              + Quality.h_measure(dim_num, n, z, ns, sample_routine, seed_init) + "");

            return;
        }

        static void test06(int dim_num, int n, double[] z)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST06 tests LAMBDA_MEASURE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("TEST06");
            Console.WriteLine("  LAMBDA_MEASURE computes the Lambda measure of quality.");
            Console.WriteLine("  The COV measure             Lambda = "
                              + Quality.lambda_measure(dim_num, n, z) + "");

            return;
        }

        static void test07(int dim_num, int n, double[] z, int ns,
                Func<int, int, int, GeometrySampleResult> sample_routine, int seed_init)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST07 tests MU_MEASURE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("TEST07");
            Console.WriteLine("  MU_MEASURE computes the Mu measure of quality.");
            Console.WriteLine("  The point distribution ratio    Mu = "
                              + Quality.mu_measure(dim_num, n, z, ns, sample_routine, seed_init) + "");

            return;
        }

        static void test08(int dim_num, int n, double[] z, int ns,
                Func<int, int, int, GeometrySampleResult> sample_routine, int seed_init)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST08 tests NU_MEASURE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("TEST08");
            Console.WriteLine("  NU_MEASURE computes the Nu measure of quality.");
            Console.WriteLine("  The cell volume deviation       Nu = "
                              + Quality.nu_measure(dim_num, n, z, ns, sample_routine, seed_init) + "");

            return;
        }

        static void test083(int n, double[] z, int nt, int[] triangle)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST083 tests Q_MEASURE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    08 November 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int triangle_order = 3;

            Console.WriteLine("");
            Console.WriteLine("TEST083");
            Console.WriteLine("  Q_MEASURE computes the Q measure of quality.");
            Console.WriteLine("  The triangle uniformity measure  Q = "
                              + Quality.q_measure(n, z, triangle_order, nt, triangle) + "");

            return;
        }

        static void test085(int dim_num, int n, double[] z)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST085 tests R0_MEASURE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("TEST085");
            Console.WriteLine("  R0_MEASURE computes the Riesz S = 0 energy measure of quality.");
            Console.WriteLine("  The R0 measure                  R0 = "
                              + Quality.r0_measure(dim_num, n, z) + "");

            return;
        }

        static void test09(int dim_num, int n, double[] z)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST09 tests SPHERE_MEASURE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("TEST09");
            Console.WriteLine("  SPHERE_MEASURE computes the sphere measure of quality.");
            Console.WriteLine("  Nonintersecting sphere volume    S = "
                              + Quality.sphere_measure(dim_num, n, z) + "");

            return;
        }

        static void test10(int dim_num, int n, double[] z, int ns,
                Func<int, int, int, GeometrySampleResult> sample_routine, int seed_init)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST10 tests TAU_MEASURE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("TEST10");
            Console.WriteLine("  TAU_MEASURE computes the Tau measure of quality.");
            Console.WriteLine("  2nd moment trace measure       Tau = "
                              + Quality.tau_measure(dim_num, n, z, ns, sample_routine, seed_init) + "");

            return;
        }

        static void test11(int dim_num, int n, double[] z)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST11 tests POINTSET_SPACING.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 July 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] gamma;
            double gamma_ave;
            double gamma_max;
            double gamma_min;
            double gamma_std;
            int i;

            Console.WriteLine("");
            Console.WriteLine("TEST11");
            Console.WriteLine("  POINTSET_SPACING computes pointset spacing parameters.");

            gamma = Spacing.pointset_spacing(dim_num, n, z);

            gamma_min = typeMethods.r8vec_min(n, gamma);
            gamma_max = typeMethods.r8vec_max(n, gamma);

            gamma_ave = 0.0;
            for (i = 0; i < n; i++)
            {
                gamma_ave = gamma_ave + gamma[i];
            }

            gamma_ave = gamma_ave / (double) (n);

            if (1 < n)
            {
                gamma_std = 0.0;
                for (i = 0; i < n; i++)
                {
                    gamma_std = gamma_std + Math.Pow(gamma[i] - gamma_ave, 2);
                }

                gamma_std = Math.Sqrt(gamma_std / (double) (n - 1));
            }
            else
            {
                gamma_std = 0.0;
            }

            Console.WriteLine("");
            Console.WriteLine("  Minimum spacing          GAMMA_MIN = " + gamma_min + "");
            Console.WriteLine("  Average spacing          GAMMA_AVE = " + gamma_ave + "");
            Console.WriteLine("  Maximum spacing          GAMMA_MAX = " + gamma_max + "");
            Console.WriteLine("  Spacing standard dev     GAMMA_STD = " + gamma_std + "");

        }

        public static double sphere_measure(int dim_num, int n, double[] z)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_MEASURE determines the pointset quality measure S.
            //
            //  Discussion:
            //
            //    This routine computes a measure of even spacing for a set of N points
            //    in the DIM_NUM-dimensional unit hypercube.  We will discuss the program
            //    as though the space is 2-dimensional and the spheres are circles, but
            //    the program may be used for general DIM_NUM-dimensional data.
            //
            //    The points are assumed to lie in the unit square.
            //
            //    The program makes a circle-packing measurement on the points
            //    by assuming that, at each point, a circle is centered; all
            //    the circles start out with zero radius, and then expand
            //    together at the same rate.  A circle stops expanding as soon
            //    as it touches any other circle.
            //
            //    The amount of area covered by the circles is compared to the
            //    area of the unit square.  This measurement has a certain amount
            //    of boundary effect: some circles will naturally extend outside
            //    the unit hypercube.  If this is a concern, is possible to restrict
            //    the circles to remain inside the unit hypercube.  In any case,
            //    this problem generally goes away as the number of points increases.
            //
            //    Since we are interested in the coverage of the unit hypercube,
            //    it is probably best if the circles are restricted.  This way,
            //    computing the area of the circles gives a measure of the even
            //    coverage of the region, relative to the presumably best possible
            //    covering, by the same number of circles, but of equal radius.
            //
            //    In the limit, the maximum relative packing density of a 2D
            //    region with equal-sized circles is 0.9069.  In 3D, a density
            //    of at least 0.74 can be achieved, and it is known that no
            //    greater than 0.7796 is possible.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 October 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, double Z[DIM_NUM*N], the points.
            //
            //    Output, double SPHERE_MEASURE, the amount of volume taken up
            //    by the nonintersecting spheres of maximum radius around each
            //    point.  Ignoring boundary effects, the "ideal" value would be
            //    1 (achievable only in 1 dimension), and the maximum value
            //    possible is the sphere packing density in the given spatial
            //    dimension.  If boundary effects can be ignored, the value of
            //    SPHERE_VOLUME reports how closely the given set of points
            //    behaves like a set of close-packed spheres.
            //
            //  Local Parameters:
            //
            //    Local, logical WALLS, is TRUE if the spheres are restricted
            //    to lie within the unit hypercube.
            //
        {
            int i;
            int j;
            double[] radius;
            double radius_ave;
            double radius_max;
            double radius_min;
            double sphere;
            bool verbose = false;
            double volume;
            bool walls = true;

            if (!typeMethods.r8mat_in_01(dim_num, n, z))
            {
                Console.WriteLine("");
                Console.WriteLine("SPHERE_MEASURE - Fatal error!");
                Console.WriteLine("  Some of the data is not inside the unit hypercube.");
                return typeMethods.r8_huge();
            }

            radius = Quality.radius_maximus(dim_num, n, z, walls);

            sphere = 0.0;
            for (i = 0; i < n; i++)
            {
                volume = Sphere.sphere_volume_nd(dim_num, radius[i]);
                sphere = sphere + volume;
            }

            if (verbose)
            {
                radius_ave = 0.0;
                radius_min = typeMethods.r8_huge();
                radius_max = 0.0;
                for (j = 0; j < n; j++)
                {
                    radius_ave = radius_ave + radius[j];
                    radius_min = Math.Min(radius_min, radius[j]);
                    radius_max = Math.Max(radius_max, radius[j]);
                }

                Console.WriteLine("");
                Console.WriteLine("  Number of dimensions is " + dim_num + "");
                Console.WriteLine("  Number of points is " + n + "");
                if (walls)
                {
                    Console.WriteLine("  Spheres are required to stay in the unit hypercube.");
                }
                else
                {
                    Console.WriteLine("  Spheres are NOT required to stay in the unit hypercube.");
                }

                Console.WriteLine("");
                Console.WriteLine("  Average radius = " + radius_ave + "");
                Console.WriteLine("  Minimum radius = " + radius_min + "");
                Console.WriteLine("  Maximum radius = " + radius_max + "");
                Console.WriteLine("  Sphere volume =  " + sphere + "");
            }

            return sphere;
        }


    }
}