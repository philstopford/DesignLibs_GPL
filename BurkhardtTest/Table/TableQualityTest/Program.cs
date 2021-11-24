using System;
using Burkardt.Pointset;
using Burkardt.Sampling;
using Burkardt.Table;
using Burkardt.TriangulationNS;
using Burkardt.Types;

namespace TableQualityTest;

internal static class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TABLE_QUALITY.
        //
        //  Discussion:
        //
        //    TABLE_QUALITY determines quality measures for a given set of points.
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
        //  Reference:
        //
        //    Max Gunzburger, John Burkardt,
        //    Uniformity Measures for Point Samples in Hypercubes.
        //
        //  Local parameters:
        //
        //    Local, int DIM_NUM, the spatial dimension of the point set.
        //
        //    Local, int N, the number of points.
        //
        //    Local, double Z[DIM_NUM*N], the point set.
        //
        //    Local, int NS, the number of sample points.
        //
    {
        Console.WriteLine("");

        Console.WriteLine("");
        Console.WriteLine("TABLE_QUALITY");
        Console.WriteLine("  Compute measures of uniform dispersion for a pointset.");
        Console.WriteLine("");
        Console.WriteLine("  Note: this routine assumes that the input pointset");
        Console.WriteLine("  is contained in the unit hypercube.  If this is not");
        Console.WriteLine("  the case, then you must rewrite the routine");
        Console.WriteLine("    SAMPLE_ROUTINE");
        Console.WriteLine("  so that it properly returns sample points in your");
        Console.WriteLine("  region, with a uniform density, or a probability");
        Console.WriteLine("  density of your choosing.");
        //
        //  If the input file was not specified, get it now.
        //
        try
        {
            int i;
            for (i = 1; i < args.Length; i++)
            {
                handle(args[i], Burkardt.HyperGeometry.Hypercube.Sample.sample_hypercube_uniform);
            }
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("TABLE_QUALITY:");
            Console.WriteLine("  Please enter the name of a file to be analyzed.");

            string input_filename = Console.ReadLine();

            handle(input_filename, Burkardt.HyperGeometry.Hypercube.Sample.sample_hypercube_uniform);

        }

        Console.WriteLine("");
        Console.WriteLine("TABLE_QUALITY:");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }


    private static void handle(string input_filename,
            Func<int, int, int, GeometrySampleResult> sample_routine)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HANDLE handles a single file.
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
        //  Reference:
        //
        //    Max Gunzburger, John Burkardt,
        //    Uniformity Measures for Point Samples in Hypercubes.
        //
        //  Parameters:
        //
        //    Input, char *INPUT_FILENAME, the name of the input file.
        //
        //  Local parameters:
        //
        //    Local, int DIM_NUM, the spatial dimension of the point set.
        //
        //    Local, int N, the number of points.
        //
        //    Local, double Z[DIM_NUM*N], the point set.
        //
        //    Local, int NS, the number of sample points.
        //
        //    Input, double *SAMPLE_ROUTINE, the name of a routine which
        //    is used to produce an DIM_NUM by N array of sample points in the region, 
        //    of the form:
        //      double *sample_routine ( int dim_num, int n, int *seed )
        //
    {
        double gamma_std;
        int i;
        const int ns = 100000;
        int nt = 0;
        const int seed_init = 123456789;
        int[] triangle = null;
        int triangle_order;

        TableHeader h = typeMethods.dtable_header_read(input_filename);
        int dim_num = h.m;
        int n = h.n;

        // 
        //  Read the point set.
        //
        double[] z = typeMethods.dtable_data_read(input_filename, dim_num, n);
        switch (dim_num)
        {
            //
            //  For 2D datasets, compute the Delaunay triangulation.
            //
            case 2:
                triangle = new int[3 * 2 * n];
                int[] triangle_neighbor = new int[3 * 2 * n];

                Delauney.dtris2(n, 0, ref z, ref nt, ref triangle, ref triangle_neighbor);
                Console.WriteLine("");
                Console.WriteLine("  Triangulated data generates " + nt + " triangles.");
                break;
            default:
                nt = 0;
                break;
        }

        Console.WriteLine("");
        Console.WriteLine("  Measures of uniform point dispersion.");
        Console.WriteLine("");
        Console.WriteLine("  The pointset was read from \"" + input_filename + "\".");
        Console.WriteLine("  The sample routine will be SAMPLE_HYPERCUBE_UNIFORM.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =      " + dim_num + "");
        Console.WriteLine("  The number of points N =         " + n + "");
        Console.WriteLine("  The number of sample points NS = " + ns + "");
        Console.WriteLine("  The random number SEED_INIT =    " + seed_init + "");
        Console.WriteLine("");

        switch (dim_num)
        {
            case 2:
                triangle_order = 3;
                Console.WriteLine("  The minimum angle measure    Alpha = "
                                  + Quality.alpha_measure(n, z, triangle_order, nt, triangle) + "");
                break;
            default:
                Console.WriteLine("  The minimum angle measure    Alpha = (omitted)");
                break;
        }

        Console.WriteLine("  Relative spacing deviation    Beta = "
                          + Quality.beta_measure(dim_num, n, z) + "");
        Console.WriteLine("  The regularity measure         Chi = "
                          + Quality.chi_measure(dim_num, n, z, ns,
                              Burkardt.HyperGeometry.Hypercube.Sample.sample_hypercube_uniform,
                              seed_init) + "");
        Console.WriteLine("  2nd moment determinant measure   D = "
                          + Quality.d_measure(dim_num, n, z, ns,
                              Burkardt.HyperGeometry.Hypercube.Sample.sample_hypercube_uniform,
                              seed_init) + "");
        Console.WriteLine("  The Voronoi energy measure       E = "
                          + Quality.e_measure(dim_num, n, z, ns,
                              Burkardt.HyperGeometry.Hypercube.Sample.sample_hypercube_uniform,
                              seed_init) + "");
        Console.WriteLine("  The mesh ratio               Gamma = "
                          + Quality.gamma_measure(dim_num, n, z) + "");
        Console.WriteLine("  The point distribution norm      H = "
                          + Quality.h_measure(dim_num, n, z, ns,
                              Burkardt.HyperGeometry.Hypercube.Sample.sample_hypercube_uniform,
                              seed_init) + "");
        Console.WriteLine("  The COV measure             Lambda = "
                          + Quality.lambda_measure(dim_num, n, z) + "");
        Console.WriteLine("  The point distribution ratio    Mu = "
                          + Quality.mu_measure(dim_num, n, z, ns,
                              Burkardt.HyperGeometry.Hypercube.Sample.sample_hypercube_uniform,
                              seed_init) + "");
        Console.WriteLine("  The cell volume deviation       Nu = "
                          + Quality.nu_measure(dim_num, n, z, ns,
                              Burkardt.HyperGeometry.Hypercube.Sample.sample_hypercube_uniform,
                              seed_init) + "");
        switch (dim_num)
        {
            case 2:
                triangle_order = 2;
                Console.WriteLine("  The triangle uniformity measure  Q = "
                                  + Quality.q_measure(n, z, triangle_order, nt, triangle) + "");
                break;
            default:
                Console.WriteLine("  The triangle uniformity measure  Q = (omitted)");
                break;
        }

        Console.WriteLine("  The Riesz S = 0 energy measure  R0 = "
                          + Quality.r0_measure(dim_num, n, z) + "");
        if (typeMethods.r8mat_in_01(dim_num, n, z))
        {
            Console.WriteLine("  Nonintersecting sphere volume    S = "
                              + Quality.sphere_measure(dim_num, n, z) + "");
        }
        else
        {
            Console.WriteLine("  Nonintersecting sphere volume    S = (omitted)");
        }

        Console.WriteLine("  2nd moment trace measure       Tau = "
                          + Quality.tau_measure(dim_num, n, z, ns,
                              Burkardt.HyperGeometry.Hypercube.Sample.sample_hypercube_uniform,
                              seed_init) + "");

        double[] gamma = Spacing.pointset_spacing(dim_num, n, z);

        double gamma_ave = 0.0;
        for (i = 0; i < n; i++)
        {
            gamma_ave += gamma[i];
        }

        gamma_ave /= n;

        double gamma_min = typeMethods.r8vec_min(n, gamma);
        double gamma_max = typeMethods.r8vec_max(n, gamma);

        switch (n)
        {
            case > 1:
            {
                gamma_std = 0.0;
                for (i = 0; i < n; i++)
                {
                    gamma_std += Math.Pow(gamma[i] - gamma_ave, 2);
                }

                gamma_std = Math.Sqrt(gamma_std / (n - 1));
                break;
            }
            default:
                gamma_std = 0.0;
                break;
        }

        Console.WriteLine("");
        Console.WriteLine("  Minimum spacing          Gamma_min = " + gamma_min + "");
        Console.WriteLine("  Average spacing          Gamma_ave = " + gamma_ave + "");
        Console.WriteLine("  Maximum spacing          Gamma_max = " + gamma_max + "");
        Console.WriteLine("  Spacing standard deviate Gamma_std = " + gamma_std + "");

    }
}