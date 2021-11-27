using System;
using System.Globalization;
using Burkardt.Composition;
using Burkardt.Interpolation;
using Burkardt.Lagrange;
using Burkardt.Types;
using Burkardt.Uniform;

namespace SparseInterpnDTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SPARSE_INTERP_ND_TEST.
        //
        //  Discussion:
        //
        //    SPARSE_INTERP_ND_TEST tests the SPARSE_INTERP_ND library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine(" ");
        Console.WriteLine("SPARSE_INTERP_ND_TEST");
        Console.WriteLine("  Test the SPARSE_INTERP_ND library.");
        Console.WriteLine("  The R8LIB library is also required.");

        int m = 1;
        int sparse_max = 9;
        test01(m, sparse_max);

        m = 2;
        sparse_max = 9;
        test01(m, sparse_max);

        m = 3;
        sparse_max = 9;
        test01(m, sparse_max);

        m = 4;
        sparse_max = 7;
        test01(m, sparse_max);
        //
        //  Terminate.
        //
        Console.WriteLine(" ");
        Console.WriteLine("SPARSE_INTERP_ND_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine(" ");
    }

    private static void test01(int m, int sparse_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01: sequence of sparse interpolants to an M-dimensional function.
        //
        //  Discussion:
        //
        //    We have functions that can generate a Lagrange interpolant to data
        //    in M dimensions, with specified order or level in each dimension.
        //
        //    We use the Lagrange function as the inner evaluator for a sparse
        //    grid procedure. 
        //
        //    The procedure computes sparse interpolants of levels 0 to SPARSE_MAX
        //    to a given function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Local, int M, the spatial dimension.
        //
        //    Input, int SPARSE_MAX, the maximum sparse grid level to try.
        //
        //  Local Parameters:
        //
        //    Local, double A[M], B[M], the upper and lower variable limits 
        //    in each dimension.
        //
        //    Local, double APP_ERROR, the averaged Euclidean norm of the 
        //    difference between the sparse interpolant and the exact function at 
        //    the interpolation points.
        //
        //    Local, int C[L_MAX+1], the sparse grid coefficient vector.
        //    Results at level L are weighted by C(L).
        //
        //    Local, int IND[M], the 1D indices defining a Lagrange grid.
        //    Each entry is a 1d "level" that specifies the order of a 
        //    Clenshaw Curtis 1D grid.
        //
        //    Local, int L, the current Lagrange grid level.
        //
        //    Local, int L_MAX, the current sparse grid level.
        //
        //    Local, int MORE, used to control the enumeration of all the
        //    Lagrange grids at a current grid level.
        //
        //    Local, int ND, the number of points used in a Lagrange grid.
        //
        //    Local, int ND_TOTAL, the total number of points used in all the
        //    Lagrange interpolants for a given sparse interpolant points that occur
        //    in more than one Lagrange grid are counted multiple times.
        //
        //    Local, int NI, the number of interpolant evaluation points.
        //
        //    Local, int SPARSE_MIN, the minimum sparse grid level to try.
        //
        //    Local, double XD[M*ND], the data points for a Lagrange grid.
        //
        //    Local, double XI[M*NI], the interpolant evaluation points.
        //
        //    Local, double ZD[ND], the data values for a Lagrange grid.
        //
        //    Local, double ZE[NI], the exact function values at XI.
        //
        //    Local, double ZI[NI], the sparse interpolant values at XI.
        //
        //    Local, double ZPI[NI], one set of Lagrange interpolant values at XI.
        //
    {
        int h = 0;
        int i;
        int l_max;
        int t = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  Sparse interpolation for a function f(x) of M-dimensional argument.");
        Console.WriteLine("  Use a sequence of sparse grids of levels 0 through SPARSE_MAX.");
        Console.WriteLine("  Invoke a general Lagrange interpolant function to do this.");
        Console.WriteLine("");
        Console.WriteLine("  Compare the exact function and the interpolants at a grid of points.");
        Console.WriteLine("");
        Console.WriteLine("  The \"order\" is the sum of the orders of all the product grids");
        Console.WriteLine("  used to make a particular sparse grid.");
        //
        //  User input.
        //
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension M = " + m + "");
        Console.WriteLine("  Maximum sparse grid level = " + sparse_max + "");
        //
        //  Define the region.
        //
        double[] a = new double[m];
        double[] b = new double[m];

        for (i = 0; i < m; i++)
        {
            a[i] = 0.0;
            b[i] = 1.0;
        }

        //
        //  Define the interpolation evaluation information.
        //
        const int ni = 100;
        int seed = 123456789;
        double[] xi = UniformRNG.r8mat_uniform_abvec_new(m, ni, a, b, ref seed);

        Console.WriteLine("  Number of interpolation points is NI = " + ni + "");

        double[] ze = f_sinr(m, ni, xi);
        //
        //  Compute a sequence of sparse grid interpolants of increasing level.
        //
        Console.WriteLine("");
        Console.WriteLine("   L     Order    ApproxError");
        Console.WriteLine("");

        int[] ind = new int[m];
        double[] zi = new double[ni];

        const int sparse_min = 0;

        for (l_max = sparse_min; l_max <= sparse_max; l_max++)
        {
            int[] c = new int[l_max + 1];
            int[] w = new int[l_max + 1];
            Smolyak.smolyak_coefficients(l_max, m, ref c, ref w);

            for (i = 0; i < ni; i++)
            {
                zi[i] = 0.0;
            }

            int nd_total = 0;

            int l_min = Math.Max(l_max + 1 - m, 0);

            int l;
            for (l = l_min; l <= l_max; l++)
            {
                bool more = false;
                while (true)
                {
                    //
                    //  Define the next product grid.
                    //
                    Comp.comp_next(l, m, ref ind, ref more, ref h, ref t);
                    //
                    //  Count the grid, find its points, evaluate the data there.
                    //
                    int nd = LagrangenD.lagrange_interp_nd_size2(m, ind);
                    double[] xd = LagrangenD.lagrange_interp_nd_grid2(m, ind, a, b, nd);
                    double[] zd = f_sinr(m, nd, xd);
                    //
                    //  Use the grid to evaluate the interpolant.
                    //
                    double[] zpi = LagrangenD.lagrange_interp_nd_value2(m, ind, a, b, nd, zd, ni, xi);
                    //
                    //  Weighted the interpolant values and add to the sparse grid interpolant.
                    //
                    nd_total += nd;
                    for (i = 0; i < ni; i++)
                    {
                        zi[i] += c[l] * zpi[i];
                    }

                    if (!more)
                    {
                        break;
                    }
                }
            }

            //
            //  Compare sparse interpolant and exact function at interpolation points.
            //
            double app_error = typeMethods.r8vec_norm_affine(ni, zi, ze) / ni;

            Console.WriteLine("  " + l.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + nd_total.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + app_error.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");

        }
    }

    private static double[] f_sinr(int m, int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F_SINR is a scalar function of an M-dimensional argument, to be interpolated.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, int N, the number of points.
        //
        //    Input, double X(M,N), the points.
        //
        //    Output, double F_SINR[N], the value of the function at each point.
        //
    {
        int j;

        double[] z = new double[n];

        for (j = 0; j < n; j++)
        {
            double r = 0.0;
            int i;
            for (i = 0; i < m; i++)
            {
                r += x[i + j * m] * x[i + j * m];
            }

            r = Math.Sqrt(r);
            z[j] = Math.Sin(r);
        }

        return z;
    }
}