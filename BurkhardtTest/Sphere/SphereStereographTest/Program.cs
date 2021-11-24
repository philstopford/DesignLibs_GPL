using System;
using Burkardt.PlaneNS;
using Burkardt.SphereNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace SphereStereographTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SPHERE_STEREOGRAPH_TEST.
        //
        //  Discussion:
        //
        //    SPHERE_STEREOGRAPH_TEST tests the SPHERE_STEREOGRAPH library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 November 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("SPHERE_STEREOGRAPH_TEST");
            
        Console.WriteLine("  Test the SPHERE_STEREOGRAPH library.");

        test01();
        test02();
        test03();

        Console.WriteLine("");
        Console.WriteLine("SPHERE_STEREOGRAPH_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 checks that the two functions are inverses.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 November 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        typeMethods.r8vecNormalData data = new();

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  SPHERE_STEREOGRAPH maps from sphere to plane.");
        Console.WriteLine("  SPHERE_STEREOGRAPH_INVERSE is the inverse map.");
        Console.WriteLine("  Check that these two functions are inverses.");

        int m = 3;
        int n = 100;
        //
        //  Check #1.
        //
        int seed = 123456789;
        double[] p1 = Burkardt.Uniform.Sphere.uniform_on_sphere01_map(m, n, ref data, ref seed);
        double[] q = Stereograph.sphere_stereograph(m, n, p1);
        double[] p2 = Stereograph.sphere_stereograph_inverse(m, n, q);
        double dif = 0.0;
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                dif += Math.Pow(p1[i + j * m] - p2[i + j * m], 2);
            }
        }

        dif = Math.Sqrt(dif) / n;
        Console.WriteLine("");
        Console.WriteLine("  Map points from sphere to plane to sphere.");
        Console.WriteLine("  RMS difference for " + n + " points was " + dif + "");

        //
        //  Check #2.
        //
        double[] q1 = UniformRNG.r8mat_uniform_01_new(m, n, ref seed);
        for (j = 0; j < n; j++)
        {
            q1[m - 1 + j * m] = 1.0;
        }

        double[] p = Stereograph.sphere_stereograph_inverse(m, n, q1);
        double[] q2 = Stereograph.sphere_stereograph(m, n, p);

        dif = 0.0;
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                dif += Math.Pow(q1[i + j * m] - q2[i + j * m], 2);
            }
        }

        dif = Math.Sqrt(dif) / n;
        Console.WriteLine("");
        Console.WriteLine("  Map points from plane to sphere to plane.");
        Console.WriteLine("  RMS difference for " + n + " points was " + dif + "");

    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 checks the generalized mapping.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 November 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        typeMethods.r8vecNormalData data = new();

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  SPHERE_STEREOGRAPH standard mapping from sphere to plane.");
        Console.WriteLine("  SPHERE_STEREOGRAPH2 generalized mapping:");
        Console.WriteLine("  (focus and center are arbitrary)");
        Console.WriteLine("  Check that these two functions can agree.");

        int m = 3;
        int n = 100;

        double[] focus = new double[m];
        for (i = 0; i < m - 1; i++)
        {
            focus[i] = 0.0;
        }

        focus[m - 1] = -1.0;

        double[] center = new double[m];
        for (i = 0; i < m; i++)
        {
            center[i] = 0.0;
        }

        //
        //  Check #1.
        //
        int seed = 123456789;
        double[] p1 = Burkardt.Uniform.Sphere.uniform_on_sphere01_map(m, n, ref data, ref seed);

        double[] q1 = Stereograph.sphere_stereograph(m, n, p1);

        double[] q2 = Stereograph.sphere_stereograph2(m, n, p1, focus, center);

        double dif = 0.0;
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                dif += Math.Pow(q1[i + j * m] - q2[i + j * m], 2);
            }
        }

        dif = Math.Sqrt(dif) / n;
        Console.WriteLine("");
        Console.WriteLine("  Map points from sphere to plane.");
        Console.WriteLine("  RMS difference for " + n + " points was " + dif + "");

        //
        //  Check #2.
        //
        q1 = UniformRNG.r8mat_uniform_01_new(m, n, ref seed);
        for (j = 0; j < n; j++)
        {
            q1[m - 1 + j * m] = 1.0;
        }

        p1 = Stereograph.sphere_stereograph_inverse(m, n, q1);

        double[] p2 = Stereograph.sphere_stereograph2_inverse(m, n, q1, focus, center);

        dif = 0.0;
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                dif += Math.Pow(p1[i + j * m] - p2[i + j * m], 2);
            }
        }

        dif = Math.Sqrt(dif) / n;
        Console.WriteLine("");
        Console.WriteLine("  Map points from plane to sphere.");
        Console.WriteLine("  RMS difference for " + n + " points was " + dif + "");

    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 checks that the two functions are inverses.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 November 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        typeMethods.r8vecNormalData data = new();

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  SPHERE_STEREOGRAPH2 maps from sphere to plane.");
        Console.WriteLine("  SPHERE_STEREOGRAPH2_INVERSE is the inverse map.");
        Console.WriteLine("  Check that these two functions are inverses.");

        const int m = 3;
        const int n = 100;
        int seed = 123456789;

        double[] focus = UniformRNG.r8vec_uniform_01_new(m, ref seed);
        double[] center = UniformRNG.r8vec_uniform_01_new(m, ref seed);
        double r = typeMethods.r8vec_norm_affine(m, focus, center);

        Console.WriteLine("");
        Console.WriteLine("  Using radius = " + r + "");
        typeMethods.r8vec_transpose_print(m, center, "  Center:");
        typeMethods.r8vec_transpose_print(m, focus, "  Focus:");
        //
        //  Check #1.
        //
        double[] p1 = Burkardt.Uniform.Sphere.uniform_on_sphere01_map(m, n, ref data, ref seed);

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                p1[i + j * m] = center[i] + r * p1[i + j * m];
            }
        }

        double[] q = Stereograph.sphere_stereograph2(m, n, p1, focus, center);

        double[] p2 = Stereograph.sphere_stereograph2_inverse(m, n, q, focus, center);

        double dif = 0.0;
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                dif += Math.Pow(p1[i + j * m] - p2[i + j * m], 2);
            }
        }

        dif = Math.Sqrt(dif) / n;
        Console.WriteLine("");
        Console.WriteLine("  Map points from sphere to plane to sphere.");
        Console.WriteLine("  RMS difference for " + n + " points was " + dif + "");

        //
        //  Check #2.
        //  We have to work hard to get random points on the plane, since
        //  all we know to begin with is the point of tangency and the normal.
        //
        double[] tang = new double[m];
        for (i = 0; i < m; i++)
        {
            tang[i] = 2.0 * center[i] - focus[i];
        }

        double[] normal = new double[m];
        for (i = 0; i < m; i++)
        {
            normal[i] = center[i] - focus[i];
        }

        double[] pr = new double[m];
        double[] pq = new double[m];
        Plane.plane_normal_basis_3d(tang, ref normal, ref pr, ref pq);

        double[] q1 = new double[m * n];

        double[] alpha = UniformRNG.r8vec_uniform_01_new(n, ref seed);
        double[] beta = UniformRNG.r8vec_uniform_01_new(n, ref seed);

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                q1[i + j * m] = tang[i] + pr[i] * alpha[j] + pq[i] * beta[j];
            }
        }

        double[] p = Stereograph.sphere_stereograph2_inverse(m, n, q1, focus, center);
        double[] q2 = Stereograph.sphere_stereograph2(m, n, p, focus, center);

        dif = 0.0;
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                dif += Math.Pow(q1[i + j * m] - q2[i + j * m], 2);
            }
        }

        dif = Math.Sqrt(dif) / n;
        Console.WriteLine("");
        Console.WriteLine("  Map points from plane to sphere to plane.");
        Console.WriteLine("  RMS difference for " + n + " points was " + dif + "");

    }
}