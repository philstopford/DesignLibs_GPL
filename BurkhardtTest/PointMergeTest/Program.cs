using System;
using System.Globalization;
using Burkardt.Types;
using Burkardt.Uniform;

namespace PointMergeTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for POINT_MERGE_TEST.
        //
        //  Discussion:
        //
        //    POINT_MERGE_TEST tests the POINT_MERGE library.
        //
        //    Compare correctness of the codes.
        //
        //    Compare speed of the codes.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int m;
        int n;
        int n_unique;
        int seed;
        double tol;

        Console.WriteLine(" ");
        Console.WriteLine("POINT_MERGE_TEST");
        Console.WriteLine("  Test the POINT_MERGE library.");
        //
        //  TEST01 gives me some confidence that, at least for zero-tolerance,
        //  the radial approach is accurate, as compared to the "Unique count"
        //  (which cannot be extended to a tolerance version in multidimensions)
        //  and the "Tol Unique Count", which is an O(N^2) algorithm.
        //
        m = 3;
        n = 10;
        n_unique = 7;
        seed = 123456789;
        test01(m, n, n_unique, seed);

        m = 4;
        n = 20;
        n_unique = 11;
        seed = 987654321;
        test01(m, n, n_unique, seed);
        //
        //  In TEST02, I want to compute the same data, but with "blurred"
        //  duplicates, and a tolerance version of the radial approach,
        //  compared to "Tol Unique Count".
        //
        m = 3;
        n = 10;
        n_unique = 7;
        tol = 0.00001;
        seed = 123456789;
        test02(m, n, n_unique, tol, seed);

        m = 4;
        n = 20;
        n_unique = 11;
        tol = 0.00001;
        seed = 987654321;
        test02(m, n, n_unique, tol, seed);
        //
        //  In TEST03, I want to measure the time required for a sequence
        //  of increasingly hard problems.
        //
        m = 3;
        n = 100;
        n_unique = n / 2;
        tol = 0.00001;
        seed = 123456789;
        test03(m, n, n_unique, tol, seed);

        m = 3;
        n = 1000;
        n_unique = n / 2;
        tol = 0.00001;
        seed = 123456789;
        test03(m, n, n_unique, tol, seed);

        m = 3;
        n = 10000;
        n_unique = n / 2;
        tol = 0.00001;
        seed = 123456789;
        test03(m, n, n_unique, tol, seed);

        switch (false)
        {
            case true:
                m = 3;
                n = 100000;
                n_unique = n / 2;
                tol = 0.00001;
                seed = 123456789;
                test03(m, n, n_unique, tol, seed);
                break;
        }

        //
        //  In TEST04, repeat TEST02, but now compute the index vector.
        //
        m = 3;
        n = 10;
        n_unique = 7;
        tol = 0.00001;
        seed = 123456789;
        test04(m, n, n_unique, tol, seed);

        m = 4;
        n = 20;
        n_unique = 11;
        tol = 0.00001;
        seed = 987654321;
        test04(m, n, n_unique, tol, seed);
        //
        //  In TEST05, I want to measure the time required for a sequence
        //  of increasingly hard problems.
        //
        m = 3;
        n = 100;
        n_unique = n / 2;
        tol = 0.00001;
        seed = 123456789;
        test05(m, n, n_unique, tol, seed);

        m = 3;
        n = 1000;
        n_unique = n / 2;
        tol = 0.00001;
        seed = 123456789;
        test05(m, n, n_unique, tol, seed);

        m = 3;
        n = 10000;
        n_unique = n / 2;
        tol = 0.00001;
        seed = 123456789;
        test05(m, n, n_unique, tol, seed);

        switch (false)
        {
            case true:
                m = 3;
                n = 100000;
                n_unique = n / 2;
                tol = 0.00001;
                seed = 123456789;
                test05(m, n, n_unique, tol, seed);
                break;
        }

        Console.WriteLine(" ");
        Console.WriteLine("POINT_MERGE_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine(" ");
    }

    private static void test01(int m, int n, int n_unique, int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests uniqueness counting with no tolerance.
        //
        //  Discussion:
        //
        //    POINT_UNIQUE_COUNT uses an O(N) algorithm.
        //    POINT_RADIAL_UNIQUE_COUNT uses an algorithm that should be,
        //      in general, O(N);
        //    POINT_TOL_UNIQUE_COUNT uses an O(N^2) algorithm.
        //
        //    For this test, we just want to make sure the algorithms agree
        //    in the counting.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a;
        double tol;
        int unique_num;

        Console.WriteLine(" ");
        Console.WriteLine("TEST01");
        Console.WriteLine("  To count the unique columns in an typeMethods.r8COL, we call");
        Console.WriteLine("  POINT_UNIQUE_COUNT,");
        Console.WriteLine("  POINT_RADIAL_UNIQUE_COUNT, (with random center)");
        Console.WriteLine("  POINT_TOL_UNIQUE_COUNT, (with zero tolerance)");
        Console.WriteLine(" ");
        Console.WriteLine("  M =     " + m + "");
        Console.WriteLine("  N =     " + n + "");
        Console.WriteLine("  SEED =  " + seed + "");

        a = typeMethods.r8col_duplicates(m, n, n_unique, ref seed);

        typeMethods.r8mat_transpose_print(m, n, a, "  Matrix with N_UNIQUE unique columns:");

        Console.WriteLine(" ");
        Console.WriteLine("  N_UNIQUE =                  " + n_unique + "");

        unique_num = typeMethods.point_unique_count(m, n, a);
        Console.WriteLine("  POINT_UNIQUE_COUNT =        " + unique_num + "");

        unique_num = typeMethods.point_radial_unique_count(m, n, a, ref seed);
        Console.WriteLine("  POINT_RADIAL_UNIQUE_COUNT = " + unique_num + "");

        tol = 0.0;
        unique_num = typeMethods.point_tol_unique_count(m, n, a, tol);
        Console.WriteLine("  POINT_TOL_UNIQUE_COUNT =    " + unique_num + "");
    }

    private static void test02(int m, int n, int n_unique, double tol, int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests uniqueness counting with a tolerance.
        //
        //  Discussion:
        //
        //    POINT_RADIAL_TOL_UNIQUE_COUNT uses an algorithm that should be,
        //      in general, O(N);
        //    POINT_TOL_UNIQUE_COUNT uses an O(N^2) algorithm.
        //
        //    For this test, we just want to make sure the algorithms agree
        //    in the counting.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a;
        int i;
        int j;
        double[] r;
        double r_norm;
        int unique_num;

        Console.WriteLine(" ");
        Console.WriteLine("TEST02");
        Console.WriteLine("  To count the unique columns in an typeMethods.r8COL, we call");
        Console.WriteLine("  POINT_RADIAL_TOL_UNIQUE_COUNT, (with random center)");
        Console.WriteLine("  POINT_TOL_UNIQUE_COUNT, (with zero tolerance)");
        Console.WriteLine(" ");
        Console.WriteLine("  M =     " + m + "");
        Console.WriteLine("  N =     " + n + "");
        Console.WriteLine("  TOL =  " + tol + "");
        Console.WriteLine("  SEED =  " + seed + "");

        a = typeMethods.r8col_duplicates(m, n, n_unique, ref seed);

        typeMethods.r8mat_transpose_print(m, n, a, "  Matrix with N_UNIQUE unique columns:");
        //
        //  The form of the tolerance test means that if two vectors are initially
        //  equal, they remain "tolerably equal" after the addition of random
        //  perturbation vectors whose 2-norm is no greater than TOL/2.
        //
        r = new double[m];

        for (j = 0; j < n; j++)
        {
            UniformRNG.r8vec_uniform_01(m, ref seed, ref r);
            r_norm = typeMethods.r8vec_norm_l2(m, r);
            for (i = 0; i < m; i++)
            {
                a[i + j * m] += 0.5 * tol * r[i] / r_norm;
            }
        }

        typeMethods.r8mat_transpose_print(m, n, a, "  Blurred matrix:");

        Console.WriteLine(" ");
        Console.WriteLine("  N_UNIQUE =                      " + n_unique + "");

        unique_num = typeMethods.point_radial_tol_unique_count(m, n, a, tol, ref seed);
        Console.WriteLine("  POINT_RADIAL_TOL_UNIQUE_COUNT = " + unique_num + "");

        unique_num = typeMethods.point_tol_unique_count(m, n, a, tol);
        Console.WriteLine("  POINT_TOL_UNIQUE_COUNT =        " + unique_num + "");

    }

    private static void test03(int m, int n, int n_unique, double tol, int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 compares timings for two uniqueness counters.
        //
        //  Discussion:
        //
        //    POINT_RADIAL_TOL_UNIQUE_COUNT uses an algorithm that should be,
        //      in general, O(N);
        //    POINT_TOL_UNIQUE_COUNT uses an O(N^2) algorithm.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a;
        DateTime ctime;
        int i;
        int j;
        double[] r;
        double r_norm;
        int unique_num;

        Console.WriteLine(" ");
        Console.WriteLine("TEST03");
        Console.WriteLine("  To count the unique columns in an typeMethods.r8COL, we call");
        Console.WriteLine("  POINT_RADIAL_TOL_UNIQUE_COUNT, (with random center)");
        Console.WriteLine("  POINT_TOL_UNIQUE_COUNT, (with zero tolerance)");
        Console.WriteLine(" ");
        Console.WriteLine("  M =     " + m + "");
        Console.WriteLine("  N =     " + n + "");
        Console.WriteLine("  TOL =  " + tol + "");
        Console.WriteLine("  SEED =  " + seed + "");

        a = typeMethods.r8col_duplicates(m, n, n_unique, ref seed);
        //
        //  The form of the tolerance test means that if two vectors are initially
        //  equal, they remain "tolerably equal" after the addition of random
        //  perturbation vectors whose 2-norm is no greater than TOL/2.
        //
        r = new double[m];

        for (j = 0; j < n; j++)
        {
            UniformRNG.r8vec_uniform_01(m, ref seed, ref r);
            r_norm = typeMethods.r8vec_norm_l2(m, r);
            for (i = 0; i < m; i++)
            {
                a[i + j * m] += 0.5 * tol * r[i] / r_norm;
            }
        }

        Console.WriteLine(" ");
        Console.WriteLine("  N_UNIQUE =                      " + n_unique + "");

        ctime = DateTime.Now;
        unique_num = typeMethods.point_radial_tol_unique_count(m, n, a, tol, ref seed);
        Console.WriteLine(" ");
        Console.WriteLine("  POINT_RADIAL_TOL_UNIQUE_COUNT = " + unique_num + "");
        Console.WriteLine("  Time = " + (DateTime.Now - ctime).TotalSeconds + "");

        ctime = DateTime.Now;
        unique_num = typeMethods.point_tol_unique_count(m, n, a, tol);
        Console.WriteLine(" ");
        Console.WriteLine("  POINT_TOL_UNIQUE_COUNT =        " + unique_num + "");
        Console.WriteLine("  Time = " + (DateTime.Now - ctime).TotalSeconds + "");
    }

    private static void test04(int m, int n, int n_unique, double tol, int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests uniqueness indexing with a tolerance.
        //
        //  Discussion:
        //
        //    POINT_RADIAL_TOL_UNIQUE_COUNT uses an algorithm that should be,
        //      in general, O(N);
        //    POINT_TOL_UNIQUE_COUNT uses an O(N^2) algorithm.
        //
        //    For this test, we just want to make sure the algorithms agree
        //    in the counting.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a;
        double dist;
        int i;
        int j;
        int k;
        double[] r;
        double r_norm;
        int[] undx;
        int unique_num;
        int[] xdnu;
        string cout;

        Console.WriteLine(" ");
        Console.WriteLine("TEST04");
        Console.WriteLine("  To index the unique columns in an typeMethods.r8COL, we call");
        Console.WriteLine("  POINT_RADIAL_TOL_UNIQUE_COUNT, (with random center)");
        Console.WriteLine("  POINT_TOL_UNIQUE_COUNT, (with zero tolerance)");
        Console.WriteLine(" ");
        Console.WriteLine("  M =     " + m + "");
        Console.WriteLine("  N =     " + n + "");
        Console.WriteLine("  TOL =  " + tol + "");
        Console.WriteLine("  SEED =  " + seed + "");

        a = typeMethods.r8col_duplicates(m, n, n_unique, ref seed);

        typeMethods.r8mat_transpose_print(m, n, a, "  Matrix with N_UNIQUE unique columns:");
        //
        //  The form of the tolerance test means that if two vectors are initially
        //  equal, they remain "tolerably equal" after the addition of random
        //  perturbation vectors whose 2-norm is no greater than TOL/2.
        //
        r = new double[m];

        for (j = 0; j < n; j++)
        {
            UniformRNG.r8vec_uniform_01(m, ref seed, ref r);
            r_norm = typeMethods.r8vec_norm_l2(m, r);
            for (i = 0; i < m; i++)
            {
                a[i + j * m] += 0.5 * tol * r[i] / r_norm;
            }
        }

        typeMethods.r8mat_transpose_print(m, n, a, "  Blurred matrix:");

        Console.WriteLine(" ");
        Console.WriteLine("  N_UNIQUE =                      " + n_unique + "");

        undx = new int[n];
        xdnu = new int[n];

        unique_num = typeMethods.point_radial_tol_unique_index(m, n, a, tol, ref seed, ref undx,
            ref xdnu);

        Console.WriteLine(" ");
        Console.WriteLine("  POINT_RADIAL_TOL_UNIQUE_INDEX");
        Console.WriteLine("  Unique_num = " + unique_num + "");

        typeMethods.i4vec_print(unique_num, undx, "  UNDX:");

        typeMethods.i4vec_print(n, xdnu, "  XDNU:");

        Console.WriteLine(" ");
        Console.WriteLine("  List of nonunique points P(J), represented by");
        Console.WriteLine("  point with index I(J).");
        Console.WriteLine(" ");
        Console.WriteLine("  J, P(J)");
        Console.WriteLine("  I(J), P(I(J))");
        Console.WriteLine("  || P(J) - P(I(J)) || (should be <= TOL)");
        Console.WriteLine(" ");
        for (j = 0; j < n; j++)
        {
            k = undx[xdnu[j]];
            if (j != k)
            {
                Console.WriteLine(" ");
                cout = "  " + j.ToString().PadLeft(4);
                for (i = 0; i < m; i++)
                {
                    cout += "  " + a[i + j * m];
                }

                Console.WriteLine(cout);
                cout = "  " + k.ToString().PadLeft(4);
                for (i = 0; i < m; i++)
                {
                    cout += "  " + a[i + k * m];
                }

                Console.WriteLine(cout);
                dist = 0.0;
                for (i = 0; i < m; i++)
                {
                    dist += Math.Pow(a[i + j * m] - a[i + k * m], 2);
                }

                dist = Math.Sqrt(dist);
                Console.WriteLine("          " + dist.ToString().PadLeft(10) + "");
            }
        }

        //
        //  The interpretation of XDNU is simpler for POINT_TOL_UNIQUE_INDEX.
        //
        unique_num = typeMethods.point_tol_unique_index(m, n, a, tol, ref xdnu);

        Console.WriteLine(" ");
        Console.WriteLine("  POINT_TOL_UNIQUE_INDEX");
        Console.WriteLine("  Unique_num = " + unique_num + "");
        Console.WriteLine(" ");
        Console.WriteLine("  List of nonunique points P(J), represented by");
        Console.WriteLine("  point with index I(J).");
        Console.WriteLine(" ");
        Console.WriteLine("  J, P(J)");
        Console.WriteLine("  I(J), P(I(J))");
        Console.WriteLine("  || P(J) - P(I(J)) || (should be <= TOL)");
        Console.WriteLine(" ");
        for (j = 0; j < n; j++)
        {
            k = xdnu[j];
            if (j != k)
            {
                Console.WriteLine(" ");
                cout = "  " + j.ToString().PadLeft(4);
                for (i = 0; i < m; i++)
                {
                    cout += "  " + a[i + j * m];
                }

                Console.WriteLine(cout);
                cout = "  " + k.ToString().PadLeft(4);
                for (i = 0; i < m; i++)
                {
                    cout += "  " + a[i + k * m];
                }

                Console.WriteLine(cout);
                dist = 0.0;
                for (i = 0; i < m; i++)
                {
                    dist += Math.Pow(a[i + j * m] - a[i + k * m], 2);
                }

                dist = Math.Sqrt(dist);
                Console.WriteLine("          " + dist.ToString().PadLeft(10) + "");
            }
        }
    }

    private static void test05(int m, int n, int n_unique, double tol, int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 times uniqueness indexing with a tolerance.
        //
        //  Discussion:
        //
        //    POINT_RADIAL_TOL_UNIQUE_COUNT uses an algorithm that should be,
        //      in general, O(N);
        //    POINT_TOL_UNIQUE_COUNT uses an O(N^2) algorithm.
        //
        //    For this test, we just want to make sure the algorithms agree
        //    in the counting.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a;
        DateTime ctime;
        int i;
        int j;
        double[] r;
        double r_norm;
        int[] undx;
        int unique_num;
        int[] xdnu;

        Console.WriteLine(" ");
        Console.WriteLine("TEST05");
        Console.WriteLine("  We time the computations in TEST04, calling");
        Console.WriteLine("  POINT_RADIAL_TOL_UNIQUE_COUNT, (with random center)");
        Console.WriteLine("  POINT_TOL_UNIQUE_COUNT, (with zero tolerance)");
        Console.WriteLine(" ");
        Console.WriteLine("  M =     " + m + "");
        Console.WriteLine("  N =     " + n + "");
        Console.WriteLine("  TOL =  " + tol + "");
        Console.WriteLine("  SEED =  " + seed + "");

        a = typeMethods.r8col_duplicates(m, n, n_unique, ref seed);
        //
        //  The form of the tolerance test means that if two vectors are initially
        //  equal, they remain "tolerably equal" after the addition of random
        //  perturbation vectors whose 2-norm is no greater than TOL/2.
        //
        r = new double[m];

        for (j = 0; j < n; j++)
        {
            UniformRNG.r8vec_uniform_01(m, ref seed, ref r);
            r_norm = typeMethods.r8vec_norm_l2(m, r);
            for (i = 0; i < m; i++)
            {
                a[i + j * m] += 0.5 * tol * r[i] / r_norm;
            }
        }

        Console.WriteLine(" ");
        Console.WriteLine("  N_UNIQUE =                      " + n_unique + "");

        undx = new int[n];
        xdnu = new int[n];

        ctime = DateTime.Now;
        unique_num = typeMethods.point_radial_tol_unique_index(m, n, a, tol, ref seed, ref undx,
            ref xdnu);

        Console.WriteLine(" ");
        Console.WriteLine("  POINT_RADIAL_TOL_UNIQUE_INDEX");
        Console.WriteLine("  Unique_num = " + unique_num + "");
        Console.WriteLine("  Time = " + (DateTime.Now - ctime).TotalSeconds + "");

        ctime = DateTime.Now;
        unique_num = typeMethods.point_tol_unique_index(m, n, a, tol, ref xdnu);

        Console.WriteLine(" ");
        Console.WriteLine("  POINT_TOL_UNIQUE_INDEX");
        Console.WriteLine("  Unique_num = " + unique_num + "");
        Console.WriteLine("  Time = " + (DateTime.Now - ctime).TotalSeconds + "");

    }
}