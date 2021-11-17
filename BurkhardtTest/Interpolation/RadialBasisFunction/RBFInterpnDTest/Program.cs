using System;
using Burkardt.Interpolation;
using Burkardt.Types;
using Burkardt.Uniform;

namespace RBFInterpnDTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for RBF_INTERP_ND_TEST.
        //
        //  Discussion:
        //
        //    RBF_INTERP_ND_TEST tests the RBF_INTERP_ND library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("RBF_INTERP_ND_TEST:");
        Console.WriteLine("  Test the RBF_INTERP_ND library.");
        Console.WriteLine("  The R8LIB library is also needed.");

        rbf_interp_nd_test01();
        rbf_interp_nd_test02();
        rbf_interp_nd_test03();
        rbf_interp_nd_test04();

        Console.WriteLine("");
        Console.WriteLine("RBF_INTERP_ND_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");

    }

    private static void rbf_interp_nd_test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RBF_INTERP_ND_TEST01 tests RBF_WEIGHTS and RBF_INTERP_ND with PHI1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 July 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a;
        double app_error;
        double b;
        double[] fd;
        double[] fe;
        double[] fi;
        int i;
        double int_error;
        int j;
        int m = 2;
        int n1d = 5;
        int nd;
        int ni;
        double r0;
        int seed;
        double[] w;
        double[] x1d;
        double[] xd;
        double[] xi;
        typeMethods.r8vecDPData data = new();

        Console.WriteLine("");
        Console.WriteLine("RBF_INTERP_ND_TEST01:");
        Console.WriteLine("  RBF_WEIGHT computes weights for RBF interpolation.");
        Console.WriteLine("  RBF_INTERP_ND evaluates the RBF interpolant.");
        Console.WriteLine("  Use the multiquadratic basis function PHI1(R).");

        a = 0.0;
        b = 2.0;

        x1d = typeMethods.r8vec_linspace_new(n1d, a, b);
        nd = (int) Math.Pow(n1d, m);
        xd = new double[m * nd];

        for (i = 0; i < m; i++)
        {
            typeMethods.r8vec_direct_product(ref data, i, n1d, x1d, m, nd, ref xd);
        }

        typeMethods.r8mat_transpose_print(m, nd, xd, "  The product points:");

        r0 = (b - a) / n1d;

        Console.WriteLine("");
        Console.WriteLine("  Scale factor R0 = " + r0 + "");

        fd = new double[nd];
        for (j = 0; j < nd; j++)
        {
            fd[j] = xd[0 + j * m] * xd[1 + j * m] * Math.Exp(-xd[0 + j * m] * xd[1 + j * m]);
        }

        typeMethods.r8vec_print(nd, fd, "  Function data:");

        w = RadialBasisFunctions.rbf_weight(m, nd, xd, r0, RadialBasisFunctions.phi1, fd);

        typeMethods.r8vec_print(nd, w, "  Weight vector:");
        //
        //  #1: Interpolation test.  Does interpolant match function at interpolation points?
        //
        ni = nd;
        xi = typeMethods.r8mat_copy_new(m, ni, xd);

        fi = RadialBasisFunctions.rbf_interp_nd(m, nd, xd, r0, RadialBasisFunctions.phi1, w, ni, xi);

        int_error = typeMethods.r8vec_norm_affine(nd, fd, fi) / nd;

        Console.WriteLine("");
        Console.WriteLine("  L2 interpolation error averaged per interpolant node = " + int_error + "");

        //
        //  #2: Approximation test.  Estimate the integral (f-interp(f))^2.
        //
        ni = 1000;
        seed = 123456789;

        xi = UniformRNG.r8mat_uniform_ab_new(m, ni, a, b, ref seed);

        fi = RadialBasisFunctions.rbf_interp_nd(m, nd, xd, r0, RadialBasisFunctions.phi1, w, ni, xi);

        fe = new double[ni];
        for (j = 0; j < ni; j++)
        {
            fe[j] = xi[0 + j * m] * xi[1 + j * m] * Math.Exp(-xi[0 + j * m] * xi[1 + j * m]);
        }

        app_error = Math.Pow(b - a, m) * typeMethods.r8vec_norm_affine(ni, fi, fe) / ni;

        Console.WriteLine("");
        Console.WriteLine("  L2 approximation error averaged per 1000 samples = " + app_error + "");

    }

    private static void rbf_interp_nd_test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RBF_INTERP_ND_TEST02 tests RBF_WEIGHTS and RBF_INTERP_ND with PHI2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 July 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a;
        double app_error;
        double b;
        double[] fd;
        double[] fe;
        double[] fi;
        int i;
        double int_error;
        int j;
        int m = 2;
        int n1d = 5;
        int nd;
        int ni;
        double r0;
        int seed;
        double[] w;
        double[] x1d;
        double[] xd;
        double[] xi;
        typeMethods.r8vecDPData data = new();

        Console.WriteLine("");
        Console.WriteLine("RBF_INTERP_ND_TEST02:");
        Console.WriteLine("  RBF_WEIGHT computes weights for RBF interpolation.");
        Console.WriteLine("  RBF_INTERP_ND evaluates the RBF interpolant.");
        Console.WriteLine("  Use the inverse multiquadratic basis function PHI2(R).");

        a = 0.0;
        b = 2.0;

        x1d = typeMethods.r8vec_linspace_new(n1d, a, b);
        nd = (int) Math.Pow(n1d, m);
        xd = new double[m * nd];

        for (i = 0; i < m; i++)
        {
            typeMethods.r8vec_direct_product(ref data, i, n1d, x1d, m, nd, ref xd);
        }

        typeMethods.r8mat_transpose_print(m, nd, xd, "  The product points:");

        r0 = (b - a) / n1d;

        Console.WriteLine("");
        Console.WriteLine("  Scale factor R0 = " + r0 + "");

        fd = new double[nd];
        for (j = 0; j < nd; j++)
        {
            fd[j] = xd[0 + j * m] * xd[1 + j * m] * Math.Exp(-xd[0 + j * m] * xd[1 + j * m]);
        }

        typeMethods.r8vec_print(nd, fd, "  Function data:");

        w = RadialBasisFunctions.rbf_weight(m, nd, xd, r0, RadialBasisFunctions.phi2, fd);

        typeMethods.r8vec_print(nd, w, "  Weight vector:");
        //
        //  #1: Interpolation test.  Does interpolant match function at interpolation points?
        //
        ni = nd;
        xi = typeMethods.r8mat_copy_new(m, ni, xd);

        fi = RadialBasisFunctions.rbf_interp_nd(m, nd, xd, r0, RadialBasisFunctions.phi2, w, ni, xi);

        int_error = typeMethods.r8vec_norm_affine(nd, fd, fi) / nd;

        Console.WriteLine("");
        Console.WriteLine("  L2 interpolation error averaged per interpolant node = " + int_error + "");

        //
        //  #2: Approximation test.  Estimate the integral (f-interp(f))^2.
        //
        ni = 1000;
        seed = 123456789;

        xi = UniformRNG.r8mat_uniform_ab_new(m, ni, a, b, ref seed);

        fi = RadialBasisFunctions.rbf_interp_nd(m, nd, xd, r0, RadialBasisFunctions.phi2, w, ni, xi);

        fe = new double[ni];
        for (j = 0; j < ni; j++)
        {
            fe[j] = xi[0 + j * m] * xi[1 + j * m] * Math.Exp(-xi[0 + j * m] * xi[1 + j * m]);
        }

        app_error = Math.Pow(b - a, m) * typeMethods.r8vec_norm_affine(ni, fi, fe) / ni;

        Console.WriteLine("");
        Console.WriteLine("  L2 approximation error averaged per 1000 samples = " + app_error + "");

    }

    private static void rbf_interp_nd_test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RBF_INTERP_ND_TEST03 tests RBF_WEIGHTS and RBF_INTERP_ND with PHI3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 July 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a;
        double app_error;
        double b;
        double[] fd;
        double[] fe;
        double[] fi;
        int i;
        double int_error;
        int j;
        int m = 2;
        int n1d = 5;
        int nd;
        int ni;
        double r0;
        int seed;
        double[] w;
        double[] x1d;
        double[] xd;
        double[] xi;
        typeMethods.r8vecDPData data = new();

        Console.WriteLine("");
        Console.WriteLine("RBF_INTERP_ND_TEST03:");
        Console.WriteLine("  RBF_WEIGHT computes weights for RBF interpolation.");
        Console.WriteLine("  RBF_INTERP_ND evaluates the RBF interpolant.");
        Console.WriteLine("  Use the thin-plate spline basis function PHI3(R).");

        a = 0.0;
        b = 2.0;

        x1d = typeMethods.r8vec_linspace_new(n1d, a, b);
        nd = (int) Math.Pow(n1d, m);
        xd = new double[m * nd];

        for (i = 0; i < m; i++)
        {
            typeMethods.r8vec_direct_product(ref data, i, n1d, x1d, m, nd, ref xd);
        }

        typeMethods.r8mat_transpose_print(m, nd, xd, "  The product points:");

        r0 = (b - a) / n1d;

        Console.WriteLine("");
        Console.WriteLine("  Scale factor R0 = " + r0 + "");

        fd = new double[nd];
        for (j = 0; j < nd; j++)
        {
            fd[j] = xd[0 + j * m] * xd[1 + j * m] * Math.Exp(-xd[0 + j * m] * xd[1 + j * m]);
        }

        typeMethods.r8vec_print(nd, fd, "  Function data:");

        w = RadialBasisFunctions.rbf_weight(m, nd, xd, r0, RadialBasisFunctions.phi3, fd);

        typeMethods.r8vec_print(nd, w, "  Weight vector:");
        //
        //  #1: Interpolation test.  Does interpolant match function at interpolation points?
        //
        ni = nd;
        xi = typeMethods.r8mat_copy_new(m, ni, xd);

        fi = RadialBasisFunctions.rbf_interp_nd(m, nd, xd, r0, RadialBasisFunctions.phi3, w, ni, xi);

        int_error = typeMethods.r8vec_norm_affine(nd, fd, fi) / nd;

        Console.WriteLine("");
        Console.WriteLine("  L2 interpolation error averaged per interpolant node = " + int_error + "");

        //
        //  #2: Approximation test.  Estimate the integral (f-interp(f))^2.
        //
        ni = 1000;
        seed = 123456789;

        xi = UniformRNG.r8mat_uniform_ab_new(m, ni, a, b, ref seed);

        fi = RadialBasisFunctions.rbf_interp_nd(m, nd, xd, r0, RadialBasisFunctions.phi3, w, ni, xi);

        fe = new double[ni];
        for (j = 0; j < ni; j++)
        {
            fe[j] = xi[0 + j * m] * xi[1 + j * m] * Math.Exp(-xi[0 + j * m] * xi[1 + j * m]);
        }

        app_error = Math.Pow(b - a, m) * typeMethods.r8vec_norm_affine(ni, fi, fe) / ni;

        Console.WriteLine("");
        Console.WriteLine("  L2 approximation error averaged per 1000 samples = " + app_error + "");

    }

    private static void rbf_interp_nd_test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RBF_INTERP_ND_TEST04 tests RBF_WEIGHTS and RBF_INTERP_ND with PHI4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 July 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a;
        double app_error;
        double b;
        double[] fd;
        double[] fe;
        double[] fi;
        int i;
        double int_error;
        int j;
        int m = 2;
        int n1d = 5;
        int nd;
        int ni;
        double r0;
        int seed;
        double[] w;
        double[] x1d;
        double[] xd;
        double[] xi;
        typeMethods.r8vecDPData data = new();

        Console.WriteLine("");
        Console.WriteLine("RBF_INTERP_ND_TEST04:");
        Console.WriteLine("  RBF_WEIGHT computes weights for RBF interpolation.");
        Console.WriteLine("  RBF_INTERP_ND evaluates the RBF interpolant.");
        Console.WriteLine("  Use the gaussian basis function PHI4(R).");

        a = 0.0;
        b = 2.0;

        x1d = typeMethods.r8vec_linspace_new(n1d, a, b);
        nd = (int) Math.Pow(n1d, m);
        xd = new double[m * nd];

        for (i = 0; i < m; i++)
        {
            typeMethods.r8vec_direct_product(ref data, i, n1d, x1d, m, nd, ref xd);
        }

        typeMethods.r8mat_transpose_print(m, nd, xd, "  The product points:");

        r0 = (b - a) / n1d;

        Console.WriteLine("");
        Console.WriteLine("  Scale factor R0 = " + r0 + "");

        fd = new double[nd];
        for (j = 0; j < nd; j++)
        {
            fd[j] = xd[0 + j * m] * xd[1 + j * m] * Math.Exp(-xd[0 + j * m] * xd[1 + j * m]);
        }

        typeMethods.r8vec_print(nd, fd, "  Function data:");

        w = RadialBasisFunctions.rbf_weight(m, nd, xd, r0, RadialBasisFunctions.phi4, fd);

        typeMethods.r8vec_print(nd, w, "  Weight vector:");
        //
        //  #1: Interpolation test.  Does interpolant match function at interpolation points?
        //
        ni = nd;
        xi = typeMethods.r8mat_copy_new(m, ni, xd);

        fi = RadialBasisFunctions.rbf_interp_nd(m, nd, xd, r0, RadialBasisFunctions.phi4, w, ni, xi);

        int_error = typeMethods.r8vec_norm_affine(nd, fd, fi) / nd;

        Console.WriteLine("");
        Console.WriteLine("  L2 interpolation error averaged per interpolant node = " + int_error + "");

        //
        //  #2: Approximation test.  Estimate the integral (f-interp(f))^2.
        //
        ni = 1000;
        seed = 123456789;

        xi = UniformRNG.r8mat_uniform_ab_new(m, ni, a, b, ref seed);

        fi = RadialBasisFunctions.rbf_interp_nd(m, nd, xd, r0, RadialBasisFunctions.phi4, w, ni, xi);

        fe = new double[ni];
        for (j = 0; j < ni; j++)
        {
            fe[j] = xi[0 + j * m] * xi[1 + j * m] * Math.Exp(-xi[0 + j * m] * xi[1 + j * m]);
        }

        app_error = Math.Pow(b - a, m) * typeMethods.r8vec_norm_affine(ni, fi, fe) / ni;

        Console.WriteLine("");
        Console.WriteLine("  L2 approximation error averaged per 1000 samples = " + app_error + "");

    }
}