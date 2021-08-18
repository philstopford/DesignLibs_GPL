using System;
using Burkardt.ClenshawCurtisNS;
using Burkardt.Composition;
using Burkardt.IntegralNS;
using Burkardt.Quadrature;
using Burkardt.Sparse;
using Burkardt.Types;
using Burkardt.Uniform;

namespace SparseGridHeissWinschelTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for SPARSE_GRID_HW_TEST.
            //
            //  Discussion:
            //
            //    SPARSE_GRID_HW_TEST tests the SPARSE_GRID_HW library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 February 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("SPARSE_GRID_HW_TEST");
            Console.WriteLine("  Test the SPARSE_GRID_HW library.");

            ccl_test();
            ccl_sparse_test();

            ccs_test();
            ccs_sparse_test();

            cce_test();
            cce_sparse_test();

            get_seq_test();

            gqn_test();
            gqn_sparse_test();
            gqn2_sparse_test();

            gqu_test();
            gqu_sparse_test();

            kpn_test();
            kpn_sparse_test();

            kpu_test();
            kpu_sparse_test();

            nwspgr_size_test();
            nwspgr_time_test();
            nwspgr_test();

            order_report();

            symmetric_sparse_size_test();

            tensor_product_test();
            tensor_product_cell_test();

            Console.WriteLine("");
            Console.WriteLine("SPARSE_GRID_HW_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void ccl_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CCL_TEST uses CCL_ORDER + CC for 1D quadrature over [0,1].
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 February 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int d;
            double e;
            double exact;
            double[] fx;
            int l;
            int n;
            double q;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("CCL_TEST:");
            Console.WriteLine("  CCL_ORDER + CC");
            Console.WriteLine("  Clenshaw Curtis Linear (CCL) quadrature over [0,1]:");
            Console.WriteLine("");
            Console.WriteLine("   Level   Nodes    Estimate  Error");
            Console.WriteLine("");

            d = 1;
            exact = Integral.fu_integral(d);

            for (l = 1; l <= 5; l++)
            {
                n = ClenshawCurtis.ccl_order(l);

                x = new double[n];
                w = new double[n];

                ClenshawCurtis.cc(n, x, w);

                fx = Integral.fu_value(d, n, x);

                q = typeMethods.r8vec_dot_product(n, w, fx);

                e = Math.Abs(q - exact) / exact;

                Console.WriteLine("  " + l.ToString().PadLeft(2)
                                       + "  " + n.ToString().PadLeft(6)
                                       + "  " + q.ToString().PadLeft(14)
                                       + "  " + e.ToString().PadLeft(14) + "");

            }
        }

        static void ccl_sparse_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CCL_SPARSE_TEST uses CCL_ORDER + CC for a sparse grid.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 February 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Local parameters:
            //
            //    Local, int D, the spatial dimension.
            //
            //    Local, int MAXK, the maximum level to check.
            //
        {
            int d;
            double error_mc;
            double error_sg;
            double estimate;
            double[] fx;
            int k;
            int maxk;
            int n;
            int n2 = 0;
            int r;
            double[] s;
            int s_num;
            int seed;
            double trueval;
            double[] w;
            double[] x;

            d = 10;
            maxk = 7;

            trueval = Integral.fu_integral(d);

            Console.WriteLine("");
            Console.WriteLine("CCL_SPARSE_TEST:");
            Console.WriteLine("  CCL_ORDER + CC");
            Console.WriteLine("  Sparse Clenshaw Curtis Linear quadrature over [0,1].");
            Console.WriteLine("");
            Console.WriteLine("   D  Level   Nodes    SG error    MC error");
            Console.WriteLine("");

            for (k = 2; k <= maxk; k++)
            {
                //
                //  Compute sparse grid estimate.
                //
                n = Grid_NodesWeights.nwspgr_size(ClenshawCurtis.ccl_order, d, k);

                x = new double[d * n];
                w = new double[n];

                Grid_NodesWeights.nwspgr(ClenshawCurtis.cc, ClenshawCurtis.ccl_order, d, k, n, ref n2, ref x, ref w);

                fx = Integral.fu_value(d, n2, x);
                estimate = typeMethods.r8vec_dot_product(n2, w, fx);
                error_sg = Math.Abs((estimate - trueval) / trueval);
                //
                //  Compute 1000 Monte Carlo estimates with same number of points, and average.
                //
                s_num = 1000;
                s = new double[s_num];
                seed = 123456789;
                for (r = 0; r < 1000; r++)
                {
                    x = UniformRNG.r8mat_uniform_01_new(d, n2, ref seed);
                    fx = Integral.fu_value(d, n2, x);
                    s[r] = typeMethods.r8vec_sum(n2, fx) / (double)(n2);
                }

                error_mc = 0.0;
                for (r = 0; r < s_num; r++)
                {
                    error_mc = error_mc + Math.Pow(s[r] - trueval, 2);
                }

                error_mc = Math.Sqrt(error_mc / (double)(s_num)) / trueval;

                Console.WriteLine("  " + d.ToString().PadLeft(2)
                                       + "  " + k.ToString().PadLeft(5)
                                       + "  " + n2.ToString().PadLeft(6)
                                       + "  " + error_sg.ToString().PadLeft(10)
                                       + "  " + error_mc.ToString().PadLeft(10) + "");

            }
        }

        static void ccs_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CCS_TEST uses CCS_ORDER + CC for 1D quadrature over [0,1].
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 February 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int d;
            double e;
            double exact;
            double[] fx;
            int l;
            int n;
            double q;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("CCS_TEST:");
            Console.WriteLine("  CCS_ORDER + CC");
            Console.WriteLine("  Clenshaw Curtis Slow quadrature over [0,1]:");
            Console.WriteLine("");
            Console.WriteLine("   Level   Nodes    Estimate  Error");
            Console.WriteLine("");

            d = 1;
            exact = Integral.fu_integral(d);

            for (l = 1; l <= 5; l++)
            {
                n = ClenshawCurtis.ccs_order(l);

                x = new double[n];
                w = new double[n];

                ClenshawCurtis.cc(n, x, w);

                fx = Integral.fu_value(d, n, x);

                q = typeMethods.r8vec_dot_product(n, w, fx);

                e = Math.Abs(q - exact) / exact;

                Console.WriteLine("  " + l.ToString().PadLeft(2)
                                       + "  " + n.ToString().PadLeft(6)
                                       + "  " + q.ToString().PadLeft(14)
                                       + "  " + e.ToString().PadLeft(14) + "");

            }
        }

        static void ccs_sparse_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CCS_SPARSE_TEST uses CCS_ORDER + CC for a sparse grid.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 February 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Local parameters:
            //
            //    Local, int D, the spatial dimension.
            //
            //    Local, int MAXK, the maximum level to check.
            //
        {
            int d;
            double error_mc;
            double error_sg;
            double estimate;
            double[] fx;
            int k;
            int maxk;
            int n;
            int n2 = 0;
            int r;
            double[] s;
            int s_num;
            int seed;
            double trueval;
            double[] w;
            double[] x;

            d = 10;
            maxk = 7;

            trueval = Integral.fu_integral(d);

            Console.WriteLine("");
            Console.WriteLine("CCS_SPARSE_TEST:");
            Console.WriteLine("  CCS_ORDER + CC");
            Console.WriteLine("  Sparse Clenshaw Curtis Slow quadrature over [0,1].");
            Console.WriteLine("");
            Console.WriteLine("   D  Level   Nodes    SG error    MC error");
            Console.WriteLine("");

            for (k = 2; k <= maxk; k++)
            {
                //
                //  Compute sparse grid estimate.
                //
                n = Grid_NodesWeights.nwspgr_size(ClenshawCurtis.ccs_order, d, k);

                x = new double[d * n];
                w = new double[n];

                Grid_NodesWeights.nwspgr(ClenshawCurtis.cc, ClenshawCurtis.ccs_order, d, k, n, ref n2, ref x, ref w);

                fx = Integral.fu_value(d, n2, x);
                estimate = typeMethods.r8vec_dot_product(n2, w, fx);
                error_sg = Math.Abs((estimate - trueval) / trueval);
                //
                //  Compute 1000 Monte Carlo estimates with same number of points, and average.
                //
                s_num = 1000;
                s = new double[s_num];
                seed = 123456789;
                for (r = 0; r < 1000; r++)
                {
                    x = UniformRNG.r8mat_uniform_01_new(d, n2, ref seed);
                    fx = Integral.fu_value(d, n2, x);
                    s[r] = typeMethods.r8vec_sum(n2, fx) / (double)(n2);
                }

                error_mc = 0.0;
                for (r = 0; r < s_num; r++)
                {
                    error_mc = error_mc + Math.Pow(s[r] - trueval, 2);
                }

                error_mc = Math.Sqrt(error_mc / (double)(s_num)) / trueval;

                Console.WriteLine("  " + d.ToString().PadLeft(2)
                                       + "  " + k.ToString().PadLeft(5)
                                       + "  " + n2.ToString().PadLeft(6)
                                       + "  " + error_sg.ToString().PadLeft(10)
                                       + "  " + error_mc.ToString().PadLeft(10) + "");

            }

            return;
        }

        static void cce_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CCE_TEST uses CCE_ORDER + CC for 1D quadrature over [0,1].
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 February 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int d;
            double e;
            double exact;
            double[] fx;
            int l;
            int n;
            double q;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("CCE_TEST:");
            Console.WriteLine("  CCE_ORDER + CC");
            Console.WriteLine("  Clenshaw Curtis Exponential quadrature over [0,1]:");
            Console.WriteLine("");
            Console.WriteLine("   Level   Nodes    Estimate  Error");
            Console.WriteLine("");

            d = 1;
            exact = Integral.fu_integral(d);

            for (l = 1; l <= 5; l++)
            {
                n = ClenshawCurtis.cce_order(l);

                x = new double[n];
                w = new double[n];

                ClenshawCurtis.cc(n, x, w);

                fx = Integral.fu_value(d, n, x);

                q = typeMethods.r8vec_dot_product(n, w, fx);

                e = Math.Abs(q - exact) / exact;

                Console.WriteLine("  " + l.ToString().PadLeft(2)
                                       + "  " + n.ToString().PadLeft(6)
                                       + "  " + q.ToString().PadLeft(14)
                                       + "  " + e.ToString().PadLeft(14) + "");

            }
        }

        static void cce_sparse_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CCE_SPARSE_TEST uses CCE_ORDER + CC for a sparse grid.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 February 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Local parameters:
            //
            //    Local, int D, the spatial dimension.
            //
            //    Local, int MAXK, the maximum level to check.
            //
        {
            int d;
            double error_mc;
            double error_sg;
            double estimate;
            double[] fx;
            int k;
            int maxk;
            int n;
            int n2 = 0;
            int r;
            double[] s;
            int s_num;
            int seed;
            double trueval;
            double[] w;
            double[] x;

            d = 10;
            maxk = 7;

            trueval = Integral.fu_integral(d);

            Console.WriteLine("");
            Console.WriteLine("CCE_SPARSE_TEST:");
            Console.WriteLine("  CCE_ORDER + CC");
            Console.WriteLine("  Sparse Clenshaw Curtis Exponential quadrature over [0,1].");
            Console.WriteLine("");
            Console.WriteLine("   D  Level   Nodes    SG error    MC error");
            Console.WriteLine("");

            for (k = 2; k <= maxk; k++)
            {
                //
                //  Compute sparse grid estimate.
                //
                n = Grid_NodesWeights.nwspgr_size(ClenshawCurtis.cce_order, d, k);

                x = new double[d * n];
                w = new double[n];

                Grid_NodesWeights.nwspgr(ClenshawCurtis.cc, ClenshawCurtis.cce_order, d, k, n, ref n2, ref x, ref w);

                fx = Integral.fu_value(d, n2, x);
                estimate = typeMethods.r8vec_dot_product(n2, w, fx);
                error_sg = Math.Abs((estimate - trueval) / trueval);
                //
                //  Compute 1000 Monte Carlo estimates with same number of points, and average.
                //
                s_num = 1000;
                s = new double[s_num];
                seed = 123456789;
                for (r = 0; r < 1000; r++)
                {
                    x = UniformRNG.r8mat_uniform_01_new(d, n2, ref seed);
                    fx = Integral.fu_value(d, n2, x);
                    s[r] = typeMethods.r8vec_sum(n2, fx) / (double)(n2);
                }

                error_mc = 0.0;
                for (r = 0; r < s_num; r++)
                {
                    error_mc = error_mc + Math.Pow(s[r] - trueval, 2);
                }

                error_mc = Math.Sqrt(error_mc / (double)(s_num)) / trueval;

                Console.WriteLine("  " + d.ToString().PadLeft(2)
                                       + "  " + k.ToString().PadLeft(5)
                                       + "  " + n2.ToString().PadLeft(6)
                                       + "  " + error_sg.ToString().PadLeft(10)
                                       + "  " + error_mc.ToString().PadLeft(10) + "");

            }

            return;
        }

        static void get_seq_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GET_SEQ_TEST tests GET_SEQ.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 January 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int d;

            int[] fs;
            int norm;
            int seq_num;

            Console.WriteLine("");
            Console.WriteLine("GET_SEQ_TEST");
            Console.WriteLine("  GET_SEQ returns all D-dimensional vectors that sum to NORM.");

            d = 3;
            norm = 6;

            Console.WriteLine("");
            Console.WriteLine("  D = " + d + "");
            Console.WriteLine("  NORM = " + norm + "");

            seq_num = Comp.num_seq(norm - d, d);

            fs = Comp.get_seq(d, norm, seq_num);

            typeMethods.i4mat_print(seq_num, d, fs, "  The compositions");

        }

        static void gqn_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GQN_TEST uses the GQN function for 1D quadrature over (-oo,+oo).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 January 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int d;
            double e;
            double exact;
            double[] fx;
            int l;
            int n;
            double q;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("");
            Console.WriteLine("GQN_TEST:");
            Console.WriteLine("  Gauss-Hermite quadrature over (-oo,+oo):");
            Console.WriteLine("");
            Console.WriteLine("   Level   Nodes    Estimate  Error");
            Console.WriteLine("");

            d = 1;
            exact = Integral.fn_integral(d);

            for (l = 1; l <= 5; l++)
            {
                n = l;
                x = new double[n];
                w = new double[n];

                GaussQuadrature.gqn(n, x, w);

                fx = Integral.fn_value(d, n, x);

                q = typeMethods.r8vec_dot_product(n, w, fx);

                e = Math.Abs(q - exact) / exact;

                Console.WriteLine("  " + l.ToString().PadLeft(2)
                                       + "  " + n.ToString().PadLeft(6)
                                       + "  " + q.ToString().PadLeft(14)
                                       + "  " + e.ToString().PadLeft(14) + "");


            }

            return;
        }

        static void gqn_sparse_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GQN_SPARSE_TEST uses the GQN function to build a sparse grid.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 January 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Local parameters:
            //
            //    Local, int D, the spatial dimension.
            //
            //    Local, int MAXK, the maximum level to check.
            //
        {
            int d;
            double error_mc;
            double error_sg;
            double estimate;
            double[] fx;
            int k;
            int maxk;
            int n;
            int n2 = 0;
            int r;
            double[] s;
            int s_num;
            int seed;
            double trueval;
            double[] w;
            double[] x;
            typeMethods.r8vecNormalData data = new typeMethods.r8vecNormalData();

            d = 10;
            maxk = 7;

            trueval = Integral.fn_integral(d);

            Console.WriteLine("");
            Console.WriteLine("GQN_SPARSE_TEST:");
            Console.WriteLine("  GQN sparse grid:");
            Console.WriteLine("  Sparse Gaussian quadrature with Hermite weight over (-oo,+oo).");
            Console.WriteLine("");
            Console.WriteLine("   D  Level   Nodes    SG error    MC error");
            Console.WriteLine("");

            for (k = 2; k <= maxk; k++)
            {
                //
                //  Compute sparse grid estimate.
                //
                n = Grid_NodesWeights.nwspgr_size(GaussQuadrature.gqn_order, d, k);
                x = new double[d * n];
                w = new double[n];
                Grid_NodesWeights.nwspgr(GaussQuadrature.gqn, GaussQuadrature.gqn_order, d, k, n, ref n2, ref x, ref w);
                fx = Integral.fn_value(d, n2, x);
                estimate = typeMethods.r8vec_dot_product(n2, w, fx);

                error_sg = Math.Abs((estimate - trueval) / trueval);

                //
                //  Compute 1000 Monte Carlo estimates with same number of points, and average.
                //
                s_num = 1000;
                s = new double[s_num];
                seed = 123456789;
                for (r = 0; r < 1000; r++)
                {
                    x = typeMethods.r8mat_normal_01_new(d, n2, ref data, ref seed);
                    fx = Integral.fn_value(d, n2, x);
                    s[r] = typeMethods.r8vec_sum(n2, fx) / (double)(n2);
                }

                error_mc = 0.0;
                for (r = 0; r < s_num; r++)
                {
                    error_mc = error_mc + Math.Pow(s[r] - trueval, 2);
                }

                error_mc = Math.Sqrt(error_mc / (double)(s_num)) / trueval;

                Console.WriteLine("  " + d.ToString().PadLeft(2)
                                       + "  " + k.ToString().PadLeft(5)
                                       + "  " + n2.ToString().PadLeft(6)
                                       + "  " + error_sg.ToString().PadLeft(10)
                                       + "  " + error_mc.ToString().PadLeft(10) + "");

            }

            return;
        }

        static void gqn2_sparse_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GQN2_SPARSE_TEST uses the GQN and GQN2_ORDER functions.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 February 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Local parameters:
            //
            //    Local, int D, the spatial dimension.
            //
            //    Local, int MAXK, the maximum level to check.
            //
        {
            int d;
            int j;
            int k;
            int maxk;
            int n;
            int n2 = 0;
            double[] w;
            double[] x;

            d = 2;
            maxk = 4;

            Console.WriteLine("");
            Console.WriteLine("GQN2_SPARSE_TEST:");
            Console.WriteLine("  GQN sparse grid:");
            Console.WriteLine("  Gauss-Hermite sparse grids over (-oo,+oo).");
            Console.WriteLine("  Use GQN2_ORDER, the growth rule N = 2 * L - 1.");

            for (k = 2; k <= maxk; k++)
            {
                Console.WriteLine("");
                Console.WriteLine("     J      W                X               Y");
                Console.WriteLine("");

                n = Grid_NodesWeights.nwspgr_size(GaussQuadrature.gqn2_order, d, k);

                x = new double[d * n];
                w = new double[n];
                Grid_NodesWeights.nwspgr(GaussQuadrature.gqn, GaussQuadrature.gqn2_order, d, k, n, ref n2, ref x,
                    ref w);

                for (j = 0; j < n2; j++)
                {
                    Console.WriteLine("  " + j.ToString().PadLeft(4)
                                           + "  " + w[j].ToString().PadLeft(14)
                                           + "  " + x[0 + j * d].ToString().PadLeft(14)
                                           + "  " + x[1 + j * d].ToString().PadLeft(14) + "");
                }

            }

            return;
        }

        static void gqu_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GQU_TEST uses the GQU function for 1D quadrature over [0,1].
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 January 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int d;
            double e;
            double exact;
            double[] fx;
            int l;
            int n;
            double q;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("GQU_TEST:");
            Console.WriteLine("  Gauss-Legendre quadrature over [0,1]:");
            Console.WriteLine("");
            Console.WriteLine("   Level   Nodes    Estimate  Error");
            Console.WriteLine("");

            d = 1;
            exact = Integral.fu_integral(d);

            for (l = 1; l <= 5; l++)
            {
                n = l;
                x = new double[n];
                w = new double[n];

                GaussQuadrature.gqu(n, x, w);

                fx = Integral.fu_value(d, n, x);

                q = typeMethods.r8vec_dot_product(n, w, fx);

                e = Math.Abs(q - exact) / exact;

                Console.WriteLine("  " + l.ToString().PadLeft(2)
                                       + "  " + n.ToString().PadLeft(6)
                                       + "  " + q.ToString().PadLeft(14)
                                       + "  " + e.ToString().PadLeft(14) + "");

            }

            return;
        }

        static void gqu_sparse_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GQU_SPARSE_TEST uses the GQU function to build a sparse grid.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 January 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Local parameters:
            //
            //    Local, int D, the spatial dimension.
            //
            //    Local, int MAXK, the maximum level to check.
            //
        {
            int d;
            double error_mc;
            double error_sg;
            double estimate;
            double[] fx;
            int k;
            int maxk;
            int n;
            int n2 = 0;
            int r;
            double[] s;
            int s_num;
            int seed;
            double trueval;
            double[] w;
            double[] x;

            d = 10;
            maxk = 7;

            trueval = Integral.fu_integral(d);

            Console.WriteLine("");
            Console.WriteLine("GQU_SPARSE_TEST:");
            Console.WriteLine("  GQU sparse grid:");
            Console.WriteLine("  Sparse Gaussian unweighted quadrature over [0,1].");
            Console.WriteLine("");
            Console.WriteLine("   D  Level   Nodes    SG error    MC error");
            Console.WriteLine("");

            for (k = 2; k <= maxk; k++)
            {
                //
                //  Compute sparse grid estimate.
                //
                n = Grid_NodesWeights.nwspgr_size(GaussQuadrature.gqu_order, d, k);
                x = new double[d * n];
                w = new double[n];
                Grid_NodesWeights.nwspgr(GaussQuadrature.gqu, GaussQuadrature.gqu_order, d, k, n, ref n2, ref x, ref w);
                fx = Integral.fu_value(d, n2, x);
                estimate = typeMethods.r8vec_dot_product(n2, w, fx);
                error_sg = Math.Abs((estimate - trueval) / trueval);
                //
                //  Compute 1000 Monte Carlo estimates with same number of points, and average.
                //
                s_num = 1000;
                s = new double[s_num];
                seed = 123456789;
                for (r = 0; r < 1000; r++)
                {
                    x = UniformRNG.r8mat_uniform_01_new(d, n2, ref seed);
                    fx = Integral.fu_value(d, n2, x);
                    s[r] = typeMethods.r8vec_sum(n2, fx) / (double)(n2);
                }

                error_mc = 0.0;
                for (r = 0; r < s_num; r++)
                {
                    error_mc = error_mc + Math.Pow(s[r] - trueval, 2);
                }

                error_mc = Math.Sqrt(error_mc / (double)(s_num)) / trueval;

                Console.WriteLine("  " + d.ToString().PadLeft(2)
                                       + "  " + k.ToString().PadLeft(5)
                                       + "  " + n2.ToString().PadLeft(6)
                                       + "  " + error_sg.ToString().PadLeft(14)
                                       + "  " + error_mc.ToString().PadLeft(14) + "");

            }

            return;
        }

        static void kpn_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KPN_TEST uses the KPN function for 1D quadrature over (-oo,+oo).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 January 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int d;
            double e;
            double exact;
            double[] fx;
            int l;
            int n;
            double q;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("KPN_TEST:");
            Console.WriteLine("  Kronrod-Patterson-Hermite quadrature over (-oo,+oo):");
            Console.WriteLine("");
            Console.WriteLine("   Level   Nodes    Estimate  Error");
            Console.WriteLine("");

            d = 1;
            exact = Integral.fn_integral(d);

            for (l = 1; l <= 5; l++)
            {
                n = KronrodPatterson.kpn_order(l);
                x = new double[n];
                w = new double[n];
                KronrodPatterson.kpn(n, x, w);

                fx = Integral.fn_value(d, n, x);

                q = typeMethods.r8vec_dot_product(n, w, fx);

                e = Math.Abs(q - exact) / exact;

                Console.WriteLine("  " + l.ToString().PadLeft(2)
                                       + "  " + n.ToString().PadLeft(6)
                                       + "  " + q.ToString().PadLeft(14)
                                       + "  " + e.ToString().PadLeft(14) + "");

            }

            return;
        }

        static void kpn_sparse_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KPN_SPARSE_TEST uses the KPN function to build a sparse grid.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 January 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Local parameters:
            //
            //    Local, int D, the spatial dimension.
            //
            //    Local, int MAXK, the maximum level to check.
            //
        {
            int d;
            double error_mc;
            double error_sg;
            double estimate;
            double[] fx;
            int k;
            int maxk;
            int n;
            int n2 = 0;
            int r;
            double[] s;
            int s_num;
            int seed;
            double trueval;
            double[] w;
            double[] x;
            typeMethods.r8vecNormalData data = new typeMethods.r8vecNormalData();

            d = 10;
            maxk = 7;

            trueval = Integral.fn_integral(d);

            Console.WriteLine("");
            Console.WriteLine("KPN_SPARSE_TEST:");
            Console.WriteLine("  KPN sparse grid:");
            Console.WriteLine("  Sparse Kronrod-Patterson quadrature with Hermite weight over (-oo,+oo).");
            Console.WriteLine("");
            Console.WriteLine("   D  Level   Nodes    SG error    MC error");
            Console.WriteLine("");

            for (k = 2; k <= maxk; k++)
            {
                //
                //  Compute sparse grid estimate.
                //
                n = Grid_NodesWeights.nwspgr_size(KronrodPatterson.kpn_order, d, k);
                x = new double[d * n];
                w = new double[n];
                Grid_NodesWeights.nwspgr(KronrodPatterson.kpn, KronrodPatterson.kpn_order, d, k, n, ref n2, ref x,
                    ref w);
                fx = Integral.fn_value(d, n2, x);
                estimate = typeMethods.r8vec_dot_product(n2, w, fx);

                error_sg = Math.Abs((estimate - trueval) / trueval);

                //
                //  Compute 1000 Monte Carlo estimates with same number of points, and average.
                //
                s_num = 1000;
                s = new double[s_num];
                seed = 123456789;
                for (r = 0; r < 1000; r++)
                {
                    x = typeMethods.r8mat_normal_01_new(d, n2, ref data, ref seed);
                    fx = Integral.fn_value(d, n2, x);
                    s[r] = typeMethods.r8vec_sum(n2, fx) / (double)(n2);
                }

                error_mc = 0.0;
                for (r = 0; r < s_num; r++)
                {
                    error_mc = error_mc + Math.Pow(s[r] - trueval, 2);
                }

                error_mc = Math.Sqrt(error_mc / (double)(s_num)) / trueval;

                Console.WriteLine("  " + d.ToString().PadLeft(2)
                                       + "  " + k.ToString().PadLeft(5)
                                       + "  " + n2.ToString().PadLeft(6)
                                       + "  " + error_sg.ToString().PadLeft(10)
                                       + "  " + error_mc.ToString().PadLeft(10) + "");

            }

            return;
        }

        static void kpu_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KPU_TEST uses the KPU function for 1D quadrature over [0,1].
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 January 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int d;
            double e;
            double exact;
            double[] fx;
            int l;
            int n;
            double q;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("KPU_TEST:");
            Console.WriteLine("  Kronrod-Patterson quadrature over [0,1]:");
            Console.WriteLine("");
            Console.WriteLine("   Level   Nodes    Estimate  Error");
            Console.WriteLine("");

            d = 1;
            exact = Integral.fu_integral(d);

            for (l = 1; l <= 5; l++)
            {
                n = KronrodPatterson.kpu_order(l);
                x = new double[n];
                w = new double[n];
                KronrodPatterson.kpu(n, x, w);

                fx = Integral.fu_value(d, n, x);

                q = typeMethods.r8vec_dot_product(n, w, fx);

                e = Math.Abs(q - exact) / exact;

                Console.WriteLine("  " + l.ToString().PadLeft(2)
                                       + "  " + n.ToString().PadLeft(6)
                                       + "  " + q.ToString().PadLeft(14)
                                       + "  " + e.ToString().PadLeft(14) + "");

            }

            return;
        }

        static void kpu_sparse_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KPU_SPARSE_TEST uses the KPU function to build a sparse grid.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 January 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Local parameters:
            //
            //    Local, int D, the spatial dimension.
            //
            //    Local, int MAXK, the maximum level to check.
            //
        {
            int d;
            double error_mc;
            double error_sg;
            double estimate;
            double[] fx;
            int k;
            int maxk;
            int n;
            int n2 = 0;
            int r;
            double[] s;
            int s_num;
            int seed;
            double trueval;
            double[] w;
            double[] x;

            d = 10;
            maxk = 7;

            trueval = Integral.fu_integral(d);

            Console.WriteLine("");
            Console.WriteLine("KPU_SPARSE_TEST:");
            Console.WriteLine("  KPU sparse grid:");
            Console.WriteLine("  Sparse Kronrod-Patterson unweighted quadrature over [0,1].");
            Console.WriteLine("");
            Console.WriteLine("   D  Level   Nodes    SG error    MC error");
            Console.WriteLine("");

            for (k = 2; k <= maxk; k++)
            {
                //
                //  Compute sparse grid estimate.
                //
                n = Grid_NodesWeights.nwspgr_size(KronrodPatterson.kpu_order, d, k);
                x = new double[d * n];
                w = new double[n];
                Grid_NodesWeights.nwspgr(KronrodPatterson.kpu, KronrodPatterson.kpu_order, d, k, n, ref n2, ref x,
                    ref w);
                fx = Integral.fu_value(d, n2, x);
                estimate = typeMethods.r8vec_dot_product(n2, w, fx);
                error_sg = Math.Abs((estimate - trueval) / trueval);
                //
                //  Compute 1000 Monte Carlo estimates with same number of points, and average.
                //
                s_num = 1000;
                s = new double[s_num];
                seed = 123456789;
                for (r = 0; r < 1000; r++)
                {
                    x = UniformRNG.r8mat_uniform_01_new(d, n2, ref seed);
                    fx = Integral.fu_value(d, n2, x);
                    s[r] = typeMethods.r8vec_sum(n2, fx) / (double)(n2);
                }

                error_mc = 0.0;
                for (r = 0; r < s_num; r++)
                {
                    error_mc = error_mc + Math.Pow(s[r] - trueval, 2);
                }

                error_mc = Math.Sqrt(error_mc / (double)(s_num)) / trueval;

                Console.WriteLine("  " + d.ToString().PadLeft(2)
                                       + "  " + k.ToString().PadLeft(5)
                                       + "  " + n2.ToString().PadLeft(6)
                                       + "  " + error_sg.ToString().PadLeft(10)
                                       + "  " + error_mc.ToString().PadLeft(10) + "");

            }

            return;
        }

        static void nwspgr_size_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NWSPGR_SIZE_TEST tests NWSPGR_SIZE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    09 January 2013
            //
            //  Author:
            //
            //    John Burkardt.
            //
        {
            int dim;
            int k;
            int r_size;

            Console.WriteLine("");
            Console.WriteLine("NWSPGR_SIZE_TEST:");
            Console.WriteLine("  NWSPGR_SIZE returns the size of a sparse grid, based on either:");
            Console.WriteLine("  one of the built-in 1D rules, or a family of 1D rules");
            Console.WriteLine("  supplied by the user.");

            dim = 2;
            k = 3;
            Console.WriteLine("");
            Console.WriteLine("  Kronrod-Patterson, [0,1], Dim " + dim + ", Level " + k + ", Symmetric");
            Console.WriteLine("");
            r_size = Grid_NodesWeights.nwspgr_size(KronrodPatterson.kpu_order, dim, k);
            Console.WriteLine("  Full          " + r_size + "");

            dim = 2;
            k = 3;
            Console.WriteLine("");
            Console.WriteLine("  Kronrod-Patterson, (-oo,+oo), Dim " + dim + ", Level " + k + ", Symmetric");
            Console.WriteLine("");
            r_size = Grid_NodesWeights.nwspgr_size(KronrodPatterson.kpn_order, dim, k);
            Console.WriteLine("  Full          " + r_size + "");

            dim = 2;
            k = 3;
            Console.WriteLine("");
            Console.WriteLine("  Gauss-Legendre, [0,1], Dim " + dim + ", Level " + k + ", Symmetric");
            Console.WriteLine("");
            r_size = Grid_NodesWeights.nwspgr_size(GaussQuadrature.gqu_order, dim, k);
            Console.WriteLine("  Full          " + r_size + "");

            dim = 2;
            k = 3;
            Console.WriteLine("");
            Console.WriteLine("  Gauss Hermite, (-oo,+oo), [0,1], Dim " + dim + ", Level " + k + ", Symmetric");
            Console.WriteLine("");
            r_size = Grid_NodesWeights.nwspgr_size(GaussQuadrature.gqn_order, dim, k);
            Console.WriteLine("  Full          " + r_size + "");

            dim = 2;
            k = 3;
            Console.WriteLine("");
            Console.WriteLine("  Clenshaw Curtis Exponential, [-1,+1], [0,1], Dim " + dim + ", Level " + k +
                              ", Unsymmetric");
            Console.WriteLine("");
            r_size = Grid_NodesWeights.nwspgr_size(ClenshawCurtis.cce_order, dim, k);
            Console.WriteLine("  Full          " + r_size + "");
            //
            //  Do a table.
            //
            Console.WriteLine("");
            Console.WriteLine("  Dimension / Level table for Clenshaw Curtis Exponential");
            Console.WriteLine("");
            string cout = " Dim: ";
            for (dim = 1; dim <= 10; dim++)
            {
                cout += "  " + dim.ToString().PadLeft(6);
            }

            Console.WriteLine(cout);
            Console.WriteLine("Level");
            for (k = 1; k <= 5; k++)
            {
                cout = "  " + k.ToString().PadLeft(2) + "  ";
                for (dim = 1; dim <= 10; dim++)
                {
                    r_size = Grid_NodesWeights.nwspgr_size(ClenshawCurtis.cce_order, dim, k);
                    cout += "  " + r_size.ToString().PadLeft(6);
                }

                Console.WriteLine(cout);
            }

            return;
        }

        static void nwspgr_time_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NWSPGR_TIME_TEST times NWSPGR.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 January 2013
            //
            //  Author:
            //
            //    John Burkardt.
            //
        {
            int dim;
            int k;
            double[] nodes;
            int r_size;
            int s_size = 0;
            DateTime t1;
            double[] weights;

            Console.WriteLine("");
            Console.WriteLine("  This function measures the time in seconds required by NWSPGR");
            Console.WriteLine("  to compute a sparse grid, based on either:");
            Console.WriteLine("  one of the built-in 1D rules, or a family of 1D rules");
            Console.WriteLine("  supplied by the user.");

            dim = 20;
            k = 5;
            Console.WriteLine("");
            Console.WriteLine("  Kronrod-Patterson, [0,1], Dim " + dim + ", Level " + k + ", Symmetric");
            Console.WriteLine("");
            r_size = Grid_NodesWeights.nwspgr_size(KronrodPatterson.kpu_order, dim, k);
            nodes = new double[dim * r_size];
            weights = new double[r_size];
            t1 = DateTime.Now;
            Grid_NodesWeights.nwspgr(KronrodPatterson.kpu, KronrodPatterson.kpu_order, dim, k, r_size, ref s_size,
                ref nodes, ref weights);
            Console.WriteLine("  Full          " + (DateTime.Now - t1).TotalSeconds + "");

            dim = 20;
            k = 5;
            Console.WriteLine("");
            Console.WriteLine("  Kronrod-Patterson, (-oo,+oo), Dim " + dim + ", Level " + k + ", Symmetric");
            Console.WriteLine("");
            r_size = Grid_NodesWeights.nwspgr_size(KronrodPatterson.kpn_order, dim, k);
            nodes = new double[dim * r_size];
            weights = new double[r_size];
            t1 = DateTime.Now;
            Grid_NodesWeights.nwspgr(KronrodPatterson.kpn, KronrodPatterson.kpn_order, dim, k, r_size, ref s_size,
                ref nodes, ref weights);
            Console.WriteLine("  Full          " + (DateTime.Now - t1).TotalSeconds + "");

            dim = 20;
            k = 5;
            Console.WriteLine("");
            Console.WriteLine("  Gauss-Legendre, [0,1], Dim " + dim + ", Level " + k + ", Symmetric");
            Console.WriteLine("");
            r_size = Grid_NodesWeights.nwspgr_size(GaussQuadrature.gqu_order, dim, k);
            nodes = new double[dim * r_size];
            weights = new double[r_size];
            t1 = DateTime.Now;
            Grid_NodesWeights.nwspgr(GaussQuadrature.gqu, GaussQuadrature.gqu_order, dim, k, r_size, ref s_size,
                ref nodes, ref weights);
            Console.WriteLine("  Full          " + (DateTime.Now - t1).TotalSeconds + "");

            dim = 20;
            k = 5;
            Console.WriteLine("");
            Console.WriteLine("  Gauss Hermite, (-oo,+oo), [0,1], Dim " + dim + ", Level " + k + ", Symmetric");
            Console.WriteLine("");
            r_size = Grid_NodesWeights.nwspgr_size(GaussQuadrature.gqn_order, dim, k);
            nodes = new double[dim * r_size];
            weights = new double[r_size];
            t1 = DateTime.Now;
            Grid_NodesWeights.nwspgr(GaussQuadrature.gqn, GaussQuadrature.gqn_order, dim, k, r_size, ref s_size,
                ref nodes, ref weights);
            Console.WriteLine("  Full          " + (DateTime.Now - t1).TotalSeconds + "");

            dim = 20;
            k = 5;
            Console.WriteLine("");
            Console.WriteLine("  Clenshaw Curtis Exponential, [-1,+1], [0,1], Dim " + dim + ", Level " + k +
                              ", Unsymmetric");
            Console.WriteLine("");
            r_size = Grid_NodesWeights.nwspgr_size(ClenshawCurtis.cce_order, dim, k);
            nodes = new double[dim * r_size];
            weights = new double[r_size];
            t1 = DateTime.Now;
            Grid_NodesWeights.nwspgr(ClenshawCurtis.cc, ClenshawCurtis.cce_order, dim, k, r_size, ref s_size, ref nodes,
                ref weights);
            Console.WriteLine("  Full          " + (DateTime.Now - t1).TotalSeconds + "");
            /*
            Do a table.
            */
            Console.WriteLine("");
            Console.WriteLine("  Dimension / Level table for Clenshaw Curtis Exponential");
            Console.WriteLine("");
            string cout = " Dim: ";
            for (dim = 1; dim <= 10; dim++)
            {
                cout += "  " + dim.ToString().PadLeft(6);
            }

            Console.WriteLine(cout);
            Console.WriteLine("Level");
            for (k = 1; k <= 5; k++)
            {
                cout = "  " + k.ToString().PadLeft(2) + "  ";
                for (dim = 1; dim <= 10; dim++)
                {
                    r_size = Grid_NodesWeights.nwspgr_size(ClenshawCurtis.cce_order, dim, k);
                    nodes = new double[dim * r_size];
                    weights = new double[r_size];
                    t1 = DateTime.Now;
                    Grid_NodesWeights.nwspgr(ClenshawCurtis.cc, ClenshawCurtis.cce_order, dim, k, r_size, ref s_size,
                        ref nodes, ref weights);
                    cout += "  " + (DateTime.Now - t1).TotalSeconds.ToString().PadLeft(10);
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
            cout = " Dim: ";
            for (dim = 11; dim <= 20; dim++)
            {
                cout += "  " + dim.ToString().PadLeft(6);
            }

            Console.WriteLine(cout);
            Console.WriteLine("Level");
            for (k = 1; k <= 5; k++)
            {
                cout = "  " + k.ToString().PadLeft(2) + "  ";
                for (dim = 11; dim <= 20; dim++)
                {
                    r_size = Grid_NodesWeights.nwspgr_size(ClenshawCurtis.cce_order, dim, k);
                    nodes = new double[dim * r_size];
                    weights = new double[r_size];
                    t1 = DateTime.Now;
                    Grid_NodesWeights.nwspgr(ClenshawCurtis.cc, ClenshawCurtis.cce_order, dim, k, r_size, ref s_size,
                        ref nodes, ref weights);
                    cout += "  " + (DateTime.Now - t1).TotalSeconds.ToString().PadLeft(10);
                }

                Console.WriteLine(cout);
            }

            return;
        }

        static void nwspgr_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NWSPGR_TEST tests NWSPGR.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 January 2013
            //
            //  Author:
            //
            //    John Burkardt.
            //
        {
            int dim;
            int k;
            double[] nodes;
            int r_size;
            int s_size = 0;
            double[] weights;

            Console.WriteLine("");
            Console.WriteLine("NWSPGR_TEST:");
            Console.WriteLine("  NWSPGR generates a sparse grid, based on either:");
            Console.WriteLine("  one of the built-in 1D rules, or a family of 1D rules");
            Console.WriteLine("  supplied by the user.");

            dim = 2;
            k = 3;
            r_size = Grid_NodesWeights.nwspgr_size(KronrodPatterson.kpu_order, dim, k);
            nodes = new double[dim * r_size];
            weights = new double[r_size];
            Grid_NodesWeights.nwspgr(KronrodPatterson.kpu, KronrodPatterson.kpu_order, dim, k, r_size, ref s_size,
                ref nodes, ref weights);
            QuadratureRule.quad_rule_print(dim, s_size, nodes, weights, "  Kronrod-Patterson, [0,1], Dim 2, Level 3");

            dim = 2;
            k = 3;
            r_size = Grid_NodesWeights.nwspgr_size(KronrodPatterson.kpn_order, dim, k);
            nodes = new double[dim * r_size];
            weights = new double[r_size];
            Grid_NodesWeights.nwspgr(KronrodPatterson.kpn, KronrodPatterson.kpn_order, dim, k, r_size, ref s_size,
                ref nodes, ref weights);
            QuadratureRule.quad_rule_print(dim, s_size, nodes, weights,
                "  Kronrod-Patterson, (-oo,+oo), Dim 2, Level 3");

            dim = 2;
            k = 3;
            r_size = Grid_NodesWeights.nwspgr_size(GaussQuadrature.gqu_order, dim, k);
            nodes = new double[dim * r_size];
            weights = new double[r_size];
            Grid_NodesWeights.nwspgr(GaussQuadrature.gqu, GaussQuadrature.gqu_order, dim, k, r_size, ref s_size,
                ref nodes, ref weights);
            QuadratureRule.quad_rule_print(dim, s_size, nodes, weights, "  Gauss-Legendre, [0,1], Dim 2, Level 3");

            dim = 2;
            k = 3;
            r_size = Grid_NodesWeights.nwspgr_size(GaussQuadrature.gqn_order, dim, k);
            nodes = new double[dim * r_size];
            weights = new double[r_size];
            Grid_NodesWeights.nwspgr(GaussQuadrature.gqn, GaussQuadrature.gqn_order, dim, k, r_size, ref s_size,
                ref nodes, ref weights);
            QuadratureRule.quad_rule_print(dim, s_size, nodes, weights,
                "  Gauss Hermite, (-oo,+oo), Dim 2, Level 3");

            dim = 2;
            k = 3;
            r_size = Grid_NodesWeights.nwspgr_size(ClenshawCurtis.cce_order, dim, k);
            nodes = new double[dim * r_size];
            weights = new double[r_size];
            Grid_NodesWeights.nwspgr(ClenshawCurtis.cc, ClenshawCurtis.cce_order, dim, k, r_size, ref s_size, ref nodes,
                ref weights);
            QuadratureRule.quad_rule_print(dim, s_size, nodes, weights,
                "  Clenshaw Curtis Exponential, [-1,+1], Dim 2, Level 3");

            return;
        }

        static void order_report()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ORDER_REPORT reports on the order of each family of rules.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 January 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int ap;
            int k;
            int[] kpn_order =
            {
                1, 3, 9, 19, 35
            };
            int l;
            int o;
            int rp;

            Console.WriteLine("");
            Console.WriteLine("ORDER_REPORT");
            Console.WriteLine("  For each family of rules, report:");
            Console.WriteLine("  L,  the level index,");
            Console.WriteLine("  RP, the required polynomial precision,");
            Console.WriteLine("  AP, the actual polynomial precision,");
            Console.WriteLine("  O,  the rule order (number of points).");

            Console.WriteLine("");
            Console.WriteLine("  GQN family");
            Console.WriteLine("  Gauss quadrature, exponential weight, (-oo,+oo)");
            Console.WriteLine("");
            Console.WriteLine("   L  RP  AP   O");
            Console.WriteLine("");

            for (l = 1; l <= 25; l++)
            {
                rp = 2 * l - 1;
                o = l;
                ap = 2 * o - 1;
                Console.WriteLine("  " + l.ToString().PadLeft(2)
                                       + "  " + rp.ToString().PadLeft(2)
                                       + "  " + ap.ToString().PadLeft(2)
                                       + "  " + o.ToString().PadLeft(2) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("  GQU family");
            Console.WriteLine("  Gauss quadrature, unit weight, [0,1]");
            Console.WriteLine("");
            Console.WriteLine("   L  RP  AP   O");
            Console.WriteLine("");

            for (l = 1; l <= 25; l++)
            {
                rp = 2 * l - 1;
                o = l;
                ap = 2 * o - 1;
                Console.WriteLine("  " + l.ToString().PadLeft(2)
                                       + "  " + rp.ToString().PadLeft(2)
                                       + "  " + ap.ToString().PadLeft(2)
                                       + "  " + o.ToString().PadLeft(2) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("  KPN family");
            Console.WriteLine("  Gauss-Kronrod-Patterson quadrature, exponential weight, (-oo,+oo)");
            Console.WriteLine("");
            Console.WriteLine("   L  RP  AP   O");
            Console.WriteLine("");

            k = 1;
            o = 1;
            ap = 1;

            for (l = 1; l <= 25; l++)
            {
                rp = 2 * l - 1;

                while (ap < rp)
                {
                    if (k == 5)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("  No higher order rule is available!");
                        break;
                    }

                    //
                    //  Can we use a simple rule?
                    //
                    if (rp < kpn_order[k])
                    {
                        o = rp;
                        ap = rp;
                    }
                    //
                    //  Otherwise, move to next higher rule.
                    //
                    else
                    {
                        k = k + 1;
                        ap = 2 * kpn_order[k - 1] - kpn_order[k - 2];
                        o = kpn_order[k - 1];
                    }
                }

                Console.WriteLine("  " + l.ToString().PadLeft(2)
                                       + "  " + rp.ToString().PadLeft(2)
                                       + "  " + ap.ToString().PadLeft(2)
                                       + "  " + o.ToString().PadLeft(2) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("  KPU family");
            Console.WriteLine("  Gauss-Kronrod-Patterson quadrature, unit weight, [0,1]");
            Console.WriteLine("");
            Console.WriteLine("   L  RP  AP   O");
            Console.WriteLine("");

            for (l = 1; l <= 25; l++)
            {
                rp = 2 * l - 1;
                o = 1;
                ap = 1;
                while (ap < rp)
                {
                    o = 2 * (o + 1) - 1;
                    ap = (3 * o + 1) / 2;
                }

                Console.WriteLine("  " + l.ToString().PadLeft(2)
                                       + "  " + rp.ToString().PadLeft(2)
                                       + "  " + ap.ToString().PadLeft(2)
                                       + "  " + o.ToString().PadLeft(2) + "");
            }

            return;
        }

        static void symmetric_sparse_size_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SYMMETRIC_SPARSE_SIZE_TEST tests SYMMETRIC_SPARSE_SIZE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 January 2013
            //
            //  Author:
            //
            //    John Burkardt.
            //
            //  Local parameters:
            //
            //    Local, int D, the spatial dimension.
            //
            //    Local, int MAXK, the maximum level to check.
            //
        {
            int test_num = 3;

            int dim;
            int[] dim_test = { 5, 5, 3 };
            double[] nodes1 =
            {
                0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
                0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                0.0, 1.0, 0.0, 0.0, 0.0, 0.0
            };
            double[] nodes2 =
            {
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                1.73205,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.73205, 0.0, 0.0, 0.0, 0.0, 1.0,
                0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.73205, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                0.0,
                0.0, 0.0, 0.0, 1.0, 1.0, 1.73205, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                0.0,
                0.0, 1.0, 1.73205, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                0.0
            };
            double[] nodes3 =
            {
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.741964, 1.0,
                1.0, 1.0, 1.0, 1.0, 1.0, 1.73205, 1.73205, 1.73205, 2.33441,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.741964, 1.0, 1.0, 1.0, 1.73205, 1.73205, 2.33441,
                0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.73205, 0.0, 0.0, 1.0, 0.0,
                0.0, 0.741964, 1.0, 1.73205, 2.33441, 0.0, 0.0, 1.0, 1.73205, 0.0, 1.0, 0.0,
                0.0, 0.0, 1.0, 1.73205, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0
            };
            int r;
            int[] r_test = { 6, 21, 23 };
            int r2 = 0;
            int test;
            double x0;

            Console.WriteLine("");
            Console.WriteLine("SYMMETRIC_SPARSE_SIZE_TEST");
            Console.WriteLine("  Given a symmetric sparse grid rule represented only by");
            Console.WriteLine("  the points with positive values, determine the total number");
            Console.WriteLine("  of points in the grid.");
            Console.WriteLine("");
            Console.WriteLine("  For dimension DIM, we report");
            Console.WriteLine("  R, the number of points in the positive orthant, and");
            Console.WriteLine("  R2, the total number of points.");
            Console.WriteLine("");
            Console.WriteLine("       DIM         R        R2");
            Console.WriteLine("");

            x0 = 0.0;

            for (test = 0; test < test_num; test++)
            {
                r = r_test[test];
                dim = dim_test[test];
                if (test == 0)
                {
                    r2 = SparseRule.symmetric_sparse_size(r, dim, nodes1, x0);
                }
                else if (test == 1)
                {
                    r2 = SparseRule.symmetric_sparse_size(r, dim, nodes2, x0);
                }
                else if (test == 2)
                {
                    r2 = SparseRule.symmetric_sparse_size(r, dim, nodes3, x0);
                }

                Console.WriteLine("  " + dim.ToString().PadLeft(8)
                                       + "  " + r.ToString().PadLeft(8)
                                       + "  " + r2.ToString().PadLeft(8) + "");
            }

            return;
        }

        static void tensor_product_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TENSOR_PRODUCT_TEST tests TENSOR_PRODUCT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 January 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int d;
            int i;
            int j;
            int n;
            int n1d;
            int order1 = 2;
            int order2 = 3;
            int order3 = 2;
            int[] order1d;
            double[] w1d;
            double[] wnd;
            double[] w1_1d = { 1.0, 1.0 };
            double[] w2_1d = { 0.25, 0.50, 0.25 };
            double[] w3_1d = { 2.50, 2.50 };
            double[] x1_1d = { -1.0, +1.0 };
            double[] x2_1d = { 2.0, 2.5, 3.0 };
            double[] x3_1d = { 10.0, 15.0 };
            double[] x1d;
            double[] xnd;

            Console.WriteLine("");
            Console.WriteLine("TENSOR_PRODUCT_TEST:");
            Console.WriteLine("  Given a sequence of 1D quadrature rules, construct the");
            Console.WriteLine("  tensor product rule.");
            //
            //  1D rule.
            //
            d = 1;
            order1d = new int[d];

            order1d[0] = order1;

            n1d = typeMethods.i4vec_sum(d, order1d);
            x1d = new double[n1d];
            w1d = new double[n1d];

            n = typeMethods.i4vec_product(d, order1d);
            xnd = new double[d * n];
            wnd = new double[n];

            j = 0;
            for (i = 0; i < order1; i++)
            {
                x1d[j] = x1_1d[i];
                w1d[j] = w1_1d[i];
                j = j + 1;
            }

            typeMethods.tensor_product(d, order1d, n1d, x1d, w1d, n, ref xnd, ref wnd);

            QuadratureRule.quad_rule_print(d, n, xnd, wnd, "  A 1D rule over [-1,+1]:");

            //
            //  2D rule.
            //
            d = 2;
            order1d = new int[d];

            order1d[0] = order1;
            order1d[1] = order2;

            n1d = typeMethods.i4vec_sum(d, order1d);
            x1d = new double[n1d];
            w1d = new double[n1d];

            n = typeMethods.i4vec_product(d, order1d);
            xnd = new double[d * n];
            wnd = new double[n];

            j = 0;
            for (i = 0; i < order1; i++)
            {
                x1d[j] = x1_1d[i];
                w1d[j] = w1_1d[i];
                j = j + 1;
            }

            for (i = 0; i < order2; i++)
            {
                x1d[j] = x2_1d[i];
                w1d[j] = w2_1d[i];
                j = j + 1;
            }

            typeMethods.tensor_product(d, order1d, n1d, x1d, w1d, n, ref xnd, ref wnd);

            QuadratureRule.quad_rule_print(d, n, xnd, wnd, "  A 2D rule over [-1,+1] x [2.0,3.0]:");

            //
            //  3D rule.
            //
            d = 3;
            order1d = new int[d];

            order1d[0] = order1;
            order1d[1] = order2;
            order1d[2] = order3;

            n1d = typeMethods.i4vec_sum(d, order1d);
            x1d = new double[n1d];
            w1d = new double[n1d];

            n = typeMethods.i4vec_product(d, order1d);
            xnd = new double[d * n];
            wnd = new double[n];

            j = 0;
            for (i = 0; i < order1; i++)
            {
                x1d[j] = x1_1d[i];
                w1d[j] = w1_1d[i];
                j = j + 1;
            }

            for (i = 0; i < order2; i++)
            {
                x1d[j] = x2_1d[i];
                w1d[j] = w2_1d[i];
                j = j + 1;
            }

            for (i = 0; i < order3; i++)
            {
                x1d[j] = x3_1d[i];
                w1d[j] = w3_1d[i];
                j = j + 1;
            }

            typeMethods.tensor_product(d, order1d, n1d, x1d, w1d, n, ref xnd, ref wnd);

            QuadratureRule.quad_rule_print(d, n, xnd, wnd,
                "  A 3D rule over [-1,+1] x [2.0,3.0] x [10.0,15.0]:");

        }

        static void tensor_product_cell_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TENSOR_PRODUCT_CELL_TEST tests TENSOR_PRODUCT_CELL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 January 2013
            //
            //  Author:
            //
            //   John Burkardt
            //
        {
            int d;
            int nc;
            int np;
            int[] nr = { 2, 3, 2 };
            int[] roff;
            double[] wc;
            double[] wp;
            double[] w1_1d = { 1.0, 1.0 };
            double[] w2_1d = { 0.25, 0.50, 0.25 };
            double[] w3_1d = { 2.50, 2.50 };
            double[] x1_1d = { -1.0, +1.0 };
            double[] x2_1d = { 2.0, 2.5, 3.0 };
            double[] x3_1d = { 10.0, 15.0 };
            double[] xc;
            double[] xp;

            Console.WriteLine("");
            Console.WriteLine("TENSOR_PRODUCT_TEST_CELL:");
            Console.WriteLine("  Given a set of 1D quadrature rules stored in a cell array,");
            Console.WriteLine("  construct the tensor product rule.");
            //
            //  We can construct ROFF once and for all.
            //
            roff = typeMethods.r8cvv_offset(3, nr);
            //
            //  1D rule.
            //
            d = 1;
            nc = typeMethods.i4vec_sum(d, nr);
            xc = new double[nc];
            typeMethods.r8cvv_rset(nc, xc, d, roff, 0, x1_1d);
            wc = new double[nc];
            typeMethods.r8cvv_rset(nc, wc, d, roff, 0, w1_1d);
            np = typeMethods.i4vec_product(d, nr);
            xp = new double[d * np];
            wp = new double[np];
            typeMethods.tensor_product_cell(nc, xc, wc, d, nr, roff, np, ref xp, ref wp);
            QuadratureRule.quad_rule_print(d, np, xp, wp, "  A 1D rule over [-1,+1]:");

            //
            //  2D rule.
            //
            d = 2;
            nc = typeMethods.i4vec_sum(d, nr);
            xc = new double[nc];
            typeMethods.r8cvv_rset(nc, xc, d, roff, 0, x1_1d);
            typeMethods.r8cvv_rset(nc, xc, d, roff, 1, x2_1d);
            wc = new double[nc];
            typeMethods.r8cvv_rset(nc, wc, d, roff, 0, w1_1d);
            typeMethods.r8cvv_rset(nc, wc, d, roff, 1, w2_1d);
            np = typeMethods.i4vec_product(d, nr);
            xp = new double[d * np];
            wp = new double[np];

            typeMethods.tensor_product_cell(nc, xc, wc, d, nr, roff, np, ref xp, ref wp);

            QuadratureRule.quad_rule_print(d, np, xp, wp, "  A 1D rule over [-1,+1]:");

            //
            //  3D rule.
            //
            d = 3;
            nc = typeMethods.i4vec_sum(d, nr);
            xc = new double[nc];
            typeMethods.r8cvv_rset(nc, xc, d, roff, 0, x1_1d);
            typeMethods.r8cvv_rset(nc, xc, d, roff, 1, x2_1d);
            typeMethods.r8cvv_rset(nc, xc, d, roff, 2, x3_1d);
            wc = new double[nc];
            typeMethods.r8cvv_rset(nc, wc, d, roff, 0, w1_1d);
            typeMethods.r8cvv_rset(nc, wc, d, roff, 1, w2_1d);
            typeMethods.r8cvv_rset(nc, wc, d, roff, 2, w3_1d);
            np = typeMethods.i4vec_product(d, nr);
            xp = new double[d * np];
            wp = new double[np];

            typeMethods.tensor_product_cell(nc, xc, wc, d, nr, roff, np, ref xp, ref wp);

            QuadratureRule.quad_rule_print(d, np, xp, wp, "  A 1D rule over [-1,+1]:");

        }
    }
}