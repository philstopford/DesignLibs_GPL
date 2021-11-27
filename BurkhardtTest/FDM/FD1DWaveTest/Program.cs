﻿using System;
using Burkardt.FDM;
using Burkardt.Function;
using Burkardt.Types;

namespace FD1DWaveTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for FD1D_WAVE_TEST.
        //
        //  Discussion:
        //
        //    FD1D_WAVE_TEST tests the FD1D_WAVE library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 January 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("FD1D_WAVE_TEST");
            
        Console.WriteLine("  Test the FD1D_WAVE library.");

        fd1d_wave_test01();
        fd1d_wave_test02();

        Console.WriteLine("");
        Console.WriteLine("FD1D_WAVE_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void fd1d_wave_test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FD1D_WAVE_TEST_01 tests the FD1D finite difference wave computation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 January 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        int t_num = 41;
        int x_num = 16;

        Console.WriteLine("");
        Console.WriteLine("FD1D_WAVE_TEST01");
        Console.WriteLine("  Try the \"shark\" wave.");

        double x1 = 0.0;
        double x2 = 1.5;
        double[] x_vec = typeMethods.r8vec_linspace_new(x_num, x1, x2);

        double t1 = 0.0;
        double t2 = 4.0;
        double[] t_vec = typeMethods.r8vec_linspace_new(t_num, t1, t2);
        double t_delta = (t2 - t1) / (t_num - 1);

        double c = 1.0;
        double alpha = FD1D_Wave.fd1d_wave_alpha(x_num, x1, x2, t_num, t1, t2, c);
        //
        //  Load the initial condition.
        //
        double[] u = new double[t_num * x_num];

        double[] u1 = u_t1_01(x_num, x_vec);

        for (j = 0; j < x_num; j++)
        {
            u[0 + j * t_num] = u1[j];
        }

        //
        //  Take the first step.
        //
        double t = t_vec[1];
        double[] u2 = FD1D_Wave.fd1d_wave_start(x_num, x_vec, t, t_delta, alpha, u_x1_01, u_x2_01,
            ut_t1_01, u1);

        for (j = 0; j < x_num; j++)
        {
            u[1 + j * t_num] = u2[j];
        }

        //
        //  Take all the other steps.
        //
        for (i = 2; i < t_num; i++)
        {
            t = t_vec[i];
            double[] u3 = FD1D_Wave.fd1d_wave_step(x_num, t, alpha, u_x1_01, u_x2_01, u1, u2);
            for (j = 0; j < x_num; j++)
            {
                u[i + j * t_num] = u3[j];
                u1[j] = u2[j];
                u2[j] = u3[j];
            }

        }

        //
        //  Write the solution to a file.
        //
        typeMethods.r8mat_write("test01_data.txt", t_num, x_num, u);

        Console.WriteLine("");
        Console.WriteLine("  Plot data written to \"test01_data.txt\".");

    }

    private static double u_x1_01(double t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    U_X1_01 evaluates U at the boundary X1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 January 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T, the time.
        //
        //    Output, double U, the value of U(T,X1).
        //
    {
        int nd = 6;
        double[] td =  {
                0.0, 0.10, 0.20, 0.30, 0.40, 0.50
            }
            ;
        double[] tv = new double[1];
        double[] ud =  {
                0.0, 2.0, 10.0, 8.0, 5.0, 0.0
            }
            ;

        tv[0] = t;

        double[] uv = Piecewise.piecewise_linear(nd, td, ud, 1, tv);

        double u = uv[0];

        return u;
    }

    private static double u_x2_01(double t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    U_X2_01 evaluates U at the boundary X2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 January 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T, the time.
        //
        //    Output, double U, the value of U(T,X2).
        //
    {
        double u = 0.0;

        return u;
    }

    private static double[] u_t1_01(int x_num, double[] x_vec)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    U_T1_01 evaluates U at the initial time T1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 January 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int X_NUM, the number of nodes.
        //
        //    Input, double X_VEC[X_NUM], the coordinates of the nodes.
        //
        //    Output, double U_T1_01[X_NUM], the value of U at the initial time. 
        //
    {
        int j;

        double[] u = new double[x_num];

        for (j = 0; j < x_num; j++)
        {
            u[j] = 0.0;
        }

        return u;
    }

    private static double[] ut_t1_01(int x_num, double[] x_vec)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UT_T1_01 evaluates dUdT at the initial time T1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 January 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int X_NUM, the number of nodes.
        //
        //    Input, double X_VEC[X_NUM], the coordinates of the nodes.
        //
        //    Output, double UT_T1_01[X_NUM], the value of dUdT at the initial time. 
        //
    {
        int j;

        double[] ut = new double[x_num];

        for (j = 0; j < x_num; j++)
        {
            ut[j] = 0.0;
        }

        return ut;
    }

    private static void fd1d_wave_test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FD1D_WAVE_TEST_02 tests the FD1D finite difference wave computation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 January 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        int t_num = 41;
        int x_num = 16;

        Console.WriteLine("");
        Console.WriteLine("FD1D_WAVE_TEST02");
        Console.WriteLine("  Try a sine curve.");

        double x1 = 0.0;
        double x2 = 1.5;
        double[] x_vec = typeMethods.r8vec_linspace_new(x_num, x1, x2);

        double t1 = 0.0;
        double t2 = 4.0;
        double[] t_vec = typeMethods.r8vec_linspace_new(t_num, t1, t2);
        double t_delta = (t2 - t1) / (t_num - 1);
        //
        //  Changing T2 to 4.5 is enough to push the algorithm into instability.
        //
        //  t2 = 4.5;
        //
        double c = 1.0;
        double alpha = FD1D_Wave.fd1d_wave_alpha(x_num, x1, x2, t_num, t1, t2, c);
        //
        //  Load the initial condition.
        //
        double[] u = new double[t_num * x_num];

        double[] u1 = u_t1_02(x_num, x_vec);

        for (j = 0; j < x_num; j++)
        {
            u[0 + j * t_num] = u1[j];
        }

        //
        //  Take the first step.
        //
        double t = t_vec[1];
        double[] u2 = FD1D_Wave.fd1d_wave_start(x_num, x_vec, t, t_delta, alpha, u_x1_02, u_x2_02,
            ut_t1_02, u1);

        for (j = 0; j < x_num; j++)
        {
            u[1 + j * t_num] = u2[j];
        }

        //
        //  Take all the other steps.
        //
        for (i = 2; i < t_num; i++)
        {
            t = t_vec[i];
            double[] u3 = FD1D_Wave.fd1d_wave_step(x_num, t, alpha, u_x1_02, u_x2_02, u1, u2);
            for (j = 0; j < x_num; j++)
            {
                u[i + j * t_num] = u3[j];
                u1[j] = u2[j];
                u2[j] = u3[j];
            }
        }

        //
        //  Write the solution to a file.
        //
        typeMethods.r8mat_write("test02_data.txt", t_num, x_num, u);

        Console.WriteLine("");
        Console.WriteLine("  Plot data written to \"test02_data.txt\".");

    }

    private static double u_x1_02(double t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    U_X1_02 evaluates U at the boundary X1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 January 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T, the time.
        //
        //    Output, double U_X1_02, the value of U(T,X1).
        //
    {
        double u = 0.0;

        return u;
    }

    private static double u_x2_02(double t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    U_X2_02 evaluates U at the boundary X2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 January 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T, the time.
        //
        //    Output, double U, the value of U(T,X2).
        //
    {
        double u = 0.0;

        return u;
    }

    private static double[] u_t1_02(int x_num, double[] x_vec)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    U_T1_02 evaluates U at the initial time T1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 January 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int X_NUM, the number of nodes.
        //
        //    Input, double X_VEC[X_NUM], the spatial node coordinates.
        //
        //    Output, double U[X_NUM], the value of U at the initial time,
        //    and every node.
        //
    {
        int j;

        double[] u = new double[x_num];

        for (j = 0; j < x_num; j++)
        {
            u[j] = Math.Sin(2.0 * Math.PI * x_vec[j]);
        }

        return u;
    }

    private static double[] ut_t1_02(int x_num, double[] x_vec)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UT_T1_02 evaluates dUdT at the initial time T1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 January 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int X_NUM, the number of spatial intervals.
        //
        //    Input, double X_VEC[X_NUM], the spatial node coordinates.
        //
        //    Output, double UT_T1_02[X_NUM], the value of dUdT at the initial time,
        //    and every node.
        //
    {
        int j;

        double[] ut = new double[x_num];

        for (j = 0; j < x_num; j++)
        {
            ut[j] = 0.0;
        }

        return ut;
    }
}