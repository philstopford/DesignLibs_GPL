using System;
using Burkardt.ODENS.RungeKuttaFehlberg;
using Burkardt.Types;

namespace RungeKuttaFehlbergTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for RKF45_TEST.
            //
            //  Discussion:
            //
            //    RKF45_TEST tests the RKF45 library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 June 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("RKF45_TEST");
            Console.WriteLine("  C++ version");
            Console.WriteLine("  Test the RKF45 library.");

            test01();
            test02();
            test03();
            test04();
            test05();
            test06();
            //
            //  Terminate.
            //
            Console.WriteLine("");
            Console.WriteLine("RKF45_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");

        }

        static void test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 solves a scalar ODE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 June 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            float abserr;
            int flag;
            int i_step;
            int n_step;
            int neqn;
            float relerr;
            float t;
            float t_out;
            float t_start;
            float t_stop;
            float[] y;
            float[] yp;
            RungeKuttaFehlberg.r4RKFData data = new RungeKuttaFehlberg.r4RKFData();

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  Solve a scalar equation using R4_RKF:");
            Console.WriteLine("");
            Console.WriteLine("  Y' = 0.25 * Y * ( 1 - Y / 20 )");
            Console.WriteLine("");

            neqn = 1;

            y = new float[neqn];
            yp = new float[neqn];

            abserr = (float) Math.Sqrt(typeMethods.r4_epsilon());
            relerr = (float) Math.Sqrt(typeMethods.r4_epsilon());

            flag = 1;

            t_start = 0.0f;
            t_stop = 20.0f;

            n_step = 5;

            t = 0.0f;
            t_out = 0.0f;
            y[0] = 1.0f;
            yp = r4_f1(t, y, yp);

            Console.WriteLine("");
            Console.WriteLine("FLAG             T          Y         Y'          Y_Exact         Error");
            Console.WriteLine("");

            Console.WriteLine(flag.ToString().PadLeft(4) + "  "
                                                         + t.ToString().PadLeft(12) + "  "
                                                         + y[0].ToString().PadLeft(12) + "  "
                                                         + yp[0].ToString().PadLeft(12) + "  "
                                                         + r4_y1x(t).ToString().PadLeft(12) + "  "
                                                         + (y[0] - r4_y1x(t)).ToString().PadLeft(12) + "");
            ;

            for (i_step = 1; i_step <= n_step; i_step++)
            {
                t = ((float) (n_step - i_step + 1) * t_start
                     + (float) (i_step - 1) * t_stop)
                    / (float) (n_step);

                t_out = ((float) (n_step - i_step) * t_start
                         + (float) (i_step) * t_stop)
                        / (float) (n_step);

                flag = RungeKuttaFehlberg.r4_rkf45(ref data, r4_f1, neqn, ref y, ref yp, ref t, t_out, ref relerr, abserr,
                    flag);

                Console.WriteLine(flag.ToString().PadLeft(4) + "  "
                                                             + t.ToString().PadLeft(12) + "  "
                                                             + y[0].ToString().PadLeft(12) + "  "
                                                             + yp[0].ToString().PadLeft(12) + "  "
                                                             + r4_y1x(t).ToString().PadLeft(12) + "  "
                                                             + (y[0] - r4_y1x(t)).ToString().PadLeft(12) + "");
            }
        }

        static void test02()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 solves a vector ODE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 June 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            float abserr;
            int flag;
            int i_step;
            int n_step;
            int neqn;
            float relerr;
            float t;
            float t_out;
            float t_start;
            float t_stop;
            float[] y;
            float[] yp;
            RungeKuttaFehlberg.r4RKFData data = new RungeKuttaFehlberg.r4RKFData();

            Console.WriteLine("");
            Console.WriteLine("TEST02");
            Console.WriteLine("  Solve a vector equation using R4_RKF:");
            Console.WriteLine("");
            Console.WriteLine("  Y'(1) =  Y(2)");
            Console.WriteLine("  Y'(2) = -Y(1)");
            Console.WriteLine("");
            Console.WriteLine("");
            Console.WriteLine("  This system is equivalent to the following");
            Console.WriteLine("  second order system:");
            Console.WriteLine("");
            Console.WriteLine("  Z\" = - Z.");

            neqn = 2;

            y = new float[neqn];
            yp = new float[neqn];

            abserr = (float) Math.Sqrt(typeMethods.r4_epsilon());
            relerr = (float) Math.Sqrt(typeMethods.r4_epsilon());

            flag = 1;

            t_start = 0.0f;
            t_stop = 2.0f * 3.14159265f;

            n_step = 12;

            t = 0.0f;
            t_out = 0.0f;

            y[0] = 1.0f;
            y[1] = 0.0f;

            Console.WriteLine("");
            Console.WriteLine("FLAG             T          Y(1)       Y(2)");
            Console.WriteLine("");

            Console.WriteLine(flag.ToString().PadLeft(4) + "  "
                                                         + t.ToString().PadLeft(12) + "  "
                                                         + y[0].ToString().PadLeft(12) + "  "
                                                         + y[1].ToString().PadLeft(12) + "");

            for (i_step = 1; i_step <= n_step; i_step++)
            {
                t = ((float) (n_step - i_step + 1) * t_start
                     + (float) (i_step - 1) * t_stop)
                    / (float) (n_step);

                t_out = ((float) (n_step - i_step) * t_start
                         + (float) (i_step) * t_stop)
                        / (float) (n_step);

                flag = RungeKuttaFehlberg.r4_rkf45(ref data, r4_f2, neqn, ref y, ref yp, ref t, t_out, ref relerr, abserr,
                    flag);

                Console.WriteLine(flag.ToString().PadLeft(4) + "  "
                                                             + t.ToString().PadLeft(12) + "  "
                                                             + y[0].ToString().PadLeft(12) + "  "
                                                             + y[1].ToString().PadLeft(12) + "");
            }
        }

        static void test03()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST03 solves a scalar ODE using single step mode.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 June 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int NEQN = 1;

            float abserr;
            int flag;
            int i_step;
            int n_step;
            float relerr;
            float t;
            float t_out;
            float t_start;
            float t_stop;
            float[] y = new float[NEQN];
            float[] yp = new float[NEQN];
            RungeKuttaFehlberg.r4RKFData data = new RungeKuttaFehlberg.r4RKFData();

            Console.WriteLine("");
            Console.WriteLine("TEST03");
            Console.WriteLine("  Solve a scalar equation using R4_RKF:");
            Console.WriteLine("");
            Console.WriteLine("  Y' = 0.25 * Y * ( 1 - Y / 20 )");
            Console.WriteLine("");
            Console.WriteLine("  This routine uses the SINGLE STEP mode.");
            Console.WriteLine("");

            abserr = (float) Math.Sqrt(typeMethods.r4_epsilon());
            relerr = (float) Math.Sqrt(typeMethods.r4_epsilon());

            flag = -1;

            t_start = 0.0f;
            t_stop = 20.0f;

            n_step = 5;

            t = 0.0f;
            t_out = 0.0f;
            y[0] = 1.0f;
            yp = r4_f1(t, y, yp);

            Console.WriteLine("");
            Console.WriteLine("FLAG             T          Y         Y'        Y_Exact         Error");
            Console.WriteLine("");

            Console.WriteLine(flag.ToString().PadLeft(4) + "  "
                                                         + t.ToString().PadLeft(12) + "  "
                                                         + y[0].ToString().PadLeft(12) + "  "
                                                         + yp[0].ToString().PadLeft(12) + "  "
                                                         + r4_y1x(t).ToString().PadLeft(12) + "  "
                                                         + (y[0] - r4_y1x(t)).ToString().PadLeft(12) + "");
            ;

            for (i_step = 1; i_step <= n_step; i_step++)
            {
                t = ((float) (n_step - i_step + 1) * t_start
                     + (float) (i_step - 1) * t_stop)
                    / (float) (n_step);

                t_out = ((float) (n_step - i_step) * t_start
                         + (float) (i_step) * t_stop)
                        / (float) (n_step);

                while (flag < 0)
                {
                    flag = RungeKuttaFehlberg.r4_rkf45(ref data, r4_f1, NEQN, ref y, ref yp, ref t, t_out, ref relerr, abserr,
                        flag);

                    Console.WriteLine(flag.ToString().PadLeft(4) + "  "
                                                                 + t.ToString().PadLeft(12) + "  "
                                                                 + y[0].ToString().PadLeft(12) + "  "
                                                                 + yp[0].ToString().PadLeft(12) + "  "
                                                                 + r4_y1x(t).ToString().PadLeft(12) + "  "
                                                                 + (y[0] - r4_y1x(t)).ToString().PadLeft(12) + "");
                }

                flag = -2;
            }

            return;
# undef NEQN
        }

        static void test04()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST04 solves a scalar ODE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 June 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int NEQN = 1;

            double abserr;
            int flag;
            int i_step;
            int n_step;
            double relerr;
            double t;
            double t_out;
            double t_start;
            double t_stop;
            double[] y = new double[NEQN];
            double[] yp = new double[NEQN];
            RungeKuttaFehlberg.r8RKFData data = new RungeKuttaFehlberg.r8RKFData();

            Console.WriteLine("");
            Console.WriteLine("TEST04");
            Console.WriteLine("  Solve a scalar equation using R8_RKF:");
            Console.WriteLine("");
            Console.WriteLine("  Y' = 0.25 * Y * ( 1 - Y / 20 )");
            Console.WriteLine("");

            abserr = Math.Sqrt(typeMethods.r8_epsilon());
            relerr = Math.Sqrt(typeMethods.r8_epsilon());

            flag = 1;

            t_start = 0.0;
            t_stop = 20.0;

            n_step = 5;

            t = 0.0;
            t_out = 0.0;
            y[0] = 1.0;
            yp = r8_f1(t, y, yp);

            Console.WriteLine("");
            Console.WriteLine("FLAG             T          Y         Y'          Y_Exact         Error");
            Console.WriteLine("");

            Console.WriteLine(flag.ToString().PadLeft(4) + "  "
                                                         + t.ToString().PadLeft(12) + "  "
                                                         + y[0].ToString().PadLeft(12) + "  "
                                                         + yp[0].ToString().PadLeft(12) + "  "
                                                         + r8_y1x(t).ToString().PadLeft(12) + "  "
                                                         + (y[0] - r8_y1x(t)).ToString().PadLeft(12) + "");
            ;

            for (i_step = 1; i_step <= n_step; i_step++)
            {
                t = ((double) (n_step - i_step + 1) * t_start
                     + (double) (i_step - 1) * t_stop)
                    / (double) (n_step);

                t_out = ((double) (n_step - i_step) * t_start
                         + (double) (i_step) * t_stop)
                        / (double) (n_step);

                flag = RungeKuttaFehlberg.r8_rkf45(ref data, r8_f1, NEQN, ref y, ref yp, ref t, t_out, ref relerr, abserr,
                    flag);

                Console.WriteLine(flag.ToString().PadLeft(4) + "  "
                                                             + t.ToString().PadLeft(12) + "  "
                                                             + y[0].ToString().PadLeft(12) + "  "
                                                             + yp[0].ToString().PadLeft(12) + "  "
                                                             + r8_y1x(t).ToString().PadLeft(12) + "  "
                                                             + (y[0] - r8_y1x(t)).ToString().PadLeft(12) + "");
            }

        }

        static void test05()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST05 solves a vector ODE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 June 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int NEQN = 2;

            double abserr;
            int flag;
            int i_step;
            int n_step;
            double relerr;
            double t;
            double t_out;
            double t_start;
            double t_stop;
            double[] y = new double[NEQN];
            double[] yp = new double[NEQN];
            RungeKuttaFehlberg.r8RKFData data = new RungeKuttaFehlberg.r8RKFData();

            Console.WriteLine("");
            Console.WriteLine("TEST05");
            Console.WriteLine("  Solve a vector equation using R8_RKF:");
            Console.WriteLine("");
            Console.WriteLine("  Y'(1) =  Y(2)");
            Console.WriteLine("  Y'(2) = -Y(1)");
            Console.WriteLine("");

            abserr = Math.Sqrt(typeMethods.r8_epsilon());
            relerr = Math.Sqrt(typeMethods.r8_epsilon());

            flag = 1;

            t_start = 0.0;
            t_stop = 2.0 * 3.14159265;

            n_step = 12;

            t = 0.0;
            t_out = 0.0;

            y[0] = 1.0;
            y[1] = 0.0;
            yp = r8_f2(t, y, yp);

            Console.WriteLine("");
            Console.WriteLine("FLAG             T          Y(1)       Y(2)");
            Console.WriteLine("");

            Console.WriteLine(flag.ToString().PadLeft(4) + "  "
                                                         + t.ToString().PadLeft(12) + "  "
                                                         + y[0].ToString().PadLeft(12) + "  "
                                                         + y[1].ToString().PadLeft(12) + "");

            for (i_step = 1; i_step <= n_step; i_step++)
            {
                t = ((double) (n_step - i_step + 1) * t_start
                     + (double) (i_step - 1) * t_stop)
                    / (double) (n_step);

                t_out = ((double) (n_step - i_step) * t_start
                         + (double) (i_step) * t_stop)
                        / (double) (n_step);

                flag = RungeKuttaFehlberg.r8_rkf45(ref data, r8_f2, NEQN, ref y, ref yp, ref t, t_out, ref relerr, abserr,
                    flag);

                Console.WriteLine(flag.ToString().PadLeft(4) + "  "
                                                             + t.ToString().PadLeft(12) + "  "
                                                             + y[0].ToString().PadLeft(12) + "  "
                                                             + y[1].ToString().PadLeft(12) + "");
            }
        }

        static void test06()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST06 solves a scalar ODE using single step mode.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 June 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int NEQN = 1;

            double abserr;
            int flag;
            int i_step;
            int n_step;
            double relerr;
            double t;
            double t_out;
            double t_start;
            double t_stop;
            double[] y = new double[NEQN];
            double[] yp = new double[NEQN];
            RungeKuttaFehlberg.r8RKFData data = new RungeKuttaFehlberg.r8RKFData();

            Console.WriteLine("");
            Console.WriteLine("TEST06");
            Console.WriteLine("  Solve a scalar equation using R8_RKF:");
            Console.WriteLine("");
            Console.WriteLine("  Y' = 0.25 * Y * ( 1 - Y / 20 )");
            Console.WriteLine("");
            Console.WriteLine("  This routine uses the SINGLE STEP mode.");
            Console.WriteLine("");

            abserr = Math.Sqrt(typeMethods.r8_epsilon());
            relerr = Math.Sqrt(typeMethods.r8_epsilon());

            flag = -1;

            t_start = 0.0;
            t_stop = 20.0;

            n_step = 5;

            t = 0.0;
            t_out = 0.0;
            y[0] = 1.0;
            yp = r8_f1(t, y, yp);

            Console.WriteLine("");
            Console.WriteLine("FLAG             T          Y         Y'        Y_Exact         Error");
            Console.WriteLine("");

            Console.WriteLine(flag.ToString().PadLeft(4) + "  "
                                                         + t.ToString().PadLeft(12) + "  "
                                                         + y[0].ToString().PadLeft(12) + "  "
                                                         + yp[0].ToString().PadLeft(12) + "  "
                                                         + r8_y1x(t).ToString().PadLeft(12) + "  "
                                                         + (y[0] - r8_y1x(t)).ToString().PadLeft(12) + "");
            ;

            for (i_step = 1; i_step <= n_step; i_step++)
            {
                t = ((double) (n_step - i_step + 1) * t_start
                     + (double) (i_step - 1) * t_stop)
                    / (double) (n_step);

                t_out = ((double) (n_step - i_step) * t_start
                         + (double) (i_step) * t_stop)
                        / (double) (n_step);

                while (flag < 0)
                {
                    flag = RungeKuttaFehlberg.r8_rkf45(ref data, r8_f1, NEQN, ref y, ref yp, ref t, t_out, ref relerr, abserr,
                        flag);

                    Console.WriteLine(flag + "  "
                                           + t.ToString().PadLeft(12) + "  "
                                           + y[0].ToString().PadLeft(12) + "  "
                                           + yp[0].ToString().PadLeft(12) + "  "
                                           + r8_y1x(t).ToString().PadLeft(12) + "  "
                                           + (y[0] - r8_y1x(t)).ToString().PadLeft(12) + "");
                }

                flag = -2;
            }
        }

        static float[] r4_f1(float t, float[] y, float[] yp)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R4_F1 evaluates the derivative for the ODE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    26 March 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, float T, the value of the independent variable.
            //
            //    Input, float Y[NEQN], the value of the dependent variable.
            //
            //    Output, float YP[NEQN], the value of the derivative dY(1:NEQN)/dT.
            //
        {
            yp[0] = (float) (0.25 * y[0] * (1.0 - y[0] / 20.0));

            return yp;
        }

        static float r4_y1x(float t)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R4_Y1X evaluates the exact solution of the ODE.
            //
            //  Modified:
            //
            //    26 March 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, float T, the value of the independent variable.
            //
            //    Output, float R4_Y1X, the exact solution.
            //
        {
            float value;

            value = (float) (20.0 / (1.0 + 19.0 * Math.Exp(-0.25 * t)));

            return value;
        }

        static float[] r4_f2(float t, float[] y, float[] yp)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R4_F2 evaluates the derivative for the ODE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    26 March 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, float T, the value of the independent variable.
            //
            //    Input, float Y(NEQN), the value of the dependent variable.
            //
            //    Output float YP(NEQN), the value of the derivative dY(1:NEQN)/dT.
            //
        {
            yp[0] = y[1];
            yp[1] = -y[0];

            return yp;
        }

        static double[] r8_f1(double t, double[] y, double[] yp)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_F1 evaluates the derivative for the ODE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    26 March 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double T, the value of the independent variable.
            //
            //    Input, double Y[NEQN], the value of the dependent variable.
            //
            //    Output, double YP[NEQN], the value of the derivative dY(1:NEQN)/dT.
            //
        {
            yp[0] = 0.25 * y[0] * (1.0 - y[0] / 20.0);

            return yp;
        }

        static double r8_y1x(double t)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_Y1X evaluates the exact solution of the ODE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    26 March 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double T, the value of the independent variable.
            //
            //    Output, double R8_Y1X, the exact solution.
            //
        {
            double value;

            value = 20.0 / (1.0 + 19.0 * Math.Exp(-0.25 * t));

            return value;
        }

        static double[] r8_f2(double t, double[] y, double[] yp)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_F2 evaluates the derivative for the ODE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    26 March 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double T, the value of the independent variable.
            //
            //    Input, double Y(NEQN), the value of the dependent variable.
            //
            //    Output double YP(NEQN), the value of the derivative dY(1:NEQN)/dT.
            //
        {
            yp[0] = y[1];
            yp[1] = -y[0];

            return yp;
        }
    }
}