using System;
using Burkardt.Uniform;

namespace FeynmanKac1DTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for FEYNMAN_KAC_1D.
        //
        //  Discussion:
        //
        //    This program is derived from section 2.5, exercise 2.2 of Petersen 
        //    and Arbenz.
        //
        //    The problem is to determine the solution U(X) of the following 
        //    partial differential equation:
        //
        //      (1/2) Laplacian U - V(X) * U = 0,
        //
        //    inside the domain D:
        //
        //      D = { X | (X/A)^2 <= 1 }
        // 
        //    with the boundary condition U(boundary(D)) = 1.
        //
        //    V(X) is the potential function:
        //
        //      V = 2 * ( (X/A^2)^2 + 1/A^2.
        //
        //    The analytic solution of this problem is already known:
        //
        //      U(X) = exp ( (X/A)^2 - 1 ).
        //
        //    Our method is via the Feynman-Kac Formula.
        //
        //    The idea is to start from any x in D, and
        //    compute x+Wx(t) where 1D Brownian motion
        //    Wx is updated each step by sqrt(h)*z,
        //    each z is an independent approximately Gaussian 
        //    random variable with zero mean and variance 1. 
        //
        //    Each x1(t) is advanced until x1(t) exits the domain D.  
        //
        //    Upon its first exit from D, the sample path x1 is stopped and a 
        //    new sample path at x is started until N such paths are completed.
        //
        //    The Feynman-Kac formula gives the solution here as
        //
        //      U(X) = (1/N) sum(1 <= I <= N) Y(tau_i),
        //
        //    where
        //
        //      Y(tau) = exp( -int(s=0..tau) v(x1(s)) ds),
        //
        //    and tau = first exit time for path x1. 
        //
        //    The integration procedure is a second order weak accurate method:
        //
        //      X(t+h)  = x1(t) + sqrt ( h ) * z
        //
        //    Here Z is an approximately normal univariate Gaussian. 
        //
        //    An Euler predictor approximates Y at the end of the step
        //
        //      Y_e     = (1 - h*v(X(t)) * Y(t), 
        //
        //    A trapezoidal rule completes the step:
        //
        //      Y(t+h)  = Y(t) - (h/2)*[v(X(t+h))*Y_e + v(X(t))*Y(t)].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 May 2012
        //
        //  Author:
        //
        //    Original C 3D version by Wesley Petersen.
        //    C++ 1D version by John Burkardt.
        //
        //  Reference:
        //
        //    Peter Arbenz, Wesley Petersen,
        //    Introduction to Parallel Computing:
        //    A Practical Guide with Examples in C,
        //    Oxford, 2004,
        //    ISBN: 0-19-851577-4,
        //    LC: QA76.59.P47.
        //
    {
        double a = 2.0;
        double chk;
        double dx;
        double err;
        double h = 0.0001;
        int i;
        int it;
        int k;
        int n = 1000;
        int n_int;
        int ni;
        double rth;
        int seed = 123456789;
        int steps;
        int steps_ave;
        double test;
        double us;
        double vh;
        double vs;
        double x;
        double x1;
        double w;
        double w_exact;
        double we;
        double wt;

        Console.WriteLine("");
        Console.WriteLine("FEYNMAN_KAC_1D:");
            
        Console.WriteLine("");
        Console.WriteLine("  Program parameters:");
        Console.WriteLine("");
        Console.WriteLine("  The calculation takes place inside an interval.");
        Console.WriteLine("  The solution will be estimated at points");
        Console.WriteLine("  on a regular spaced grid within the interval.");
        Console.WriteLine("  Each solution will be estimated by computing " + n + " trajectories");
        Console.WriteLine("  from the point to the boundary.");
        Console.WriteLine("");
        Console.WriteLine("    (X/A)^2 = 1");
        Console.WriteLine("");
        Console.WriteLine("  The interval parameter A is:");
        Console.WriteLine("");
        Console.WriteLine("    A = " + a + "");
        Console.WriteLine("");
        Console.WriteLine("  Path stepsize H = " + h + "");
        //
        //  Choose the spacing so we have about ni points on or in the interval.
        //
        ni = 21;

        Console.WriteLine("");
        Console.WriteLine("  X coordinate discretized by " + ni + 2 + " points");
        //
        //  RTH is the scaled stepsize.
        //
        rth = Math.Sqrt(h);

        err = 0.0;
        //
        //  Loop over the points.
        //
        Console.WriteLine("");
        Console.WriteLine("     I     K       X           W exact" +
                          "      W Approx        Error      Ave Steps  Test");
        Console.WriteLine("");

        k = 0;
        n_int = 0;

        for (i = 0; i <= ni + 1; i++)
        {
            x = ((ni - i) * -a
                 + (i - 1) * a)
                / (ni - 1);

            k += 1;

            test = a * a - x * x;

            switch (test)
            {
                case < 0.0:
                    w_exact = 1.0;
                    wt = 1.0;
                    steps_ave = 0;
                    Console.WriteLine("  " + i.ToString().PadLeft(4)
                                           + "  " + k.ToString().PadLeft(4)
                                           + "  " + x.ToString().PadLeft(12)
                                           + "  " + w_exact.ToString().PadLeft(12)
                                           + "  " + wt.ToString().PadLeft(12)
                                           + "  " + Math.Abs(w_exact - wt).ToString().PadLeft(12)
                                           + "  " + steps_ave.ToString().PadLeft(8)
                                           + "  " + test.ToString().PadLeft(8) + "");
                    continue;
            }

            n_int += 1;
            //
            //  Compute the exact solution at this point (x,y,z).
            //
            w_exact = Math.Exp(Math.Pow(x / a, 2) - 1.0);
            //
            //  Now try to estimate the solution at this point.
            //
            wt = 0.0;
            steps = 0;

            for (it = 1; it <= n; it++)
            {

                x1 = x;
                // 
                //  W = exp(-int(s=0..t) v(X)ds) 
                //
                w = 1.0;
                //
                //  CHK is < 1.0 while the point is inside the interval.
                //
                chk = 0.0;

                while (chk < 1.0)
                {
                    //
                    //  Determine DX.
                    //
                    us = UniformRNG.r8_uniform_01(ref seed) - 0.5;
                    dx = us switch
                    {
                        < 0.0 => -rth,
                        _ => +rth
                    };

                    vs = potential(a, x1);
                    //
                    //  Move to the new point.
                    //
                    x1 += dx;

                    steps += 1;

                    vh = potential(a, x1);

                    we = (1.0 - h * vs) * w;
                    w -= 0.5 * h * (vh * we + vs * w);

                    chk = Math.Pow(x1 / a, 2);
                }

                wt += w;
            }

            //
            //   WT is the average of the sum of the different trials.
            //
            wt /= n;
            steps_ave = steps / n;
            //
            //  Add error in WT to the running L2 error in the solution.
            //
            err += Math.Pow(w_exact - wt, 2);

            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + k.ToString().PadLeft(4)
                                   + "  " + x.ToString().PadLeft(12)
                                   + "  " + w_exact.ToString().PadLeft(12)
                                   + "  " + wt.ToString().PadLeft(12)
                                   + "  " + Math.Abs(w_exact - wt).ToString().PadLeft(12)
                                   + "  " + steps_ave.ToString().PadLeft(8)
                                   + "  " + test.ToString().PadLeft(8) + "");
        }

        //
        //  Compute the RMS error for all the points.
        //
        err = Math.Sqrt(err / n_int);

        Console.WriteLine("");
        Console.WriteLine("  RMS absolute error in solution = " + err + "");
        Console.WriteLine("");
        Console.WriteLine("FEYNMAN_KAC_1D:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static double potential(double a, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POTENTIAL evaluates the potential function V(X).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 May 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, the parameter that defines the region.
        //
        //    Input, double X, the coordinate of the point.
        //
        //    Output, double POTENTIAL, the value of the potential function.
        //
    {
        double value = 0;

        value = 2.0 * Math.Pow(x / a / a, 2) + 1.0 / a / a;

        return value;
    }
}