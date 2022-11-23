using System;
using System.Globalization;
using Burkardt.Uniform;

namespace FeynmanKac3DTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for FEYNMAN_KAC_3D.
        //
        //  Discussion:
        //
        //    This program is derived from section 2.5, exercise 2.2 of Petersen and Arbenz.
        //
        //    The problem is to determine the solution U(X,Y,Z) of the following 
        //    partial differential equation:
        //
        //      (1/2) Laplacian U - V(X,Y,Z) * U = 0,
        //
        //    inside the elliptic domain D:
        // 
        //      D = { (X,Y,Z) | (X/A)^2+(Y/B)^2+(Z/C)^2 <= 1 }
        //   
        //    with the boundary condition U(boundary(D)) = 1.
        //
        //    The V(X,Y,Z) is the potential function:
        //
        //      V = 2 * ( (X/A^2)^2 + (Y/B^2)^2 + (Z/C^2)^2 ) + 1/A^2 + 1/B^2 + 1/C^2.
        //
        //    The analytic solution of this problem is already known:
        //
        //      U(X,Y,Z) = exp ( (X/A)^2 + (Y/B)^2 + (Z/C)^2 - 1 ).
        //
        //    Our method is via the Feynman-Kac Formula.
        //
        //    The idea is to start from any (x,y,z) in D, and
        //    compute (x+Wx(t),y+Wy(t),z+Wz(t)) where 3-D Brownian motion
        //    (Wx,Wy,Wz) is updated each step by sqrt(h)*(z1,z2,z3),
        //    each z1,z2,z3 are independent approximately Gaussian 
        //    random variables with zero mean and variance 1. 
        //
        //    Each (x1(t),x2(t),x3(t)) is advanced until (x1,x2,x3) exits 
        //    the domain D.  
        //
        //    Upon its first exit from D, the sample path (x1,x2,x3) is stopped and a 
        //    new sample path at (x,y,z) is started until N such paths are completed.
        // 
        //    The Feynman-Kac formula gives the solution here as
        //
        //      U(X,Y,Z) = (1/N) sum(1 <= I <= N) Y(tau_i),
        //
        //    where
        //
        //      Y(tau) = exp( -int(s=0..tau) v(x1(s),x2(s),x3(s)) ds),
        //
        //    and tau = first exit time for path (x1,x2,x3). 
        //
        //    The integration procedure is a second order weak accurate method:
        //
        //      X(t+h)  = [ x1(t) + sqrt ( h ) * z1 ]
        //                [ x2(t) + sqrt ( h ) * z2 ]
        //                [ x3(t) + sqrt ( h ) * z3 ]
        //
        //    Here Z1, Z2, and Z3 are approximately normal univariate Gaussians. 
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
        //    Original C version by Wesley Petersen.
        //    C++ version by John Burkardt.
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
        const double a = 3.0;
        const double b = 2.0;
        const double c = 1.0;
        const int dim = 3;
        const double h = 0.001;
        int i;
        const int N = 10000;
        int ni;
        int nj;
        int nk;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("FEYNMAN-KAC_3D:");
            
        Console.WriteLine("");
        Console.WriteLine("  Program parameters:");
        Console.WriteLine("");
        Console.WriteLine("  The calculation takes place inside a 3D ellipsoid.");
        Console.WriteLine("  A rectangular grid of points will be defined.");
        Console.WriteLine("  The solution will be estimated for those grid points");
        Console.WriteLine("  that lie inside the ellipsoid.");
        Console.WriteLine("");
        Console.WriteLine("  Each solution will be estimated by computing " + N
                                                                            + " trajectories");
        Console.WriteLine("  from the point to the boundary.");
        Console.WriteLine("");
        Console.WriteLine("    (X/A)^2 + (Y/B)^2 + (Z/C)^2 = 1");
        Console.WriteLine("");
        Console.WriteLine("  The ellipsoid parameters A, B, C are set to:");
        Console.WriteLine("");
        Console.WriteLine("    A = " + a + "");
        Console.WriteLine("    B = " + b + "");
        Console.WriteLine("    C = " + c + "");
        Console.WriteLine("  Stepsize H = " + h + "");
        //
        //  RTH is the scaled stepsize.
        //
        double rth = Math.Sqrt(dim * h);
        //
        //  Choose the spacing so we have about 10 points in the shortest direction.
        //
        if (Math.Abs(a - Math.Min(Math.Min(a, b), c)) <= typeMethods.r8_epsilon())
        {
            ni = 11;
            nj = 1 + (int) Math.Ceiling(b / a) * (ni - 1);
            nk = 1 + (int) Math.Ceiling(c / a) * (ni - 1);
        }
        else if (Math.Abs(b - Math.Min(Math.Min(a, b), c)) <= typeMethods.r8_epsilon())
        {
            nj = 11;
            ni = 1 + (int) Math.Ceiling(a / b) * (nj - 1);
            nk = 1 + (int) Math.Ceiling(c / b) * (nj - 1);
        }
        else
        {
            nk = 11;
            ni = 1 + (int) Math.Ceiling(a / c) * (nk - 1);
            nj = 1 + (int) Math.Ceiling(b / c) * (nk - 1);
        }

        Console.WriteLine("");
        Console.WriteLine("  X coordinate marked by " + ni + " points");
        Console.WriteLine("  Y coordinate marked by " + nj + " points");
        Console.WriteLine("  Z coordinate marked by " + nk + " points");
        //
        //  Loop over the grid points.
        //
        double err = 0.0;
        int n_inside = 0;

        Console.WriteLine("");
        Console.WriteLine("     X           Y           Z          W approx" +
                          "    W exact     Error   Ave Steps");
        Console.WriteLine("");

        for (i = 1; i <= ni; i++)
        {
            double x = ((ni - i) * -a
                        + (i - 1) * a)
                       / (ni - 1);

            int j;
            for (j = 1; j <= nj; j++)
            {
                double y = ((nj - j) * -b
                            + (j - 1) * b)
                           / (nj - 1);

                int k;
                for (k = 1; k <= nk; k++)
                {
                    double z = ((nk - k) * -c
                                + (k - 1) * c)
                               / (nk - 1);

                    double chk = Math.Pow(x / a, 2) + Math.Pow(y / b, 2) + Math.Pow(z / c, 2);

                    double w_exact;
                    double wt;
                    int steps_ave;
                    switch (chk)
                    {
                        case > 1.0:
                            w_exact = 1.0;
                            wt = 1.0;
                            steps_ave = 0;
                            Console.WriteLine("  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                                   + "  " + y.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                                   + "  " + z.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                                   + "  " + wt.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                                   + "  " + w_exact.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                                   + "  " + Math.Abs(w_exact - wt).ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                                   + "  " + steps_ave.ToString().PadLeft(8) + "");
                            continue;
                    }

                    n_inside += 1;
                    //
                    //  Compute the exact solution at this point (x,y,z).
                    //
                    w_exact = Math.Exp(Math.Pow(x / a, 2)
                        + Math.Pow(y / b, 2)
                        + Math.Pow(z / c, 2) - 1.0);
                    //
                    //  Now try to estimate the solution at this point.
                    //
                    wt = 0.0;

                    int steps = 0;
                    // 
                    //  Do N paths
                    //
                    int trial;
                    for (trial = 0; trial < N; trial++)
                    {
                        double x1 = x;
                        double x2 = y;
                        double x3 = z;
                        // 
                        //  W = exp(-int(s=0..t) v(X)ds) 
                        //
                        double w = 1.0;
                        //
                        //  CHK is < 1.0 while the point is inside the ellipsoid.
                        //
                        chk = 0.0;
                        while (chk < 1.0)
                        {
                            //
                            //  Determine DX, DY, DZ.
                            //
                            double ut = UniformRNG.r8_uniform_01(ref seed);
                            double dx;
                            double us;
                            switch (ut)
                            {
                                case < 1.0 / 3.0:
                                {
                                    us = UniformRNG.r8_uniform_01(ref seed) - 0.5;
                                    dx = us switch
                                    {
                                        < 0.0 => -rth,
                                        _ => rth
                                    };

                                    break;
                                }
                                default:
                                    dx = 0.0;
                                    break;
                            }

                            ut = UniformRNG.r8_uniform_01(ref seed);
                            double dy;
                            switch (ut)
                            {
                                case < 1.0 / 3.0:
                                {
                                    us = UniformRNG.r8_uniform_01(ref seed) - 0.5;
                                    dy = us switch
                                    {
                                        < 0.0 => -rth,
                                        _ => rth
                                    };

                                    break;
                                }
                                default:
                                    dy = 0.0;
                                    break;
                            }

                            ut = UniformRNG.r8_uniform_01(ref seed);
                            double dz;
                            switch (ut)
                            {
                                case < 1.0 / 3.0:
                                {
                                    us = UniformRNG.r8_uniform_01(ref seed) - 0.5;
                                    dz = us switch
                                    {
                                        < 0.0 => -rth,
                                        _ => rth
                                    };

                                    break;
                                }
                                default:
                                    dz = 0.0;
                                    break;
                            }

                            double vs = potential(a, b, c, x1, x2, x3);
                            //
                            //  Move to the new point.
                            //
                            x1 += dx;
                            x2 += dy;
                            x3 += dz;

                            steps += 1;

                            double vh = potential(a, b, c, x1, x2, x3);

                            double we = (1.0 - h * vs) * w;
                            w -= 0.5 * h * (vh * we + vs * w);

                            chk = Math.Pow(x1 / a, 2)
                                  + Math.Pow(x2 / b, 2)
                                  + Math.Pow(x3 / c, 2);
                        }

                        wt += w;
                    }

                    //
                    //  WT is the average of the sum of the different trials.
                    //
                    wt /= N;
                    steps_ave = steps / N;
                    //
                    //  Add error in WT to the running L2 error in the solution.
                    //
                    err += Math.Pow(w_exact - wt, 2);

                    Console.WriteLine("  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                           + "  " + y.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                           + "  " + z.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                           + "  " + wt.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                           + "  " + w_exact.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                           + "  " + Math.Abs(w_exact - wt).ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                           + "  " + steps_ave.ToString().PadLeft(8) + "");
                }
            }
        }

        //
        //  Compute the RMS error for all the points.
        //
        err = Math.Sqrt(err / n_inside);

        Console.WriteLine("");
        Console.WriteLine("  RMS absolute error in solution = " + err + "");
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("FEYNMAN_KAC_3D:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static double potential(double a, double b, double c, double x, double y, double z)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POTENTIAL evaluates the potential function V(X,Y,Z).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, C, the parameters that define the ellipse.
        //
        //    Input, double X, Y, Z, the coordinates of the point.
        //
        //    Output, double POTENTIAL, the value of the potential function at (X,Y,Z).
        //
    {
        double value = 0;

        value = 2.0 * (Math.Pow(x / a / a, 2)
                       + Math.Pow(y / b / b, 2)
                       + Math.Pow(z / c / c, 2))
                + 1.0 / a / a
                + 1.0 / b / b
                + 1.0 / c / c;

        return value;
    }


}