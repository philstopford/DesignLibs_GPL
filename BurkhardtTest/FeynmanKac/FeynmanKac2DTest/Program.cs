using System;
using Burkardt.Uniform;

namespace FeynmanKac2DTest;

internal class Program
{
      private static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for FEYNMAN_KAC_2D.
            //
            //  Discussion:
            //
            //    This program is derived from section 2.5, exercise 2.2 of Petersen and Arbenz.
            //
            //    The problem is to determine the solution U(X,Y) of the following 
            //    partial differential equation:
            //
            //      (1/2) Laplacian U - V(X,Y) * U = 0,
            //
            //    inside the elliptic domain D:
            // 
            //      D = { (X,Y) | (X/A)^2+(Y/B)^2 <= 1 }
            //   
            //    with the boundary condition U(boundary(D)) = 1.
            //
            //    The V(X,Y) is the potential function:
            //
            //      V = 2 * ( (X/A^2)^2 + (Y/B^2)^2 ) + 1/A^2 + 1/B^2.
            //
            //    The analytic solution of this problem is already known:
            //
            //      U(X,Y) = exp ( (X/A)^2 + (Y/B)^2 - 1 ).
            //
            //    Our method is via the Feynman-Kac Formula.
            //
            //    The idea is to start from any (x,y) in D, and
            //    compute (x+Wx(t),y+Wy(t)) where 2D Brownian motion
            //    (Wx,Wy) is updated each step by sqrt(h)*(z1,z2),
            //    each z1,z2 are independent approximately Gaussian 
            //    random variables with zero mean and variance 1. 
            //
            //    Each (x1(t),x2(t)) is advanced until (x1,x2) exits 
            //    the domain D.  
            //
            //    Upon its first exit from D, the sample path (x1,x2) is stopped and a 
            //    new sample path at (x,y) is started until N such paths are completed.
            // 
            //    The Feynman-Kac formula gives the solution here as
            //
            //      U(X,Y) = (1/N) sum(1 <= I <= N) Y(tau_i),
            //
            //    where
            //
            //      Y(tau) = exp( -int(s=0..tau) v(x1(s),x2(s)) ds),
            //
            //    and tau = first exit time for path (x1,x2). 
            //
            //    The integration procedure is a second order weak accurate method:
            //
            //      X(t+h)  = [ x1(t) + sqrt ( h ) * z1 ]
            //                [ x2(t) + sqrt ( h ) * z2 ]
            //
            //    Here Z1, Z2 are approximately normal univariate Gaussians. 
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
            //    31 May 2012
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
            double a = 2.0;
            double b = 1.0;
            double chk;
            int dim = 2;
            double dx;
            double dy;
            double err;
            double h = 0.0001;
            int i;
            int j;
            int k;
            int N = 10000;
            int n_inside;
            int ni;
            int nj;
            double rth;
            int seed = 123456789;
            int steps;
            int steps_ave;
            double us;
            double ut;
            double vh;
            double vs;
            double x;
            double x1;
            double x2;
            double y;
            double w;
            double w_exact;
            double we;
            double wt;

            Console.WriteLine("");
            Console.WriteLine("FEYNMAN-KAC_2D:");
                  
            Console.WriteLine("");
            Console.WriteLine("  Program parameters:");
            Console.WriteLine("");
            Console.WriteLine("  The calculation takes place inside a 2D ellipse.");
            Console.WriteLine("  A rectangular grid of points will be defined.");
            Console.WriteLine("  The solution will be estimated for those grid points");
            Console.WriteLine("  that lie inside the ellipse.");
            Console.WriteLine("");
            Console.WriteLine("  Each solution will be estimated by computing " + N
                  + " trajectories");
            Console.WriteLine("  from the point to the boundary.");
            Console.WriteLine("");
            Console.WriteLine("    (X/A)^2 + (Y/B)^2 = 1");
            Console.WriteLine("");
            Console.WriteLine("  The ellipsoid parameters A, B are set to:");
            Console.WriteLine("");
            Console.WriteLine("    A = " + a + "");
            Console.WriteLine("    B = " + b + "");
            Console.WriteLine("  Stepsize H = " + h + "");
            //
            //  RTH is the scaled stepsize.
            //
            rth = Math.Sqrt(dim * h);
            //
            //  Choose the spacing so we have about 10 points in the shorter direction.
            //
            if (a < b)
            {
                  ni = 11;
                  nj = 1 + (int) Math.Ceiling(b / a) * (ni - 1);
            }
            else
            {
                  nj = 11;
                  ni = 1 + (int) Math.Ceiling(a / b) * (nj - 1);
            }

            Console.WriteLine("");
            Console.WriteLine("  X coordinate marked by %d points");
            Console.WriteLine(ni + "");
            Console.WriteLine("  Y coordinate marked by %d points");
            Console.WriteLine(nj + "");
            //
            //  Loop over the grid points.
            //
            err = 0.0;
            n_inside = 0;

            Console.WriteLine("");
            Console.WriteLine("     X        Y        W Approx     W Exact    Error    Ave Steps");
            Console.WriteLine("");

            for (j = 1; j <= nj; j++)
            {
                  x = ((nj - j) * -a
                       + (j - 1) * a)
                      / (nj - 1);

                  for (i = 1; i <= ni; i++)
                  {
                        y = ((ni - i) * -b
                             + (i - 1) * b)
                            / (ni - 1);

                        chk = Math.Pow(x / a, 2) + Math.Pow(y / b, 2);

                        switch (chk)
                        {
                              case > 1.0:
                                    w_exact = 1.0;
                                    wt = 1.0;
                                    steps_ave = 0;
                                    Console.WriteLine("  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                                           + "  " + y.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                                           + "  " + wt.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                                           + "  " + w_exact.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                                           + "  " + Math.Abs(w_exact - wt).ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                                           + "  " + steps_ave.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
                                    continue;
                        }

                        n_inside += 1;
                        //
                        //  Compute the exact solution at this point (x,y).
                        //
                        w_exact = Math.Exp(Math.Pow(x / a, 2)
                              + Math.Pow(y / b, 2) - 1.0);
                        //
                        //  Now try to estimate the solution at this point.
                        //
                        wt = 0.0;

                        steps = 0;
                        // 
                        //  Do N paths
                        //
                        for (k = 0; k < N; k++)
                        {
                              x1 = x;
                              x2 = y;
                              // 
                              //  W = exp(-int(s=0..t) v(X)ds) 
                              //
                              w = 1.0;
                              //
                              //  CHK is < 1.0 while the point is inside the ellipse.
                              //
                              chk = 0.0;
                              while (chk < 1.0)
                              {
                                    //
                                    //  Determine DX, DY.
                                    //
                                    ut = UniformRNG.r8_uniform_01(ref seed);
                                    switch (ut)
                                    {
                                          case < 1.0 / 2.0:
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
                                    switch (ut)
                                    {
                                          case < 1.0 / 2.0:
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

                                    vs = potential(a, b, x1, x2);
                                    //
                                    //  Move to the new point.
                                    //
                                    x1 += dx;
                                    x2 += dy;

                                    steps += 1;

                                    vh = potential(a, b, x1, x2);

                                    we = (1.0 - h * vs) * w;
                                    w -= 0.5 * h * (vh * we + vs * w);

                                    chk = Math.Pow(x1 / a, 2)
                                          + Math.Pow(x2 / b, 2);
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
                                               + "  " + wt.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                               + "  " + w_exact.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                               + "  " + Math.Abs(w_exact - wt).ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                               + "  " + steps_ave.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
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
            Console.WriteLine("FEYNMAN_KAC_2D:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
      }

      private static double potential ( double a, double b, double x, double y )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POTENTIAL evaluates the potential function.
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
            //    Input, double A, B, the parameters that define the ellipse.
            //
            //    Input, double X, Y, the coordinates of the point.
            //
            //    Output, double POTENTIAL, the value of the potential function.
            //
      {
            double value = 0;

            value = 2.0 * ( Math.Pow ( x / a / a, 2 )
                            + Math.Pow ( y / b / b, 2 ) )
                    + 1.0 / a / a
                    + 1.0 / b / b;

            return value;
      }
}