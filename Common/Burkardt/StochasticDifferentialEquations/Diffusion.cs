using System;
using Burkardt.SolveNS;

namespace Burkardt.StochasticDifferentialEquations
{
    public static class Diffusion
    {
        public static double[] diffusivity_1d_xk(double dc0, int m, double[] omega, int n,
                double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DIFFUSIVITY_1D_XK evaluates a 1D stochastic diffusivity function.
            //
            //  Discussion:
            //
            //    The 1D diffusion equation has the form
            //
            //      - d/dx ( DC(X) Del U(X) ) = F(X)
            //
            //    where DC(X) is a function called the diffusivity.
            //
            //    In the stochastic version of the problem, the diffusivity function
            //    includes the influence of stochastic parameters:
            //
            //      - d/dx ( DC(X;OMEGA) d/dx U(X) ) = F(X).
            //
            //    In this function, the domain is assumed to be the unit interval [0.1].
            //
            //
            //    For DC0 = 1 and F(X) = 0, with boundary conditions U(0:OMEGA) = 0,
            //    U(1;OMEGA) = 1, the exact solution is
            //
            //    If OMEGA ~= 0:
            //
            //      U(X;OMEGA) = log ( 1 + OMEGA * X ) / log ( 1 + OMEGA )
            //
            //    If OMEGA = 0:
            //
            //      U(X;OMEGA) = X
            //
            //    In the numerical experiments described in the paper, OMEGA was taken
            //    to be a random variable with a Beta, or Uniform, or Gaussian or 
            //    Poisson or Binomial distribution.
            //
            //    For the Gaussian and Poisson distributions, the positivity requirement 
            //    could not be guaranteed, and the experiments were simply made with a 
            //    "small" variance of 0.1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 August 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Dongbin Xiu, George Karniadakis,
            //    Modeling uncertainty in steady state diffusion problems via
            //    generalized polynomial chaos,
            //    Computer Methods in Applied Mechanics and Engineering,
            //    Volume 191, 2002, pages 4927-4948.
            //
            //  Parameters:
            //
            //    Input, double DC0, the constant term in the expansion of the 
            //    diffusion coefficient.
            //
            //    Input, int M, the number of stochastic parameters.
            //
            //    Input, double OMEGA[M], the stochastic parameters.  
            //
            //    Input, int N, the number of evaluation points.
            //
            //    Input, double X[N], the point where the diffusion coefficient 
            //    is to be evaluated.
            //
            //    Output, double DIFFUSIVITY_1D_XK[N], the value of the diffusion coefficient 
            //    at X.
            //
        {
            double[] dc;
            int j;
            int k;
            double w;

            k = 0;
            w = 1.0;

            dc = new double[n];

            for (j = 0; j < n; j++)
            {
                dc[j] = 0.0;
            }

            while (k < m)
            {
                if (k < m)
                {
                    k = k + 1;
                    for (j = 0; j < n; j++)
                    {
                        dc[j] = dc[j] + omega[k - 1] * Math.Sin(w * Math.PI * x[j]);
                    }
                }

                if (k < m)
                {
                    k = k + 1;
                    for (j = 0; j < n; j++)
                    {
                        dc[j] = dc[j] + omega[k - 1] * Math.Cos(w * Math.PI * x[j]);
                    }
                }

                w = w + 1.0;

            }

            for (j = 0; j < n; j++)
            {
                dc[j] = Math.Exp(-0.125) * dc[j];
            }

            for (j = 0; j < n; j++)
            {
                dc[j] = dc0 + Math.Exp(dc[j]);
            }

            return dc;
        }

        public static double[] diffusivity_2d_bnt(double dc0, double[] omega, int n, double[] x,
                double[] y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DIFFUSIVITY_2D_BNT evaluates a 2D stochastic diffusivity function.
            //
            //  Discussion:
            //
            //    The 2D diffusion equation has the form
            //
            //      - Del ( DC(X,Y) Del U(X,Y) ) = F(X,Y)
            //
            //    where DC(X,Y) is a function called the diffusivity.
            //
            //    In the stochastic version of the problem, the diffusivity function
            //    includes the influence of stochastic parameters:
            //
            //      - Del ( DC(X,Y;OMEGA) Del U(X,Y;OMEGA) ) = F(X,Y).
            //
            //    In this function, the domain is the rectangle [-1.5,0]x[-0.4,0.8].
            //
            //    The four stochastic parameters OMEGA(1:4) are assumed to be independent
            //    identically distributed random variables with mean value zero and 
            //    variance 1.  The distribution is typically taken to be Gaussian or
            //    uniform.
            //
            //    A collocation approach to this problem would then use the roots of
            //    Hermite or Legendre polynomials.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 August 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Ivo Babuska, Fabio Nobile, Raul Tempone,
            //    A stochastic collocation method for elliptic partial differential equations
            //    with random input data,
            //    SIAM Journal on Numerical Analysis,
            //    Volume 45, Number 3, 2007, pages 1005-1034.
            //
            //  Parameters:
            //
            //    Input, double DC0, the constant term in the expansion of the 
            //    diffusion coefficient.  Take DC0 = 10.
            //
            //    Input, double OMEGA[4], the stochastic parameters.
            //
            //    Input, int N, the number of evaluation points.
            //
            //    Input, double X[N], Y[N], the points where the diffusion 
            //    coefficient is to be evaluated.
            //
            //    Output, double DIFFUSIVITY_2D_BNT[N], the value of the diffusion
            //    coefficient at (X,Y).
            //
        {
            double[] arg;
            double[] dc;
            int j;

            arg = new double[n];

            for (j = 0; j < n; j++)
            {
                arg[j] = omega[0] * Math.Cos(Math.PI * x[j])
                         + omega[1] * Math.Sin(Math.PI * x[j])
                         + omega[2] * Math.Cos(Math.PI * y[j])
                         + omega[3] * Math.Sin(Math.PI * y[j]);
            }

            for (j = 0; j < n; j++)
            {
                arg[j] = Math.Exp(-0.125) * arg[j];
            }

            dc = new double[n];
            for (j = 0; j < n; j++)
            {
                dc[j] = dc0 + Math.Exp(arg[j]);
            }

            return dc;
        }

        public static double[] diffusivity_2d_elman(double a, double cl, double dc0, int m_1d,
                double[] omega, int n1, int n2, double[] x, double[] y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DIFFUSIVITY_2D_ELMAN evaluates a 2D stochastic diffusivity function.
            //
            //  Discussion:
            //
            //    The 2D diffusion equation has the form
            //
            //      - Del ( DC(X,Y) Del U(X,Y) ) = F(X,Y)
            //
            //    where DC(X,Y) is a function called the diffusivity.
            //
            //    In the stochastic version of the problem, the diffusivity function
            //    includes the influence of stochastic parameters:
            //
            //      - Del ( DC(X,Y;OMEGA) Del U(X,Y;OMEGA) ) = F(X,Y).
            //
            //    In this function, the domain is assumed to be the square [-A,+A]x[-A,+A].
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 August 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Howard Elman, Darran Furnaval,
            //    Solving the stochastic steady-state diffusion problem using multigrid,
            //    IMA Journal on Numerical Analysis,
            //    Volume 27, Number 4, 2007, pages 675-688.
            //
            //    Roger Ghanem, Pol Spanos,
            //    Stochastic Finite Elements: A Spectral Approach,
            //    Revised Edition,
            //    Dover, 2003,
            //    ISBN: 0486428184,
            //    LC: TA347.F5.G56.
            //
            //  Parameters:
            //
            //    Input, double A, the "radius" of the square region.  The region
            //    is assumed to be [-A,+A]x[-A,+A].
            //    0 < A.
            //
            //    Input, double CL, the correlation length.
            //    0 < CL.
            //
            //    Input, double DC0, the constant term in the expansion of the 
            //    diffusion coefficient.  Take DC0 = 10.
            //
            //    Input, int M_1D, the first and second dimensions of the
            //    stochastic parameter array.
            //
            //    Input, double OMEGA[M_1D*M_1D], the stochastic parameters.
            //
            //    Input, int N1, N2, the dimensions of the X and Y arrays.
            //
            //    Input, double X[N1*N2], Y[N1*N2], the points where the diffusion 
            //    coefficient is to be evaluated.
            //
            //    Output, double DIFFUSIVITY_2D_ELMAN[N1*N2], the value of the diffusion 
            //    coefficient at X.
            //
        {
            double[] c_1dx;
            double[] c_1dy;
            double[] dc;
            int i;
            int i1;
            int i2;
            int j;
            int k;
            double[] lambda_1d;
            double[] theta_1d;
            //
            //  Compute THETA.
            //
            theta_1d = Theta.theta_solve(a, cl, m_1d);
            //
            //  Compute LAMBDA_1D.
            //
            lambda_1d = new double[m_1d];

            for (i = 0; i < m_1d; i++)
            {
                lambda_1d[i] = 2.0 * cl / (1.0 + cl * cl * theta_1d[i] * theta_1d[i]);
            }

            //
            //  Compute C_1DX(1:M1D) and C_1DY(1:M1D) at (X,Y).
            //
            c_1dx = new double[m_1d * n1 * n2];
            c_1dy = new double[m_1d * n1 * n2];

            for (k = 0; k < n2; k++)
            {
                for (j = 0; j < n1; j++)
                {
                    for (i = 0; i < m_1d; i++)
                    {
                        c_1dx[i + j * m_1d + k * m_1d * n1] = 0.0;
                        c_1dy[i + j * m_1d + k * m_1d * n1] = 0.0;
                    }
                }
            }

            i = 0;

            for (;;)
            {
                if (m_1d <= i)
                {
                    break;
                }

                for (k = 0; k < n2; k++)
                {
                    for (j = 0; j < n1; j++)
                    {
                        c_1dx[i + j * m_1d + k * m_1d * n1] = Math.Cos(theta_1d[i] * a * x[j + k * n1])
                                                              / Math.Sqrt(a + Math.Sin(2.0 * theta_1d[i] * a)
                                                                  / (2.0 * theta_1d[i]));

                        c_1dy[i + j * m_1d + k * m_1d * n1] = Math.Cos(theta_1d[i] * a * y[j + k * n1])
                                                              / Math.Sqrt(a + Math.Sin(2.0 * theta_1d[i] * a)
                                                                  / (2.0 * theta_1d[i]));
                    }
                }

                i = i + 1;

                if (m_1d <= i)
                {
                    break;
                }

                for (k = 0; k < n2; k++)
                {
                    for (j = 0; j < n1; j++)
                    {
                        c_1dx[i + j * m_1d + k * m_1d * n1] = Math.Sin(theta_1d[i] * a * x[j + k * n1])
                                                              / Math.Sqrt(a - Math.Sin(2.0 * theta_1d[i] * a)
                                                                  / (2.0 * theta_1d[i]));

                        c_1dy[i + j * m_1d + k * m_1d * n1] = Math.Sin(theta_1d[i] * a * y[j + k * n1])
                                                              / Math.Sqrt(a - Math.Sin(2.0 * theta_1d[i] * a)
                                                                  / (2.0 * theta_1d[i]));
                    }
                }

                i = i + 1;
            }

            //
            //  Evaluate the diffusion coefficient DC at (X,Y).
            //
            dc = new double[n1 * n2];

            for (k = 0; k < n2; k++)
            {
                for (j = 0; j < n1; j++)
                {
                    dc[j + k * n1] = dc0;
                    for (i2 = 0; i2 < m_1d; i2++)
                    {
                        for (i1 = 0; i1 < m_1d; i1++)
                        {
                            dc[j + k * n1] = dc[j + k * n1] + Math.Sqrt(lambda_1d[i1] * lambda_1d[i2])
                                * c_1dx[i1 + j * m_1d + k * m_1d * n1] * c_1dy[i2 + j * m_1d + k * m_1d * n1]
                                * omega[i1 + i2 * m_1d];
                        }
                    }
                }
            }

            return dc;
        }

        public static double[] diffusivity_2d_ntw(double cl, double dc0, int m, double[] omega,
                int n, double[] x, double[] y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DIFFUSIVITY_2D_NTW evaluates a 2D stochastic diffusivity function.
            //
            //  Discussion:
            //
            //    The 2D diffusion equation has the form
            //
            //      - Del ( DC(X,Y) Del U(X,Y) ) = F(X,Y)
            //
            //    where DC(X,Y) is a function called the diffusivity.
            //
            //    In the stochastic version of the problem, the diffusivity function
            //    includes the influence of stochastic parameters:
            //
            //      - Del ( DC(X,Y;OMEGA) Del U(X,Y;OMEGA) ) = F(X,Y).
            //
            //    In this function, the domain is the rectangle [0,D]x[0,D] where D = 1.
            //
            //    Note that in this problem the diffusivity has a one-dimensional
            //    spatial dependence on X, but not on Y
            //
            //    The random variables OMEGA are independent, have zero mean and unit
            //    variance, and are uniformly distributed in [-sqrt(3),+sqrt(3)].
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 August 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Xiang Ma, Nicholas Zabaras,
            //    An adaptive hierarchical sparse grid collocation algorithm for the solution
            //    of stochastic differential equations,
            //    Journal of Computational Physics,
            //    Volume 228, pages 3084-3113, 2009.
            //
            //    Fabio Nobile, Raul Tempone, Clayton Webster,
            //    A Sparse Grid Stochastic Collocation Method for Partial Differential
            //    Equations with Random Input Data,
            //    SIAM Journal on Numerical Analysis,
            //    Volume 46, Number 5, 2008, pages 2309-2345.
            //
            //  Parameters:
            //
            //    Input, double CL, the desired physical correlation length for 
            //    the coefficient.
            //
            //    Input, double DC0, the constant term in the expansion of the 
            //    diffusion coefficient.  Take DC0 = 0.5.
            //
            //    Input, int M, the number of terms in the expansion.
            //
            //    Input, double OMEGA[M], the stochastic parameters.
            //
            //    Input, int N, the number of evaluation points.
            //
            //    Input, double X[N], Y[N], the points where the diffusion 
            //    coefficient is to be evaluated.
            //
            //    Output, double DIFFUSIVITY_2D_NTW[N], the value of the diffusion coefficient
            //    at (X,Y).
            //
        {
            double d;
            double[] dc;
            double[] dc_arg;
            int i;
            double ihalf_r8;
            int j;
            double l;
            double lp;
            double[] phi;
            double zeta;
            double zeta_arg;

            d = 1.0;
            lp = Math.Max(d, 2.0 * cl);
            l = cl / lp;

            dc_arg = new double[n];

            for (j = 0; j < n; j++)
            {
                dc_arg[j] = 1.0 + omega[0] * Math.Sqrt(Math.Sqrt(Math.PI) * l / 2.0);
            }

            dc = new double[n];
            phi = new double[n];

            for (i = 2; i <= m; i++)
            {
                ihalf_r8 = (double) (i / 2);
                zeta_arg = -Math.Pow(ihalf_r8 * Math.PI * l, 2) / 8.0;
                zeta = Math.Sqrt(Math.Sqrt(Math.PI) * l) * Math.Exp(zeta_arg);

                if ((i % 2) == 0)
                {
                    for (j = 0; j < n; j++)
                    {
                        phi[j] = Math.Sin(ihalf_r8 * Math.PI * x[j] / lp);
                    }
                }
                else
                {
                    for (j = 0; j < n; j++)
                    {
                        phi[j] = Math.Cos(ihalf_r8 * Math.PI * x[j] / lp);
                    }
                }

                for (j = 0; j < n; j++)
                {
                    dc_arg[j] = dc_arg[j] + zeta * phi[j] * omega[i - 1];
                }
            }

            for (j = 0; j < n; j++)
            {
                dc[j] = dc0 + Math.Exp(dc_arg[j]);
            }

            return dc;
        }
    }
}