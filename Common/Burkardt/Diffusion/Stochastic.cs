using System;
using Burkardt.Types;

namespace Burkardt.Diffusion
{
    public static class Stochastic
    {
        public static double[] diffusivity_1d_pwc(int nc, double[] xc, double[] vc, int np,
        double[] xp )
//****************************************************************************80
//
//  Purpose:
//
//    DIFFUSIVITY_1D_PWC: piecewise constant diffusivity function in 1D.
//
//  Discussion:
//
//    A piecewise constant function is defined over NC intervals, 
//    with interval IC associated with constant value VC(I),
//    separated by NC-1 ascending sorted breakpoints.
//
//    The function is to be evaluated at NP points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 March 2019
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NC, the number of function values.
//
//    Input, double XC[NC-1], the breakpoints, in ascending order.
//
//    Input, double VC[NC], the function values over each interval.
//
//    Input, int NP, the number of evaluation points.
//
//    Input, double XP[NP], the evaluation points.
//
//    Output, double DIFFUSIVITY_1D_PWC[NP], the function value at the 
//    evaluation points.
//
        {
            double[] vp = new double[np];

            for (int ip = 0; ip < np; ip++)
            {
                int kc = 0;
                for (int ic = 0; ic < nc - 1; ic++)
                {
                    if (xp[ip] < xc[ic])
                    {
                        break;
                    }

                    kc = kc + 1;
                }

                vp[ip] = vc[kc];
            }

            return vp;
        }

        public static double[] diffusivity_1d_xk(double dc0, int m, double[] omega, int n,
        double[] x )
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
            int k = 0;
            double w = 1.0;

            double[] dc = new double[n];

            for (int j = 0; j < n; j++)
            {
                dc[j] = 0.0;
            }

            while (k < m)
            {
                if (k < m)
                {
                    k = k + 1;
                    for (int j = 0; j < n; j++)
                    {
                        dc[j] = dc[j] + omega[k - 1] * Math.Sin(w * Math.PI * x[j]);
                    }
                }

                if (k < m)
                {
                    k = k + 1;
                    for (int j = 0; j < n; j++)
                    {
                        dc[j] = dc[j] + omega[k - 1] * Math.Cos(w * Math.PI * x[j]);
                    }
                }

                w = w + 1.0;

            }

            for (int j = 0; j < n; j++)
            {
                dc[j] = Math.Exp(-0.125) * dc[j];
            }

            for (int j = 0; j < n; j++)
            {
                dc[j] = dc0 + Math.Exp(dc[j]);
            }

            return dc;
        }

        public static double[] diffusivity_2d_bnt(double dc0, double[] omega, int n, double[] x,
        double[] y )
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
            double[] arg = new double[n];

            for (int j = 0; j < n; j++)
            {
                arg[j] = omega[0] * Math.Cos(Math.PI * x[j])
                         + omega[1] * Math.Sin(Math.PI * x[j])
                         + omega[2] * Math.Cos(Math.PI * y[j])
                         + omega[3] * Math.Sin(Math.PI * y[j]);
            }

            for (int j = 0; j < n; j++)
            {
                arg[j] = Math.Exp(-0.125) * arg[j];
            }

            double[] dc = new double[n];
            for (int j = 0; j < n; j++)
            {
                dc[j] = dc0 + Math.Exp(arg[j]);
            }

            return dc;
        }

        public static double[] diffusivity_2d_elman(double a, double cl, double dc0, int m_1d,
            double[] omega, int n1, int n2, double[] x, double[] y )
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
//
//  Compute THETA.
//
            double[] theta_1d = theta_solve(a, cl, m_1d);
//
//  Compute LAMBDA_1D.
//
            double[] lambda_1d = new double[m_1d];

            int i = 0;

            for (i = 0; i < m_1d; i++)
            {
                lambda_1d[i] = 2.0 * cl / (1.0 + cl * cl * theta_1d[i] * theta_1d[i]);
            }

//
//  Compute C_1DX(1:M1D) and C_1DY(1:M1D) at (X,Y).
//
            double[] c_1dx = new double[m_1d * n1 * n2];
            double[] c_1dy = new double[m_1d * n1 * n2];

            for (int k = 0; k < n2; k++)
            {
                for (int j = 0; j < n1; j++)
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

                for (int k = 0; k < n2; k++)
                {
                    for (int j = 0; j < n1; j++)
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

                for (int k = 0; k < n2; k++)
                {
                    for (int j = 0; j < n1; j++)
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
            double[] dc = new double[n1 * n2];

            for (int k = 0; k < n2; k++)
            {
                for (int j = 0; j < n1; j++)
                {
                    dc[j + k * n1] = dc0;
                    for (int i2 = 0; i2 < m_1d; i2++)
                    {
                        for (int i1 = 0; i1 < m_1d; i1++)
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
        int n, double[] x, double[] y )
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
            double d = 1.0;
            double lp = Math.Max(d, 2.0 * cl);
            double l = cl / lp;

            double[] dc_arg = new double[n];

            for (int j = 0; j < n; j++)
            {
                dc_arg[j] = 1.0 + omega[0] * Math.Sqrt(Math.Sqrt(Math.PI) * l / 2.0);
            }

            double[] dc = new double[n];
            double[] phi = new double[n];

            for (int i = 2; i <= m; i++)
            {
                double ihalf_r8 = (double) (i / 2);
                double zeta_arg = -Math.Pow(ihalf_r8 * Math.PI * l, 2) / 8.0;
                double zeta = Math.Sqrt(Math.Sqrt(Math.PI) * l) * Math.Exp(zeta_arg);

                if ((i % 2) == 0)
                {
                    for (int j = 0; j < n; j++)
                    {
                        phi[j] = Math.Sin(ihalf_r8 * Math.PI * x[j] / lp);
                    }
                }
                else
                {
                    for (int j = 0; j < n; j++)
                    {
                        phi[j] = Math.Cos(ihalf_r8 * Math.PI * x[j] / lp);
                    }
                }

                for (int j = 0; j < n; j++)
                {
                    dc_arg[j] = dc_arg[j] + zeta * phi[j] * omega[i - 1];
                }
            }

            for (int j = 0; j < n; j++)
            {
                dc[j] = dc0 + Math.Exp(dc_arg[j]);
            }
            return dc;
        }

        public static double[] diffusivity_2d_pwc(int h, int w, double a, double b, double c,
            double d, double[] omega, int n, double[] x, double[] y )
//****************************************************************************80
//
//  Purpose:
//
//    DIFFUSIVITY_2D_PWC: piecewise constant diffusivity function.
//
//  Discussion:
//
//    The 2D stochastic diffusion equation has the form
//
//      - Del ( RHO(X,Y;OMEGA) Del U(X,Y;OMEGA) ) = F(X,Y).
//
//    Here, the diffusivity RHO is assumed to be a piecewise constant function,
//    defined on rectangular grid over [A,B]x[C,D] that is H cells high
//    and W cells wide.  The parameters OMEGA(H*W) are assumed to be given 
//    positive values that are bounded away from 0.
//
//    The underlying grid is assumed to be equally spaced, so each cell
//    is (D-C)/M units high and (B-A)/N wide.
//
//       ^  ^  +-------+-------+-------+-------+
//       D  H  | (H,1) | (H,2) | ----- | (H,W) |
//       |  |  +-------+-------+-------+-------+
//             | ----- | ----- | (I,J) | ----- |
//       Y  I  +-------+-------+-------+-------+
//             | (2,1) | (2,2) | ----- | (2,W) |
//       |  |  +-------+-------+-------+-------+
//       |  |  | (1,1) | (1,2) | ----- | (1,W) |
//       C  1  +-------+-------+-------+-------+
//       |  |
//       +--+--1-------------- J --------------W->
//       +--+--A-------------- X --------------B->
//
//    Indexing is tricky, since our (X,Y) coordinate system indexing and
//    conventions for indexing an array and for ordering indices do not
//    match up.  An arbitrary choice must be made, and here we associate
//    "width" and "X" and "J", as well as "height" and "Y" and "I".
//    We start in the lower left corner, and proceed right first, and
//    then move up one row.
//
//    The following coordinate systems are available:
//
//    * cell index (I,J):  (1,1) <= (I,J) <= (H,W). 
//      I refers to "height" and J to "width".
//
//    * cell counter K:    1     <= K     <= H*W
//      Count from lower left to lower right, then
//      go up one row.
//
//    * physical cell center coordinates (XC,YC)
//      the (X,Y) coordinates of the centers of a cell
//
//    * physical coordinates (X,Y) in [A,B]x[C,D].
//      the (X,Y) coordinates of any point
//      
//    * normalized coordinates (X01,Y01) in [0,1]x[0,1].
//   
//    Transformations from one coordinate system to another:
//
//    * (I,J) --> (XC,YC)
//                X01 = (2*J-1)/2/W
//                Y01 = (2*I-1)/2/H
//                XC = A + (B-A) * X01
//                YC = C + (D-C) * Y01
//
//    * (I,J) --> K
//                K = (I-1)*W+J
//
//    * K     --> (I,J)
//                I = floor ( K - 1 ) / W ) + 1
//                J = mod ( K - 1, W ) + 1
//
//    * K     --> (XC,YC)
//                I = floor ( K - 1 ) / W ) + 1
//                J = mod ( K - 1, W ) + 1
//                X01 = (2*J-1)/2/W
//                Y01 = (2*I-1)/2/H
//                XC = A + (B-A) * X01
//                YC = C + (D-C) * Y01
//
//    * (X,Y) --> (I,J)
//                X01 = ( X - A ) / ( B - A )
//                Y01 = ( Y - C ) / ( D - C )
//                I = round ( ( 2 * h * Y01 + 1 ) / 2 );
//                J = round ( ( 2 * w * X01 + 1 ) / 2 );
//
//    * (X,Y) --> K
//                X01 = ( X - A ) / ( B - A )
//                Y01 = ( Y - C ) / ( D - C )
//                I = round ( ( 2 * h * Y01 + 1 ) / 2 );
//                J = round ( ( 2 * w * X01 + 1 ) / 2 );
//                K = (I-1)*W+J
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2019
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int H, W, specifies the number of elements.
//
//    Input, double A, B, the lower and upper limits of X.
//
//    Input, double C, D, the lower and upper limits of Y.
//
//    Input, double OMEGA[H*W], the parameters.
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the points where the diffusion coefficient is to 
//    be evaluated.
//
//    Output, double DIFFUSIVITY_2D_PWC[N], the value of the diffusion 
//    coefficient at each point.
//
        {
            double[] rho = new double[n];

            for (int ij = 0; ij < n; ij++)
            {
                double x01 = (x[ij] - a) / (b - a);
                double y01 = (y[ij] - c) / (d - c);

                int i = (int) Math.Round((2 * h * y01 - 1) / 2);
                i = Math.Max(i, 0);
                i = Math.Min(i, h - 1);

                int j = (int) Math.Round((2 * w * x01 - 1) / 2);
                j = Math.Max(j, 0);
                j = Math.Min(j, w - 1);

                int k = i * w + j;
                rho[ij] = omega[k];
            }

            return rho;
        }



        public static double[] theta_solve(double a, double cl, int m)
//****************************************************************************80
//
//  Purpose:
//
//    THETA_SOLVE solves a pair of transcendental equations.
//
//  Discussion:
//
//    The vector THETA returned by this function is needed in order to define
//    the terms in a Karhunen-Loeve expansion of a diffusion coefficient.
//
//    The two equations are:
//
//      1/CL - THETA * TAN ( A * THETA ) = 0
//      THETA - 1/CL * TAN ( A * THETA ) = 0
//
//    A and CL are taken to be positive.  Over each open interval 
//
//      ( n - 1/2 pi, n + 1/2 pi ) / A, for N = 0, 1, ...
//
//    the function TAN ( A * THETA ) monotonically rises from -oo to +00; 
//    therefore, it can be shown that there is one root of each equation 
//    in every interval of this form.  Moreover, because of the positivity
//    of A and CL, we can restrict our search to the interval 
//
//      [ n pi, n + 1/2 pi ) / A, for N = 0, 1, ...
//
//    This function computes K such roots, starting in the first interval,
//    finding those two roots, moving to the next interval, and so on, until
//    the requested number of roots have been found.  Odd index roots will
//    correspond to the first equation, and even index roots to the second.
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
//    Howard Elman, Darran Furnival,
//    Solving the Stochastic Steady-State Diffusion Problem Using Multigrid,
//    University of Maryland Department of Computer Science,
//    Technical Report TR-4786.
//
//  Parameters:
//
//    Input, double A, the "radius" of the domain, D = (-A,A)x(-A,A).
//    0 < A.
//
//    Input, double CL, the correlation length.
//    0 < CL.
//
//    Input, int M, the number of values to compute.
//
//    Output, double THETA_SOLVE[M], the values of Theta.
//
        {
            double xc = 0;

            int k;

            double[] theta = new double[m];
            for (k = 0; k < m; k++)
            {
                theta[k] = 0.0;
            }

//
//  [ XA_INIT, XB_INIT] = [ n * pi, n+1/2 pi ] / a, n = 0, 1, 2, ...
//
            double xa_init = 0.0;
            double xb_init = (Math.PI / 2.0) / a;

            k = 0;
            for (;;)
            {
//
//  Seek root of equation 1 in interval.
//
                if (m <= k)
                {
                    break;
                }

                k = k + 1;
                double xa = xa_init;
                double fa = 1.0 / cl - xa * Math.Tan(a * xa);
                double ftol = typeMethods.r8_epsilon() * (Math.Abs(fa) + 1.0);
                double xb = xb_init;
                double fc = fa;
                double bmatol = 100.0 * typeMethods.r8_epsilon() * (Math.Abs(xa) + Math.Abs(xb));

                while (bmatol < xb - xa)
                {
                    xc = (xa + xb) / 2.0;
                    fc = 1.0 / cl - xc * Math.Tan(a * xc);

                    if (Math.Abs(fc) <= ftol)
                    {
                        break;
                    }
                    else if (0.0 < fc)
                    {
                        xa = xc;
                    }
                    else
                    {
                        xb = xc;
                    }
                }

                theta[k - 1] = xc;
//
//  Seek root of equation 2 in interval.
//
                if (m <= k)
                {
                    break;
                }

                k = k + 1;
//
//  In the first interval, we need to skip the zero root of equation 2.
//
                if (k == 2)
                {
                    k = k - 1;
                }
                else
                {
                    xa = xa_init;
                    fa = xa - Math.Tan(a * xa) / cl;
                    ftol = typeMethods.r8_epsilon() * (Math.Abs(fa) + 1.0);
                    xb = xb_init;

                    while (bmatol < xb - xa)
                    {
                        xc = (xa + xb) / 2.0;
                        fc = xc - Math.Tan(a * xc) / cl;

                        if (Math.Abs(fc) <= ftol)
                        {
                            break;
                        }
                        else if (0.0 < fc)
                        {
                            xa = xc;
                        }
                        else
                        {
                            xb = xc;
                        }
                    }

                    theta[k - 1] = xc;
                }

//
//  Advance the interval.
//
                xa_init = xa_init + Math.PI / a;
                xb_init = xb_init + Math.PI / a;
            }

            return theta;
        }
    }
}