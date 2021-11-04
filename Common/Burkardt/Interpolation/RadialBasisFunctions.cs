using System;
using Burkardt.Types;

namespace Burkardt.Interpolation
{
    public static class RadialBasisFunctions
    {
        public static double[] phi1(int n, double[] r, double r0, double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PHI1 evaluates the multiquadric radial basis function.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 June 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
            //    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
            //    Third Edition,
            //    Cambridge University Press, 2007,
            //    ISBN13: 978-0-521-88068-8,
            //    LC: QA297.N866.
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input, double R[N], the radial separation.
            //    0 < R.
            //
            //    Input, double R0, a scale factor.
            //
            //    Output, double V[N], the value of the radial basis function.
            //
        {
            int i;

            for (i = 0; i < n; i++)
            {
                v[i] = Math.Sqrt(r[i] * r[i] + r0 * r0);
            }

            return v;
        }

        public static double[] phi2(int n, double[] r, double r0, double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PHI2 evaluates the inverse multiquadric radial basis function.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 June 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
            //    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
            //    Third Edition,
            //    Cambridge University Press, 2007,
            //    ISBN13: 978-0-521-88068-8,
            //    LC: QA297.N866.
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input, double R[N], the radial separation.
            //    0 < R.
            //
            //    Input, double R0, a scale factor.
            //
            //    Output, double V[N], the value of the radial basis function.
            //
        {
            int i;

            for (i = 0; i < n; i++)
            {
                v[i] = 1.0 / Math.Sqrt(r[i] * r[i] + r0 * r0);
            }

            return v;
        }

        public static double[] phi3(int n, double[] r, double r0, double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PHI3 evaluates the thin-plate spline radial basis function.
            //
            //  Discussion:
            //
            //    Note that PHI3(R,R0) is negative if R < R0.  Thus, for this basis function,
            //    it may be desirable to choose a value of R0 smaller than any possible R.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 June 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
            //    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
            //    Third Edition,
            //    Cambridge University Press, 2007,
            //    ISBN13: 978-0-521-88068-8,
            //    LC: QA297.N866.
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input, double R[N], the radial separation.
            //    0 < R.
            //
            //    Input, double R0, a scale factor.
            //
            //    Output, double V[N], the value of the radial basis function.
            //
        {
            int i;

            for (i = 0; i < n; i++)
            {
                if (r[i] <= 0.0)
                {
                    v[i] = 0.0;
                }
                else
                {
                    v[i] = r[i] * r[i] * Math.Log(r[i] / r0);
                }
            }

            return v;
        }

        public static double[] phi4(int n, double[] r, double r0, double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PHI4 evaluates the gaussian radial basis function.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 June 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
            //    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
            //    Third Edition,
            //    Cambridge University Press, 2007,
            //    ISBN13: 978-0-521-88068-8,
            //    LC: QA297.N866.
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input, double R[N], the radial separation.
            //    0 < R.
            //
            //    Input, double R0, a scale factor.
            //
            //    Output, double V[N], the value of the radial basis function.
            //
        {
            int i;

            for (i = 0; i < n; i++)
            {
                v[i] = Math.Exp(-0.5 * r[i] * r[i] / r0 / r0);
            }

            return v;
        }

        public static double[] rbf_interp(int m, int nd, double[] xd, double r0,
                Func<int, double[], double, double[], double[]> phi, double[] w,
                int ni, double[] xi)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    RBF_INTERP evaluates a radial basis function interpolant.
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
            //  Reference:
            //
            //    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
            //    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
            //    Third Edition,
            //    Cambridge University Press, 2007,
            //    ISBN13: 978-0-521-88068-8,
            //    LC: QA297.N866.
            //
            //  Parameters:
            //
            //    Input, int M, the spatial dimension.
            //
            //    Input, int ND, the number of data points.
            //
            //    Input, double XD[M*ND], the data points.
            //
            //    Input, double R0, a scale factor.  R0 should be larger than the typical
            //    separation between points, but smaller than the maximum separation.
            //    The value of R0 has a significant effect on the resulting interpolant.
            //
            //    Input, void PHI ( int N, double R[], double R0, double V[] ), a 
            //    function to evaluate the radial basis functions.
            //
            //    Input, double W[ND], the weights, as computed by RBF_WEIGHTS.
            //
            //    Input, int NI, the number of interpolation points.
            //
            //    Input, double XI[M*NI], the interpolation points.
            //
            //    Output, double RBF_INTERP[NI], the interpolated values.
            //
        {
            double[] fi;
            int i;
            int j;
            int k;
            double[] r;
            double[] v;

            fi = new double[ni];
            r = new double[nd];
            v = new double[nd];

            for (i = 0; i < ni; i++)
            {
                for (j = 0; j < nd; j++)
                {
                    r[j] = 0.0;
                    for (k = 0; k < m; k++)
                    {
                        r[j] = r[j] + Math.Pow(xi[k + i * m] - xd[k + j * m], 2);
                    }

                    r[j] = Math.Sqrt(r[j]);
                }

                v = phi(nd, r, r0, v);

                fi[i] = typeMethods.r8vec_dot_product(nd, v, w);
            }


            return fi;
        }

        public static double[] rbf_interp_nd(int m, int nd, double[] xd, double r0,
                Func<int, double[], double, double[], double[]> phi, double[] w,
                int ni, double[] xi)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    RBF_INTERP_ND evaluates a radial basis function interpolant.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 June 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
            //    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
            //    Third Edition,
            //    Cambridge University Press, 2007,
            //    ISBN13: 978-0-521-88068-8,
            //    LC: QA297.N866.
            //
            //  Parameters:
            //
            //    Input, int M, the spatial dimension.
            //
            //    Input, int ND, the number of data points.
            //
            //    Input, double XD[M*ND], the data points.
            //
            //    Input, double R0, a scale factor.  R0 should be larger than the typical
            //    separation between points, but smaller than the maximum separation.
            //    The value of R0 has a significant effect on the resulting interpolant.
            //
            //    Input, void PHI ( int N, double R[], double R0, double V[] ), a 
            //    function to evaluate the radial basis functions.
            //
            //    Input, double W[ND], the weights, as computed by RBF_WEIGHTS.
            //
            //    Input, int NI, the number of interpolation points.
            //
            //    Input, double XI[M*NI], the interpolation points.
            //
            //    Output, double RBF_INTERP_ND[NI], the interpolated values.
            //
        {
            double[] fi;
            int i;
            int j;
            int k;
            double[] r;
            double[] v;

            fi = new double[ni];
            r = new double[nd];
            v = new double[nd];

            for (i = 0; i < ni; i++)
            {
                for (j = 0; j < nd; j++)
                {
                    r[j] = 0.0;
                    for (k = 0; k < m; k++)
                    {
                        r[j] = r[j] + Math.Pow(xi[k + i * m] - xd[k + j * m], 2);
                    }

                    r[j] = Math.Sqrt(r[j]);
                }

                phi(nd, r, r0, v);

                fi[i] = typeMethods.r8vec_dot_product(nd, v, w);
            }

            return fi;
        }

        public static double[] rbf_weight(int m, int nd, double[] xd, double r0,
                Func<int, double[], double, double[], double[]> phi,
                double[] fd)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    RBF_WEIGHT computes weights for radial basis function interpolation.
            //
            //  Discussion:
            //
            //    We assume that there are N (nonsingular) equations in N unknowns.
            //
            //    However, it should be clear that, if we are willing to do some kind
            //    of least squares calculation, we could allow for singularity,
            //    inconsistency, or underdetermine systems.  This could be associated
            //    with data points that are very close or repeated, a smaller number
            //    of data points than function values, or some other ill-conditioning
            //    of the system arising from a peculiarity in the point spacing.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 June 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
            //    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
            //    Third Edition,
            //    Cambridge University Press, 2007,
            //    ISBN13: 978-0-521-88068-8,
            //    LC: QA297.N866.
            //
            //  Parameters:
            //
            //    Input, int M, the spatial dimension.
            //
            //    Input, int ND, the number of data points.
            //
            //    Input, double XD[M*ND], the data points.
            //
            //    Input, double R0, a scale factor.  R0 should be larger than the typical
            //    separation between points, but smaller than the maximum separation.
            //    The value of R0 has a significant effect on the resulting interpolant.
            //
            //    Input, void PHI ( int N, double R[], double R0, double V[] ), a 
            //    function to evaluate the radial basis functions.
            //
            //    Input, double FD[ND], the function values at the data points.
            //
            //    Output, double RBF_WEIGHT[ND], the weights.
            //
        {
            double[] a;
            int i;
            int j;
            int k;
            double[] r;
            double[] v;
            double[] w;

            a = new double[nd * nd];
            r = new double[nd];
            v = new double[nd];

            for (i = 0; i < nd; i++)
            {
                for (j = 0; j < nd; j++)
                {
                    r[j] = 0.0;
                    for (k = 0; k < m; k++)
                    {
                        r[j] = r[j] + Math.Pow(xd[k + i * m] - xd[k + j * m], 2);
                    }

                    r[j] = Math.Sqrt(r[j]);
                }

                v = phi(nd, r, r0, v);

                for (j = 0; j < nd; j++)
                {
                    a[i + j * nd] = v[j];
                }
            }

            //
            //  Solve for the weights.
            //
            w = typeMethods.r8mat_solve_svd(nd, nd, a, fd);

            return w;
        }
    }
}