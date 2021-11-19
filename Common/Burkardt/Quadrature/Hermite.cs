using System;
using Burkardt.ClenshawCurtisNS;
using Burkardt.IntegralNS;
using Burkardt.MatrixNS;
using Burkardt.Types;

namespace Burkardt.Quadrature;

using Monomial = Burkardt.MonomialNS.Monomial;
public static class HermiteQuadrature
{
    public static void gen_hermite_compute_points ( int order, double alpha, ref double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEN_HERMITE_COMPUTE_POINTS computes Generalized Hermite quadrature points.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order.
        //
        //    Input, double ALPHA, the exponent of the X factor.
        //    -1.0 < ALPHA.
        //
        //    Output, double X[ORDER], the abscissas.
        //
    {
        double[] w = new double[order];

        gen_hermite_compute ( order, alpha, ref x, ref w );
    }

    public static double[] gen_hermite_compute_points_np ( int order, int np, double[] p, double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEN_HERMITE_COMPUTE_POINTS_NP: Generalized Hermite quadrature points.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order.
        //
        //    Input, int NP, the number of parameters.
        //
        //    Input, double P[NP], contains parameters.
        //    P[0] = ALPHA, the exponent of the X factor. -1.0 < ALPHA.
        //
        //    Output, double X[ORDER], the abscissas.
        //
    {
        double alpha = p[0];

        gen_hermite_compute_points ( order, alpha, ref x );

        return x;
    }
    public static void gen_hermite_compute_weights ( int order, double alpha, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEN_HERMITE_COMPUTE_WEIGHTS computes Generalized Hermite quadrature weights.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order.
        //
        //    Input, double ALPHA, the exponent of the X factor.
        //    -1.0 < ALPHA.
        //
        //    Output, double W[ORDER], the weights.
        //
    {
        double[] x = new double[order];

        gen_hermite_compute ( order, alpha, ref x, ref w );
    }

    public static double[] gen_hermite_compute_weights_np ( int order, int np, double[] p,
            double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEN_HERMITE_COMPUTE_WEIGHTS_NP: Generalized Hermite quadrature weights.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order.
        //
        //    Input, int NP, the number of parameters.
        //
        //    Input, double P[NP], contains parameters.
        //    P[0] = ALPHA, the exponent of the X factor. -1.0 < ALPHA.
        //
        //    Output, double W[ORDER], the weights.
        //
    {
        double alpha = p[0];

        gen_hermite_compute_weights ( order, alpha, ref w );

        return w;
    }
        
    public static void gen_hermite_compute_np ( int order, int np, double[] p, ref double[] x,
            ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEN_HERMITE_COMPUTE_NP computes a Generalized Hermite quadrature rule.
        //
        //  Discussion:
        //
        //    The integral:
        //
        //      Integral ( -oo < x < +oo ) |x|^ALPHA exp(-x^2) f(x) dx
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Philip Davis, Philip Rabinowitz,
        //    Methods of Numerical Integration,
        //    Second Edition,
        //    Dover, 2007,
        //    ISBN: 0486453391,
        //    LC: QA299.3.D28.
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order.
        //    1 <= ORDER.
        //
        //    Input, int NP, the number of parameters.
        //
        //    Input, double P[NP], contains parameters.
        //    P[0] = ALPHA, the exponent of the X factor. -1.0 < ALPHA.
        //
        //    Output, double X[ORDER], the abscissas.
        //
        //    Output, double W[ORDER], the weights.
        //
    {
        double alpha = p[0];

        gen_hermite_compute ( order, alpha, ref x, ref w );
    }
        
    public static void hermite_compute_np ( int order, int np, double[] p, ref double[] x,
            ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_COMPUTE_NP computes a Hermite quadrature rule.
        //
        //  Discussion:
        //
        //    The abscissas are the zeros of the N-th order Hermite polynomial.
        //
        //    The integral:
        //
        //      Integral ( -oo < X < +oo ) exp ( - X * X ) * F(X) dX
        //
        //    The quadrature rule:
        //
        //      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 June 2009
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Arthur Stroud, Don Secrest,
        //    Gaussian Quadrature Formulas,
        //    Prentice Hall, 1966,
        //    LC: QA299.4G3S7.
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order.
        //    1 <= ORDER.
        //
        //    Input, int NP, the number of parameters.
        //
        //    Input, double P[NP], parameters which are not needed by this function.
        //
        //    Output, double X[ORDER], the abscissas.
        //
        //    Output, double W[ORDER], the weights.
        //
    {
        hermite_compute ( order, ref x, ref w );
    }
        
    public static void hermite_compute_weights ( int order, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_COMPUTE_WEIGHTS computes Hermite quadrature weights.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order.
        //
        //    Output, double W[ORDER], the weights.
        //
    {
        double[] x = new double[order];

        hermite_compute ( order, ref x, ref w );
    }

    public static double[] hermite_compute_weights_np ( int order, int np, double[] p, double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_COMPUTE_WEIGHTS_NP computes Hermite quadrature weights.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order.
        //
        //    Input, int NP, the number of parameters.
        //
        //    Input, double P[NP], parameters which are not needed by this function.
        //
        //    Output, double W[ORDER], the weights.
        //
    {
        hermite_compute_weights ( order, ref w );

        return w;
    }
        
    public static void hermite_compute ( int n, ref double[] x, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_COMPUTE computes a Gauss-Hermite quadrature rule.
        //
        //  Discussion:
        //
        //    The code uses an algorithm by Elhay and Kautsky.
        //
        //    The abscissas are the zeros of the N-th order Hermite polynomial.
        //
        //    The integral:
        //
        //      integral ( -oo < x < +oo ) exp ( - x * x ) * f(x) dx
        //
        //    The quadrature rule:
        //
        //      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 May 2012
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Sylvan Elhay, Jaroslav Kautsky,
        //    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
        //    Interpolatory Quadrature,
        //    ACM Transactions on Mathematical Software,
        //    Volume 13, Number 4, December 1987, pages 399-415.
        //
        //  Parameters:
        //
        //    Input, int N, the number of abscissas.
        //
        //    Output, double X[N], the abscissas.
        //
        //    Output, double W[N], the weights.
        //
    {
        int i;
        //
        //  Define the zero-th moment.
        //
        const double arg = 0.5;
        double zemu = typeMethods.r8_gamma ( arg );
        //
        //  Define the Jacobi matrix.
        //
        double[] bj = new double[n];

        for ( i = 0; i < n; i++ )
        {
            bj[i] = Math.Sqrt ( (i + 1) / 2.0 );
        }

        for ( i = 0; i < n; i++ )
        {
            x[i] = 0.0;
        }

        w[0] = Math.Sqrt ( zemu );
        for ( i = 1; i < n; i++ )
        {
            w[i] = 0.0;
        }
        //
        //  Diagonalize the Jacobi matrix.
        //
        IMTQLX.imtqlx ( n, ref x, ref bj, ref w );
        x[(n - 1) / 2] = (n % 2) switch
        {
            //
            //  If N is odd, force the middle X to be exactly zero.
            //
            1 => 0.0,
            _ => x[(n - 1) / 2]
        };

        for ( i = 0; i < n; i++ )
        {
            w[i] *= w[i];
        }
    }
        
    public static void gen_hermite_compute ( int n, double alpha, ref double[] x, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEN_HERMITE_COMPUTE computes a generalized Gauss-Hermite quadrature rule.
        //
        //  Discussion:
        //
        //    The code uses an algorithm by Elhay and Kautsky.
        //
        //    The integral:
        //
        //      integral ( -oo < x < +oo ) |x|^alpha exp(-x^2) f(x) dx
        //
        //    The quadrature rule:
        //
        //      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 April 2011
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Sylvan Elhay, Jaroslav Kautsky,
        //    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
        //    Interpolatory Quadrature,
        //    ACM Transactions on Mathematical Software,
        //    Volume 13, Number 4, December 1987, pages 399-415.
        //
        //  Parameters:
        //
        //    Input, int N, the number of abscissas.
        //
        //    Input, double ALPHA, the parameter.
        //    -1.0 < ALPHA.
        //
        //    Output, double X[N], the abscissas.
        //
        //    Output, double W[N], the weights.
        //
    {
        int i;
        //
        //  Define the zero-th moment.
        //
        double zemu = typeMethods.r8_gamma ( ( alpha + 1.0 ) / 2.0 );
        //
        //  Define the Jacobi matrix.
        //
        double[] bj = new double[n];

        for ( i = 0; i < n; i++ )
        {
            double i_r8 = i + 1;
            bj[i] = (i % 2) switch
            {
                0 => (i_r8 + alpha) / 2.0,
                _ => i_r8 / 2.0
            };
        }

        for ( i = 0; i < n; i++ )
        {
            bj[i] = Math.Sqrt ( bj[i] );
        }

        for ( i = 0; i < n; i++ )
        {
            x[i] = 0.0;
        }

        w[0] = Math.Sqrt ( zemu );
        for ( i = 1; i < n; i++ )
        {
            w[i] = 0.0;
        }
        //
        //  Diagonalize the Jacobi matrix.
        //
        IMTQLX.imtqlx ( n, ref x, ref bj, ref w );

        for ( i = 0; i < n; i++ )
        {
            w[i] *= w[i];
        }
    }
        
    public static double monomial_quadrature_gen_hermite ( int expon, double alpha, int order, 
            int option, double[] w, double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MONOMIAL_QUADRATURE_GEN_HERMITE applies a quadrature rule to a monomial.
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
        //    Input, int EXPON, the exponent.
        //
        //    Input, double ALPHA, the power of |X| in the weighting function.
        //
        //    Input, int ORDER, the number of points in the rule.
        //
        //    Input, int OPTION, indicates standard or modified rule.
        //    0, standard Gauss-Hermite rule for integrand |x|^alpha exp(-x^2)*f(x).
        //    1, modified Gauss-Hermite rule for integrand                     f(x).
        //
        //    Input, double W[ORDER], the quadrature weights.
        //
        //    Input, double X[ORDER], the quadrature points.
        //
        //    Output, double MONOMIAL_QUADRATURE_GEN_HERMITE, the quadrature error.
        //
    {
        int i;
        //
        //  Get the exact value of the integral of the monomial.
        //
        double exact = Integral.gen_hermite_integral ( expon, alpha );
        //
        //  Evaluate the unweighted monomial at the quadrature points.
        //
        double quad = 0.0;
        switch (option)
        {
            case 0:
            {
                for ( i = 0; i < order; i++ )
                {
                    quad += w[i] * Math.Pow ( x[i], expon );
                }

                break;
            }
            default:
            {
                for ( i = 0; i < order; i++ )
                {
                    quad += w[i] * Math.Pow ( Math.Abs ( x[i] ), alpha )
                                 * Math.Exp ( - x[i] * x[i] ) * Math.Pow ( x[i], expon );
                }

                break;
            }
        }

        double quad_error = exact switch
        {
            //
            //  Error:
            //
            0.0 => Math.Abs(quad),
            _ => Math.Abs((quad - exact) / exact)
        };
        return quad_error;
    }
        
    public static void hermite_ek_compute ( int n, ref double[] x, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_EK_COMPUTE computes a Gauss-Hermite quadrature rule.
        //
        //  Discussion:
        //
        //    The code uses an algorithm by Elhay and Kautsky.
        //
        //    The abscissas are the zeros of the N-th order Hermite polynomial.
        //
        //    The integral:
        //
        //      integral ( -oo < x < +oo ) exp ( - x * x ) * f(x) dx
        //
        //    The quadrature rule:
        //
        //      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 April 2011
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Sylvan Elhay, Jaroslav Kautsky,
        //    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
        //    Interpolatory Quadrature,
        //    ACM Transactions on Mathematical Software,
        //    Volume 13, Number 4, December 1987, pages 399-415.
        //
        //  Parameters:
        //
        //    Input, int N, the number of abscissas.
        //
        //    Output, double X[N], the abscissas.
        //
        //    Output, double W[N], the weights.
        //
    {
        int i;
        //
        //  Define the zero-th moment.
        //
        double arg = 0.5;
        double zemu = Helpers.Gamma ( arg );
        //
        //  Define the Jacobi matrix.
        //
        double[] bj = new double[n];

        for ( i = 0; i < n; i++ )
        {
            bj[i] = Math.Sqrt ( (i + 1) / 2.0 );
        }

        for ( i = 0; i < n; i++ )
        {
            x[i] = 0.0;
        }

        w[0] = Math.Sqrt ( zemu );
        for ( i = 1; i < n; i++ )
        {
            w[i] = 0.0;
        }
        //
        //  Diagonalize the Jacobi matrix.
        //
        IMTQLX.imtqlx ( n, ref x, ref bj, ref w );

        for ( i = 0; i < n; i++ )
        {
            w[i] *= w[i];
        }
    }
        
    public static void hermite_abscissa(int dim_num, int point_num, int[] grid_index,
            int[] grid_base, ref double[] grid_point, int gridIndex = 0, int gridBaseIndex = 0, int gridPointIndex = 0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_ABSCISSA sets abscissas for multidimensional Gauss-Hermite quadrature.
        //
        //  Discussion:
        //
        //    The "nesting" as it occurs for Gauss-Hermite sparse grids simply
        //    involves the use of a specified set of permissible orders for the
        //    rule.  
        //
        //    The X array lists the (complete) Gauss-Legendre abscissas for rules 
        //    of order 1, 3, 7, 15, 31, 63 or 127, in order. 
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 October 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int POINT_NUM, the number of points.
        //
        //    Input, int GRID_INDEX[DIM_NUM*POINT_NUM], the index of the abscissa
        //    from the rule, for each dimension and point.
        //
        //    Input, int GRID_BASE[DIM_NUM], the number of points used in the 
        //    rule for a given dimension.
        //
        //    Output, double GRID_POINT[DIM_NUM*POINT_NUM], the grid points of abscissas.
        //
    {
        int dim;
        int point;
        int[] skip =  {
                0, 1, 4, 11, 26, 57, 120, 247
            }
            ;
        double[] x =  {
                0.0E+00,
                -0.122474487139158904909864203735E+01,
                0.0E+00,
                0.122474487139158904909864203735E+01,
                -0.265196135683523349244708200652E+01,
                -0.167355162876747144503180139830E+01,
                -0.816287882858964663038710959027E+00,
                0.0E+00,
                0.816287882858964663038710959027E+00,
                0.167355162876747144503180139830E+01,
                0.265196135683523349244708200652E+01,
                -0.449999070730939155366438053053E+01,
                -0.366995037340445253472922383312E+01,
                -0.296716692790560324848896036355E+01,
                -0.232573248617385774545404479449E+01,
                -0.171999257518648893241583152515E+01,
                -0.113611558521092066631913490556E+01,
                -0.565069583255575748526020337198E+00,
                0.0E+00,
                0.565069583255575748526020337198E+00,
                0.113611558521092066631913490556E+01,
                0.171999257518648893241583152515E+01,
                0.232573248617385774545404479449E+01,
                0.296716692790560324848896036355E+01,
                0.366995037340445253472922383312E+01,
                0.449999070730939155366438053053E+01,
                -6.9956801237185402753248521473232E+00,
                -6.2750787049428601427036567812530E+00,
                -5.6739614446185883296332558789276E+00,
                -5.1335955771123807045862968913996E+00,
                -4.6315595063128599420667997654336E+00,
                -4.1562717558181451724831352315314E+00,
                -3.7007434032314694224497164589673E+00,
                -3.2603207323135408104645401509648E+00,
                -2.8316804533902054557015640151425E+00,
                -2.4123177054804201051740184582119E+00,
                -2.0002585489356389657975562598571E+00,
                -1.5938858604721398261388419455550E+00,
                -1.1918269983500464260821358649242E+00,
                -0.79287697691530893968593032998830E+00,
                -0.39594273647142311094670041663436E+00,
                0.0000000000000000000000000000000E+00,
                0.39594273647142311094670041663436E+00,
                0.79287697691530893968593032998830E+00,
                1.1918269983500464260821358649242E+00,
                1.5938858604721398261388419455550E+00,
                2.0002585489356389657975562598571E+00,
                2.4123177054804201051740184582119E+00,
                2.8316804533902054557015640151425E+00,
                3.2603207323135408104645401509648E+00,
                3.7007434032314694224497164589673E+00,
                4.1562717558181451724831352315314E+00,
                4.6315595063128599420667997654336E+00,
                5.1335955771123807045862968913996E+00,
                5.6739614446185883296332558789276E+00,
                6.2750787049428601427036567812530E+00,
                6.9956801237185402753248521473232E+00,
                -10.435499877854168053468115427285E+00,
                -9.8028759912974963635223935286507E+00,
                -9.2792019543050391319404745506496E+00,
                -8.8118581437284546442526628275570E+00,
                -8.3807683451863219343010651043788E+00,
                -7.9755950801420373181541806298501E+00,
                -7.5901395198641066762479783194468E+00,
                -7.2203167078889678461161324222529E+00,
                -6.8632544331795368527353285876066E+00,
                -6.5168348106821160605273395854042E+00,
                -6.1794379922705969862418461787263E+00,
                -5.8497884000810673462526582961482E+00,
                -5.5268572526403031425047575122840E+00,
                -5.2097979830408354861575136416263E+00,
                -4.8979018644975742350745099214868E+00,
                -4.5905665744435190229271294569091E+00,
                -4.2872733352824404031727616199454E+00,
                -3.9875699104197157485227052068068E+00,
                -3.6910577000963465117322810559754E+00,
                -3.3973817713303911852755941806287E+00,
                -3.1062230279282566329138616746036E+00,
                -2.8172919672837977750747135657355E+00,
                -2.5303236304712010926855221718499E+00,
                -2.2450734604812066298995918179330E+00,
                -1.9613138583081485293922008411321E+00,
                -1.6788312791720137520802800622638E+00,
                -1.3974237486049625107570752063702E+00,
                -1.1168987050996462690510970277840E+00,
                -0.83707109558947615977737795461293E+00,
                -0.55776166427908221668763665253822E+00,
                -0.27879538567115223986687628627202E+00,
                0.00000000000000000000000000000000E+00,
                0.27879538567115223986687628627202E+00,
                0.55776166427908221668763665253822E+00,
                0.83707109558947615977737795461293E+00,
                1.1168987050996462690510970277840E+00,
                1.3974237486049625107570752063702E+00,
                1.6788312791720137520802800622638E+00,
                1.9613138583081485293922008411321E+00,
                2.2450734604812066298995918179330E+00,
                2.5303236304712010926855221718499E+00,
                2.8172919672837977750747135657355E+00,
                3.1062230279282566329138616746036E+00,
                3.3973817713303911852755941806287E+00,
                3.6910577000963465117322810559754E+00,
                3.9875699104197157485227052068068E+00,
                4.2872733352824404031727616199454E+00,
                4.5905665744435190229271294569091E+00,
                4.8979018644975742350745099214868E+00,
                5.2097979830408354861575136416263E+00,
                5.5268572526403031425047575122840E+00,
                5.8497884000810673462526582961482E+00,
                6.1794379922705969862418461787263E+00,
                6.5168348106821160605273395854042E+00,
                6.8632544331795368527353285876066E+00,
                7.2203167078889678461161324222529E+00,
                7.5901395198641066762479783194468E+00,
                7.9755950801420373181541806298501E+00,
                8.3807683451863219343010651043788E+00,
                8.8118581437284546442526628275570E+00,
                9.2792019543050391319404745506496E+00,
                9.8028759912974963635223935286507E+00,
                10.435499877854168053468115427285E+00,
                -15.228338148167350978246954433464E+00,
                -14.669595158833972632746354112896E+00,
                -14.209085995284870755168244250887E+00,
                -13.799722290211676634645246746673E+00,
                -13.423518590070950062438258321855E+00,
                -13.071208660474601901583995439649E+00,
                -12.737235652415686338138003924072E+00,
                -12.417939378869715805445879624069E+00,
                -12.110749020947747600132123508132E+00,
                -11.813772198267727195134584136191E+00,
                -11.525565112572696599167888588564E+00,
                -11.244994583785543445194384194300E+00,
                -10.971150569840247423423040263881E+00,
                -10.703288201027481347670940744690E+00,
                -10.440787957772772867742591798027E+00,
                -10.183127473450343888624126450357E+00,
                -9.9298610495114250736847004273684E+00,
                -9.6806044412474728038150712732737E+00,
                -9.4350233389881650135019598506287E+00,
                -9.1928244988460305715774195052527E+00,
                -8.9537488108565404323807890169970E+00,
                -8.7175658087076307363833999548548E+00,
                -8.4840692689832473326097180339984E+00,
                -8.2530736454457156579694124243888E+00,
                -8.0244111514703375578594739796798E+00,
                -7.7979293513870105420829120455591E+00,
                -7.5734891556083454022834960763301E+00,
                -7.3509631392269052701961258043733E+00,
                -7.1302341220350710668064025713431E+00,
                -6.9111939615465713197465633109366E+00,
                -6.6937425208758294190074417381666E+00,
                -6.4777867811645365448144903821487E+00,
                -6.2632400742737354345609723857092E+00,
                -6.0500214161419845694465474482388E+00,
                -5.8380549248774187386601690807757E+00,
                -5.6272693105464816659423455794909E+00,
                -5.4175974259243240722848425872924E+00,
                -5.2089758693153983587570258372239E+00,
                -5.0013446320386360038520809107373E+00,
                -4.7946467843764925009748509930857E+00,
                -4.5888281947698372951606485031212E+00,
                -4.3838372778464736294253744407459E+00,
                -4.1796247675352031349421189892408E+00,
                -3.9761435120673355916035814195920E+00,
                -3.7733482881250526721004678400057E+00,
                -3.5711956317782180447199756485249E+00,
                -3.3696436841717397896643629240035E+00,
                -3.1686520501953630191857798261495E+00,
                -2.9681816685955910267761649521505E+00,
                -2.7681946921824058801226545958892E+00,
                -2.5686543769473501723144013022363E+00,
                -2.3695249790490401080012474645702E+00,
                -2.1707716587411506879498498083695E+00,
                -1.9723603904195020079324743227565E+00,
                -1.7742578780516791584676442103681E+00,
                -1.5764314753267801315519597621879E+00,
                -1.3788491099261778091441557053728E+00,
                -1.1814792113700685848678583598423E+00,
                -0.98429064194027277726568984213773E+00,
                -0.78725263021825034151596831878971E+00,
                -0.59033470680942102142230439346102E+00,
                -0.39350664185130136568037826200185E+00,
                -0.19673838392423251964272239737078E+00,
                0.0000000000000000000000000000000E+00,
                0.19673838392423251964272239737078E+00,
                0.39350664185130136568037826200185E+00,
                0.59033470680942102142230439346102E+00,
                0.78725263021825034151596831878971E+00,
                0.98429064194027277726568984213773E+00,
                1.1814792113700685848678583598423E+00,
                1.3788491099261778091441557053728E+00,
                1.5764314753267801315519597621879E+00,
                1.7742578780516791584676442103681E+00,
                1.9723603904195020079324743227565E+00,
                2.1707716587411506879498498083695E+00,
                2.3695249790490401080012474645702E+00,
                2.5686543769473501723144013022363E+00,
                2.7681946921824058801226545958892E+00,
                2.9681816685955910267761649521505E+00,
                3.1686520501953630191857798261495E+00,
                3.3696436841717397896643629240035E+00,
                3.5711956317782180447199756485249E+00,
                3.7733482881250526721004678400057E+00,
                3.9761435120673355916035814195920E+00,
                4.1796247675352031349421189892408E+00,
                4.3838372778464736294253744407459E+00,
                4.5888281947698372951606485031212E+00,
                4.7946467843764925009748509930857E+00,
                5.0013446320386360038520809107373E+00,
                5.2089758693153983587570258372239E+00,
                5.4175974259243240722848425872924E+00,
                5.6272693105464816659423455794909E+00,
                5.8380549248774187386601690807757E+00,
                6.0500214161419845694465474482388E+00,
                6.2632400742737354345609723857092E+00,
                6.4777867811645365448144903821487E+00,
                6.6937425208758294190074417381666E+00,
                6.9111939615465713197465633109366E+00,
                7.1302341220350710668064025713431E+00,
                7.3509631392269052701961258043733E+00,
                7.5734891556083454022834960763301E+00,
                7.7979293513870105420829120455591E+00,
                8.0244111514703375578594739796798E+00,
                8.2530736454457156579694124243888E+00,
                8.4840692689832473326097180339984E+00,
                8.7175658087076307363833999548548E+00,
                8.9537488108565404323807890169970E+00,
                9.1928244988460305715774195052527E+00,
                9.4350233389881650135019598506287E+00,
                9.6806044412474728038150712732737E+00,
                9.9298610495114250736847004273684E+00,
                10.183127473450343888624126450357E+00,
                10.440787957772772867742591798027E+00,
                10.703288201027481347670940744690E+00,
                10.971150569840247423423040263881E+00,
                11.244994583785543445194384194300E+00,
                11.525565112572696599167888588564E+00,
                11.813772198267727195134584136191E+00,
                12.110749020947747600132123508132E+00,
                12.417939378869715805445879624069E+00,
                12.737235652415686338138003924072E+00,
                13.071208660474601901583995439649E+00,
                13.423518590070950062438258321855E+00,
                13.799722290211676634645246746673E+00,
                14.209085995284870755168244250887E+00,
                14.669595158833972632746354112896E+00,
                15.228338148167350978246954433464E+00
            }
            ;

        for (dim = 0; dim < dim_num; dim++)
        {
            switch (grid_base[gridBaseIndex + dim])
            {
                case < 0:
                    Console.WriteLine("");
                    Console.WriteLine("HERMITE_ABSCISSA - Fatal error!");
                    Console.WriteLine("  Some base values are less than 0.");
                    return;
            }
        }

        for (dim = 0; dim < dim_num; dim++)
        {
            switch (grid_base[gridBaseIndex + dim])
            {
                case > 63:
                    Console.WriteLine("");
                    Console.WriteLine("HERMITE_ABSCISSA - Fatal error!");
                    Console.WriteLine("  Some base values are greater than 63.");
                    return;
            }
        }

        for (point = 0; point < point_num; point++)
        {
            for (dim = 0; dim < dim_num; dim++)
            {
                int level = (int)Math.Log2(grid_base[gridBaseIndex + dim] + 1);

                int pointer = skip[level] + grid_index[gridIndex + dim + point * dim_num] + grid_base[gridBaseIndex + dim];

                grid_point[gridPointIndex + dim + point * dim_num] = x[pointer];
            }
        }

    }

    public static void hermite_genz_keister_lookup_points(int n, ref double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_GENZ_KEISTER_LOOKUP_POINTS looks up Genz-Keister Hermite abscissas.
        //
        //  Discussion:
        //
        //    The integral:
        //
        //      integral ( -oo <= x <= +oo ) f(x) exp ( - x * x ) dx
        //
        //    The quadrature rule:
        //
        //      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
        //
        //    A nested family of rules for the Hermite integration problem
        //    was produced by Genz and Keister.  The structure of the nested
        //    family was denoted by 1+2+6+10+?, that is, it comprised rules
        //    of successive orders O = 1, 3, 9, 19, and a final rule of order
        //    35, 37, 41 or 43.
        //
        //    The precisions of these rules are P = 1, 5, 15, 29, 
        //    with the final rule of precision 51, 55, 63 or 67.
        //
        //    Three related families begin the same way, but end with a different final
        //    rule.  As a convenience, this function includes these final rules as well:
        //
        //    Designation  Orders       Precisions
        //
        //    1+2+6+10+16, 1,3,9,19,35  1,5,15,29,51
        //    1+2+6+10+18  1,3,9,19,37  1,5,15,29,55
        //    1+2+6+10+22  1,3,9,19,41  1,5,15,29,63
        //    1+2+6+10+24  1,3,9,19,43  1,5,15,29,67
        //
        //    Some of the data in this function was kindly supplied directly by
        //    Alan Genz on 24 April 2011.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Alan Genz, Bradley Keister,
        //    Fully symmetric interpolatory rules for multiple integrals
        //    over infinite regions with Gaussian weight,
        //    Journal of Computational and Applied Mathematics,
        //    Volume 71, 1996, pages 299-309
        //
        //    Florian Heiss, Viktor Winschel,
        //    Likelihood approximation by numerical integration on sparse grids,
        //    Journal of Econometrics,
        //    Volume 144, 2008, pages 62-80.
        //
        //    Thomas Patterson,
        //    The Optimal Addition of Points to Quadrature Formulae,
        //    Mathematics of Computation,
        //    Volume 22, Number 104, October 1968, pages 847-856.
        //
        //  Parameters:
        //
        //    Input, int N, the order.
        //    N must be 1, 3, 9, 19, 35, 37, 41, or 43.
        //
        //    Output, double X[N], the abscissas.
        //
    {
        switch (n)
        {
            case 1:
                x[0] = 0.0000000000000000E+00;
                break;
            case 3:
                x[0] = -1.2247448713915889E+00;
                x[1] = 0.0000000000000000E+00;
                x[2] = 1.2247448713915889E+00;
                break;
            case 9:
                x[0] = -2.9592107790638380E+00;
                x[1] = -2.0232301911005157E+00;
                x[2] = -1.2247448713915889E+00;
                x[3] = -5.2403354748695763E-01;
                x[4] = 0.0000000000000000E+00;
                x[5] = 5.2403354748695763E-01;
                x[6] = 1.2247448713915889E+00;
                x[7] = 2.0232301911005157E+00;
                x[8] = 2.9592107790638380E+00;
                break;
            case 19:
                x[0] = -4.4995993983103881E+00;
                x[1] = -3.6677742159463378E+00;
                x[2] = -2.9592107790638380E+00;
                x[3] = -2.2665132620567876E+00;
                x[4] = -2.0232301911005157E+00;
                x[5] = -1.8357079751751868E+00;
                x[6] = -1.2247448713915889E+00;
                x[7] = -8.7004089535290285E-01;
                x[8] = -5.2403354748695763E-01;
                x[9] = 0.0000000000000000E+00;
                x[10] = 5.2403354748695763E-01;
                x[11] = 8.7004089535290285E-01;
                x[12] = 1.2247448713915889E+00;
                x[13] = 1.8357079751751868E+00;
                x[14] = 2.0232301911005157E+00;
                x[15] = 2.2665132620567876E+00;
                x[16] = 2.9592107790638380E+00;
                x[17] = 3.6677742159463378E+00;
                x[18] = 4.4995993983103881E+00;
                break;
            case 35:
                x[0] = -6.3759392709822356E+00;
                x[1] = -5.6432578578857449E+00;
                x[2] = -5.0360899444730940E+00;
                x[3] = -4.4995993983103881E+00;
                x[4] = -4.0292201405043713E+00;
                x[5] = -3.6677742159463378E+00;
                x[6] = -3.3491639537131945E+00;
                x[7] = -2.9592107790638380E+00;
                x[8] = -2.5705583765842968E+00;
                x[9] = -2.2665132620567876E+00;
                x[10] = -2.0232301911005157E+00;
                x[11] = -1.8357079751751868E+00;
                x[12] = -1.5794121348467671E+00;
                x[13] = -1.2247448713915889E+00;
                x[14] = -8.7004089535290285E-01;
                x[15] = -5.2403354748695763E-01;
                x[16] = -1.7606414208200893E-01;
                x[17] = 0.0000000000000000E+00;
                x[18] = 1.7606414208200893E-01;
                x[19] = 5.2403354748695763E-01;
                x[20] = 8.7004089535290285E-01;
                x[21] = 1.2247448713915889E+00;
                x[22] = 1.5794121348467671E+00;
                x[23] = 1.8357079751751868E+00;
                x[24] = 2.0232301911005157E+00;
                x[25] = 2.2665132620567876E+00;
                x[26] = 2.5705583765842968E+00;
                x[27] = 2.9592107790638380E+00;
                x[28] = 3.3491639537131945E+00;
                x[29] = 3.6677742159463378E+00;
                x[30] = 4.0292201405043713E+00;
                x[31] = 4.4995993983103881E+00;
                x[32] = 5.0360899444730940E+00;
                x[33] = 5.6432578578857449E+00;
                x[34] = 6.3759392709822356E+00;
                break;
            case 37:
                x[0] = -6.853200069757519;
                x[1] = -6.124527854622158;
                x[2] = -5.521865209868350;
                x[3] = -4.986551454150765;
                x[4] = -4.499599398310388;
                x[5] = -4.057956316089741;
                x[6] = -3.667774215946338;
                x[7] = -3.315584617593290;
                x[8] = -2.959210779063838;
                x[9] = -2.597288631188366;
                x[10] = -2.266513262056788;
                x[11] = -2.023230191100516;
                x[12] = -1.835707975175187;
                x[13] = -1.561553427651873;
                x[14] = -1.224744871391589;
                x[15] = -0.870040895352903;
                x[16] = -0.524033547486958;
                x[17] = -0.214618180588171;
                x[18] = 0.000000000000000;
                x[19] = 0.214618180588171;
                x[20] = 0.524033547486958;
                x[21] = 0.870040895352903;
                x[22] = 1.224744871391589;
                x[23] = 1.561553427651873;
                x[24] = 1.835707975175187;
                x[25] = 2.023230191100516;
                x[26] = 2.266513262056788;
                x[27] = 2.597288631188366;
                x[28] = 2.959210779063838;
                x[29] = 3.315584617593290;
                x[30] = 3.667774215946338;
                x[31] = 4.057956316089741;
                x[32] = 4.499599398310388;
                x[33] = 4.986551454150765;
                x[34] = 5.521865209868350;
                x[35] = 6.124527854622158;
                x[36] = 6.853200069757519;
                break;
            case 41:
                x[0] = -7.251792998192644;
                x[1] = -6.547083258397540;
                x[2] = -5.961461043404500;
                x[3] = -5.437443360177798;
                x[4] = -4.953574342912980;
                x[5] = -4.4995993983103881;
                x[6] = -4.070919267883068;
                x[7] = -3.6677742159463378;
                x[8] = -3.296114596212218;
                x[9] = -2.9592107790638380;
                x[10] = -2.630415236459871;
                x[11] = -2.2665132620567876;
                x[12] = -2.043834754429505;
                x[13] = -2.0232301911005157;
                x[14] = -1.8357079751751868;
                x[15] = -1.585873011819188;
                x[16] = -1.2247448713915889;
                x[17] = -0.87004089535290285;
                x[18] = -0.52403354748695763;
                x[19] = -0.195324784415805;
                x[20] = 0.0000000000000000;
                x[21] = 0.195324784415805;
                x[22] = 0.52403354748695763;
                x[23] = 0.87004089535290285;
                x[24] = 1.2247448713915889;
                x[25] = 1.585873011819188;
                x[26] = 1.8357079751751868;
                x[27] = 2.0232301911005157;
                x[28] = 2.043834754429505;
                x[29] = 2.2665132620567876;
                x[30] = 2.630415236459871;
                x[31] = 2.9592107790638380;
                x[32] = 3.296114596212218;
                x[33] = 3.6677742159463378;
                x[34] = 4.070919267883068;
                x[35] = 4.4995993983103881;
                x[36] = 4.953574342912980;
                x[37] = 5.437443360177798;
                x[38] = 5.961461043404500;
                x[39] = 6.547083258397540;
                x[40] = 7.251792998192644;
                break;
            case 43:
                x[0] = -10.167574994881873;
                x[1] = -7.231746029072501;
                x[2] = -6.535398426382995;
                x[3] = -5.954781975039809;
                x[4] = -5.434053000365068;
                x[5] = -4.952329763008589;
                x[6] = -4.4995993983103881;
                x[7] = -4.071335874253583;
                x[8] = -3.6677742159463378;
                x[9] = -3.295265921534226;
                x[10] = -2.9592107790638380;
                x[11] = -2.633356763661946;
                x[12] = -2.2665132620567876;
                x[13] = -2.089340389294661;
                x[14] = -2.0232301911005157;
                x[15] = -1.8357079751751868;
                x[16] = -1.583643465293944;
                x[17] = -1.2247448713915889;
                x[18] = -0.87004089535290285;
                x[19] = -0.52403354748695763;
                x[20] = -0.196029453662011;
                x[21] = 0.0000000000000000;
                x[22] = 0.196029453662011;
                x[23] = 0.52403354748695763;
                x[24] = 0.87004089535290285;
                x[25] = 1.2247448713915889;
                x[26] = 1.583643465293944;
                x[27] = 1.8357079751751868;
                x[28] = 2.0232301911005157;
                x[29] = 2.089340389294661;
                x[30] = 2.2665132620567876;
                x[31] = 2.633356763661946;
                x[32] = 2.9592107790638380;
                x[33] = 3.295265921534226;
                x[34] = 3.6677742159463378;
                x[35] = 4.071335874253583;
                x[36] = 4.4995993983103881;
                x[37] = 4.952329763008589;
                x[38] = 5.434053000365068;
                x[39] = 5.954781975039809;
                x[40] = 6.535398426382995;
                x[41] = 7.231746029072501;
                x[42] = 10.167574994881873;
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("HERMITE_GENZ_KEISTER_LOOKUP_POINTS - Fatal error!");
                Console.WriteLine("  Illegal input value of N.");
                Console.WriteLine("  N must be 1, 3, 9, 19, 35, 37, 41 or 43.");
                break;
        }
    }

    public static double[] hermite_genz_keister_lookup_points_np(int n, int np, double[] p,
            double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_GENZ_KEISTER_LOOKUP_POINTS_NP looks up Genz-Keister Hermite abscissas.
        //
        //  Discussion:
        //
        //    The integral:
        //
        //      integral ( -oo <= x <= +oo ) f(x) exp ( - x * x ) dx
        //
        //    The quadrature rule:
        //
        //      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
        //
        //    A nested family of rules for the Hermite integration problem
        //    was produced by Genz and Keister.  The structure of the nested
        //    family was denoted by 1+2+6+10+?, that is, it comprised rules
        //    of successive orders O = 1, 3, 9, 19, and a final rule of order
        //    35, 37, 41 or 43.
        //
        //    The precisions of these rules are P = 1, 5, 15, 29, 
        //    with the final rule of precision 51, 55, 63 or 67.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Alan Genz, Bradley Keister,
        //    Fully symmetric interpolatory rules for multiple integrals
        //    over infinite regions with Gaussian weight,
        //    Journal of Computational and Applied Mathematics,
        //    Volume 71, 1996, pages 299-309
        //
        //    Florian Heiss, Viktor Winschel,
        //    Likelihood approximation by numerical integration on sparse grids,
        //    Journal of Econometrics,
        //    Volume 144, 2008, pages 62-80.
        //
        //    Thomas Patterson,
        //    The Optimal Addition of Points to Quadrature Formulae,
        //    Mathematics of Computation,
        //    Volume 22, Number 104, October 1968, pages 847-856.
        //
        //  Parameters:
        //
        //    Input, int N, the order.
        //    N must be 1, 3, 9, 19, 35, 37, 41 or 43.
        //
        //    Input, int NP, the number of parameters.
        //
        //    Input, double P[NP], parameters which are not needed by this function.
        //
        //    Output, double X[N], the abscissas.
        //
    {
        hermite_genz_keister_lookup_points(n, ref x);

        return x;
    }

    public static void hermite_genz_keister_lookup_weights(int n, ref double[] w)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_GENZ_KEISTER_LOOKUP_WEIGHTS looks up Genz-Keister Hermite weights.
        //
        //  Discussion:
        //
        //    The integral:
        //
        //      integral ( -oo <= x <= +oo ) f(x) exp ( - x * x ) dx
        //
        //    The quadrature rule:
        //
        //      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
        //
        //    A nested family of rules for the Hermite integration problem
        //    was produced by Genz and Keister.  The structure of the nested
        //    family was denoted by 1+2+6+10+?, that is, it comprised rules
        //    of successive orders O = 1, 3, 9, 19, and a final rule of order
        //    35, 37, 41 or 43.
        //
        //    The precisions of these rules are P = 1, 5, 15, 29, 
        //    with the final rule of precision 51, 55, 63 or 67.
        //
        //    Three related families begin the same way, but end with a different final
        //    rule.  As a convenience, this function includes these final rules as well:
        //
        //    Designation  Orders       Precisions
        //
        //    1+2+6+10+16, 1,3,9,19,35  1,5,15,29,51
        //    1+2+6+10+18  1,3,9,19,37  1,5,15,29,55
        //    1+2+6+10+22  1,3,9,19,41  1,5,15,29,63
        //    1+2+6+10+24  1,3,9,19,43  1,5,15,29,67
        //
        //    Some of the data in this function was kindly supplied directly by
        //    Alan Genz on 24 April 2011.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Alan Genz, Bradley Keister,
        //    Fully symmetric interpolatory rules for multiple integrals
        //    over infinite regions with Gaussian weight,
        //    Journal of Computational and Applied Mathematics,
        //    Volume 71, 1996, pages 299-309
        //
        //    Florian Heiss, Viktor Winschel,
        //    Likelihood approximation by numerical integration on sparse grids,
        //    Journal of Econometrics,
        //    Volume 144, 2008, pages 62-80.
        //
        //    Thomas Patterson,
        //    The Optimal Addition of Points to Quadrature Formulae,
        //    Mathematics of Computation,
        //    Volume 22, Number 104, October 1968, pages 847-856.
        //
        //  Parameters:
        //
        //    Input, int N, the order.
        //    N must be 1, 3, 9, 19, 35, 37, 41, or 43.
        //
        //    Output, double W[N], the weights.
        //
    {
        switch (n)
        {
            case 1:
                w[0] = 1.7724538509055159E+00;
                break;
            case 3:
                w[0] = 2.9540897515091930E-01;
                w[1] = 1.1816359006036772E+00;
                w[2] = 2.9540897515091930E-01;
                break;
            case 9:
                w[0] = 1.6708826306882348E-04;
                w[1] = 1.4173117873979098E-02;
                w[2] = 1.6811892894767771E-01;
                w[3] = 4.7869428549114124E-01;
                w[4] = 4.5014700975378197E-01;
                w[5] = 4.7869428549114124E-01;
                w[6] = 1.6811892894767771E-01;
                w[7] = 1.4173117873979098E-02;
                w[8] = 1.6708826306882348E-04;
                break;
            case 19:
                w[0] = 1.5295717705322357E-09;
                w[1] = 1.0802767206624762E-06;
                w[2] = 1.0656589772852267E-04;
                w[3] = 5.1133174390883855E-03;
                w[4] = -1.1232438489069229E-02;
                w[5] = 3.2055243099445879E-02;
                w[6] = 1.1360729895748269E-01;
                w[7] = 1.0838861955003017E-01;
                w[8] = 3.6924643368920851E-01;
                w[9] = 5.3788160700510168E-01;
                w[10] = 3.6924643368920851E-01;
                w[11] = 1.0838861955003017E-01;
                w[12] = 1.1360729895748269E-01;
                w[13] = 3.2055243099445879E-02;
                w[14] = -1.1232438489069229E-02;
                w[15] = 5.1133174390883855E-03;
                w[16] = 1.0656589772852267E-04;
                w[17] = 1.0802767206624762E-06;
                w[18] = 1.5295717705322357E-09;
                break;
            case 35:
                w[0] = 1.8684014894510604E-18;
                w[1] = 9.6599466278563243E-15;
                w[2] = 5.4896836948499462E-12;
                w[3] = 8.1553721816916897E-10;
                w[4] = 3.7920222392319532E-08;
                w[5] = 4.3737818040926989E-07;
                w[6] = 4.8462799737020461E-06;
                w[7] = 6.3328620805617891E-05;
                w[8] = 4.8785399304443770E-04;
                w[9] = 1.4515580425155904E-03;
                w[10] = 4.0967527720344047E-03;
                w[11] = 5.5928828911469180E-03;
                w[12] = 2.7780508908535097E-02;
                w[13] = 8.0245518147390893E-02;
                w[14] = 1.6371221555735804E-01;
                w[15] = 2.6244871488784277E-01;
                w[16] = 3.3988595585585218E-01;
                w[17] = 9.1262675363737921E-04;
                w[18] = 3.3988595585585218E-01;
                w[19] = 2.6244871488784277E-01;
                w[20] = 1.6371221555735804E-01;
                w[21] = 8.0245518147390893E-02;
                w[22] = 2.7780508908535097E-02;
                w[23] = 5.5928828911469180E-03;
                w[24] = 4.0967527720344047E-03;
                w[25] = 1.4515580425155904E-03;
                w[26] = 4.8785399304443770E-04;
                w[27] = 6.3328620805617891E-05;
                w[28] = 4.8462799737020461E-06;
                w[29] = 4.3737818040926989E-07;
                w[30] = 3.7920222392319532E-08;
                w[31] = 8.1553721816916897E-10;
                w[32] = 5.4896836948499462E-12;
                w[33] = 9.6599466278563243E-15;
                w[34] = 1.8684014894510604E-18;
                break;
            case 37:
                w[0] = 0.337304188079177058E-20;
                w[1] = 0.332834739632930463E-16;
                w[2] = 0.323016866782871498E-13;
                w[3] = 0.809333688669950037E-11;
                w[4] = 0.748907559239519284E-09;
                w[5] = 0.294146671497083432E-07;
                w[6] = 0.524482423744884136E-06;
                w[7] = 0.586639457073896277E-05;
                w[8] = 0.571885531470621903E-04;
                w[9] = 0.41642095727577091E-03;
                w[10] = 0.174733389581099482E-02;
                w[11] = 0.313373786000304381E-02;
                w[12] = 0.768092665770660459E-02;
                w[13] = 0.274962713372148476E-01;
                w[14] = 0.783630990508037449E-01;
                w[15] = 0.16611584261479281E+00;
                w[16] = 0.253636910481387185E+00;
                w[17] = 0.261712932511430884E+00;
                w[18] = 0.171719680968980257E+00;
                w[19] = 0.261712932511430884E+00;
                w[20] = 0.253636910481387185E+00;
                w[21] = 0.16611584261479281E+00;
                w[22] = 0.783630990508037449E-01;
                w[23] = 0.274962713372148476E-01;
                w[24] = 0.768092665770660459E-02;
                w[25] = 0.313373786000304381E-02;
                w[26] = 0.174733389581099482E-02;
                w[27] = 0.41642095727577091E-03;
                w[28] = 0.571885531470621903E-04;
                w[29] = 0.586639457073896277E-05;
                w[30] = 0.524482423744884136E-06;
                w[31] = 0.294146671497083432E-07;
                w[32] = 0.748907559239519284E-09;
                w[33] = 0.809333688669950037E-11;
                w[34] = 0.323016866782871498E-13;
                w[35] = 0.332834739632930463E-16;
                w[36] = 0.337304188079177058E-20;
                break;
            case 41:
                w[0] = 0.117725656974405367E-22;
                w[1] = 0.152506745534300636E-18;
                w[2] = 0.202183949965101288E-15;
                w[3] = 0.724614869051195508E-13;
                w[4] = 0.103121966469463034E-10;
                w[5] = 0.710371395169350952E-09;
                w[6] = 0.264376044449260516E-07;
                w[7] = 0.558982787078644997E-06;
                w[8] = 0.675628907134744976E-05;
                w[9] = 0.512198007019776873E-04;
                w[10] = 0.335013114947200879E-03;
                w[11] = 0.249379691096933139E-02;
                w[12] = -0.25616995850607458E-01;
                w[13] = 0.317007878644325588E-01;
                w[14] = 0.125041498584003435E-02;
                w[15] = 0.293244560924894295E-01;
                w[16] = 0.799536390803302298E-01;
                w[17] = 0.164543666806555251E+00;
                w[18] = 0.258718519718241095E+00;
                w[19] = 0.293588795735908566E+00;
                w[20] = 0.997525375254611951E-01;
                w[21] = 0.293588795735908566E+00;
                w[22] = 0.258718519718241095E+00;
                w[23] = 0.164543666806555251E+00;
                w[24] = 0.799536390803302298E-01;
                w[25] = 0.293244560924894295E-01;
                w[26] = 0.125041498584003435E-02;
                w[27] = 0.317007878644325588E-01;
                w[28] = -0.25616995850607458E-01;
                w[29] = 0.249379691096933139E-02;
                w[30] = 0.335013114947200879E-03;
                w[31] = 0.512198007019776873E-04;
                w[32] = 0.675628907134744976E-05;
                w[33] = 0.558982787078644997E-06;
                w[34] = 0.264376044449260516E-07;
                w[35] = 0.710371395169350952E-09;
                w[36] = 0.103121966469463034E-10;
                w[37] = 0.724614869051195508E-13;
                w[38] = 0.202183949965101288E-15;
                w[39] = 0.152506745534300636E-18;
                w[40] = 0.117725656974405367E-22;
                break;
            case 43:
                w[0] = 0.968100020641528185E-37;
                w[1] = 0.15516931262860431E-22;
                w[2] = 0.175937309107750992E-18;
                w[3] = 0.217337608710893738E-15;
                w[4] = 0.747837010380540069E-13;
                w[5] = 0.104028132097205732E-10;
                w[6] = 0.70903573389336778E-09;
                w[7] = 0.263481722999966618E-07;
                w[8] = 0.560127964848432175E-06;
                w[9] = 0.680410934802210232E-05;
                w[10] = 0.508343873102544037E-04;
                w[11] = 0.32753080006610181E-03;
                w[12] = 0.267479828788552937E-02;
                w[13] = -0.687704270963253854E-02;
                w[14] = 0.119383201790913588E-01;
                w[15] = 0.248083722871002796E-02;
                w[16] = 0.29000335749726387E-01;
                w[17] = 0.798689557875757008E-01;
                w[18] = 0.164609842422580606E+00;
                w[19] = 0.258535954731607738E+00;
                w[20] = 0.292243810406117141E+00;
                w[21] = 0.102730713753441829E+00;
                w[22] = 0.292243810406117141E+00;
                w[23] = 0.258535954731607738E+00;
                w[24] = 0.164609842422580606E+00;
                w[25] = 0.798689557875757008E-01;
                w[26] = 0.29000335749726387E-01;
                w[27] = 0.248083722871002796E-02;
                w[28] = 0.119383201790913588E-01;
                w[29] = -0.687704270963253854E-02;
                w[30] = 0.267479828788552937E-02;
                w[31] = 0.32753080006610181E-03;
                w[32] = 0.508343873102544037E-04;
                w[33] = 0.680410934802210232E-05;
                w[34] = 0.560127964848432175E-06;
                w[35] = 0.263481722999966618E-07;
                w[36] = 0.70903573389336778E-09;
                w[37] = 0.104028132097205732E-10;
                w[38] = 0.747837010380540069E-13;
                w[39] = 0.217337608710893738E-15;
                w[40] = 0.175937309107750992E-18;
                w[41] = 0.15516931262860431E-22;
                w[42] = 0.968100020641528185E-37;
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("HERMITE_GENZ_KEISTER_LOOKUP_WEIGHTS - Fatal error!");
                Console.WriteLine("  Illegal input value of N.");
                Console.WriteLine("  N must be 1, 3, 9, 19, 35, 37, 41 or 43.");
                break;
        }
    }

    public static double[] hermite_genz_keister_lookup_weights_np ( int n, int np, double[] p,
            double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_GENZ_KEISTER_LOOKUP_WEIGHTS_NP looks up Genz-Keister Hermite weights.
        //
        //  Discussion:
        //
        //    The integral:
        //
        //      integral ( -oo <= x <= +oo ) f(x) exp ( - x * x ) dx
        //
        //    The quadrature rule:
        //
        //      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
        //
        //    A nested family of rules for the Hermite integration problem
        //    was produced by Genz and Keister.  The structure of the nested
        //    family was denoted by 1+2+6+10+?, that is, it comprised rules
        //    of successive orders O = 1, 3, 9, 19, and a final rule of order
        //    35, 37, 41 or 43.
        //
        //    The precisions of these rules are P = 1, 5, 15, 29, 
        //    with the final rule of precision 51, 55, 63 or 67.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Alan Genz, Bradley Keister,
        //    Fully symmetric interpolatory rules for multiple integrals
        //    over infinite regions with Gaussian weight,
        //    Journal of Computational and Applied Mathematics,
        //    Volume 71, 1996, pages 299-309
        //
        //    Florian Heiss, Viktor Winschel,
        //    Likelihood approximation by numerical integration on sparse grids,
        //    Journal of Econometrics,
        //    Volume 144, 2008, pages 62-80.
        //
        //    Thomas Patterson,
        //    The Optimal Addition of Points to Quadrature Formulae,
        //    Mathematics of Computation,
        //    Volume 22, Number 104, October 1968, pages 847-856.
        //
        //  Parameters:
        //
        //    Input, int N, the order.
        //    N must be 1, 3, 9, 19, 35, 37, 41 or 43.
        //
        //    Input, int NP, the number of parameters.
        //
        //    Input, double P[NP], parameters which are not needed by this function.
        //
        //    Output, double W[N], the weights.
        //
    {
        hermite_genz_keister_lookup_weights ( n, ref w );

        return w;
    }
        

    public static void hermite_weights(int order, ref double[] weight)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_WEIGHTS returns weights for certain Gauss-Hermite quadrature rules.
        //
        //  Discussion:
        //
        //    The allowed orders are 1, 3, 7, 15, 31, 63 and 127.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 October 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //    Arthur Stroud, Don Secrest,
        //    Gaussian Quadrature Formulas,
        //    Prentice Hall, 1966,
        //    LC: QA299.4G3S7.
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order of the rule.
        //    ORDER must be 1, 3, 7, 15, 31, 63 or 127.
        //
        //    Output, double WEIGHT[ORDER], the weights.
        //    The weights are positive, symmetric and should sum to SQRT(PI).
        //
    {
        switch (order)
        {
            case 1:
                weight[1 - 1] = 1.77245385090551602729816748334E+00;
                break;
            case 3:
                weight[1 - 1] = 0.295408975150919337883027913890E+00;
                weight[2 - 1] = 0.118163590060367735153211165556E+01;
                weight[3 - 1] = 0.295408975150919337883027913890E+00;
                break;
            case 7:
                weight[1 - 1] = 0.971781245099519154149424255939E-03;
                weight[2 - 1] = 0.545155828191270305921785688417E-01;
                weight[3 - 1] = 0.425607252610127800520317466666E+00;
                weight[4 - 1] = 0.810264617556807326764876563813E+00;
                weight[5 - 1] = 0.425607252610127800520317466666E+00;
                weight[6 - 1] = 0.545155828191270305921785688417E-01;
                weight[7 - 1] = 0.971781245099519154149424255939E-03;
                break;
            case 15:
                weight[1 - 1] = 0.152247580425351702016062666965E-08;
                weight[2 - 1] = 0.105911554771106663577520791055E-05;
                weight[3 - 1] = 0.100004441232499868127296736177E-03;
                weight[4 - 1] = 0.277806884291277589607887049229E-02;
                weight[5 - 1] = 0.307800338725460822286814158758E-01;
                weight[6 - 1] = 0.158488915795935746883839384960E+00;
                weight[7 - 1] = 0.412028687498898627025891079568E+00;
                weight[8 - 1] = 0.564100308726417532852625797340E+00;
                weight[9 - 1] = 0.412028687498898627025891079568E+00;
                weight[10 - 1] = 0.158488915795935746883839384960E+00;
                weight[11 - 1] = 0.307800338725460822286814158758E-01;
                weight[12 - 1] = 0.277806884291277589607887049229E-02;
                weight[13 - 1] = 0.100004441232499868127296736177E-03;
                weight[14 - 1] = 0.105911554771106663577520791055E-05;
                weight[15 - 1] = 0.152247580425351702016062666965E-08;
                break;
            case 31:
                weight[1 - 1] = 0.46189683944498305857470556847735E-21;
                weight[2 - 1] = 0.51106090079112519643027197715274E-17;
                weight[3 - 1] = 0.58995564987355133075257722133966E-14;
                weight[4 - 1] = 0.18603735214463569590294465062239E-11;
                weight[5 - 1] = 0.23524920032013205739850619940094E-09;
                weight[6 - 1] = 0.14611988344865057576066495091513E-07;
                weight[7 - 1] = 0.50437125589241034841778074689627E-06;
                weight[8 - 1] = 0.10498602757642934202945441341697E-04;
                weight[9 - 1] = 0.13952090395003623854995664958146E-03;
                weight[10 - 1] = 0.12336833073030489880608311394968E-02;
                weight[11 - 1] = 0.74827999140119116765002499116934E-02;
                weight[12 - 1] = 0.31847230731201222775249585776902E-01;
                weight[13 - 1] = 0.96717948160569462991143316029341E-01;
                weight[14 - 1] = 0.21213278866810461318136114862419E+00;
                weight[15 - 1] = 0.33877265789305344906000174083214E+00;
                weight[16 - 1] = 0.39577855609737786462923720809676E+00;
                weight[17 - 1] = 0.33877265789305344906000174083214E+00;
                weight[18 - 1] = 0.21213278866810461318136114862419E+00;
                weight[19 - 1] = 0.96717948160569462991143316029341E-01;
                weight[20 - 1] = 0.31847230731201222775249585776902E-01;
                weight[21 - 1] = 0.74827999140119116765002499116934E-02;
                weight[22 - 1] = 0.12336833073030489880608311394968E-02;
                weight[23 - 1] = 0.13952090395003623854995664958146E-03;
                weight[24 - 1] = 0.10498602757642934202945441341697E-04;
                weight[25 - 1] = 0.50437125589241034841778074689627E-06;
                weight[26 - 1] = 0.14611988344865057576066495091513E-07;
                weight[27 - 1] = 0.23524920032013205739850619940094E-09;
                weight[28 - 1] = 0.18603735214463569590294465062239E-11;
                weight[29 - 1] = 0.58995564987355133075257722133966E-14;
                weight[30 - 1] = 0.51106090079112519643027197715274E-17;
                weight[31 - 1] = 0.46189683944498305857470556847735E-21;
                break;
            case 63:
                weight[1 - 1] = 0.37099206434787551197827130470031E-47;
                weight[2 - 1] = 0.10400778615192299534481914814892E-41;
                weight[3 - 1] = 0.19796804708258311251124226474396E-37;
                weight[4 - 1] = 0.84687478191640015120141181138947E-34;
                weight[5 - 1] = 0.13071305930779945903630127634063E-30;
                weight[6 - 1] = 0.93437837175367456929765381518998E-28;
                weight[7 - 1] = 0.36027426635173044862245783257252E-25;
                weight[8 - 1] = 0.82963863115951789374753323156164E-23;
                weight[9 - 1] = 0.12266629909105281472971700203949E-20;
                weight[10 - 1] = 0.12288435628797061539461585325494E-18;
                weight[11 - 1] = 0.86925536958188009075932426691516E-17;
                weight[12 - 1] = 0.44857058689176221240330804981619E-15;
                weight[13 - 1] = 0.17335817955735154599902643794700E-13;
                weight[14 - 1] = 0.51265062385038307838565047455223E-12;
                weight[15 - 1] = 0.11808921844532942490513037158404E-10;
                weight[16 - 1] = 0.21508698297808025739828859845140E-09;
                weight[17 - 1] = 0.31371929535285447801497640621672E-08;
                weight[18 - 1] = 0.37041625984781705796752840204084E-07;
                weight[19 - 1] = 0.35734732949879669663960738150956E-06;
                weight[20 - 1] = 0.28393114498380927832990899215541E-05;
                weight[21 - 1] = 0.18709113003730498008961134765721E-04;
                weight[22 - 1] = 0.10284880800653635546698378640623E-03;
                weight[23 - 1] = 0.47411702610173128107201781718693E-03;
                weight[24 - 1] = 0.18409222622384813438539657470055E-02;
                weight[25 - 1] = 0.60436044551187631655712178246467E-02;
                weight[26 - 1] = 0.16829299199599730926458559757600E-01;
                weight[27 - 1] = 0.39858264027692992170237391875317E-01;
                weight[28 - 1] = 0.80467087993950415219587554532823E-01;
                weight[29 - 1] = 0.13871950817615293377792092082674E+00;
                weight[30 - 1] = 0.20448695346833761570957197160475E+00;
                weight[31 - 1] = 0.25799889943058042204920467417642E+00;
                weight[32 - 1] = 0.27876694884838411919175686949858E+00;
                weight[33 - 1] = 0.25799889943058042204920467417642E+00;
                weight[34 - 1] = 0.20448695346833761570957197160475E+00;
                weight[35 - 1] = 0.13871950817615293377792092082674E+00;
                weight[36 - 1] = 0.80467087993950415219587554532823E-01;
                weight[37 - 1] = 0.39858264027692992170237391875317E-01;
                weight[38 - 1] = 0.16829299199599730926458559757600E-01;
                weight[39 - 1] = 0.60436044551187631655712178246467E-02;
                weight[40 - 1] = 0.18409222622384813438539657470055E-02;
                weight[41 - 1] = 0.47411702610173128107201781718693E-03;
                weight[42 - 1] = 0.10284880800653635546698378640623E-03;
                weight[43 - 1] = 0.18709113003730498008961134765721E-04;
                weight[44 - 1] = 0.28393114498380927832990899215541E-05;
                weight[45 - 1] = 0.35734732949879669663960738150956E-06;
                weight[46 - 1] = 0.37041625984781705796752840204084E-07;
                weight[47 - 1] = 0.31371929535285447801497640621672E-08;
                weight[48 - 1] = 0.21508698297808025739828859845140E-09;
                weight[49 - 1] = 0.11808921844532942490513037158404E-10;
                weight[50 - 1] = 0.51265062385038307838565047455223E-12;
                weight[51 - 1] = 0.17335817955735154599902643794700E-13;
                weight[52 - 1] = 0.44857058689176221240330804981619E-15;
                weight[53 - 1] = 0.86925536958188009075932426691516E-17;
                weight[54 - 1] = 0.12288435628797061539461585325494E-18;
                weight[55 - 1] = 0.12266629909105281472971700203949E-20;
                weight[56 - 1] = 0.82963863115951789374753323156164E-23;
                weight[57 - 1] = 0.36027426635173044862245783257252E-25;
                weight[58 - 1] = 0.93437837175367456929765381518998E-28;
                weight[59 - 1] = 0.13071305930779945903630127634063E-30;
                weight[60 - 1] = 0.84687478191640015120141181138947E-34;
                weight[61 - 1] = 0.19796804708258311251124226474396E-37;
                weight[62 - 1] = 0.10400778615192299534481914814892E-41;
                weight[63 - 1] = 0.37099206434787551197827130470031E-47;
                break;
            case 127:
                weight[1 - 1] = 0.12504497577050595552677230002883E-100;
                weight[2 - 1] = 0.17272798059419131415318615789672E-93;
                weight[3 - 1] = 0.89321681571986548608031150791499E-88;
                weight[4 - 1] = 0.77306185240893578449625186483810E-83;
                weight[5 - 1] = 0.20143957652648255497735460506196E-78;
                weight[6 - 1] = 0.21503714733610239701351039429345E-74;
                weight[7 - 1] = 0.11341924208594594813715533569504E-70;
                weight[8 - 1] = 0.33489139011795051950683388483136E-67;
                weight[9 - 1] = 0.60486548964016681064424451668405E-64;
                weight[10 - 1] = 0.71375092946352177824971347343892E-61;
                weight[11 - 1] = 0.57884563374885556636801095624030E-58;
                weight[12 - 1] = 0.33581166223858230300409326551248E-55;
                weight[13 - 1] = 0.14394641949253923568603163698953E-52;
                weight[14 - 1] = 0.46821808383216117724080263903889E-50;
                weight[15 - 1] = 0.11817054440684264071348471955361E-47;
                weight[16 - 1] = 0.23581659156008927203181682045005E-45;
                weight[17 - 1] = 0.37814427940797540210712758405540E-43;
                weight[18 - 1] = 0.49411031115771638145610738414006E-41;
                weight[19 - 1] = 0.53255303775425059266087298458297E-39;
                weight[20 - 1] = 0.47854390680131484999315199332765E-37;
                weight[21 - 1] = 0.36191883445952356128627543209554E-35;
                weight[22 - 1] = 0.23232083386343554805352497446119E-33;
                weight[23 - 1] = 0.12753331411008716683688974281454E-31;
                weight[24 - 1] = 0.60277753850758742112436095241270E-30;
                weight[25 - 1] = 0.24679773241777200207460855084439E-28;
                weight[26 - 1] = 0.88019567691698482573264198727415E-27;
                weight[27 - 1] = 0.27482489212040561315005725890593E-25;
                weight[28 - 1] = 0.75468218903085486125222816438456E-24;
                weight[29 - 1] = 0.18303134636280466270545996891835E-22;
                weight[30 - 1] = 0.39355990860860813085582448449811E-21;
                weight[31 - 1] = 0.75293161638581191068419292570042E-20;
                weight[32 - 1] = 0.12857997786722855037584105682618E-18;
                weight[33 - 1] = 0.19659326888445857792541925311450E-17;
                weight[34 - 1] = 0.26986511907214101894995783364250E-16;
                weight[35 - 1] = 0.33344414303198856330118301113874E-15;
                weight[36 - 1] = 0.37173303125150639885726463109574E-14;
                weight[37 - 1] = 0.37473954472839737091885387788983E-13;
                weight[38 - 1] = 0.34230094493397259538669512076007E-12;
                weight[39 - 1] = 0.28385303724993373166810860630552E-11;
                weight[40 - 1] = 0.21406920290454669208938772802828E-10;
                weight[41 - 1] = 0.14706331273431716244229273183839E-09;
                weight[42 - 1] = 0.92173940967434659264335883218167E-09;
                weight[43 - 1] = 0.52781663936972714041837056042506E-08;
                weight[44 - 1] = 0.27650497044951117835905283127679E-07;
                weight[45 - 1] = 0.13267855842539464770913063113371E-06;
                weight[46 - 1] = 0.58380944276113062188573331195042E-06;
                weight[47 - 1] = 0.23581561724775629112332165335800E-05;
                weight[48 - 1] = 0.87524468034280444703919485644809E-05;
                weight[49 - 1] = 0.29876790535909012274846532159647E-04;
                weight[50 - 1] = 0.93874435720072545206729594267039E-04;
                weight[51 - 1] = 0.27170762627931172053444716883938E-03;
                weight[52 - 1] = 0.72493929742498358979684249380921E-03;
                weight[53 - 1] = 0.17841208326763432884316727108264E-02;
                weight[54 - 1] = 0.40524855186046131499765636276283E-02;
                weight[55 - 1] = 0.85000263041544110385806705526917E-02;
                weight[56 - 1] = 0.16471142241609687824005585301760E-01;
                weight[57 - 1] = 0.29499296248213632269675010319119E-01;
                weight[58 - 1] = 0.48847387114300011006959603975676E-01;
                weight[59 - 1] = 0.74807989768583731416517226905270E-01;
                weight[60 - 1] = 0.10598520508090929403834368934301E+00;
                weight[61 - 1] = 0.13893945309051540832066283010510E+00;
                weight[62 - 1] = 0.16856236074207929740526975049765E+00;
                weight[63 - 1] = 0.18927849580120432177170145550076E+00;
                weight[64 - 1] = 0.19673340688823289786163676995151E+00;
                weight[65 - 1] = 0.18927849580120432177170145550076E+00;
                weight[66 - 1] = 0.16856236074207929740526975049765E+00;
                weight[67 - 1] = 0.13893945309051540832066283010510E+00;
                weight[68 - 1] = 0.10598520508090929403834368934301E+00;
                weight[69 - 1] = 0.74807989768583731416517226905270E-01;
                weight[70 - 1] = 0.48847387114300011006959603975676E-01;
                weight[71 - 1] = 0.29499296248213632269675010319119E-01;
                weight[72 - 1] = 0.16471142241609687824005585301760E-01;
                weight[73 - 1] = 0.85000263041544110385806705526917E-02;
                weight[74 - 1] = 0.40524855186046131499765636276283E-02;
                weight[75 - 1] = 0.17841208326763432884316727108264E-02;
                weight[76 - 1] = 0.72493929742498358979684249380921E-03;
                weight[77 - 1] = 0.27170762627931172053444716883938E-03;
                weight[78 - 1] = 0.93874435720072545206729594267039E-04;
                weight[79 - 1] = 0.29876790535909012274846532159647E-04;
                weight[80 - 1] = 0.87524468034280444703919485644809E-05;
                weight[81 - 1] = 0.23581561724775629112332165335800E-05;
                weight[82 - 1] = 0.58380944276113062188573331195042E-06;
                weight[83 - 1] = 0.13267855842539464770913063113371E-06;
                weight[84 - 1] = 0.27650497044951117835905283127679E-07;
                weight[85 - 1] = 0.52781663936972714041837056042506E-08;
                weight[86 - 1] = 0.92173940967434659264335883218167E-09;
                weight[87 - 1] = 0.14706331273431716244229273183839E-09;
                weight[88 - 1] = 0.21406920290454669208938772802828E-10;
                weight[89 - 1] = 0.28385303724993373166810860630552E-11;
                weight[90 - 1] = 0.34230094493397259538669512076007E-12;
                weight[91 - 1] = 0.37473954472839737091885387788983E-13;
                weight[92 - 1] = 0.37173303125150639885726463109574E-14;
                weight[93 - 1] = 0.33344414303198856330118301113874E-15;
                weight[94 - 1] = 0.26986511907214101894995783364250E-16;
                weight[95 - 1] = 0.19659326888445857792541925311450E-17;
                weight[96 - 1] = 0.12857997786722855037584105682618E-18;
                weight[97 - 1] = 0.75293161638581191068419292570042E-20;
                weight[98 - 1] = 0.39355990860860813085582448449811E-21;
                weight[99 - 1] = 0.18303134636280466270545996891835E-22;
                weight[100 - 1] = 0.75468218903085486125222816438456E-24;
                weight[101 - 1] = 0.27482489212040561315005725890593E-25;
                weight[102 - 1] = 0.88019567691698482573264198727415E-27;
                weight[103 - 1] = 0.24679773241777200207460855084439E-28;
                weight[104 - 1] = 0.60277753850758742112436095241270E-30;
                weight[105 - 1] = 0.12753331411008716683688974281454E-31;
                weight[106 - 1] = 0.23232083386343554805352497446119E-33;
                weight[107 - 1] = 0.36191883445952356128627543209554E-35;
                weight[108 - 1] = 0.47854390680131484999315199332765E-37;
                weight[109 - 1] = 0.53255303775425059266087298458297E-39;
                weight[110 - 1] = 0.49411031115771638145610738414006E-41;
                weight[111 - 1] = 0.37814427940797540210712758405540E-43;
                weight[112 - 1] = 0.23581659156008927203181682045005E-45;
                weight[113 - 1] = 0.11817054440684264071348471955361E-47;
                weight[114 - 1] = 0.46821808383216117724080263903889E-50;
                weight[115 - 1] = 0.14394641949253923568603163698953E-52;
                weight[116 - 1] = 0.33581166223858230300409326551248E-55;
                weight[117 - 1] = 0.57884563374885556636801095624030E-58;
                weight[118 - 1] = 0.71375092946352177824971347343892E-61;
                weight[119 - 1] = 0.60486548964016681064424451668405E-64;
                weight[120 - 1] = 0.33489139011795051950683388483136E-67;
                weight[121 - 1] = 0.11341924208594594813715533569504E-70;
                weight[122 - 1] = 0.21503714733610239701351039429345E-74;
                weight[123 - 1] = 0.20143957652648255497735460506196E-78;
                weight[124 - 1] = 0.77306185240893578449625186483810E-83;
                weight[125 - 1] = 0.89321681571986548608031150791499E-88;
                weight[126 - 1] = 0.17272798059419131415318615789672E-93;
                weight[127 - 1] = 0.12504497577050595552677230002883E-100;
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("HERMITE_WEIGHTS - Fatal error!");
                Console.WriteLine("  Illegal value of ORDER = " + order + "");
                Console.WriteLine("  Legal values are 1, 3, 7, 15, 31, 63 and 127.");
                break;
        }
    }

    public static int[] index_level_hermite(int level, int level_max, int dim_num, int point_num,
            int[] grid_index, int[] grid_base)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX_LEVEL_HERMITE: determine first level at which given index is generated.
        //
        //  Discussion:
        //
        //    We are constructing a sparse grid of Gauss-Hermite points.  The grid
        //    is built up of product grids, with a characteristic LEVEL.  
        //
        //    We are concerned with identifying points in this product grid which
        //    have actually been generated previously, on a lower value of LEVEL.
        //
        //    This routine determines the lowest value of LEVEL at which each of
        //    the input points would be generated.
        //
        //    In 1D, given LEVEL, the number of points is ORDER = 2**(LEVEL+1) + 1,
        //    (except that LEVEL = 0 implies ORDER = 1), the BASE is (ORDER-1)/2, 
        //    and the point INDEX values range from -BASE to +BASE.
        //
        //    The values of INDEX and BASE allow us to determine the abstract
        //    properties of the point.  In particular, if INDEX is 0, the corresponding
        //    Gauss-Legendre abscissa is 0, the special "nested" value we need
        //    to take care of.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 October 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    A Sparse Grid Stochastic Collocation Method for Partial Differential
        //    Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2309-2345.
        //
        //  Parameters:
        //
        //    Input, int LEVEL, the level at which these points were 
        //    generated.  LEVEL_MIN <= LEVEL <= LEVEL_MAX.
        //
        //    Input, int LEVEL_MAX, the maximum level.
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int POINT_NUM, the number of points to be tested.
        //
        //    Input, int GRID_INDEX[DIM_NUM*POINT_NUM], the indices of the 
        //    points to be tested.
        //
        //    Input, int GRID_BASE[DIM_NUM], the "base", which is essentially
        //    the denominator of the index.
        //
        //    Output, int INDEX_LEVEL_HERM[POINT_NUM], the value of LEVEL at 
        //    which the point would first be generated.  This will be the same as
        //    the input value of LEVEL, unless the point has an INDEX of 0 and
        //    a corresponding BASE that is NOT zero.
        //
    {
        int point;

        int[] grid_level = new int[point_num];

        int level_min = Math.Max(0, level_max + 1 - dim_num);
        //
        //  If a point has a DIM-th component whose INDEX is 0, then the 
        //  value of LEVEL at which this point would first be generated is
        //  less than LEVEL, unless the DIM-th component of GRID_BASE is 0.
        //
        for (point = 0; point < point_num; point++)
        {
            grid_level[point] = Math.Max(level, level_min);

            int dim;
            for (dim = 0; dim < dim_num; dim++)
            {
                grid_level[point] = grid_index[dim + point * dim_num] switch
                {
                    0 => Math.Max(grid_level[point] - grid_base[dim], level_min),
                    _ => grid_level[point]
                };
            }
        }

        return grid_level;
    }

    public static double monomial_quadrature_hermite ( int dim_num, int[] expon, int point_num, 
            double[] weight, double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MONOMIAL_QUADRATURE_HERMITE applies Hermite quadrature to a monomial.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 October 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int EXPON[DIM_NUM], the exponents.
        //
        //    Input, int POINT_NUM, the number of points in the rule.
        //
        //    Input, double WEIGHT[POINT_NUM], the quadrature weights.
        //
        //    Input, double X[DIM_NUM*POINT_NUM], the quadrature points.
        //
        //    Output, double MONOMIAL_QUADRATURE, the quadrature error.
        //
    {
        int point;
        //
        //  Get the exact value of the integral of the unscaled monomial.
        //
        double exact = Integral.hermite_integral_nd ( dim_num, expon );
        //
        //  Evaluate the monomial at the quadrature points.
        //
        double[] value = Monomial.monomial_value ( dim_num, point_num, x, expon );
        //
        //  Compute the weighted sum.
        //
        double quad = 0.0;
        for ( point = 0; point < point_num; point++ )
        {
            quad += weight[point] * value[point];
        }

        double quad_error = exact switch
        {
            //
            //  If the exact value is nonzero, use it to scale the data.
            //
            0.0 => Math.Abs(quad),
            _ => Math.Abs((quad - exact) / exact)
        };

        return quad_error;
    }
        
    public static double[] product_weight_hermite ( int dim_num, int[] order_1d, int order_nd )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PRODUCT_WEIGHT_HERMITE: weights for a product Gauss-Hermite rule.
        //
        //  Discussion:
        //
        //    This routine computes the weights for a quadrature rule which is
        //    a product of 1D Gauss-Hermite rules of varying order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 October 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int ORDER_1D[DIM_NUM], the order of the 1D rules.
        //
        //    Input, int ORDER_ND, the order of the product rule.
        //
        //    Output, double PRODUCT_WEIGHT_HERM[ORDER_ND], the product rule weights.
        //
    {
        int dim;
        int order;
        typeMethods.r8vecDPData data = new();

        double[] w_nd = new double[order_nd];
  
        for ( order = 0; order < order_nd; order++ )
        {
            w_nd[order] = 1.0;
        }

        for ( dim = 0; dim < dim_num; dim++ )
        {
            double[] w_1d = new double[order_1d[dim]];
    
            hermite_weights ( order_1d[dim], ref w_1d );

            typeMethods.r8vec_direct_product2 ( ref data, dim, order_1d[dim], w_1d, dim_num, 
                order_nd, ref w_nd );

        }
        return w_nd;
    }
        
    public static void h_quadrature_rule(int nt, ref double[] t, ref double[] wts)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    H_QUADRATURE_RULE: quadrature for H(i,x).
        //
        //  Discussion:
        //
        //    H(i,x) is the physicist's Hermite polynomial of degree I.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NT, the order of the rule.
        //
        //    Output, double T[NT], WTS[NT], the points and weights 
        //    of the rule.
        //
    {
        int i;
            

        for (i = 0; i < nt; i++)
        {
            t[i] = 0.0;
        }

        double[] bj = new double[nt];

        for (i = 0; i < nt; i++)
        {
            bj[i] = Math.Sqrt((i + 1) / 2.0);
        }

        for (i = 0; i < nt; i++)
        {
            wts[i] = 0.0;
        }

        wts[0] = Math.Sqrt(Math.Sqrt(Math.PI));

        IMTQLX.imtqlx(nt, ref t, ref bj, ref wts);

        for (i = 0; i < nt; i++)
        {
            wts[i] *= wts[i];
        }
    }
        
    public static void hc_compute_weights_from_points ( int nhalf, double[] xhalf, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HC_COMPUTE_WEIGHTS_FROM_POINTS: Hermite-Cubic weights, user-supplied points.
        //
        //  Discussion:
        //
        //    An interval [A,B] has been divided by NHALF points X; at each
        //    point both function and derivative information is available.
        //
        //    The piecewise cubic Hermite interpolant is constructed for this data.
        //
        //    A quadrature rule is determined for the interpolant.
        //
        //    There will be N=2*NHALF weights.  If the quadrature rule is to be written 
        //    out, one would normally list each point twice, so that the number of points 
        //    and weights are equal.  The listing of the same point value twice is an
        //    implicit indication that both function and derivative values should be
        //    used.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 March 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NHALF, the number of points, not counting repetitions.
        //
        //    Input, double XHALF[NHALF], the points, without repetition.
        //
        //    Output, double W[2*NHALF], the weights.  The first two weights are 
        //    associated with the first point, and so on.
        //
    {
        int j;

        w[0+0*2] =    0.5 * ( xhalf[1] - xhalf[0] );
        w[1+0*2] = Math.Pow ( xhalf[1] - xhalf[0], 2 ) / 12.0;

        for ( j = 1; j < nhalf - 1; j++ )
        {
            w[0+j*2] = 0.5 * ( xhalf[j+1] - xhalf[j-1] );
            w[1+j*2] =       ( xhalf[j+1] - xhalf[j-1] ) 
                * ( xhalf[j+1] - 2.0 * xhalf[j] + xhalf[j-1] ) / 12.0;
        }

        w[0+(nhalf-1)*2] =      0.5 * ( xhalf[nhalf-1] - xhalf[nhalf-2]   );
        w[1+(nhalf-1)*2] = - Math.Pow ( xhalf[nhalf-2] - xhalf[nhalf-1], 2 ) / 12.0;
    }
        
    public static void hcc_compute ( int n, ref double[] x, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HCC_COMPUTE computes a Hermite-Cubic-Chebyshev-Spacing quadrature rule.
        //
        //  Discussion:
        //
        //    For the HCE rule, we assume that an interval has been divided by
        //    M nodes X into Chebyshev-spaced subintervals, and that at each
        //    abscissa both function and derivative information is available.
        //    The piecewise cubic Hermite interpolant is constructed for this data.
        //    The quadrature rule uses N = 2 * M abscissas, where each node is
        //    listed twice, and the weights occur in pairs, with the first multiplying
        //    the function value and the second the derivative.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 March 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order.
        //    1 <= N.
        //
        //    Output, double X[N], the abscissas.
        //
        //    Output, double W[N], the weights.
        //
    {
        int nhalf = n / 2;
        double[] xhalf = new double[nhalf];

        ClenshawCurtis.clenshaw_curtis_compute_points ( nhalf, ref xhalf );
        typeMethods.r8vec_stutter ( nhalf, xhalf, 2, ref x );
        hc_compute_weights_from_points ( nhalf, xhalf, ref w );
    }
        
    public static void hermite_compute_points ( int order, ref double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_COMPUTE_POINTS computes Hermite quadrature points.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order.
        //
        //    Output, double X[ORDER], the abscissas.
        //
    {
        double[] w = new double[order];

        hermite_compute ( order, ref x, ref w );
    }

    public static double[] hermite_compute_points_np ( int order, int np, double[] p, double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_COMPUTE_POINTS_NP computes Hermite quadrature points.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order.
        //
        //    Input, int NP, the number of parameters.
        //
        //    Input, double P[NP], parameters which are not needed by this function.
        //
        //    Output, double X[ORDER], the abscissas.
        //
    {
        hermite_compute_points ( order, ref x );

        return x;
    }
        
    public static void hcc_compute_weights ( int n, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HCC_COMPUTE_WEIGHTS: Hermite-Cubic-Chebyshev-Spacing quadrature weights.
        //
        //  Discussion:
        //
        //    For the HCE rule, we assume that an interval has been divided by
        //    M nodes X into Chebyshev-spaced subintervals, and that at each
        //    abscissa both function and derivative information is available.
        //    The piecewise cubic Hermite interpolant is constructed for this data.
        //    The quadrature rule uses N = 2 * M abscissas, where each node is
        //    listed twice, and the weights occur in pairs, with the first multiplying
        //    the function value and the second the derivative.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 March 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order.
        //
        //    Output, double W[N], the weights.
        //
    {
        if ( n % 2 != 0 )
        {
            Console.WriteLine("");
            Console.WriteLine("HCC_COMPUTE_WEIGHTS - Fatal error!");
            Console.WriteLine("  Order of rule N is not even.");
            return;
        }

        int nhalf = n / 2;
        double[] xhalf = new double[nhalf];

        hc_compute_weights_from_points ( nhalf, xhalf, ref w );

    }
        
    public static void hce_compute ( int n, ref double[] x, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HCE_COMPUTE computes a Hermite-Cubic-Equal-Spacing quadrature rule.
        //
        //  Discussion:
        //
        //    For the HCE rule, we assume that an interval has been divided by
        //    M nodes X into equally spaced subintervals, and that at each
        //    abscissa both function and derivative information is available.
        //    The piecewise cubic Hermite interpolant is constructed for this data.
        //    The quadrature rule uses N = 2 * M abscissas, where each node is
        //    listed twice, and the weights occur in pairs, with the first multiplying
        //    the function value and the second the derivative.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 March 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order.
        //    1 <= N.
        //
        //    Output, double X[N], the abscissas.
        //
        //    Output, double W[N], the weights.
        //
    {
        const double a_low = 0.0;
        const double a_high = 1.0;

        int nhalf = n / 2;

        double[] xhalf = typeMethods.r8vec_linspace_new ( nhalf, a_low, a_high );
        typeMethods.r8vec_stutter ( nhalf, xhalf, 2, ref x );
        hc_compute_weights_from_points ( nhalf, xhalf, ref w );
    }

    public static void he_quadrature_rule(int nt, ref double[] t, ref double[] wts)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HE_QUADRATURE_RULE: quadrature for He(i,x).
        //
        //  Discussion:
        //
        //    He(i,x) represents the probabilist's Hermite polynomial.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NT, the order of the rule.
        //
        //    Output, double T[NT], WTS[NT], the points and weights 
        //    of the rule.
        //
    {
        int i;
            

        for (i = 0; i < nt; i++)
        {
            t[i] = 0.0;
        }

        double[] bj = new double[nt];

        for (i = 0; i < nt; i++)
        {
            bj[i] = Math.Sqrt((i + 1) / 2.0);
        }

        for (i = 0; i < nt; i++)
        {
            wts[i] = 0.0;
        }

        wts[0] = Math.Sqrt(Math.Sqrt(Math.PI));

        IMTQLX.imtqlx(nt, ref t, ref bj, ref wts);

        for (i = 0; i < nt; i++)
        {
            t[i] *= Math.Sqrt(2.0);
        }

        for (i = 0; i < nt; i++)
        {
            wts[i] = wts[i] * wts[i] * Math.Sqrt(2.0);
        }
    }

    public static void hf_quadrature_rule(int nt, ref double[] t, ref double[] wts)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HF_QUADRATURE_RULE: quadrature for Hf(i,x).
        //
        //  Discussion:
        //
        //    Hf(i,x) represents the Hermite function of "degree" I.   
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NT, the order of the rule.
        //
        //    Output, double T[NT], WTS[NT], the points and weights 
        //    of the rule.
        //
    {
        int i;
            

        for (i = 0; i < nt; i++)
        {
            t[i] = 0.0;
        }

        double[] bj = new double[nt];

        for (i = 0; i < nt; i++)
        {
            bj[i] = Math.Sqrt((i + 1) / 2.0);
        }

        for (i = 0; i < nt; i++)
        {
            wts[i] = 0.0;
        }

        wts[0] = Math.Sqrt(Math.Sqrt(Math.PI));

        IMTQLX.imtqlx(nt, ref t, ref bj, ref wts);

        for (i = 0; i < nt; i++)
        {
            wts[i] = wts[i] * wts[i] * Math.Exp(t[i] * t[i]);
        }
    }
}