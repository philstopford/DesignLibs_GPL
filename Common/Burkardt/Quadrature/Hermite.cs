using System;
using Burkardt.ClenshawCurtisNS;
using Burkardt.IntegralNS;
using Burkardt.MatrixNS;
using Burkardt.Types;

namespace Burkardt.Quadrature
{
    using Monomial = Burkardt.MonomialNS.Monomial;
    public static class HermiteQuadrature
    {
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
            double alpha;

            alpha = p[0];

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
            double arg;
            double[] bj;
            int i;
            double zemu;
            //
            //  Define the zero-th moment.
            //
            arg = 0.5;
            zemu = typeMethods.r8_gamma ( arg );
            //
            //  Define the Jacobi matrix.
            //
            bj = new double[n];

            for ( i = 0; i < n; i++ )
            {
                bj[i] = Math.Sqrt ( ( double ) ( i + 1 ) / 2.0 );
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
            //
            //  If N is odd, force the middle X to be exactly zero.
            //
            if ( ( n % 2 ) == 1 )
            {
                x[(n-1)/2] = 0.0;
            }

            for ( i = 0; i < n; i++ )
            {
                w[i] = w[i] * w[i];
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
            double[] bj;
            int i;
            double i_r8;
            double zemu;
            //
            //  Define the zero-th moment.
            //
            zemu = typeMethods.r8_gamma ( ( alpha + 1.0 ) / 2.0 );
            //
            //  Define the Jacobi matrix.
            //
            bj = new double[n];

            for ( i = 0; i < n; i++ )
            {
                i_r8 = ( double ) ( i + 1 );
                if ( ( i % 2 ) == 0 )
                {
                    bj[i] = ( i_r8 + alpha ) / 2.0;
                }
                else
                {
                    bj[i] = i_r8 / 2.0;
                }
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
            MatrixNS.IMTQLX.imtqlx ( n, ref x, ref bj, ref w );

            for ( i = 0; i < n; i++ )
            {
                w[i] = w[i] * w[i];
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
            double exact;
            int i;
            double quad;
            double quad_error;
            //
            //  Get the exact value of the integral of the monomial.
            //
            exact = Integral.gen_hermite_integral ( expon, alpha );
            //
            //  Evaluate the unweighted monomial at the quadrature points.
            //
            quad = 0.0;
            if ( option == 0 )
            {
                for ( i = 0; i < order; i++ )
                {
                    quad = quad + w[i] * Math.Pow ( x[i], expon );
                }
            }
            else
            {
                for ( i = 0; i < order; i++ )
                {
                    quad = quad + w[i] * Math.Pow ( Math.Abs ( x[i] ), alpha )
                                       * Math.Exp ( - x[i] * x[i] ) * Math.Pow ( x[i], expon );
                }
            }
            //
            //  Error:
            //
            if ( exact == 0.0 )
            {
                quad_error = Math.Abs ( quad );
            }
            else
            {
                quad_error = Math.Abs ( ( quad - exact ) / exact );
            }
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
            double arg;
            double[] bj;
            int i;
            double zemu;
            //
            //  Define the zero-th moment.
            //
            arg = 0.5;
            zemu = Helpers.Gamma ( arg );
            //
            //  Define the Jacobi matrix.
            //
            bj = new double[n];

            for ( i = 0; i < n; i++ )
            {
                bj[i] = Math.Sqrt ( ( double ) ( i + 1 ) / 2.0 );
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
                w[i] = w[i] * w[i];
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
            int level;
            int point;
            int pointer;
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
                if (grid_base[gridBaseIndex + dim] < 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("HERMITE_ABSCISSA - Fatal error!");
                    Console.WriteLine("  Some base values are less than 0.");
                    return;
                }
            }

            for (dim = 0; dim < dim_num; dim++)
            {
                if (63 < grid_base[gridBaseIndex + dim])
                {
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
                    level = (int)Math.Log2(grid_base[gridBaseIndex + dim] + 1);

                    pointer = skip[level] + (grid_index[gridIndex + (dim + point * dim_num)] + grid_base[gridBaseIndex + dim]);

                    grid_point[gridPointIndex + (dim + point * dim_num)] = x[pointer];
                }
            }

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
            if (order == 1)
            {
                weight[1 - 1] = 1.77245385090551602729816748334E+00;
            }
            else if (order == 3)
            {
                weight[1 - 1] = 0.295408975150919337883027913890E+00;
                weight[2 - 1] = 0.118163590060367735153211165556E+01;
                weight[3 - 1] = 0.295408975150919337883027913890E+00;
            }
            else if (order == 7)
            {
                weight[1 - 1] = 0.971781245099519154149424255939E-03;
                weight[2 - 1] = 0.545155828191270305921785688417E-01;
                weight[3 - 1] = 0.425607252610127800520317466666E+00;
                weight[4 - 1] = 0.810264617556807326764876563813E+00;
                weight[5 - 1] = 0.425607252610127800520317466666E+00;
                weight[6 - 1] = 0.545155828191270305921785688417E-01;
                weight[7 - 1] = 0.971781245099519154149424255939E-03;
            }
            else if (order == 15)
            {
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
            }
            else if (order == 31)
            {
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
            }
            else if (order == 63)
            {
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
            }
            else if (order == 127)
            {
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
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("HERMITE_WEIGHTS - Fatal error!");
                Console.WriteLine("  Illegal value of ORDER = " + order + "");
                Console.WriteLine("  Legal values are 1, 3, 7, 15, 31, 63 and 127.");
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
            int dim;
            int[] grid_level;
            int level_min;
            int point;

            grid_level = new int[point_num];

            level_min = Math.Max(0, level_max + 1 - dim_num);
            //
            //  If a point has a DIM-th component whose INDEX is 0, then the 
            //  value of LEVEL at which this point would first be generated is
            //  less than LEVEL, unless the DIM-th component of GRID_BASE is 0.
            //
            for (point = 0; point < point_num; point++)
            {
                grid_level[point] = Math.Max(level, level_min);

                for (dim = 0; dim < dim_num; dim++)
                {
                    if (grid_index[dim + point * dim_num] == 0)
                    {
                        grid_level[point] = Math.Max(grid_level[point] - grid_base[dim], level_min);
                    }
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
            double exact;
            int point;
            double quad;
            double quad_error;
            double[] value;
            //
            //  Get the exact value of the integral of the unscaled monomial.
            //
            exact = Integral.hermite_integral_nd ( dim_num, expon );
            //
            //  Evaluate the monomial at the quadrature points.
            //
            value = Monomial.monomial_value ( dim_num, point_num, x, expon );
            //
            //  Compute the weighted sum.
            //
            quad = 0.0;
            for ( point = 0; point < point_num; point++ )
            {
                quad = quad + weight[point] * value[point];
            }
            //
            //  If the exact value is nonzero, use it to scale the data.
            //
            if ( exact == 0.0 )
            {
                quad_error = Math.Abs ( quad );
            }
            else
            {
                quad_error = Math.Abs ( ( quad - exact ) / exact );
            }

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
            double[] w_1d;
            double[] w_nd;
            typeMethods.r8vecDPData data = new typeMethods.r8vecDPData();

            w_nd = new double[order_nd];
  
            for ( order = 0; order < order_nd; order++ )
            {
                w_nd[order] = 1.0;
            }

            for ( dim = 0; dim < dim_num; dim++ )
            {
                w_1d = new double[order_1d[dim]];
    
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
            double[] bj;
            int i;
            const double r8_pi = 3.141592653589793;

            for (i = 0; i < nt; i++)
            {
                t[i] = 0.0;
            }

            bj = new double[nt];

            for (i = 0; i < nt; i++)
            {
                bj[i] = Math.Sqrt((double)(i + 1) / 2.0);
            }

            for (i = 0; i < nt; i++)
            {
                wts[i] = 0.0;
            }

            wts[0] = Math.Sqrt(Math.Sqrt(r8_pi));

            IMTQLX.imtqlx(nt, ref t, ref bj, ref wts);

            for (i = 0; i < nt; i++)
            {
                wts[i] = wts[i] * wts[i];
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
            int nhalf;
            double[] xhalf;

            nhalf = n / 2;
            xhalf = new double[nhalf];

            ClenshawCurtis.clenshaw_curtis_compute_points ( nhalf, ref xhalf );
            typeMethods.r8vec_stutter ( nhalf, xhalf, 2, ref x );
            hc_compute_weights_from_points ( nhalf, xhalf, ref w );
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
            int nhalf;
            double[] xhalf;

            if ( ( n % 2 ) != 0 )
            {
                Console.WriteLine("");
                Console.WriteLine("HCC_COMPUTE_WEIGHTS - Fatal error!");
                Console.WriteLine("  Order of rule N is not even.");
                return;
            }

            nhalf = n / 2;
            xhalf = new double[nhalf];

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
            double a_high = 1.0;
            double a_low = 0.0;
            int nhalf;
            double[] xhalf;

            a_low = 0.0;
            a_high = 1.0;

            nhalf = n / 2;

            xhalf = typeMethods.r8vec_linspace_new ( nhalf, a_low, a_high );
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
            double[] bj;
            int i;
            const double r8_pi = 3.141592653589793;

            for (i = 0; i < nt; i++)
            {
                t[i] = 0.0;
            }

            bj = new double[nt];

            for (i = 0; i < nt; i++)
            {
                bj[i] = Math.Sqrt((double)(i + 1) / 2.0);
            }

            for (i = 0; i < nt; i++)
            {
                wts[i] = 0.0;
            }

            wts[0] = Math.Sqrt(Math.Sqrt(r8_pi));

            IMTQLX.imtqlx(nt, ref t, ref bj, ref wts);

            for (i = 0; i < nt; i++)
            {
                t[i] = t[i] * Math.Sqrt(2.0);
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
            double[] bj;
            int i;
            const double r8_pi = 3.141592653589793;

            for (i = 0; i < nt; i++)
            {
                t[i] = 0.0;
            }

            bj = new double[nt];

            for (i = 0; i < nt; i++)
            {
                bj[i] = Math.Sqrt((double)(i + 1) / 2.0);
            }

            for (i = 0; i < nt; i++)
            {
                wts[i] = 0.0;
            }

            wts[0] = Math.Sqrt(Math.Sqrt(r8_pi));

            IMTQLX.imtqlx(nt, ref t, ref bj, ref wts);

            for (i = 0; i < nt; i++)
            {
                wts[i] = wts[i] * wts[i] * Math.Exp(t[i] * t[i]);
            }
        }
    }
}