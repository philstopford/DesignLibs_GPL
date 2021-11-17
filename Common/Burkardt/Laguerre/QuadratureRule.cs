using System;
using Burkardt.IntegralNS;
using Burkardt.MatrixNS;
using Burkardt.Types;

namespace Burkardt.Laguerre;

public static partial class QuadratureRule
{
    public static void laguerre_compute_np ( int order, int np, double[] p, ref double[] x,
            ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGUERRE_COMPUTE_NP computes a Laguerre quadrature rule.
        //
        //  Discussion:
        //
        //    The integral:
        //
        //      Integral ( 0 <= X < +oo ) exp ( - X ) * F(X) dX
        //
        //    The quadrature rule:
        //
        //      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
        //
        //    The integral:
        //
        //      Integral ( A <= X < +oo ) F(X) dX
        //
        //    The quadrature rule:
        //
        //      Sum ( 1 <= I <= ORDER ) W(I) * exp ( X(I) ) * F ( X(I) )
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
        laguerre_compute ( order, ref x, ref w );
    }
        
    public static void laguerre_compute ( int n, ref double[] x, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGUERRE_COMPUTE: Laguerre quadrature rule by the Elhay-Kautsky method.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 April 2011
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
        //    Input, int N, the order.
        //
        //    Output, double X[N], the abscissas.
        //
        //    Output, double W[N], the weights.
        //
    {
        double[] bj;
        int i;
        double zemu;
        //
        //  Define the zero-th moment.
        //
        zemu = 1.0;
        //
        //  Define the Jacobi matrix.
        //
        bj = new double[n];

        for ( i = 0; i < n; i++ )
        {
            bj[i] = i + 1;
        }

        for ( i = 0; i < n; i++ )
        {
            x[i] = 2 * i + 1;
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
        
    public static void gen_laguerre_compute_np ( int order, int np, double[] p, ref double[] x,
            ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEN_LAGUERRE_COMPUTE_NP computes a Generalized Laguerre quadrature rule.
        //
        //  Discussion:
        //
        //    In the simplest case, ALPHA is 0, and we are approximating the
        //    integral from 0 to +oo of exp(-X) * F(X).  When this is so,
        //    it is easy to modify the rule to approximate the integral from
        //    A to +oo as well.
        //
        //    If ALPHA is nonzero, then there is no simple way to extend the
        //    rule to approximate the integral from A to +oo.  The simplest
        //    procedures would be to approximate the integral from 0 to A.
        //
        //    If the integral to approximate is:
        //
        //        Integral ( A <= X < +oo ) exp ( - X ) * F(X) dX
        //      or
        //        Integral ( 0 <= X < +oo ) exp ( - X ) * X^ALPHA * F(X) dX
        //
        //    then the quadrature rule is:
        //
        //      exp ( - A ) * Sum ( 1 <= I <= ORDER ) W(I) * F ( A+X(I) )
        //    or
        //      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
        //
        //
        //    If the integral to approximate is:
        //
        //        Integral ( A <= X < +oo ) F(X) dX
        //      or
        //        Integral ( 0 <= X < +oo ) X^ALPHA * F(X) dX
        //
        //    then the quadrature rule is:
        //
        //      exp ( - A ) * Sum ( 1 <= I <= ORDER )
        //        W(I) * exp(A+X(I)) * F ( A+X(I) )
        //    or
        //      Sum ( 1 <= I <= ORDER ) W(I) * exp(X(I)) * F ( X(I) )
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
        //    Input, double P[1], contains parameters.
        //    P[0] = ALPHA, the exponent of the X factor.
        //    Set ALPHA = 0.0 for the simplest rule.
        //    ALPHA must be nonnegative.
        //
        //    Output, double X[ORDER], the abscissas.
        //
        //    Output, double W[ORDER], the weights.
        //
    {
        double alpha;

        alpha = p[0];

        gen_laguerre_compute(order, alpha, ref x, ref w);
    }
        
    public static void gen_laguerre_compute_weights ( int order, double alpha, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEN_LAGUERRE_COMPUTE_WEIGHTS: Generalized Laguerre quadrature weights.
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
        //    Set ALPHA = 0.0 for the simplest rule.
        //    ALPHA must be nonnegative.
        //
        //    Output, double W[ORDER], the weights.
        //
    {
        double[] x;

        x = new double[order];

        gen_laguerre_compute ( order, alpha, ref x, ref w );
    }

    public static double[] gen_laguerre_compute_weights_np ( int order, int np, double[] p,
            double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEN_LAGUERRE_COMPUTE_WEIGHTS_NP: Generalized Laguerre quadrature weights.
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
        //    P[0] = ALPHA, the exponent of the X factor.
        //    Set ALPHA = 0.0 for the simplest rule.
        //    ALPHA must be nonnegative.
        //
        //    Output, double W[ORDER], the weights.
        //
    {
        double alpha;

        alpha = p[0];

        gen_laguerre_compute_weights ( order, alpha, ref w );

        return w;
    }

    public static void gen_laguerre_compute(int n, double alpha, ref double[] x, ref double[] w)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEN_LAGUERRE_COMPUTE: generalized Gauss-Laguerre quadrature rule.
        //
        //  Discussion:
        //
        //    The integral:
        //
        //      integral ( 0 <= x < +oo ) exp ( - x ) * x^alpha * f(x) dx
        //
        //    The quadrature rule:
        //
        //      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
        //
        //    The integral:
        //
        //      integral ( 0 <= x < +oo ) x^alpha * f(x) dx
        //
        //    The quadrature rule:
        //
        //      sum ( 1 <= i <= n ) w(i) * exp ( x(i) ) * f ( x(i) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 April 2011
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
        //    Input, int N, the order.
        //
        //    Input, double ALPHA, the exponent of the X factor.
        //    ALPHA must be nonnegative.
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
        zemu = typeMethods.r8_gamma(alpha + 1.0);
        //
        //  Define the Jacobi matrix.
        //
        bj = new double[n];

        for (i = 0; i < n; i++)
        {
            i_r8 = i + 1;
            bj[i] = Math.Sqrt(i_r8 * (i_r8 + alpha));
        }

        for (i = 0; i < n; i++)
        {
            i_r8 = i + 1;
            x[i] = 2.0 * i_r8 - 1.0 + alpha;
        }

        w[0] = Math.Sqrt(zemu);

        for (i = 1; i < n; i++)
        {
            w[i] = 0.0;
        }

        //
        //  Diagonalize the Jacobi matrix.
        //
        IMTQLX.imtqlx(n, ref x, ref bj, ref w);

        for (i = 0; i < n; i++)
        {
            w[i] *= w[i];
        }
    }

    public static double monomial_quadrature_gen_laguerre ( int expon, double alpha, int order, 
            int option, double[] w, double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MONOMIAL_QUADRATURE_GEN_LAGUERRE applies a quadrature rule to a monomial.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int EXPON, the exponent.
        //
        //    Input, double ALPHA, the exponent of X in the weight factor.
        //
        //    Input, int ORDER, the number of points in the rule.
        //
        //    Input, int OPTION, indicates standard or modified rule.
        //    0, standard generalized Gauss-Laguerre rule for 
        //       integrand x^alpha*exp(-x)*f(x).
        //    1, modified generalized Gauss-Laguerre rule for 
        //       integrand                 f(x).
        //
        //    Input, double W[ORDER], the quadrature weights.
        //
        //    Input, double X[ORDER], the quadrature points.
        //
        //    Output, double MONOMIAL_QUADRATURE_GEN_LAGUERRE, the quadrature error.
        //
    {
        double exact;
        int i;
        double quad;
        double quad_error;
        //
        //  Get the exact value of the integral.
        //
        exact = Integral.gen_laguerre_integral ( expon, alpha );
        //
        //  Evaluate the unweighted monomial at the quadrature points.
        //
        quad = 0.0;

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
                    quad += w[i] * Math.Pow ( x[i], alpha ) * Math.Exp ( - x[i] ) * Math.Pow ( x[i], expon );
                }

                break;
            }
        }
        //
        //  Error:
        //
        quad_error = Math.Abs ( quad - exact ) / exact;

        return quad_error;
    }
        
    public static double[] product_weight_laguerre ( int dim_num, int[] order_1d, int order_nd )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PRODUCT_WEIGHT_LAGUERRE: weights for a product Gauss-Laguerre rule.
        //
        //  Discussion:
        //
        //    This routine computes the weights for a quadrature rule which is
        //    a product of 1D Gauss-Laguerre rules of varying order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 October 2007
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
        typeMethods.r8vecDPData data = new();

        w_nd = new double[order_nd];
  
        for ( order = 0; order < order_nd; order++ )
        {
            w_nd[order] = 1.0;
        }

        for ( dim = 0; dim < dim_num; dim++ )
        {
            w_1d = new double[order_1d[dim]];
    
            laguerre_weights ( order_1d[dim], ref w_1d );

            typeMethods.r8vec_direct_product2 (ref data, dim, order_1d[dim], w_1d, dim_num, 
                order_nd, ref w_nd );

        }
        return w_nd;
    }
    public static void laguerre_abscissa(int dim_num, int point_num, int[] grid_index,
            int[] grid_base, ref double[] grid_point, int gridIndex = 0, int gridBaseIndex = 0, int gridPointIndex = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGUERRE_ABSCISSA sets abscissas for multidimensional Gauss-Laguerre quadrature.
        //
        //  Discussion:
        //
        //    The "nesting" as it occurs for Gauss-Laguerre sparse grids simply
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
        //    10 October 2007
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
        //    Output, double GRID_POINT[DIM_NUM], the grid points of abscissas.
        //
    {
        int dim;
        int level;
        int point;
        int pointer;
        int[] skip = { 0, 1, 4, 11, 26, 57, 120, 247 };
        double[] x =
        {
            1.0E+00,
            0.415774556783479083311533873128E+00,
            0.229428036027904171982205036136E+01,
            0.628994508293747919686641576551E+01,
            0.193043676560362413838247885004E+00,
            0.102666489533919195034519944317E+01,
            0.256787674495074620690778622666E+01,
            0.490035308452648456810171437810E+01,
            0.818215344456286079108182755123E+01,
            0.127341802917978137580126424582E+02,
            0.193957278622625403117125820576E+02,
            0.933078120172818047629030383672E-01,
            0.492691740301883908960101791412E+00,
            0.121559541207094946372992716488E+01,
            0.226994952620374320247421741375E+01,
            0.366762272175143727724905959436E+01,
            0.542533662741355316534358132596E+01,
            0.756591622661306786049739555812E+01,
            0.101202285680191127347927394568E+02,
            0.131302824821757235640991204176E+02,
            0.166544077083299578225202408430E+02,
            0.207764788994487667729157175676E+02,
            0.256238942267287801445868285977E+02,
            0.314075191697539385152432196202E+02,
            0.385306833064860094162515167595E+02,
            0.480260855726857943465734308508E+02,
            0.45901947621108290743496080275224E-01,
            0.24198016382477204890408974151714E+00,
            0.59525389422235073707330165005414E+00,
            1.1066894995329987162111308789792E+00,
            1.7775956928747727211593727482675E+00,
            2.6097034152566806503893375925315E+00,
            3.6051968023400442698805817554243E+00,
            4.7667470844717611313629127271123E+00,
            6.0975545671817409269925429328463E+00,
            7.6014009492331374229360106942867E+00,
            9.2827143134708894182536695297710E+00,
            11.146649755619291358993815629587E+00,
            13.199189576244998522464925028637E+00,
            15.447268315549310075809325891801E+00,
            17.898929826644757646725793817752E+00,
            20.563526336715822170743048968779E+00,
            23.451973482011858591050255575933E+00,
            26.577081352118260459975876986478E+00,
            29.953990872346445506951917840024E+00,
            33.600759532902202735410313885784E+00,
            37.539164407330440882887902558001E+00,
            41.795830870182219981347945853330E+00,
            46.403866806411123136029227604386E+00,
            51.405314476797755161861461088395E+00,
            56.854992868715843620511922055660E+00,
            62.826855908786321453677523304806E+00,
            69.425277191080345623322251656443E+00,
            76.807047763862732837609972285484E+00,
            85.230358607545669169387065607043E+00,
            95.188939891525629981308606853957E+00,
            107.95224382757871475002440117666E+00,
            0.22768893732576153785994330248562E-01,
            0.11998325242727824715771416426383E+00,
            0.29494185444770149577427738517405E+00,
            0.54779087896237725363865073775856E+00,
            0.87869061179931901673895567052285E+00,
            1.2878464335919706302309207788611E+00,
            1.7755123815388553763979463268728E+00,
            2.3419925567085989256055628337716E+00,
            2.9876423223246473939976731053629E+00,
            3.7128695992018000346299637413422E+00,
            4.5181363349503584391105568561550E+00,
            5.4039601781825946286902599782736E+00,
            6.3709163787865330220392250891777E+00,
            7.4196399339311711154888493199004E+00,
            8.5508280008403328312589048722235E+00,
            9.7652425999245366807004592977996E+00,
            11.063713635140661736220550410604E+00,
            12.447142262356492749798687569289E+00,
            13.916504641057818562912967008183E+00,
            15.472856110036296424777143607779E+00,
            17.117335833863588753116900303886E+00,
            18.851171974154856850873483787506E+00,
            20.675687448056515660377265667433E+00,
            22.592306346311528381292277759986E+00,
            24.602561094972638883700642760037E+00,
            26.708100458737343969779087998829E+00,
            28.910698500451382640177718103234E+00,
            31.212264631175912885477773820802E+00,
            33.614854909101154836598842888345E+00,
            36.120684774484823056306328740825E+00,
            38.732143442933582145626041607663E+00,
            41.451810222318741191114726181363E+00,
            44.282473071479233839358857134636E+00,
            47.227149784295686898935095231536E+00,
            50.289112264240695761749021839419E+00,
            53.471914456788652808348280619542E+00,
            56.779424636342062213099781057119E+00,
            60.215862909019862886417550114424E+00,
            63.785845004235974631701139601836E+00,
            67.494433702293885830374325695045E+00,
            71.347199604295266286654803376075E+00,
            75.350293425653234254290504744279E+00,
            79.510532629986309149555391354778E+00,
            83.835506080872257843339817658508E+00,
            88.333701570354369086112766326498E+00,
            93.014662728558547405303399037100E+00,
            97.889184147578140043386727677112E+00,
            102.96955690741381650783952746778E+00,
            108.26988161961595392226350967206E+00,
            113.80647350287462738934485955901E+00,
            119.59839538830458666962452963285E+00,
            125.66817255856119431291196303280E+00,
            132.04277272091165746585590583045E+00,
            138.75498418103789078167590567526E+00,
            145.84541318313540358283994248439E+00,
            153.36548459497863623710815962660E+00,
            161.38215194813761243562172669592E+00,
            169.98570600665839438795175301156E+00,
            179.30366247401580910251827858515E+00,
            189.52789596532475473668721332981E+00,
            200.97521159924656741628671841018E+00,
            214.25368536638788642698056296400E+00,
            230.93465747089703971246562985079E+00,
            0.11339635298518611691893169631306E-01,
            0.59749753435726620281348237057387E-01,
            0.14685098690746167612388223687431E+00,
            0.27267590735859553131378008278900E+00,
            0.43724600644192665554577035869932E+00,
            0.64058688222566929533576416399983E+00,
            0.88272968639058364481487653650042E+00,
            1.1637114160166537661560584700951E+00,
            1.4835750152834613891313584861012E+00,
            1.8423694351613565380686320809853E+00,
            2.2401496839579024244513315656522E+00,
            2.6769768780141303692167869961238E+00,
            3.1529182957082825565771508308846E+00,
            3.6680474360304752540226339926515E+00,
            4.2224440823301888455977876667425E+00,
            4.8161943715870502475665535087286E+00,
            5.4493908694559416755862178908416E+00,
            6.1221326512997254193944584763155E+00,
            6.8345253894122668112237994973336E+00,
            7.5866814466367472174205986836847E+00,
            8.3787199765932725254842120659452E+00,
            9.2107670307426558777922506102445E+00,
            10.082955672528643809166439353647E+00,
            10.995426098858125429803147358780E+00,
            11.948325769197725997610605127857E+00,
            12.941809542585531053723381098192E+00,
            13.976039822878506520014405668679E+00,
            15.051186712579523631574796365435E+00,
            16.167428175612852922977395051768E+00,
            17.324950209443673446561163712616E+00,
            18.523947026965688560811711309349E+00,
            19.764621248611504104071669386884E+00,
            21.047184105173183606877044020054E+00,
            22.371855651855542817648123918101E+00,
            23.738864994122497183652313788712E+00,
            25.148450525937368234077278385644E+00,
            26.600860181041749607253384279755E+00,
            28.096351697964619201753961292129E+00,
            29.635192899504178910610227138642E+00,
            31.217661987479759144214467152615E+00,
            32.844047853610430460522951341338E+00,
            34.514650407441149149105635947422E+00,
            36.229780922306804019615388508885E+00,
            37.989762400399956435968780140278E+00,
            39.794929958089961778396437141707E+00,
            41.645631232730180705153990897484E+00,
            43.542226812286859549950892993822E+00,
            45.485090689228791137996151336673E+00,
            47.474610740231964719468766599146E+00,
            49.511189233379087716728884584381E+00,
            51.595243364671244443182771266934E+00,
            53.727205825819316758288140069145E+00,
            55.907525405447553305830605991732E+00,
            58.136667626022439197077526025660E+00,
            60.415115419018590295707192053805E+00,
            62.743369841051809700207126742685E+00,
            65.121950833949996311956025417139E+00,
            67.551398031997886314411872443149E+00,
            70.032271619884584511229871192030E+00,
            72.565153245206849090888669416801E+00,
            75.150646989739935299354362325096E+00,
            77.789380404085816000647405462136E+00,
            80.482005610750729205803962926758E+00,
            83.229200481195914886796120019048E+00,
            86.031669892953582966798238732643E+00,
            88.890147073512051099652518544282E+00,
            91.805395038358177994971250170499E+00,
            94.778208131331583205387031034825E+00,
            97.809413676305116411054110115424E+00,
            100.89987375017285940371939762172E+00,
            104.05048708821598934704076845022E+00,
            107.26219113414600428423116401414E+00,
            110.53596424851500530602771351277E+00,
            113.87282809075839485348376187652E+00,
            117.27385019192517774095477886379E+00,
            120.74014673718880106173978002719E+00,
            124.27288557955698354259506446928E+00,
            127.87328950885942645093841745425E+00,
            131.54263980314366921809377742137E+00,
            135.28228009311836970132738106369E+00,
            139.09362057432970013964422086977E+00,
            142.97814260643601776808227753574E+00,
            146.93740374437366549441080969072E+00,
            150.97304325252187127492511437460E+00,
            155.08678816034612572229641420609E+00,
            159.28045992663288235401956989889E+00,
            163.55598178957571104015967182053E+00,
            167.91538689194360134245547184721E+00,
            172.36082728473812536838156191681E+00,
            176.89458392960192176311674993508E+00,
            181.51907784036813069227528834025E+00,
            186.23688252828112373861202530357E+00,
            191.05073794450929196790836610789E+00,
            195.96356614879879837839002542988E+00,
            200.97848897600025153696475526130E+00,
            206.09884802468871112127283042753E+00,
            211.32822735671655260572377256981E+00,
            216.67047937658230323477089465777E+00,
            222.12975445929687246267304963754E+00,
            227.71053502072232419089132431317E+00,
            233.41767488282602453367775322563E+00,
            239.25644498830308620018749667089E+00,
            245.23258677871567172531254018984E+00,
            251.35237488718128030005500991754E+00,
            257.62269123792061413076191882313E+00,
            264.05111322908240551754377241831E+00,
            270.64601945722796749299111718606E+00,
            277.41671750163651071798388218104E+00,
            284.37359974220870326674402873120E+00,
            291.52833521346495719581282021650E+00,
            298.89410837028248600878895615414E+00,
            306.48591978262611320418112423947E+00,
            314.32096986471177487400007507615E+00,
            322.41915589128679683349440361344E+00,
            330.80372663802405651933847334878E+00,
            339.50216127832433747735367595958E+00,
            348.54737559472697355480761787441E+00,
            357.97942028029845454049007443090E+00,
            367.84794520076004578858341422871E+00,
            378.21590623135532818332979188889E+00,
            389.16539141251004101579475325153E+00,
            400.80729331451702589996361286427E+00,
            413.29853681779384418008260081859E+00,
            426.87579153663675538288509017051E+00,
            441.93085485310841412460309271842E+00,
            459.21804639888429981971267313224E+00,
            480.69378263388373859704269229304E+00
        };

        for (dim = 0; dim < dim_num; dim++)
        {
            switch (grid_base[gridBaseIndex + dim])
            {
                case < 1:
                    Console.WriteLine("");
                    Console.WriteLine("LAGUERRE_ABSCISSA - Fatal error!");
                    Console.WriteLine("  Some base values are less than 1.");
                    return;
            }
        }

        for (dim = 0; dim < dim_num; dim++)
        {
            switch (grid_base[gridBaseIndex + dim])
            {
                case > 127:
                    Console.WriteLine("");
                    Console.WriteLine("LAGUERRE_ABSCISSA - Fatal error!");
                    Console.WriteLine("  Some base values are greater than 127.");
                    return;
            }
        }

        for (point = 0; point < point_num; point++)
        {
            for (dim = 0; dim < dim_num; dim++)
            {
                level = (int)Math.Log2(grid_base[gridBaseIndex + dim] + 1) - 1;

                pointer = skip[level] + grid_index[gridIndex + dim + point * dim_num];

                switch (pointer)
                {
                    case < 1:
                    case > 247:
                        Console.WriteLine("");
                        Console.WriteLine("LAGUERRE_ABSCISSA - Fatal error!");
                        Console.WriteLine("  POINTER out of bounds.");
                        Console.WriteLine("  POINTER    = " + pointer + "");
                        Console.WriteLine("  POINT      = " + point + "");
                        Console.WriteLine("  DIM        = " + dim + "");
                        Console.WriteLine("  GRID_BASE  = " + grid_base[gridBaseIndex + dim] + "");
                        Console.WriteLine("  LEVEL      = " + level + "");
                        Console.WriteLine("  GRID_INDEX = " + grid_index[gridIndex + dim + point * dim_num] + "");
                        return;
                    default:
                        grid_point[gridPointIndex + dim + point * dim_num] = x[pointer - 1];
                        break;
                }
            }
        }

    }

    public static void laguerre_compute(int order, ref double[] xtab, ref double[] weight,
            double alpha)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGUERRE_COMPUTE computes a Gauss-Laguerre quadrature rule.
        //
        //  Discussion:
        //
        //    In the simplest case, ALPHA is 0, and we are approximating the
        //    integral from 0 to +oo of EXP(-X) * F(X).  When this is so,
        //    it is easy to modify the rule to approximate the integral from
        //    A to +oo as well.
        //
        //    If ALPHA is nonzero, then there is no simple way to extend the
        //    rule to approximate the integral from A to +oo.  The simplest
        //    procedures would be to approximate the integral from 0 to A.
        //
        //    The integration interval is [ A, +oo ) or [ 0, +oo ).
        //
        //    The weight function is w(x) = exp ( -x ) or exp ( -x ) * x^alpha.
        //
        //
        //    If the integral to approximate is:
        //
        //        Integral ( A <= X < +oo ) EXP ( - X ) * F(X) dX
        //      or
        //        Integral ( 0 <= X < +oo ) EXP ( - X ) * X^ALPHA * F(X) dX
        //
        //    then the quadrature rule is:
        //
        //      EXP ( - A ) * Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( A+XTAB(I) )
        //    or
        //      sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
        //
        //    If the integral to approximate is:
        //
        //        Integral ( A <= X < +oo ) F(X) dX
        //      or
        //        Integral ( 0 <= X < +oo ) X^ALPHA * F(X) dX
        //
        //    then the quadrature rule is:
        //
        //      EXP ( - A ) * Sum ( 1 <= I <= ORDER ) 
        //        WEIGHT(I) * EXP(A+XTAB(I)) * F ( A+XTAB(I) )
        //    or
        //      sum ( 1 <= I <= ORDER ) WEIGHT(I) * EXP(XTAB(I)) * F ( XTAB(I) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 May 2006
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
        //    Input, int ORDER, the order of the quadrature rule to be computed.
        //    ORDER must be at least 1.
        //
        //    Output, double XTAB[ORDER], the Gauss-Laguerre abscissas.
        //
        //    Output, double WEIGHT[ORDER], the Gauss-Laguerre weights.
        //
        //    Input, double ALPHA, the exponent of the X factor.
        //    Set ALPHA = 0.0 for the simplest rule.
        //    ALPHA must be nonnegative.
        //
    {
        double[] b;
        double[] c;
        double cc;
        double dp2 = 0;
        int i;
        double p1 = 0;
        double prod;
        double r1;
        double r2;
        double ratio;
        double x = 0;

        b = new double[order];
        c = new double[order];
        //
        //  Set the recursion coefficients.
        //
        for (i = 0; i < order; i++)
        {
            b[i] = alpha + (2 * i + 1);
        }

        for (i = 0; i < order; i++)
        {
            c[i] = i * (alpha + i);
        }

        prod = 1.0;
        for (i = 1; i < order; i++)
        {
            prod *= c[i];
        }

        cc = typeMethods.r8_gamma(alpha + 1.0) * prod;

        for (i = 0; i < order; i++)
        {
            switch (i)
            {
                //
                //  Compute an estimate for the root.
                //
                case 0:
                    x = (1.0 + alpha) * (3.0 + 0.92 * alpha) /
                        (1.0 + 2.4 * order + 1.8 * alpha);
                    break;
                case 1:
                    x += (15.0 + 6.25 * alpha) /
                         (1.0 + 0.9 * alpha + 2.5 * order);
                    break;
                default:
                    r1 = (1.0 + 2.55 * (i - 1))
                         / (1.9 * (i - 1));

                    r2 = 1.26 * (i - 1) * alpha /
                         (1.0 + 3.5 * (i - 1));

                    ratio = (r1 + r2) / (1.0 + 0.3 * alpha);

                    x += ratio * (x - xtab[i - 2]);
                    break;
            }

            //
            //  Use iteration to find the root.
            //
            PolynomialNS.Laguerre.laguerre_root(ref x, order, alpha, ref dp2, ref p1, b, c);
            //
            //  Set the abscissa and weight.
            //
            xtab[i] = x;
            weight[i] = cc / dp2 / p1;
        }
    }
        
    public static void l_quadrature_rule ( int n, ref double[] x, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    L_QUADRATURE_RULE: Gauss-Laguerre quadrature based on L(n,x).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 March 2012
        //
        //  Author:
        //
        //    John Burkardt.
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
        //    Input, int N, the order.
        //
        //    Output, double X[N], the abscissas.
        //
        //    Output, double W[N], the weights.
        //
    {
        double[] bj;
        int i;
        double zemu;
        //
        //  Define the zero-th moment.
        //
        zemu = 1.0;
        //
        //  Define the Jacobi matrix.
        //
        bj = new double[n];
        for ( i = 0; i < n; i++ )
        {
            bj[i] = i + 1;
        }

        for ( i = 0; i < n; i++ )
        {
            x[i] = 2 * i + 1;
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
        
    public static void gen_laguerre_compute_points ( int order, double alpha, ref double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEN_LAGUERRE_COMPUTE_POINTS: Generalized Laguerre quadrature points.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 March 2009
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
        //    Set ALPHA = 0.0 for the simplest rule.
        //    ALPHA must be nonnegative.
        //
        //    Output, double X[ORDER], the abscissas.
        //
    {
        double[] w;

        w = new double[order];

        gen_laguerre_compute ( order, alpha, ref x, ref w );
    }

    public static double[] gen_laguerre_compute_points_np ( int order, int np, double[] p,
            double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEN_LAGUERRE_COMPUTE_POINTS_NP: Generalized Laguerre quadrature points.
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
        //    P[0] = ALPHA, the exponent of the X factor.
        //    Set ALPHA = 0.0 for the simplest rule.
        //    ALPHA must be nonnegative.
        //
        //    Output, double X[ORDER], the abscissas.
        //
    {
        double alpha;

        alpha = p[0];

        gen_laguerre_compute_points ( order, alpha, ref x );

        return x;
    }
        
    public static void laguerre_compute_points ( int order, ref double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGUERRE_COMPUTE_POINTS computes Laguerre quadrature points.
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
        double[] w;

        w = new double[order];

        laguerre_compute ( order, ref x, ref w );
    }

    public static double[] laguerre_compute_points_np ( int order, int np, double[] p, double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGUERRE_COMPUTE_POINTS_NP computes Laguerre quadrature points.
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
        laguerre_compute_points ( order, ref x );

        return x;
    }

    public static void laguerre_compute_weights ( int order, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGUERRE_COMPUTE_WEIGHTS computes Laguerre quadrature weights.
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
        double[] x;

        x = new double[order];

        laguerre_compute ( order, ref x, ref w );
    }

    public static double[] laguerre_compute_weights_np ( int order, int np, double[] p, double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGUERRE_COMPUTE_WEIGHTS_NP computes Laguerre quadrature weights.
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
        laguerre_compute_weights ( order, ref w );

        return w;
    }
        
    public static void laguerre_handle(int order, double a, int option, string output)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGUERRE_HANDLE computes the requested Gauss-Laguerre rule and outputs it.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order of the rule.
        //
        //    Input, double A, the left endpoint of the interval.
        //
        //    Input, int OPTION.
        //    * 0, the integral has the form 
        //      Integral ( A <= x < +oo ) exp(-x) f(x) dx
        //    * 1, the integral has the form 
        //      Integral ( A <= x < +oo )         f(x) dx.
        //
        //    Input, STING OUTPUT_FILE, specifies the output.
        //    * "C++'", print as C++ code.
        //    * "F77", print as FORTRAN77 code.
        //    * "F90", print as FORTRAN90 code.
        //    * "MAT", print as MATLAB code.
        //    * file,  write files "file_w.txt", "file_x.txt", "file_r.txt" 
        //      defining weights, abscissas, and region.
        // 
    {
        int i;
        string output_r;
        string output_w;
        string output_x;
        double[] r;
        double[] w;
        double[] x;

        r = new double[2];
        w = new double[order];
        x = new double[order];

        r[0] = a;
        r[1] = typeMethods.r8_huge();

        laguerre_compute(order, ref x, ref w);
        switch (option)
        {
            //
            //  Modify weights if requested.
            //
            case 1:
            {
                for (i = 0; i < order; i++)
                {
                    w[i] = Math.Exp(x[i]) * w[i];
                }

                break;
            }
        }

        switch (output)
        {
            case "C++":
            {
                Console.WriteLine("//");
                Console.WriteLine("//  Weights W, abscissas X and range R");
                Console.WriteLine("//  for a Gauss-Laguerre quadrature rule");
                Console.WriteLine("//  ORDER = " + order + "");
                Console.WriteLine("//  A = " + a + "");
                Console.WriteLine("//");

                switch (option)
                {
                    case 0:
                        Console.WriteLine("//  OPTION = 0, Standard rule:");
                        Console.WriteLine("//    Integral ( A <= x < +oo ) exp(-x) f(x) dx");
                        Console.WriteLine("//    is to be approximated by");
                        Console.WriteLine("//    sum ( 1 <= I <= ORDER ) w(i) * f(x(i)).");
                        break;
                    default:
                        Console.WriteLine("//  OPTION = 1, modified rule:");
                        Console.WriteLine("//    Integral ( A <= x < +oo ) f(x) dx");
                        Console.WriteLine("//    is to be approximated by");
                        Console.WriteLine("//    sum ( 1 <= I <= ORDER ) w(i) * f(x(i)).");
                        break;
                }

                Console.WriteLine("//");

                for (i = 0; i < order; i++)
                {
                    Console.WriteLine("  w[" + i + "] = "
                                      + w[i].ToString("0.################") + ";");
                }

                Console.WriteLine("");
                for (i = 0; i < order; i++)
                {
                    Console.WriteLine("  x[" + i + "] = "
                                      + x[i].ToString("0.################") + ";");
                }

                Console.WriteLine("");
                for (i = 0; i < 2; i++)
                {
                    Console.WriteLine("  r[" + i + "] = " + r[i] + ";");
                }

                break;
            }
            case "F77":
            {
                Console.WriteLine("c");
                Console.WriteLine("c  Weights W, abscissas X and range R");
                Console.WriteLine("c  for a Gauss-Laguerre quadrature rule");
                Console.WriteLine("c  ORDER = " + order + "");
                Console.WriteLine("c  A = " + a + "");
                Console.WriteLine("c");
                switch (option)
                {
                    case 0:
                        Console.WriteLine("c  OPTION = 0, Standard rule:");
                        Console.WriteLine("c    Integral ( A <= x < +oo ) exp(-x) f(x) dx");
                        Console.WriteLine("c    is to be approximated by");
                        Console.WriteLine("c    sum ( 1 <= I <= ORDER ) w(i) * f(x(i)).");
                        break;
                    default:
                        Console.WriteLine("c  OPTION = 1, modified rule:");
                        Console.WriteLine("c    Integral ( A <= x < +oo ) f(x) dx");
                        Console.WriteLine("c    is to be approximated by");
                        Console.WriteLine("c    sum ( 1 <= I <= ORDER ) w(i) * f(x(i)).");
                        break;
                }

                Console.WriteLine("c");

                for (i = 0; i < order; i++)
                {
                    Console.WriteLine("      w(" + i + 1 + ") = "
                                      + w[i].ToString("0.################") + "");
                }

                Console.WriteLine("");
                for (i = 0; i < order; i++)
                {
                    Console.WriteLine("      x(" + i + 1 + ") = "
                                      + x[i].ToString("0.################") + "");
                }

                Console.WriteLine("");
                for (i = 0; i < 2; i++)
                {
                    Console.WriteLine("      r(" + i + 1 + ") = " + r[i] + "");
                }

                break;
            }
            case "F90":
            {
                Console.WriteLine("!");
                Console.WriteLine("!  Weights W, abscissas X and range R");
                Console.WriteLine("!  for a Gauss-Laguerre quadrature rule");
                Console.WriteLine("!  ORDER = " + order + "");
                Console.WriteLine("!  A = " + a + "");
                Console.WriteLine("!");
                switch (option)
                {
                    case 0:
                        Console.WriteLine("!  OPTION = 0, Standard rule:");
                        Console.WriteLine("!    Integral ( A <= x < +oo ) exp(-x) f(x) dx");
                        Console.WriteLine("!    is to be approximated by");
                        Console.WriteLine("!    sum ( 1 <= I <= ORDER ) w(i) * f(x(i)).");
                        break;
                    default:
                        Console.WriteLine("!  OPTION = 1, modified rule:");
                        Console.WriteLine("!    Integral ( A <= x < +oo ) f(x) dx");
                        Console.WriteLine("!    is to be approximated by");
                        Console.WriteLine("!    sum ( 1 <= I <= ORDER ) w(i) * f(x(i)).");
                        break;
                }

                Console.WriteLine("!");

                for (i = 0; i < order; i++)
                {
                    Console.WriteLine("  w(" + i + 1 + ") = "
                                      + w[i].ToString("0.################") + "");
                }

                Console.WriteLine("");
                for (i = 0; i < order; i++)
                {
                    Console.WriteLine("  x(" + i + 1 + ") = "
                                      + x[i].ToString("0.################") + "");
                }

                Console.WriteLine("");
                for (i = 0; i < 2; i++)
                {
                    Console.WriteLine("  r(" + i + 1 + ") = " + r[i] + "");
                }

                break;
            }
            case "MAT":
            {
                Console.WriteLine("%");
                Console.WriteLine("%  Weights W, abscissas X and range R");
                Console.WriteLine("%  for a Gauss-Laguerre quadrature rule");
                Console.WriteLine("%  ORDER = " + order + "");
                Console.WriteLine("%  A = " + a + "");
                Console.WriteLine("%");
                switch (option)
                {
                    case 0:
                        Console.WriteLine("%  OPTION = 0, Standard rule:");
                        Console.WriteLine("%   Integral ( A <= x < +oo ) exp(-x) f(x) dx");
                        Console.WriteLine("%    is to be approximated by");
                        Console.WriteLine("%    sum ( 1 <= I <= ORDER ) w(i) * f(x(i)).");
                        break;
                    default:
                        Console.WriteLine("%  OPTION = 1, modified rule:");
                        Console.WriteLine("%    Integral ( A <= x < +oo ) f(x) dx");
                        Console.WriteLine("%    is to be approximated by");
                        Console.WriteLine("%    sum ( 1 <= I <= ORDER ) w(i) * f(x(i)).");
                        break;
                }

                Console.WriteLine("%");

                for (i = 0; i < order; i++)
                {
                    Console.WriteLine("  w(" + i + 1 + ") = "
                                      + w[i].ToString("0.################") + ";");
                }

                Console.WriteLine("");
                for (i = 0; i < order; i++)
                {
                    Console.WriteLine("  x(" + i + 1 + ") = "
                                      + x[i].ToString("0.################") + ";");
                }

                Console.WriteLine("");
                for (i = 0; i < 2; i++)
                {
                    Console.WriteLine("  r(" + i + 1 + ") = " + r[i] + ";");
                }

                break;
            }
            default:
            {
                switch (option)
                {
                    case 0:
                        output_w = output + "%s_w.txt";
                        output_x = output + "%s_x.txt";
                        output_r = output + "%s_r.txt";
                        break;
                    default:
                        output_w = output + "%s_modified_w.txt";
                        output_x = output + "%s_modified_x.txt";
                        output_r = output + "%s_modified_r.txt";
                        break;
                }

                Console.WriteLine("");
                Console.WriteLine("  Creating quadrature files.");
                Console.WriteLine("");
                Console.WriteLine("  Root file name is     \"" + output + "\".");
                Console.WriteLine("");
                Console.WriteLine("  Weight file will be   \"" + output_w + "\".");
                Console.WriteLine("  Abscissa file will be \"" + output_x + "\".");
                Console.WriteLine("  Region file will be   \"" + output_r + "\".");

                typeMethods.r8mat_write(output_w, 1, order, w);
                typeMethods.r8mat_write(output_x, 1, order, x);
                typeMethods.r8mat_write(output_r, 1, 2, r);
                break;
            }
        }
    }

    public static void laguerre_weights(int order, ref double[] weight)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGUERRE_WEIGHTS returns weights for certain Gauss-Laguerre quadrature rules.
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
        //    10 October 2007
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
        //    The weights are positive, symmetric and should sum to 1.
        //
    {
        switch (order)
        {
            case 1:
                weight[1 - 1] = 1.0E+00;
                break;
            case 3:
                weight[1 - 1] = 0.711093009929173015449590191143E+00;
                weight[2 - 1] = 0.278517733569240848801444888457E+00;
                weight[3 - 1] = 0.103892565015861357489649204007E-01;
                break;
            case 7:
                weight[1 - 1] = 0.409318951701273902130432880018E+00;
                weight[2 - 1] = 0.421831277861719779929281005417E+00;
                weight[3 - 1] = 0.147126348657505278395374184637E+00;
                weight[4 - 1] = 0.206335144687169398657056149642E-01;
                weight[5 - 1] = 0.107401014328074552213195962843E-02;
                weight[6 - 1] = 0.158654643485642012687326223234E-04;
                weight[7 - 1] = 0.317031547899558056227132215385E-07;
                break;
            case 15:
                weight[1 - 1] = 0.218234885940086889856413236448E+00;
                weight[2 - 1] = 0.342210177922883329638948956807E+00;
                weight[3 - 1] = 0.263027577941680097414812275022E+00;
                weight[4 - 1] = 0.126425818105930535843030549378E+00;
                weight[5 - 1] = 0.402068649210009148415854789871E-01;
                weight[6 - 1] = 0.856387780361183836391575987649E-02;
                weight[7 - 1] = 0.121243614721425207621920522467E-02;
                weight[8 - 1] = 0.111674392344251941992578595518E-03;
                weight[9 - 1] = 0.645992676202290092465319025312E-05;
                weight[10 - 1] = 0.222631690709627263033182809179E-06;
                weight[11 - 1] = 0.422743038497936500735127949331E-08;
                weight[12 - 1] = 0.392189726704108929038460981949E-10;
                weight[13 - 1] = 0.145651526407312640633273963455E-12;
                weight[14 - 1] = 0.148302705111330133546164737187E-15;
                weight[15 - 1] = 0.160059490621113323104997812370E-19;
                break;
            case 31:
                weight[1 - 1] = 0.11252789550372583820847728082801E+00;
                weight[2 - 1] = 0.21552760818089123795222505285045E+00;
                weight[3 - 1] = 0.23830825164569654731905788089234E+00;
                weight[4 - 1] = 0.19538830929790229249915303390711E+00;
                weight[5 - 1] = 0.12698283289306190143635272904602E+00;
                weight[6 - 1] = 0.67186168923899300670929441993508E-01;
                weight[7 - 1] = 0.29303224993879487404888669311974E-01;
                weight[8 - 1] = 0.10597569915295736089529380314433E-01;
                weight[9 - 1] = 0.31851272582386980320974842433019E-02;
                weight[10 - 1] = 0.79549548307940382922092149012477E-03;
                weight[11 - 1] = 0.16480052126636687317862967116412E-03;
                weight[12 - 1] = 0.28229237864310816393860971468993E-04;
                weight[13 - 1] = 0.39802902551008580387116174900106E-05;
                weight[14 - 1] = 0.45931839841801061673729694510289E-06;
                weight[15 - 1] = 0.43075545187731100930131457465897E-07;
                weight[16 - 1] = 0.32551249938271570855175749257884E-08;
                weight[17 - 1] = 0.19620246675410594996247151593142E-09;
                weight[18 - 1] = 0.93190499086617587129534716431331E-11;
                weight[19 - 1] = 0.34377541819411620520312597898311E-12;
                weight[20 - 1] = 0.96795247130446716997405035776206E-14;
                weight[21 - 1] = 0.20368066110115247398010624219291E-15;
                weight[22 - 1] = 0.31212687280713526831765358632585E-17;
                weight[23 - 1] = 0.33729581704161052453395678308350E-19;
                weight[24 - 1] = 0.24672796386616696011038363242541E-21;
                weight[25 - 1] = 0.11582201904525643634834564576593E-23;
                weight[26 - 1] = 0.32472922591425422434798022809020E-26;
                weight[27 - 1] = 0.49143017308057432740820076259666E-29;
                weight[28 - 1] = 0.34500071104808394132223135953806E-32;
                weight[29 - 1] = 0.87663710117162041472932760732881E-36;
                weight[30 - 1] = 0.50363643921161490411297172316582E-40;
                weight[31 - 1] = 0.19909984582531456482439549080330E-45;
                break;
            case 63:
                weight[1 - 1] = 0.57118633213868979811587283390476E-01;
                weight[2 - 1] = 0.12067476090640395283319932036351E+00;
                weight[3 - 1] = 0.15925001096581873723870561096472E+00;
                weight[4 - 1] = 0.16875178327560799234596192963585E+00;
                weight[5 - 1] = 0.15366641977668956696193711310131E+00;
                weight[6 - 1] = 0.12368770614716481641086652261948E+00;
                weight[7 - 1] = 0.89275098854848671545279150057422E-01;
                weight[8 - 1] = 0.58258485446105944957571825725160E-01;
                weight[9 - 1] = 0.34546657545992580874717085812508E-01;
                weight[10 - 1] = 0.18675685985714656798286552591203E-01;
                weight[11 - 1] = 0.92233449044093536528490075241649E-02;
                weight[12 - 1] = 0.41671250684839592762582663470209E-02;
                weight[13 - 1] = 0.17238120299900582715386728541955E-02;
                weight[14 - 1] = 0.65320845029716311169340559359043E-03;
                weight[15 - 1] = 0.22677644670909586952405173207471E-03;
                weight[16 - 1] = 0.72127674154810668410750270234861E-04;
                weight[17 - 1] = 0.21011261180466484598811536851241E-04;
                weight[18 - 1] = 0.56035500893357212749181536071292E-05;
                weight[19 - 1] = 0.13673642785604888017836641282292E-05;
                weight[20 - 1] = 0.30507263930195817240736097189550E-06;
                weight[21 - 1] = 0.62180061839309763559981775409241E-07;
                weight[22 - 1] = 0.11566529551931711260022448996296E-07;
                weight[23 - 1] = 0.19614588267565478081534781863335E-08;
                weight[24 - 1] = 0.30286171195709411244334756404054E-09;
                weight[25 - 1] = 0.42521344539400686769012963452599E-10;
                weight[26 - 1] = 0.54202220578073819334698791381873E-11;
                weight[27 - 1] = 0.62627306838597672554166850420603E-12;
                weight[28 - 1] = 0.65474443156573322992307089591924E-13;
                weight[29 - 1] = 0.61815575808729181846302500000047E-14;
                weight[30 - 1] = 0.52592721363507381404263991342633E-15;
                weight[31 - 1] = 0.40230920092646484015391506025408E-16;
                weight[32 - 1] = 0.27600740511819536505013824207729E-17;
                weight[33 - 1] = 0.16936946756968296053322009855265E-18;
                weight[34 - 1] = 0.92689146872177087314963772462726E-20;
                weight[35 - 1] = 0.45093739060365632939780140603959E-21;
                weight[36 - 1] = 0.19435162876132376573629962695374E-22;
                weight[37 - 1] = 0.73926270895169207037999639194513E-24;
                weight[38 - 1] = 0.24714364154434632615980126000066E-25;
                weight[39 - 1] = 0.72288649446741597655145390616476E-27;
                weight[40 - 1] = 0.18407617292614039362985209905608E-28;
                weight[41 - 1] = 0.40583498566841960105759537058880E-30;
                weight[42 - 1] = 0.77000496416438368114463925286343E-32;
                weight[43 - 1] = 0.12488505764999334328843314866038E-33;
                weight[44 - 1] = 0.17185000226767010697663950619912E-35;
                weight[45 - 1] = 0.19896372636672396938013975755522E-37;
                weight[46 - 1] = 0.19199671378804058267713164416870E-39;
                weight[47 - 1] = 0.15278588285522166920459714708240E-41;
                weight[48 - 1] = 0.99054752688842142955854138884590E-44;
                weight[49 - 1] = 0.51597523673029211884228858692990E-46;
                weight[50 - 1] = 0.21249846664084111245693912887783E-48;
                weight[51 - 1] = 0.67903852766852910591172042494884E-51;
                weight[52 - 1] = 0.16466654148296177467908300517887E-53;
                weight[53 - 1] = 0.29509065402691055027053659375033E-56;
                weight[54 - 1] = 0.37838420647571051984882241014675E-59;
                weight[55 - 1] = 0.33358130068542431878174667995217E-62;
                weight[56 - 1] = 0.19223461022273880981363303073329E-65;
                weight[57 - 1] = 0.67812696961083016872779388922288E-69;
                weight[58 - 1] = 0.13404752802440604607620468935693E-72;
                weight[59 - 1] = 0.13109745101805029757648048223928E-76;
                weight[60 - 1] = 0.52624863881401787388694579143866E-81;
                weight[61 - 1] = 0.63780013856587414257760666006511E-86;
                weight[62 - 1] = 0.12997078942372924566347473916943E-91;
                weight[63 - 1] = 0.10008511496968754063443740168421E-98;
                break;
            case 127:
                weight[1 - 1] = 0.28773246692000124355770010301506E-01;
                weight[2 - 1] = 0.63817468175134649363480949265236E-01;
                weight[3 - 1] = 0.91919669721570571389864194652717E-01;
                weight[4 - 1] = 0.11054167914413766381245463002967E+00;
                weight[5 - 1] = 0.11879771633375850188328329422643E+00;
                weight[6 - 1] = 0.11737818530052695148804451630074E+00;
                weight[7 - 1] = 0.10819305984180551488335145581193E+00;
                weight[8 - 1] = 0.93827075290489628080377261401107E-01;
                weight[9 - 1] = 0.76966450960588843995822485928431E-01;
                weight[10 - 1] = 0.59934903912939714332570730063476E-01;
                weight[11 - 1] = 0.44417742073889001371708316272923E-01;
                weight[12 - 1] = 0.31385080966252320983009372215062E-01;
                weight[13 - 1] = 0.21172316041924506411370709025015E-01;
                weight[14 - 1] = 0.13650145364230541652171185564626E-01;
                weight[15 - 1] = 0.84172852710599172279366657385445E-02;
                weight[16 - 1] = 0.49674990059882760515912858620175E-02;
                weight[17 - 1] = 0.28069903895001884631961957446400E-02;
                weight[18 - 1] = 0.15192951003941952460445341057817E-02;
                weight[19 - 1] = 0.78789028751796084086217287140548E-03;
                weight[20 - 1] = 0.39156751064868450584507324648999E-03;
                weight[21 - 1] = 0.18652434268825860550093566260060E-03;
                weight[22 - 1] = 0.85173160415576621908809828160247E-04;
                weight[23 - 1] = 0.37285639197853037712145321577724E-04;
                weight[24 - 1] = 0.15648416791712993947447805296768E-04;
                weight[25 - 1] = 0.62964340695224829035692735524979E-05;
                weight[26 - 1] = 0.24288929711328724574541379938222E-05;
                weight[27 - 1] = 0.89824607890051007201922871545035E-06;
                weight[28 - 1] = 0.31844174740760353710742966328091E-06;
                weight[29 - 1] = 0.10821272905566839211861807542741E-06;
                weight[30 - 1] = 0.35245076750635536015902779085340E-07;
                weight[31 - 1] = 0.11001224365719347407063839761738E-07;
                weight[32 - 1] = 0.32904079616717932125329343003261E-08;
                weight[33 - 1] = 0.94289145237889976419772700772988E-09;
                weight[34 - 1] = 0.25882578904668318184050195309296E-09;
                weight[35 - 1] = 0.68047437103370762630942259017560E-10;
                weight[36 - 1] = 0.17131398805120837835399564475632E-10;
                weight[37 - 1] = 0.41291744524052865469443922304935E-11;
                weight[38 - 1] = 0.95264189718807273220707664873469E-12;
                weight[39 - 1] = 0.21032604432442425932962942047474E-12;
                weight[40 - 1] = 0.44427151938729352860940434285789E-13;
                weight[41 - 1] = 0.89760500362833703323319846405449E-14;
                weight[42 - 1] = 0.17341511407769287074627948346848E-14;
                weight[43 - 1] = 0.32028099548988356631494379835210E-15;
                weight[44 - 1] = 0.56531388950793682022660742095189E-16;
                weight[45 - 1] = 0.95329672799026591234588044025896E-17;
                weight[46 - 1] = 0.15353453477310142565288509437552E-17;
                weight[47 - 1] = 0.23608962179467365686057842132176E-18;
                weight[48 - 1] = 0.34648742794456611332193876653230E-19;
                weight[49 - 1] = 0.48515241897086461320126957663545E-20;
                weight[50 - 1] = 0.64786228633519813428137373790678E-21;
                weight[51 - 1] = 0.82476020965403242936448553126316E-22;
                weight[52 - 1] = 0.10005361880214719793491658282977E-22;
                weight[53 - 1] = 0.11561395116207304954233181263632E-23;
                weight[54 - 1] = 0.12719342731167922655612134264961E-24;
                weight[55 - 1] = 0.13316584714165372967340004160814E-25;
                weight[56 - 1] = 0.13261218454678944033646108509198E-26;
                weight[57 - 1] = 0.12554995447643949807286074138324E-27;
                weight[58 - 1] = 0.11294412178579462703240913107219E-28;
                weight[59 - 1] = 0.96491020279562119228500608131696E-30;
                weight[60 - 1] = 0.78241846768302099396733076955632E-31;
                weight[61 - 1] = 0.60181503542219626658249939076636E-32;
                weight[62 - 1] = 0.43882482704961741551510518054138E-33;
                weight[63 - 1] = 0.30314137647517256304035802501863E-34;
                weight[64 - 1] = 0.19826016543944539545224676057020E-35;
                weight[65 - 1] = 0.12267623373665926559013654872402E-36;
                weight[66 - 1] = 0.71763931692508888943812834967620E-38;
                weight[67 - 1] = 0.39659378833836963584113716149270E-39;
                weight[68 - 1] = 0.20688970553868040099581951696677E-40;
                weight[69 - 1] = 0.10179587017979517245268418427523E-41;
                weight[70 - 1] = 0.47200827745986374625714293679649E-43;
                weight[71 - 1] = 0.20606828985553374825744353490744E-44;
                weight[72 - 1] = 0.84627575907305987245899032156188E-46;
                weight[73 - 1] = 0.32661123687088798658026998931647E-47;
                weight[74 - 1] = 0.11833939207883162380564134612682E-48;
                weight[75 - 1] = 0.40211209123895013807243250164050E-50;
                weight[76 - 1] = 0.12799824394111125389430292847476E-51;
                weight[77 - 1] = 0.38123877747548846504399051365162E-53;
                weight[78 - 1] = 0.10612057542701156767898551949650E-54;
                weight[79 - 1] = 0.27571446947200403594113572720812E-56;
                weight[80 - 1] = 0.66772544240928492881306904862856E-58;
                weight[81 - 1] = 0.15052438383868234954068178600268E-59;
                weight[82 - 1] = 0.31538986800113758526689068500772E-61;
                weight[83 - 1] = 0.61326614299483180785237418887960E-63;
                weight[84 - 1] = 0.11048510030324810567549119229368E-64;
                weight[85 - 1] = 0.18410563538091348076979665543900E-66;
                weight[86 - 1] = 0.28323926570052832195543883237652E-68;
                weight[87 - 1] = 0.40154409843763655508670978777418E-70;
                weight[88 - 1] = 0.52351530215683708779772201956106E-72;
                weight[89 - 1] = 0.62634476665005100555787696642851E-74;
                weight[90 - 1] = 0.68612210535666530365348093803922E-76;
                weight[91 - 1] = 0.68651298840956019297134099761855E-78;
                weight[92 - 1] = 0.62581388433728084867318704240915E-80;
                weight[93 - 1] = 0.51833271237514904046803469968027E-82;
                weight[94 - 1] = 0.38893621571918443533108973497673E-84;
                weight[95 - 1] = 0.26357711379476932781525533730623E-86;
                weight[96 - 1] = 0.16078851293917979699005509638883E-88;
                weight[97 - 1] = 0.87978042070968939637972577886624E-91;
                weight[98 - 1] = 0.43013405077495109903408697802188E-93;
                weight[99 - 1] = 0.18713435881342838527144321803729E-95;
                weight[100 - 1] = 0.72125744708060471675805761366523E-98;
                weight[101 - 1] = 0.24508746062177874383231742333023E-100;
                weight[102 - 1] = 0.73042094619470875777647865078327E-103;
                weight[103 - 1] = 0.18983290818383463537886818579820E-105;
                weight[104 - 1] = 0.42757400244246684123093264825902E-108;
                weight[105 - 1] = 0.82894681420515755691423485228897E-111;
                weight[106 - 1] = 0.13729432219324400013067050156048E-113;
                weight[107 - 1] = 0.19265464126404973222043166489406E-116;
                weight[108 - 1] = 0.22693344503301354826140809941334E-119;
                weight[109 - 1] = 0.22209290603717355061909071271535E-122;
                weight[110 - 1] = 0.17851087685544512662856555121755E-125;
                weight[111 - 1] = 0.11630931990387164467431190485525E-128;
                weight[112 - 1] = 0.60524443584652392290952805077893E-132;
                weight[113 - 1] = 0.24729569115063528647628375096400E-135;
                weight[114 - 1] = 0.77789065006489410364997205809045E-139;
                weight[115 - 1] = 0.18409738662712607039570678274636E-142;
                weight[116 - 1] = 0.31900921131079114970179071968597E-146;
                weight[117 - 1] = 0.39179487139174199737617666077555E-150;
                weight[118 - 1] = 0.32782158394188697053774429820559E-154;
                weight[119 - 1] = 0.17793590713138888062819640128739E-158;
                weight[120 - 1] = 0.58882353408932623157467835381214E-163;
                weight[121 - 1] = 0.10957236509071169877747203273886E-167;
                weight[122 - 1] = 0.10281621114867000898285076975760E-172;
                weight[123 - 1] = 0.41704725557697758145816510853967E-178;
                weight[124 - 1] = 0.58002877720316101774638319601971E-184;
                weight[125 - 1] = 0.18873507745825517106171619101120E-190;
                weight[126 - 1] = 0.69106601826730911682786705950895E-198;
                weight[127 - 1] = 0.43506813201105855628383313334402E-207;
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("LAGUERRE_WEIGHTS - Fatal error!");
                Console.WriteLine("  Illegal value of ORDER = " + order + "");
                Console.WriteLine("  Legal values are 1, 3, 7, 15, 31, 63 and 127.");
                break;
        }
    }

}