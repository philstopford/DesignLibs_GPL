namespace Burkardt.CDFLib;

public static partial class CDF
{
    public static void cumchi ( double x, double df, ref double cum, ref double ccum )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CUMCHI evaluates the cumulative chi-square distribution.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 February 2021
        //
        //  Author:
        //
        //    Barry Brown, James Lovato, Kathy Russell.
        //
        //  Parameters:
        //
        //    Input, double *X, the upper limit of integration.
        //
        //    Input, double *DF, the degrees of freedom of the
        //    chi-square distribution.
        //
        //    Output, double *CUM, the cumulative chi-square distribution.
        //
        //    Output, double *CCUM, the complement of the cumulative
        //    chi-square distribution.
        //
    {
        double a = df * 0.5;
        double xx = x * 0.5;
        cumgam ( xx, a, ref cum, ref ccum );

    }
}