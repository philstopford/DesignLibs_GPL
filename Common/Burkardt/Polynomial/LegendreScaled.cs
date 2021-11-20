namespace Burkardt.PolynomialNS;

public static class LegendreScaled
{
    public static double[] klegeypols(double x, double y, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    KLEGEYPOLS evaluates scaled Legendre polynomials.
        //
        //  Discussion:
        //
        //    This routine evaluate a sequence of scaled Legendre polynomials
        //    P_n(x/y) y^n, with the parameter y in [0,1].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU GPL license.
        //
        //  Modified:
        //
        //    30 June 2014
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Hong Xiao, Zydrunas Gimbutas,
        //    A numerical algorithm for the construction of efficient quadrature
        //    rules in two and higher dimensions,
        //    Computers and Mathematics with Applications,
        //    Volume 59, 2010, pages 663-676.
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Input, double Y, the parameter.
        //
        //    Input, int N, the highest degree to be evaluated.
        //
        //    Output, double KLEGEYPOLS[N+1], the polynomial values.
        //
    {
        int k;

        double[] pols = new double[n + 1];

        double pkp1 = 1.0;
        pols[0] = pkp1;
        switch (n)
        {
            case 0:
                return pols;
        }

        double pk = pkp1;
        pkp1 = x;
        pols[1] = pkp1;
        switch (n)
        {
            case 1:
                return pols;
        }

        for (k = 1; k < n; k++)
        {
            double pkm1 = pk;
            pk = pkp1;
            pkp1 = ((2.0 * k + 1.0) * x * pk - k * pkm1 * y * y) / (k + 1.0);
            pols[k + 1] = pkp1;
        }

        return pols;
    }

    public static void klegeypols3(double x, double y, int n, ref double[] pols, ref double[] dersx,
            ref double[] dersy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    KLEGEYPOLS3 evaluates scaled Legendre polynomials and derivatives.
        //
        //  Discussion:
        //
        //    This routine evaluates a sequence of scaled Legendre polynomials
        //    P_n(x/y) y^n, with the parameter y in [0,1], together with their
        //    derivatives with respect to the parameters x and y.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU GPL license.
        //
        //  Modified:
        //
        //    30 June 2014
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Hong Xiao, Zydrunas Gimbutas,
        //    A numerical algorithm for the construction of efficient quadrature
        //    rules in two and higher dimensions,
        //    Computers and Mathematics with Applications,
        //    Volume 59, 2010, pages 663-676.
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Input, double Y, the parameter value.
        //
        //    Input, int N, the highest degree to be evaluated.
        //
        //    Output, double POLS[N+1], the polynomial values.
        //
        //    Output, double DERSX[N+1], the derivatives with respect to X.
        //
        //    Output, double DERSY[N+1], the derivatives with respect to Y.
        //
    {
        int k;

        double pkp1 = 1.0;
        pols[0] = pkp1;
        double dkp1 = 0.0;
        dersx[0] = dkp1;
        double ykp1 = 0.0;
        dersy[0] = ykp1;

        switch (n)
        {
            case 0:
                return;
        }

        double pk = pkp1;
        pkp1 = x;
        pols[1] = pkp1;
        double dk = dkp1;
        dkp1 = 1.0;
        dersx[1] = dkp1;
        double yk = ykp1;
        ykp1 = 0.0;
        dersy[1] = ykp1;

        switch (n)
        {
            case 1:
                return;
        }

        for (k = 1; k <= n - 1; k++)
        {
            double pkm1 = pk;
            pk = pkp1;
            double dkm1 = dk;
            dk = dkp1;
            double ykm1 = yk;
            yk = ykp1;
            pkp1 = ((2.0 * k + 1.0) * x * pk - k * pkm1 * y * y) / (k + 1.0);
            dkp1 = ((2.0 * k + 1.0) * (x * dk + pk)
                    - k * dkm1 * y * y) / (k + 1.0);
            ykp1 = ((2.0 * k + 1.0) * (x * yk)
                    - k * (pkm1 * 2.0 * y + ykm1 * y * y))
                   / (k + 1.0);
            pols[k + 1] = pkp1;
            dersx[k + 1] = dkp1;
            dersy[k + 1] = ykp1;
        }
    }
}