namespace Burkardt.Probability;

public static class Multinoulli
{
    public static double multinoulli_pdf ( int x, int n, double[] theta )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MULTINOULLI_PDF evaluates the Multinoulli PDF.
        //
        //  Discussion:
        //
        //    PDF(X) = THETA(X) for 0 <= X < N.
        //           = 0 otherwise
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 September 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int X, the index of the outcome.
        //    0 <= X < N.
        //
        //    Input, int N, the number of legal outcomes.
        //
        //    Input, double THETA[N], the probability of each outcome.
        //
        //    Output, double MULTINOULLI_PDF, the probability of outcome X.
        //
    {
        double value = x switch
        {
            >= 0 when x < n => theta[x],
            _ => 0.0
        };

        return value;
    }
}