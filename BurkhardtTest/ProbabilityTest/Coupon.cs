namespace Burkardt.ProbabilityTest
{
    partial class Program
    {
static void coupon_complete_pdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    COUPON_COMPLETE_PDF_TEST tests COUPON_COMPLETE_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2016
//
//  Author:
//
//    John Burkardt
//
{
int box_num;
double cdf;
double pdf;
int type_num;

Console.WriteLine("");
Console.WriteLine("COUPON_COMPLETE_PDF_TEST");
Console.WriteLine("  COUPON_COMPLETE_PDF evaluates the coupon collector's");
Console.WriteLine("  complete collection pdf.");
Console.WriteLine("");

for (type_num = 2; type_num <= 4; type_num++)
{
Console.WriteLine("");
Console.WriteLine("  Number of coupon types is " + type_num + "");
Console.WriteLine("");
Console.WriteLine("   BOX_NUM      PDF             CDF");
Console.WriteLine("");
cdf = 0.0;
for (box_num = 1; box_num <= 20; box_num++)
{
pdf = coupon_complete_pdf(type_num, box_num);
cdf = cdf + pdf;
Console.WriteLine("  " + setw(8) + box_num
  + "  " + setw(14) + pdf
  + "  " + setw(14) + cdf + "");
}
}

return;
}

static void coupon_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    COUPON_SAMPLE_TEST tests COUPON_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2016
//
//  Author:
//
//    John Burkardt
//
{
# define N_TRIAL 10
# define MAX_TYPE 25

double average;
int coupon[MAX_TYPE];
double expect;
int i;
int n_coupon;
int n_type;
int seed = 123456789;

Console.WriteLine("");
Console.WriteLine("COUPON_SAMPLE_TEST:");
Console.WriteLine("  COUPON_SAMPLE samples the coupon PDF.");
Console.WriteLine("");

for (n_type = 5; n_type <= MAX_TYPE; n_type = n_type + 5)
{
Console.WriteLine("");
Console.WriteLine("  Number of coupon types is " + n_type + "");
expect = (double) (n_type) * log((double) (n_type));
Console.WriteLine("  Expected wait is about " + expect + "");
Console.WriteLine("");

average = 0.0;
for (i = 1; i <= N_TRIAL; i++)
{
coupon_sample(n_type, seed, coupon, &n_coupon);

Console.WriteLine("  "
  + setw(6) + i + "  "
  + setw(6) + n_coupon + "");

average = average + (double) (n_coupon);
}

average = average / (double) (N_TRIAL);
Console.WriteLine("");
Console.WriteLine("  Average wait was " + average + "");
}

return;
# undef N_TRIAL
# undef MAX_TRIAL
}
        
    }
}