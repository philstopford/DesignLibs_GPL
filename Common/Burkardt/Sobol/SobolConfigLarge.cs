namespace Burkardt.Sobol;

public static partial class SobolSampler
{
    public partial class SobolConfigLarge
    {
        public static int DIM_MAX = 40;
        public static int DIM_MAX2 = 1111;
        public static int LOG_MAX = 62;
            
        public long seed { get; set; }
        public long atmost { get; set; }
        public bool initialized { get; set; }
        public long dim_num_save { get; set; }
        public long[] lastq { get; set; }
        public long maxcol { get; set; }

        public double recipd { get; set; }
        public long seed_save { get; set; }

        public long[,] v { get; set; }

        public double[] quasi { get; set; }

        public SobolConfigLarge(int maxdim)
        {
            initialized = false;
            dim_num_save = 0;
            lastq = new long[DIM_MAX2];
            seed_save = 1;
            v = new long[DIM_MAX2, LOG_MAX];
            quasi = new double[maxdim];
            maxcol = LOG_MAX;
        }
    }
}