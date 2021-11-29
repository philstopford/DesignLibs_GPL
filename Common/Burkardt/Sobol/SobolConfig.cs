namespace Burkardt.Sobol;

public static partial class SobolSampler
{
    public partial class SobolConfig
    {
        public const int DIM_MAX = 40;
        public const int DIM_MAX2 = 1111;
        public const int LOG_MAX = 30;
            
        public int seed { get; set; }
        public int atmost { get; set; }
        public bool initialized { get; set; }
        public int dim_num_save { get; set; }
        public int[] lastq { get; set; }
        public int maxcol { get; set; }
            
        public float recipd { get; set; }
        public int seed_save { get; set; }

        public int[,] v { get; set; }
            
        public float[] quasi { get; set; }

        public SobolConfig(int maxdim)
        {
            initialized = false;
            dim_num_save = 0;
            lastq = new int[DIM_MAX2];
            seed_save = 1;
            v = new int[DIM_MAX2, LOG_MAX];
            quasi = new float[maxdim];
        }
    }
}