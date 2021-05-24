using System;
using entropyRNG;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static int get_seed()
        {
            return RNG.nextint(1, Int32.MaxValue);
        }

    }
}