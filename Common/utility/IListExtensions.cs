using System.Collections.Generic;
using entropyRNG;

namespace utility
{
    public static class IListExtensions {
        /// <summary>
        /// Shuffles the element order of the specified list.
        /// </summary>
        public static void Shuffle<T>(this IList<T> ts) {
            var count = ts.Count;
            var last = count - 1;

            for (var i = 0; i < last; ++i) {
                // var r = rnd.Next(i, count);
                var r = RNG.nextint(i, count);
                var tmp = ts[i];
                ts[i] = ts[r];
                ts[r] = tmp;
            }
        }
    }

    public static class testing
    {
        static void testShuffle()
        {
            List<int> t = new List<int>() {1, 2, 3, 4, 5, 6};
            t.Shuffle();
            int[] t2 = {1, 2, 3, 4, 5, 6};
            t2.Shuffle();

        }
        
    }
}