using System.Collections.Generic;
using entropyRNG;

namespace utility;

public static class IListExtensions {
    /// <summary>
    /// Shuffles the element order of the specified list.
    /// </summary>
    public static void Shuffle<T>(this IList<T> ts) {
        int count = ts.Count;
        int last = count - 1;

        for (int i = 0; i < last; ++i) {
            // var r = rnd.Next(i, count);
            int r = RNG.nextint(i, count);
            (ts[i], ts[r]) = (ts[r], ts[i]);
        }
    }
}
