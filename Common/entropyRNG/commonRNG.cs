using System.Collections.Generic;

namespace entropyRNG;

public static class commonRNG
{
    public static List<string> rngTypes = new()
    { "System.Random", "Mersenne Twister", "Crypto"/*, "SIMD Fast Mersenne Twister",
										   "Rei SIMD Fast Mersenne Twister", "Rei Mersenne Twister", "Rei Linear Congruential Generator"*/
    };
    public enum rngIndex { system_random, mtwister, crypto/*, sfmt, rei_sfmt, rei_mt, rei_lcg*/ };

}