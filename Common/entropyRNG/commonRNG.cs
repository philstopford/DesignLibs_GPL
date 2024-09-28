using System.Collections.Generic;

namespace entropyRNG;

public static class commonRNG
{
    public static readonly List<string> rngTypes = ["System.Random", "Mersenne Twister", "Crypto"];
    public enum rngIndex { system_random, mtwister, crypto/*, sfmt, rei_sfmt, rei_mt, rei_lcg*/ }

}