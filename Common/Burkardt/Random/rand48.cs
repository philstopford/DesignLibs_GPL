namespace Burkardt.RandomNS
{
    /// The RAND48 linear congruential algorithm (lrand48). This is a widely used but legacy
    /// generator with 48 bits of internal state. It is not suitable to cryptographic applications.
    /// On construction, the initial state is populated with a fixed arbitrary value.
    public class Rand48
    {
        
        public double drand48() {
            _s = (0x5DEECE66DL * _s + 0xBL) & ((1L << 48) - 1);
            return (double)_s / (1L << 48);
        }
        
        private ulong _s;

        public Rand48()
        {
            SetSeed(0xD0C3F8CBU);
        }

        public string AlgorithmName { get; } = "RAND48";

        public ulong MaxNext { get; } = int.MaxValue;

        public int SeedLength { get; } = 6;

        public ulong Next()
        {
            // Modulus of 2^48
            _s = (25214903917UL * _s + 11UL) & 0x0000FFFFFFFFFFFFUL;
            return _s >> 17;
        }

        public void SetSeedB(byte[] seed)
        {
            ulong s = seed[0];
            s = (s << 8) | seed[1];
            s = (s << 8) | seed[2];
            s = (s << 8) | seed[3];
            s = (s << 8) | seed[4];

            SetSeed((s << 8) | seed[5]);
        }

        public void SetSeed(ulong seed)
        {
            _s = (seed << 16) | 0x0000330EU;
        }

        public Rand48 Clone()
        {
            var clone = new Rand48();
            clone._s = _s;
            return clone;
        }
        
    }
}

