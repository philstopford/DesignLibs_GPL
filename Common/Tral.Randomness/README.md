# TRAL.RANDOMNESS
TRAL.RANDOMNESS is a random number class library for C#. It provides common algorithm implementations and
an extended set of random number routines which can be used with any algorithm.

Target Framework: .NET Standard 2.1

There are no dependencies.

Project Page: https://kuiper.zone/tral-randomness-csharp-prng/

## QUICK START
Install the NuGet Release or download and build the source. Also available from NuGet.

### Basic Routines
Here, we create a default generator and demonstrate a selection of routines:

    using Tral.Randomness;

    ...

    // This will initialize generator
    var rand = new RandomGenerator();

    // Integers
    uint u32 = rand.Next32();
    ulong u64 = rand.Next64();
    int i32 = rand.NextInt(5, 10);
    long i64 = rand.NextLong(-5, 10);

    // Doubles
    var d1 = rand.NextDouble();
    var d2 = rand.NextDouble(0, 10);
    var d3 = rand.NextOpenDouble();
    var d4 = rand.NextStdNormal();

    // Bytes
    var b = rand.NextBytes(64);

    // Shuffle
    var shuffled = new int[]{0, 1, 2, 3, 4};
    rand.Shuffle(shuffled);

    // Strings
    var s1 = rand.NextString(3, "TheAlphabet");
    var s2 = rand.NextString(10, IRandomRoutines.AlphaNumericMixed);

### Global Routines
For convenience, we can also use a singleton provided:

    int i32 = RandomGenerator.Global.NextInt(5, 10);

The Global instance is local to the thread.

### Algorithms
Declare the "internal algorithm" to be used as follows:

    using Tral.Randomness;
    using Tral.Randomness.Algorithms;

    ...

    // Xoshiro256++
    var xopp256 = new RandomGenerator<Xoshiro256pp>();

    // Xoshiro256**
    var xoss256 = new RandomGenerator<Xoshiro256ss>();

    // WELL512
    var well512 = new RandomGenerator<Well512>();

    // MT19937-32
    var mt32 = new RandomGenerator<MersenneTwister32>();

    // MT19937-64
    var mt64 = new RandomGenerator<MersenneTwister64>();

    // ISAAC-32
    var isaac32 = new RandomGenerator<Isaac32>();

    // ISAAC-64
    var isaac64 = new RandomGenerator<Isaac64>();

    // SplitMix64
    var split64 = new RandomGenerator<SplitMix64>();

    // RAND48
    var rand48 = new RandomGenerator<Rand48>();

We can pass references to any of the above, without having to specify the algorithm, using the interface
types "IRandomGenerator" or "IRandomRoutines".

### Seeding and Jumping
As follows:

    // Not randomized (false)
    var rand = new RandomGenerator<Xoshiro256pp>(false);

    // Initialize
    rand.Randomize();

    // How many bytes do we need to seed?
    Console.WriteLine(rand.SeedLength);

    // Xoshiro256++ requires 32 bytes
    rand.SetSeed(seedBytes);

    // We can also seed with a simple integer, although how
    // this is done by RandomGenerator is implementation specific.
    rand.SetSeed(98765);

    // And, where supported:
    if (rand.IsJumpable)
    {
        rand.Jump();
    }

It is also possible to **create your own algorithm** by implementing the minimal interface
"IRandomAlgorithm" or, more commonly: "ISeedableAlgorithm" or "IJumpableAlgorithm". This can then be
used with RandomGenerator<TAlgo> which provides the extended routines. See SplitMix64 and Xoshiro256pp
for straightforward implementation examples.

## More Information
You should find TRAL.RANDOMNESS meaningfully documented by code comments.

Project Page: https://kuiper.zone/tral-randomness-csharp-prng/

Andy Thomas / https://kuiper.zone


# NOTICE

Copyright 2020 Andy Thomas

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.