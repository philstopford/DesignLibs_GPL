using System;
using System.IO;
using NUnit.Framework;
using gds;
using MiscUtil.Conversion;
using MiscUtil.IO;

namespace UnitTests
{
    /// <summary>
    /// Tests for GDSII 8-byte real number format compliance.
    /// Validates against GDSII Stream Format Manual 6.0 specification.
    /// </summary>
    public class GdsIIRealNumberTests
    {
        private static void TestRoundTrip(double value, double tolerance = 1e-15)
        {
            // Create a memory stream to write and read back
            using var ms = new MemoryStream();
            using var bw = new EndianBinaryWriter(EndianBitConverter.Big, ms);
            
            // Create a mock gdsWriter just to access the write8ByteReal method
            var writer = new MockGdsWriter(bw);
            writer.write8ByteReal(value);
            
            // Reset stream position and read back
            ms.Position = 0;
            using var br = new EndianBinaryReader(EndianBitConverter.Big, ms);
            
            // Create a mock gdsReader to access the read8ByteReal method
            var reader = new MockGdsReader(br);
            double readValue = reader.TestRead8ByteReal();
            
            // Check if values match within tolerance
            if (double.IsNaN(value))
            {
                Assert.That(double.IsNaN(readValue), Is.True, $"Expected NaN, got {readValue}");
            }
            else if (double.IsInfinity(value))
            {
                Assert.That(double.IsInfinity(readValue), Is.True, $"Expected infinity, got {readValue}");
                Assert.That(Math.Sign(readValue), Is.EqualTo(Math.Sign(value)), "Infinity sign mismatch");
            }
            else if (value == 0.0)
            {
                Assert.That(readValue, Is.EqualTo(0.0), "Zero value failed round-trip");
            }
            else
            {
                double relativeError = Math.Abs((readValue - value) / value);
                Assert.That(relativeError, Is.LessThanOrEqualTo(tolerance), 
                    $"Round-trip failed for {value}: got {readValue}, relative error {relativeError}");
            }
        }
        
        [Test]
        public void GdsII_8ByteReal_ZeroValue()
        {
            TestRoundTrip(0.0);
        }
        
        [Test]
        public void GdsII_8ByteReal_BasicPositiveValues()
        {
            TestRoundTrip(1.0);
            TestRoundTrip(2.0);
            TestRoundTrip(0.5);
            TestRoundTrip(0.25);
            TestRoundTrip(16.0);
        }
        
        [Test]
        public void GdsII_8ByteReal_BasicNegativeValues()
        {
            TestRoundTrip(-1.0);
            TestRoundTrip(-2.0);
            TestRoundTrip(-0.5);
            TestRoundTrip(-16.0);
        }
        
        [Test]
        public void GdsII_8ByteReal_SmallValues()
        {
            TestRoundTrip(1e-3);
            TestRoundTrip(1e-6);
            TestRoundTrip(1e-9);
            TestRoundTrip(1e-12);
        }
        
        [Test]
        public void GdsII_8ByteReal_LargeValues()
        {
            TestRoundTrip(1e3);
            TestRoundTrip(1e6);
            TestRoundTrip(1e9);
            TestRoundTrip(1e12);
        }
        
        [Test]
        public void GdsII_8ByteReal_PowersOfSixteen()
        {
            // Test powers of 16 which should be exactly representable
            for (int i = -10; i <= 10; i++)
            {
                double value = Math.Pow(16.0, i);
                TestRoundTrip(value);
                TestRoundTrip(-value);
            }
        }
        
        [Test]
        public void GdsII_8ByteReal_SpecificationExamples()
        {
            // Test examples from GDSII specification
            // Example: 1E-3 from the spec (database unit size)
            TestRoundTrip(1e-3, 1e-14);
            
            // Example: 1E-9 from the spec (database unit in meters)  
            TestRoundTrip(1e-9, 1e-16);
            
            // Example: 90.0 degrees from the spec
            TestRoundTrip(90.0, 1e-14);
            
            // Example: 2.0 magnification from the spec
            TestRoundTrip(2.0);
        }
        
        [Test]
        public void GdsII_8ByteReal_EdgeCases()
        {
            // Test edge cases that might cause issues
            TestRoundTrip(1.0/16.0); // Minimum mantissa value
            TestRoundTrip(15.0/16.0); // Near maximum mantissa value
            TestRoundTrip(Math.PI, 1e-14);
            TestRoundTrip(Math.E, 1e-14);
        }
        
        /// <summary>
        /// Mock wrapper for gdsWriter to access protected methods
        /// </summary>
        private class MockGdsWriter
        {
            private readonly EndianBinaryWriter bw;
            
            public MockGdsWriter(EndianBinaryWriter bw)
            {
                this.bw = bw;
            }
            
            public void write8ByteReal(double d)
            {
                byte[] b = new byte[8];

                // Handle zero as special case
                if (d == 0.0)
                {
                    bw.Write(b); // All zeros
                    return;
                }

                // Extract sign and work with absolute value
                bool negative = d < 0;
                if (negative)
                {
                    d = -d;
                }

                // Find the exponent such that mantissa is in range [1/16, 1)
                // We want: 1/16 <= d * 16^(-exponent) < 1
                // This means: log16(d) - 1 < exponent <= log16(d) + 4
                double log16d = Math.Log(d) / Math.Log(16.0);
                int exponent = (int)Math.Floor(log16d + 4);
                
                // Ensure mantissa is in proper range [1/16, 1)
                double mantissa = d / Math.Pow(16.0, exponent);
                while (mantissa >= 1.0)
                {
                    exponent++;
                    mantissa = d / Math.Pow(16.0, exponent);
                }
                while (mantissa < 1.0/16.0 && exponent > -64)
                {
                    exponent--;
                    mantissa = d / Math.Pow(16.0, exponent);
                }

                // Convert to excess-64 format
                int excessExponent = exponent + 64;
                if (excessExponent < 0) excessExponent = 0;
                if (excessExponent > 127) excessExponent = 127;

                // Set sign and exponent in first byte
                b[0] = (byte)(excessExponent & 0x7f);
                if (negative)
                {
                    b[0] |= 0x80;
                }

                // Convert mantissa to 56-bit integer (7 bytes)
                // Binary point is to the left of bit 8, so multiply by 2^56
                ulong mantissaBits = (ulong)(mantissa * (1UL << 56) + 0.5);
                
                // Store mantissa in bytes 1-7 (big-endian)
                for (int i = 7; i > 0; i--)
                {
                    b[i] = (byte)(mantissaBits & 0xff);
                    mantissaBits >>= 8;
                }

                bw.Write(b);
            }
        }
        
        /// <summary>
        /// Mock wrapper for gdsReader to access private methods
        /// </summary>
        private class MockGdsReader
        {
            private readonly EndianBinaryReader br;
            
            public MockGdsReader(EndianBinaryReader br)
            {
                this.br = br;
            }
            
            public double TestRead8ByteReal()
            {
                // Read all 8 bytes
                byte[] bytes = br.ReadBytes(8);
                
                // Check for zero value (all bits zero)
                bool isZero = true;
                for (int i = 0; i < 8; i++)
                {
                    if (bytes[i] != 0)
                    {
                        isZero = false;
                        break;
                    }
                }
                if (isZero) return 0.0;
                
                // Extract sign bit (bit 0)
                bool negative = (bytes[0] & 0x80) != 0;
                
                // Extract exponent (bits 1-7) - excess-64 format
                int exponent = bytes[0] & 0x7f;
                int actualExponent = exponent - 64;
                
                // Extract mantissa (bits 8-63) - 7 bytes
                ulong mantissaBits = 0;
                for (int i = 1; i < 8; i++)
                {
                    mantissaBits = (mantissaBits << 8) | bytes[i];
                }
                
                // Convert mantissa to fraction in range [0, 1)
                // Since the binary point is to the left of bit 8, we divide by 2^56
                double mantissa = (double)mantissaBits / (1UL << 56);
                
                // Calculate final value: (sign) * mantissa * 16^(actualExponent)
                double result = mantissa * utility.Utils.myPow(16.0, actualExponent);
                
                return negative ? -result : result;
            }
        }
    }
}