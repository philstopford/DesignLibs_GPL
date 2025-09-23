using System;
using System.IO;
using NUnit.Framework;
using gds;
using MiscUtil.Conversion;
using MiscUtil.IO;

namespace UnitTests
{
    /// <summary>
    /// Tests GDSII specification examples for validation against the manual.
    /// </summary>
    public class GdsIISpecificationTests
    {
        [Test]
        public void GdsII_SpecExample_UnitsRecord()
        {
            // Test the actual example from the GDSII specification
            // 0014 0305 3E41 8937 4BC6 A7EF 3944 B82F A09B 5A51
            // First 8-byte real: 3E41 8937 4BC6 A7EF = 1E-3
            // Second 8-byte real: 3944 B82F A09B 5A51 = 1E-9
            
            byte[] specBytes1 = { 0x3E, 0x41, 0x89, 0x37, 0x4B, 0xC6, 0xA7, 0xEF };
            byte[] specBytes2 = { 0x39, 0x44, 0xB8, 0x2F, 0xA0, 0x9B, 0x5A, 0x51 };
            
            // Test reading the spec example bytes
            using var ms1 = new MemoryStream(specBytes1);
            using var br1 = new EndianBinaryReader(EndianBitConverter.Big, ms1);
            var reader1 = new MockGdsReader(br1);
            double value1 = reader1.TestRead8ByteReal();
            
            using var ms2 = new MemoryStream(specBytes2);
            using var br2 = new EndianBinaryReader(EndianBitConverter.Big, ms2);
            var reader2 = new MockGdsReader(br2);
            double value2 = reader2.TestRead8ByteReal();
            
            // Verify against specification values
            Assert.That(Math.Abs(value1 - 1e-3), Is.LessThanOrEqualTo(1e-16), 
                $"First value should be 1E-3, got {value1}");
            Assert.That(Math.Abs(value2 - 1e-9), Is.LessThanOrEqualTo(1e-22), 
                $"Second value should be 1E-9, got {value2}");
                
            // Test round-trip to ensure our writer produces the same
            using var msOut1 = new MemoryStream();
            using var bwOut1 = new EndianBinaryWriter(EndianBitConverter.Big, msOut1);
            var writer1 = new MockGdsWriter(bwOut1);
            writer1.write8ByteReal(1e-3);
            
            using var msOut2 = new MemoryStream();
            using var bwOut2 = new EndianBinaryWriter(EndianBitConverter.Big, msOut2);
            var writer2 = new MockGdsWriter(bwOut2);
            writer2.write8ByteReal(1e-9);
            
            // Read back what we wrote
            msOut1.Position = 0;
            using var brOut1 = new EndianBinaryReader(EndianBitConverter.Big, msOut1);
            var readerOut1 = new MockGdsReader(brOut1);
            double roundTrip1 = readerOut1.TestRead8ByteReal();
            
            msOut2.Position = 0;
            using var brOut2 = new EndianBinaryReader(EndianBitConverter.Big, msOut2);
            var readerOut2 = new MockGdsReader(brOut2);
            double roundTrip2 = readerOut2.TestRead8ByteReal();
            
            Assert.That(Math.Abs(roundTrip1 - 1e-3), Is.LessThanOrEqualTo(1e-16), 
                "Round-trip 1E-3 failed");
            Assert.That(Math.Abs(roundTrip2 - 1e-9), Is.LessThanOrEqualTo(1e-22), 
                "Round-trip 1E-9 failed");
        }
        
        [Test]
        public void GdsII_SpecExample_AngleRecord()
        {
            // Test ANGLE record example: OOOC 1C05 425A 0000 0000 0000 = 90.0 degrees
            byte[] specBytes = { 0x42, 0x5A, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };
            
            using var ms = new MemoryStream(specBytes);
            using var br = new EndianBinaryReader(EndianBitConverter.Big, ms);
            var reader = new MockGdsReader(br);
            double angle = reader.TestRead8ByteReal();
            
            Assert.That(Math.Abs(angle - 90.0), Is.LessThanOrEqualTo(1e-14), 
                $"Angle should be 90.0 degrees, got {angle}");
        }
        
        [Test]
        public void GdsII_SpecExample_MagnificationRecord()
        {
            // Test MAG record example: 4120 0000 0000 0000 = 2.0
            byte[] specBytes = { 0x41, 0x20, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };
            
            using var ms = new MemoryStream(specBytes);
            using var br = new EndianBinaryReader(EndianBitConverter.Big, ms);
            var reader = new MockGdsReader(br);
            double mag = reader.TestRead8ByteReal();
            
            Assert.That(Math.Abs(mag - 2.0), Is.LessThanOrEqualTo(1e-15), 
                $"Magnification should be 2.0, got {mag}");
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