using System;
using System.Globalization;
using System.IO;
using System.IO.Compression;
using System.Runtime.Serialization;
using System.Security.Cryptography;

namespace utility
{
    public static class Utils
    {
        public static string friendlyNumber(Int32 number)
        {
            return pFriendlyNumber(number);
        }

        public static string friendlyNumber(Int64 number)
        {
            return pFriendlyNumber(number);
        }

        public static string friendlyNumber(double number)
        {
            return pFriendlyNumber(number);
        }

        static string pFriendlyNumber(double number)
        {
            string ret;
            if (number >= 1E12)
            {
                ret = (number / 1E12).ToString("0.##") + " trillion";
            }
            else if (number >= 1E9)
            {
                ret = (number / 1E9).ToString("0.##") + " billion";
            }
            else if (number >= 1E6)
            {
                ret = (number / 1E6).ToString("0.##") + " million";
            }
            else if (number >= 1E3)
            {
                ret = (number / 1E3).ToString("0.##") + " thousand";
            }
            else if (number >= 1E2)
            {
                ret = (number / 1E2).ToString("0.##") + " hundred";
            }
            else
            {
                ret = number.ToString(CultureInfo.InvariantCulture);
            }

            return ret;
        }

        public static double myPow(double num, int exp)
        {
            // Micro-optimization
            if (exp == 2)
            {
                return num * num;
            }

            double result = 1.0;
            bool invertResult = false;

            if (exp < 0)
            {
                invertResult = true;
                exp = Math.Abs(exp);
            }

            while (exp > 0)
            {
                if (exp % 2 == 1)
                    result *= num;
                exp >>= 1;
                num *= num;
            }

            if (invertResult)
            {
                result = 1.0f / result;
            }

            return result;
        }

        public static double toRadians(double deg)
        {
            return deg * (Math.PI / 180.0f);
        }

        public static double toDegrees(double rad)
        {
            return rad / (Math.PI / 180.0f);
        }

        public static byte[] compress(byte[] data)
        {
            MemoryStream compressedStream = new MemoryStream();
            GZipStream zipStream = new GZipStream(compressedStream, CompressionMode.Compress);
            zipStream.Write(data, 0, data.Length);
            zipStream.Close();
            return compressedStream.ToArray();
        }

        public static byte[] decompress(byte[] data)
        {
            MemoryStream compressedStream = new MemoryStream(data);
            GZipStream zipStream = new GZipStream(compressedStream, CompressionMode.Decompress);
            MemoryStream resultStream = new MemoryStream();
            zipStream.CopyTo(resultStream);
            return resultStream.ToArray();
        }

        // Below is from : https://alexmg.com/compute-any-hash-for-any-object-in-c/

        /// <summary>
        ///     Gets a hash of the current instance.
        /// </summary>
        /// <typeparam name="T">
        ///     The type of the Cryptographic Service Provider to use.
        /// </typeparam>
        /// <param name="instance">
        ///     The instance being extended.
        /// </param>
        /// <returns>
        ///     A base 64 encoded string representation of the hash.
        /// </returns>
        public static string GetHash<T>(this object instance) where T : HashAlgorithm, new()
        {
            T cryptoServiceProvider = new T();
            return computeHash(instance, cryptoServiceProvider);
        }

        /// <summary>
        ///     Gets a key based hash of the current instance.
        /// </summary>
        /// <typeparam name="T">
        ///     The type of the Cryptographic Service Provider to use.
        /// </typeparam>
        /// <param name="instance">
        ///     The instance being extended.
        /// </param>
        /// <param name="key">
        ///     The key passed into the Cryptographic Service Provider algorithm.
        /// </param>
        /// <returns>
        ///     A base 64 encoded string representation of the hash.
        /// </returns>
        public static string GetKeyedHash<T>(this object instance, byte[] key) where T : KeyedHashAlgorithm, new()
        {
            T cryptoServiceProvider = new T { Key = key };
            return computeHash(instance, cryptoServiceProvider);
        }

        /// <summary>
        ///     Gets a MD5 hash of the current instance.
        /// </summary>
        /// <param name="instance">
        ///     The instance being extended.
        /// </param>
        /// <returns>
        ///     A base 64 encoded string representation of the hash.
        /// </returns>
        public static string GetMD5Hash(this object instance)
        {
            return instance.GetHash<MD5CryptoServiceProvider>();
        }

        /// <summary>
        ///     Gets a SHA1 hash of the current instance.
        /// </summary>
        /// <param name="instance">
        ///     The instance being extended.
        /// </param>
        /// <returns>
        ///     A base 64 encoded string representation of the hash.
        /// </returns>
        public static string GetSHA1Hash(this object instance)
        {
            return instance.GetHash<SHA1CryptoServiceProvider>();
        }

        private static string computeHash<T>(object instance, T cryptoServiceProvider) where T : HashAlgorithm, new()
        {
            DataContractSerializer serializer = new DataContractSerializer(instance.GetType());
            using (MemoryStream memoryStream = new MemoryStream())
            {
                serializer.WriteObject(memoryStream, instance);
                cryptoServiceProvider.ComputeHash(memoryStream.ToArray());
                return Convert.ToBase64String(cryptoServiceProvider.Hash);
            }
        }
    }
}
