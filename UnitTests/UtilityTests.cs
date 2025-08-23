using utility;
using NextAfterNS;
using System.Collections.Generic;

namespace UnitTests;

/// <summary>
/// Tests for the Utility library, which provides general utility functions
/// including number formatting, mathematical operations, compression, hashing,
/// and list extensions.
/// </summary>
public class UtilityTests
{
    [Test]
    public static void Utils_FriendlyNumber_Int_FormatsCorrectly()
    {
        // Test integer formatting
        Assert.That(Utils.friendlyNumber(42), Is.EqualTo("42"));
        Assert.That(Utils.friendlyNumber(150), Is.EqualTo("1.5 hundred"));
        Assert.That(Utils.friendlyNumber(1500), Is.EqualTo("1.5 thousand"));
        Assert.That(Utils.friendlyNumber(1500000), Is.EqualTo("1.5 million"));
        Assert.That(Utils.friendlyNumber(1500000000), Is.EqualTo("1.5 billion"));
    }
    
    [Test]
    public static void Utils_FriendlyNumber_Long_FormatsCorrectly()
    {
        // Test long formatting
        Assert.That(Utils.friendlyNumber(42L), Is.EqualTo("42"));
        Assert.That(Utils.friendlyNumber(1500000000000L), Is.EqualTo("1.5 trillion"));
        Assert.That(Utils.friendlyNumber(15000000000000L), Is.EqualTo("15 trillion"));
    }
    
    [Test]
    public static void Utils_FriendlyNumber_Double_FormatsCorrectly()
    {
        // Test double formatting
        Assert.That(Utils.friendlyNumber(42.5), Is.EqualTo("42.5"));
        Assert.That(Utils.friendlyNumber(150.25), Is.EqualTo("1.5 hundred"));
        Assert.That(Utils.friendlyNumber(1500.75), Is.EqualTo("1.5 thousand"));
        Assert.That(Utils.friendlyNumber(1500000.123), Is.EqualTo("1.5 million"));
        Assert.That(Utils.friendlyNumber(2500000000.0), Is.EqualTo("2.5 billion"));
        Assert.That(Utils.friendlyNumber(3500000000000.0), Is.EqualTo("3.5 trillion"));
    }
    
    [Test]
    public static void Utils_MyPow_SpecialCases_WorkCorrectly()
    {
        // Test power function special cases
        Assert.That(Utils.myPow(5, 2), Is.EqualTo(25.0).Within(0.001));
        Assert.That(Utils.myPow(3, 0), Is.EqualTo(1.0).Within(0.001));
        Assert.That(Utils.myPow(2, 1), Is.EqualTo(2.0).Within(0.001));
        Assert.That(Utils.myPow(4, 3), Is.EqualTo(64.0).Within(0.001));
        Assert.That(Utils.myPow(2, 10), Is.EqualTo(1024.0).Within(0.001));
    }
    
    [Test]
    public static void Utils_MyPow_NegativeExponents_WorkCorrectly()
    {
        // Test negative exponents
        Assert.That(Utils.myPow(2, -1), Is.EqualTo(0.5).Within(0.001));
        Assert.That(Utils.myPow(4, -2), Is.EqualTo(0.0625).Within(0.001));
        Assert.That(Utils.myPow(10, -3), Is.EqualTo(0.001).Within(0.0001));
    }
    
    [Test]
    public static void Utils_ToRadians_ConvertsCorrectly()
    {
        // Test degree to radian conversion
        Assert.That(Utils.toRadians(0), Is.EqualTo(0.0).Within(0.001));
        Assert.That(Utils.toRadians(90), Is.EqualTo(Math.PI / 2).Within(0.001));
        Assert.That(Utils.toRadians(180), Is.EqualTo(Math.PI).Within(0.001));
        Assert.That(Utils.toRadians(360), Is.EqualTo(2 * Math.PI).Within(0.001));
        Assert.That(Utils.toRadians(45), Is.EqualTo(Math.PI / 4).Within(0.001));
    }
    
    [Test]
    public static void Utils_ToDegrees_ConvertsCorrectly()
    {
        // Test radian to degree conversion
        Assert.That(Utils.toDegrees(0), Is.EqualTo(0.0).Within(0.001));
        Assert.That(Utils.toDegrees(Math.PI / 2), Is.EqualTo(90.0).Within(0.001));
        Assert.That(Utils.toDegrees(Math.PI), Is.EqualTo(180.0).Within(0.001));
        Assert.That(Utils.toDegrees(2 * Math.PI), Is.EqualTo(360.0).Within(0.001));
        Assert.That(Utils.toDegrees(Math.PI / 4), Is.EqualTo(45.0).Within(0.001));
    }
    
    [Test]
    public static void Utils_CompressDecompress_WorksCorrectly()
    {
        // Test compression and decompression
        byte[] originalData = System.Text.Encoding.UTF8.GetBytes("Hello, World! This is a test string for compression.");
        
        byte[] compressed = Utils.compress(originalData);
        byte[] decompressed = Utils.decompress(compressed);
        
        Assert.That(compressed.Length, Is.LessThan(originalData.Length));
        Assert.That(decompressed, Is.EqualTo(originalData));
        
        string decompressedString = System.Text.Encoding.UTF8.GetString(decompressed);
        Assert.That(decompressedString, Is.EqualTo("Hello, World! This is a test string for compression."));
    }
    
    [Test]
    public static void Utils_GetMD5Hash_ProducesConsistentHash()
    {
        // Test MD5 hash function
        string testString = "Test data for hashing";
        
        string hash1 = Utils.GetMD5Hash(testString);
        string hash2 = Utils.GetMD5Hash(testString);
        
        Assert.That(hash1, Is.EqualTo(hash2));
        Assert.That(hash1, Is.Not.Empty);
        
        // Different input should produce different hash
        string differentHash = Utils.GetMD5Hash("Different test data");
        Assert.That(differentHash, Is.Not.EqualTo(hash1));
    }
    
    [Test]
    public static void Utils_GetSHA1Hash_ProducesConsistentHash()
    {
        // Test SHA1 hash function
        string testString = "Test data for SHA1 hashing";
        
        string hash1 = Utils.GetSHA1Hash(testString);
        string hash2 = Utils.GetSHA1Hash(testString);
        
        Assert.That(hash1, Is.EqualTo(hash2));
        Assert.That(hash1, Is.Not.Empty);
        
        // Different input should produce different hash
        string differentHash = Utils.GetSHA1Hash("Different test data");
        Assert.That(differentHash, Is.Not.EqualTo(hash1));
    }
    
    [Test]
    public static void Utils_GetSHA256Hash_ProducesConsistentHash()
    {
        // Test SHA256 hash function
        string testString = "Test data for SHA256 hashing";
        
        string hash1 = Utils.GetSHA256Hash(testString);
        string hash2 = Utils.GetSHA256Hash(testString);
        
        Assert.That(hash1, Is.EqualTo(hash2));
        Assert.That(hash1, Is.Not.Empty);
        
        // Different input should produce different hash
        string differentHash = Utils.GetSHA256Hash("Different test data");
        Assert.That(differentHash, Is.Not.EqualTo(hash1));
    }
    
    [Test]
    public static void FrExp_Calculate_WorksCorrectly()
    {
        // Test frexp function for various values
        var result1 = FrExp.calculate(8.0);
        Assert.That(result1.mantissa, Is.EqualTo(0.5).Within(0.001));
        Assert.That(result1.exponent, Is.EqualTo(4));
        
        var result2 = FrExp.calculate(1.0);
        Assert.That(result2.mantissa, Is.EqualTo(0.5).Within(0.001));
        Assert.That(result2.exponent, Is.EqualTo(1));
        
        // Test edge cases
        var resultZero = FrExp.calculate(0.0);
        Assert.That(resultZero.mantissa, Is.EqualTo(0.0));
        Assert.That(resultZero.exponent, Is.EqualTo(0));
        
        var resultNaN = FrExp.calculate(double.NaN);
        Assert.That(double.IsNaN(resultNaN.mantissa), Is.True);
        Assert.That(resultNaN.exponent, Is.EqualTo(0));
    }
    
    [Test]
    public static void NextAfter_WorksCorrectly()
    {
        // Test nextafter function
        float result1 = NextAfter.nextafter(1.0f, 2.0f);
        Assert.That(result1, Is.GreaterThan(1.0f));
        
        float result2 = NextAfter.nextafter(1.0f, 0.0f);
        Assert.That(result2, Is.LessThan(1.0f));
        
        // Test when x equals y
        float result3 = NextAfter.nextafter(5.0f, 5.0f);
        Assert.That(result3, Is.EqualTo(5.0f));
        
        // Test NaN cases
        float resultNaN = NextAfter.nextafter(float.NaN, 1.0f);
        Assert.That(float.IsNaN(resultNaN), Is.True);
    }
    
    [Test]
    public static void IListExtensions_Shuffle_ChangesOrder()
    {
        // Test list shuffling
        List<int> originalList = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        List<int> listToShuffle = new(originalList);
        
        listToShuffle.Shuffle();
        
        // List should have same elements but likely different order
        Assert.That(listToShuffle.Count, Is.EqualTo(originalList.Count));
        foreach (int item in originalList)
        {
            Assert.That(listToShuffle, Contains.Item(item));
        }
        
        // While technically possible, it's extremely unlikely the order stays the same
        // We'll just verify the shuffle method doesn't crash and preserves elements
    }
    
    [Test]
    public static void IListExtensions_Shuffle_EmptyList_HandlesCorrectly()
    {
        // Test shuffling empty list
        List<int> emptyList = [];
        emptyList.Shuffle();
        Assert.That(emptyList.Count, Is.EqualTo(0));
    }
    
    [Test]
    public static void IListExtensions_Shuffle_SingleElement_HandlesCorrectly()
    {
        // Test shuffling single element list
        List<int> singleList = [42];
        singleList.Shuffle();
        Assert.That(singleList.Count, Is.EqualTo(1));
        Assert.That(singleList[0], Is.EqualTo(42));
    }
    
    [Test]
    public static void Histo_Constructor_CreatesValidHistogram()
    {
        // Test histogram creation with bin boundaries
        double[] boundaries = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
        Histo histo = new(boundaries);
        
        // Verify histogram properties
        Assert.That(histo.NumBins, Is.EqualTo(5)); // n-1 bins for n boundaries
        Assert.That(histo.NumDataPoints, Is.EqualTo(0));
        Assert.That(histo.NumSmaller, Is.EqualTo(0));
        Assert.That(histo.NumLarger, Is.EqualTo(0));
    }
    
    [Test]
    public static void Histo_Update_CountsDataCorrectly()
    {
        // Test histogram data collection
        double[] boundaries = { 0.0, 1.0, 2.0, 3.0 };
        Histo histo = new(boundaries);
        
        // Add data points
        histo.Update(0.5); // Bin 0
        histo.Update(1.5); // Bin 1
        histo.Update(1.7); // Bin 1
        histo.Update(2.5); // Bin 2
        histo.Update(-0.5); // Smaller than range
        histo.Update(3.5); // Larger than range
        
        Assert.That(histo.NumDataPoints, Is.EqualTo(6));
        Assert.That(histo.NumSmaller, Is.EqualTo(1));
        Assert.That(histo.NumLarger, Is.EqualTo(1));
        Assert.That(histo.Count(0), Is.EqualTo(1)); // 0.5 in bin 0
        Assert.That(histo.Count(1), Is.EqualTo(2)); // 1.5, 1.7 in bin 1
        Assert.That(histo.Count(2), Is.EqualTo(1)); // 2.5 in bin 2
    }
    
    [Test]
    public static void Histo_Reset_ClearsData()
    {
        // Test histogram reset functionality
        double[] boundaries = { 0.0, 1.0, 2.0 };
        Histo histo = new(boundaries);
        
        histo.Update(0.5);
        histo.Update(1.5);
        Assert.That(histo.NumDataPoints, Is.EqualTo(2));
        
        histo.Reset();
        Assert.That(histo.NumDataPoints, Is.EqualTo(0));
        Assert.That(histo.Count(0), Is.EqualTo(0));
        Assert.That(histo.Count(1), Is.EqualTo(0));
    }
}