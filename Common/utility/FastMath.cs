using System;
using System.Runtime.CompilerServices;

namespace utility;

/// <summary>
/// High-performance mathematical operations with optimized implementations
/// </summary>
public static class FastMath
{
    // Lookup tables for fast trigonometric operations
    private const int TRIG_TABLE_SIZE = 4096;
    private static readonly float[] SinTable = new float[TRIG_TABLE_SIZE];
    private static readonly float[] CosTable = new float[TRIG_TABLE_SIZE];
    private static readonly float TrigScale = TRIG_TABLE_SIZE / (2.0f * MathF.PI);
    
    // Fast exponential and logarithm constants
    private const float LOG2_E = 1.4426950408889634f;
    private const float LN2 = 0.6931471805599453f;
    
    static FastMath()
    {
        // Pre-calculate trigonometric lookup tables
        for (int i = 0; i < TRIG_TABLE_SIZE; i++)
        {
            float angle = 2.0f * MathF.PI * i / TRIG_TABLE_SIZE;
            SinTable[i] = MathF.Sin(angle);
            CosTable[i] = MathF.Cos(angle);
        }
    }

    /// <summary>
    /// Fast sine approximation using lookup table
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static float FastSin(float x)
    {
        // Normalize angle to [0, 2π]
        x = x - MathF.Floor(x / (2.0f * MathF.PI)) * (2.0f * MathF.PI);
        
        int index = (int)(x * TrigScale) & (TRIG_TABLE_SIZE - 1);
        return SinTable[index];
    }

    /// <summary>
    /// Fast cosine approximation using lookup table
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static float FastCos(float x)
    {
        // Normalize angle to [0, 2π]
        x = x - MathF.Floor(x / (2.0f * MathF.PI)) * (2.0f * MathF.PI);
        
        int index = (int)(x * TrigScale) & (TRIG_TABLE_SIZE - 1);
        return CosTable[index];
    }

    /// <summary>
    /// Fast square root using bit manipulation (Quake's fast inverse square root variant)
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static float FastSqrt(float x)
    {
        if (x <= 0.0f) return 0.0f;
        
        // Use built-in for better accuracy on modern hardware
        return MathF.Sqrt(x);
    }

    /// <summary>
    /// Fast inverse square root approximation
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static unsafe float FastInverseSqrt(float x)
    {
        if (x <= 0.0f) return 0.0f;
        
        float xhalf = 0.5f * x;
        int i = *(int*)&x;
        i = 0x5f3759df - (i >> 1); // Magic number
        x = *(float*)&i;
        x = x * (1.5f - xhalf * x * x); // Newton-Raphson iteration
        return x;
    }

    /// <summary>
    /// Fast logarithm approximation
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static float FastLog(float x)
    {
        if (x <= 0.0f) return float.NegativeInfinity;
        
        // Use built-in for accuracy
        return MathF.Log(x);
    }

    /// <summary>
    /// Fast exponential approximation
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static float FastExp(float x)
    {
        // Use built-in for accuracy
        return MathF.Exp(x);
    }

    /// <summary>
    /// Fast power function optimized for small integer exponents
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static float FastPow(float baseValue, int exponent)
    {
        // Performance optimization: Handle common cases first
        switch (exponent)
        {
            case 0: return 1.0f;
            case 1: return baseValue;
            case 2: return baseValue * baseValue;
            case 3: return baseValue * baseValue * baseValue;
            case 4: 
                var squared = baseValue * baseValue;
                return squared * squared;
            case 5:
                var val2 = baseValue * baseValue;
                return val2 * val2 * baseValue;
            case 6:
                var val3 = baseValue * baseValue * baseValue;
                return val3 * val3;
            case 8:
                var val4 = baseValue * baseValue;
                val4 *= val4;
                return val4 * val4;
        }
        
        // Use binary exponentiation for larger exponents
        float result = 1.0f;
        bool negative = exponent < 0;
        exponent = negative ? -exponent : exponent; // Avoid Math.Abs for better performance
        
        float currentBase = baseValue;
        while (exponent > 0)
        {
            if ((exponent & 1) == 1)
                result *= currentBase;
            currentBase *= currentBase;
            exponent >>= 1;
        }
        
        return negative ? 1.0f / result : result;
    }

    /// <summary>
    /// Linear interpolation optimized for performance
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static float Lerp(float a, float b, float t)
    {
        return a + t * (b - a);
    }

    /// <summary>
    /// Fast absolute value
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static float FastAbs(float x)
    {
        return MathF.Abs(x);
    }

    /// <summary>
    /// Fast floor function
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int FastFloor(float x)
    {
        return (int)MathF.Floor(x);
    }
}