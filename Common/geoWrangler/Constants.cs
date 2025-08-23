namespace geoWrangler;

/// <summary>
/// Constants used throughout the geoWrangler library for geometric operations.
/// Provides precision values, tolerance settings, and common scaling factors
/// used in polygon processing and coordinate transformations.
/// </summary>
public static class Constants
{
    /// <summary>
    /// Decimal precision used for rounding operations in geometric calculations.
    /// Set to 4 decimal places to balance precision with computational efficiency.
    /// </summary>
    public const int roundingDecimalPrecision = 4;
    
    /// <summary>
    /// Default tolerance value for geometric comparisons and floating-point operations.
    /// Used to handle numerical precision issues in geometric calculations.
    /// </summary>
    public const double tolerance = 0.001;

    /// <summary>
    /// Scaling factor of 100 (1E2) commonly used for coordinate transformations.
    /// Often used to convert between different coordinate systems or units.
    /// </summary>
    public const double scalar_1E2 = 100;
    
    /// <summary>
    /// Inverse of the 1E2 scaling factor (1/100 = 0.01).
    /// Used for reverse transformations when scaling down coordinates.
    /// </summary>
    public const double scalar_1E2_inv = 1.0 / scalar_1E2;
    
    /// <summary>
    /// Scaling factor of 10000 (1E4) used for high-precision coordinate transformations.
    /// Provides greater precision when working with fine geometric details.
    /// </summary>
    public const double scalar_1E4 = 10000;
    
    /// <summary>
    /// Inverse of the 1E4 scaling factor (1/10000 = 0.0001).
    /// Used for reverse transformations when scaling down high-precision coordinates.
    /// </summary>
    public const double scalar_1E4_inv = 1.0 / scalar_1E4;
}
