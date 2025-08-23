using geoLib;
using Clipper2Lib;

namespace UnitTests;

/// <summary>
/// Tests for the GeoLib library, which provides basic geometric primitives
/// including points, arrays, matrices, and rectangles. These were designed
/// for multi-platform operation where System.Drawing was not available.
/// </summary>
public class GeoLibTests
{
    [Test]
    public static void MyVertex_Constructor_CreatesValidVertex()
    {
        // Test basic vertex creation
        MyVertex vertex = new(10.5, 20.3, typeDirection.left1, true, false, typeVertex.corner);
        
        Assert.That(vertex.X, Is.EqualTo(10.5).Within(0.001));
        Assert.That(vertex.Y, Is.EqualTo(20.3).Within(0.001));
        Assert.That(vertex.direction, Is.EqualTo(typeDirection.left1));
        Assert.That(vertex.vertical, Is.True);
        Assert.That(vertex.inner, Is.False);
        Assert.That(vertex.type, Is.EqualTo(typeVertex.corner));
        Assert.That(vertex.xBiasApplied, Is.False);
        Assert.That(vertex.yBiasApplied, Is.False);
    }
    
    [Test]
    public static void MyVertex_CopyConstructor_CreatesIdenticalVertex()
    {
        // Test copy constructor
        MyVertex original = new(15.0, 25.0, typeDirection.right1, false, true, typeVertex.center);
        original.xBiasApplied = true;
        original.yBiasApplied = true;
        
        MyVertex copy = new(original);
        
        Assert.That(copy.X, Is.EqualTo(original.X));
        Assert.That(copy.Y, Is.EqualTo(original.Y));
        Assert.That(copy.direction, Is.EqualTo(original.direction));
        Assert.That(copy.vertical, Is.EqualTo(original.vertical));
        Assert.That(copy.inner, Is.EqualTo(original.inner));
        Assert.That(copy.type, Is.EqualTo(original.type));
        // Note: bias properties are reset in copy constructor
        Assert.That(copy.xBiasApplied, Is.False);
        Assert.That(copy.yBiasApplied, Is.False);
    }
    
    [Test]
    public static void MyVertex_PropertyModification_WorksCorrectly()
    {
        MyVertex vertex = new(0, 0, typeDirection.up1, false, false, typeVertex.corner);
        
        // Test property modification
        vertex.X = 100.5;
        vertex.Y = 200.7;
        vertex.xBiasApplied = true;
        vertex.yBiasApplied = true;
        
        Assert.That(vertex.X, Is.EqualTo(100.5).Within(0.001));
        Assert.That(vertex.Y, Is.EqualTo(200.7).Within(0.001));
        Assert.That(vertex.xBiasApplied, Is.True);
        Assert.That(vertex.yBiasApplied, Is.True);
    }
    
    [Test]
    public static void MyRound_Properties_WorkCorrectly()
    {
        MyRound round = new();
        
        // Test property setting
        round.index = 5;
        round.verFace = 10;
        round.horFace = 15;
        round.MaxRadius = 25.5;
        round.direction = typeRound.inner;
        
        Assert.That(round.index, Is.EqualTo(5));
        Assert.That(round.verFace, Is.EqualTo(10));
        Assert.That(round.horFace, Is.EqualTo(15));
        Assert.That(round.MaxRadius, Is.EqualTo(25.5).Within(0.001));
        Assert.That(round.direction, Is.EqualTo(typeRound.inner));
    }
    
    [Test]
    public static void GeoLibArray_Properties_WorkCorrectly()
    {
        GeoLibArray array = new();
        
        // Test property setting
        array.point = new Point64(10, 20);
        
        Assert.That(array.point.X, Is.EqualTo(10));
        Assert.That(array.point.Y, Is.EqualTo(20));
        // Note: pitch and count are init-only properties, tested via constructor
    }
    
    [Test]
    public static void GeoLibArrayF_Properties_WorkCorrectly()
    {
        GeoLibArrayF arrayF = new();
        
        // Test property setting
        arrayF.point = new PointD(10.5, 20.7);
        arrayF.pitch = new PointD(5.5, 7.3);
        arrayF.count = new Point64(3, 4);
        
        Assert.That(arrayF.point.x, Is.EqualTo(10.5).Within(0.001));
        Assert.That(arrayF.point.y, Is.EqualTo(20.7).Within(0.001));
        Assert.That(arrayF.pitch.x, Is.EqualTo(5.5).Within(0.001));
        Assert.That(arrayF.pitch.y, Is.EqualTo(7.3).Within(0.001));
        Assert.That(arrayF.count.X, Is.EqualTo(3));
        Assert.That(arrayF.count.Y, Is.EqualTo(4));
    }
    
    [Test]
    public static void GeoLibRectangle_DefaultConstructor_CreatesEmptyRectangle()
    {
        GeoLibRectangle rect = new();
        
        // Default constructor should create zeroed rectangle
        Assert.That(rect.Location.X, Is.EqualTo(0));
        Assert.That(rect.Location.Y, Is.EqualTo(0));
        Assert.That(rect.Width, Is.EqualTo(0));
        Assert.That(rect.Height, Is.EqualTo(0));
    }
    
    [Test]
    public static void GeoLibRectangle_ParameterizedConstructor_CreatesCorrectRectangle()
    {
        GeoLibRectangle rect = new(10, 20, 100, 200);
        
        Assert.That(rect.Location.X, Is.EqualTo(10));
        Assert.That(rect.Location.Y, Is.EqualTo(20));
        Assert.That(rect.Width, Is.EqualTo(100));
        Assert.That(rect.Height, Is.EqualTo(200));
    }
    
    [Test]
    public static void GeoLibRectangle_CopyConstructor_CreatesIdenticalRectangle()
    {
        GeoLibRectangle original = new(15, 25, 150, 250);
        GeoLibRectangle copy = new(original);
        
        Assert.That(copy.Location.X, Is.EqualTo(original.Location.X));
        Assert.That(copy.Location.Y, Is.EqualTo(original.Location.Y));
        Assert.That(copy.Width, Is.EqualTo(original.Width));
        Assert.That(copy.Height, Is.EqualTo(original.Height));
    }
    
    [Test]
    public static void GeoLibRectangle_Offset_MovesRectangleCorrectly()
    {
        GeoLibRectangle rect = new(10, 20, 100, 200);
        Point64 offset = new(5, -3);
        
        rect.Offset(offset);
        
        Assert.That(rect.Location.X, Is.EqualTo(15)); // 10 + 5
        Assert.That(rect.Location.Y, Is.EqualTo(17)); // 20 + (-3)
        Assert.That(rect.Width, Is.EqualTo(100)); // Width unchanged
        Assert.That(rect.Height, Is.EqualTo(200)); // Height unchanged
    }
    
    [Test]
    public static void GeoLibVector3_IntConstructor_CreatesValidVector()
    {
        GeoLibVector3 vector = new(10, 20, 30);
        
        Assert.That(vector.x, Is.EqualTo(10.0).Within(0.001));
        Assert.That(vector.y, Is.EqualTo(20.0).Within(0.001));
        Assert.That(vector.z, Is.EqualTo(30.0).Within(0.001));
    }
    
    [Test]
    public static void GeoLibVector3_DoubleConstructor_CreatesValidVector()
    {
        GeoLibVector3 vector = new(10.5, 20.7, 30.3);
        
        Assert.That(vector.x, Is.EqualTo(10.5).Within(0.001));
        Assert.That(vector.y, Is.EqualTo(20.7).Within(0.001));
        Assert.That(vector.z, Is.EqualTo(30.3).Within(0.001));
    }
    
    [Test]
    public static void GeoLibVector3_CopyConstructor_CreatesIdenticalVector()
    {
        GeoLibVector3 original = new(15.5, 25.7, 35.3);
        GeoLibVector3 copy = new(original);
        
        Assert.That(copy.x, Is.EqualTo(original.x).Within(0.001));
        Assert.That(copy.y, Is.EqualTo(original.y).Within(0.001));
        Assert.That(copy.z, Is.EqualTo(original.z).Within(0.001));
    }
    
    [Test]
    public static void TypeDirectionEnum_HasAllExpectedValues()
    {
        // Test that all expected enum values exist
        Assert.That(Enum.IsDefined(typeof(typeDirection), typeDirection.left1), Is.True);
        Assert.That(Enum.IsDefined(typeof(typeDirection), typeDirection.right1), Is.True);
        Assert.That(Enum.IsDefined(typeof(typeDirection), typeDirection.up1), Is.True);
        Assert.That(Enum.IsDefined(typeof(typeDirection), typeDirection.down1), Is.True);
        Assert.That(Enum.IsDefined(typeof(typeDirection), typeDirection.tilt1), Is.True);
    }
    
    [Test]
    public static void TypeVertexEnum_HasAllExpectedValues()
    {
        // Test that all expected enum values exist
        Assert.That(Enum.IsDefined(typeof(typeVertex), typeVertex.corner), Is.True);
        Assert.That(Enum.IsDefined(typeof(typeVertex), typeVertex.center), Is.True);
    }
    
    [Test]
    public static void TypeRoundEnum_HasAllExpectedValues()
    {
        // Test that all expected enum values exist
        Assert.That(Enum.IsDefined(typeof(typeRound), typeRound.inner), Is.True);
        Assert.That(Enum.IsDefined(typeof(typeRound), typeRound.exter), Is.True);
    }
    
    [Test]
    public static void GeoLibRectangle_EdgeCases_HandleCorrectly()
    {
        // Test negative values
        GeoLibRectangle rect = new(-10, -20, 50, 60);
        Assert.That(rect.Location.X, Is.EqualTo(-10));
        Assert.That(rect.Location.Y, Is.EqualTo(-20));
        
        // Test zero dimensions
        GeoLibRectangle zeroRect = new(0, 0, 0, 0);
        Assert.That(zeroRect.Width, Is.EqualTo(0));
        Assert.That(zeroRect.Height, Is.EqualTo(0));
        
        // Test large offset
        Point64 largeOffset = new(1000000, -1000000);
        rect.Offset(largeOffset);
        Assert.That(rect.Location.X, Is.EqualTo(999990)); // -10 + 1000000
        Assert.That(rect.Location.Y, Is.EqualTo(-1000020)); // -20 + (-1000000)
    }
}