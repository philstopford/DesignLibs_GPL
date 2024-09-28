using Eto.Drawing;
using Eto.Forms;
using Eto.Veldrid;
using System.Numerics;
using Veldrid;

namespace VeldridEto;

public partial class VeldridDriver
{
	// Core Veldrid stuff. This shouldn't need to be messed with
	private VeldridSurface? _surface;

	public VeldridSurface? Surface
	{
		get => _surface;
		init
		{
			_surface = value;

			Surface!.Draw += (sender, e) => Draw();
		}
	}

	private uint[]? gridIndices;
	private uint[]? axesIndices;

	private uint[]? linesIndices;
	private uint[]? tessIndices;
	private uint[]? polyIndices;
	private uint[]? pointsIndices;

	private float axisZ;
	private float gridZ;

	private DeviceBuffer? GridVertexBuffer;
	private DeviceBuffer? GridIndexBuffer;
	private DeviceBuffer? AxesVertexBuffer;
	private DeviceBuffer? AxesIndexBuffer;

	private DeviceBuffer? VertexBuffer { get; set; }
	private DeviceBuffer? IndexBuffer { get; set; }
	private DeviceBuffer? ModelBuffer { get; set; }

	private DeviceBuffer? ViewBuffer;

	private DeviceBuffer? LinesVertexBuffer;
	private DeviceBuffer? PointsVertexBuffer;
	private DeviceBuffer? PolysVertexBuffer;
	private DeviceBuffer? TessVertexBuffer;
	private DeviceBuffer? LinesIndexBuffer;
	private DeviceBuffer? PointsIndexBuffer;
	private DeviceBuffer? PolysIndexBuffer;
	private DeviceBuffer? TessIndexBuffer;

	private Pipeline? PointsPipeline;
	private Pipeline? LinePipeline;
	private Pipeline? LinesPipeline;
	private Pipeline? FilledPipeline;

	private Matrix4x4 ModelMatrix { get; set; } = Matrix4x4.Identity;
	private Matrix4x4 ViewMatrix;
	private ResourceSet? ModelMatrixSet { get; set; }

	private ResourceSet? ViewMatrixSet;

	private CommandList? CommandList { get; set; }

	/*
	private Shader VertexShader { get; set; }
	private Shader FragmentShader { get; set; }
	*/

	private bool Ready = false;
	private PointF savedLocation;

	private const float pointWidth = 0.50f;
	private bool hasFocus;
	private bool keyHandlerApplied;

	// Use for drag handling.
	public bool dragging { get; set; }
	private float x_orig;
	private float y_orig;

	private DateTime CurrentTime;
	private DateTime PreviousTime = DateTime.Now;

	private ContextMenu? menu;
}
