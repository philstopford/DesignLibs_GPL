using System;
using System.Collections.Generic;
using System.Numerics;
using Veldrid;

namespace VeldridEto;

public class ProgressivePolygonRenderer
{
	private readonly OVPSettings _settings;
	
	// Batch tracking for incremental updates
	private readonly List<VertexPositionColor> _progressivePolyList = new();
	private readonly List<VertexPositionColor> _progressivePointsList = new();
	private readonly List<VertexPositionColor> _progressiveTessList = new();
	private readonly List<uint> _progressivePolyIndices = new();
	private readonly List<uint> _progressivePointsIndices = new();
	private readonly List<uint> _progressiveTessIndices = new();
	
	// Track processed counts for each type
	private int _processedForegroundPolys = 0;
	private int _processedBackgroundPolys = 0;
	private int _processedTessellatedPolys = 0;
	
	public event Action<RenderUpdateEventArgs>? RenderUpdateReady;
	
	public ProgressivePolygonRenderer(OVPSettings settings)
	{
		_settings = settings ?? throw new ArgumentNullException(nameof(settings));
	}
	
	public void Reset()
	{
		_progressivePolyList.Clear();
		_progressivePointsList.Clear();
		_progressiveTessList.Clear();
		_progressivePolyIndices.Clear();
		_progressivePointsIndices.Clear();
		_progressiveTessIndices.Clear();
		
		_processedForegroundPolys = 0;
		_processedBackgroundPolys = 0;
		_processedTessellatedPolys = 0;
	}
	
	public RenderBatchData ProcessBatch(PolygonBatch batch)
	{
		var batchData = new RenderBatchData();
		
		switch (batch.PolygonType)
		{
			case PolygonType.Foreground:
				ProcessForegroundBatch(batch, batchData);
				break;
			case PolygonType.Background:
				ProcessBackgroundBatch(batch, batchData);
				break;
			case PolygonType.Tessellated:
				ProcessTessellatedBatch(batch, batchData);
				break;
		}
		
		// Trigger render update
		RenderUpdateReady?.Invoke(new RenderUpdateEventArgs
		{
			BatchData = batchData,
			TotalVertices = _progressivePolyList.Count + _progressivePointsList.Count + _progressiveTessList.Count,
			BatchType = batch.PolygonType,
			IsLastBatch = batch.IsLastBatch
		});
		
		return batchData;
	}
	
	private void ProcessForegroundBatch(PolygonBatch batch, RenderBatchData batchData)
	{
		var tessPolyListCount = _settings.tessPolyList?.Count ?? 0;
		var polyZStep = 1.0f / Math.Max(1, tessPolyListCount + (_settings.polyList?.Count ?? 0) + 1);
		
		var batchPolyVertices = new List<VertexPositionColor>();
		var batchPointVertices = new List<VertexPositionColor>();
		var batchPolyIndices = new List<uint>();
		var batchPointIndices = new List<uint>();
		
		foreach (var poly in batch.Polygons)
		{
			float alpha = poly.alpha;
			if (_settings.drawFilled())
			{
				alpha = 1.0f;
			}
			
			float polyZ = (tessPolyListCount * polyZStep) + (_processedForegroundPolys * polyZStep);
			int polyLength = poly.poly.Length - 1;
			
			// Calculate offsets for this polygon
			uint polyIndexOffset = (uint)_progressivePolyList.Count + (uint)batchPolyVertices.Count;
			uint pointIndexOffset = (uint)_progressivePointsList.Count + (uint)batchPointVertices.Count;
			
			// Process polygon vertices
			for (int pt = 0; pt < polyLength; pt++)
			{
				// Add polygon line segments
				batchPolyVertices.Add(new VertexPositionColor(
					new Vector3(poly.poly[pt].X, poly.poly[pt].Y, polyZ),
					new RgbaFloat(poly.color.R, poly.color.G, poly.color.B, alpha)));
				
				batchPolyVertices.Add(new VertexPositionColor(
					new Vector3(poly.poly[pt + 1].X, poly.poly[pt + 1].Y, polyZ),
					new RgbaFloat(poly.color.R, poly.color.G, poly.color.B, alpha)));
				
				// Add indices for the line segment
				batchPolyIndices.Add(polyIndexOffset + (uint)(pt * 2));
				batchPolyIndices.Add(polyIndexOffset + (uint)(pt * 2 + 1));
				
				// Add point vertices if enabled
				if (_settings.drawPoints())
				{
					const float pointWidth = 2.0f;
					// Create 6 vertices for 2 triangles (quad) for each point
					AddPointQuad(batchPointVertices, batchPointIndices, 
						poly.poly[pt], pointWidth, polyZ, poly.color, alpha, pointIndexOffset + (uint)(pt * 6));
				}
			}
			
			_processedForegroundPolys++;
		}
		
		// Add to progressive lists
		_progressivePolyList.AddRange(batchPolyVertices);
		_progressivePointsList.AddRange(batchPointVertices);
		_progressivePolyIndices.AddRange(batchPolyIndices);
		_progressivePointsIndices.AddRange(batchPointIndices);
		
		// Update batch data
		batchData.PolyVertices = _progressivePolyList.ToArray();
		batchData.PolyIndices = _progressivePolyIndices.ToArray();
		batchData.PointVertices = _progressivePointsList.ToArray();
		batchData.PointIndices = _progressivePointsIndices.ToArray();
	}
	
	private void ProcessBackgroundBatch(PolygonBatch batch, RenderBatchData batchData)
	{
		var polyZStep = 1.0f / Math.Max(1, (_settings.bgPolyList?.Count ?? 0) + 1);
		
		var batchPolyVertices = new List<VertexPositionColor>();
		var batchPolyIndices = new List<uint>();
		
		foreach (var poly in batch.Polygons)
		{
			float alpha = poly.alpha;
			float polyZ = _processedBackgroundPolys * polyZStep;
			
			int polyLength = poly.poly.Length - 1;
			uint indexOffset = (uint)_progressivePolyList.Count + (uint)batchPolyVertices.Count;
			
			for (int pt = 0; pt < polyLength; pt++)
			{
				batchPolyVertices.Add(new VertexPositionColor(
					new Vector3(poly.poly[pt].X, poly.poly[pt].Y, polyZ),
					new RgbaFloat(poly.color.R, poly.color.G, poly.color.B, alpha)));
				
				batchPolyVertices.Add(new VertexPositionColor(
					new Vector3(poly.poly[pt + 1].X, poly.poly[pt + 1].Y, polyZ),
					new RgbaFloat(poly.color.R, poly.color.G, poly.color.B, alpha)));
				
				batchPolyIndices.Add(indexOffset + (uint)(pt * 2));
				batchPolyIndices.Add(indexOffset + (uint)(pt * 2 + 1));
			}
			
			_processedBackgroundPolys++;
		}
		
		_progressivePolyList.AddRange(batchPolyVertices);
		_progressivePolyIndices.AddRange(batchPolyIndices);
		
		batchData.PolyVertices = _progressivePolyList.ToArray();
		batchData.PolyIndices = _progressivePolyIndices.ToArray();
	}
	
	private void ProcessTessellatedBatch(PolygonBatch batch, RenderBatchData batchData)
	{
		var batchTessVertices = new List<VertexPositionColor>();
		var batchTessIndices = new List<uint>();
		
		var polyZStep = 1.0f / Math.Max(1, (_settings.tessPolyList?.Count ?? 0) + 1);
		
		foreach (var poly in batch.Polygons)
		{
			float alpha = poly.alpha;
			float polyZ = _processedTessellatedPolys * polyZStep;
			
			uint indexOffset = (uint)_progressiveTessList.Count + (uint)batchTessVertices.Count;
			
			// Each tessellated polygon should have exactly 3 vertices (triangle)
			for (int pt = 0; pt < 3; pt++)
			{
				batchTessVertices.Add(new VertexPositionColor(
					new Vector3(poly.poly[pt].X, poly.poly[pt].Y, polyZ),
					new RgbaFloat(poly.color.R, poly.color.G, poly.color.B, alpha)));
			}
			
			// Add triangle indices
			batchTessIndices.Add(indexOffset);
			batchTessIndices.Add(indexOffset + 1);
			batchTessIndices.Add(indexOffset + 2);
			
			_processedTessellatedPolys++;
		}
		
		_progressiveTessList.AddRange(batchTessVertices);
		_progressiveTessIndices.AddRange(batchTessIndices);
		
		batchData.TessVertices = _progressiveTessList.ToArray();
		batchData.TessIndices = _progressiveTessIndices.ToArray();
	}
	
	private static void AddPointQuad(List<VertexPositionColor> vertices, List<uint> indices,
		Eto.Drawing.PointF point, float pointWidth, float z, Eto.Drawing.Color color, float alpha, uint startIndex)
	{
		var rgbaColor = new RgbaFloat(color.R, color.G, color.B, alpha);
		
		// Create quad vertices (2 triangles)
		vertices.Add(new VertexPositionColor(new Vector3(point.X - pointWidth / 2.0f, point.Y - pointWidth / 2.0f, z), rgbaColor));
		vertices.Add(new VertexPositionColor(new Vector3(point.X - pointWidth / 2.0f, point.Y + pointWidth / 2.0f, z), rgbaColor));
		vertices.Add(new VertexPositionColor(new Vector3(point.X + pointWidth / 2.0f, point.Y - pointWidth / 2.0f, z), rgbaColor));
		vertices.Add(new VertexPositionColor(new Vector3(point.X + pointWidth / 2.0f, point.Y - pointWidth / 2.0f, z), rgbaColor));
		vertices.Add(new VertexPositionColor(new Vector3(point.X - pointWidth / 2.0f, point.Y + pointWidth / 2.0f, z), rgbaColor));
		vertices.Add(new VertexPositionColor(new Vector3(point.X + pointWidth / 2.0f, point.Y + pointWidth / 2.0f, z), rgbaColor));
		
		// Add triangle indices for the quad
		for (uint i = 0; i < 6; i++)
		{
			indices.Add(startIndex + i);
		}
	}
	
	public RenderBatchData GetCurrentRenderData()
	{
		return new RenderBatchData
		{
			PolyVertices = _progressivePolyList.ToArray(),
			PolyIndices = _progressivePolyIndices.ToArray(),
			PointVertices = _progressivePointsList.ToArray(),
			PointIndices = _progressivePointsIndices.ToArray(),
			TessVertices = _progressiveTessList.ToArray(),
			TessIndices = _progressiveTessIndices.ToArray()
		};
	}
}

public class RenderBatchData
{
	public VertexPositionColor[]? PolyVertices { get; set; }
	public uint[]? PolyIndices { get; set; }
	public VertexPositionColor[]? PointVertices { get; set; }
	public uint[]? PointIndices { get; set; }
	public VertexPositionColor[]? TessVertices { get; set; }
	public uint[]? TessIndices { get; set; }
}

public class RenderUpdateEventArgs : EventArgs
{
	public RenderBatchData? BatchData { get; set; }
	public int TotalVertices { get; set; }
	public PolygonType BatchType { get; set; }
	public bool IsLastBatch { get; set; }
}