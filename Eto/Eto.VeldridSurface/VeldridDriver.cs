using Eto.Veldrid;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Text;
using Veldrid;
using Veldrid.SPIRV;
using VeldridEto;

namespace VeldridEto;
	public struct VertexPositionColor
	{
		public static uint SizeInBytes = (uint)Marshal.SizeOf(typeof(VertexPositionColor));

		public Vector3 Position;
		public RgbaFloat Color;

		public VertexPositionColor(Vector3 position, RgbaFloat color)
		{
			Position = position;
			Color = color;
		}
	}

	/// <summary>
	/// A class that controls rendering to a VeldridSurface.
	/// </summary>
	/// <remarks>
	/// VeldridSurface is only a basic control that lets you render to the screen
	/// using Veldrid. How exactly to do that is up to you; this driver class is
	/// only one possible approach, and in all likelihood not the most efficient.
	/// </remarks>
	public partial class VeldridDriver
	{

		public VeldridDriver(ref OVPSettings settings, ref VeldridSurface surface)
		{
			ovpSettings = settings;
			Surface = surface;
			addKeyHandlers();
			Clock.Interval = 1.0f / 60.0f;
			Clock.Elapsed += Clock_Elapsed!;
		}

		public void SetUpVeldrid()
		{
			CreateResources();

			Ready = true;
		}

		private void CreateResources()
		{
			// Check if swapchain is available - if not, defer resource creation
			if (Surface?.Swapchain == null)
			{
				Console.WriteLine("[DEBUG] Swapchain not available yet, deferring resource creation");
				// For Wayland systems, the swapchain might not be available immediately
				// This is normal and expected - we'll retry when VeldridInitialized is called again
				return;
			}
			
			// Add debugging for graphics device backend
			Console.WriteLine($"[DEBUG] Graphics device backend type: {Surface.GraphicsDevice?.BackendType}");
			Console.WriteLine($"[DEBUG] Expected backend type: {Surface.Backend}");
			
			// Check for backend mismatch
			if (Surface.GraphicsDevice?.BackendType != Surface.Backend)
			{
				Console.WriteLine($"[ERROR] Backend mismatch! Expected: {Surface.Backend}, Got: {Surface.GraphicsDevice?.BackendType}");
				throw new InvalidOperationException($"Graphics device backend mismatch. Expected: {Surface.Backend}, Got: {Surface.GraphicsDevice?.BackendType}");
			}
			
			// Veldrid.SPIRV is an additional library that complements Veldrid
			// by simplifying the development of cross-backend shaders, and is
			// currently the recommended approach to doing so:
			//
			//   https://veldrid.dev/articles/portable-shaders.html
			//
			// If you decide against using it, you can try out Veldrid developer
			// mellinoe's other project, ShaderGen, or drive yourself crazy by
			// writing and maintaining custom shader code for each platform.
			// Load GLSL source code from embedded resources
			string vertexShaderSource = LoadGlslSource(ShaderStages.Vertex);
			string fragmentShaderSource = LoadGlslSource(ShaderStages.Fragment);

			Console.WriteLine($"[DEBUG] Vertex shader source loaded: {vertexShaderSource?.Length ?? 0} characters");
			Console.WriteLine($"[DEBUG] Fragment shader source loaded: {fragmentShaderSource?.Length ?? 0} characters");

			if (string.IsNullOrEmpty(vertexShaderSource) || string.IsNullOrEmpty(fragmentShaderSource))
			{
				throw new InvalidOperationException("Failed to load shader source code");
			}

			CrossCompileOptions? options = new();
			switch (Surface!.GraphicsDevice!.BackendType)
			{
				// InvertVertexOutputY and FixClipSpaceZ address two major
				// differences between Veldrid's various graphics APIs, as
				// discussed here:
				//
				//   https://veldrid.dev/articles/backend-differences.html
				//
				// Note that the only reason those options are useful in this
				// example project is that the vertices being drawn are stored
				// the way Vulkan stores vertex data. The options will therefore
				// properly convert from the Vulkan style to whatever's used by
				// the destination backend. If you store vertices in a different
				// coordinate system, these may not do anything for you, and
				// you'll need to handle the difference in your shader code.
				case GraphicsBackend.Metal:
					options.InvertVertexOutputY = true;
					break;
				case GraphicsBackend.Direct3D11:
					options.InvertVertexOutputY = true;
					break;
				case GraphicsBackend.OpenGL:
					options.FixClipSpaceZ = true;
					options.InvertVertexOutputY = true;
					break;
				default:
					break;
			}

			ResourceFactory factory = Surface.GraphicsDevice.ResourceFactory;

			ResourceLayout viewMatrixLayout = factory.CreateResourceLayout(
				new ResourceLayoutDescription(
					new ResourceLayoutElementDescription(
						"ViewMatrix",
						ResourceKind.UniformBuffer,
						ShaderStages.Vertex)));

			ViewBuffer = factory.CreateBuffer(
				new BufferDescription(64, BufferUsage.UniformBuffer));

			ViewMatrixSet = factory.CreateResourceSet(new ResourceSetDescription(
				viewMatrixLayout, ViewBuffer));

			ShaderDescription vertex = new(ShaderStages.Vertex, Encoding.UTF8.GetBytes(vertexShaderSource), "main");
			ShaderDescription fragment = new(ShaderStages.Fragment, Encoding.UTF8.GetBytes(fragmentShaderSource), "main");
			
			Console.WriteLine($"[DEBUG] About to create shaders using CreateFromSpirv for backend: {Surface.GraphicsDevice.BackendType}");
			Console.WriteLine($"[DEBUG] ResourceFactory type: {factory.GetType().FullName}");
			
			Shader[] shaders;
			try
			{
				// Test if SPIRV compilation is available by attempting a simple operation
				Console.WriteLine("[DEBUG] Testing SPIRV library availability...");
				
				shaders = factory.CreateFromSpirv(vertex, fragment, options);
				Console.WriteLine($"[DEBUG] Successfully created {shaders.Length} shaders using SPIRV");
			}
			catch (DllNotFoundException dllEx)
			{
				Console.WriteLine($"[ERROR] SPIRV native library not found: {dllEx.Message}");
				Console.WriteLine("[ERROR] This is likely due to missing libveldrid-spirv.so");
				Console.WriteLine("[ERROR] Falling back to alternative shader creation method...");
				
				// For now, throw a more descriptive error
				throw new InvalidOperationException(
					"SPIRV shader compilation is not available. This is likely due to missing native dependencies (libveldrid-spirv.so). " +
					"Please ensure the Veldrid.SPIRV native library is properly installed.", dllEx);
			}
			catch (ArgumentNullException argEx)
			{
				Console.WriteLine($"[ERROR] Null argument in shader creation: {argEx.Message}");
				Console.WriteLine($"[ERROR] This might be due to backend mismatch or SPIRV compilation failure");
				
				// Check if the source is actually null
				Console.WriteLine($"[DEBUG] Vertex shader source null: {string.IsNullOrEmpty(vertexShaderSource)}");
				Console.WriteLine($"[DEBUG] Fragment shader source null: {string.IsNullOrEmpty(fragmentShaderSource)}");
				
				throw new InvalidOperationException(
					$"Shader creation failed with null bytes. Backend: {Surface.GraphicsDevice.BackendType}, " +
					$"Vertex source: {vertexShaderSource?.Length ?? 0}, Fragment source: {fragmentShaderSource?.Length ?? 0}", argEx);
			}
			catch (Exception ex)
			{
				Console.WriteLine($"[ERROR] Failed to create shaders: {ex.GetType().Name}: {ex.Message}");
				Console.WriteLine($"[ERROR] Stack trace: {ex.StackTrace}");
				
				// Provide more context in the error message
				throw new InvalidOperationException(
					$"Shader creation failed for {Surface.GraphicsDevice.BackendType} backend. " +
					$"This might be due to missing native dependencies or backend compatibility issues. " +
					$"Original error: {ex.GetType().Name}: {ex.Message}", ex);
			}

			ResourceLayout modelMatrixLayout = factory.CreateResourceLayout(
				new ResourceLayoutDescription(
					new ResourceLayoutElementDescription(
						"ModelMatrix",
						ResourceKind.UniformBuffer,
						ShaderStages.Vertex)));

			ModelBuffer = factory.CreateBuffer(
				new BufferDescription(64, BufferUsage.UniformBuffer));

			ModelMatrixSet = factory.CreateResourceSet(new ResourceSetDescription(
				modelMatrixLayout, ModelBuffer));

			VertexBuffer =
				factory.CreateBuffer(new BufferDescription(4 * VertexPositionColor.SizeInBytes,
					BufferUsage.VertexBuffer));
			IndexBuffer = factory.CreateBuffer(new BufferDescription(4 * sizeof(ushort), BufferUsage.IndexBuffer));

			// Veldrid.SPIRV, when cross-compiling to HLSL, will always produce
			// TEXCOORD semantics; VertexElementSemantic.TextureCoordinate thus
			// becomes necessary to let D3D11 work alongside Vulkan and OpenGL.
			//
			//   https://github.com/mellinoe/veldrid/issues/121
			//
			VertexLayoutDescription vertexLayout = new(
				new VertexElementDescription("Position", VertexElementSemantic.TextureCoordinate,
					VertexElementFormat.Float3),
				new VertexElementDescription("Color", VertexElementSemantic.TextureCoordinate,
					VertexElementFormat.Float4));

			create_pipelines(ref factory, ref viewMatrixLayout, ref modelMatrixLayout, ref shaders, ref vertexLayout);

			CommandList = factory.CreateCommandList();
			
			// VeldridDriver resources successfully created - no debug message needed
		}

		public void CompleteResourceCreation()
		{
			// This method can be called when the swapchain becomes available
			// to complete resource creation that was deferred
			if (Surface?.Swapchain != null && !Ready)
			{
				Console.WriteLine("[DEBUG] Completing deferred VeldridDriver resource creation");
				CreateResources();
				if (CommandList != null)
				{
					Ready = true;
				}
			}
		}

		private void create_pipelines(ref ResourceFactory factory, ref ResourceLayout viewMatrixLayout,
			ref ResourceLayout modelMatrixLayout, ref Shader[] shaders, ref VertexLayoutDescription vertexLayout)
		{
			LinePipeline = factory.CreateGraphicsPipeline(new GraphicsPipelineDescription
			{
				BlendState = BlendStateDescription.SINGLE_OVERRIDE_BLEND,
				DepthStencilState = new DepthStencilStateDescription(
					depthTestEnabled: true,
					depthWriteEnabled: true,
					comparisonKind: ComparisonKind.LessEqual),
				RasterizerState = new RasterizerStateDescription(
					cullMode: FaceCullMode.Back,
					fillMode: PolygonFillMode.Solid,
					frontFace: FrontFace.Clockwise,
					depthClipEnabled: true,
					scissorTestEnabled: false),
				PrimitiveTopology = PrimitiveTopology.LineList,
				ResourceLayouts = [viewMatrixLayout, modelMatrixLayout],
				ShaderSet = new ShaderSetDescription(
					vertexLayouts: [vertexLayout],
					shaders: shaders),
				Outputs = Surface!.Swapchain!.Framebuffer.OutputDescription
			});

			LinesPipeline = factory.CreateGraphicsPipeline(new GraphicsPipelineDescription
			{
				BlendState = BlendStateDescription.SINGLE_ALPHA_BLEND,
				DepthStencilState = new DepthStencilStateDescription(
					depthTestEnabled: false,
					depthWriteEnabled: false,
					comparisonKind: ComparisonKind.LessEqual),
				RasterizerState = new RasterizerStateDescription(
					cullMode: FaceCullMode.Back,
					fillMode: PolygonFillMode.Solid,
					frontFace: FrontFace.Clockwise,
					depthClipEnabled: false,
					scissorTestEnabled: false),
				PrimitiveTopology = PrimitiveTopology.LineList,
				ResourceLayouts = [viewMatrixLayout, modelMatrixLayout],
				ShaderSet = new ShaderSetDescription(
					vertexLayouts: [vertexLayout],
					shaders: shaders),
				Outputs = Surface.Swapchain.Framebuffer.OutputDescription
			});

			FilledPipeline = factory.CreateGraphicsPipeline(new GraphicsPipelineDescription
			{
				BlendState = BlendStateDescription.SINGLE_ALPHA_BLEND,
				DepthStencilState = new DepthStencilStateDescription(
					depthTestEnabled: false,
					depthWriteEnabled: false,
					comparisonKind: ComparisonKind.LessEqual),
				RasterizerState = new RasterizerStateDescription(
					cullMode: FaceCullMode.None,
					fillMode: PolygonFillMode.Solid,
					frontFace: FrontFace.CounterClockwise,
					depthClipEnabled: false,
					scissorTestEnabled: false),
				PrimitiveTopology = PrimitiveTopology.TriangleList,
				ResourceLayouts = [viewMatrixLayout, modelMatrixLayout],
				ShaderSet = new ShaderSetDescription(
					vertexLayouts: [vertexLayout],
					shaders: shaders),
				Outputs = Surface.Swapchain.Framebuffer.OutputDescription
			});
		}

		private string LoadGlslSource(ShaderStages stage)
		{
			string name = $"VertexColor-{stage.ToString().ToLowerInvariant()}.450.glsl";
			string full = $"Eto.VeldridSurface.shaders.{name}";

			Console.WriteLine($"[DEBUG] Loading GLSL source: {name} from resource: {full}");

			using (Stream? stream = GetType().Assembly.GetManifestResourceStream(full))
			{
				if (stream == null)
				{
					Console.WriteLine($"[ERROR] Could not find shader resource: {full}");
					throw new InvalidOperationException($"Could not find shader resource: {full}");
				}
				
				Console.WriteLine($"[DEBUG] Shader resource found, length: {stream.Length} bytes");
				
				using (StreamReader reader = new(stream))
				{
					string result = reader.ReadToEnd();
					Console.WriteLine($"[DEBUG] Successfully loaded {result.Length} characters of GLSL source for {stage} shader");
					return result;
				}
			}
		}
	}
