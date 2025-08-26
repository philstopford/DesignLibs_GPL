using Eto.Veldrid;
using System.Numerics;
using System.Runtime.InteropServices;
using System.IO;
using System.Linq;
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
			// Check if libveldrid-spirv.so is available for SPIRV compilation
			bool nativeLibExists = File.Exists("libveldrid-spirv.so");
			Console.WriteLine($"[DEBUG] libveldrid-spirv.so exists in current directory: {nativeLibExists}");
			if (!nativeLibExists)
			{
				Console.WriteLine($"[DEBUG] Current directory: {Directory.GetCurrentDirectory()}");
				var files = Directory.GetFiles(".", "*.so").Take(5);
				Console.WriteLine($"[DEBUG] Available .so files: {string.Join(", ", files)}");
			}
			
			// Veldrid.SPIRV is an additional library that complements Veldrid
			// by simplifying the development of cross-backend shaders, and is
			// currently the recommended approach to doing so:
			//
			//   https://veldrid.dev/articles/portable-shaders.html
			//
			// Load GLSL source files and compile them to SPIRV
			string vertexShaderSource = LoadGlslSource(ShaderStages.Vertex);
			string fragmentShaderSource = LoadGlslSource(ShaderStages.Fragment);

			Console.WriteLine($"[DEBUG] Vertex shader source loaded: {vertexShaderSource.Length} characters");
			Console.WriteLine($"[DEBUG] Fragment shader source loaded: {fragmentShaderSource.Length} characters");

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

			Console.WriteLine($"[DEBUG] About to compile GLSL to SPIRV for backend: {Surface.GraphicsDevice.BackendType}");
			Console.WriteLine($"[DEBUG] ResourceFactory type: {factory.GetType().Name}");
			
			Shader[] shaders;
			try 
			{
				// Use SpirvCompilation to compile GLSL source to SPIRV and then create shaders
				var compilation = SpirvCompilation.CompileGlslToSpirv(
					vertexShaderSource,
					"main", 
					ShaderStages.Vertex,
					new GlslCompileOptions(false));
				
				var vertexSpirvBytes = compilation.SpirvBytes;
				Console.WriteLine($"[DEBUG] Vertex GLSL compiled to SPIRV: {vertexSpirvBytes.Length} bytes");
				
				compilation = SpirvCompilation.CompileGlslToSpirv(
					fragmentShaderSource,
					"main", 
					ShaderStages.Fragment,
					new GlslCompileOptions(false));
				
				var fragmentSpirvBytes = compilation.SpirvBytes;
				Console.WriteLine($"[DEBUG] Fragment GLSL compiled to SPIRV: {fragmentSpirvBytes.Length} bytes");
				
				// Now create shader descriptions with actual SPIRV bytecode
				ShaderDescription vertex = new(ShaderStages.Vertex, vertexSpirvBytes, "main", true);
				ShaderDescription fragment = new(ShaderStages.Fragment, fragmentSpirvBytes, "main", true);
				
				shaders = factory.CreateFromSpirv(vertex, fragment, options);
				Console.WriteLine($"[DEBUG] Successfully created {shaders.Length} shaders using GLSL->SPIRV compilation");
			}
			catch (Exception ex)
			{
				Console.WriteLine($"[DEBUG] GLSL->SPIRV shader compilation failed: {ex.GetType().Name}: {ex.Message}");
				Console.WriteLine($"[DEBUG] Exception stack trace: {ex.StackTrace}");
				
				// If GLSL compilation fails, provide specific guidance
				Console.WriteLine($"[DEBUG] This error suggests either:");
				Console.WriteLine($"[DEBUG] 1. libveldrid-spirv.so is not available or not loadable");
				Console.WriteLine($"[DEBUG] 2. The GLSL source contains syntax errors");
				Console.WriteLine($"[DEBUG] 3. The SPIRV compilation toolchain is not working properly");
				
				throw;
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
					Console.WriteLine($"[DEBUG] Shader resource not found: {full}");
					throw new InvalidOperationException($"Shader resource not found: {full}");
				}
				
				Console.WriteLine($"[DEBUG] Shader resource found, length: {stream.Length} bytes");
				using (StreamReader reader = new(stream))
				{
					string source = reader.ReadToEnd();
					Console.WriteLine($"[DEBUG] Successfully loaded {source.Length} characters of GLSL source for {stage} shader");
					return source;
				}
			}
		}
	}
