using Eto.Veldrid;
using System.Numerics;
using System.Runtime.InteropServices;
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
			Console.WriteLine("VeldridDriver: SetUpVeldrid called");
			
			try
			{
				Console.WriteLine($"VeldridDriver: Surface GraphicsDevice: {Surface?.GraphicsDevice}");
				Console.WriteLine($"VeldridDriver: Surface Swapchain: {Surface?.Swapchain}");
				
				CreateResources();
				Console.WriteLine("VeldridDriver: Resources created successfully");

				Ready = true;
				Console.WriteLine("VeldridDriver: Ready set to true");
			}
			catch (Exception ex)
			{
				Console.WriteLine($"VeldridDriver: Error in SetUpVeldrid: {ex.Message}");
				Console.WriteLine($"VeldridDriver: Stack trace: {ex.StackTrace}");
				throw;
			}
		}

		private void CreateResources()
		{
			Console.WriteLine("VeldridDriver: CreateResources called");
			
			// Load GLSL shader source code
			byte[] vertexShaderBytes = LoadShaderBytes(ShaderStages.Vertex);
			byte[] fragmentShaderBytes = LoadShaderBytes(ShaderStages.Fragment);

			Console.WriteLine($"VeldridDriver: Loaded vertex shader: {vertexShaderBytes.Length} bytes");
			Console.WriteLine($"VeldridDriver: Loaded fragment shader: {fragmentShaderBytes.Length} bytes");

			// Convert byte arrays to strings for GLSL compilation
			string vertexShaderCode = System.Text.Encoding.UTF8.GetString(vertexShaderBytes);
			string fragmentShaderCode = System.Text.Encoding.UTF8.GetString(fragmentShaderBytes);

			Console.WriteLine($"VeldridDriver: Using graphics backend: {Surface!.GraphicsDevice!.BackendType}");

			// Veldrid.SPIRV is used to cross-compile GLSL to the target backend
			CrossCompileOptions? options = new();
			switch (Surface.GraphicsDevice.BackendType)
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

			// Compile GLSL source code to SPIR-V bytecode
			Console.WriteLine("VeldridDriver: Compiling GLSL to SPIRV...");
			SpirvCompilationResult vertexResult = SpirvCompilation.CompileGlslToSpirv(
				vertexShaderCode, 
				"vertex.glsl", 
				ShaderStages.Vertex, 
				new GlslCompileOptions());
			
			SpirvCompilationResult fragmentResult = SpirvCompilation.CompileGlslToSpirv(
				fragmentShaderCode, 
				"fragment.glsl", 
				ShaderStages.Fragment, 
				new GlslCompileOptions());

			Console.WriteLine($"VeldridDriver: Vertex SPIRV compilation succeeded: {vertexResult.SpirvBytes.Length} bytes");
			Console.WriteLine($"VeldridDriver: Fragment SPIRV compilation succeeded: {fragmentResult.SpirvBytes.Length} bytes");

			// Use the compiled SPIR-V bytecode to create shaders
			ShaderDescription vertex = new(ShaderStages.Vertex, vertexResult.SpirvBytes, "main");
			ShaderDescription fragment = new(ShaderStages.Fragment, fragmentResult.SpirvBytes, "main");
			
			Console.WriteLine("VeldridDriver: Creating shaders from SPIRV...");
			Shader[] shaders = factory.CreateFromSpirv(vertex, fragment, options);
			Console.WriteLine($"VeldridDriver: Successfully created {shaders.Length} shaders");

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

			Console.WriteLine("VeldridDriver: Creating pipelines...");
			create_pipelines(ref factory, ref viewMatrixLayout, ref modelMatrixLayout, ref shaders, ref vertexLayout);

			CommandList = factory.CreateCommandList();
			Console.WriteLine("VeldridDriver: CreateResources completed successfully");
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

		private byte[] LoadShaderBytes(ShaderStages stage)
		{
			string name = $"VertexColor-{stage.ToString().ToLowerInvariant()}.450.glsl";
			string full = $"Eto.VeldridSurface.shaders.{name}";

			Console.WriteLine($"VeldridDriver: Loading shader resource: {full}");

			using (Stream? stream = GetType().Assembly.GetManifestResourceStream(full))
			{
				if (stream == null)
				{
					Console.WriteLine($"VeldridDriver: Failed to load shader resource: {full}");
					throw new InvalidOperationException($"Could not load shader resource: {full}");
				}

				Console.WriteLine($"VeldridDriver: Successfully loaded shader resource, length: {stream.Length}");
				
				using (BinaryReader reader = new(stream))
				{
					return reader.ReadBytes((int)stream.Length);
				}
			}
		}
	}
