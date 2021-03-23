namespace geoWrangler
{
	// Experimental - not really used for anything at the moment, but it might be useful in the future.

	// A C# program to find maximal 
	// Bipartite matching. 
	// This code is contributed by nitin mittal. }
	// https://www.geeksforgeeks.org/maximum-bipartite-matching/

	class GFG
	{
		// M is number of applicants 
		// and N is number of jobs 
		static int M = 6;
		static int N = 6;

		// A DFS based recursive function 
		// that returns true if a matching 
		// for vertex u is possible 
		bool bpm(bool[,] bpGraph, int u,
				bool[] seen, int[] matchR)
		{
			// Try every job one by one 
			for (int v = 0; v < N; v++)
			{
				// If applicant u is interested 
				// in job v and v is not visited 
				if (bpGraph[u, v] && !seen[v])
				{
					// Mark v as visited 
					seen[v] = true;

					// If job 'v' is not assigned to 
					// an applicant OR previously assigned 
					// applicant for job v (which is matchR[v]) 
					// has an alternate job available. 
					// Since v is marked as visited in the above 
					// line, matchR[v] in the following recursive 
					// call will not get job 'v' again 
					if (matchR[v] < 0 || bpm(bpGraph, matchR[v],
											seen, matchR))
					{
						matchR[v] = u;
						return true;
					}
				}
			}
			return false;
		}

		bool bpm(bool[][] bpGraph, int u,
		bool[] seen, int[] matchR)
		{
			// Try every job one by one 
			for (int v = 0; v < N; v++)
			{
				// If applicant u is interested 
				// in job v and v is not visited 
				if (bpGraph[u][v] && !seen[v])
				{
					// Mark v as visited 
					seen[v] = true;

					// If job 'v' is not assigned to 
					// an applicant OR previously assigned 
					// applicant for job v (which is matchR[v]) 
					// has an alternate job available. 
					// Since v is marked as visited in the above 
					// line, matchR[v] in the following recursive 
					// call will not get job 'v' again 
					if (matchR[v] < 0 || bpm(bpGraph, matchR[v],
											seen, matchR))
					{
						matchR[v] = u;
						return true;
					}
				}
			}
			return false;
		}

		// Returns maximum number of 
		// matching from M to N 
		int maxBPM(bool[,] bpGraph)
		{
			// An array to keep track of the 
			// applicants assigned to jobs. 
			// The value of matchR[i] is the 
			// applicant number assigned to job i, 
			// the value -1 indicates nobody is assigned. 
			int[] matchR = new int[N];

			// Initially all jobs are available 
			for (int i = 0; i < N; ++i)
				matchR[i] = -1;

			// Count of jobs assigned to applicants 
			int result = 0;
			for (int u = 0; u < M; u++)
			{
				// Mark all jobs as not 
				// seen for next applicant. 
				bool[] seen = new bool[N];
				for (int i = 0; i < N; ++i)
					seen[i] = false;

				// Find if the applicant 
				// 'u' can get a job 
				if (bpm(bpGraph, u, seen, matchR))
					result++;
			}
			return result;
		}

		int maxBPM(bool[][] bpGraph)
		{
			// An array to keep track of the 
			// applicants assigned to jobs. 
			// The value of matchR[i] is the 
			// applicant number assigned to job i, 
			// the value -1 indicates nobody is assigned. 
			int[] matchR = new int[N];

			// Initially all jobs are available 
			for (int i = 0; i < N; ++i)
				matchR[i] = -1;

			// Count of jobs assigned to applicants 
			int result = 0;
			for (int u = 0; u < M; u++)
			{
				// Mark all jobs as not 
				// seen for next applicant. 
				bool[] seen = new bool[N];
				for (int i = 0; i < N; ++i)
					seen[i] = false;

				// Find if the applicant 
				// 'u' can get a job 
				if (bpm(bpGraph, u, seen, matchR))
					result++;
			}
			return result;
		}

		/*
		// Driver Code 
		public static void Main()
		{
			// Let us create a bpGraph shown  
			// in the above example 
			bool[,] bpGraph = new bool[,]
							  {{false, true, true,
							false, false, false},
						   {true, false, false,
							true, false, false},
						   {false, false, true,
							false, false, false},
						   {false, false, true,
							true, false, false},
						   {false, false, false,
							false, false, false},
						   {false, false, false,
							false, false, true}};
			GFG m = new GFG();
			Console.Write("Maximum number of applicants that can" +
									" get job is " + m.maxBPM(bpGraph));
		}
		*/
	}
}
