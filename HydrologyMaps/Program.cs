// See https://aka.ms/new-console-template for more information

using System.Diagnostics;
using RancHydrologyMaps;

namespace HydrologyMaps;

public static class Program
{
    public static void Main(string[] _)
    {
        int width = 512; // only square heightmaps supported

        HydrologyParameters parameters = new HydrologyParameters();
        parameters.SpaceBetweenRiverMouthCandidates = 100;
        parameters.MinDistanceBetweenRiverNodes = 25;
        parameters.InterNodeDist = 15;
        parameters.NodeExpansionMaxTries = 10;
        parameters.MaxNodePriority = 10;
        parameters.MinAngleBetweenBranches = 40;
        parameters.ContinueProbability = 0.4;
        parameters.AsymmetricBranchProbability = 0.3;
        parameters.RiverCarveBankSlope = 3f;
        parameters.SymmetricBranchProbability = 0.3;
        parameters.RiverCarveDepth = 0.3f;
        parameters.IslandShapeChaos = 0.5f;
        HydrologyMapGen hydrologyMapGen = new HydrologyMapGen(parameters);

        int seed = 5; // If seed is -1 then a random map is generated AND the program will continuously loop through maps
        bool stepByStep = false;
        bool saveMultipleFiles = false;
        for (int i = 0; i < (seed == -1 || stepByStep ? 200 : 1); i++)
        {
            var map = hydrologyMapGen.GenerateIsland(width, seed, (stepByStep ? i : 300));

            string filename = "./heightmap" + (saveMultipleFiles ? i : "") + ".png";
            HydrologyRenderer.SaveMapAsPNG(map, filename, 
                false, true, true, true);

            // open heightmap in preferred image viewer if not already open
            string preferredImageViewer = "ImageGlass";
            if (!Process.GetProcessesByName(preferredImageViewer).Any())
            {
                ProcessStartInfo startInfo = new ProcessStartInfo("cmd.exe", "/C start \"\" ./heightmap.png");
                Process.Start(startInfo);
            }

            Thread.Sleep(1000);
        }
    }
}