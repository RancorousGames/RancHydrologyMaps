// See https://aka.ms/new-console-template for more information

using System.Diagnostics;

namespace HydrologyMaps;

public static class Program
{
    public static void Main(string[] args)
    {
        int width = 512; // only square heightmaps supported

        HydrologyParameters parameters = new HydrologyParameters();
        parameters.SpaceBetweenRiverMouthCandidates = 100;
        parameters.MinDistanceBetweenRiverNodes = 25;
        parameters.InterNodeDist = 15;
        parameters.NodeExpansionMaxTries = 10;
        parameters.MaxNodePriority = 10;
        parameters.ContinueProbability = 0.1;
        parameters.MinAngleBetweenBranches = 40;
        parameters.AsymmetricBranchProbability = 0.3;
        parameters.SymmetricBranchProbability = 0.6;
        HydrologyMapGen hydrologyMapGen = new HydrologyMapGen(parameters);
        HydrologyRenderer hydrologyMapRender = new HydrologyRenderer();

        int seed = 1562471381;// 178145830;
        bool stepByStep = false;
        for (int i = 0; i < (seed == -1 || stepByStep ? 200 : 1); i++)
        {
            var map = hydrologyMapGen.GenerateIsland(width, seed, (stepByStep ? i : 300));

            hydrologyMapRender.SaveMapAsPNG(map, "./heightmap.png", false, true,  true, finalHeightmapRender: true, parameters);

            // open heightmap
            var ps =  Process.GetProcesses();

            var ps1 = ps.Where(x => x.ProcessName.Contains("mage") || x.ProcessName.Contains("heightmap")).ToList();
            
            if (!Process.GetProcessesByName("ImageGlass").Any())
            {
                ProcessStartInfo startInfo = new ProcessStartInfo("cmd.exe", "/C start \"\" ./heightmap.png");
                Process.Start(startInfo);
            }

            Thread.Sleep(300);
        }
    }
}