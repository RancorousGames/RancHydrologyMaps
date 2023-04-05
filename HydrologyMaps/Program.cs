// See https://aka.ms/new-console-template for more information

using System.Diagnostics;

namespace HydrologyMaps;

public static class Program
{
    public static void Main(string[] args)
    {
        int width = 512;
        int height = 512;

        HydrologyParameters parameters = new HydrologyParameters();
        parameters.SpaceBetweenRiverMouthCandidates = 65;
        parameters.MinDistanceBetweenRiverNodes = 10;
        parameters.InterNodeDist = 15;
        parameters.NodeExpansionMaxTries = 10;
        parameters.MaxNodePriority = 10;
        parameters.ContinueProbability = 0.1;
        parameters.MinAngleBetweenBranches = 40;
        parameters.AsymmetricBranchProbability = 0.3;
        parameters.SymmetricBranchProbability = 0.6;
        HydrologyMapGen hydrologyMapGen = new HydrologyMapGen(parameters);
        HydrologyRenderer hydrologyMapRender = new HydrologyRenderer();

        int seed = 5; // 326765686; //120087508 441937323;// 200756544; //1291693340;  // 341179573
        //    for (int i = 0; i < 100; i++)// (seed == -1 ? 100 : 1); i++)
        bool stepByStep = false;
        for (int i = 0; i < (seed == -1 || stepByStep ? 200 : 1); i++)
        {
            var map = hydrologyMapGen.GenerateIsland(width, height, seed, (stepByStep ? i : 300));

            hydrologyMapRender.SaveMapAsPNG(map, "./heightmap.png", false, true,  false, parameters);

            // open heightmap
            var ps =  Process.GetProcesses();

            var ps1 = ps.Where(x => x.ProcessName.Contains("mage") || x.ProcessName.Contains("heightmap")).ToList();
            
            if (!Process.GetProcessesByName("ImageGlass").Any())
            {
                ProcessStartInfo startInfo = new ProcessStartInfo("cmd.exe", "/C start \"\" ./heightmap.png");
                Process.Start(startInfo);
            }

            Thread.Sleep(1000);
        }
    }
}