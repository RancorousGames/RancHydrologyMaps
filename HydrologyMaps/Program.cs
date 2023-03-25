// See https://aka.ms/new-console-template for more information

namespace HydrologyMaps;

public static class Program
{
    public static void Main(string[] args)
    {
        int width = 512;
        int height = 512;

        HydrologyParameters parameters = new HydrologyParameters();
        parameters.SpaceBetweenRiverMouthCandidates = 35;
        HydrologyMapGen hydrologyMapGen = new HydrologyMapGen(parameters);
        HydrologyRenderer hydrologyMapRender = new HydrologyRenderer();

        int seed = 326765686; //120087508 441937323;// 200756544; //1291693340;  // 341179573
        for (int i = 0; i < 200;i++)//(seed == -1 ? 100 : 1); i++)
        {
            Console.WriteLine(i);
            var map = hydrologyMapGen.GenerateIsland(width, height, seed, i);

            hydrologyMapRender.SaveMapAsPNG(map, ".heightmap.png", true, true);
            Thread.Sleep(250);
        }
    }



}