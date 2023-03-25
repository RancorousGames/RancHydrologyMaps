// See https://aka.ms/new-console-template for more information

using System.Drawing;
using System.Drawing.Imaging;
using HydrologyMaps;


public class Program
{
    public static void Main(string[] args)
    {
        int width = 512;
        int height = 512;

        HydrologyParameters parameters = new HydrologyParameters();
        parameters.SpaceBetweenRiverMouthCandidates = 110;
        HydrologyMapGen hydrologyMapGen = new HydrologyMapGen(parameters);
        HydrologyRenderer hydrologyMapRender = new HydrologyRenderer();

        int seed = 200756544; //1291693340;  // 341179573
        for (int i = 0; i < (seed == -1 ? 100 : 1); i++)
        {
            var map = hydrologyMapGen.GenerateIsland(width, height, seed);

            hydrologyMapRender.SaveMapAsPNG(map, ".heightmap.png");
            Thread.Sleep(1500);
        }
    }



}

