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

        HydrologyMapGen hydrologyMapGen = new HydrologyMapGen();

        for (int i = 0; i < 1000; i++)
        {
            var map = hydrologyMapGen.GenerateIsland(width, height);

            hydrologyMapGen.SaveMapAsPNG(map, ".heightmap.png");
            Thread.Sleep(1500);
        }
    }



}

