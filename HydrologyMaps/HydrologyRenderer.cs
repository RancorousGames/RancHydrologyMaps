using System.Drawing;
using System.Drawing.Imaging;

namespace HydrologyMaps;

public class HydrologyRenderer
{
    public void SaveMapAsPNG(HydrologyMap map, string outputPath, bool renderVoronoi, bool renderExtraPoints)
    {
        var heightmap = map.HeightMap;
        int width = heightmap.GetLength(0);
        int height = heightmap.GetLength(1);


        using (Bitmap bitmap = new Bitmap(width, height))
        {
            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    int grayValue = (int)(heightmap[x, y] * 255);
                    bitmap.SetPixel(x, y, Color.FromArgb(grayValue, grayValue, grayValue));
                }
            }

            
            if (renderVoronoi)
            {
                foreach (var edge in map.VoronoiEdges)
                {
                    DrawLineOnBitmap(bitmap, new Point((int)edge.X1, (int)edge.Y1),
                        new Point((int)edge.X2, (int)edge.Y2),
                        Color.Yellow, 1);

                    //Thread.Sleep(500);
                    //bitmap.Save(outputPath, ImageFormat.Png);
                }
            }

            foreach (var edge in map.RiverEdges)
            {
                DrawLineOnBitmap(bitmap, edge.P1.Point, edge.P2.Point, Color.Blue, 1);
            }

            if (renderExtraPoints)
            {
                foreach (var voronoiNode in map.VoronoiNodes)
                {
                    bitmap.SetPixel(voronoiNode.Point.X, voronoiNode.Point.Y, Color.Fuchsia);
                }
            }
            
            foreach (var riverMouth in map.AllRiverNodes)
            {
                bitmap.SetPixel(riverMouth.X, riverMouth.Y, Color.Red);
                if (riverMouth.Type == NodeType.Border)
                    bitmap.SetPixel(riverMouth.X + (int)(riverMouth.Direction.X*5), riverMouth.Y + (int)(riverMouth.Direction.Y*5), Color.Coral);

                bitmap.SetPixel(riverMouth.X + 1, riverMouth.Y - 1, Color.Red);
                bitmap.SetPixel(riverMouth.X + 1, riverMouth.Y + 1, Color.Red);

                bitmap.SetPixel(riverMouth.X - 1, riverMouth.Y - 1, Color.Red);

                bitmap.SetPixel(riverMouth.X - 1, riverMouth.Y + 1, Color.Red);
            }





            bitmap.Save(outputPath, ImageFormat.Png);
        }
    }


    static void DrawLineOnBitmap(Bitmap bitmap, Point coord1, Point coord2, Color lineColor, int lineWidth)
    {
        using (Graphics graphics = Graphics.FromImage(bitmap))
        {
            using (Pen pen = new Pen(lineColor, lineWidth))
            {
                graphics.DrawLine(pen, coord1.X, coord1.Y, coord2.X, coord2.Y);
            }
        }
    }
}