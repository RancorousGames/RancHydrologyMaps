using System.Drawing;
using System.Drawing.Imaging;

namespace HydrologyMaps;

public class HydrologyRenderer
{
    public void SaveMapAsPNG(HydrologyMap map, string outputPath, bool renderVoronoi, bool renderExtraPoints, bool drawRiverEdges,
        HydrologyParameters parameters)
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

               //    if (heightmap[x, y] == 0)
               //    {
               //        bitmap.SetPixel(x, y, Color.Black);
               //        continue;
               //    }
               //    
               //    var cell = map.GridCells[x, y];
               //    var val = Math.Max(0, cell.NearestNode.FlowRate - cell.DistToNode * 3);


               //    int col = Math.Max(0, (int)(255 - val * 3));
               //    bitmap.SetPixel(x, y, Color.FromArgb(col,col,col));
                }
            }


            if (renderVoronoi)
            {
                foreach (var edge in map.VoronoiEdges)
                {
                    DrawLineOnBitmap(bitmap, new Point((int)edge.X1, (int)edge.Y1),
                        new Point((int)edge.X2, (int)edge.Y2),
                        Color.Yellow, 1);
                }
            }

            if (drawRiverEdges)
            {
                foreach (var edge in map.RiverEdges)
                {
                    int prioBasedColor = (int)((166 / edge.FlowRate));
                    var color = Color.FromArgb(prioBasedColor, prioBasedColor, 255);
                    DrawLineOnBitmap(bitmap, edge.P1.Point, edge.P2.Point, color, 1);
                }
            }


            if (renderExtraPoints)
            {
                foreach (var voronoiNode in map.VoronoiNodes)
                {
                    bitmap.SetPixel(voronoiNode.Point.X, voronoiNode.Point.Y, Color.Fuchsia);
                }
            }

            
            foreach (var riverNode in map.AllRiverNodes)
            {
                var flowRateNormalized = (int)Math.Min(255, riverNode.FlowRate / 100 * 255);
                var col = Color.FromArgb(255, flowRateNormalized, flowRateNormalized);
                bitmap.SetPixel(riverNode.X, riverNode.Y,col);
                if (riverNode.Type == NodeType.Border)
                    bitmap.SetPixel(riverNode.X + (int)(riverNode.Direction.X * 5),
                        riverNode.Y + (int)(riverNode.Direction.Y * 5), Color.Coral);

                bitmap.SetPixel(riverNode.X + 1, riverNode.Y - 1, col);
                bitmap.SetPixel(riverNode.X + 1, riverNode.Y + 1, col);

                bitmap.SetPixel(riverNode.X - 1, riverNode.Y - 1, col);

                bitmap.SetPixel(riverNode.X - 1, riverNode.Y + 1, col);
            }

            foreach (var edge in HydrologyMapGen.DebugGraphEdges1)
            {
                DrawLineOnBitmap(bitmap, new Point((int)edge.X1, (int)edge.Y1),
                    new Point((int)edge.X2, (int)edge.Y2),
                    Color.Fuchsia, 1);
            }

            foreach (var edge in HydrologyMapGen.DebugGraphEdges2)
            {
                DrawLineOnBitmap(bitmap, new Point((int)edge.X1, (int)edge.Y1),
                    new Point((int)edge.X2, (int)edge.Y2),
                    Color.Chocolate, 1);
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