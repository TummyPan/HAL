package ModelFinal;

import HAL.Tools.FileIO;

import java.util.Arrays;

public class InitialVesselsUniform extends ModelFinal.ExampleCell {
    public static int[][] vesselPosition(int x, int y, int nV) {
        // square space per vessel
        double squareSpaceVessel = x * y / nV;
        int spaceVessel = (int) Math.sqrt(squareSpaceVessel);
        int xVessels = x / spaceVessel;
        int yVessels = y / spaceVessel;
        int noVessels = xVessels * yVessels;
        int xRemainder = x - (spaceVessel * xVessels);
        int yRemainder = y - (spaceVessel * yVessels);
        int xBuffer = (spaceVessel + xRemainder) / 2;
        int yBuffer = (spaceVessel + yRemainder) / 2;

        // create vessel site storage matrix
        int[][] vesselSites = new int[noVessels][2];

        // populate initial vessel site
        //vesselSites[0][0] = xBuffer;
        //vesselSites[0][1] = yBuffer;

// STILL FIXING THIS BIT!!!!
        // populate vessel site storage matrix
        for (int i = 0; i < noVessels; i++) {
            vesselSites[i][0] = xBuffer + (i / yVessels) * spaceVessel;
            vesselSites[i][1] = yBuffer + (i % xVessels) * spaceVessel;
        }
        return vesselSites;
    }

    public static void writeVesselLocations(int[][] vS, FileIO vesselLocations) {
        for (int i = 0; i < vS.length; i++) {
            vesselLocations.Write(vS[i][0] + "," + vS[i][1] + "\n");
        }
    }

    public static void main(String[] args) {
        // declare parameters

        // set up data recording
        FileIO initVessels = new FileIO("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master" +
                "\\vesselLocationsUniform.csv", "w");

        // determine vessel locations
        int[][] vS = InitialVesselsUniform.vesselPosition(x, y, nV);
        System.out.println(Arrays.deepToString(vS));

        // record data
        InitialVesselsUniform.writeVesselLocations(vS, initVessels);
        initVessels.Close();
    }
}
