package ModelFinal;

import HAL.Tools.FileIO;

import java.util.Arrays;

public class InitialVesselsRandom extends ModelFinal.ExampleCell {
    public static int[][] vesselPosition(int x, int y, double sigmaMin, int nV) {
        // create vessel site storage matrix
        int[][] vesselSites = new int[nV][2];
        // create available positions matrix
        int[][] positions = new int[x * y][3];
        int ticker = 0;
        int j = 0;
        for (int i = 0; i < x * y; i++) {
            positions[i][0] = j % x;
            positions[i][1] = i % x;
            if (ticker % x == x - 1) {
                j = j + 1;
            }
            ticker = ticker + 1;
        }
        for (int k = 0; k < nV; k++) {
            // choose index for location of vessel
            double rN = Math.random() * (positions.length);
            int location = (int) rN;
            // create vessel site
            vesselSites[k][0] = positions[location][0];
            vesselSites[k][1] = positions[location][1];

            // determine points more than sigmaMin distance from vessels
            for (int i = 0; i < positions.length; i++) {
                if (Math.sqrt(Math.pow((positions[location][0] - positions[i][0]), 2)
                        + Math.pow((positions[location][1] - positions[i][1]), 2)) < sigmaMin) {
                    positions[i][2] = 0;
                } else {
                    positions[i][2] = 1;
                }
            }
            int counter = 0;
            for (int i = 0; i < positions.length; i++) {
                if (positions[i][2] == 1) {
                    counter++;
                }
            }
            int[][] newPositions = new int[counter][3];

            int idx = 0;

            for (int i = 0; i < positions.length; i++) {
                if (positions[i][2] == 1) {
                    newPositions[idx][0] = positions[i][0];
                    newPositions[idx][1] = positions[i][1];
                    idx++;
                }
            }
            positions = newPositions;
            if (positions.length == 0) {
                System.out.println("no positions");
                break;
            }
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
        FileIO initVessels = new FileIO("C:\\MyFiles\\PHD Data\\Data" +
                "\\vesselLocations.csv", "w");

        // determine vessel locations
        int[][] vS = InitialVesselsRandom.vesselPosition(x, y, sigmaMin, nV);
        System.out.println(Arrays.deepToString(vS));

        // record data
        InitialVesselsRandom.writeVesselLocations(vS, initVessels);
        initVessels.Close();
    }
}
