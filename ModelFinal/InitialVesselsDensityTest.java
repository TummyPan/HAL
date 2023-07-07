package ModelFinal;

import HAL.Tools.FileIO;

import java.util.Arrays;

public class InitialVesselsDensityTest extends ModelFinal.ExampleCell {
    public static int[][] vesselPosition(int x, int y, double sigmaMin, int nV) {
        // create vessel site storage matrix
        int[][] vesselSites = new int[nV][2];
        // create available positions matrix
        int[][] positions = new int[x * y][3];
        int ticker = 0;
        int j = 0;
        for (int i = 0; i < x * y; i++) {
            positions[i][0] = j % x;
            positions[i][1] = i % y;
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
/*        FileIO initVessels = new FileIO("C:\\MyFiles\\PHD Data\\Data" +
                "\\vesselLocationsDensityTest.csv", "w");
*/
        FileIO initVesselsQ1 = new FileIO("C:\\MyFiles\\PHD Data\\VesselLocations" +
                "\\vesselLocationsDensityTest_Q1.csv", "w");
        FileIO initVesselsQ2 = new FileIO("C:\\MyFiles\\PHD Data\\VesselLocations" +
                "\\vesselLocationsDensityTest_Q2.csv", "w");
        FileIO initVesselsQ3 = new FileIO("C:\\MyFiles\\PHD Data\\VesselLocations" +
                "\\vesselLocationsDensityTest_Q3.csv", "w");
        FileIO initVesselsQ4 = new FileIO("C:\\MyFiles\\PHD Data\\VesselLocations" +
                "\\vesselLocationsDensityTest_Q4.csv", "w");

        int nVQ1 = 600;
        int nVQ2 = 200;
        int nVQ3 = 250;
        int nVQ4 = 300;

        double nVMultiplier = 830;

        double sigmaMin_nVQ1 = nVMultiplier * Math.sqrt((xDomain * yDomain) / nVQ1);
        double sigmaMin_nVQ2 = nVMultiplier * Math.sqrt((xDomain * yDomain) / nVQ2);
        double sigmaMin_nVQ3 = nVMultiplier * Math.sqrt((xDomain * yDomain) / nVQ3);
        double sigmaMin_nVQ4 = nVMultiplier * Math.sqrt((xDomain * yDomain) / nVQ4);

        // determine vessel locations
        int[][] vSQ1 = InitialVesselsDensityTest.vesselPosition(x, y, sigmaMin_nVQ1, nVQ1);
        //int[][] vSQ2 = InitialVesselsDensityTest.vesselPosition(x, y, sigmaMin_nVQ2, nVQ2);
        //int[][] vSQ3 = InitialVesselsDensityTest.vesselPosition(x, y, sigmaMin_nVQ3, nVQ3);
        //int[][] vSQ4 = InitialVesselsDensityTest.vesselPosition(x, y, sigmaMin_nVQ4, nVQ4);
        System.out.println(Arrays.deepToString(vSQ1));
        //System.out.println(Arrays.deepToString(vSQ2));
        //System.out.println(Arrays.deepToString(vSQ3));
        //System.out.println(Arrays.deepToString(vSQ4));



        /*int[][] vS = new int[nVLow + nVHigh][2];
        System.out.println(vS.length);
        for (int i = 0; i < vS.length; i++) {
            if (i < nVLow) {
                vS[i][0] = vSLow[i][0];
                vS[i][1] = vSLow[i][1];
            } else {
                int indexFix = i - nVLow;
                vS[i][0] = vSHigh[indexFix][0] + x / 2;
                vS[i][1] = vSHigh[indexFix][1];
            }
        }*/

        //System.out.println(Arrays.deepToString(vS));

        // record data
        InitialVesselsDensityTest.writeVesselLocations(vSQ1, initVesselsQ1);
        //InitialVesselsDensityTest.writeVesselLocations(vSQ2, initVesselsQ2);
        //InitialVesselsDensityTest.writeVesselLocations(vSQ3, initVesselsQ3);
        //InitialVesselsDensityTest.writeVesselLocations(vSQ4, initVesselsQ4);
        initVesselsQ1.Close();
        //initVesselsQ2.Close();
        //initVesselsQ3.Close();
        //initVesselsQ4.Close();
    }
}
