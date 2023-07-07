package ModelFinal;

import HAL.Gui.GridWindow;
import HAL.Tools.FileIO;

import java.io.IOException;

public class InitialStroma extends ExampleCell {
    public static void main(String[]args) throws IOException {

        // log start time of algorithm
        long startTime = System.nanoTime();

        // number of iterations
        int timesteps = 30000;

        // declare time and stroma population storage arrays
        double [] timeVector = new double[timesteps];
        double [] stromaPopulation = new double[timesteps];
        double [] meanStromaAge = new double[timesteps];

        // set up animation
        GridWindow win = new GridWindow(x, y, 1);
        ExampleGrid model = new ExampleGrid(x, y);

        // name input folder
        String inputFolder = "C:\\MyFiles\\PHD Data\\";

        // name results folder
        String resultsFolder = "C:\\MyFiles\\PHD Data\\Data\\";

        int nVtemp = 600;

        // set up data recording
        FileIO homeostaticStromaCount = new FileIO(resultsFolder +
                "stromaHomeostasisCountDensityTest_" + nVtemp + ".csv", "w");
        FileIO homeostaticStroma = new FileIO(resultsFolder +
                "stromaHomeostasisDensityTest_" + nVtemp + ".csv", "w");

        // read vessel initialisation
        int [][] vesselSetUp = model.readVesselInitialisation(inputFolder +
                "\\VesselLocations\\vesselLocationsDensityTest_Q1.csv", nVtemp);

        // initialise vessel cells
        model.initialiseVesselCells(vesselSetUp, nVtemp, VESSEL);

        // fill region with stroma cells exclude vessel sites
        for (int i = 0; i < x * y; i++) {
            if (model.GetAgent(i) != null) {
                // vessel site
            } else {
                model.initialStromaCell(i, STROMA, timestepPerDay);
            }
        }

        model.drawModelAgents(win);

        // save initial to png
        //win.ToPNG("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\homeostasisinit.png");
        //win.ToPNG("C:\\Users\\2146974\\OneDrive - Swansea University\\Code\\HAL-master\\homeostasisinit.png");

        for (int i = 0; i < timesteps; i++) {
            win.TickPause(10);

            // create time counter
            int timeCounter = i;


            // collect stroma population data
            double timeCount = i * 1/timestepPerDay;
            double stromaCn = model.Pop() - nV;
            //double meanStroma = model.averageStromaAge(STROMA);
            timeVector[i] = timeCount;
            stromaPopulation[i] = stromaCn;
            //meanStromaAge[i] = (10 * meanStroma);

            // print out message every thousand iterations
            if (i % 1000 == 0) {
                System.out.println("iteration: " + i);
            }

            // age cells
            model.ageCells(STROMA, ACTIVATEDSTROMA, CANCER, UNREACTIVESTROMA);

            // divide cells
            model.divideCells(stromaNCI, STROMA, ACTIVATEDSTROMA, stromaIMT, CANCER, cancerIMT, proliferationCancerSite, UNREACTIVESTROMA, vesselSetUp);

            // stroma cells dying
            model.stromaDieCells(stromaDieProb, STROMA, ACTIVATEDSTROMA, UNREACTIVESTROMA);

            // increment timestep
            model.IncTick();

            // draw
            model.drawModelAgents(win);

            // save final set up as png
 /*           if (i == timesteps - 1) {
                win.ToPNG("C:\\Users\\amymm\\Documents\\PHD Data\\Data\\homeostasisfin.png");
            }*/
        }

        // record count, posistion and age of each stroma cell in homeostatic tissue
        int stromaHomeostasisCn = (int)model.populationType(STROMA);
        homeostaticStromaCount.Write(Integer.toString(stromaHomeostasisCn));
        model.stromaLocations(STROMA, UNREACTIVESTROMA, homeostaticStroma);
        homeostaticStromaCount.Close();
        homeostaticStroma.Close();

        // plot stroma population
        model.plotStromaPopulation(timeVector, stromaPopulation, meanStromaAge, timesteps);

        // log end time of algorithm
        long endTime = System.nanoTime();
        // time elapsed
        long timeElapsed = endTime - startTime;
        long timeSeconds = timeElapsed / 1000000000;
        long timeMinutes = timeSeconds / 60;

        System.out.println("Simulation finished!");
        System.out.println("Time taken: " + timeSeconds + " seconds");
        System.out.println("Time taken: " + timeMinutes + " minutes");


    }
}
