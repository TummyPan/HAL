package ModelFinal;

import HAL.Gui.GridWindow;

import java.io.IOException;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;

public class HomeostaticStromaTest extends ExampleCell {
    public static void main(String[]args) throws IOException {

        // log start time of algorithm
        long startTime = System.nanoTime();

        // log algorithm date and time of run
        LocalDateTime rawDateTime = LocalDateTime.now();
        DateTimeFormatter formatDateTime = DateTimeFormatter.ofPattern("dd-MM-yyyy_HH-mm-ss");
        String runDateTime = rawDateTime.format(formatDateTime);

        // name input folder
        String inputFolder = "C:\\MyFiles\\PHD Data\\";

        // name results folder
        String resultsFolder = "C:\\MyFiles\\PHD Data\\ThesisFigures\\";

        // number of iterations
        int timesteps = 30000;

        // declare time and stroma population storage arrays
        double [] timeVector = new double[timesteps];
        double [] stromaPopulation = new double[timesteps];

        // set up animation
        GridWindow win = new GridWindow(x, y, 1);
        ExampleGrid model = new ExampleGrid(x, y);

        // read vessel initialisation
        int [][] vesselSetUp = model.readVesselInitialisation(inputFolder + "\\VesselLocations\\vesselLocations.csv", nV);

        // initialise vessel cells
        model.initialiseVesselCells(vesselSetUp, nV, VESSEL);

        // read homeostatic stroma set up
        // number of initial stroma cells from homeostatic run
        int homeostaticStromaCount = model.countHomeostaticStroma(inputFolder + "StromaHomeostasis\\stromaHomeostasisCount.csv");
        // determine position and age of initial stroma population
        double[][] homeostaticStroma = model.homeostaticStromaSetUp(inputFolder + "StromaHomeostasis\\stromaHomeostasis.csv", homeostaticStromaCount);

        // initialise stroma cells
        model.initialiseStroma(homeostaticStroma, STROMA, UNREACTIVESTROMA, 1);

        // create two boxes one with no stroma the other all stroma
        model.stromaBoxes(x, y, STROMA);

        model.drawModelAgents(win);

        // save image as png
        //win.ToPNG("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\homeostasisBoxesInit.png");
        //win.ToPNG("C:\\Users\\2146974\\OneDrive - Swansea University\\Code\\HAL-master\\homeostasisBoxesInit.png");

        model.positionData(-1, 0, 0, x, y, runDateTime, resultsFolder);

        for (int i = 0; i < timesteps; i++) {
            win.TickPause(1);

            // collect stroma population data
            double timeCount = i * 1/timestepPerDay;
            double stromaCn = model.Pop() - nV;
            timeVector[i] = timeCount;
            stromaPopulation[i] = stromaCn;

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

            // collect data
            if (i % 100 == 0) {
                model.positionData(i, 0, 0, x, y, runDateTime, resultsFolder);
            }

        }

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
