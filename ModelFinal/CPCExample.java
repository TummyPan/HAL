package ModelFinal;

import HAL.Gui.GridWindow;
import HAL.Tools.FileIO;

import java.io.IOException;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;

public class CPCExample extends ExampleCell {
    public static void main(String[]args) throws IOException {

        // log start time of algorithm
        long startTime = System.nanoTime();

        // log algorithm date and time of run
        LocalDateTime rawDateTime = LocalDateTime.now();
        DateTimeFormatter formatDateTime = DateTimeFormatter.ofPattern("dd-MM-yyyy_HH-mm-ss");
        String runDateTime = rawDateTime.format(formatDateTime);

        // name results folder
        //String inputFolder = "C:\\MyFiles\\PHD Data\\";

        // name results folder
        String resultsFolder = "C:\\MyFiles\\PHD Data\\NeighbourhoodInvestigation\\CPCExample\\";

        for (int k = 0; k < 1; k++) {
            System.out.println("Run " + k);
            for (int i = 0; i < 1; i++) {

                // print simulation number
                System.out.println("    Simulation " + i);

                // number of iterations
                int timesteps = 10; // number of iterations of model

                int xVal = 100;
                int yVal = 100;

                // set up animation
                GridWindow agentLayer = new GridWindow(xVal, yVal, 3);
                ExampleGrid model = new ExampleGrid(xVal, yVal);

                // initialisation
                for (int j = 0; j < 4; j++) {
                    for (int l = 0; l < 4; l++) {
                        model.initialCellCPCExample((j + 1) * 20, (l + 1) * 20, QUIESCENTCANCER);
                    }
                }

                model.drawModelAgents(agentLayer);
                model.positionDataCPCExample(100, xVal, yVal, resultsFolder);

                // main loop
                for (int j = 0; j < timesteps; j++) {
                    // time between animation frames
                    agentLayer.TickPause(1);

                    System.out.println("iteration: " + j);

                    // count agents
                    double agentCn = model.populationType(QUIESCENTCANCER);

                    // age cells
                    model.CPCExampleDivideCells(agentCn);

                    // increment timestep
                    model.IncTick();

                    // draw
                    model.drawModelAgents(agentLayer);

                    model.positionDataCPCExample(j, xVal, yVal, resultsFolder);

                }
            }
        }

        // log end time of algorithm
        long endTime = System.nanoTime();
        // time elapsed
        long timeElapsed = endTime - startTime;
        long timeSeconds = timeElapsed / 1000000000;
        long timeMinutes = timeSeconds / 60;

        // indicate simulation has finished
        System.out.println("Simulation finished!");
        System.out.println("Time taken: " + timeSeconds + " seconds");
        System.out.println("Time taken: " + timeMinutes + " minutes");
    }
}


