package ModelFinal;

import HAL.Gui.GridWindow;
import HAL.Tools.FileIO;

import java.io.IOException;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;

public class NeighboursInformationSameInitialCondition extends ExampleCell {
    public static void main(String[] args) throws IOException {

        // log start time of algorithm
        long startTime = System.nanoTime();

        // log algorithm date and time of run
        LocalDateTime rawDateTime = LocalDateTime.now();
        DateTimeFormatter formatDateTime = DateTimeFormatter.ofPattern("dd-MM-yyyy_HH-mm-ss");
        String runDateTime = rawDateTime.format(formatDateTime);

        // name results folder
        String inputFolder = "C:\\MyFiles\\PHD Data\\";

        // name results folder
        String resultsFolder = "C:\\MyFiles\\PHD Data\\GrowCancer\\";

        // set up animation
        GridWindow agentLayer = new GridWindow(x, y, 1);
        ExampleGrid model = new ExampleGrid(x, y);

        // declare passive and reactive stroma recorder
        FileIO locationsPassiveReactiveStroma = new FileIO(resultsFolder + "\\locationsPassiveReactiveStroma.csv", "w");

        // read vessel initialisation
        int[][] vesselSetUp = model.readVesselInitialisation(inputFolder + "\\VesselLocations\\vesselLocations.csv", nV);

        // initialise vessel cells
        model.initialiseVesselCells(vesselSetUp, nV, VESSEL);

        // determine position and age of initial stroma population
        int initialStromaCount = model.countInitialStroma(inputFolder + "\\GrowCancer\\stromaInitialCount.csv");
        double[][] initialStroma = model.initialStromaSetUp(inputFolder + "\\GrowCancer\\locationsStroma.csv", initialStromaCount);

        // determine position and age of initial cancer population
        int initialCancerCount = model.countInitialCancer(inputFolder + "\\GrowCancer\\cancerInitialCount.csv");
        double[][] initialCancer = model.initialCancerSetUp(inputFolder + "\\GrowCancer\\locationsCancer.csv", initialCancerCount);

        // initialise stroma cells
        model.initialiseStroma(initialStroma, STROMA, UNREACTIVESTROMA, reactiveStromaProportion);

        // initialise cancer
        model.initialiseCancer(initialCancer, CANCER, QUIESCENTCANCER, cellIDCount, vesselSetUp, clusterBoxSize, VESSEL, x, y);

        // record passive and reactive stroma initial conditions
        model.stromaPassiveReactiveLocations(STROMA, UNREACTIVESTROMA, locationsPassiveReactiveStroma);

        locationsPassiveReactiveStroma.Close();

    }
}


