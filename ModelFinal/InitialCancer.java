package ModelFinal;

import HAL.Gui.GifMaker;
import HAL.Gui.GridWindow;
import HAL.Tools.FileIO;

import java.io.IOException;
import java.util.Arrays;

public class InitialCancer extends ExampleCell {
    public static void main(String[] args) throws IOException {

        // log start time of algorithm
        long startTime = System.nanoTime();

        // set up animation
        GridWindow agentLayer = new GridWindow(x, y, 1);
        GridWindow drugLayer = new GridWindow(x, y, 1);
        GridWindow proliferationLayer = new GridWindow(x, y, 1);
        ExampleGrid model = new ExampleGrid(x, y);

        // declare animation as gif
        /*GifMaker growCancerGif = new GifMaker("C:\\Users\\amymm\\OneDrive - Swansea University\\Code" +
                "\\HAL-master\\GrowCancer.gif", 10, false);

*/

        // name input folder
        String inputFolder = "C:\\MyFiles\\PHD Data\\";

        // name results folder
        String resultsFolder = "C:\\MyFiles\\PHD Data\\Data\\";

        int nVtemp = 600;

        /*GifMaker growCancerProliferationGif = new GifMaker("C:\\Users\\amymm\\Documents\\PHD Data" +
                "\\Data\\GrowCancerProliferation.gif", 1, false);*/

        // declare data storage
        FileIO locationsStroma = new FileIO(resultsFolder + "locationsStromaDensityTest_" + nVtemp + ".csv", "w");
        FileIO countInitialStroma = new FileIO(resultsFolder + "stromaInitialCountDensityTest_" + nVtemp + ".csv", "w");
        FileIO locationsCancer = new FileIO(resultsFolder + "locationsCancerDensityTest_" + nVtemp + ".csv", "w");
        FileIO countInitialCancer = new FileIO(resultsFolder + "cancerInitialCountDensityTest_" + nVtemp + ".csv", "w");
        FileIO proliferationSignalValues = new FileIO(resultsFolder + "valuesProliferationSignalDensityTest_" + nVtemp + ".csv", "w");

        // read vessel initialisation
        int[][] vesselSetUp = model.readVesselInitialisation(inputFolder + "VesselLocations" +
                "\\vesselLocationsDensityTest_Q1.csv", nVtemp);

        // initialise vessel cells
        model.initialiseVesselCells(vesselSetUp, nVtemp, VESSEL);

        // read homeostatic stroma set up
        // number of initial stroma cells from homeostatic run
        int homeostaticStromaCount = model.countHomeostaticStroma(inputFolder +
                "StromaHomeostasis\\stromaHomeostasisCountDensityTest_" + nVtemp + ".csv");
        // determine position and age of initial stroma population
        double [][] homeostaticStroma = model.homeostaticStromaSetUp(inputFolder +
                "StromaHomeostasis\\stromaHomeostasisDensityTest_" + nVtemp + ".csv", homeostaticStromaCount);

        // initialise stroma cells
        model.initialiseStroma(homeostaticStroma, STROMA, UNREACTIVESTROMA, reactiveStromaProportion);

        // initialise cancer cell
        model.initialCancerEvent(x, y, nVtemp, vesselSetUp, QUIESCENTCANCER, cancerIMT, proliferationCancerSite);

        // draw model
        model.drawModelLayers(agentLayer, drugLayer, proliferationLayer);

        // store frame in gif
//        growCancerGif.AddFrame(agentLayer);
//        growCancerProliferationGif.AddFrame(proliferationLayer);

//        agentLayer.ToPNG("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\growTumourInit.png");

        // declare time counter
        int j = 0;

        // main loop
        while (cancerDetected == false && cancerEliminated == false) {
            // time between animation frames
            agentLayer.TickPause(1);
            drugLayer.TickPause(1);
            proliferationLayer.TickPause(1);

            // create time counter
            int timeCounter = j;

            // collect cancer data
            double cancerCn = model.populationType(CANCER);

            // print out message every thousand iterations
            if (j % 1000 == 0) {
                System.out.println("iteration: " + j);
                System.out.println("cancer population: " + cancerCn);
            }

            // age cells
            model.ageCells(STROMA, ACTIVATEDSTROMA, CANCER, UNREACTIVESTROMA);

            // divide cells
            model.divideCells(stromaNCI, STROMA, ACTIVATEDSTROMA, stromaIMT, CANCER, cancerIMT, proliferationCancerSite, UNREACTIVESTROMA, vesselSetUp);

            int [] treatmentTimestep = model.treatmentScheduleTimesteps(timestepPerDay, treatmentDaysDrug, holidayDaysDrug);

            // stroma cells dying
            model.stromaDieCells(stromaDieProb, STROMA, ACTIVATEDSTROMA, UNREACTIVESTROMA);

            // cancer cell proliferation update and age cancer cells
            model.cancerDieProliferationUpdate(CANCER, QUIESCENTCANCER, cancerEliminated);

            // update PDEGrid values
            model.updatePDEValues(x, y, timeCounter, treatmentTimestep[0], treatmentTimestep[1], VESSEL,
                    drugDiffusionCoefficientTimestep, drugConcentrationVessel, drugRemovalRateVesselTimestep,
                    proliferationDiffusionCoefficientTimestep, CANCER, QUIESCENTCANCER, autocrineProliferationSignalProductionTimestep,
                    STROMA, ACTIVATEDSTROMA, paracrineProliferationSignalProductionTimestep,
                    proliferationDegradationDrugTimestep, initialTreatmentDays, timestepPerDay);

            // increment timestep
            model.IncTick();

            // check for cancer elimination
/*            if (cancerEliminated == true) {
                break;
            }
*/
            // check if cancer cell population has reached detection threshold
            double cancerPop = model.populationType(CANCER);
            double quiescentCancerPop = model.populationType(QUIESCENTCANCER);
            double totalCancerPop = cancerPop + quiescentCancerPop;
            if (cancerDetected == false) {
                if (totalCancerPop >= cancerDetectionThreshold) {
                    cancerDetected = true;
                    System.out.println("cancer detected iteration: " + j);
                }
            }

            // increment timestep
            model.IncTick();

            // draw
            model.drawModelLayers(agentLayer, drugLayer, proliferationLayer);

            // store frame in gif
//            growCancerGif.AddFrame(agentLayer);
            //growCancerProliferationGif.AddFrame(proliferationLayer);

            // save image as png
/*            if (j == 2000) {
                agentLayer.ToPNG("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\growTumour2000.png");
            }

            // save image as png
            if (j == 4000) {
                agentLayer.ToPNG("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\growTumour4000.png");
            }

            // save image as png
            if (j == 6000) {
                agentLayer.ToPNG("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\growTumour6000.png");
            }
*/
            // increase time counter
            j++;
        }

        //growCancerProliferationGif.Close();

        // store stroma and cancer location data
        int stromaInitCn = (int)model.populationType(STROMA);
        int unreactiveStromaInitCn = (int)model.populationType(UNREACTIVESTROMA);
        int totalStromaInitCn = stromaInitCn + unreactiveStromaInitCn;
        countInitialStroma.Write(Integer.toString(totalStromaInitCn));
        countInitialStroma.Close();
        model.stromaLocations(STROMA, UNREACTIVESTROMA, locationsStroma);
        locationsStroma.Close();
        int cancerInitCn = (int)model.populationType(CANCER);
        int quiescentCancerInitCn = (int)model.populationType(QUIESCENTCANCER);
        int totalCancerInitCn = cancerInitCn + quiescentCancerInitCn;
        countInitialCancer.Write(Integer.toString(totalCancerInitCn));
        countInitialCancer.Close();
        model.cancerLocations(CANCER, QUIESCENTCANCER, locationsCancer);
        locationsCancer.Close();
        model.finalProliferationSignal(proliferationSignalValues);
        proliferationSignalValues.Close();

        //agentLayer.ToPNG("C:\\Users\\amymm\\Documents\\PHD Data\\GrowCancer\\growTumourAgentsFin.png");
        //drugLayer.ToPNG("C:\\Users\\amymm\\Documents\\PHD Data\\GrowCancer\\growTumourDrugFin.png");
        //proliferationLayer.ToPNG("C:\\Users\\amymm\\Documents\\PHD Data\\GrowCancer\\growTumourProliferationFin.png");

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
