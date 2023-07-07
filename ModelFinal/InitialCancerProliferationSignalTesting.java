package ModelFinal;

import HAL.Gui.GifMaker;
import HAL.Gui.GridWindow;
import HAL.Tools.FileIO;

import java.io.IOException;

public class InitialCancerProliferationSignalTesting extends ExampleCell {
    public static void main(String[] args) throws IOException {

        // log start time of algorithm
        long startTime = System.nanoTime();

        // set length of the extended time
        int extendedTimeLength = 10;

        // set up animation
        GridWindow agentLayer = new GridWindow(x, y, 8);
        GridWindow drugLayer = new GridWindow(x, y, 8);
        GridWindow proliferationLayer = new GridWindow(x, y, 8);
        ExampleGrid model = new ExampleGrid(x, y);

        // declare animation as gif
        /*GifMaker growCancerGif = new GifMaker("C:\\Users\\amymm\\OneDrive - Swansea University\\Code" +
                "\\HAL-master\\GrowCancer.gif", 10, false);

*/
        GifMaker growCancerProliferationGif = new GifMaker("C:\\Users\\amymm\\Documents\\PHD Data" +
                "\\Data\\GrowCancerProliferation.gif", 1, false);

        // declare data storage
        FileIO locationsStroma = new FileIO("C:\\Users\\amymm\\Documents\\PHD Data" +
                "\\Data\\locationsStroma.csv", "w");
        FileIO countInitialStroma = new FileIO("C:\\Users\\amymm\\Documents\\PHD Data" +
                "\\Data\\stromaInitialCount.csv", "w");
        FileIO locationsCancer = new FileIO("C:\\Users\\amymm\\Documents\\PHD Data" +
                "\\Data\\locationsCancer.csv", "w");
        FileIO countInitialCancer = new FileIO("C:\\Users\\amymm\\Documents\\PHD Data" +
                "\\Data\\cancerInitialCount.csv", "w");
        FileIO proliferationSignalValues = new FileIO("C:\\Users\\amymm\\Documents\\PHD Data" +
                "\\Data\\valuesProliferationSignal.csv", "w");

        // declare data storage
        FileIO locationsStromaEXTENDED = new FileIO("C:\\Users\\amymm\\Documents\\PHD Data" +
                "\\Data\\locationsStromaEXTENDED.csv", "w");
        FileIO countInitialStromaEXTENDED = new FileIO("C:\\Users\\amymm\\Documents\\PHD Data" +
                "\\Data\\stromaInitialCountEXTENDED.csv", "w");
        FileIO locationsCancerEXTENDED = new FileIO("C:\\Users\\amymm\\Documents\\PHD Data" +
                "\\Data\\locationsCancerEXTENDED.csv", "w");
        FileIO countInitialCancerEXTENDED = new FileIO("C:\\Users\\amymm\\Documents\\PHD Data" +
                "\\Data\\cancerInitialCountEXTENDED.csv", "w");
        FileIO proliferationSignalValuesEXTENDED = new FileIO("C:\\Users\\amymm\\Documents\\PHD Data" +
                "\\Data\\valuesProliferationSignalEXTENDED.csv", "w");

        // read vessel initialisation
        int[][] vesselSetUp = model.readVesselInitialisation("C:\\Users\\amymm\\Documents\\PHD Data\\VesselLocations" +
                "\\vesselLocationsSMALL.csv", nV);

        // initialise vessel cells
        model.initialiseVesselCells(vesselSetUp, nV, VESSEL);

        // read homeostatic stroma set up
        // number of initial stroma cells from homeostatic run
        int homeostaticStromaCount = model.countHomeostaticStroma("C:\\Users\\amymm\\Documents\\PHD Data\\StromaHomeostasis" +
                "\\stromaHomeostasisCountSMALL.csv");
        // determine position and age of initial stroma population
        double [][] homeostaticStroma = model.homeostaticStromaSetUp("C:\\Users\\amymm\\Documents\\PHD Data\\StromaHomeostasis" +
                "\\stromaHomeostasisSMALL.csv", homeostaticStromaCount);

        // initialise stroma cells
        model.initialiseStroma(homeostaticStroma, STROMA, UNREACTIVESTROMA, reactiveStromaProportion);

        // initialise cancer cell
        model.initialCancerEvent(x, y, nV, vesselSetUp, QUIESCENTCANCER, cancerIMT, proliferationCancerSite);

        // draw model
        model.drawModelLayers(agentLayer, drugLayer, proliferationLayer);

        // store frame in gif
//        growCancerGif.AddFrame(agentLayer);
        growCancerProliferationGif.AddFrame(proliferationLayer);

//        agentLayer.ToPNG("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\growTumourInit.png");

        // declare time counter
        int j = 0;
        int extendedTime = -1;

        // main loop
        while (cancerDetectedEXTENDED == false) {
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
            if (cancerEliminated == true) {
                break;
            }

            // check if cancer cell population has reached detection threshold
            double cancerPop = model.populationType(CANCER);
            double quiescentCancerPop = model.populationType(QUIESCENTCANCER);
            double totalCancerPop = cancerPop + quiescentCancerPop;
            if (cancerDetected == false) {
                if (totalCancerPop >= cancerDetectionThreshold) {
                    cancerDetected = true;
                    extendedTime = extendedTimeLength;

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

                    agentLayer.ToPNG("C:\\Users\\amymm\\Documents\\PHD Data\\Data\\growTumourAgentsFin.png");
                    //drugLayer.ToPNG("C:\\Users\\amymm\\Documents\\PHD Data\\Data\\growTumourDrugFin.png");
                    proliferationLayer.ToPNG("C:\\Users\\amymm\\Documents\\PHD Data\\Data\\growTumourProliferationFin.png");

                    System.out.println("cancer detected iteration: " + j);
                }
            }
            if (cancerDetected == true) {
                int imageNumber = extendedTimeLength - extendedTime;
                agentLayer.ToPNG("C:\\Users\\amymm\\Documents\\PHD Data\\Data\\growTumourAgentsFin_" + imageNumber + ".png");
                //drugLayer.ToPNG("C:\\Users\\amymm\\Documents\\PHD Data\\Data\\growTumourDrugFin" + imageNumber + ".png");
                proliferationLayer.ToPNG("C:\\Users\\amymm\\Documents\\PHD Data\\Data\\growTumourProliferationFin" + imageNumber + ".png");

                extendedTime--;
            }
            if (cancerDetectedEXTENDED == false) {
                if (extendedTime == 0) {
                    cancerDetectedEXTENDED = true;
                    //System.out.println("cancer detected iteration: " + j);
                }
            }



            // increment timestep
            model.IncTick();

            // draw
            model.drawModelLayers(agentLayer, drugLayer, proliferationLayer);

            // store frame in gif
//            growCancerGif.AddFrame(agentLayer);
            growCancerProliferationGif.AddFrame(proliferationLayer);

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

        growCancerProliferationGif.Close();

        // store stroma and cancer location data
        int stromaInitCnEXTENDED = (int)model.populationType(STROMA);
        int unreactiveStromaInitCnEXTENDED = (int)model.populationType(UNREACTIVESTROMA);
        int totalStromaInitCnEXTENDED = stromaInitCnEXTENDED + unreactiveStromaInitCnEXTENDED;
        countInitialStromaEXTENDED.Write(Integer.toString(totalStromaInitCnEXTENDED));
        countInitialStromaEXTENDED.Close();
        model.stromaLocations(STROMA, UNREACTIVESTROMA, locationsStromaEXTENDED);
        locationsStromaEXTENDED.Close();
        int cancerInitCnEXTENDED = (int)model.populationType(CANCER);
        int quiescentCancerInitCnEXTENDED = (int)model.populationType(QUIESCENTCANCER);
        int totalCancerInitCnEXTENDED = cancerInitCnEXTENDED + quiescentCancerInitCnEXTENDED;
        countInitialCancerEXTENDED.Write(Integer.toString(totalCancerInitCnEXTENDED));
        countInitialCancerEXTENDED.Close();
        model.cancerLocations(CANCER, QUIESCENTCANCER, locationsCancerEXTENDED);
        locationsCancerEXTENDED.Close();
        model.finalProliferationSignal(proliferationSignalValuesEXTENDED);
        proliferationSignalValuesEXTENDED.Close();

        agentLayer.ToPNG("C:\\Users\\amymm\\Documents\\PHD Data\\Data\\growTumourAgentsFinEXTENDED.png");
        //drugLayer.ToPNG("C:\\Users\\amymm\\Documents\\PHD Data\\Data\\growTumourDrugFinEXTENDED.png");
        proliferationLayer.ToPNG("C:\\Users\\amymm\\Documents\\PHD Data\\Data\\growTumourProliferationFinEXTENDED.png");

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
