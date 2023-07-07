package ModelFinal;

import HAL.Gui.GifMaker;
import HAL.Gui.GridWindow;
import HAL.Tools.FileIO;

import java.io.IOException;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;

public class LatentActivationDeactivation extends ExampleCell {
    public static void main(String[]args) throws IOException {

        // log start time of algorithm
        long startTime = System.nanoTime();

        // log algorithm date and time of run
        LocalDateTime rawDateTime = LocalDateTime.now();
        DateTimeFormatter formatDateTime = DateTimeFormatter.ofPattern("dd-MM-yyyy_HH-mm-ss");
        String runDateTime = rawDateTime.format(formatDateTime);

        for (int i = 0; i < 1; i++) {

            // print simulation number
            System.out.println("Simulation " + i);

            // number of iterations
            int timesteps = 4000; // number of iterations of model
            double noDays = (int) Math.round(timesteps / timestepPerDay);
            System.out.println("Number of days: " + noDays);

            // declare time and population storage arrays
            double[] timeVector = new double[timesteps];
            double[] stromaPopulation = new double[timesteps];
            double[] cancerPopulation = new double[timesteps];
            double[] activatedStromaPopulation = new double[timesteps];
            double[] contributingActivatedStromaPopulation = new double[timesteps];
            double[] quiescentCancerPopulation = new double[timesteps];
            double[] meanFieldDrug = new double[timesteps];
            ArrayList activationDelayTimes = new ArrayList<Double>();
            ArrayList deactivationDelayTimes = new ArrayList<Double>();

            // declare time and population storage arrays Small
            double[] stromaPopulationSmall = new double[timesteps];
            double[] cancerPopulationSmall = new double[timesteps];
            double[] activatedStromaPopulationSmall = new double[timesteps];
            double[] contributingActivatedStromaPopulationSmall = new double[timesteps];
            double[] quiescentCancerPopulationSmall = new double[timesteps];

            // set up animation
            GridWindow agentLayer = new GridWindow(x, y, 1);
            GridWindow drugLayer = new GridWindow(x, y, 1);
            GridWindow proliferationLayer = new GridWindow(x, y, 1);
            GridWindow agentLayerSmall = new GridWindow(xSmall, ySmall, 3);
            GridWindow drugLayerSmall = new GridWindow(xSmall, ySmall, 3);
            GridWindow proliferationLayerSmall = new GridWindow(xSmall, ySmall, 3);
            ExampleGrid model = new ExampleGrid(x, y);
            // declare animation as gif
            /*GifMaker agentGif = new GifMaker("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\Agents_" + runDateTime + "_sim_" + i + ".gif", 10, false);
            GifMaker drugGif = new GifMaker("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\Drug_" + runDateTime + "_sim_" + i + ".gif", 10, false);
            GifMaker proliferationGif = new GifMaker("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\Proliferation_" + runDateTime + "_sim_" + i + ".gif", 10, false);
            GifMaker agentSmallGif = new GifMaker("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\AgentsSmall_" + runDateTime + "_sim_" + i + ".gif", 10, false);
            GifMaker drugSmallGif = new GifMaker("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\DrugSmall_" + runDateTime + "_sim_" + i + ".gif", 10, false);
            GifMaker proliferationSmallGif = new GifMaker("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\ProliferationSmall_" + runDateTime + "_sim_" + i + ".gif", 10, false);
            */
            // declare parameter recorder
            FileIO parameters = new FileIO("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\Parameter Values_" + runDateTime + "_sim_"+ i + ".csv", "w");

            // declare drug mean field storage
            //FileIO MFDrug = new FileIO("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\Drug Mean Field" + runDateTime + " sim "+ i + ".csv", "w");
            //FileIO populations = new FileIO("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\Populations_" + runDateTime + "_sim_" + i + ".csv", "w");
            //FileIO populationsSmall = new FileIO("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\PopulationsSmall_" + runDateTime + "_sim_" + i + ".csv", "w");
            //FileIO DelayTimesActivation = new FileIO("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\DelayTimesActivation_" + runDateTime + "_sim_" + i + ".csv", "w");
            //FileIO DelayTimesDeactivation = new FileIO("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\DelayTimesDeactivation_" + runDateTime + "_sim_" + i + ".csv", "w");


            model.getParameters(parameters, x, y, xSmall, ySmall, xDomain, yDomain, timestepPerDay, STROMA, VESSEL, CANCER, ACTIVATEDSTROMA, QUIESCENTCANCER, UNREACTIVESTROMA, stromaNCI, stromaDieProb, stromaIMTDays
                    , stromaDrugThresholdReactivity, stromaActivationNormalisationProbability, reactiveStromaProportion
                    , cancerIMT, drugDiffusionCoefficientDays, drugRemovalRateVesselDays, drugConcentrationVessel
                    , proliferationDiffusionCoefficientDays, autocrineProliferationSignalProductionDays
                    , paracrineProliferationSignalProductionDays, proliferationDegradationDrugDays
                    , proliferationCancerSite, activationDistMeanDays, activationDistStdDevDays, deactivationDistMeanDays, deactivationDistStdDevDays, initialTreatmentDays, treatmentDaysDrug, holidayDaysDrug);


            // read vessel initialisation
            int[][] vesselSetUp = model.readVesselInitialisation("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\Results\\vesselLocations" +
                    "\\vesselLocations.csv", nV);

            // initialise vessel cells
            model.initialiseVesselCells(vesselSetUp, nV, VESSEL);

            // determine position and age of initial stroma population
            int initialStromaCount = model.countInitialStroma("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\Results\\GrowCancer" +
                    "\\stromaInitialCount.csv");
            double[][] initialStroma = model.initialStromaSetUp("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\Results\\GrowCancer" +
                    "\\locationsStroma.csv", initialStromaCount);

            // determine position and age of initial cancer population
            int initialCancerCount = model.countInitialCancer("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\Results\\GrowCancer" +
                    "\\cancerInitialCount.csv");
            double[][] initialCancer = model.initialCancerSetUp("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\Results\\GrowCancer" +
                    "\\locationsCancer.csv", initialCancerCount);

            // determine intial proliferation signal over domain
            double[] proliferationSignalInit = model.initialProliferationSignal("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\Results\\GrowCancer" +
                    "\\valuesProliferationSignal.csv", x, y);

            // initialise stroma cells
            model.initialiseStroma(initialStroma, STROMA, UNREACTIVESTROMA, reactiveStromaProportion);

            // initialise proliferation signal
            model.initialiseProliferationSignal(proliferationSignalInit);

            // initialise cancer
            model.initialiseCancer(initialCancer, CANCER, QUIESCENTCANCER, cellIDCount, vesselSetUp, clusterBoxSize, VESSEL, x, y);

            model.drawModelLayers(agentLayer, drugLayer, proliferationLayer);
            model.drawModelLayersSmall(x, y, xSmall, ySmall, agentLayerSmall, drugLayerSmall, proliferationLayerSmall);

            //agentLayer.ToPNG("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\Agents_Init_" + runDateTime + "_sim_" + i + ".png");
            //drugLayer.ToPNG("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\Drug_Init_" + runDateTime + "_sim_" + i + ".png");
            //proliferationLayer.ToPNG("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\Proliferation_Init_" + runDateTime + "_sim_" + i + ".png");
            //agentLayerSmall.ToPNG("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\AgentsSmall_Init_" + runDateTime + "_sim_" + i + ".png");
            //drugLayerSmall.ToPNG("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\DrugSmall_Init_" + runDateTime + "_sim_" + i + ".png");
            //proliferationLayerSmall.ToPNG("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\ProliferationSmall_Init_" + runDateTime + "_sim_" + i + ".png");

            // main loop
            for (int j = 0; j < timesteps; j++) {
                // time between animation frames
                agentLayer.TickPause(1);
                drugLayer.TickPause(1);
                proliferationLayer.TickPause(1);
                agentLayerSmall.TickPause(1);
                drugLayerSmall.TickPause(1);
                proliferationLayerSmall.TickPause(1);

                // create time counter
                int timeCounter = j;


                // collect population data
                /*double timeCount = j * 1 / timestepPerDay;
                double stromaCn = model.populationType(STROMA);
                double cancerCn = model.populationType(CANCER);
                double activatedStromaCn = model.populationType(ACTIVATEDSTROMA);
                double contributingActivatedStromaCn = model.contributingActivatedStromaPopulationType(STROMA, ACTIVATEDSTROMA);
                double quiescentCancerCn = model.populationType(QUIESCENTCANCER);
                //double drugMF = model.drugMeanField();
                double stromaSmallCn = model.populationTypeSmall(STROMA);
                double cancerSmallCn = model.populationTypeSmall(CANCER);
                double activatedStromaSmallCn = model.populationTypeSmall(ACTIVATEDSTROMA);
                double contributingActivatedStromaSmallCn = model.contributingActivatedStromaPopulationTypeSmall(STROMA, ACTIVATEDSTROMA);
                double quiescentCancerSmallCn = model.populationTypeSmall(QUIESCENTCANCER);
                timeVector[j] = timeCount;
                stromaPopulation[j] = stromaCn;
                cancerPopulation[j] = cancerCn;
                activatedStromaPopulation[j] = activatedStromaCn;
                contributingActivatedStromaPopulation[j] = contributingActivatedStromaCn;
                quiescentCancerPopulation[j] = quiescentCancerCn;
                stromaPopulationSmall[j] = stromaSmallCn;
                cancerPopulationSmall[j] = cancerSmallCn;
                activatedStromaPopulationSmall[j] = activatedStromaSmallCn;
                contributingActivatedStromaPopulationSmall[j] = contributingActivatedStromaSmallCn;
                quiescentCancerPopulationSmall[j] = quiescentCancerSmallCn;
                //meanFieldDrug[j] = drugMF;
*/
                // print out message every 500 iterations
                if (j % 500 == 0) {
                    System.out.println("iteration: " + j);
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

                // stroma cells activated or normalised
                model.stromaStatusUpdate(STROMA, ACTIVATEDSTROMA, stromaDrugThresholdReactivity, CANCER, QUIESCENTCANCER, stromaActivationNormalisationProbability, activationDistMean, activationDistStdDev, activationDelayTimes, deactivationDistMean, deactivationDistStdDev, deactivationDelayTimes);

                // update activated stroma neighbours
                model.neighboursActivatedStroma(CANCER, STROMA, ACTIVATEDSTROMA);

                // update PDEGrid values
                model.updatePDEValues(x, y, timeCounter, treatmentTimestep[0], treatmentTimestep[1], VESSEL,
                        drugDiffusionCoefficientTimestep, drugConcentrationVessel, drugRemovalRateVesselTimestep,
                        proliferationDiffusionCoefficientTimestep, CANCER, QUIESCENTCANCER, autocrineProliferationSignalProductionTimestep,
                        STROMA, ACTIVATEDSTROMA, paracrineProliferationSignalProductionTimestep,
                        proliferationDegradationDrugTimestep, initialTreatmentDays, timestepPerDay);

                // increment timestep
                model.IncTick();

                // draw
                model.drawModelLayers(agentLayer, drugLayer, proliferationLayer);
                model.drawModelLayersSmall(x, y, xSmall, ySmall, agentLayerSmall, drugLayerSmall, proliferationLayerSmall);

                // store frame in gif
  /*              if (j % 10 == 0) {
                    agentGif.AddFrame(agentLayer);
                    drugGif.AddFrame(drugLayer);
                    proliferationGif.AddFrame(proliferationLayer);
                    agentSmallGif.AddFrame(agentLayerSmall);
                    drugSmallGif.AddFrame(drugLayerSmall);
                    proliferationSmallGif.AddFrame(proliferationLayerSmall);
                }
*/
                // store data
                if (j % 100 == 0) {
                    // record timestep data
                    //model.recordAgentData(CANCER, QUIESCENTCANCER, STROMA, ACTIVATEDSTROMA, i, j, runDateTime);
                    // record each position on domain data
                    //model.positionData(i, j, x, y, runDateTime);
                }


                // save image as png
                /*if (j % 1000 == 0) {
                    agentLayer.ToPNG("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\Agents_" + j + "_" + runDateTime + "_sim_"  + i + ".png");
                    drugLayer.ToPNG("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\Drug_" + j + "_" + runDateTime + "_sim_"  + i + ".png");
                    proliferationLayer.ToPNG("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\Proliferation_" + j + "_" + runDateTime + "_sim_"  + i + ".png");
                    agentLayerSmall.ToPNG("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\AgentsSmall_" + j + "_" + runDateTime + "_sim_"  + i + ".png");
                    drugLayerSmall.ToPNG("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\DrugSmall_" + j + "_" + runDateTime + "_sim_"  + i + ".png");
                    proliferationLayerSmall.ToPNG("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\ProliferationSmall_" + j + "_" + runDateTime + "_sim_"  + i + ".png");
                }*/

                // check for cancer elimination to break loop
                cancerEliminated = model.testElimination(CANCER, QUIESCENTCANCER);

                // check for cancer elimination
                if (cancerEliminated == true) {
                    break;
                }

            }

            agentLayer.ToPNG("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\Agents_fin_" + runDateTime + "_sim_"  + i + ".png");
            drugLayer.ToPNG("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\Drug_fin_" + runDateTime + "_sim_"  + i + ".png");
            proliferationLayer.ToPNG("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\Proliferation_fin_" + runDateTime + "_sim_"  + i + ".png");
            agentLayerSmall.ToPNG("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\AgentsSmall_fin_" + runDateTime + "_sim_"  + i + ".png");
            drugLayerSmall.ToPNG("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\DrugSmall_fin_" + runDateTime + "_sim_"  + i + ".png");
            proliferationLayerSmall.ToPNG("C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\ProliferationSmall_fin_" + runDateTime + "_sim_"  + i + ".png");

            // end gif storage
            /*agentGif.Close();
            drugGif.Close();
            proliferationGif.Close();
            agentSmallGif.Close();
            drugSmallGif.Close();
            proliferationSmallGif.Close();
*/

            // record drug mean field as csv
            //model.meanFieldDrug(meanFieldDrug, MFDrug);
            //MFDrug.Close();

            // record population data as csv
            //model.populationData(timeVector, stromaPopulation, cancerPopulation, activatedStromaPopulation, contributingActivatedStromaPopulation, quiescentCancerPopulation, populations);
            //model.populationData(timeVector, stromaPopulationSmall, cancerPopulationSmall, activatedStromaPopulationSmall, contributingActivatedStromaPopulationSmall, quiescentCancerPopulationSmall, populationsSmall);
            //model.meanDelayTimeData(activationDelayTimes, DelayTimesActivation);
            //model.meanDelayTimeData(deactivationDelayTimes, DelayTimesDeactivation);
            //populations.Close();
            //populationsSmall.Close();
            //DelayTimesActivation.Close();
            //DelayTimesDeactivation.Close();
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


