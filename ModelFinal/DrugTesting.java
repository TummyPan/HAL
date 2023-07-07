package ModelFinal;

import HAL.Gui.GridWindow;
import HAL.Tools.FileIO;

import java.io.IOException;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;

public class DrugTesting extends ExampleCell {
    public static void main(String[]args) throws IOException {

        // log start time of algorithm
        long startTime = System.nanoTime();

        // log algorithm date and time of run
        LocalDateTime rawDateTime = LocalDateTime.now();
        DateTimeFormatter formatDateTime = DateTimeFormatter.ofPattern("dd-MM-yyyy_HH-mm-ss");
        String runDateTime = rawDateTime.format(formatDateTime);

        // name results folder
        String inputFolder = "C:\\MyFiles\\PHD Data\\";

        // name results folder
        String resultsFolder = "C:\\MyFiles\\PHD Data\\DrugMF\\";

        //double [] drugRemovalRateVesselDaysVary = {5, 10, 15}; // (days) removal rate of drugA at vessel site

        for (int k = 0; k < 1; k++) {
            System.out.println("Run " + k);
            for (int i = 0; i < 1; i++) {

                // print simulation number
                System.out.println("    Simulation " + i);

                // number of iterations
                int timesteps = 30000; // number of iterations of model
                double noDays = (int) Math.round(timesteps / timestepPerDay);
                System.out.println("Number of days: " + noDays);
/*                System.out.println("stroma: " + STROMA);
                System.out.println("activated stroma: " + ACTIVATEDSTROMA);
                System.out.println("unreactive stroma: " + UNREACTIVESTROMA);
                System.out.println("cancer: " + CANCER);
                System.out.println("quiescent cancer: " + QUIESCENTCANCER);
                System.out.println("vessel: " + VESSEL);
*/
                // declare time and population storage arrays
                double[] timeVector = new double[timesteps];
                //double[] stromaPopulation = new double[timesteps];
                //double[] cancerPopulation = new double[timesteps];
                //double[] activatedStromaPopulation = new double[timesteps];
                //double[] contributingActivatedStromaPopulation = new double[timesteps];
                //double[] quiescentCancerPopulation = new double[timesteps];
                double[] meanFieldDrug = new double[timesteps];
                //ArrayList activationDelayTimes = new ArrayList<Double>();
                //ArrayList deactivationDelayTimes = new ArrayList<Double>();

                // declare time and population storage arrays Small
                //double[] stromaPopulationSmall = new double[timesteps];
                //double[] cancerPopulationSmall = new double[timesteps];
                //double[] activatedStromaPopulationSmall = new double[timesteps];
                //double[] contributingActivatedStromaPopulationSmall = new double[timesteps];
                //double[] quiescentCancerPopulationSmall = new double[timesteps];

                // set up animation
                GridWindow agentLayer = new GridWindow(x, y, 1);
                GridWindow drugLayer = new GridWindow(x, y, 1);
                GridWindow proliferationLayer = new GridWindow(x, y, 1);
                //GridWindow agentLayerSmall = new GridWindow(xSmall, ySmall, 3);
                //GridWindow drugLayerSmall = new GridWindow(xSmall, ySmall, 3);
                //GridWindow proliferationLayerSmall = new GridWindow(xSmall, ySmall, 3);
                ExampleGrid model = new ExampleGrid(x, y);
                // declare animation as gif
/*                GifMaker agentGif = new GifMaker(resultsFolder + "\\Movies\\Agents_" + runDateTime + "_run_" + k + "_sim_" + i + ".gif", 10, false);
                GifMaker drugGif = new GifMaker(resultsFolder + "\\Movies\\Drug_" + runDateTime + "_run_" + k + "_sim_" + i + ".gif", 10, false);
                GifMaker proliferationGif = new GifMaker(resultsFolder + "\\Movies\\Proliferation_" + runDateTime + "_run_" + k + "_sim_" + i + ".gif", 10, false);
                GifMaker agentSmallGif = new GifMaker(resultsFolder + "\\Movies\\AgentsSmall_" + runDateTime + "_run_" + k + "_sim_" + i + ".gif", 10, false);
                GifMaker drugSmallGif = new GifMaker(resultsFolder + "\\Movies\\DrugSmall_" + runDateTime + "_run_" + k + "_sim_" + i + ".gif", 10, false);
                GifMaker proliferationSmallGif = new GifMaker(resultsFolder + "\\Movies\\ProliferationSmall_" + runDateTime + "_run_" + k + "_sim_" + i + ".gif", 10, false);
*/
                // declare parameter recorder
                FileIO parameters = new FileIO(resultsFolder + "\\Data\\Parameter Values_" + runDateTime + "_run_" + k + "_sim_" + i + ".csv", "w");

                // declare storage
                //FileIO populations = new FileIO(resultsFolder + "\\Data\\Populations_" + runDateTime + "_run_" + k + "_sim_" + i + ".csv", "w");
                FileIO drugMeanFieldValues = new FileIO(resultsFolder + "\\Data\\DrugMF_" + runDateTime + "_run_" + k + "_sim_" + i + ".csv", "w");
                //FileIO populationsSmall = new FileIO(resultsFolder + "\\Data\\PopulationsSmall_" + runDateTime + "_run_" + k + "_sim_" + i + ".csv", "w");
                //FileIO DelayTimesActivation = new FileIO(resultsFolder + "\\Data\\DelayTimesActivation_" + runDateTime + "_run_" + k + "_sim_" + i + ".csv", "w");
                //FileIO DelayTimesDeactivation = new FileIO(resultsFolder + "\\Data\\DelayTimesDeactivation_" + runDateTime + "_run_" + k + "_sim_" + i + ".csv", "w");

                //double drugRemovalRateVesselTimestep = drugRemovalRateVesselDaysVary[i] / timestepPerDay;

                model.getParameters(parameters, x, y, xSmall, ySmall, xDomain, yDomain, timestepPerDay, STROMA, VESSEL, CANCER, ACTIVATEDSTROMA, QUIESCENTCANCER, UNREACTIVESTROMA, stromaNCI, stromaDieProb, stromaIMTDays
                        , stromaDrugThresholdReactivity, stromaActivationNormalisationProbability, reactiveStromaProportion
                        , cancerIMT, drugDiffusionCoefficientDays, drugRemovalRateVesselDays, drugConcentrationVessel
                        , proliferationDiffusionCoefficientDays, autocrineProliferationSignalProductionDays
                        , paracrineProliferationSignalProductionDays, proliferationDegradationDrugDays
                        , proliferationCancerSite, activationDistMeanDays, activationDistStdDevDays, deactivationDistMeanDays, deactivationDistStdDevDays, initialTreatmentDays, treatmentDaysDrug, holidayDaysDrug);


                // read vessel initialisation
                int[][] vesselSetUp = model.readVesselInitialisation(inputFolder + "\\VesselLocations\\vesselLocations.csv", nV);

                // initialise vessel cells
                model.initialiseVesselCells(vesselSetUp, nV, VESSEL);

                // determine position and age of initial stroma population
                //int initialStromaCount = model.countInitialStroma(inputFolder + "\\GrowCancer\\stromaInitialCount.csv");
                //double[][] initialStroma = model.initialStromaSetUp(inputFolder + "\\GrowCancer\\locationsStroma.csv", initialStromaCount);
                //double[][] initialStroma = model.initialPassiveReactiveStromaSetUp(inputFolder + "\\GrowCancer\\locationsPassiveReactiveStroma.csv", initialStromaCount);

                // determine position and age of initial cancer population
                //int initialCancerCount = model.countInitialCancer(inputFolder + "\\GrowCancer\\cancerInitialCount.csv");
                //double[][] initialCancer = model.initialCancerSetUp(inputFolder + "\\GrowCancer\\locationsCancer.csv", initialCancerCount);

                // determine intial proliferation signal over domain
                //double[] proliferationSignalInit = model.initialProliferationSignal(inputFolder + "\\GrowCancer\\valuesProliferationSignal.csv", x, y);

                // initialise stroma cells
                //model.initialisePassiveReactiveStroma(initialStroma, STROMA, UNREACTIVESTROMA, reactiveStromaProportion);

                // initialise proliferation signal
                //model.initialiseProliferationSignal(proliferationSignalInit);

                // initialise cancer
                //model.initialiseCancer(initialCancer, CANCER, QUIESCENTCANCER, cellIDCount, vesselSetUp, clusterBoxSize, VESSEL, x, y);

                model.drawModelLayers(agentLayer, drugLayer, proliferationLayer);
                //model.drawModelLayersSmall(x, y, xSmall, ySmall, agentLayerSmall, drugLayerSmall, proliferationLayerSmall);

/*                agentLayer.ToPNG(resultsFolder + "\\Images\\Agents_Init_" + runDateTime + "_run_" + k + "_sim_" + i + ".png");
                drugLayer.ToPNG(resultsFolder + "\\Images\\Drug_Init_" + runDateTime + "_run_" + k + "_sim_" + i + ".png");
                proliferationLayer.ToPNG(resultsFolder + "\\Images\\Proliferation_Init_" + runDateTime + "_run_" + k + "_sim_" + i + ".png");
                agentLayerSmall.ToPNG(resultsFolder + "\\Images\\AgentsSmall_Init_" + runDateTime + "_run_" + k + "_sim_" + i + ".png");
                drugLayerSmall.ToPNG(resultsFolder + "\\Images\\DrugSmall_Init_" + runDateTime + "_run_" + k + "_sim_" + i + ".png");
                proliferationLayerSmall.ToPNG(resultsFolder + "\\Images\\ProliferationSmall_Init_" + runDateTime + "_run_" + k + "_sim_" + i + ".png");
*/
                // main loop
                for (int j = 0; j < timesteps; j++) {
                    // time between animation frames
                    agentLayer.TickPause(1);
                    drugLayer.TickPause(1);
                    proliferationLayer.TickPause(1);
                    //agentLayerSmall.TickPause(1);
                    //drugLayerSmall.TickPause(1);
                    //proliferationLayerSmall.TickPause(1);

                    // create time counter
                    int timeCounter = j;


                    // collect population data
                    double timeCount = j * 1 / timestepPerDay;
                    //double stromaCn = model.populationType(STROMA);
                    //double cancerCn = model.populationType(CANCER);
                    //double activatedStromaCn = model.populationType(ACTIVATEDSTROMA);
                    //double contributingActivatedStromaCn = model.contributingActivatedStromaPopulationType(STROMA, ACTIVATEDSTROMA);
                    //double quiescentCancerCn = model.populationType(QUIESCENTCANCER);
                    double drugMF = model.drugMeanField();
                    //double stromaSmallCn = model.populationTypeSmall(STROMA);
                    //double cancerSmallCn = model.populationTypeSmall(CANCER);
                    //double activatedStromaSmallCn = model.populationTypeSmall(ACTIVATEDSTROMA);
                    //double contributingActivatedStromaSmallCn = model.contributingActivatedStromaPopulationTypeSmall(STROMA, ACTIVATEDSTROMA);
                    //double quiescentCancerSmallCn = model.populationTypeSmall(QUIESCENTCANCER);
                    timeVector[j] = timeCount;
                    meanFieldDrug[j] = drugMF;
                    //stromaPopulation[j] = stromaCn;
                    //cancerPopulation[j] = cancerCn;
                    //activatedStromaPopulation[j] = activatedStromaCn;
                    //contributingActivatedStromaPopulation[j] = contributingActivatedStromaCn;
                    //quiescentCancerPopulation[j] = quiescentCancerCn;
                    //stromaPopulationSmall[j] = stromaSmallCn;
                    //cancerPopulationSmall[j] = cancerSmallCn;
                    //activatedStromaPopulationSmall[j] = activatedStromaSmallCn;
                    //contributingActivatedStromaPopulationSmall[j] = contributingActivatedStromaSmallCn;
                    //quiescentCancerPopulationSmall[j] = quiescentCancerSmallCn;

                    // print out message every 500 iterations
                    if (j % 500 == 0) {
                        System.out.println("        iteration: " + j);
                    }

                    // age cells
                    //model.ageCells(STROMA, ACTIVATEDSTROMA, CANCER, UNREACTIVESTROMA);

                    // divide cells
                    //model.divideCells(stromaNCI, STROMA, ACTIVATEDSTROMA, stromaIMT, CANCER, cancerIMT, proliferationCancerSite, UNREACTIVESTROMA, vesselSetUp);

                    int[] treatmentTimestep = model.treatmentScheduleTimesteps(timestepPerDay, treatmentDaysDrug, holidayDaysDrug);

                    // stroma cells dying
                    //model.stromaDieCells(stromaDieProb, STROMA, ACTIVATEDSTROMA, UNREACTIVESTROMA);

                    // cancer cell proliferation update and age cancer cells
                    //model.cancerDieProliferationUpdate(CANCER, QUIESCENTCANCER, cancerEliminated);

                    // stroma cells activated or normalised
                    //model.stromaStatusUpdate(STROMA, ACTIVATEDSTROMA, stromaDrugThresholdReactivity, CANCER, QUIESCENTCANCER, stromaActivationNormalisationProbability, activationDistMean, activationDistStdDev, activationDelayTimes, deactivationDistMean, deactivationDistStdDev, deactivationDelayTimes);

                    // update activated stroma neighbours
                    //model.neighboursActivatedStroma(CANCER, STROMA, ACTIVATEDSTROMA);

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
                    //model.drawModelLayersSmall(x, y, xSmall, ySmall, agentLayerSmall, drugLayerSmall, proliferationLayerSmall);

                    // store frame in gif
/*                    if (j % 10 == 0) {
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
                        //model.recordAgentData(CANCER, QUIESCENTCANCER, STROMA, ACTIVATEDSTROMA, i, j, k, runDateTime);
                        // record each position on domain data
                        //model.positionData(i, j, k, x, y, runDateTime, resultsFolder);
                        //model.recordNeighbours(x, y, i, j, k, STROMA, ACTIVATEDSTROMA, UNREACTIVESTROMA, CANCER, QUIESCENTCANCER, VESSEL, resultsFolder, runDateTime);
                        //model.recordNeighbourInformation(i, j, k, x, y, xSmall, ySmall, runDateTime, CANCER, ACTIVATEDSTROMA, UNREACTIVESTROMA, resultsFolder);
                    }


                    // save image as png
/*                    if (j % 1000 == 0) {
                        agentLayer.ToPNG(resultsFolder + "\\Images\\Agents_" + j + "_" + runDateTime + "_run_" + k + "_sim_" + i + ".png");
                        drugLayer.ToPNG(resultsFolder + "\\Images\\Drug_" + j + "_" + runDateTime + "_run_" + k + "_sim_" + i + ".png");
                        proliferationLayer.ToPNG(resultsFolder + "\\Images\\Proliferation_" + j + "_" + runDateTime + "_run_" + k + "_sim_" + i + ".png");
                        agentLayerSmall.ToPNG(resultsFolder + "\\Images\\AgentsSmall_" + j + "_" + runDateTime + "_run_" + k + "_sim_" + i + ".png");
                        drugLayerSmall.ToPNG(resultsFolder + "\\Images\\DrugSmall_" + j + "_" + runDateTime + "_run_" + k + "_sim_" + i + ".png");
                        proliferationLayerSmall.ToPNG(resultsFolder + "\\Images\\ProliferationSmall_" + j + "_" + runDateTime + "_run_" + k + "_sim_" + i + ".png");
                    }
*/
                    // check for cancer elimination to break loop
                    //cancerEliminated = model.testElimination(CANCER, QUIESCENTCANCER);

                    // check for cancer elimination
                    //if (cancerEliminated == true) {
                    //    break;
                    //}

                }

/*                agentLayer.ToPNG(resultsFolder + "\\Images\\Agents_fin_" + runDateTime + "_run_" + k + "_sim_" + i + ".png");
                drugLayer.ToPNG(resultsFolder + "\\Images\\Drug_fin_" + runDateTime + "_run_" + k + "_sim_" + i + ".png");
                proliferationLayer.ToPNG(resultsFolder + "\\Images\\Proliferation_fin_" + runDateTime + "_run_" + k + "_sim_" + i + ".png");
                agentLayerSmall.ToPNG(resultsFolder + "\\Images\\AgentsSmall_fin_" + runDateTime + "_run_" + k + "_sim_" + i + ".png");
                drugLayerSmall.ToPNG(resultsFolder + "\\Images\\DrugSmall_fin_" + runDateTime + "_run_" + k + "_sim_" + i + ".png");
                proliferationLayerSmall.ToPNG(resultsFolder + "\\Images\\ProliferationSmall_fin_" + runDateTime + "_run_" + k + "_sim_" + i + ".png");
*/
                // end gif storage
/*                agentGif.Close();
                drugGif.Close();
                proliferationGif.Close();
                agentSmallGif.Close();
                drugSmallGif.Close();
                proliferationSmallGif.Close();
*/

                // record population data as csv
                //model.populationData(timeVector, stromaPopulation, cancerPopulation, activatedStromaPopulation, contributingActivatedStromaPopulation, quiescentCancerPopulation, populations);
                model.meanFieldDrug(timeVector, meanFieldDrug, drugMeanFieldValues);
                //model.populationData(timeVector, stromaPopulationSmall, cancerPopulationSmall, activatedStromaPopulationSmall, contributingActivatedStromaPopulationSmall, quiescentCancerPopulationSmall, populationsSmall);
                //model.meanDelayTimeData(activationDelayTimes, DelayTimesActivation);
                //model.meanDelayTimeData(deactivationDelayTimes, DelayTimesDeactivation);
                //populations.Close();
                drugMeanFieldValues.Close();
                //populationsSmall.Close();
                //DelayTimesActivation.Close();
                //DelayTimesDeactivation.Close();
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


