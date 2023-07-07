package ModelFinal;

import HAL.GridsAndAgents.AgentGrid2D;
import HAL.GridsAndAgents.AgentSQ2Dunstackable;
import HAL.GridsAndAgents.PDEGrid2D;
import HAL.Gui.GridWindow;
import HAL.Rand;
import HAL.Tools.FileIO;
import HAL.Tools.Internal.Gaussian;
import HAL.Util;
import org.knowm.xchart.BitmapEncoder;
import org.knowm.xchart.SwingWrapper;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYSeries;
import org.knowm.xchart.style.Styler;
import org.knowm.xchart.style.markers.SeriesMarkers;
import org.lwjgl.Sys;
import org.w3c.dom.ls.LSOutput;

import java.awt.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.StringTokenizer;
import java.time.LocalDateTime; // import the LocalDateTime class

import static HAL.Util.*;

class ExampleCell extends AgentSQ2Dunstackable<ExampleGrid> {
    // declare individual cell attributes
    public int cellID;
    public int type;
    public double IMTAge;
    public double cancerVariedIMT;
    public double hd;
    public double hp;
    public double activationProbability;
    public double distanceNearestVessel;
    public double vesselClusterValue;
    public int activatedStromaNeighbours;
    public double activationDelay;
    public double deactivationDelay;

    // declare parameters

    // grid parameters
    static int x = 300; // number of cells horizontal

    static int y = 300; // number of cells vertical
    static int xSmall = 100;
    static int ySmall = 100;
    static double xDomain = 0.3; // (cm) horizontal length of sample space
    static double yDomain = 0.3; // (cm) vertical length of sample space
    static double deltaX = xDomain / x; // (cm) length of horizontal gridpoint
    static double deltaY = yDomain / y; // (cm) length of vertical gridpoint
    static double spaceConversion = (1 / deltaX) * (1 / deltaY); // space conversion parameter (used for PDEGrids)

    // time parameters
    static double timestepHour = 0.576;//0.6; // (hours) length of timestep
    // number of timesteps in a day (each timestep is 0.6 hours that is 40 timesteps in a day!)
    static double timestepPerDay = 24 / timestepHour; // number of timesteps in day

    // colour sets
    // tol
    Color D1 = new Color(51, 34, 136, 255); // set colour for cancer population
    Color D2 = new Color(17,119,51, 255);
    Color D3 = new Color(68, 170, 153, 255);
    Color D4 = new Color(136,204,238, 255); // set colour for quiescent cancer population
    Color D5 = new Color(221, 204, 119, 255); // set colour for stroma population
    Color D6 = new Color(204, 102, 119, 255); // set colour for activated stroma population
    Color D7 = new Color(170,68,153, 255);
    Color D8 = new Color(136, 34, 85, 255); // set colour for vessel population


    // Noemis scheme
    Color E1 = new Color(246, 189, 76, 255); // set colour for stroma population
    Color E2 = new Color(0,60,115, 255); // set colour for activated stroma population
    Color E3 = new Color(204, 204, 204, 255); // set colour for inactive stroma population
    Color E4 = new Color(203,0,93, 255); // set colour for cancer population
    Color E5 = new Color(255, 102, 102, 255); // set colour for quiescent cancer population
    Color E6 = new Color(93, 14, 40, 255); // set colour for vessel population
    Color E7 = new Color(255,255,255, 255); // set colour for empty


    // agent paramenters Noemis colours
/*    static int STROMA = RGB256(246, 189, 76); // D5 set colour for stroma cells
    static int VESSEL = RGB256(93, 14, 40); // D8 set colour for vessel cells
    static int CANCER = RGB256(203,0,93); // D7 set colour for cancer cells
    static int ACTIVATEDSTROMA = RGB256(0, 60, 115); // D2 set colour for activated stroma cells
    static int QUIESCENTCANCER = RGB256(255, 102, 102); // D6 set colour for quiescent cancer cells
    static int UNREACTIVESTROMA = RGB256(204, 204, 204); // set colour for unreactive stroma
*/

    // agent paramenters Noemis colours
    static int STROMA = -1; // set code for stroma cells
    static int VESSEL = -2; // set code for vessel cells
    static int CANCER = -3; // set code for cancer cells
    static int ACTIVATEDSTROMA = -4; // set code for activated stroma cells
    static int QUIESCENTCANCER = -5; // set code for quiescent cancer cells
    static int UNREACTIVESTROMA = -6; // set code for unreactive stroma


    // agent paramenters my colours
/*    static int STROMA = RGB256(221, 204, 119); // D5 set colour for stroma cells
    static int VESSEL = RGB256(136, 34, 85); // D8 set colour for vessel cells
    static int CANCER = RGB256(170,68,153); // D7 set colour for cancer cells
    static int ACTIVATEDSTROMA = RGB256(17, 119, 51); // D2 set colour for activated stroma cells
    static int QUIESCENTCANCER = RGB256(204, 102, 119); // D6 set colour for quiescent cancer cells
    static int UNREACTIVESTROMA = RGB256(204, 204, 204); // set colour for unreactive stroma
*/
    // vessel parameters
    static double sigmaMean = 0.016; // (cm) average distance between vessels
    static double sigmaMinRaw = 0.008; // (cm) minimum distance between vessels
    static double sigmaMin = sigmaMinRaw / deltaX; // gridpoint minimum distance between vessels
    static int nV = (int) ((xDomain * yDomain) / Math.pow(sigmaMean, 2)); // number of vessels in tissue
    static int nVU = 324;
    static int clusterBoxSize = 3; // size of box to get vessel density measure


    // stroma parameters
    static int stromaNCI = 6; // stroma contact inhibition (number of empty neighbour cells required for division)
    static double stromaDieProb = 0.0001; // stroma turnover probability
    static int stromaIMTDays = 1; // (days) stroma cell inter mitotic time
    static int stromaIMT = stromaIMTDays * (int)timestepPerDay; // (timesteps) stroma cell inter mitotic time
    static double stromaDrugThresholdReactivity = 0.93;//0.9; // stroma drug threshold reactivity
    static double stromaActivationNormalisationProbability = 0.001; // 10 x stromaDieProb;// stroma activation normalisation probability
    static double reactiveStromaProportion = 0.5; //proportion of stroma cell that can be activated

    static double activationDistMeanDays = 0 / 24; // (days) mean for activation distribution
    static double activationDistStdDevDays = 0 / 24; // (days) standard deviation for activation distribution
    static double activationDistMean = activationDistMeanDays * (int)timestepPerDay; // (timesteps) mean for activation distribution
    static double activationDistStdDev = activationDistStdDevDays * (int)timestepPerDay; // (timesteps) standard deviation for activation distribution
    static double deactivationDistMeanDays = 0 / 24; // (days) mean for activation distribution
    static double deactivationDistStdDevDays = 0 / 24; // (days) standard deviation for activation distribution
    static double deactivationDistMean = deactivationDistMeanDays * (int)timestepPerDay; // (timesteps) mean for activation distribution
    static double deactivationDistStdDev = deactivationDistStdDevDays * (int)timestepPerDay; // (timesteps) standard deviation for activation distribution

    // cancer parameters
    static int cellIDCount = 0;
    static double cancerDetectionThreshold = 1E4; // (Picco) number of cancer cells required for detection
    static boolean cancerDetected = false; // cancer has reached detectable size
    static boolean cancerDetectedEXTENDED = false; // cancer has reached detectable size

    // cancer proliferation signal threshold proliferation
    static int cancerIMTDays = 1; // (days) cancer cell inter mitotic time
    static int cancerIMT = cancerIMTDays *(int)timestepPerDay; // (timesteps) cancer cell inter mitotic time
    static boolean cancerEliminated = false; // boolean for cancer eliminated

    // cancer resistance parameters Cmax = 2 and Cmin = 0.5
    static double DcancerResistanceCost = 0.996099; // value determined by Mathematica - CostParameters or Python - CostHeatMap
    static double AcancerResistanceCost = 0.4892; // value determined by Mathematica - CostParameters or Python - CostHeatMap
    static double BcancerResistanceCost = 0.645299; // value determined by Mathematica - CostParameters or Python - CostHeatMap

    // drug parameters
    static double drugDiffusionCoefficientDays = 1E-5; // (cm^2 per day) diffusion coefficient drug
    // (cm^2 per day) diffusion coefficient proliferation signal
    static double drugDiffusionCoefficientTimestep = spaceConversion * drugDiffusionCoefficientDays / timestepPerDay;
    // (cm^2 per timestep) diffusion coefficient proliferation signal
    static double drugRemovalRateVesselDays = 10; // (days) removal rate of drugA at vessel site
    static double drugRemovalRateVesselTimestep = drugRemovalRateVesselDays / timestepPerDay;
    static double drugConcentrationVessel = 1; // drug concentration on delivery

    // proliferation signal parameters
    static double proliferationDiffusionCoefficientDays = 1E-7;
    // (cell space per timestep) diffusion coefficient drugA
    static double proliferationDiffusionCoefficientTimestep = proliferationDiffusionCoefficientDays / timestepPerDay;
    static double autocrineProliferationSignalProductionDays = 0.65;//0.77;//0.7;//0.72;//0.75;//0.8;//0.65; // (days) autocrine proliferation signal rate
    static double paracrineProliferationSignalProductionDays = 1.54;//3.15;//1.54; // (days) paracrine proliferation signal rate
    static double proliferationDegradationDrugDays = 0.79; // (days) proliferation signal degradation rate
    // (timesteps) removal rate of drugA at vessel site
    static double autocrineProliferationSignalProductionTimestep = autocrineProliferationSignalProductionDays
            / timestepPerDay; // (timesteps) autocrine proliferation signal rate
    static double paracrineProliferationSignalProductionTimestep = paracrineProliferationSignalProductionDays
            / timestepPerDay; // (timesteps) paracrine proliferation signal rate
    // (timesteps) proliferation signal degradation rate
    static double proliferationDegradationDrugTimestep = proliferationDegradationDrugDays / timestepPerDay;
    static double proliferationCancerSite = 0.15; // proliferation signal at initial cancer cell site


    // treatment parameters
    static int treatmentDaysDrug = 40; // number of treatment days for drugA
    static int holidayDaysDrug = 10; // number of holiday days for drugA
    static int initialTreatmentDays = 0; // number of days for initial treatment to get drug to max concentration


    public void stromaDieCell(double stromaDieProb) {
        // turnover of stroma cells with given probability
        if (G.rng.Double() < stromaDieProb) {
            //cell will die
            Dispose();
            return;
        }
    }

    public void stromaDivide(int position, int coordI) {
        // dividing cell dies and new cell of same type is created in same position
        int cellType = this.type;
        Dispose();
        ExampleCell resetCell = G.NewAgentSQ(coordI);
        resetCell.type = cellType;
        resetCell.IMTAge = 0;
        // create daughter cell and assign attributes
        ExampleCell daughterCell = G.NewAgentSQ(position);
        daughterCell.type = cellType;
        daughterCell.IMTAge = 0;
    }

    public void CPCExampleCellDivide(int position) {
        // dividing cell dies and new cell of same type is created in same position
        int cellType = this.type;
        // create daughter cell and assign attributes
        ExampleCell daughterCell = G.NewAgentSQ(position);
        daughterCell.type = cellType;
    }

    public void divideStromaCell(int stromaNCI) {
        // get x and y coordinates of cell
        int coordI = this.Isq();
        int options = MapEmptyHood(G.divHood);
        if (options >= stromaNCI) {
            int position = G.divHood[G.rng.Int(options)];
            //System.out.println("position: " + position);
            stromaDivide(position, coordI);
        }
    }

    public void divideCPCExampleCell() {
        // get x and y coordinates of cell
        int options = MapEmptyHood(G.CPCExampleDivHood);
        if (options != 0) {
            int position = G.CPCExampleDivHood[G.rng.Int(options)];
            CPCExampleCellDivide(position);
        }
    }



    public boolean seedCancer(int QUIESCENTCANCER, boolean cancerSeeded, double cancerIMT, double proliferationCancerSite) {
        // check neighbourhood of vessel cell for empty positions
        int vesselEmptyNeighbours = MapEmptyHood(G.divHood);
        // seeds cancer cell in random empty neighbour of vessel
        if (vesselEmptyNeighbours > 0) {
            int position = G.divHood[G.rng.Int(vesselEmptyNeighbours)];
            ExampleCell cancerSeed = G.NewAgentSQ(position);
            cancerSeed.type = QUIESCENTCANCER;
            cancerSeed.IMTAge = 0;
            // thresholds for death and proliferation
            cancerSeed.hd = 0.1;
            cancerSeed.hp = 0.8;
            // random variation of each cancer cells IMT
            double cancerSeedCancerResistanceCost = DcancerResistanceCost / ((cancerSeed.hd + AcancerResistanceCost) * (cancerSeed.hp + BcancerResistanceCost));
            cancerSeed.cancerVariedIMT = cancerIMT * cancerSeedCancerResistanceCost;
            G.proliferation.Add(cancerSeed.Xsq(), cancerSeed.Ysq(), proliferationCancerSite);
            G.proliferation.Update();
            cancerSeeded = true;
        }
        return cancerSeeded;
    }

    public boolean checkCancerNeighbour(int CANCER, int QUIESCENTCANCER) {
        boolean cancerPresent = false;
        int [] neighbourhood = G.divHood;
        int hoodsize = MapOccupiedHood(neighbourhood);
        for (int i = 0; i < hoodsize; i++) {
            ExampleCell checkCell = G.GetAgent(neighbourhood[i]);
            if (checkCell.type == CANCER || checkCell.type == QUIESCENTCANCER) {
                cancerPresent = true;
            }
        }
        return cancerPresent;
    }

    public int checkActivatedStromaNeighbour(int STROMA, int ACTIVATEDSTROMA) {
        int activatedStromaCn = 0;
        int [] neighbourhood = G.divHood;
        int hoodsize = MapOccupiedHood(neighbourhood);
        for (int i = 0; i < hoodsize; i++) {
            ExampleCell checkCell = G.GetAgent(neighbourhood[i]);
            if ((checkCell.type == ACTIVATEDSTROMA && checkCell.activationDelay < 0) || (checkCell.type == STROMA && checkCell.deactivationDelay > 0)) {
                activatedStromaCn++;
            }
        }
        return activatedStromaCn;
    }

    public int checkCancerNeighbourCount(int CANCER, int QUIESCENTCANCER) {
        int cancerCn = 0;
        int [] neighbourhood = G.divHood;
        int hoodsize = MapOccupiedHood(neighbourhood);
        for (int i = 0; i < hoodsize; i++) {
            ExampleCell checkCell = G.GetAgent(neighbourhood[i]);
            if ((checkCell.type == CANCER) || (checkCell.type == QUIESCENTCANCER)) {
                cancerCn++;
            }
        }
        return cancerCn;
    }

    public int checkVesselNeighbour(int VESSEL) {
        int vesselCn = 0;
        int [] neighbourhood = G.divHood;
        int hoodsize = MapOccupiedHood(neighbourhood);
        for (int i = 0; i < hoodsize; i++) {
            ExampleCell checkCell = G.GetAgent(neighbourhood[i]);
            if (checkCell.type == VESSEL) {
                vesselCn++;
            }
        }
        return vesselCn;
    }

    public int checkStromaNeighbour(int STROMA, int ACTIVATEDSTROMA, int UNREACTIVESTROMA) {
        int stromaCn = 0;
        int [] neighbourhood = G.divHood;
        int hoodsize = MapOccupiedHood(neighbourhood);
        for (int i = 0; i < hoodsize; i++) {
            ExampleCell checkCell = G.GetAgent(neighbourhood[i]);
            if ((checkCell.type == ACTIVATEDSTROMA && checkCell.activationDelay > 0) || (checkCell.type == STROMA && checkCell.deactivationDelay < 0) || (checkCell.type == UNREACTIVESTROMA)) {
                stromaCn++;
            }
        }
        return stromaCn;
    }

    public int checkEmptyNeighbour() {
        int [] neighbourhood = G.divHood;
        int emptyCn = MapEmptyHood(neighbourhood);
        return emptyCn;
    }

    public void cancerDieCell(int CANCER, int QUIESCENTCANCER, boolean cancerEliminated) {
        //cell will die
        Dispose();
        //System.out.println("cancer dies");
        G.checkElimination(CANCER, QUIESCENTCANCER, cancerEliminated);
        return;
    }

    public void cancerDivide(int position, int coordI, int cancerIMT, double proliferationCancerSite, int[][] vesselSetUp, double DcancerResistanceCost, double AcancerResistanceCost, double BcancerResistanceCost) {
        // dividing cell dies and new cell of same type is created in same position
        int cellType = this.type;
        double hdParent = this.hd;
        double hpParent = this.hp;
        Dispose();
        ExampleCell resetCell = G.NewAgentSQ(coordI);
        resetCell.cellID = cellIDCount;
        cellIDCount++;
        resetCell.type = cellType;
        resetCell.IMTAge = 0;

        // no mutation
        resetCell.hd = 0.2;
        resetCell.hp = 0.8;
        resetCell.cancerVariedIMT =  cancerIMT *  + (2*Math.random() - 1)*0.1;

        // mutation
        /*
        double hdProbabilityReset = G.rng.Gaussian(0, 0.25);
        double hpProbabilityReset = G.rng.Gaussian(0, 0.25);
        correctProbability(hdProbabilityReset);
        correctProbability(hpProbabilityReset);
        double hdResetNew = hdParent + 0.1 * Math.round(hdProbabilityReset);
        double hpResetNew = hpParent + 0.1 * Math.round(hpProbabilityReset);
        boolean testResetCellHd = false;
        double [] resetAdj = testProbability(hdResetNew, hpResetNew);
        double hdResetAdj = resetAdj[0];
        double hpResetAdj = resetAdj[1];
        resetCell.hd = hdResetAdj;
        resetCell.hp = hpResetAdj;
        double resetCellCancerResistanceCost = DcancerResistanceCost / ((resetCell.hd + AcancerResistanceCost) * (resetCell.hp + BcancerResistanceCost));
        // random variation of each cancer cells IMT
        resetCell.cancerVariedIMT =  cancerIMT * resetCellCancerResistanceCost;
        */

        // vessel information
        resetCell.distanceNearestVessel = G.nearestVesselDist(vesselSetUp, resetCell.Xsq(), resetCell.Ysq());
        resetCell.vesselClusterValue = G.clusterValue(clusterBoxSize, resetCell.Xsq(), resetCell.Ysq(), VESSEL, x, y);
        // create daughter cell and assign attributes
        ExampleCell daughterCell = G.NewAgentSQ(position);
        daughterCell.cellID = cellIDCount;
        cellIDCount++;
        daughterCell.type = QUIESCENTCANCER;
        daughterCell.IMTAge = 0;

        // no mutation
        daughterCell.hd = 0.2;
        daughterCell.hp = 0.8;
        daughterCell.cancerVariedIMT =  cancerIMT *  + (2*Math.random() - 1)*0.1;

        // mutation
        /*
        double hdProbabilityDaughter = G.rng.Gaussian(0, 0.25);
        double hpProbabilityDaughter = G.rng.Gaussian(0, 0.25);
        correctProbability(hdProbabilityDaughter);
        correctProbability(hpProbabilityDaughter);
        double hdDaughterNew = hdParent + 0.1 * Math.round(hdProbabilityDaughter);
        double hpDaughterNew = hpParent + 0.1 * Math.round(hpProbabilityDaughter);
        boolean testDaughterCellHd = false;
        double [] daughterAdj = testProbability(hdDaughterNew, hpDaughterNew);
        double hdDaughterAdj = daughterAdj[0];
        double hpDaughterAdj = daughterAdj[1];
        daughterCell.hd = hdDaughterAdj;
        daughterCell.hp = hpDaughterAdj;
        double daughterCellCancerResistanceCost = DcancerResistanceCost / ((daughterCell.hd + AcancerResistanceCost) * (daughterCell.hp + BcancerResistanceCost));
        // random variation of each cancer cells IMT
        daughterCell.cancerVariedIMT =  cancerIMT * daughterCellCancerResistanceCost;
        */

        // vessel information
        daughterCell.distanceNearestVessel = G.nearestVesselDist(vesselSetUp, daughterCell.Xsq(), daughterCell.Ysq());
        daughterCell.vesselClusterValue = G.clusterValue(clusterBoxSize, daughterCell.Xsq(), daughterCell.Ysq(), VESSEL, x, y);
        // update proliferation signal PDEGrid
        G.proliferation.Add(daughterCell.Xsq(), daughterCell.Ysq(), proliferationCancerSite);
        G.proliferation.Update();
    }

    public double [] testProbability(double valueHd, double valueHp) {
        double [] returnValues = new double [2];
        if (valueHd >= valueHp) {
            valueHp = valueHd + 0.1;
        }
        if (valueHd < 0.1) {
            valueHd = 0.1;
        } else if (valueHd > 0.8) {
            valueHd = 0.8;
        }
        if (valueHp < 0.2) {
            valueHp = 0.2;
        } else if (valueHp > 0.9) {
            valueHp = 0.9;
        }
        returnValues[0] = valueHd;
        returnValues[1] = valueHp;
        return returnValues;
    }

    public void correctProbability (double probabilityValue) {
        if (probabilityValue < -1.4) {
            probabilityValue = -1;
        }
        else if (probabilityValue > 1.4) {
            probabilityValue = 1;
        }
    }
    public void divideCancerCell(int cancerIMT, double proliferationCancerSite, int[][] vesselSetUp) {
        // declare neighbourhood integers
        int options;
        // get x and y coordinates of cell
        int coordI = this.Isq();
        options = MapEmptyHood(G.divHood);
        if (options > 0) {
            int position = G.divHood[G.rng.Int(options)];
            cancerDivide(position, coordI, cancerIMT, proliferationCancerSite, vesselSetUp, DcancerResistanceCost, AcancerResistanceCost, BcancerResistanceCost);
            //System.out.println("cancer divides");
        }
    }
}

public class ExampleGrid extends AgentGrid2D<ExampleCell> {
    /*public ExampleGrid(int x, int y) {
        super(x, y, ExampleCell.class);
    }*/
    // random number
    Rand rng = new Rand();
    // neighbourhood function
    int [] divHood = Util.MooreHood(false);
    int [] CPCExampleDivHood = Util.VonNeumannHood(false);
    // declare PDE grids
    PDEGrid2D drug;
    PDEGrid2D proliferation;

    public ExampleGrid(int x, int y) {
        // link agents and grid
        super(x, y, ExampleCell.class);
        // create PDE grid for drug
        drug = new PDEGrid2D(x, y);
        // create PDE grid for proliferation signal
        proliferation = new PDEGrid2D(x, y);
    }

    public void printProliferationSignal() {
        int belowThresholdDeath = 0;
        //System.out.println("proliferation.Length: " + proliferation.Length());
        for (int i = 0; i < proliferation.Length(); i++) {
            if ((proliferation.Get(i) > 0) && (proliferation.Get(i) < 0.1)) {
                belowThresholdDeath++;
            }
        }
        System.out.println("cancer to die: " + belowThresholdDeath);
        System.out.println("-----------------");
        System.out.println(Arrays.toString(proliferation.GetField()));
        System.out.println("-----------------");
        //System.out.println(Arrays.toString(proliferation.GetDeltas()));
    }

    public void getParameters(FileIO parameters
            , int x
            , int y
            , int xSmall
            , int ySmall
            , double xDomain
            , double yDomain
            , double timestepPerDay
            , int STROMA
            , int VESSEL
            , int CANCER
            , int ACTIVATEDSTROMA
            , int QUIESCENTSTROMA
            , int UNREACTIVESTROMA
            , int stromaNCI
            , double stromaDieProb
            , int stromaIMTDays
            , double stromaDrugThresholdReactivity
            , double stromaActivationNormalisationProbability
            , double reactiveStromaProportion
            , int cancerIMT
            , double drugDiffusionCoefficientDays
            , double drugRemovalRateVesselDays
            , double drugConcentrationVessel
            , double proliferationDiffusionCoefficientDays
            , double autocrineProliferationSignalProductionDays
            , double paracrineProliferationSignalProductionDays
            , double proliferationDegradationDrugDays
            , double proliferationCancerSite
            , double activationDistMeanDays
            , double activationDistStdDevDays
            , double deactivationDistMeanDays
            , double deactivationDistStdDevDays
            , int initialTreatmentDays
            , int treatmentDaysDrugA
            , int holidayDaysDrugA) {
        parameters.Write("x, " + (double)x + "\n");
        parameters.Write("y, " + (double)y + "\n");
        parameters.Write("xSmall, " + (double)xSmall + "\n");
        parameters.Write("ySmall, " + (double)ySmall + "\n");
        parameters.Write("xDomain, " + xDomain + "\n");
        parameters.Write("yDomain, " + yDomain + "\n");
        parameters.Write("timestepPerDay, " + timestepPerDay + "\n");
        parameters.Write("STROMA, " + STROMA + "\n");
        parameters.Write("VESSEL, " + VESSEL + "\n");
        parameters.Write("CANCER, " + CANCER + "\n");
        parameters.Write("ACTIVATEDSTROMA, " + ACTIVATEDSTROMA + "\n");
        parameters.Write("QUIESCENTCANCER, " + QUIESCENTSTROMA + "\n");
        parameters.Write("UNREACTIVESTROMA, " + UNREACTIVESTROMA + "\n");
        parameters.Write("stromaNCI, " + (double)stromaNCI + "\n");
        parameters.Write("stromaDieProb, " + stromaDieProb + "\n");
        parameters.Write("stromaIMTDays, " + (double)stromaIMTDays + "\n");
        parameters.Write("stromaDrugThresholdReactivity, " + stromaDrugThresholdReactivity + "\n");
        parameters.Write("stromaActivationNormalisationProbability, " + stromaActivationNormalisationProbability + "\n");
        parameters.Write("reactiveStromaProportion, " + reactiveStromaProportion + "\n");
        parameters.Write("cancerIMT, " + (double)cancerIMT + "\n");
        parameters.Write("drugDiffusionCoefficientDays, " + drugDiffusionCoefficientDays + "\n");
        parameters.Write("drugRemovalRateVesselDays, " + drugRemovalRateVesselDays + "\n");
        parameters.Write("drugConcentrationVessel, " + drugConcentrationVessel + "\n");
        parameters.Write("proliferationDiffusionCoefficientDays, " + proliferationDiffusionCoefficientDays + "\n");
        parameters.Write("autocrineProliferationSignalProductionDays, " + autocrineProliferationSignalProductionDays + "\n");
        parameters.Write("paracrineProliferationSignalProductionDays, " + paracrineProliferationSignalProductionDays + "\n");
        parameters.Write("proliferationDegradationDrugDays, " + proliferationDegradationDrugDays + "\n");
        parameters.Write("proliferationCancerSite, " + proliferationCancerSite + "\n");
        parameters.Write("activationDistMeanDays, " + activationDistMeanDays + "\n");
        parameters.Write("activationDistStdDevDays, " + activationDistStdDevDays + "\n");
        parameters.Write("deactivationDistMeanDays, " + deactivationDistMeanDays + "\n");
        parameters.Write("deactivationDistStdDevDays, " + deactivationDistStdDevDays + "\n");
        parameters.Write("initialTreatmentDays, " + initialTreatmentDays + "\n");
        parameters.Write("treatmentDaysDrugA, " + (double)treatmentDaysDrugA + "\n");
        parameters.Write("holidayDaysDrugA, " + (double)holidayDaysDrugA + "\n");
        parameters.Close();
    }

    public int [] treatmentScheduleTimesteps(double timestepPerDay, int treatmentDaysDrug, int holidayDaysDrug) {
        int treatmentTimesteps = (int) (treatmentDaysDrug * timestepPerDay); // number of treatment timesteps for drug
        int holidayTimesteps = (int) (holidayDaysDrug* timestepPerDay); // number of holiday timesteps for drug
        int [] schedule = {treatmentTimesteps, holidayTimesteps};
        return schedule;
    }

    public void initialStromaCell(int iCord, int STROMA, double timestepPerDay) {
        // create initial stroma cell
        ExampleCell newCell = NewAgentSQ(iCord);
        newCell.type = STROMA;
        newCell.IMTAge = (int) (timestepPerDay*Math.random());
        newCell.activationProbability = Math.random();
    }

    public void initialCellCPCExample(int xCord, int yCord, int CELLTYPE) {
        // create initial stroma cell
        ExampleCell newCell = NewAgentSQ(xCord, yCord);
        newCell.type = CELLTYPE;
    }

    public void ageCells(int STROMA, int ACTIVATEDSTROMA, int CANCER, int UNREACTIVESTROMA) {
        ShuffleAgents(rng);
        // age stroma cell
        for (ExampleCell cell : this) {
            if (cell.type == ACTIVATEDSTROMA) {
                double currentAge = cell.IMTAge;
                cell.IMTAge = currentAge + 1;
                double currentDelay = cell.activationDelay;
                cell.activationDelay = currentDelay - 1;
            } else if (cell.type == STROMA || cell.type == CANCER || cell.type == UNREACTIVESTROMA) {
                double currentAge = cell.IMTAge;
                cell.IMTAge = currentAge + 1;
                double currentDelay = cell.deactivationDelay;
                cell.deactivationDelay = currentDelay - 1;
            }
        }
    }

    public void stromaDieCells(double stromaDieProb, int STROMA, int ACTIVATEDSTROMA, int UNREACTIVESTROMA){
        ShuffleAgents(rng);
        // check for stroma cell type and turnover
        for (ExampleCell cell:this) {
            if ((cell.type == STROMA || cell.type == ACTIVATEDSTROMA || cell.type == UNREACTIVESTROMA)) {
                cell.stromaDieCell(stromaDieProb);
            }
        }
    }

    public void stromaStatusUpdate(int STROMA, int ACTIVATEDSTROMA, double stromaDrugThresholdReactivity,
                                                      int CANCER, int QUIESCENTCANCER, double stromaActivationNormalisationProbability,
                                                            double activationDistMean, double activationDistStdDev, ArrayList activationDelayTime,
                                   double deactivationDistMean, double deactivationDistStdDev, ArrayList deactivationDelayTime) {
        ShuffleAgents(rng);
        // check for stroma cell neighbouring cancer cells and drug threshold
        for (ExampleCell cell:this) {
            boolean cancerNeighbour = false;
            // activate stroma if it has cancer neighbour and drug is equal or above threshold
            cancerNeighbour = cell.checkCancerNeighbour(CANCER, QUIESCENTCANCER);
            if (cell.type == STROMA && cancerNeighbour == true) {
                // remove test for drug concentration so that it is alway met.
                double drugAmount = drug.Get(cell.Isq());
                if (drugAmount >= stromaDrugThresholdReactivity && rng.Double() < stromaActivationNormalisationProbability) {
                    cell.type = ACTIVATEDSTROMA;
                    cell.activationDelay = rng.Gaussian(activationDistMean, activationDistStdDev);
                    activationDelayTime.add(cell.activationDelay);
                //System.out.println("stroma activated!");
                }
                // deactivate stroma if drug is under threshold
            } else if (cell.type == ACTIVATEDSTROMA) {
                double drugAmount = drug.Get(cell.Isq());
                if (drugAmount < stromaDrugThresholdReactivity) {
                    cell.type = STROMA;
                    cell.deactivationDelay = rng.Gaussian(deactivationDistMean, deactivationDistStdDev);
                    deactivationDelayTime.add(cell.deactivationDelay);
                }
            }
        }
    }

    public void initialiseDrugValues(int VESSEL, double drugConcentrationVessel) {
        // initialise drug from vessel sites
        for (ExampleCell cell:this) {
            if (cell.type == VESSEL) {
                drug.Set(cell.Xsq(), cell.Ysq(), drugConcentrationVessel);
            }
        }
        drug.Update();
    }

    public void updatePDEValues(int x, int y, int timeCounter, int treatmentTimestepDrug, int holidayTimestepDrug,
                                int VESSEL, double drugDiffusionCoefficientTimestep, double drugConcentrationVessel,
                                double drugRemovalRateVesselTimestep, double proliferationDiffusionCoefficientTimestep,
                                int CANCER, int QUIESCENTCANCER, double autocrineProliferationSignalProductionTimestep, int STROMA, int ACTIVATEDSTROMA,
                                double paracrineProliferationSignalProductionTimestep,
                                double proliferationDegradationDrugTimestep, int initialTreatmentDays, double timestepPerDay) {

        boolean treatmentOn = false;
        // initial treatment period
        int initTreatmentTimesteps = initialTreatmentTimesteps(initialTreatmentDays, timestepPerDay);
        if (timeCounter < initTreatmentTimesteps) {
            treatmentOn = true;
            // check if treatment is on
        } else if ((timeCounter - initTreatmentTimesteps) % (treatmentTimestepDrug + holidayTimestepDrug) < treatmentTimestepDrug) {
            //System.out.println("treatment on timestep: " + timeCounter);
            treatmentOn = true;
        }


        // apply diffusion to drug according to treatment schedule
        if (treatmentOn == true) {
            initialiseDrugValues(VESSEL, drugConcentrationVessel);
            //drug.Diffusion(drugDiffusionCoefficientTimestep);
            drug.DiffusionADI(drugDiffusionCoefficientTimestep);
        } else {
            //drug.Diffusion(drugDiffusionCoefficientTimestep);
            drug.DiffusionADI(drugDiffusionCoefficientTimestep);
        }

        // remove drug at vessel sites
        for (ExampleCell cell : this) {
            if (cell.type == VESSEL) {
                drug.Add(cell.Isq(), -drugRemovalRateVesselTimestep);
            }
        }

        // apply diffusion to proliferation signal
        //proliferation.Diffusion(proliferationDiffusionCoefficientTimestep);
        proliferation.DiffusionADI(proliferationDiffusionCoefficientTimestep);
        // autocrine increase
        for (ExampleCell cell:this) {
            if (cell.type == CANCER || cell.type == QUIESCENTCANCER) {
                proliferation.Add(cell.Isq(), autocrineProliferationSignalProductionTimestep);
            }
        }
        // paracrine increase
        for (ExampleCell cell:this) {
            if (cell.type == CANCER || cell.type == QUIESCENTCANCER) {
                // check for number of activated stroma neighbours
                int activatedStromaNeighbourCn = cell.checkActivatedStromaNeighbour(STROMA, ACTIVATEDSTROMA);
                double paracrineInput = activatedStromaNeighbourCn * paracrineProliferationSignalProductionTimestep;
                proliferation.Add(cell.Isq(), paracrineInput);
            }
        }
        // degradation due to drug
        for (int i = 0; i < x * y; i++) {
            double drugDegradationValue = -proliferationDegradationDrugTimestep*drug.Get(i);
            proliferation.Add(i, drugDegradationValue);
        }
        // create arrays of PDE values for bounds test
        double[] drugCurrentValues = drug.GetField();
        double[] drugDeltaValues = drug.GetDeltas();
        double[] proliferationCurrentValues = proliferation.GetField();
        double[] proliferationDeltaValues = proliferation.GetDeltas();

        /* check for negative drug values and set to zero and check for negative proliferation values and set to zero
            and check for proliferation values greater than one and set to one*/
        for (int i = 0; i < x * y; i++) {
            if (drugCurrentValues[i] + drugDeltaValues[i] < 0) {
                drug.Set(i, 0);
            } else if (drugCurrentValues[i] + drugDeltaValues[i] > 1) {
                System.out.println("Drug value exceeds 1!");
            }
            if (proliferationCurrentValues[i] + proliferationDeltaValues[i] < 0) {
                proliferation.Set(i, 0);
            }
            if (proliferationCurrentValues[i] + proliferationDeltaValues[i] > 1) {
                proliferation.Set(i, 1);
            }
        }
        // update drug grid
        drug.Update();
        // update proliferation grid
        proliferation.Update();
    }

    public int initialTreatmentTimesteps(int initialTreatmentDays, double timestepPerDay) {
        //double multiplier = initialTreatmentDays;
        int timestepInitialTreatment = (int) Math.round(initialTreatmentDays * timestepPerDay);
        //System.out.println("Initial treatment timesteps: " + timestepInitialTreatment);
        return timestepInitialTreatment;
    }

    public double drugMeanField() {
        double mFValue = drug.GetAvg();
        return mFValue;
    }

    public void cancerDieProliferationUpdate(int CANCER, int QUIESCENTCANCER, boolean cancerEliminated) {
        ShuffleAgents(rng);
        for (ExampleCell cell:this) {
            if (cell.type == CANCER || cell.type == QUIESCENTCANCER) {
                double proliferationAmount = proliferation.Get(cell.Isq());
                if (proliferationAmount < cell.hd) {
                    // cancer cell dies
                    cell.cancerDieCell(CANCER, QUIESCENTCANCER, cancerEliminated);
                    //System.out.println("Cancer Dies");
                } else if (proliferationAmount < cell.hp) {
                    // cancer cell is quiescent
                    cell.type = QUIESCENTCANCER;
                } else if (proliferationAmount >= cell.hp) {
                    // cancer cells is not quiescent
                    cell.type = CANCER;
                }
            }
        }
    }

    public void divideCells(int stromaNCI, int STROMA, int ACTIVATEDSTROMA, int stromaIMT, int CANCER, int cancerIMT, double proliferationCancerSite, int UNREACTIVESTROMA, int[][] vesselSetUp){
        ShuffleAgents(rng);
        // check cell type and its IMTAge and divide
        for (ExampleCell cell:this) {
            if ((cell.type == STROMA || cell.type == ACTIVATEDSTROMA || cell.type == UNREACTIVESTROMA) && cell.IMTAge >= stromaIMT) {
                cell.divideStromaCell(stromaNCI);
            } else if (cell.type == CANCER && cell.IMTAge >= cell.cancerVariedIMT) {
                cell.divideCancerCell(cancerIMT, proliferationCancerSite, vesselSetUp);
            }
        }
    }

    public void CPCExampleDivideCells(double agentCn){
        ShuffleAgents(rng);
        int countA = (int)(agentCn);
        // check cell type and its IMTAge and divide
        for (int i = 0; i < countA; i++) {
            int randomAgent = (int) (countA * Math.random());
            int iterCount = 0;
            for (ExampleCell cell:this) {
                if (iterCount == randomAgent) {
                    cell.divideCPCExampleCell();
                }
                iterCount++;
            }
        }
    }
    public void checkElimination(int CANCER, int QUIESCENTCANCER, boolean cancerEliminated) {
        int cancerCn = 0;
        for (ExampleCell cell:this) {
            if (cell.type == CANCER || cell.type == QUIESCENTCANCER) {
                cancerCn++;
            }
        }
        if (cancerCn == 0) {
            System.out.println("cancer eliminated!");
            cancerEliminated = true;
            //System.out.println("In checkElimination if statement - " + cancerEliminated);
        }
        //System.out.println("In checkElimination outside if statment - " + cancerEliminated);
    }

    public boolean testElimination(int CANCER, int QUIESCENTCANCER) {
        int cancerCn = 0;
        boolean cancerEliminated = false;
        for (ExampleCell cell:this) {
            if (cell.type == CANCER || cell.type == QUIESCENTCANCER) {
                cancerCn++;
            }
        }
        if (cancerCn == 0) {
            System.out.println("cancer eliminated!");
            cancerEliminated = true;
            //System.out.println("In checkElimination if statement - " + cancerEliminated);
        }
        return cancerEliminated;
        //System.out.println("In checkElimination outside if statment - " + cancerEliminated);
    }


/*    static int STROMA = RGB256(246, 189, 76); // D5 set colour for stroma cells
    static int VESSEL = RGB256(93, 14, 40); // D8 set colour for vessel cells
    static int CANCER = RGB256(203,0,93); // D7 set colour for cancer cells
    static int ACTIVATEDSTROMA = RGB256(0, 60, 115); // D2 set colour for activated stroma cells
    static int QUIESCENTCANCER = RGB256(255, 102, 102); // D6 set colour for quiescent cancer cells
    static int UNREACTIVESTROMA = RGB256(204, 204, 204); // set colour for unreactive stroma
*/


    public void drawModelLayers(GridWindow agentLayer, GridWindow drugLayer, GridWindow proliferationLayer) {
        for (int i = 0; i < length; i++) {
            // set background colour
            int agentColor = Util.WHITE;
            int drugColor = Util.WHITE;
            int prolifColor = Util.WHITE;
            // colour each agent accoring to type and any vacant cell will take PDE grid colour
            ExampleCell cell=GetAgent(i);
            if (GetAgent(i)!=null) {
                // STROMA
                if (cell.type == -1) {
                    agentColor = RGB256(246, 189, 76);
                }
                // VESSEL
                else if (cell.type == -2) {
                    agentColor = RGB256(93, 14, 40);
                }
                // CANCER
                else if (cell.type == -3) {
                    agentColor = RGB256(203,0,93);
                }
                // ACTIVATEDSTROMA
                else if (cell.type == -4) {
                    agentColor = RGB256(0, 60, 115); // Blue
                    //agentColor = RGB256(0, 255, 255);
                }
                // QUIESCENTCANCER
                else if (cell.type == -5) {
                    agentColor = RGB256(255, 102, 102);
                }
                // UNREACTIVESTROMA
                else if (cell.type == -6) {
                    agentColor = RGB256(204, 204, 204);
                }
            }
            drugColor = HeatMapBGR(drug.Get(i)); // turn on for drug
            prolifColor = HeatMapRBG(proliferation.Get(i)); // turn on for proliferation signal
            //System.out.println(prolifColor);
            // set up display
            agentLayer.SetPix(i,agentColor);
            drugLayer.SetPix(i,drugColor);
            proliferationLayer.SetPix(i,prolifColor);
        }
    }

    public void drawModelAgents(GridWindow agentLayer) {
        for (int i = 0; i < length; i++) {
            // set background colour
            int agentColor = Util.WHITE;
            // colour each agent accoring to type and any vacant cell will take PDE grid colour
            ExampleCell cell=GetAgent(i);
            if (GetAgent(i)!=null) {
                // STROMA
                if (cell.type == -1) {
                    agentColor = RGB256(246, 189, 76);
                }
                // VESSEL
                else if (cell.type == -2) {
                    agentColor = RGB256(93, 14, 40);
                }
                // CANCER
                else if (cell.type == -3) {
                    agentColor = RGB256(203,0,93);
                }
                // ACTIVATEDSTROMA
                else if (cell.type == -4) {
                    agentColor = RGB256(0, 60, 115);
                }
                // QUIESCENTCANCER
                else if (cell.type == -5) {
                    agentColor = RGB256(255, 102, 102);
                }
                // UNREACTIVESTROMA
                else if (cell.type == -6) {
                    agentColor = RGB256(204, 204, 204);
                }
            }
            //System.out.println(prolifColor);
            // set up display
            agentLayer.SetPix(i,agentColor);
        }
    }


    public void drawModelLayersSmall(int x, int y, int xSmall, int ySmall, GridWindow agentLayer, GridWindow drugLayer, GridWindow proliferationLayer) {
        for (int i = 0; i < length; i++) {
            int cX = i / x;
            int cY = i % y;
            int xBd_L = (int)((x - xSmall) / 2);
            int xBd_U = xBd_L + xSmall;
            int yBd_L = (int)((y - ySmall) / 2);
            int yBd_U = yBd_L + ySmall;
            if (((cX >= xBd_L) && (cX < xBd_U)) && ((cY >= yBd_L) && (cY < yBd_U))) {
                // set background colour
                int agentColor = Util.WHITE;
                int drugColor = Util.WHITE;
                int prolifColor = Util.WHITE;
                // colour each agent accoring to type and any vacant cell will take PDE grid colour
                ExampleCell cell=GetAgent(i);
                if (GetAgent(i)!=null) {
                    // STROMA
                    if (cell.type == -1) {
                        agentColor = RGB256(246, 189, 76);
                    }
                    // VESSEL
                    else if (cell.type == -2) {
                        agentColor = RGB256(93, 14, 40);
                    }
                    // CANCER
                    else if (cell.type == -3) {
                        agentColor = RGB256(203,0,93);
                    }
                    // ACTIVATEDSTROMA
                    else if (cell.type == -4) {
                        agentColor = RGB256(0, 60, 115);
                    }
                    // QUIESCENTCANCER
                    else if (cell.type == -5) {
                        agentColor = RGB256(255, 102, 102);
                    }
                    // UNREACTIVESTROMA
                    else if (cell.type == -6) {
                        agentColor = RGB256(204, 204, 204);
                    }
                }
                drugColor = HeatMapBGR(drug.Get(i)); // turn on for drug
                prolifColor = HeatMapRBG(proliferation.Get(i)); // turn on for proliferation signal
                //System.out.println(prolifColor);
                // set up display
                int xPos = cX - xBd_L;
                int yPos = cY - yBd_L;
                agentLayer.SetPix(xPos, yPos ,agentColor);
                drugLayer.SetPix(xPos, yPos ,drugColor);
                proliferationLayer.SetPix(xPos, yPos ,prolifColor);
            }
        }
    }

    public int countInitialStroma(String initialStromaCn) throws IOException {
        // Get initial stroma count from homeostatic run
        // Reads file and stores as string
        FileReader file = new FileReader(new File(initialStromaCn));
        // declare reading variables
        BufferedReader br = new BufferedReader(file);
        String temp = br.readLine();
        // Tokenise each element in string
        StringTokenizer defaultTokenizer = new StringTokenizer(temp);
        // Store each token as element in double array
        int setUp;
        int convertedValue = Integer.parseInt(defaultTokenizer.nextToken());
        setUp = convertedValue;
        return setUp;
    }



    public double[][] initialStromaSetUp(String initialStromaPositions, int initialStromaCount) throws IOException {
        // Get initial stroma sites from homeostasis run
        // Reads file and stores as string
        FileReader file = new FileReader(new File(initialStromaPositions));
        // declare reading variables
        BufferedReader br = new BufferedReader(file);
        String temp = "";
        double[][] setUp = new double[initialStromaCount][2];
        int i = 0;
        int j = 0;
        while ((temp = br.readLine()) != null) {
            // Tokenise each element in string
            StringTokenizer defaultTokenizer = new StringTokenizer(temp, ",");
            // Store each token as element in double array
            while (defaultTokenizer.hasMoreTokens()) {
                double convertedValue = Double.parseDouble(defaultTokenizer.nextToken());
                setUp[i][j] = convertedValue;
                j++;
            }
            j = 0;
            i++;
        }
        return setUp;
    }

    public double[][] initialPassiveReactiveStromaSetUp(String initialStromaPositions, int initialStromaCount) throws IOException {
        // Get initial stroma sites from homeostasis run
        // Reads file and stores as string
        FileReader file = new FileReader(new File(initialStromaPositions));
        // declare reading variables
        BufferedReader br = new BufferedReader(file);
        String temp = "";
        double[][] setUp = new double[initialStromaCount][3];
        int i = 0;
        int j = 0;
        while ((temp = br.readLine()) != null) {
            // Tokenise each element in string
            StringTokenizer defaultTokenizer = new StringTokenizer(temp, ",");
            // Store each token as element in double array
            while (defaultTokenizer.hasMoreTokens()) {
                double convertedValue = Double.parseDouble(defaultTokenizer.nextToken());
                setUp[i][j] = convertedValue;
                j++;
            }
            j = 0;
            i++;
        }
        return setUp;
    }


    public void initialiseStroma(double[][] homeostaticStroma, int STROMA, int UNREACTIVESTROMA, double reactiveStromaProportion) {
        for (int i = 0; i < homeostaticStroma.length; i++) {
            // create initial stroma cell and assign attributes
            ExampleCell newCell = NewAgentSQ((int)homeostaticStroma[i][0]);
            newCell.activationProbability = Math.random();
            if (newCell.activationProbability <= reactiveStromaProportion) {
                newCell.type = STROMA;
            } else {
                newCell.type = UNREACTIVESTROMA;
            }
            newCell.IMTAge = homeostaticStroma[i][1];
        }
    }

    public void initialisePassiveReactiveStroma(double[][] homeostaticStroma, int STROMA, int UNREACTIVESTROMA, double reactiveStromaProportion) {
        for (int i = 0; i < homeostaticStroma.length; i++) {
            // create initial stroma cell and assign attributes
            ExampleCell newCell = NewAgentSQ((int)homeostaticStroma[i][0]);
            if (homeostaticStroma[i][2] == -1) {
                newCell.type = STROMA;
            } else {
                newCell.type = UNREACTIVESTROMA;
            }
            newCell.IMTAge = homeostaticStroma[i][1];
        }
    }


    public int countInitialCancer(String initialCancerCn) throws IOException {
        // Get initial stroma count from homeostatic run
        // Reads file and stores as string
        FileReader file = new FileReader(new File(initialCancerCn));
        // declare reading variables
        BufferedReader br = new BufferedReader(file);
        String temp = br.readLine();
        // Tokenise each element in string
        StringTokenizer defaultTokenizer = new StringTokenizer(temp);
        // Store each token as element in double array
        int setUp;
        int convertedValue = Integer.parseInt(defaultTokenizer.nextToken());
        setUp = convertedValue;
        return setUp;
    }

    public double[][] initialCancerSetUp(String initialCancerPositions, int initialCancerCount) throws IOException {
        // Get initial cancer sites from homeostasis run
        // Reads file and stores as string
        FileReader file = new FileReader(new File(initialCancerPositions));
        // declare reading variables
        BufferedReader br = new BufferedReader(file);
        String temp = "";
        double[][] setUp = new double[initialCancerCount][6];
        int i = 0;
        int j = 0;
        while ((temp = br.readLine()) != null) {
            // Tokenise each element in string
            StringTokenizer defaultTokenizer = new StringTokenizer(temp, ",");
            // Store each token as element in double array
            while (defaultTokenizer.hasMoreTokens()) {
                double convertedValue = Double.parseDouble(defaultTokenizer.nextToken());
                setUp[i][j] = convertedValue;
                j++;
            }
            j = 0;
            i++;
        }
        //System.out.println(Arrays.deepToString(setUp));
        return setUp;
    }

    public double[] initialProliferationSignal(String initialProliferationValues, int x, int y) throws IOException {
        // Get initial cancer sites from homeostasis run
        // Reads file and stores as string
        FileReader file = new FileReader(new File(initialProliferationValues));
        // declare reading variables
        BufferedReader br = new BufferedReader(file);
        String temp = "";
        double[] setUp = new double[x * y];
        int i = 0;
        while ((temp = br.readLine()) != null) {
            // Tokenise each element in string
            StringTokenizer defaultTokenizer = new StringTokenizer(temp, "/n");
            // Store each token as element in double array
            while (defaultTokenizer.hasMoreTokens()) {
                double convertedValue = Double.parseDouble(defaultTokenizer.nextToken());
                setUp[i] = convertedValue;
            }
            i++;
        }
        //System.out.println(Arrays.toString(setUp));
        return setUp;
    }

    public void meanFieldDrug(double[] timeVector,double[] meanFieldDrug, FileIO MFDrug) {
        MFDrug.Write("time,drug MF\n");
        for (int i = 0; i < meanFieldDrug.length; i++) {
            MFDrug.Write(timeVector[i] + "," + meanFieldDrug[i] + "\n");
        }
    }

    public void populationData(double[] timeVector,double[] stromaPopulation,double[] cancerPopulation,double[] activatedStromaPopulation,double[] contributingActivatedStromaPopulation,double[] quiescentCancerPopulation, FileIO populations) {
        populations.Write("time,stroma,cancer,activated stroma,contributing activated stroma,quiescent cancer\n");
        for (int i = 0; i < timeVector.length; i++) {
            populations.Write(timeVector[i] + "," + stromaPopulation[i] + "," + cancerPopulation[i] + "," + activatedStromaPopulation[i] + "," + contributingActivatedStromaPopulation[i] + "," + quiescentCancerPopulation[i] + "\n");
        }
    }

    public void quadrantPopulationData(double[] timeVector,double[] cancerPopulationQ1,double[] cancerPopulationQ2,double[] cancerPopulationQ3,double[] cancerPopulationQ4, FileIO quadrantPopulations) {
        quadrantPopulations.Write("time,Q1,Q2,Q3,Q4\n");
        for (int i = 0; i < timeVector.length; i++) {
            quadrantPopulations.Write(timeVector[i] + "," + cancerPopulationQ1[i] + "," + cancerPopulationQ2[i] + "," + cancerPopulationQ3[i] + "," + cancerPopulationQ4[i] + "\n");
        }
    }

    public void populationDataSTROMA(double[] timeVector,double[] stromaPopulation, FileIO populations) {
        populations.Write("time,stroma\n");
        for (int i = 0; i < timeVector.length; i++) {
            populations.Write(timeVector[i] + "," + stromaPopulation[i] + "\n");
        }
    }

    public void meanDelayTimeData(ArrayList activationDelayTime, FileIO ActivateDelay) {
        ActivateDelay.Write("activation delay\n");
        for (int i = 0; i < activationDelayTime.size(); i++) {
            ActivateDelay.Write(activationDelayTime.get(i) + "\n");
        }
    }

    public void initialiseCancer(double[][] initialCancer, int CANCER, int QUIESCENTCANCER, int cellIDCount, int[][]vesselSetUp, int clusterBoxSize, int VESSEL, int x, int y) {
        for (int i = 0; i < initialCancer.length; i++) {
            // create initial stroma cell and assign attributes
            ExampleCell newCell = NewAgentSQ((int)initialCancer[i][0]);
            newCell.cellID = cellIDCount;
            cellIDCount++;
            if ((int)initialCancer[i][1] == 11) {
                newCell.type = CANCER;
            } else if ((int)initialCancer[i][1] == 22) {
                newCell.type = QUIESCENTCANCER;
            }
            newCell.IMTAge = initialCancer[i][2];
            newCell.cancerVariedIMT = initialCancer[i][3];
// will need to fix this bit when growing cancer again!!!!
            newCell.hd = 0.2;
            newCell.hp = 0.8;
            newCell.distanceNearestVessel = nearestVesselDist(vesselSetUp, newCell.Xsq(), newCell.Ysq());
            newCell.vesselClusterValue = clusterValue(clusterBoxSize, newCell.Xsq(), newCell.Ysq(), VESSEL, x, y);
        }
    }

    public void neighboursActivatedStroma(int CANCER, int STROMA, int ACTIVATEDSTROMA) {
        for (ExampleCell cell:this) {
            if (cell.type == CANCER) {
                cell.activatedStromaNeighbours = cell.checkActivatedStromaNeighbour(STROMA, ACTIVATEDSTROMA);
            }
        }
    }

    public double clusterValue(int clusterBoxSize, int xPos, int yPos, int VESSEL, int x, int y) {
        int vesselCn = 0;
        for (int i = 0; i < clusterBoxSize; i++) {
            for (int j = 0; j < clusterBoxSize; j++) {
                if ((((xPos - (i+1) >= 0) && (xPos - (i+1) < x)) && ((yPos - (j+1) >= 0) && (yPos - (j+1) < y))) && (GetAgent(xPos - (i+1), yPos - (j+1)) != null)) {
                    if (GetAgent(xPos - (i + 1), yPos - (j + 1)).type == VESSEL) {
                        vesselCn++;
                    }
                }
                if ((((xPos + (i+1) >= 0) && (xPos + (i+1) < x)) && ((yPos + (j+1) >= 0) && (yPos + (j+1) < y))) && (GetAgent(xPos + (i+1), yPos + (j+1)) != null)) {
                    if (GetAgent(xPos + (i + 1), yPos + (j + 1)).type == VESSEL) {
                        vesselCn++;
                    }
                }
            }
        }
        double vesselClusterDensity = vesselCn / (2 * clusterBoxSize + 1)^2;
        return vesselClusterDensity;
    }

    public double nearestVesselDist(int[][] vesselSetUp, int xCoord, int yCoord) {
        double [] distance = new double [vesselSetUp.length];
        for (int i = 0; i < vesselSetUp.length; i++) {
            //System.out.println(Math.pow(xCoord - vesselSetUp[i][0], 2) + Math.pow(yCoord - vesselSetUp[i][1], 2));
            distance[i] = Math.sqrt(Math.pow(xCoord - vesselSetUp[i][0], 2) + Math.pow(yCoord - vesselSetUp[i][1], 2));
        }
        double minDist = getMin(distance);
        return minDist;
    }

    public static double getMin(double[] inputArray){
        double minValue = inputArray[0];
        for(int i=1;i<inputArray.length;i++){
            if(inputArray[i] < minValue){
                minValue = inputArray[i];
            }
        }
        return minValue;
    }

    public void initialiseProliferationSignal(double[] proliferationSignalInit) {
        // update proliferation signal PDEGrid
        for (int i = 0; i < proliferationSignalInit.length; i++) {
            proliferation.Add(i, proliferationSignalInit[i]);
        }
        proliferation.Update();
    }

    public void finalProliferationSignal(FileIO proliferationSignalValues) {
        double[] proliferationFinalValues = proliferation.GetField();
        for (int i = 0; i < proliferationFinalValues.length; i++) {
            proliferationSignalValues.Write(proliferationFinalValues[i] + "\n");
        }
    }

    public void initialCancerEvent(int x, int y, int nV, int [][] LocVessels, int QUIESCENTCANCER, int cancerIMT, double proliferationCancerSite) {
        // locate vessel cells in the centre of the grid
        int centreX = x / 2;
        int centreY = y / 2;
        // declare radius counter
        int centreRadius = 0;
        // declare cancer seeded boolean
        boolean cancerSeeded = false;
        // determine where to seed cancer cell
        while (cancerSeeded == false && centreRadius <= Math.min(centreX, centreY)) {
            // declare vessels in radius counter
            int centreVesselsCn = 0;
            // declare outer radius
            int newCentreRadius = centreRadius + 1;
            if (centreRadius < Math.min(centreX, centreY)) {
                for (int i = 0; i < nV; i++) {
                    // count vessels in radius
                    double distMeasure = Math.sqrt(Math.pow(LocVessels[i][0] - centreX, 2)
                            + Math.pow(LocVessels[i][1] - centreY, 2));
                    if (distMeasure > centreRadius && distMeasure <= newCentreRadius) {
                        centreVesselsCn++;
                    }
                }
                // declare storage array for vessels in radius
                int[][] centreVessels = new int[centreVesselsCn][2];
                // declare position counter
                int positionCn = 0;
                for (int i = 0; i < nV; i++) {
                    double distMeasure = Math.sqrt(Math.pow(LocVessels[i][0] - centreX, 2)
                            + Math.pow(LocVessels[i][1] - centreY, 2));
                    // store location of vessel site in radius
                    if (distMeasure > centreRadius && distMeasure <= newCentreRadius) {
                        centreVessels[positionCn][0] = LocVessels[i][0];
                        centreVessels[positionCn][1] = LocVessels[i][1];
                        positionCn++;
                    }
                }
                for (int i = 0; i < centreVesselsCn; i++) {
                    // IS THERE ANOTHER WAY TO DO THIS???? problem is cancerSeeded is false in while loop and here!
                    if (cancerSeeded == false) {
                        // seed cancer next to vessel in radius
                        ExampleCell centreVesselCell = GetAgent(centreVessels[i][0], centreVessels[i][1]);
                        cancerSeeded = centreVesselCell.seedCancer(QUIESCENTCANCER, cancerSeeded, cancerIMT, proliferationCancerSite);
                    }
                }
            } else {
                // indicate there are no valid positions to seed cancer
                System.out.println("No valid positions to seed cancer!");
            }
            // increase centre radius counter
            centreRadius++;
        }
    }


    public int [][] readVesselInitialisation(String initialVesselPositions, int nV) throws IOException {
        FileReader file = new FileReader(new File(initialVesselPositions));
        // declare reading variables
        BufferedReader br = new BufferedReader(file);
        String temp = "";
        int[][] setUp = new int[nV][2];
        int i = 0;
        int j = 0;
        while ((temp = br.readLine()) != null) {
            // Tokenise each element in string
            StringTokenizer defaultTokenizer = new StringTokenizer(temp, ",");
            // Store each token as element in double array
            while (defaultTokenizer.hasMoreTokens()) {
                int convertedValue = Integer.parseInt(defaultTokenizer.nextToken());
                setUp[i][j] = convertedValue;
                j++;
            }
            j = 0;
            i++;
        }
        return setUp;
    }

    public void initialiseVesselCells(int[][] vesselSetUp, int nV, int VESSEL) {
        for (int j = 0; j < nV; j++) {
            NewAgentSQ(vesselSetUp[j][0], vesselSetUp[j][1]).type = VESSEL;
        }
    }

    public int countHomeostaticStroma(String initialHomeostaticStromaCount) throws IOException {
        // Get initial stroma count from homeostatic run
        // Reads file and stores as string
        FileReader file = new FileReader(new File(initialHomeostaticStromaCount));
        // declare reading variables
        BufferedReader br = new BufferedReader(file);
        String temp = br.readLine();
        // Tokenise each element in string
        StringTokenizer defaultTokenizer = new StringTokenizer(temp);
        // Store each token as element in double array
        int setUp;
        int convertedValue = Integer.parseInt(defaultTokenizer.nextToken());
        setUp = convertedValue;
        return setUp;
    }

    public double[][] homeostaticStromaSetUp(String initialHomeostaticStroma, int homestaticStromaCount) throws IOException {
        // Get initial stroma sites from homeostasis run
        // Reads file and stores as string
        FileReader file = new FileReader(new File(initialHomeostaticStroma));
        // declare reading variables
        BufferedReader br = new BufferedReader(file);
        String temp = "";
        double[][] setUp = new double[homestaticStromaCount][2];
        int i = 0;
        int j = 0;
        while ((temp = br.readLine()) != null) {
            // Tokenise each element in string
            StringTokenizer defaultTokenizer = new StringTokenizer(temp, ",");
            // Store each token as element in double array
            while (defaultTokenizer.hasMoreTokens()) {
                double convertedValue = Double.parseDouble(defaultTokenizer.nextToken());
                setUp[i][j] = convertedValue;
                j++;
            }
            j = 0;
            i++;
        }
        return setUp;
    }

// need to fix this bit when start from scratch
    public void stromaBoxes(int x, int y, int STROMA) {
        // Leave two blank regions in tissue
        int lowerBoundX = x / 5;
        int upperBoundX = 4 * x / 5;
        int lowerBoundY1 = y / 5;
        int upperBoundY1 = 2 * y / 5;
        int lowerBoundY2 = 3 * y / 5;
        int upperBoundY2 = 4 * y / 5;
        for (ExampleCell cell:this) {
            if (cell.type == STROMA && ((cell.Xsq() >= lowerBoundX && cell.Xsq() < upperBoundX) &&
                    ((cell.Ysq() >= lowerBoundY1 && cell.Ysq() < upperBoundY1) || (cell.Ysq() >= lowerBoundY2
                            && cell.Ysq() < upperBoundY2)))) {
                cell.stromaDieCell(1);
            }
        }
        // Fill one region with stroma cells
        for (int j = lowerBoundX; j < upperBoundX; j++) {
            for (int k = lowerBoundY1; k < upperBoundY1; k++) {
                int xCord = j;
                int yCord = k;
                if (GetAgent(xCord, yCord) != null) {
                } else {
                    ExampleCell newCell = NewAgentSQ(xCord, yCord);
                    newCell.type = STROMA;
                    newCell.IMTAge = 600 + (2*Math.random() - 1)*300;
                }
            }
        }
    }

    public void stromaLocations(int STROMA, int UNREACTIVESTROMA, FileIO locationsStroma) {
        for (ExampleCell cell:this) {
            if (cell.type == STROMA || cell.type == UNREACTIVESTROMA) {
                locationsStroma.Write(cell.Isq() + "," + cell.IMTAge + "\n");
            }
        }
    }

    public void stromaPassiveReactiveLocations(int STROMA, int UNREACTIVESTROMA, FileIO locationsPassiveReactiveStroma) {
        for (ExampleCell cell:this) {
            if (cell.type == STROMA || cell.type == UNREACTIVESTROMA) {
                locationsPassiveReactiveStroma.Write(cell.Isq() + "," + cell.IMTAge + "," + cell.type + "\n");
            }
        }
    }


    public void cancerLocations(int CANCER, int QUIESCENTCANCER, FileIO locationsCancer) {
        for (ExampleCell cell:this) {
            if (cell.type == CANCER) {
                int cellType = 11;
                locationsCancer.Write(cell.Isq() + "," + cellType + "," + cell.IMTAge + "," + cell.cancerVariedIMT + "," + cell.hd + "," + cell.hp + "\n");
            } else if (cell.type == QUIESCENTCANCER) {
                int cellType = 22;
                locationsCancer.Write(cell.Isq() + "," + cellType + "," + cell.IMTAge + "," + cell.cancerVariedIMT + "," + cell.hd + "," + cell.hp + "\n");
            }
        }
    }


    public double populationType(int TYPE) {
        double typeCount = 0;
        for (ExampleCell cell:this) {
            if (cell.type == TYPE) {
                typeCount++;
            }
        }
    return typeCount;
    }

    public double [] cancerCountQuadrants(int x, int y, int CANCER, int QUIESCENTCANCER) {
        //double [] quadrantCounts =  new double[4];
        double countQ1 = 0;
        double countQ2 = 0;
        double countQ3 = 0;
        double countQ4 = 0;
        for (ExampleCell cell:this) {
            if ((cell.type == CANCER) || (cell.type == QUIESCENTCANCER)) {
                if ((cell.Xsq() >= 0) && (cell.Xsq() < x / 2) && (cell.Ysq() >= y / 2) && (cell.Ysq() < y)){
                    countQ1++;
                } else if ((cell.Xsq() >= x / 2) && (cell.Xsq() < x) && (cell.Ysq() >= y / 2) && (cell.Ysq() < y)) {
                    countQ2++;
                } else if ((cell.Xsq() >= x / 2) && (cell.Xsq() < x) && (cell.Ysq() >= 0) && (cell.Ysq() < y)) {
                    countQ3++;
                } else {
                    countQ4++;
                }
            }
        }
        double [] quadrantCounts = {countQ1, countQ2, countQ3, countQ4};
        return quadrantCounts;
    }

    public double populationTypeSmall(int TYPE) {
        double typeCount = 0;
        for (ExampleCell cell:this) {
            if (cell.type == TYPE && ((cell.Xsq() >= 100 && cell.Xsq() < 200) && (cell.Ysq() >= 100 && cell.Ysq() < 200))) {
                typeCount++;
            }
        }
        return typeCount;
    }

    public double contributingActivatedStromaPopulationType(int STROMA, int ACTIVATEDSTROMA) {
        double typeCount = 0;
        for (ExampleCell cell:this) {
            if ((cell.type == ACTIVATEDSTROMA && cell.activationDelay < 0) || (cell.type == STROMA && cell.deactivationDelay > 0)) {
                typeCount++;
            }
        }
        return typeCount;
    }

    public double contributingActivatedStromaPopulationTypeSmall(int STROMA, int ACTIVATEDSTROMA) {
        double typeCount = 0;
        for (ExampleCell cell:this) {
            if (((cell.type == ACTIVATEDSTROMA && cell.activationDelay < 0) || (cell.type == STROMA && cell.deactivationDelay > 0)) && ((cell.Xsq() >= 100 && cell.Xsq() < 200) && (cell.Ysq() >= 100 && cell.Ysq() < 200))) {
                typeCount++;
            }
        }
        return typeCount;
    }

    public void recordAgentData(String resultsFolder, int CANCER, int QUIESCENTCANCER, int STROMA, int ACTIVATEDSTROMA, int i, int j, int k,  String runDateTime) {
        // create csv files
        FileIO agentCancer = new FileIO(resultsFolder + "/Data/CancerData_" + runDateTime + "_run_" + k + "_sim_" + i + "_timestep_" + j + ".csv", "w");
        FileIO agentStroma = new FileIO(resultsFolder + "/Data/StromaData_" + runDateTime + "_run_" + k + "_sim_" + i + "_timestep_" + j + ".csv", "w");
        // Write data headings
        agentCancer.Write("cellID,cell type,cell position,cell varied IMT,cell age,cell hd,cell hp,distance nearest vessel,vessel density,activated stroma neighbours\n");
        agentStroma.Write("cell type,cell position,cell age,cell activation probability\n");
        for (ExampleCell cell:this) {
            if (cell.type == CANCER || cell.type == QUIESCENTCANCER) {
                agentCancer.Write(cell.cellID + "," + cell.type + "," + cell.Isq() + "," + cell.cancerVariedIMT + "," + cell.IMTAge + "," + cell.hd + "," + cell.hp + "," + cell.distanceNearestVessel + "," + cell.vesselClusterValue +  "," + cell.activatedStromaNeighbours + "\n");
            } else if (cell.type == STROMA || cell.type == ACTIVATEDSTROMA) {
                agentStroma.Write(cell.type + "," + cell.Isq() + "," + cell.IMTAge + "," + cell.activationProbability + "\n");
            }
        }
        // close timestep data files
        agentCancer.Close();
        agentStroma.Close();
    }

    public void positionData(int i, int j, int k,  int x, int y, String runDateTime, String resultsFolder) {
        // create csv files
        FileIO positionInformation = new FileIO(resultsFolder + "/Data/PositionData_" + runDateTime + "_run_" + k + "_sim_" + i + "_timestep_" + j + ".csv", "w");
        // Write data headings
        positionInformation.Write("cell position,cell type,drug value,proliferation value\n");
        for (int l = 0; l < x * y; l++) {
            //System.out.println(k);
            if (GetAgent(l) == null) {
                positionInformation.Write(l + "," + "0" + "," + drug.Get(l) + "," + proliferation.Get(l) + "\n");
            } else {
                positionInformation.Write(l + "," + GetAgent(l).type + "," + drug.Get(l) + "," + proliferation.Get(l) + "\n");
            }
        }
        // close timestep data files
        positionInformation.Close();
    }

    public void positionDataCPCExample(int j, int x, int y, String resultsFolder) {
        // create csv files
        FileIO positionInformation = new FileIO(resultsFolder + "/Data/PositionData_timestep_" + j + ".csv", "w");
        // Write data headings
        positionInformation.Write("cell position,cell type\n");
        for (int l = 0; l < x * y; l++) {
            //System.out.println(k);
            if (GetAgent(l) == null) {
                positionInformation.Write(l + "," + "0" + "\n");
            } else {
                positionInformation.Write(l + "," + GetAgent(l).type + "\n");
            }
        }
        // close timestep data files
        positionInformation.Close();
    }


    public int neighbourTypeCount(int position, int cellType) {
        // counts the number of agents of cellType in neighbourhood
        int typeCount = 0;
        int neighbourhood = MapHood(divHood, position);
        for (int i = 0; i < neighbourhood; i++) {
            if (GetAgent(divHood[i]) != null) {
                if (GetAgent(divHood[i]).type == cellType) {
                    typeCount++;
                }
            }
        }
        return typeCount;
    }

    public void recordNeighbours(int x, int y, int i, int j, int k, int STROMA, int ACTIVATEDSTROMA, int UNREACTIVESTROMA, int CANCER, int QUIESCENTCANCER, int VESSEL, String resultsFolder, String runDateTime) {
        // create csv files
        FileIO neighbourInformation = new FileIO(resultsFolder + "/Data/NeighbourData_" + runDateTime + "_run_" + k + "_sim_" + i + "_timestep_" + j + ".csv", "w");
        // Write data headings
        neighbourInformation.Write("cell position,cancer neighbours,quiescent cancer neighbours,activated stroma neighbours,vessel neighbours,stroma neighbours,unreactive stroma neighbours,empty neighbours\n");
        for (int l = 0; l < x * y; l++) {
            int cancerNeighbourCount = neighbourTypeCount(l, CANCER);
            int quiescentCancerNeighbourCount = neighbourTypeCount(l, QUIESCENTCANCER);
            int activatedStromaNeighbourCount = neighbourTypeCount(l, ACTIVATEDSTROMA);
            int vesselNeighbourCount = neighbourTypeCount(l, VESSEL);
            int stromaNeighbourCount = neighbourTypeCount(l, STROMA);
            int unreactiveStromaNeighbourCount = neighbourTypeCount(l, UNREACTIVESTROMA);
            int emptyNeighbourCount = MapEmptyHood(divHood, l);
            neighbourInformation.Write(l + "," + cancerNeighbourCount + "," + quiescentCancerNeighbourCount + "," + activatedStromaNeighbourCount + "," + vesselNeighbourCount + "," + stromaNeighbourCount + "," + unreactiveStromaNeighbourCount + "," + emptyNeighbourCount + "\n");
            /*if (GetAgent(l) != null) {
                int cancerNeighbourCount = GetAgent(l).checkCancerNeighbourCount(CANCER, QUIESCENTCANCER);
                int activatedStromaNeighbourCount = GetAgent(l).checkActivatedStromaNeighbour(STROMA, ACTIVATEDSTROMA);;
                int vesselNeighbourCount = GetAgent(l).checkVesselNeighbour(VESSEL);;
                int stromaNeighbourCount = GetAgent(l).checkStromaNeighbour(STROMA, ACTIVATEDSTROMA, UNREACTIVESTROMA);;
                int emptyNeighbourCount = GetAgent(l).checkEmptyNeighbour();;
                neighbourInformation.Write(GetAgent(l).Isq() + "," + cancerNeighbourCount + "," + activatedStromaNeighbourCount + "," + vesselNeighbourCount + "," + stromaNeighbourCount + "," + emptyNeighbourCount + "\n");
            } else {*/
        }
        neighbourInformation.Close();
    }



    public void recordNeighbourInformation(int i, int j, int k,  int x, int y, int xSmall, int ySmall, String runDateTime, int CANCER, int ACTIVATEDSTROMA, int UNREACTIVESTROMA, String resultsFolder) {
        // Note that this only works for the xSmall and ySmall
        // create csv files
        FileIO neighbourInformation = new FileIO(resultsFolder + "/Data/NeighbourData_" + runDateTime + "_run_" + k + "_sim_" + i + "_timestep_" + j + ".csv", "w");
        // Write data headings
        neighbourInformation.Write("cell position,cancer neighbours,activated stroma neighbours,unreactive stroma neighbours\n");
        for (int l = 0; l < xSmall * ySmall; l++) {
            int cancerNeighbours = 0;
            int activatedStromaNeighbours = 0;
            int unreactiveStromaNeighbours = 0;
            int xBuff = (x - xSmall) / 2;
            int yBuff = (y - ySmall) / 2;
            int cX = (int) Math.floor(l / ySmall) + xBuff;
            int cY = l % xSmall + yBuff;
            for (int m = 0; m < 3; m++) {
                if (GetAgent(cX - 1, cY - 1 + l) != null && GetAgent(cX - 1, cY - 1 + l).type == CANCER) {
                    cancerNeighbours = cancerNeighbours + 1;
                } else if (GetAgent(cX - 1, cY - 1 + l) != null && GetAgent(cX - 1, cY - 1 + l).type == ACTIVATEDSTROMA) {
                    activatedStromaNeighbours = activatedStromaNeighbours + 1;
                } else if (GetAgent(cX - 1, cY - 1 + l) != null && GetAgent(cX - 1, cY - 1 + l).type == UNREACTIVESTROMA) {
                    unreactiveStromaNeighbours = unreactiveStromaNeighbours + 1;
                }
            }
            for (int m = 0; m < 2; m++) {
                if (GetAgent(cX, cY - 1 + 2 * l) != null && GetAgent(cX, cY - 1 + 2 * l).type == CANCER) {
                    cancerNeighbours = cancerNeighbours + 1;
                } else if (GetAgent(cX, cY - 1 + 2 * l) != null && GetAgent(cX, cY - 1 + 2 * l).type == ACTIVATEDSTROMA) {
                    activatedStromaNeighbours = activatedStromaNeighbours + 1;
                } else if (GetAgent(cX, cY - 1 + 2 * l) != null && GetAgent(cX, cY - 1 + 2 * l).type == UNREACTIVESTROMA) {
                    unreactiveStromaNeighbours = unreactiveStromaNeighbours + 1;
                }
            }
            for (int m = 0; m < 3; m++) {
                if (GetAgent(cX + 1, cY - 1 + l) != null && GetAgent(cX + 1, cY - 1 + l).type == CANCER) {
                    cancerNeighbours = cancerNeighbours + 1;
                } else if (GetAgent(cX + 1, cY - 1 + l) != null && GetAgent(cX + 1, cY - 1 + l).type == ACTIVATEDSTROMA) {
                    activatedStromaNeighbours = activatedStromaNeighbours + 1;
                } else if (GetAgent(cX + 1, cY - 1 + l) != null && GetAgent(cX + 1, cY - 1 + l).type == UNREACTIVESTROMA) {
                    unreactiveStromaNeighbours = unreactiveStromaNeighbours + 1;
                }
            }
            neighbourInformation.Write(l + "," + cancerNeighbours + "," + activatedStromaNeighbours + "," + unreactiveStromaNeighbours + "\n");
        }
        // close data file
        neighbourInformation.Close();
    }

    public void plotStromaPopulation(double [] timeVector, double [] stromaPopulation, double [] meanStromaAge, int timesteps) throws IOException {

        XYChart chart = new XYChart(500, 400);

        chart.setTitle("Stroma Population");
        chart.setXAxisTitle("time(days)");
        chart.setYAxisTitle("stroma cell count");

        // Customize Chart
        //chart.getStyler().setPlotBackgroundColor(ChartColor.getAWTColor(ChartColor.GREY));
        //chart.getStyler().setPlotGridLinesColor(new Color(255, 255, 255));
        chart.getStyler().setChartBackgroundColor(Color.WHITE);
        //chart.getStyler().setLegendBackgroundColor(Color.WHITE);
        //chart.getStyler().setChartFontColor(Color.RED);
        //chart.getStyler().setChartTitleBoxBackgroundColor(new Color(0, 222, 0));
        chart.getStyler().setChartTitleBoxVisible(false);
        chart.getStyler().setChartTitleBoxBorderColor(Color.BLACK);
        chart.getStyler().setPlotGridLinesVisible(false);
        //chart.getStyler().setLegendVisible(false);


        //chart.getStyler().setAxisTickPadding(20);

        //chart.getStyler().setAxisTickMarkLength(15);

        chart.getStyler().setAxisTicksLineVisible(false);
        chart.getStyler().setPlotBorderVisible(true);
        chart.getStyler().setPlotMargin(0);
        chart.getStyler().setXAxisMin(0.);
        chart.getStyler().setYAxisMin(0.);


        chart.getStyler().setChartTitleFont(new Font(Font.SANS_SERIF, Font.BOLD, 18));
        chart.getStyler().setLegendFont(new Font(Font.SERIF, Font.PLAIN, 16));
        chart.getStyler().setLegendPosition(Styler.LegendPosition.InsideSE);
        //chart.getStyler().setLegendSeriesLineLength(12);
        chart.getStyler().setAxisTitleFont(new Font(Font.SANS_SERIF, Font.ITALIC, 16));
        chart.getStyler().setAxisTickLabelsFont(new Font(Font.SERIF, Font.PLAIN, 11));
        //chart.getStyler().setDatePattern("dd-MMM");
        //chart.getStyler().setDecimalPattern("#");
        //chart.getStyler().setXAxisTickMarkSpacingHint(120);
        //chart.getStyler().setLocale(Locale.GERMAN);



        Color testSTROMA = new Color(255, 215, 0, 255);
        Color meanSTROMA = new Color(219, 112, 147, 255);


        XYSeries series1 = chart.addSeries("Stroma Population", timeVector, stromaPopulation);
        XYSeries series2 = chart.addSeries("Mean Stroma Age (x10)", timeVector, meanStromaAge);


        series1.setLineColor(testSTROMA);
        series1.setMarker(SeriesMarkers.NONE);

        series2.setLineColor(meanSTROMA);
        series2.setMarker(SeriesMarkers.NONE);


        //BitmapEncoder.saveBitmap(chart, "C:\\Users\\amymm\\OneDrive - Swansea University\\Code\\HAL-master\\MyHomeostasis", BitmapEncoder.BitmapFormat.PNG);
        //BitmapEncoder.saveBitmap(chart, "C:\\Users\\2146974\\OneDrive - Swansea University\\Code\\HAL-master\\MyHomeostasis", BitmapEncoder.BitmapFormat.PNG);

        new SwingWrapper<XYChart>(chart).displayChart();
    }
}
