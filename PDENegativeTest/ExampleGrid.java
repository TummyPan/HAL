package PDENegativeTest;

import HAL.GridsAndAgents.AgentGrid2D;
import HAL.GridsAndAgents.AgentSQ2Dunstackable;
import HAL.GridsAndAgents.PDEGrid2D;
import HAL.Gui.GifMaker;
import HAL.Gui.GridWindow;
import HAL.Rand;
import HAL.Util;
import org.knowm.xchart.BitmapEncoder;
import org.knowm.xchart.SwingWrapper;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYSeries;
import org.knowm.xchart.style.Styler;
import org.knowm.xchart.style.markers.SeriesMarkers;

import java.awt.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.StringTokenizer;

import static HAL.Util.HeatMapBGR;
import static HAL.Util.RGB256;

class ExampleCell extends AgentSQ2Dunstackable<ExampleGrid> {
    // declare individual cell attributes
    public int type;
    public double IMTAge;
    public double cancerVariedIMT;


    public void stromaDieCell(double stromaDieProb) {
        // turnover of stroma cells with given probability
        if (G.rng.Double() < stromaDieProb) {
            //cell will die
            Dispose();
            return;
        }
    }

    public void stromaDivide(int position, int coordX, int coordY) {
        // dividing cell dies and new cell of same type is created in same position
        int cellType = this.type;
        Dispose();
        ExampleCell resetCell = G.NewAgentSQ(coordX, coordY);
        resetCell.type = cellType;
        resetCell.IMTAge = 0;
        // create daughter cell and assign attributes
        ExampleCell daughterCell = G.NewAgentSQ(position);
        daughterCell.type = cellType;
        daughterCell.IMTAge = 0;
    }

    public void divideStromaCell(int stromaNCI, int x, int y) {
        // get x and y coordinates of cell
        int coordX = this.Xsq();
        int coordY = this.Ysq();
        int options = MapEmptyHood(G.divHood);
        if (options > 0 && options >= stromaNCI) {
            int position = G.divHood[G.rng.Int(options)];
            //System.out.println("position: " + position);
            stromaDivide(position, coordX, coordY);
        }
    }

    public boolean seedCancer(int CANCER, boolean cancerSeeded, double cancerIMT, double cancerIMTVariation
            , double proliferationCancerSite) {
        // check neighbourhood of vessel cell for empty positions
        int vesselEmptyNeighbours = MapEmptyHood(G.divHood);
        // seeds cancer cell in random empty neighbour of vessel
        if (vesselEmptyNeighbours > 0) {
            int position = G.divHood[G.rng.Int(vesselEmptyNeighbours)];
            ExampleCell cancerSeed = G.NewAgentSQ(position);
            cancerSeed.type = CANCER;
            cancerSeed.IMTAge = 0;
            // random variation of each cancer cells IMT
            cancerSeed.cancerVariedIMT = cancerIMT + (2*Math.random() - 1)*cancerIMTVariation;
            G.proliferation.Add(cancerSeed.Xsq(), cancerSeed.Ysq(), proliferationCancerSite);
            G.proliferation.Update();
            cancerSeeded = true;
        }
        return cancerSeeded;
    }

    public boolean checkCancerNeighbour(int x, int y, int CANCER) {
        boolean cancerPresent = false;
        int [] neighbourhood = G.divHood;
        int hoodsize = MapOccupiedHood(neighbourhood);
        for (int i = 0; i < hoodsize; i++) {
            ExampleCell checkCell = G.GetAgent(neighbourhood[i]);
            if (checkCell.type == CANCER) {
                cancerPresent = true;
            }
        }
        return cancerPresent;
    }

    public void cancerDieCell(int CANCER) {
        //cell will die
        Dispose();
        System.out.println("cancer dies");
        G.checkElimination(CANCER);
        return;
    }

    public void cancerDivide(int position, int coordX, int coordY, int cancerIMT, double cancerIMTVariation, double proliferationCancerSite) {
        // dividing cell dies and new cell of same type is created in same position
        int cellType = this.type;
        Dispose();
        ExampleCell resetCell = G.NewAgentSQ(coordX, coordY);
        resetCell.type = cellType;
        resetCell.IMTAge = 0;
        // random variation of each cancer cells IMT
        resetCell.cancerVariedIMT =  cancerIMT + (2*Math.random() - 1)*cancerIMTVariation;
        // create daughter cell and assign attributes
        ExampleCell daughterCell = G.NewAgentSQ(position);
        daughterCell.type = cellType;
        daughterCell.IMTAge = 0;
        // random variation of each cancer cells IMT
        daughterCell.cancerVariedIMT =  cancerIMT + (2*Math.random() - 1)*cancerIMTVariation;
        // update proliferation signal PDEGrid
        G.proliferation.Add(daughterCell.Xsq(), daughterCell.Ysq(), proliferationCancerSite);
        G.proliferation.Update();
    }

    public void divideCancerCell(int x, int y, int cancerIMT, double cancerIMTVariation, double proliferationCancerSite) {
        // declare neighbourhood integers
        int options;
        // get x and y coordinates of cell
        int coordX = this.Xsq();
        int coordY = this.Ysq();
        options = MapEmptyHood(G.divHood);
        if (options > 0) {
            int position = G.divHood[G.rng.Int(options)];
            cancerDivide(position, coordX, coordY, cancerIMT, cancerIMTVariation, proliferationCancerSite);
            //System.out.println("cancer divides");
        }
    }
}

public class ExampleGrid extends AgentGrid2D<ExampleCell> {
    // random number
    Rand rng = new Rand();
    // neighbourhood function
    int [] divHood = Util.MooreHood(false);
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

    public void ageStromaCells(int STROMA, int ACTIVATEDSTROMA) {
        ShuffleAgents(rng);
        // age stroma cell
        for (ExampleCell cell : this) {
            if ((cell.type == STROMA || cell.type == ACTIVATEDSTROMA)) {
                double currentAge = cell.IMTAge;
                cell.IMTAge = currentAge + 1;
            }
        }
    }

    public void stromaDieCells(double stromaDieProb, int STROMA, int ACTIVATEDSTROMA){
        ShuffleAgents(rng);
        // check for stroma cell type and turnover
        for (ExampleCell cell:this) {
            if ((cell.type == STROMA || cell.type == ACTIVATEDSTROMA)) {
                cell.stromaDieCell(stromaDieProb);
            }
        }
    }

    public void divideStromaCells(int stromaNCI, int STROMA, int x, int y, int ACTIVATEDSTROMA, int stromaIMT){
        ShuffleAgents(rng);
        // check for stroma cell and its IMTAge and divide
        for (ExampleCell cell:this) {
            if ((cell.type == STROMA || cell.type == ACTIVATEDSTROMA) && cell.IMTAge >= stromaIMT) {
                cell.divideStromaCell(stromaNCI, x, y);
            }
        }
    }

    public void stromaStatusUpdate(int x, int y, int STROMA, int ACTIVATEDSTROMA, double stromaDrugThresholdReactivity, int CANCER) {
        ShuffleAgents(rng);
        // check for stroma cell neighbouring cancer cells and drug threshold
        for (ExampleCell cell:this) {
            boolean cancerNeighbour = false;
            // activate stroma if it has cancer neighbour and drug is equal or above threshold
            if (cell.type == STROMA) {
                double drugAmount = drug.Get(cell.Isq());
                cancerNeighbour = cell.checkCancerNeighbour(x, y, CANCER);
                if (cancerNeighbour == true && drugAmount >= stromaDrugThresholdReactivity) {
                    cell.type = ACTIVATEDSTROMA;
                    //System.out.println("stroma activated!");
                }
            // normalise stroma if drug is under threshold
            } else if (cell.type == ACTIVATEDSTROMA) {
                double drugAmount = drug.Get(cell.Isq());
                if (drugAmount < stromaDrugThresholdReactivity) {
                    cell.type = STROMA;
                }
            }
        }
    }

    public void initialiseDrugAValues(int VESSEL, double drugConcentrationVessel) {
        // initialise drug from vessel sites
        for (ExampleCell cell:this) {
            if (cell.type == VESSEL) {
                drug.Set(cell.Xsq(), cell.Ysq(), drugConcentrationVessel);
            }
        }
        drug.Update();
    }

    public void updatePDEValues(int x, int y, int timeCounter, int treatmentTimestepDrugA, int holidayTimestepDrugA,
                                 int VESSEL, double drugDiffusionCoefficientTimestep, double drugConcentrationVessel,
                                 double drugRemovalRateVesselTimestep, double proliferationDiffusionCoefficientTimestep,
                                 int CANCER, double autocrineProliferationSignalProductionTimestep, int ACTIVATEDSTROMA,
                                 double paracrineProliferationSignalProductionTimestep,
                                 double proliferationDegradationDrugTimestep, boolean cancerDetected) {
        // start drug treatment if cancer has reached detectable size
        if (cancerDetected == false){
            // no treatment
        } else {
            // check if treatment is on
            boolean treatmentOn = false;
            if (timeCounter % (treatmentTimestepDrugA + holidayTimestepDrugA) <= treatmentTimestepDrugA - 1) {
                treatmentOn = true;
            }
            // apply diffusion to drug according to treatment schedule
            if (treatmentOn == true) {
                initialiseDrugAValues(VESSEL, drugConcentrationVessel);
                //drug.Diffusion(drugDiffusionCoefficientTimestep);
                drug.DiffusionADI(drugDiffusionCoefficientTimestep);
            } else {
                //drug.Diffusion(drugDiffusionCoefficientTimestep);
                drug.DiffusionADI(drugDiffusionCoefficientTimestep);
            }
            // remove drug at vessel sites
            for (ExampleCell cell : this) {
                if (cell.type == VESSEL) {
                    drug.Add(cell.Xsq(), cell.Ysq(), -drugRemovalRateVesselTimestep);
                }
            }
        }
        // apply diffusion to proliferation signal
        //proliferation.Diffusion(proliferationDiffusionCoefficientTimestep);
        proliferation.DiffusionADI(proliferationDiffusionCoefficientTimestep);
        // autocrine increase
        for (ExampleCell cell:this) {
            if (cell.type == CANCER) {
                proliferation.Add(cell.Xsq(), cell.Ysq(), autocrineProliferationSignalProductionTimestep);
            }
        }
        // paracrine increase
        for (ExampleCell cell:this) {
            if (cell.type == ACTIVATEDSTROMA) {
                proliferation.Add(cell.Xsq(), cell.Ysq(), paracrineProliferationSignalProductionTimestep);
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
                if (GetAgent(i).type == VESSEL) {
                    //System.out.println("vessel site!");
                } else {
                    System.out.println("negative value at non vessel site: " + i);
                }
                drug.Set(i, 0);
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

    public void initialiseCancer(int x, int y, int nV, int [][] LocVessels, int CANCER, int cancerIMT,
                                 double cancerIMTVariation, double proliferationCancerSite) {
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
                        cancerSeeded = centreVesselCell.seedCancer(CANCER, cancerSeeded, cancerIMT, cancerIMTVariation
                                , proliferationCancerSite);
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

    public void cancerDieProliferationUpdate(int CANCER, double cancerProliferationSignalThresholdDeath,
                                             double cancerProliferationSignalThresholdProliferation) {
        for (ExampleCell cell:this) {
            if (cell.type == CANCER) {
                double proliferationAmount = proliferation.Get(cell.Isq());
                if (proliferationAmount < cancerProliferationSignalThresholdDeath) {
                    // cancer cell dies
                    cell.cancerDieCell(CANCER);
                } else if (proliferationAmount < cancerProliferationSignalThresholdProliferation) {
                    // cancer cell is senescent
                } else if (proliferationAmount >= cancerProliferationSignalThresholdProliferation) {
                    // age cancer cells that are not senescent
                    double currentAge = cell.IMTAge;
                    cell.IMTAge = currentAge + 1;
                }
            }
        }
    }

    public void divideCancerCells(int x, int y, int CANCER, int cancerIMT, double cancerIMTVariation, double proliferationCancerSite) {
        ShuffleAgents(rng);
        // check for cancer cell and its cancerVariedIMT and divide
        for (ExampleCell cell:this) {
            if (cell.type == CANCER && cell.IMTAge >= cell.cancerVariedIMT) {
                cell.divideCancerCell(x, y, cancerIMT, cancerIMTVariation, proliferationCancerSite);
            }
        }
    }

    public void checkElimination(int CANCER) {
        int cancerCn = 0;
        for (ExampleCell cell:this) {
            if (cell.type == CANCER) {
                cancerCn++;
            }
        }
        if (cancerCn == 0) {
            System.out.println("cancer eliminated!");
        }
    }

    public void drawModel(GridWindow win) {
        for (int i = 0; i < length; i++) {
            // set background colour
            int color = Util.WHITE;
            // colour each agent accoring to type and any vacant cell will take PDE grid colour
            ExampleCell cell=GetAgent(i);
            if (GetAgent(i)!=null) {
                //color = cell.type;
                color = HeatMapBGR(drug.Get(i) * 1E1); // turn on for drug
                //color = HeatMapRBG(proliferation.Get(i) * 1E1); // turn on for proliferation signal
            } else{
                color = HeatMapBGR(drug.Get(i) * 1E1); // turn on for drug
                //color = HeatMapRBG(proliferation.Get(i) * 1E1); // turn on for proliferation signal
            }
            // set up display
            win.SetPix(i,color);
        }
    }

    public int countHomeostaticStroma() throws IOException {
        // Get initial stroma count from homeostatic run
        // Reads file and stores as string
        FileReader file = new FileReader(new File
                ("C:\\Users\\amymm\\OneDrive\\Documents\\PhD\\Swansea\\Code\\Results\\stromaHomeostasis" +
                        "\\100KstromaHomeostasisCount.csv"));
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


    public double[][] homeostaticStromaSetUp(int homestaticStromaCount) throws IOException {
        // Get initial stroma sites from homeostasis run
        // Reads file and stores as string
        FileReader file = new FileReader(new File
                ("C:\\Users\\amymm\\OneDrive\\Documents\\PhD\\Swansea\\Code\\Results\\stromaHomeostasis" +
                        "\\100KstromaHomeostasis.csv"));
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

    public void initialiseStroma(double[][] homeostaticStroma, int STROMA) {
        for (int i = 0; i < homeostaticStroma.length; i++) {
            // create initial stroma cell and assign attributes
            ExampleCell newCell = NewAgentSQ((int)homeostaticStroma[i][0]);
            newCell.type = STROMA;
            newCell.IMTAge = homeostaticStroma[i][1];
        }
    }

    public int [][] readVesselInitialisation(int nV) throws IOException {
        FileReader file = new FileReader(new File
                ("C:\\Users\\amymm\\OneDrive\\Documents\\PhD\\Swansea\\Code\\Results\\vesselLocations\\vesselLocations.csv"));
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

    public double [] totalStromaPopulation(double[] stromaPopulation, double[] activatedStromaPopulation) {
        double[] total = new double[stromaPopulation.length];
        for (int i = 0; i < total.length; i++) {
            total[i] = stromaPopulation[i] + activatedStromaPopulation[i];
        }
        return total;
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

    public void plotPopulations(double [] timeVector, double [] stromaPopulation, double [] cancerPopulation, double [] activatedStromaPopulation, double [] totalStromaPopulation, int timesteps)
            throws IOException {

        // declare and create plot
        XYChart chart = new XYChart(500, 400);
        // chart attributes
        chart.setTitle("Populations"); // set plot title
        chart.setXAxisTitle("time(days)"); // set plot x axis title
        chart.setYAxisTitle("Populations"); // set plot y axis title
        //chart.getStyler().setPlotBackgroundColor(ChartColor.getAWTColor(ChartColor.GREY));
        //chart.getStyler().setPlotGridLinesColor(new Color(255, 255, 255));
        chart.getStyler().setChartBackgroundColor(Color.WHITE);  // set plot background colour
        //chart.getStyler().setLegendBackgroundColor(Color.WHITE);
        //chart.getStyler().setChartFontColor(Color.RED);
        //chart.getStyler().setChartTitleBoxBackgroundColor(new Color(0, 222, 0));
        chart.getStyler().setChartTitleBoxVisible(false); // boolean for title box visible
        chart.getStyler().setChartTitleBoxBorderColor(Color.BLACK); // set plot border colour
        chart.getStyler().setPlotGridLinesVisible(false); // boolean for grid lines visible
        chart.getStyler().setLegendVisible(true); // boolean for legend visible
        //chart.getStyler().setAxisTickPadding(20);
        //chart.getStyler().setAxisTickMarkLength(15);
        chart.getStyler().setAxisTicksLineVisible(false); // boolean for tick lines visible
        chart.getStyler().setPlotBorderVisible(true); // boolean for border visible
        chart.getStyler().setPlotMargin(0); // set plot margin
        chart.getStyler().setXAxisMin(0.); // set plot x axis minimum
        chart.getStyler().setYAxisMin(0.); // set plot y axis minimum
        chart.getStyler().setChartTitleFont(new Font(Font.SANS_SERIF, Font.BOLD, 18)); // set plot title font
        chart.getStyler().setLegendFont(new Font(Font.SERIF, Font.PLAIN, 16)); // set legend font
        chart.getStyler().setLegendPosition(Styler.LegendPosition.OutsideE); // set legend position
        //chart.getStyler().setLegendSeriesLineLength(12);
        chart.getStyler().setAxisTitleFont(new Font(Font.SANS_SERIF, Font.ITALIC, 16)); // set plot axes title font
        chart.getStyler().setAxisTickLabelsFont(new Font(Font.SERIF, Font.PLAIN, 11)); // set plot axes tick font
        //chart.getStyler().setDatePattern("dd-MMM");
        //chart.getStyler().setDecimalPattern("#");
        //chart.getStyler().setXAxisTickMarkSpacingHint(120);
        //chart.getStyler().setLocale(Locale.GERMAN);
        Color STROMA = new Color(255, 222, 173, 255); // set colour for stroma population
        Color ACTIVATEDSTROMA = new Color(173, 216, 230, 255); // set colour for activated stroma population
        Color TOTALSTROMA = new Color(144, 238, 144, 175); // set colour for total stroma population
        Color CANCER = new Color(138, 43, 226, 255); // set colour for cancer population
        // add series to plot
        XYSeries series1 = chart.addSeries("Stroma", timeVector, stromaPopulation);
        XYSeries series2 = chart.addSeries("Cancer", timeVector, cancerPopulation);
        XYSeries series3 = chart.addSeries("Activated Stroma", timeVector, activatedStromaPopulation);
        XYSeries series4 = chart.addSeries("Total Stroma", timeVector, totalStromaPopulation);
        // set stroma population line colour and shape
        series1.setLineColor(STROMA);
        series1.setMarker(SeriesMarkers.NONE);
        // set cancer population line colour and shape
        series2.setLineColor(CANCER);
        series2.setMarker(SeriesMarkers.NONE);
        // set activated stroma population line colour and shape
        series3.setLineColor(ACTIVATEDSTROMA);
        series3.setMarker(SeriesMarkers.NONE);
        // set total stroma population line colour and shape
        series4.setLineColor(TOTALSTROMA);
        series4.setMarker(SeriesMarkers.NONE);
        // save plot as png
        BitmapEncoder.saveBitmap(chart, "Populations", BitmapEncoder.BitmapFormat.PNG);
        // display plot
        new SwingWrapper<XYChart>(chart).displayChart();
    }

    public static void main(String[]args) throws IOException {

        // log start time of algorithm
        long startTime = System.nanoTime();

        // declare parameters

        // grid parameters
        int x = 300; // number of cells horizontal
        int y = 300; // number of cells vertical
        double xDomain = 0.3; // (cm) horizontal length of sample space
        double yDomain = 0.3; // (cm) vertical length of sample space
        double deltaX = xDomain / x; // (cm) length of horizontal gridpoint
        double deltaY = yDomain / y; // (cm) length of vertical gridpoint
        double spaceConversion = (1 / deltaX) * (1 / deltaY); // space conversion parameter (used for PDEGrids)



        // number of iterations
        int timesteps = 500; // number of iterations of model

        // time parameters
        double timestepHour = 0.6; // (hours) length of timestep
        // number of timesteps in a day (each timestep is 0.6 hours that is 40 timesteps in a day!)
        double timestepPerDay = 24 / timestepHour; // number of timesteps in day
        double timestepDay = 1 / timestepPerDay; // (days) length of timestep
        double noDays = timesteps / timestepPerDay;

        //System.out.println("Number of days: " + noDays);


        // agent paramenters
        int STROMA = RGB256(255, 222, 173); // set colour for stroma cells
        int VESSEL = RGB256(128,0,0); // set colour for vessel cells
        int CANCER = RGB256(138,43,226); // set colour for cancer cells
        int ACTIVATEDSTROMA = RGB256(173,216,230); // set colour for activated stroma cells

        // vessel parameters
        double sigmaMean = 0.016; // (cm) average distance between vessels
        double sigmaMinRaw = 0.008; // (cm) minimum distance between vessels
        double sigmaMin = sigmaMinRaw / deltaX; // gridpoint minimum distance between vessels
        int nV = (int) ((xDomain * yDomain) / Math.pow(sigmaMean, 2)); // number of vessels in tissue

        // stroma parameters
        int stromaNCI = 6; // stroma contact inhibition (number of empty neighbour cells required for division)
        double stromaDieProb = 0.0001; // stroma turnover probability
        double stromaInitialisationProb = 0.45;//0.5449; // stroma initialisation probability (Picco)
        int stromaIMTDays = 1; // (days) stroma cell inter mitotic time
        int stromaIMT = stromaIMTDays * (int)timestepPerDay; // (timesteps) stroma cell inter mitotic time
        double stromaDrugThresholdReactivity = 0.1;//0.9; // stroma drug threshold reactivity
        double stromaActivationNormalisationProbability = 0.01; // stroma activation normalisation probability

        // cancer parameters
        double cancerDetectionThreshold = 1E4; // (Picco) number of cancer cells required for detection
        //System.out.println("cancer detection threshold: " + cancerDetectionThreshold);
        boolean cancerDetected = true;//false; // cancer has reached detectable size
        double cancerProliferationSignalThresholdDeath = 0.1; // cancer proliferation signal threshold death
        // cancer proliferation signal threshold proliferation
        double cancerProliferationSignalThresholdProliferation = 0.8;
        int cancerIMTDays = 1; // (days) cancer cell inter mitotic time
        int cancerIMT = cancerIMTDays *(int)timestepPerDay; // (timesteps) cancer cell inter mitotic time
        double cancerIMTVariation = 0.1*timestepPerDay; // variation of cancer cell IMT

        // drug parameters
        double drugDiffusionCoefficientDays = 1E-5; // (cm^2 per day) diffusion coefficient drug
        // (cm^2 per day) diffusion coefficient proliferation signal
        double drugDiffusionCoefficientTimestep = spaceConversion * drugDiffusionCoefficientDays / timestepPerDay;
        // (cm^2 per timestep) diffusion coefficient proliferation signal
        double drugRemovalRateVesselDays = 10; // (days) removal rate of drugA at vessel site
        double drugRemovalRateVesselTimestep = drugRemovalRateVesselDays / timestepPerDay;
        double drugConcentrationVessel = 1; // drug concentration on delivery

        //double diffCoef = drugDiffusionCoefficientDays * timestepDay / Math.pow(deltaX, 2); // non-dimensionalised diffusion coefficient

        //System.out.println("timesteps per day: " + timestepPerDay);
        //System.out.println("drugDiffusionCoefficientTimestep: " + drugDiffusionCoefficientTimestep);

        //System.out.println("diffCoef: " + diffCoef);


        // proliferation signal parameters
        double proliferationDiffusionCoefficientDays = 1E-7;
        // (cell space per timestep) diffusion coefficient drugA
        double proliferationDiffusionCoefficientTimestep = proliferationDiffusionCoefficientDays / timestepPerDay;
        double autocrineProliferationSignalProductionDays = 0.65; // (days) autocrine proliferation signal rate
        double paracrineProliferationSignalProductionDays = 1.54; // (days) paracrine proliferation signal rate
        double proliferationDegradationDrugDays = 0.79; // (days) proliferation signal degradation rate
        // (timesteps) removal rate of drugA at vessel site
        double autocrineProliferationSignalProductionTimestep = autocrineProliferationSignalProductionDays
                / timestepPerDay; // (timesteps) autocrine proliferation signal rate
        double paracrineProliferationSignalProductionTimestep = paracrineProliferationSignalProductionDays
                / timestepPerDay; // (timesteps) paracrine proliferation signal rate
        // (timesteps) proliferation signal degradation rate
        double proliferationDegradationDrugTimestep = proliferationDegradationDrugDays / timestepPerDay;
        double proliferationCancerSite = 1;//autocrineProliferationSignalProductionDays; // proliferation signal at initial cancer cell site


        // treatment parameters
        int treatmentDaysDrugA = 1; // number of treatment days for drugA
        // number of treatment timesteps for drugA
        int treatmentTimestepDrugA = treatmentDaysDrugA * (int)timestepPerDay;
        int holidayDaysDrugA = 5; // number of holiday days for drugA
        int holidayTimestepDrugA = holidayDaysDrugA * (int)timestepPerDay; // number of holiday timesteps for drugA



        // declare time and population storage arrays
        double [] timeVector = new double [timesteps];
        double [] stromaPopulation = new double [timesteps];
        double [] cancerPopulation = new double[timesteps];
        double [] activatedStromaPopulation = new double[timesteps];

        // set up animation
        GridWindow win = new GridWindow(x, y, 1);
        ExampleGrid model = new ExampleGrid(x, y);
        // declare animation as gif
        GifMaker testGif = new GifMaker("Animation.gif", 10, false);

        // read vessel initialisation
        int [][] vesselSetUp = model.readVesselInitialisation(nV);

        // initialise vessel cells
        model.initialiseVesselCells(vesselSetUp, nV, VESSEL);

        // read homeostatic stroma set up
        // number of initial stroma cells from homeostatic run
        int homestaticStromaCount = model.countHomeostaticStroma();
        // determine position and age of initial stroma population
        double [][] homeostaticStroma = model.homeostaticStromaSetUp(homestaticStromaCount);

        // initialise stroma cells
        model.initialiseStroma(homeostaticStroma, STROMA);

        // initialise cancer
        //model.initialiseCancer(x, y, nV, vesselSetUp, CANCER, cancerIMT, cancerIMTVariation, proliferationCancerSite);

        model.drawModel(win);

        // main loop
        for (int j = 0; j < timesteps; j++) {
            // time between animation frames
            win.TickPause(1);

            // create time counter
            int timeCounter = j;


            // collect population data
            double timeCount = j * 1/timestepPerDay;
            double stromaCn = model.populationType(STROMA);
            double cancerCn = model.populationType(CANCER);
            double activatedStromaCn = model.populationType(ACTIVATEDSTROMA);
            timeVector[j] = timeCount;
            stromaPopulation[j] = stromaCn;
            cancerPopulation[j] = cancerCn;
            activatedStromaPopulation[j] = activatedStromaCn;

            // print out message every thousand iterations
            if (j % 1000 == 0) {
                System.out.println("iteration: " + j);
                //System.out.println("cancer population: " + cancerCn);
            }

            // check if cancer cell population has reached detection threshold
            /*if (cancerDetected == false) {
                if (cancerCn >= cancerDetectionThreshold) {
                    cancerDetected = true;
                    System.out.println("cancer detected iteration: " + j);
                }
            }*/

            // stroma cells dying
            model.stromaDieCells(stromaDieProb, STROMA, ACTIVATEDSTROMA);

            // stroma cells activated or normalised
            /*if (cancerDetected == true) {
                model.stromaStatusUpdate(x, y, STROMA, ACTIVATEDSTROMA, stromaDrugThresholdReactivity, CANCER);
            }*/

            // cancer cell proliferation update and age cancer cells
            /*model.cancerDieProliferationUpdate(CANCER, cancerProliferationSignalThresholdDeath,
                    cancerProliferationSignalThresholdProliferation);*/

            // age stroma cells
            model.ageStromaCells(STROMA, ACTIVATEDSTROMA);

            // stroma cells dividing
            model.divideStromaCells(stromaNCI, STROMA, x, y, VESSEL, stromaIMT);

            // cancer cells dividing
            //model.divideCancerCells(x, y, CANCER, cancerIMT, cancerIMTVariation, proliferationCancerSite);

            // update PDEGrid values
            //System.out.println("in main : " + j);
            model.updatePDEValues(x, y, timeCounter, treatmentTimestepDrugA, holidayTimestepDrugA, VESSEL,
                    drugDiffusionCoefficientTimestep, drugConcentrationVessel, drugRemovalRateVesselTimestep,
                    proliferationDiffusionCoefficientTimestep, CANCER, autocrineProliferationSignalProductionTimestep,
                    ACTIVATEDSTROMA, paracrineProliferationSignalProductionTimestep,
                    proliferationDegradationDrugTimestep, cancerDetected);

            // increment timestep
            model.IncTick();

            // draw
            model.drawModel(win);
            // store frame in gif
            testGif.AddFrame(win);

            // save final set up as png
            if (j % 20 == 0) {
                win.ToPNG("drug" + j + ".png");
            }

        }

        //create total stroma population array
        double [] totalStromaPopulation = model.totalStromaPopulation(stromaPopulation, activatedStromaPopulation);
        // end gif storage
        testGif.Close();

        // plot populations
        model.plotPopulations(timeVector, stromaPopulation, cancerPopulation, activatedStromaPopulation, totalStromaPopulation, timesteps);

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
