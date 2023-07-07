package PDETest;

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
import org.knowm.xchart.style.lines.SeriesLines;
import org.knowm.xchart.style.markers.SeriesMarkers;

import java.awt.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.StringTokenizer;

import static HAL.Util.*;

class ExampleCell extends AgentSQ2Dunstackable<ExampleGrid> {
    public int type;
    public double IMTAge;


    public void StromaDieCell(double stromaDieProb) {
        if (G.rng.Double() < stromaDieProb) {
            //cell will die
            Dispose();
            return;
        }
    }

    public void stromaDivide(int position, int STROMA, int coordX, int coordY) {
        // reset IMT
        Dispose();
        ExampleCell resetCell = G.NewAgentSQ(coordX, coordY);
        resetCell.type = STROMA;
        // create daughter cell and assign attributes
        ExampleCell daughterCell = G.NewAgentSQ(position);
        daughterCell.type = STROMA;
    }

    public void divideCell(int stromaNCI, int STROMA, int x, int y) {
        int neighbours;
        int options;
        // change to always divide if space!
        double randomNo = G.rng.Double();
        // cell will divide
        int coordX = this.Xsq();
        int coordY = this.Ysq();
        // check for agents on boundary - assume all squares outside domain are occupied!
        // first the corners
        if (this.Xsq() == 0 && this.Ysq() == 0 || this.Xsq() == x - 1 && this.Ysq() == 0 || this.Xsq() == 0
                && this.Ysq() == y - 1 || this.Xsq() == x - 1 && this.Ysq() == y - 1) {
            neighbours = MapOccupiedHood(G.divHood) + 5;
            options = MapEmptyHood(G.divHood);
            if (options > 0 && neighbours <= 8 - stromaNCI) {
                int position = G.divHood[G.rng.Int(options)];
                stromaDivide(position, STROMA, coordX, coordY);
            } else {
                return;
            }
        } else if (this.Xsq() == 0 && this.Ysq() > 0 && this.Ysq() < y - 1 || this.Xsq() == x - 1 && this.Ysq() > 0
                && this.Ysq() < y - 1 || this.Ysq() == 0 && this.Xsq() > 0 && this.Xsq() < x - 1
                || this.Ysq() == y - 1 && this.Xsq() > 0 && this.Xsq() < x - 1) {
            // then the sides
            neighbours = MapOccupiedHood(G.divHood) + 3;
            options = MapEmptyHood(G.divHood);
            if (options > 0 && neighbours <= 8 - stromaNCI) {
                int position = G.divHood[G.rng.Int(options)];
                stromaDivide(position, STROMA, coordX, coordY);
            } else {
                return;
            }
        } else {
            // all the others
            neighbours = MapOccupiedHood(G.divHood);
            options = MapEmptyHood(G.divHood);
            if (options > 0 && neighbours <= 8 - stromaNCI) {
                int position = G.divHood[G.rng.Int(options)];
                stromaDivide(position, STROMA, coordX, coordY);
            } else {
                return;
            }
        }
    }
}


public class ExampleGrid extends AgentGrid2D<ExampleCell> {
    Rand rng = new Rand();
    int [] divHood = Util.MooreHood(false);
    PDEGrid2D drug;

    public ExampleGrid(int x, int y) {
        super(x, y, ExampleCell.class);
        drug = new PDEGrid2D(x, y);
    }

    // create initial stroma cell
    public void createInitialStromaCell(int xCord, int yCord, int STROMA, double timestepPerDay) {
        ExampleCell newCell = NewAgentSQ(xCord, yCord);
        newCell.type = STROMA;
        newCell.IMTAge = timestepPerDay*Math.random(); // use uniform distribution. Look into cell division patterns in tissues.
    }

    public void ageStromaCells(int STROMA, double timestep) {
        ShuffleAgents(rng);
        for (ExampleCell cell : this) {
            if (cell.type == STROMA) {
                double currentAge = cell.IMTAge;
                cell.IMTAge = currentAge + timestep;
            }
        }
    }

    public void StromaDieCells(double stromaDieProb, int STROMA){
        ShuffleAgents(rng);
        for (ExampleCell cell:this) {
            if (cell.type == STROMA) {
                cell.StromaDieCell(stromaDieProb);
            }
        }
    }
    public void divideCells(int stromaNCI, int STROMA, int x, int y, int VESSEL, int stromaIMT){
        ShuffleAgents(rng);
        for (ExampleCell cell:this) {
            if (cell.type != VESSEL && cell.IMTAge >= stromaIMT) {
                cell.divideCell(stromaNCI, STROMA, x, y);
            }
        }
    }

    public void initialiseDrugAValues(int VESSEL, double drugConcentrationVessel) {
        for (ExampleCell cell:this) {
            if (cell.type == VESSEL) {
                drug.Set(cell.Xsq(), cell.Ysq(), drugConcentrationVessel);
                //System.out.println(drug.Get(cell.Xsq(), cell.Ysq()));
            }
        }
    }

    public void updateDrugValues(int x, int y, int timeCounter, int treatmentTimestepDrugA, int holidayTimestepDrugA, int VESSEL, double drugDiffusionCoefficientTimestep, double drugConentrationVessel, double drugRemovalRateTimesteps) {
        boolean treatmentOn = false;
        if (timeCounter % (treatmentTimestepDrugA + holidayTimestepDrugA) <= treatmentTimestepDrugA - 1) {
            treatmentOn = true;
        }

        if (treatmentOn == true) {
            initialiseDrugAValues(VESSEL, drugConentrationVessel);
            //drug.Diffusion(drugDiffusionCoefficientTimestep);
            drug.DiffusionADI(drugDiffusionCoefficientTimestep);
        } else {
            //drug.Diffusion(drugDiffusionCoefficientTimestep);
            drug.DiffusionADI(drugDiffusionCoefficientTimestep);
        }

        // remove drug at vessel sites
        for (ExampleCell cell:this) {
            if (cell.type == VESSEL) {
                drug.Add(cell.Xsq(), cell.Ysq(), - drugRemovalRateTimesteps);
            }
        }

        drug.Update();
        for (int i = 0; i < x * y; i++) {
            if (drug.Get(i) < 0) {
                drug.Set(i, 0);
            }
        }
        drug.Update();
    }

    public void DrawModel(GridWindow win) {
        for (int i = 0; i < length; i++) {
            int color = Util.WHITE;
            ExampleCell cell=GetAgent(i);
            //System.out.println(drug.Get(i));
            if (GetAgent(i)!=null) {
                //color = cell.type;
                color = HeatMapBGR(drug.Get(i) * 1E1);
            } else{
                color = HeatMapBGR(drug.Get(i) * 1E1);
            }
            win.SetPix(i,color);
        }
    }

    public double[] PiccoSetUp(int x, int y) throws IOException {
        // Get initial stroma and vessel sites

        // Reads file and stores as string
        //FileReader file = new FileReader(new File
        // ("C:\\\\Users\\\\2146974\\\\OneDrive - Swansea University\\\\Code\\\\Noemi's code\\\\homlatt299.txt"));
        FileReader file = new FileReader(new File
                ("C:\\Users\\amymm\\Documents\\PhD\\Swansea\\Code\\Noemis code\\homlatt299.txt"));
        BufferedReader br = new BufferedReader(file);
        String temp = br.readLine();

        // Tokenise each element in string
        StringTokenizer defaultTokenizer = new StringTokenizer(temp);

        // Store each token as element in double array
        double[] setUp = new double[x * y];
        int i = 0;
        while (defaultTokenizer.hasMoreTokens()) {
            double convertedValue = Double.parseDouble(defaultTokenizer.nextToken());
            setUp[i] = convertedValue;
            i++;
        }
        return setUp;
    }

    public int[][] VesselStromaCounts(int x, int y, double [] PiccoSetUp) {
        // Determine how many vessels and stroma cells
        int vesselCn = 0;
        int stromaCn = 0;
        int [][] vesselStromaCounts = new int [1][2];
        for (int j = 0; j < x * y; j++) {
            if (PiccoSetUp[j] == 1.5) {
                stromaCn++;
            } else if (PiccoSetUp[j] == 3) {
                vesselCn++;
            }
        }
        vesselStromaCounts [0][0] = vesselCn;
        vesselStromaCounts [0][1] = stromaCn;

        return vesselStromaCounts;
    }

    public int [][] VesselLoc(int x, int y, int vesselCn, double [] setUp) {
        // Create location matrix for vessels and stroma cells
        int[][] vesselLoc = new int[vesselCn][2];
        // Store locations of initial vessels and stroma cells
        int vesselLocCn = 0;
        for (int j = 0; j < x*y; j++) {
            if (setUp[j] == 3) {
                vesselLoc[vesselLocCn][0] = Math.floorDiv(j, x);
                vesselLoc[vesselLocCn][1] = j % y;
                vesselLocCn++;
            }
        }
        return vesselLoc;
    }

    public int [][] StromaLoc(int x, int y, int stromaCn, double [] setUp) {
        // Create location matrix for vessels and stroma cells
        int[][] stromaLoc = new int[stromaCn][2];
        // Store locations of initial vessels and stroma cells
        int stromaLocCn = 0;
        for (int j = 0; j < x * y; j++) {
            if (setUp[j] == 1.5) {
                stromaLoc[stromaLocCn][0] = Math.floorDiv(j, x);
                stromaLoc[stromaLocCn][1] = j % y;
                stromaLocCn++;
            }
        }
        return stromaLoc;
    }

    public void plotStromaPopulation(double [] timeVector, double [] stromaPopulation, int countStroma, int timesteps) throws IOException {

        double [] piccoStromaPopulation = new double[timesteps];
        for (int i = 0; i < timesteps; i++) {
            piccoStromaPopulation[i] = countStroma;
        }

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
        chart.getStyler().setLegendVisible(false);


        //chart.getStyler().setAxisTickPadding(20);

        //chart.getStyler().setAxisTickMarkLength(15);

        chart.getStyler().setAxisTicksLineVisible(false);
        chart.getStyler().setPlotBorderVisible(true);
        chart.getStyler().setPlotMargin(0);
        chart.getStyler().setXAxisMin(0.);
        chart.getStyler().setYAxisMin(0.);


        chart.getStyler().setChartTitleFont(new Font(Font.SANS_SERIF, Font.BOLD, 18));
        chart.getStyler().setLegendFont(new Font(Font.SERIF, Font.PLAIN, 16));
        chart.getStyler().setLegendPosition(Styler.LegendPosition.InsideNE);
        //chart.getStyler().setLegendSeriesLineLength(12);
        chart.getStyler().setAxisTitleFont(new Font(Font.SANS_SERIF, Font.ITALIC, 16));
        chart.getStyler().setAxisTickLabelsFont(new Font(Font.SERIF, Font.PLAIN, 11));
        //chart.getStyler().setDatePattern("dd-MMM");
        //chart.getStyler().setDecimalPattern("#");
        //chart.getStyler().setXAxisTickMarkSpacingHint(120);
        //chart.getStyler().setLocale(Locale.GERMAN);



        Color testSTROMA = new Color(255, 222, 173, 255);
        Color piccoSTROMA = new Color(221, 160, 221, 255);


        XYSeries series1 = chart.addSeries("Stroma Population", timeVector, stromaPopulation);
        XYSeries series2 = chart.addSeries("Picco Stroma Population", timeVector, piccoStromaPopulation);

        series1.setLineColor(testSTROMA);
        series1.setMarker(SeriesMarkers.NONE);

        series2.setLineColor(piccoSTROMA);
        series2.setMarker(SeriesMarkers.NONE);
        series2.setLineStyle(SeriesLines.DASH_DOT);

        BitmapEncoder.saveBitmap(chart, "HomeostasisComparePicco", BitmapEncoder.BitmapFormat.PNG);

        new SwingWrapper<XYChart>(chart).displayChart();
    }


    public static void main(String[]args) throws IOException {
        int x = 300;
        int y = 300;
        int stromaNCI = 6;
        double timestep = 0.6; // (hours) length of timestep
        double timestepPerDay = 24 / timestep;
        int timesteps = 200; // each timestep is 0.6 hours that is 40 timesteps in a day!
        double stromaDieProb = 0.0001;

        double stromaInitialisationProb = 0.5449; // Picco ProbF something to do with poisson process????? was to intiate the stroma. Will speed up process.

        int STROMA = RGB256(255, 222, 173);
        int VESSEL = RGB(139,0,0);
        int CANCER = RGB(0,128,0);
        int stromaIMTDays = 1; // (days) stroma cell inter mitotic time
        int stromaIMT = stromaIMTDays * (int)timestepPerDay; // timestep stroma cell inter mitotic time
        double xDomain = 0.3; // (cm) horizontal length of sample space
        double yDomain = 0.3; // (cm) vertical length of sample space
        double deltaX = xDomain / x; // (cm) length of horizontal gridpoint
        double deltaY = yDomain / y; // (cm) length of vertical gridpoint
        double spaceConversion = (1 / deltaX) * (1 / deltaY);


        //System.out.println(spaceConversion);


        double sigmaMean = 0.016; // (cm) average distance between vessels
        double sigmaMinRaw = 0.008; // (cm) minimum distance between vessels
        //double sigmaMin = sigmaMinRaw / deltaX; // gridpoint minimum distance between vessels
        int nV = (int) ((xDomain * yDomain) / Math.pow(sigmaMean, 2)); // number of vessels in tissue
        int treatmentDaysDrugA = 1;
        int treatmentTimestepDrugA = treatmentDaysDrugA * (int)timestepPerDay;
        int holidayDaysDrugA = 1;
        int holidayTimestepDrugA = holidayDaysDrugA * (int)timestepPerDay;
        double drugDiffusionCoefficientDays = 1E-5; // (cm^2 per day) diffusion coefficient drug
        double proliferationDiffusionCoefficientDays = 1E-7; // (cm^2 per day) diffusion coefficient proliferation signal
        double drugDiffusionCoefficientTimestep = spaceConversion * drugDiffusionCoefficientDays / timestepPerDay; // (cell space per timestep) diffusion coefficient drug

        //System.out.println(drugDiffusionCoefficientTimestep);

        double proliferationDiffusionCoefficientTimestep = proliferationDiffusionCoefficientDays / timestepPerDay; // (cm^2 per timestep) diffusion coefficient proliferation signal
        double drugRemovalRateVesselDays = 10; //per day
        double autocrineProliferationSignalProductionDays = 0.65; //per day
        double paracrineProliferationSignalProductionDays = 1.54; //per day
        double proliferationDegradationDrugDays = 0.79; //per day
        double drugRemovalRateVesselTimestep = drugRemovalRateVesselDays / timestepPerDay; //per timestep
        double autocrineProliferationSignalProductionTimestep = autocrineProliferationSignalProductionDays / timestepPerDay; //per timestep
        double paracrineProliferationSignalProductionTimestep = paracrineProliferationSignalProductionDays / timestepPerDay; //per timestep
        double proliferationDegradationDrugTimestep = proliferationDegradationDrugDays / timestepPerDay; //per timestep
        double cancerProliferationSignalThresholdDeath = 0.1;
        double cancerProliferationSignalThresholdProliferation = 0.8;
        int cancerIMTDays = 1;
        int cancerIMT = cancerIMTDays *(int)timestepPerDay;
        double stromaDrugThresholdReactivity = 0.9;
        double stromaActivationNormalisationProbability = 0.01;
        double drugConcentrationVessel = 1;


        // note for cancer cells include small variation in IMT


        double [] timeVector = new double [timesteps];
        double [] stromaPopulation = new double [timesteps];


        // set up animation
        GridWindow win = new GridWindow(x, y, 2);
        ExampleGrid model = new ExampleGrid(x, y);
        GifMaker testGif = new GifMaker("testGif.gif", 10, false);



        // read Picco initial set up
        double [] piccoSetUp = model.PiccoSetUp(x, y);

        // count Picco vessels and stroma
        int [][] countsVesselStroma = model.VesselStromaCounts(x, y, piccoSetUp);
        int countVessel = countsVesselStroma[0][0];
        int countStroma = countsVesselStroma[0][1];

        // locate vessels
        int [][] LocVessels = model.VesselLoc(x, y, countVessel, piccoSetUp);

        // Initialise vessel cells
        for (int j = 0; j < countVessel; j++) {
            model.NewAgentSQ(LocVessels[j][0], LocVessels[j][1]).type = VESSEL;
        }

        // initialise PDEGrid values, here we have value of 1 at vessels sites
        model.initialiseDrugAValues(VESSEL, drugConcentrationVessel);

        // locate stroma
        int [][] LocStroma = model.StromaLoc(x, y, countStroma, piccoSetUp);

        // Initialise stroma cells

        // Fill region with stroma cells exclude vessel sites.
        for (int j = 0; j < x; j++) {
            for (int k = 0; k < y; k++) {
                int xCord = j;
                int yCord = k;
                if (model.GetAgent(xCord, yCord) != null) {
                } else {
                    double randomNo = Math.random();
                    // probability of cell existing at cite given contact inhibition
                    if (randomNo <= stromaInitialisationProb) {
                        model.createInitialStromaCell(xCord, yCord, STROMA, timestepPerDay);
                    }
                }
            }
        }



        for (int j = 0; j < timesteps; j++) {
            win.TickPause(10);

            int timeCounter = j;

            // collect stroma population data
            double timeCount = j * 1/timestepPerDay;
            double stromaCn = model.Pop() - nV;
            timeVector[j] = timeCount;
            stromaPopulation[j] = stromaCn;

            // model step cells dying
            model.StromaDieCells(stromaDieProb, STROMA);

            // age cells
            model.ageStromaCells(STROMA, timestep);

            // model step cells dividing
            model.divideCells(stromaNCI, STROMA, x, y, VESSEL, stromaIMT);

            // update PDEGrid values
            model.updateDrugValues(x, y, timeCounter, treatmentTimestepDrugA, holidayTimestepDrugA, VESSEL, drugDiffusionCoefficientTimestep, drugConcentrationVessel, drugRemovalRateVesselTimestep);

            // increment timestep
            model.IncTick();

            // draw
            model.DrawModel(win);
            testGif.AddFrame(win);

        }

        testGif.Close();

        // plot stroma population
        //model.plotStromaPopulation(timeVector, stromaPopulation, countStroma, timesteps);


        // output number of stroma in initialised sets
        //int initStromaCn = model.Pop() - nV;
        //System.out.println("Number of stroma in Picco initial set up : " + countStroma);
        //System.out.println("Number of stroma in Milne initial set up : " + initStromaCn);
        System.out.println("Simulation finished!");


    }
}
