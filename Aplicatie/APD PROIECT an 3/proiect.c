#include <ansi_c.h>
#include <cvirte.h>		
#include <userint.h>
#include <utility.h>
#include<formatio.h>
#include "proiect.h"
#include <analysis.h>
#include <string.h>
#include <math.h>


#include "histograma.h"

#define SAMPLE_RATE		0
#define NPOINTS			1
#define GRAPH_NPOINTS 					500
#define ALPHA            0.1
#define WINDOW_SIZE      100
#define M_PI 3.14

//static double rowSignal[ GRAPH_NPOINTS+1 ];
static int panelHandle;
static int histogramPanel;
static int freqPanel;
int soundHandle = -1;
int main (int argc, char *argv[])
{
	if (InitCVIRTE (0, argv, 0) == 0)
		return -1;	/* out of memory */
	if ((panelHandle = LoadPanel (0, "proiect.uir", PANEL)) < 0)
		return -1;
	if ((histogramPanel = LoadPanel (0, "histograma.uir", PANEL_HIST)) < 0)
		return -1;
		if ((freqPanel = LoadPanel (0, "proiect.uir", PANEL_FREQ)) < 0)
		return -1;
	
    
	DisplayPanel (panelHandle);
	DisplayPanel(histogramPanel);
	
	RunUserInterface ();
	DiscardPanel (panelHandle);
	DiscardPanel(histogramPanel);
	return 0;
}

int CVICALLBACK OnExit (int panel, int control, int event,
						void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			if(ConfirmPopup ("Quit", "Do you Really want to quit ?"))
			{
				QuitUserInterface (0);
			}
			
			break;
	}
	return 0;
}




//==============================================================================
// Global variables
int waveInfo[2]; //waveInfo[0] = sampleRate
				 //waveInfo[1] = number of elements

static int signalLength = 4096;      // Lungimea semnalului
static int fftWindowSize = 1024;     // Dimensiunea FFT (valoare implicit?)
static int currentStart = 0;  

double sampleRate = 0.0;
int npoints = 0;
	
double *waveData = 0;
double *convertedSpectrum  = 0;
int windowType = 0;
int filterSpectre = 0;

int fftSize = 1024;
int timerID = 0;

int line1 = -1;
int line2 = -1;
int startPoint  =0;


double windowParameter;
WindowConst* windowConstants;

double dt;
double autoSpectrum[];
double df;
double searchFrequency;
int DFT = 1024;

double *outputSpectre = 0;

ssize_t frequncySpan;
double frequencyPeak;
double powerPeak;

int currentSegmentStart = 0; 
int segmentSize = 0; 
//segmentSize = (int)sampleRate;


int CVICALLBACK OnLoadButtonCB (int panel, int control, int event,
								void *callbackData, int eventData1, int eventData2)
{
	
	switch (event)
	{
		case EVENT_COMMIT:
			FileToArray ("wafeInfo.txt", waveInfo, VAL_INTEGER, 2, 1, VAL_GROUPS_TOGETHER, VAL_GROUPS_AS_COLUMNS, VAL_ASCII);
			
			sampleRate = waveInfo[SAMPLE_RATE];
			
			segmentSize = (int)sampleRate; 
			
			npoints = waveInfo[NPOINTS];
			
			//alocare memorie pentru numarul de puncte
			waveData = (double *) calloc(npoints, sizeof(double));
			
			//incarcare din fisierul .txt in memorie (vector)
			FileToArray("waveData.txt", waveData, VAL_DOUBLE, npoints, 1, VAL_GROUPS_TOGETHER, VAL_GROUPS_AS_COLUMNS, VAL_ASCII);
			
			if (npoints > (10 * sampleRate)) { 
				int newNPoints = 6 * sampleRate; 
				double *tempData = (double *)calloc(newNPoints, sizeof(double));
				memcpy(tempData, waveData, newNPoints * sizeof(double));
				free(waveData);
				waveData = tempData;
				npoints = newNPoints; 
			}
			
			PlotY (panel, PANEL_MAIN_PANEL_IDC_GRAPH_, waveData, sampleRate, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, 1, VAL_RED);
			double skewness;
			Moment(waveData, 256, 4, &skewness);
			double kurtosis;
			Moment(waveData, 256, 3, &kurtosis);
			double maxVal = 0.0;
			double minVal = 0.0;
			int maxIndex = 0;
			int minIndex = 0;
			double mean = 0.0;
			double dispersia = 0.0;
	       	MaxMin1D(waveData,sampleRate,&maxVal,&maxIndex,&minVal,&minIndex);
			Mean(waveData, sampleRate, &mean);
	        Variance(waveData, sampleRate, &mean, &dispersia);
	        
			int zeroCrossings = 0;
			for (int i = 1; i < sampleRate; i++) {
				if ((waveData[i-1] < 0 && waveData[i] > 0) || (waveData[i-1] > 0 && waveData[i] < 0)) {
					zeroCrossings++;
				}
			}
	        
			
	        SetCtrlVal(panelHandle, PANEL_NUMERIC_3, minVal);
	        SetCtrlVal(panelHandle, PANEL_NUMERIC_2, maxVal);
	        SetCtrlVal(panelHandle, PANEL_NUMERIC, mean);
			SetCtrlVal(panelHandle, PANEL_NUMERIC_4, dispersia);
			SetCtrlVal(panelHandle, PANEL_NUMERIC_5, zeroCrossings);
			SetCtrlVal(panelHandle, PANEL_NUMERIC_6, skewness);
			SetCtrlVal(panelHandle,	PANEL_NUMERIC_7, kurtosis);
			SetCtrlVal(panelHandle, PANEL_NUMERIC_8, minIndex);
			SetCtrlVal(panelHandle, PANEL_NUMERIC_9, maxIndex);
			
			/*ScaledWindowEx(waveData, npoints, RECTANGLE_, windowParameter,windowConstants);
			AutoPowerSpectrum(waveData, npoints, dt, autoSpectrum, df);
			
			PowerFrequencyEstimate(autoSpectrum, npoints, searchFrequency, windowConstants, df, frequncySpan, frequencyPeak);
			
			SpectrumUnitConversion(autoSpectrum, npoints, RECTANGLE_, SCALING_MODE_LINEAR, DISPLAY_UNIT_VRMS, df, windowConstants, );*/
			
			int N = 0;
			
			WindowConst winConst;
            ScaledWindowEx(waveData, npoints, RECTANGLE_, 0.0, &winConst);
            
            // Calculare spectru de putere
            double *autoSpectrum = (double *) calloc(npoints / 2, sizeof(double));
            double df;
            dt = 1.0 / sampleRate;
            AutoPowerSpectrum(waveData, npoints, dt, autoSpectrum, &df);
            
            // Estimare frecven?? dominant?
            
            PowerFrequencyEstimate(autoSpectrum, npoints / 2, 0.0, winConst, df, 0, &frequencyPeak, &powerPeak);
			
            
            // Conversie spectru
            convertedSpectrum = (double *) calloc(npoints / 2, sizeof(double));
            char unitString[32];
            SpectrumUnitConversion(autoSpectrum, npoints / 2, RECTANGLE_, SCALING_MODE_LINEAR, DISPLAY_UNIT_VRMS, df, winConst, convertedSpectrum, unitString);
            
            // Afi?are rezultate
            PlotY(freqPanel, PANEL_FREQ_PANEL_FREQ, convertedSpectrum, npoints / 2, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, 1, VAL_BLUE);
            SetCtrlVal(freqPanel, PANEL_FREQ_PANEL_NUMERIC_PEAK_FR, frequencyPeak);
			SetCtrlVal(freqPanel, PANEL_FREQ_PANEL_NUMERIC_PEAK_PW, powerPeak);
            
            // Eliberare memorie
            //free(autoSpectrum);
           // free(convertedSpectrum);
			break;
	}
	return 0;
}



void FiltruOrdinI(double *signal, double *filt, int npoints, double alpha){
	filt[0] = signal[0];
	for(int i = 1; i < npoints ; i ++){
		filt[i] = (1 - alpha) * filt[i - 1] + alpha * signal[i];
	}
}

void FiltruPrinMediere(double *inputData, double *outputData, int dataLength, int windowSize) {
    double *window = (double *)malloc(windowSize * sizeof(double));
    int halfWindow = windowSize / 2;

    for (int i = halfWindow; i < dataLength - halfWindow; i++) {
        for (int j = 0; j < windowSize; j++) {
            window[j] = inputData[i - halfWindow + j];
        }

        
        for (int j = 0; j < windowSize - 1; j++) {
            for (int k = j + 1; k < windowSize; k++) {
                if (window[j] > window[k]) {
                    double temp = window[j];
                    window[j] = window[k];
                    window[k] = temp;
                }
            }
        }

        outputData[i] = window[halfWindow];
    }

    free(window);
}

int CVICALLBACK FilterCommandButton (int panel, int control, int event,
									 void *callbackData, int eventData1, int eventData2)
{
	char window_size[10];
	GetCtrlVal(panel, PANEL_STRING_3, window_size);
    int value = atoi(window_size);
	
	switch (event)
	{
		case EVENT_COMMIT:
			{
				int filterType;
				double alphaValue = ALPHA;
				double *filteredData = (double *)calloc(npoints, sizeof(double));
				double *segmentData = waveData + currentSegmentStart;
				
				GetCtrlVal(panelHandle, PANEL_FILTER_TYPE_RING, &filterType);
				
				
				if (filterType == 2) {  
					GetCtrlVal(panelHandle, PANEL_ALPHA_NUMERIC, &alphaValue);
					 if (alphaValue < 0.00 || alphaValue > 1.00) {
                    		MessagePopup("Eroare", "Valoarea lui Alpha trebuie sa fie între 0 si 1");
                    		free(filteredData); 
                    		return -1;
                	}
					FiltruOrdinI(segmentData, filteredData, sampleRate, alphaValue);
					DeleteGraphPlot(panelHandle, PANEL_GRAPH_2, -1, VAL_IMMEDIATE_DRAW);
					PlotY (panel, PANEL_GRAPH_2, filteredData, sampleRate, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, 1, VAL_RED);
					free(filteredData);
					} 
				else if (filterType == 1) {  
					if (value != 16 && value != 32) {
      					 MessagePopup("Eroare", "Dimensiunea ferestrei trebuie sa fie 16 sau 32.\n");
 	   					 return -1;
    				}	
					FiltruPrinMediere(segmentData, filteredData, sampleRate, value);
					DeleteGraphPlot(panelHandle, PANEL_GRAPH_2, -1, VAL_IMMEDIATE_DRAW);
					PlotY (panel, PANEL_GRAPH_2, filteredData, sampleRate, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, 1, VAL_RED);
					free(filteredData);
				}
			}
			break;
	}
	return 0;
}

int CVICALLBACK OnFilterTypeChanged (int panel, int control, int event,
                                     void *callbackData, int eventData1, int eventData2)
{
    int filterType;
    switch (event)
    {
        case EVENT_COMMIT:
            
            GetCtrlVal(panel, PANEL_FILTER_TYPE_RING, &filterType);

            
            if (filterType == 1) {
                SetCtrlAttribute(panel, PANEL_ALPHA_NUMERIC, ATTR_DIMMED, 0);  
            }
            else {
                SetCtrlAttribute(panel, PANEL_ALPHA_NUMERIC, ATTR_DIMMED, 1);  
            }
            break;
    }
    return 0;
}

void PlotCurrentSegment(int panel) {
    double *segmentData = waveData + currentSegmentStart;
	DeleteGraphPlot(panelHandle, PANEL_MAIN_PANEL_IDC_GRAPH_, -1, VAL_IMMEDIATE_DRAW);
	
	PlotY(panel, PANEL_MAIN_PANEL_IDC_GRAPH_, segmentData, segmentSize, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, 1, VAL_RED);
	double skewness;
	Moment(waveData, 256, 4, &skewness);
	double kurtosis;
	Moment(waveData, 256, 3, &kurtosis);
	double maxVal = 0.0;
	double minVal = 0.0;
	int maxIndex = 0;
	int minIndex = 0;
	double mean = 0.0;
	double dispersia = 0.0;
	MaxMin1D(segmentData,sampleRate,&maxVal,&maxIndex,&minVal,&minIndex);
	Mean(segmentData, sampleRate, &mean);
	Variance(segmentData, sampleRate, &mean, &dispersia);
	        
	int zeroCrossings = 0;
	for (int i = 1; i < sampleRate; i++) {
		if ((segmentData[i-1] < 0 && segmentData[i] > 0) || (segmentData[i-1] > 0 && segmentData[i] < 0)) {
				zeroCrossings++;
				}
			}
	        

			
	SetCtrlVal(panelHandle, PANEL_NUMERIC_3, minVal);
	SetCtrlVal(panelHandle, PANEL_NUMERIC_2, maxVal);
	SetCtrlVal(panelHandle, PANEL_NUMERIC, mean);
	SetCtrlVal(panelHandle, PANEL_NUMERIC_4, dispersia);
	SetCtrlVal(panelHandle, PANEL_NUMERIC_5, zeroCrossings);
	SetCtrlVal(panelHandle, PANEL_NUMERIC_6, skewness);		
	SetCtrlVal(panelHandle,	PANEL_NUMERIC_7, kurtosis);
	SetCtrlVal(panelHandle, PANEL_NUMERIC_8, minIndex);
	SetCtrlVal(panelHandle, PANEL_NUMERIC_9, maxIndex);
}

int CVICALLBACK OnPrev (int panel, int control, int event,
                        void *callbackData, int eventData1, int eventData2) {
	char data_next[50];
	char data_prev[50];
    if (event == EVENT_COMMIT) {
        if (currentSegmentStart >= segmentSize) { 
            currentSegmentStart -= segmentSize;
            PlotCurrentSegment(panel);
			GetCtrlVal(panel, PANEL_STRING_2, data_next);
            int value = atoi(data_next);
            value--;
            sprintf(data_next, "%d", value);
            SetCtrlVal(panel, PANEL_STRING_2, data_next);
			GetCtrlVal(panel, PANEL_STRING, data_prev);
            int value2 = atoi(data_prev);
            value2--;
            sprintf(data_prev, "%d", value2);
            SetCtrlVal(panel, PANEL_STRING, data_prev);
        }
    }
    return 0;
}

int CVICALLBACK OnNext (int panel, int control, int event,
                        void *callbackData, int eventData1, int eventData2) {
	char data_next[50];
	char data_prev[50];
    if (event == EVENT_COMMIT) {
        if (currentSegmentStart + segmentSize < npoints) { 
            currentSegmentStart += segmentSize;
            PlotCurrentSegment(panel);
            GetCtrlVal(panel, PANEL_STRING_2, data_next);
            int value = atoi(data_next);
            value++;
            sprintf(data_next, "%d", value);
            SetCtrlVal(panel, PANEL_STRING_2, data_next);
			GetCtrlVal(panel, PANEL_STRING, data_prev);
            int value2 = atoi(data_prev);
            value2++;
            sprintf(data_prev, "%d", value2);
            SetCtrlVal(panel, PANEL_STRING, data_prev);
        }
    }
    return 0;
}

int CVICALLBACK SaveButton (int panel, int control, int event,
							void *callbackData, int eventData1, int eventData2)
{
	char filePath[260] = {0};
	char filePath2[260] = {0};
	int bitmapID;
	switch (event)
	{
		case EVENT_COMMIT:
			if(FileSelectPopup("", "*.jpg","JPEG Files (*.jpg)","Save As", VAL_SAVE_BUTTON, 0,0,1,0,filePath)){
				GetCtrlDisplayBitmap(panel, PANEL_MAIN_PANEL_IDC_GRAPH_, 1, &bitmapID);
				SaveBitmapToJPEGFile(bitmapID, filePath, JPEG_PROGRESSIVE, 100);
				
			}
			 
			if(FileSelectPopup("", "*.jpg","JPEG Files (*.jpg)","Save filtered As", VAL_SAVE_BUTTON, 0,0,1,0,filePath2)){
				GetCtrlDisplayBitmap(panel, PANEL_GRAPH_2, 1, &bitmapID);
				SaveBitmapToJPEGFile(bitmapID, filePath2, JPEG_PROGRESSIVE, 100);
				
			}
			break;
	}
	return 0;
}

int CVICALLBACK OnEnvelopeButton (int panel, int control, int event,
								  void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			   
	
			double max, min;
			int minI ,maxI;
			double* peakL, *peakA, *peakD, *xMax, *yMax, *xMin, *yMin;
			ssize_t nrMin,nrMax;
			
			MaxMin1D(waveData, sampleRate, &max, &maxI, &min, &minI);
			PeakDetector(waveData, sampleRate, max/20, sampleRate/100*2 , 1, 1, 1, &nrMax, &peakL, &peakA, &peakD);
			xMax=(double*)malloc(sizeof(double) * nrMax);
			yMax=(double*)malloc(sizeof(double) * nrMax);
			for(int i = 0;i < nrMax;i ++)
			{
				xMax[i] = peakL[i];
				yMax[i] = peakA[i];
			}
			
			FreeAnalysisMem(peakL);
			FreeAnalysisMem(peakA);
			FreeAnalysisMem(peakD);
			
			PeakDetector(waveData,sampleRate, min/20, sampleRate/100*2,0, 1, 1, &nrMin, &peakL, &peakA, &peakD);
			xMin=(double*)malloc(sizeof(double) * nrMin);
			yMin=(double*)malloc(sizeof(double) * nrMin);
			for(int i = 0;i < nrMin;i ++)
			{
				xMin[i] = peakL[i];
				yMin[i] = peakA[i];
			}
			FreeAnalysisMem(peakL);
			FreeAnalysisMem(peakA);
			FreeAnalysisMem(peakD);
			DeleteGraphPlot(panel,PANEL_MAIN_PANEL_IDC_GRAPH_,-1,VAL_IMMEDIATE_DRAW);
			PlotXY(panel, PANEL_MAIN_PANEL_IDC_GRAPH_, xMax, yMax, nrMax, VAL_DOUBLE,VAL_DOUBLE, VAL_CONNECTED_POINTS, VAL_SOLID_CIRCLE, VAL_SOLID, 1, VAL_RED);
			PlotXY(panel, PANEL_MAIN_PANEL_IDC_GRAPH_, xMin, yMin, nrMin, VAL_DOUBLE,VAL_DOUBLE, VAL_CONNECTED_POINTS, VAL_SOLID_CIRCLE, VAL_SOLID, 1, VAL_GREEN);

			break;
	}
	return 0;
}
void calculateHistogram(void)
{
	double *segmentData = waveData + currentSegmentStart;
	double maxVal = 0.0;
	double minVal = 0.0;
	int maxIndex = 0;
	int minIndex = 0;
	
	MaxMin1D(segmentData,sampleRate,&maxVal,&maxIndex,&minVal,&minIndex);
		
	ssize_t bins = (ssize_t)sqrt(npoints);
	ssize_t histogram[bins];
	double axis[bins];
	Histogram(waveData, npoints, minVal, maxVal, histogram, axis, bins);
	DeleteGraphPlot(histogramPanel,PANEL_MAIN_PANEL_IDC_GRAPH_,-1,VAL_IMMEDIATE_DRAW);
    PlotXY (histogramPanel, PANEL_MAIN_PANEL_IDC_GRAPH_, axis,  histogram, bins, VAL_DOUBLE, VAL_SSIZE_T, VAL_VERTICAL_BAR, VAL_EMPTY_SQUARE, VAL_SOLID, 1, VAL_RED);
}
int CVICALLBACK ShowHistogramButton (int panel, int control, int event,
                                      void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			calculateHistogram();
			break;
	}
	return 0;
}


int CVICALLBACK OnSwitchPanelCB (int panel, int control, int event,
								 void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			if(panel == panelHandle)
			{
				SetCtrlVal(freqPanel, PANEL_FREQ_IDC_SWITCH_PANEL, 1);
				DisplayPanel(freqPanel);
				HidePanel(panel);
			}
			else
			{
				SetCtrlVal(panelHandle, PANEL_FREQ_IDC_SWITCH_PANEL, 0);
				DisplayPanel(panelHandle);
				HidePanel(panel);
			}
			break;
	}
	return 0;
}

int CVICALLBACK OnFFTSizeChange(int panel, int control, int event, void *callbackData, int eventData1, int eventData2) {
    if (event == EVENT_COMMIT) {
        GetCtrlVal(panel, control, &fftSize);
    }
    return 0;
}



void ApplyWindow(double* signal, int N, int windowType, double* output) {
   /* if (signal == NULL) {
        printf("Eroare: Memorie nealocat? pentru signal.\n");
        return;
    }*/
    //printf("Apelare ApplyWindow cu signal = %p\n", (void*)signal);  // Diagnosticare

    int i;
    double w = 0;

  
    for (i = 0; i < N; i++) {
        if (windowType == 1) {  // Blackman
            w = 0.42 - 0.5 * cos(2 * M_PI * i / (N - 1)) + 0.08 * cos(4 * M_PI * i / (N - 1));
        }
        else if (windowType == 2) {  // FlatTop
            w = 1 - 1.93 * cos(2 * M_PI * i / (N - 1)) + 1.29 * cos(4 * M_PI * i / (N - 1)) - 0.388 * cos(6 * M_PI * i / (N - 1)) + 0.028 * cos(8 * M_PI * i / (N - 1));
        }
        else {
            
            return;
        }
        output[i] = signal[i] * w;  
    }
 
}

int CVICALLBACK WindowCommandButton (int panel, int control, int event,
									 void *callbackData, int eventData1, int eventData2)
{
	
	int N;
	GetCtrlVal(freqPanel, PANEL_FREQ_NUMERIC, &N);
	switch (event)
	{
		case EVENT_COMMIT:
			{
				
				double *filteredData = (double *)calloc(N, sizeof(double));
				GetCtrlVal(freqPanel, PANEL_FREQ_FILTER_TYPE_RING, &windowType);	
				if (windowType == 2) {  
					//GetCtrlVal(panelHandle, PANEL_ALPHA_NUMERIC, &alphaValue);
					 ApplyWindow(convertedSpectrum, N, windowType, filteredData);
					 DeleteGraphPlot(freqPanel, PANEL_FREQ_PANEL_FREQ_MODIFY, -1, VAL_IMMEDIATE_DRAW);
					 PlotY(freqPanel, PANEL_FREQ_PANEL_FREQ_MODIFY, filteredData, N, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, 1, VAL_BLUE);
					} 
				else if (windowType == 1) {  
					ApplyWindow(convertedSpectrum, N, windowType, filteredData);
					DeleteGraphPlot(freqPanel, PANEL_FREQ_PANEL_FREQ_MODIFY, -1, VAL_IMMEDIATE_DRAW);
					PlotY(freqPanel, PANEL_FREQ_PANEL_FREQ_MODIFY, filteredData, N, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, 1, VAL_BLUE);
				}
			}
			break;
	}
	return 0;
}


int ComputePowerSpectrum(double inputArray[], unsigned int numberOfElements, double samplingRate, double convertedSpectrum[], double *frequencyInterval, double *powerPeak, double *freqPeak)
{
	WindowConst windowConstant;
	ScaledWindowEx(inputArray,numberOfElements,RECTANGLE_,0,&windowConstant);
	double* autoSpectrum=malloc(sizeof(double)*(numberOfElements/2));
	AutoPowerSpectrum(inputArray,numberOfElements,1.0/samplingRate,autoSpectrum,frequencyInterval);
	PowerFrequencyEstimate(autoSpectrum,numberOfElements/2,-1,windowConstant,*frequencyInterval,7,freqPeak,powerPeak);
	char l[32]="V";
	SpectrumUnitConversion(autoSpectrum,numberOfElements/2,SPECTRUM_POWER,SCALING_MODE_LINEAR,DISPLAY_UNIT_VRMS,*frequencyInterval,windowConstant,convertedSpectrum,l);
	free(autoSpectrum);
	return 0;
}

void DesignBesselBandpass(double* h, int N, double f1, double f2, double fs) {
    for (int n = 0; n < N; n++) {
        double omega = 2 * M_PI * n / fs;
        if (omega >= 2 * M_PI * f1 / fs && omega <= 2 * M_PI * f2 / fs) {
            h[n] = exp(-omega); 
        } else {
            h[n] = 0.0;
        }
    }
}
void DesignEquirippleBandpass(double* h, int N, double f1, double f2, double fs) {
    double delta = 2.0 / (N - 1);
    for (int n = 0; n < N; n++) {
        double omega = M_PI * n * delta;
        if (omega >= 2 * M_PI * f1 / fs && omega <= 2 * M_PI * f2 / fs) {
            h[n] = 1.0;
        } else {
            h[n] = 0.0;
        }
   }
}


int CVICALLBACK OnTimerTick (int panel, int control, int event,
								 void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_TIMER_TICK:
			DeleteGraphPlot(freqPanel,PANEL_MAIN_PANEL_IDC_GRAPH_,-1,VAL_IMMEDIATE_DRAW);
			PlotY (freqPanel, PANEL_MAIN_PANEL_IDC_GRAPH_, waveData + currentSegmentStart, sampleRate, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, 1, VAL_RED);
			//int DFT = 0;
			GetCtrlVal(freqPanel, PANEL_FREQ_NUMERIC, &DFT);
			ComputePowerSpectrum(waveData + currentSegmentStart + startPoint,DFT,sampleRate,convertedSpectrum,&df,&powerPeak,&frequencyPeak);
			
			SetCtrlVal( freqPanel, PANEL_FREQ_PANEL_NUMERIC_PEAK_FR, frequencyPeak);
			SetCtrlVal( freqPanel, PANEL_FREQ_PANEL_NUMERIC_PEAK_PW, powerPeak);
	
			DeleteGraphPlot (freqPanel, PANEL_FREQ_PANEL_FREQ, -1, VAL_IMMEDIATE_DRAW);
	
    		PlotWaveform(freqPanel, PANEL_FREQ_PANEL_FREQ, convertedSpectrum, DFT/2 ,VAL_DOUBLE, 1.0, 0.0, 0.0, df,
                                    VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID,  VAL_CONNECTED_POINTS, VAL_GREEN);
			double minY,maxY,n;
			GetAxisRange(freqPanel, PANEL_FREQ_MAIN_PANEL_IDC_GRAPH_, (int*)&n, &n, &n,(int*)&n,&minY,&maxY);
			
			/*if(line1!=-1)
			{
				DeleteGraphPlot(freqPanel,PANEL_FREQ_MAIN_PANEL_IDC_GRAPH_,line1,VAL_DELAYED_DRAW);
				DeleteGraphPlot(freqPanel,PANEL_FREQ_MAIN_PANEL_IDC_GRAPH_,line2,VAL_DELAYED_DRAW);
				//DeleteGraphPlot(panel,FREQ_PANEL_GRAPH,fill,VAL_DELAYED_DRAW);
			}*/
			line1=PlotLine(freqPanel, PANEL_FREQ_MAIN_PANEL_IDC_GRAPH_, startPoint, minY, startPoint, maxY, VAL_GREEN);
			line2=PlotLine(freqPanel, PANEL_FREQ_MAIN_PANEL_IDC_GRAPH_, startPoint+DFT, minY, startPoint+DFT, maxY, VAL_GREEN);
			//fill=PlotRectangle(panel, FREQ_PANEL_GRAPH, startPoint, minY, startPoint+DFT, maxY, VAL_WHITE,VAL_TRANSPARENT);
			DeleteGraphPlot (freqPanel, PANEL_FREQ_PANEL_FREQ_MODIFY_2, -1, VAL_IMMEDIATE_DRAW);
			PlotY(freqPanel, PANEL_FREQ_PANEL_FREQ_MODIFY_2, waveData + currentSegmentStart + startPoint, DFT, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_RED);
			double *filteredData = (double *)calloc(DFT, sizeof(double));
			double *h = (double*)malloc(DFT * sizeof(double));
			outputSpectre = (double*)malloc(DFT * sizeof(double));
			GetCtrlVal(freqPanel, PANEL_FREQ_FILTER_TYPE_RING_2, &filterSpectre);
			if (filterSpectre == 1) {  
					//GetCtrlVal(panelHandle, PANEL_ALPHA_NUMERIC, &alphaValue);
					 DesignBesselBandpass(h, DFT, sampleRate / 2, 3 * sampleRate / 4, sampleRate);
					 for (int i = 0; i < DFT; i++) {
					 	outputSpectre[i] = 0.0;
					 	for (int j = 0; j < DFT; j++) {
					        if (i - j >= 0)
           						 outputSpectre[i] += waveData[currentSegmentStart + startPoint + j] * h[i - j];
   							}
					}
			} 
			else if (filterSpectre == 2) {  
					DesignEquirippleBandpass(h, DFT, sampleRate / 4, sampleRate / 2, sampleRate);
					for (int i = 0; i < DFT; i++) {
				    	outputSpectre[i] = 0.0;
				    	for (int j = 0; j < DFT; j++) {
				        	if (i - j >= 0)
				            	outputSpectre[i] += waveData[currentSegmentStart + startPoint + j] * h[i - j];
				   		}
					}
			}
			GetCtrlVal(freqPanel, PANEL_FREQ_FILTER_TYPE_RING, &windowType);	
			if (windowType == 2) {  
					//GetCtrlVal(panelHandle, PANEL_ALPHA_NUMERIC, &alphaValue);
					 ApplyWindow(outputSpectre, DFT, windowType, filteredData);
					 DeleteGraphPlot(freqPanel, PANEL_FREQ_PANEL_FREQ_MODIFY_3, -1, VAL_IMMEDIATE_DRAW);
					 PlotY(freqPanel, PANEL_FREQ_PANEL_FREQ_MODIFY_3, filteredData, DFT, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, 1, VAL_BLUE);
					} 
			else if (windowType == 1) {  
					ApplyWindow(outputSpectre, DFT, windowType, filteredData);
					DeleteGraphPlot(freqPanel, PANEL_FREQ_PANEL_FREQ_MODIFY_3, -1, VAL_IMMEDIATE_DRAW);
					PlotY(freqPanel, PANEL_FREQ_PANEL_FREQ_MODIFY_3, filteredData, DFT, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, 1, VAL_BLUE);
				}
			ComputePowerSpectrum(filteredData,DFT,sampleRate,convertedSpectrum,&df,&powerPeak,&frequencyPeak);
			DeleteGraphPlot (freqPanel, PANEL_FREQ_PANEL_FREQ_MODIFY, -1, VAL_IMMEDIATE_DRAW);
	
    		PlotWaveform(freqPanel, PANEL_FREQ_PANEL_FREQ_MODIFY, convertedSpectrum, DFT/2 ,VAL_DOUBLE, 1.0, 0.0, 0.0, df,
                                    VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID,  VAL_CONNECTED_POINTS, VAL_GREEN);
			
		
			
			
			
			if(startPoint+2*DFT<=sampleRate)
				startPoint+=DFT;
			else if(startPoint>=sampleRate-DFT)
			{
				startPoint=0;
			}
			else
			{
				startPoint=sampleRate-DFT;
			}
			break;
			
	}
	return 0;
}

int CVICALLBACK OnSwitch (int panel, int control, int event,
						  void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			int val;
			GetCtrlVal(freqPanel, control, &val);
			SetCtrlAttribute(freqPanel, PANEL_FREQ_TIMER_2, ATTR_ENABLED, val);
			break;
	}
	return 0;
}

int CVICALLBACK SaveCB (int panel, int control, int event,
						 void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		//int DFT = 0;
		GetCtrlVal(freqPanel, PANEL_FREQ_NUMERIC, &DFT);
		case EVENT_COMMIT:
			char filePath[300]={0};
			GetProjectDir(filePath);
			strcat(filePath,"\\Spectrum\\");
			char fileName[70];
			int part=startPoint/DFT;
			time_t t;
			struct tm* timeInfo;
			char timeStr[50];
			time(&t);
			timeInfo=localtime(&t);
			strftime(timeStr,sizeof(timeStr),"%Y-%m-%d_%H-%M",timeInfo);
			sprintf(fileName,"S%d_Part%d_DFT%d_%s.jpeg",currentSegmentStart + segmentSize,part,DFT,timeStr);
			int bitmapID;
			GetCtrlDisplayBitmap(freqPanel, PANEL_FREQ_MAIN_PANEL_IDC_GRAPH_, 1, &bitmapID);
			strcat(filePath,fileName);
			SaveBitmapToJPEGFile(bitmapID,filePath,JPEG_PROGRESSIVE,100);
			filePath[strlen(filePath)-5]='\0';
			strcat(filePath,"_spectrum_unf.jpeg");
			GetCtrlDisplayBitmap(freqPanel,PANEL_FREQ_PANEL_FREQ,1,&bitmapID);
			SaveBitmapToJPEGFile(bitmapID,filePath,JPEG_PROGRESSIVE,100);
			filePath[strlen(filePath)-strlen("_spectrum_unf.jpeg")]='\0';
			strcat(filePath,"_seg_unf.jpeg");
			GetCtrlDisplayBitmap(freqPanel,PANEL_FREQ_PANEL_FREQ_MODIFY_2,1,&bitmapID);
			SaveBitmapToJPEGFile(bitmapID,filePath,JPEG_PROGRESSIVE,100);
			
			filePath[strlen(filePath)-strlen("_seg_unf.jpeg")]='\0';
			strcat(filePath,"_seg_fil.jpeg");
			GetCtrlDisplayBitmap(freqPanel,PANEL_FREQ_PANEL_FREQ_MODIFY_3,1,&bitmapID);
			SaveBitmapToJPEGFile(bitmapID,filePath,JPEG_PROGRESSIVE,100);
			
			filePath[strlen(filePath)-strlen("_seg_fil.jpeg")]='\0';
			strcat(filePath,"_spectrum_fil.jpeg");
			GetCtrlDisplayBitmap(freqPanel,PANEL_FREQ_PANEL_FREQ_MODIFY,1,&bitmapID);
			SaveBitmapToJPEGFile(bitmapID,filePath,JPEG_PROGRESSIVE,100);
			DiscardBitmap(bitmapID);
			break;
	}
	return 0;
}