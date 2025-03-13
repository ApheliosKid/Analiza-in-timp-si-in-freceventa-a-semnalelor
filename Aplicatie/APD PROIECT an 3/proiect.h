/**************************************************************************/
/* LabWindows/CVI User Interface Resource (UIR) Include File              */
/*                                                                        */
/* WARNING: Do not add to, delete from, or otherwise modify the contents  */
/*          of this include file.                                         */
/**************************************************************************/

#include <userint.h>

#ifdef __cplusplus
    extern "C" {
#endif

     /* Panels and Controls: */

#define  PANEL                            1
#define  PANEL_MAIN_PANEL_IDC_GRAPH_      2       /* control type: graph, callback function: (none) */
#define  PANEL_COMMANDBUTTON              3       /* control type: command, callback function: OnLoadButtonCB */
#define  PANEL_IDC_SWITCH_PANEL           4       /* control type: binary, callback function: OnSwitchPanelCB */
#define  PANEL_GRAPH_2                    5       /* control type: graph, callback function: (none) */
#define  PANEL_COMMANDBUTTON_2            6       /* control type: command, callback function: OnExit */
#define  PANEL_NUMERIC_7                  7       /* control type: numeric, callback function: (none) */
#define  PANEL_NUMERIC_6                  8       /* control type: numeric, callback function: (none) */
#define  PANEL_NUMERIC_9                  9       /* control type: numeric, callback function: (none) */
#define  PANEL_NUMERIC_8                  10      /* control type: numeric, callback function: (none) */
#define  PANEL_NUMERIC_5                  11      /* control type: numeric, callback function: (none) */
#define  PANEL_NUMERIC_4                  12      /* control type: numeric, callback function: (none) */
#define  PANEL_NUMERIC_3                  13      /* control type: numeric, callback function: (none) */
#define  PANEL_NUMERIC_2                  14      /* control type: numeric, callback function: (none) */
#define  PANEL_NUMERIC                    15      /* control type: numeric, callback function: (none) */
#define  PANEL_COMMANDBUTTON_4            16      /* control type: command, callback function: OnNext */
#define  PANEL_COMMANDBUTTON_3            17      /* control type: command, callback function: OnPrev */
#define  PANEL_STRING_2                   18      /* control type: string, callback function: (none) */
#define  PANEL_STRING                     19      /* control type: string, callback function: (none) */
#define  PANEL_TEXTMSG                    20      /* control type: textMsg, callback function: (none) */
#define  PANEL_TEXTMSG_2                  21      /* control type: textMsg, callback function: (none) */
#define  PANEL_STRING_3                   22      /* control type: string, callback function: (none) */
#define  PANEL_COMMANDBUTTON_7            23      /* control type: command, callback function: OnEnvelopeButton */
#define  PANEL_COMMANDBUTTON_8            24      /* control type: command, callback function: ShowHistogramButton */
#define  PANEL_COMMANDBUTTON_5            25      /* control type: command, callback function: FilterCommandButton */
#define  PANEL_FILTER_TYPE_RING           26      /* control type: ring, callback function: OnFilterTypeChanged */
#define  PANEL_ALPHA_NUMERIC              27      /* control type: numeric, callback function: (none) */
#define  PANEL_COMMANDBUTTON_6            28      /* control type: command, callback function: SaveButton */

#define  PANEL_FREQ                       2
#define  PANEL_FREQ_MAIN_PANEL_IDC_GRAPH_ 2       /* control type: graph, callback function: (none) */
#define  PANEL_FREQ_PANEL_FREQ_MODIFY_3   3       /* control type: graph, callback function: (none) */
#define  PANEL_FREQ_PANEL_FREQ_MODIFY_2   4       /* control type: graph, callback function: (none) */
#define  PANEL_FREQ_PANEL_FREQ_MODIFY     5       /* control type: graph, callback function: (none) */
#define  PANEL_FREQ_COMMANDBUTTON_2       6       /* control type: command, callback function: OnExit */
#define  PANEL_FREQ_PANEL_FREQ            7       /* control type: graph, callback function: (none) */
#define  PANEL_FREQ_NUMERIC               8       /* control type: numeric, callback function: OnFFTSizeChange */
#define  PANEL_FREQ_IDC_SWITCH_PANEL      9       /* control type: binary, callback function: OnSwitchPanelCB */
#define  PANEL_FREQ_TIMER_2               10      /* control type: timer, callback function: OnTimerTick */
#define  PANEL_FREQ_PANEL_NUMERIC_PEAK_PW 11      /* control type: numeric, callback function: (none) */
#define  PANEL_FREQ_PANEL_NUMERIC_PEAK_FR 12      /* control type: numeric, callback function: (none) */
#define  PANEL_FREQ_FILTER_TYPE_RING_2    13      /* control type: ring, callback function: OnFilterTypeChanged */
#define  PANEL_FREQ_FILTER_TYPE_RING      14      /* control type: ring, callback function: (none) */
#define  PANEL_FREQ_COMMANDBUTTON_4       15      /* control type: command, callback function: OnNext */
#define  PANEL_FREQ_COMMANDBUTTON_5       16      /* control type: command, callback function: SaveCB */
#define  PANEL_FREQ_COMMANDBUTTON_3       17      /* control type: command, callback function: OnPrev */
#define  PANEL_FREQ_STRING_2              18      /* control type: string, callback function: (none) */
#define  PANEL_FREQ_STRING                19      /* control type: string, callback function: (none) */
#define  PANEL_FREQ_TEXTMSG               20      /* control type: textMsg, callback function: (none) */
#define  PANEL_FREQ_BINARYSWITCH          21      /* control type: binary, callback function: OnSwitch */


     /* Control Arrays: */

#define  CTRLARRAY                        1

     /* Menu Bars, Menus, and Menu Items: */

          /* (no menu bars in the resource file) */


     /* Callback Prototypes: */

int  CVICALLBACK FilterCommandButton(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK OnEnvelopeButton(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK OnExit(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK OnFFTSizeChange(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK OnFilterTypeChanged(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK OnLoadButtonCB(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK OnNext(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK OnPrev(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK OnSwitch(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK OnSwitchPanelCB(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK OnTimerTick(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK SaveButton(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK SaveCB(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK ShowHistogramButton(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);


#ifdef __cplusplus
    }
#endif