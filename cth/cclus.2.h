/*
   Defs for cth clustering program.

   Modification History
   14-04-02  dale   Create

*/

#ifndef CTH_CLUSTER_H
#define CTH_CLUSTER_H

#define _FILE_OFFSET_BITS 64

#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <sstream>
#include <iomanip>
#include <map>
#include <stack>
#include <cfloat>
#include <cmath>
#include <stdlib.h>
#include <getopt.h>
#include <errno.h>
#include <ctype.h>
#include <libgen.h>

using namespace std;

extern bool Debug;

const int E_MARKER=98;
const int I_MARKER=97;
const int MIN_CELL=100;
const int MAX_CELL=999;

/*
   adt file: each tick is 0.5 ms  0.0005 sec  500 us  5,000 ticks/sec
   bdt file: each tick is 0.5 ms  0.0005 sec  500 us  5,000 ticks/sec
   edt file: each tick is 0.1 ms  0.0001 sec  100 us 10,000 ticks/sec
*/
const double ADT_TICK=0.5;  // ms
const double BDT_TICK=0.5;
const double EDT_TICK=0.1;


const double ADT_SECS=1000/ADT_TICK;
const double BDT_SECS=1000/BDT_TICK;
const double EDT_SECS=1000/EDT_TICK;


enum F_TYPE { EDT, BDT, ADT };

/* The input file list has several fields taken from a spreadsheet.
   Fields are tab-delimited, fields can be missing, but tabs need to be there.
   Format:
   path    fname    boundary code    len of ctl period (secs)   #cells    notes
*/

class FileFacts
{
   public:
      FileFacts():pmark(-1),ctl_period_len(-1),num_cells(-1) {};

      double TimeVal(long time)  // given time in file units, return milliseconds
      {
         switch (ftype)
         {
            case BDT:
               return time * BDT_TICK; break;
            case EDT:
               return time * EDT_TICK; break;
            default:
               return time * ADT_TICK; break;
         }
      }

      double TickVal()  // in seconds
      {
         switch (ftype)
         {
            case BDT:
               return BDT_TICK; break;
            case EDT:
               return EDT_TICK; break;
            default:
               return ADT_TICK; break;
         }
      }
      void print() {
                     cout << "Path:            " << this->path << endl;
                     cout << "Fname:           " << this->fname << endl;
                     cout << "Period mark:     " << this->pmark_str << endl;
                     cout << "Period mark val: " << this->pmark << endl;
                     cout << "File Type:       " << this->ftype << endl;
                     cout << "Ctl Period Len:  " << this->ctl_period_len << endl;
                     cout << "Notes:           " << this->notes << endl;
                     cout << endl;
                   };

   public:
      string   path;
      string   fname;
      string   cmplt_name;
      string   pmark_str;        // start of control or stim period
      int      pmark;
      int      ctl_period_len;
      int      num_cells;
      string   notes;
      F_TYPE   ftype;
};

typedef vector <FileFacts> VFiles;
typedef vector <FileFacts>::iterator VFilesIter;


/* One text line of data, hopefully  "[chan][time]"
*/
class DataLine
{
  public:
     DataLine(string& line)
     {
        chan_str = line.substr(0,5);
        chan = atol(chan_str.c_str());  // no error checking
        time_str = line.substr(5,10);
        time = atol(time_str.c_str());
     }
     long GetTime()  { return time; }
     long GetChan()  { return chan; }
     bool IsEPhase() { return chan == E_MARKER; }
     bool IsIPhase() { return chan == I_MARKER; }
     bool IsNeuron() { return chan >= MIN_CELL && chan <= MAX_CELL; }
     bool IsMarker(long marker) { return chan == marker; }

   private:
      DataLine() {}
      string time_str;
      string chan_str;
      long time;
      long chan;
};


/* We don't know how many of these we will need, it depends on the number of
   neuron channels in a file.  So we'll add a NeuroBin to the
   map container as we find new neurons in a file stream.
*/
class NeuroBins;
typedef map <long, NeuroBins>  FBins;
typedef FBins::iterator FBinsIter;

/* One map entry per file, the key is the path/filename
   Holds a bunch of NeuroBins for each file.
*/
typedef map <string,FBins> AllBins;
typedef AllBins::iterator AllBinsIter;

/* A single bin.  We accumulate time and spike counts,
   Algorithm is:   (n0/t + n1/t + n2/t . . . ) / num of cycles
                  That is, mean of the spikes/second
*/
class OneBin
{
   public:
      OneBin():time(0.0),spikes(0),rate(0.0), mean_rate(0.0), total_spikes(0) {}
      OneBin(double t, double sec): time(t),spikes(0),rate(0),
                                    mean_rate(0), convert(sec),total_spikes(0) {}
      virtual ~OneBin() {}

      void SetBinTime(double t, double cvt) { time = t; convert = cvt;}
      void IncSpike() { ++spikes; ++total_spikes;
                       if (Debug) cout << "got " << spikes << " spikes" << endl;
                      }
      unsigned long GetSpikes() {return spikes;}
      void ResetSpike() { spikes = 0; }
      double UpdateRate() { rate += spikes/((time*convert)/1000);  // in seconds
                             if (Debug) cout << endl << "Spikes: " << spikes << "  time " \
                             << time << " convert " << convert << " Rate: " << rate; return rate;
                          }
      void CalcMeanRate(long num_cycles) {mean_rate=rate/num_cycles;}
      double GetRate() { return rate;}
      double GetMeanRate() { return mean_rate;}
      unsigned long GetTotalSpikes() { return total_spikes; }

   private:
      double time;          // current bin width in file time units, varies
      unsigned long spikes;
      double rate;         // in spikes/second
      double mean_rate;
      double convert;      // secs/tick, e.g. edt = 0.1
      unsigned long total_spikes;

};
typedef vector <OneBin> BinArray;
typedef vector <OneBin>::iterator BinArrayIter;



/* Store data in the current phase for saving in the next
   phase's bins here.
*/
typedef stack <DataLine> DStack;

/* This is a collection of bins.  There will be one per spike channel in the files.
   The default is 10 e and 10 i bins, but can be any number (and can even
   be different.)
   The union of the e and i bins creates an n-dimensional point we save to a file
*/
class NeuroBins
{
   public:
      NeuroBins(int c_num ,size_t e_num, size_t i_num, double e_t, double i_t, double s_per_t)
      {
         chan = c_num;
         num_e_bins = e_num;
         num_i_bins = i_num;
         e_bins.assign(e_num,OneBin(e_t,s_per_t));
         i_bins.assign(i_num,OneBin(i_t,s_per_t));
      }
      virtual ~NeuroBins () {}

      static void ResetNumCycles() { num_cycles = 0; }
      static void ResetPeriod() { period_start = period_end = 0; }
      static void PeriodStart(unsigned long p_s) { period_start = p_s; if (Debug) cout << "PS: " << p_s << endl; }
      static void PeriodEnd(unsigned long p_e) { period_end = p_e; if (Debug) cout << "PE: " << p_e << endl; }
      static void IncNumCycles() { ++num_cycles; }
      static long GetNumCycles() { return num_cycles; }
      static void SetSecsTick(double s_p_t) { secs_per_tick = s_p_t; if (Debug) cout << "SPT: " << s_p_t << endl; }
      static double GetSecsTick() {return secs_per_tick;}
 
      static void SetUpIntervals(long e_s, long i_s, long n_e)
                          {
                            e_start = e_s;
                            i_start = i_s;
                            next_e  = n_e;
                            e_time = (double) (i_s - e_s) / num_e_bins;
                            i_time = (double) (n_e - i_s) / num_i_bins;
                          }
      static void SetNumEBins(int e_b)  {num_e_bins = e_b;}
      static void SetNumIBins(int i_b)  {num_i_bins = i_b;}
      static double GetETime()  {return e_time;}
      static double GetITime()  {return i_time;}


      size_t GetNumEBins() {return e_bins.size();}
      size_t GetNumIBins() {return i_bins.size();}
      void UpdateE(DataLine&);
      void UpdateI(DataLine&);
      void UpdateBinRates();
      void UpdateBinTimes();
      void MeanRate();
      int  GetChan() { return chan; }

      static void PrintPeriodTime() {cout << "Period Time: " << (period_end - period_start) * secs_per_tick/1000 << " secs" << endl; }
      void PrintEBins() ;
      void PrintIBins();

      double GetBinPeak();
      double GetBinArea();
      int GetTotalSpikes();

      void WriteOctaveHeader(ofstream &, string, int, int);
      void WriteOEBins(ofstream &, double);
      void WriteOIBins(ofstream &, double);
      void WriteCEBins(ofstream &, double);
      void WriteCIBins(ofstream &, double);
      void WriteGCHeader(ofstream &);
      void WriteGCEBins(ofstream &, double, int);
      void WriteGCIBins(ofstream &, double);


          // manage stacked events for next phase
      bool HaveE() { return !saved_e.empty(); }
      bool HaveI() { return !saved_i.empty(); }
      void SaveForE(DataLine d) { saved_e.push(d); }
      void SaveForI(DataLine d) { saved_i.push(d); }
      DataLine& GetE() { return  saved_e.top();}
      DataLine& GetI() { return  saved_i.top();}
      void RemoveE() { saved_e.pop(); }
      void RemoveI() { saved_i.pop(); }

   private:
      NeuroBins(){}

      BinArray e_bins;  // bin data managed by these
      BinArray i_bins;
      DStack saved_e;   // boundary events that belong in
      DStack saved_i;   // next phase goes into these
      int chan;

        // shared by all instances
      static long e_start;      // time of current e phase in file units
      static long i_start;      // time of current i phase
      static long next_e;       // time of next e phase
      static long period_start; // time of next e phase
      static long period_end;   // time of next e phase
      static long num_cycles;   // total number of cycles we accumulated counts for
      static unsigned int num_e_bins;
      static unsigned int num_i_bins;
      static double e_time;      // width of e bin in file time units
      static double i_time;
      static double secs_per_tick;  // file secs/tick (.1, .5)

};


#endif
