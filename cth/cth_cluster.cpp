/* Program to scan .edt, .bdt, and .adt files, accumulate firing rates for 
   control and stim periods on neuron channels into a variable number of bins,
   and generate a file with the bins representing N-dimensional points for
   subsequent processing by other software.

   Modification History
   09-May-17  dale Add in cough ctl interval/stim #cycles experiment support.
   20-Apr-17  dale Add in filtering experiments by vagotomized flag.
   05-Aug-16  dale Anniversary!  Adding in real 2nd period handling.  Markers
                   plus interval values for ctl and stim portions.  
                   Now a cell can generate 2 CTHs.
   25-Aug-15  dale Finally bite the bullet and create a holder for per-period data.
   14-Apr-15  dale Add in writing info processing by R programs.
   03-Nov-14  dale Interface to .xls files and read in stereotaxic coordinates.
   30-sep-14  dale Option to split periods into separate files.
   04-aug-14  dale Start adding in new distance calculation and creating
                   a distance matrix for octave programs.
   31-jul-14  dale I's now before E's, caculate stddev and other numbers
                   to create a distance matrix using a non-eucldean
                   distance function.
   06-may-14  dale Just use first (ctl) period, not all periods.
   02-apr-14  dale Create
*/

#include "cth_cluster.h"

// Globals
string InName;
string OutName;
int    NumEBins=50;
int    NumIBins=50;
bool   Debug = false;
bool   Oct = true;      // make a file Octave likes
bool   Csv = false;     // make a generic .cvs file 
bool   Gcsv = false;    // make a .csv file ggobi likes
NORM   Norm = MEAN;     // type of normalization 
bool   AllP  = false;   // include all periods in output
bool   DoCurves = false; // use R program to curve-fit cths
                        // (or not, this can take a looong time)
bool   WrCurves = false; // write outfile for cth curve fitting R prog
                         // but do not run the program.
bool   RdCurves = false; // read from previous run, don't redo curves
string RdFile;           // over-ride autofile name and use this one instead
                         // (If you create a curve file with all possible CTHs in it
                         // it can be used for subsets and experimenting with different
                         // bin sizes and other params without having to generate.) 
                         // a new curve file.) 
int vago_state = 0;
bool have_swall1=false;
bool have_lareflex=false;

unsigned long RejectThresh; // neurons with fewer than this many are not saved
unsigned long NumSparse = 0;  // This many of them

// This is data global to all periods in all files
PhaseHist I_Hist[MAX_PERIODS];
PhaseHist E_Hist[MAX_PERIODS];
array <double,MAX_PERIODS> maxCTH;
bool FilePass1=true;


// Class globals shared by OneBin and NeuroBins objects.
long NeuroBins::i_start;               // varies per cycle
long NeuroBins::e_start;
long NeuroBins::next_i; 
long NeuroBins::period_start;
long NeuroBins::period_end;
unsigned long NeuroBins::stim_end;

double NeuroBins::secs_per_tick;      // varies per file

unsigned int  NeuroBins::num_e_bins;
unsigned int  NeuroBins::num_i_bins;

double NeuroBins::e_time;      // width of e bin in file time units
double NeuroBins::i_time;      // width of i bin
long   NeuroBins::ep_time;     // width of e phase in file ticks
long   NeuroBins::ip_time;     // width of i bin

static void Usage(char *name)
{
   cout << endl << "Program to:" << endl;
   cout << "(1) Read .edt, .bdt, and .adt files."  << endl;
   cout << "(2) Accumulate firing rates of neuron channels in control and " << endl;
   cout << "     optionally other periods into a variable number of bins." << endl;
   cout << "(3) Generate text file(s) with CTHs and other information." << endl;
   cout << "The file(s) created by this program are inputs to cth cluster analysis software." << endl;
   cout << endl << "Usage: " << name << " -f spreadsheet.xls [-o basename] [-b number_of_bins | -e number of E bins -i number of I bins] [-np | -nu | -nr] [-c] [-g] [-pa] [-tc | -tw | -tr] [-h] -v[0 | 1 | 2]" << endl << endl;
   cout << "-f spreadsheet.xls  Spreadsheet file containing experiment info." << endl;
   cout << "                    This is a required argument." << endl;
   cout << "-o basename  Use basename for all input and output files." << endl;
   cout << "             The final output will be basename.cth." << endl;
   cout << "-c Create a .csv file.  Will create file: basename.csv " << endl;
   cout << "-g Create a .csv file that ggobi likes.  Will create file: basename_g.csv" << endl;
   cout << "-np Peak normalization. Default is mean (area) normalization." << endl;
   cout << "-nu Unit normalization." << endl; 
   cout << "-nn No normalization, output is raw spikes/sec values." << endl; 
   cout << "-pa Include all periods in files and save to one file." << endl << "    Default is to save Control and optionally Control/Stim periods in separate files." << endl << "    Note that the brainstem program will not be able to load .csv files" << endl << "    created from these files. The clustering program can load these files." << endl;
   cout << "-tc Do complete histogram curve fitting related processing using cth_curve.R." << endl;
   cout << "    By default, this is not done and the output file will contain zero values" << endl;
   cout << "    for the curve plotting vectors if no other -t option is used." << endl;
   cout << "    It can take a VERY long time to do this on the typical workstation." << endl;
   cout << "    (If you really want to do this, run this on a computational server and" << endl;
   cout << "    copy the output file to someplace visible to your workstation.)" << endl;
   cout << "-tw Create an output file that cth_curve.R can read, but do not run cth_curve.R" << endl;
   cout << "    The file will be basename_to_r.txt." << endl; 
   cout << "-tr Read the file created previously by cth_curve.R and add it to the output." << endl;
   cout << "    By default, the file will be basename_from_r.txt (see next option.)" << endl;
   cout << "-t filename   Do not automatically create a curve fitting filename from the -o option." << endl;
   cout << "              Read filename instead.  This lets you share a curve file.  Implies -tr." << endl;
   cout << "-vn, n=0 Include only non-vagotomized experiments. This is the default." << endl;
   cout << "     n=1 Include only vagotomimized experiments." << endl;
   cout << "     n=2 Include all experiments." << endl;

   cout << "-h This help." << endl;
   cout << "Default output file name: cth_[# of bins]_type.cth" << endl;
   cout << "The _type will be one of ctl, vco2, cco2, tb_cgh, lar_cgh, swall1, or allp" << endl;
   cout << "Default E bins:   50" << endl;
   cout << "Default I bins:   50" << endl;
}

static int ParseArgs(int argc, char *argv[])
{
   int cmd;
   const char* optstr = "t:n:p:v:";
   const struct option opts[] = { 
                                   {"f", required_argument, nullptr, 'f'},
                                   {"o", required_argument, nullptr, 'o'},
                                   {"b", required_argument, nullptr, 'b'},
                                   {"e", required_argument, nullptr, 'e'},
                                   {"i", required_argument, nullptr, 'i'},
                                   {"h", no_argument, nullptr, 'h'},
                                   {"d", no_argument, nullptr, 'd'},
                                   {"c", no_argument, nullptr, 'c'},
                                   {"g", no_argument, nullptr, 'g'},
                                   {"k", no_argument, nullptr, 'k'},
                                   {"r", no_argument, nullptr, 'r'},
                                   { 0,0,0,0} };

   while ((cmd = getopt_long_only(argc, argv, optstr, opts, nullptr )) != -1)
   {
      switch (cmd)
      {
         case 'f':
               InName=optarg;
               if (InName.size() == 0)
               {
                  cout << "Input file is missing, aborting. . .\n";
                  exit(1);
               }
               break;

         case 'o':
               OutName=optarg;
               if (OutName.size() == 0)
               {
                  cout << "Output file is missing, aborting. . .\n";
                  exit(1);
               }
               break;

         case 'b':
               errno = 0;
               NumEBins=strtol(optarg,nullptr,10);
               if (errno != 0)
               {
                  cout << "Number of bins is missing, aborting. . .\n";
                  exit(1);
               }
               if ((NumEBins % 2) != 0)
               {
                  cout << "Number of bins must be even, aborting. . .\n";
                  exit(1);
               }
               NumEBins /= 2;  // split in half
               NumIBins = NumEBins;
               break;

         case 'e':
               errno = 0;
               NumEBins=strtol(optarg,nullptr,10);
               if (errno != 0)
               {
                  cout << "Number of E bins is missing, aborting. . .\n";
                  exit(1);
               }
               if ((NumEBins % 2) != 0)
               {
                  cout << "Number of e bins must be even, aborting. . .\n";
                  exit(1);
               }
               break;

         case 'i':
               errno = 0;
               NumIBins=strtol(optarg,nullptr,10);
               if (errno != 0)
               {
                  cout << "Number of I bins is missing, aborting. . .\n";
                  exit(1);
               }
               if ((NumIBins % 2) != 0)
               {
                  cout << "Number of i bins must be even, aborting. . .\n";
                  exit(1);
               }
               break;

         case 'h':
            Usage(argv[0]);
            exit(1);
            break;

         case 'd':
            cout << "Debug turned on" << endl;
            Debug = true;
            break;

         case 'c':
            Csv = true;
            break;

         case 'g':
            Gcsv = true;
            break;

         case 'n':
            if (!strcmp(optarg,"p"))
               Norm = PEAK;
            else if (!strcmp(optarg,"u"))
               Norm = UNIT;
            else if (!strcmp(optarg,"r"))
               Norm = NONE;
            else
               cout << "Unrecognized -n option: " << optarg << " (ignored)" << endl;
            break;

         case 'v':
            if (!strcmp(optarg,"0"))
               vago_state = 0;
            else if (!strcmp(optarg,"1"))
               vago_state = 1;
            else if (!strcmp(optarg,"2"))
               vago_state = 2;
            else
               cout << "Unrecognized -v option: " << optarg << " (ignored)" << endl;
            break;

         case 'p':
            if (!strcmp(optarg,"a"))
               AllP = true;
            else
               cout << "Unrecognized -p option: " << optarg << " (ignored)" << endl;
            break;
         case 't':
            if (!strcmp(optarg,"c"))
               DoCurves=true;
            else if (!strcmp(optarg,"r"))
            {
               DoCurves=false;
               RdCurves=true;
            }
            else if (!strcmp(optarg,"w"))
            {
               DoCurves=false;
               WrCurves=true;
            }
            else  // assume if -t something that something is a file name
            {
               RdFile = optarg;
               DoCurves=false;
               RdCurves=true;
            }
            break;

         case '?':
         default:
           Usage(argv[0]);
           break;
      }
   }

     // use defaults for missing options
   if (InName.size() == 0)
   {
      cout << "An input .xls file that contains a list of input files and other " << endl;
      cout << "information is required." << endl << "Aborting. . ." << endl;
      Usage(argv[0]);
      exit(1);
   }

   return true;
}

void Tokenize(const string& str, vector<string>& tokens, const string& delimiters = " ")
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}



// Update some running tallies at the end of a cycle
void OneBin::UpdateBinStats(PeriodInfo& p_data) 
{
   double nrate = spikes/((time*convert)/1000);  // in seconds
   rate += nrate;
   if (Debug) cout << "  " << "Spikes: " << spikes << " total spikes " << total_spikes << "  time " << time << " Rate: " << rate << " nrate: " << nrate  << " spikes/sec" << endl; 

   // running totals for SD
   double delta = nrate-mean;
   mean += delta/p_data.num_cycles;
   mean_sqr += delta*(nrate-mean);
}


// Finished, caculate some interesting stats.
void OneBin::CalcBinStats(PeriodInfo& p_data) 
{
   mean_rate=rate/p_data.num_cycles; 
   sd = sqrt(mean_sqr/p_data.num_cycles); 
   std_err = sd/sqrt(p_data.num_cycles);
}

/* Have a new neuron event for the E bins, 
   update running sums.
*/
void NeuroBins::UpdateE(DataLine& d_line)
{
   unsigned int idx = (d_line.GetTime() - e_start) / e_time;

   if (idx < 0)
      cout << "E:  We are skipping cycles, need to deal with this case." << endl;

      // If a neuron event happens at the same time as the next I marker and
      // occurs in the file before the marker, put it on the I phase's stack
      // for inclusion in the next I phase.
   if (idx >= num_e_bins)
   {
      SaveForI(d_line);
      if (Debug) cout << "time " << d_line.GetTime() << " idx " << idx << " E SAVED FOR NEXT I" << endl;
   }
   else
   {
      e_bins[idx].IncSpike();
      if (Debug) cout << "update e bin # " << idx << endl;
   }
}

/* Have a new neuron event for the I bins, 
   update running sums.
*/
void NeuroBins::UpdateI(DataLine& d_line)
{
   unsigned int idx = (d_line.GetTime() - i_start) / i_time;
      // If a neuron event happens at the same time as the next E marker and
      // occurs in the file before the marker, put it on the E phase's stack
      // for inclusion in the next E phase in the next cycle.
   if (idx < 0)
      cout << "I:  We are skipping cycles, need to deal with this case." << endl;
   if (idx >= num_i_bins)
   {
      SaveForE(d_line);
      if (Debug) cout << "time " << d_line.GetTime() << " idx " << idx << " I SAVED FOR NEXT E" << endl;
   }
   else
   {
      i_bins[idx].IncSpike();
      if (Debug) cout << "update i bin # " << idx << endl;
   }
}

/* Have a new neuron event for the swallow1 bins, 
   update running sums.
*/
void NeuroBins::UpdateS(DataLine& d_line)
{
   unsigned int idx = (d_line.GetTime() - i_start) / i_time;

     // if we hit the end of a swallow, toss event, these are not
     // contiguous periods of time, there is no carry over
   if (idx < 0)
   {
      cout << "I:  We are skipping cycles, need to deal with this case." << endl;
      return; // odd
   }
   if (idx < num_i_bins)
      i_bins[idx].IncSpike();
   else if (idx < num_i_bins+num_e_bins)
      e_bins[idx-num_e_bins].IncSpike();

}

// Update global phase duration histogram counters for I and E for current period/file
void NeuroBins::UpdatePHist(int period)
{
   PhaseHist::iterator ph;

   ph = I_Hist[period].find(ip_time);
   if (ph == I_Hist[period].end())
      ph = I_Hist[period].insert(make_pair(ip_time,1)).first;
   else
      ph->second++;

   ph = E_Hist[period].find(ep_time);
   if (ph == E_Hist[period].end())
      ph = E_Hist[period].insert(make_pair(ep_time,1)).first;
   else
      ph->second++;
}

void PrintIHist()
{
   for (int  per=CTL_PERIOD; per < MAX_PERIODS; ++per)
   {
      cout << "I hist " << per << " len: " << I_Hist[per].size() << endl;
      for (PhaseHistIter it = I_Hist[per].begin(); it != I_Hist[per].end(); ++it)
         cout << (*it).first << "  " << (*it).second << endl;
   }
}

void PrintEHist()
{
   for (int per=CTL_PERIOD; per < MAX_PERIODS; ++per)
   {
      cout << "E hist " << per << " len: " << E_Hist[per].size() << endl;
      for (PhaseHistIter it = E_Hist[per].begin(); it != E_Hist[per].end(); ++it)
         cout << (*it).first << "  " << (*it).second << endl;
   }
}


// Scan one tick histograms and calculate info for second pass for a period
void PeriodCalcs(PeriodInfo& p_data)
{
   size_t count;
   size_t maxcnt;
   int max_at = 1;

   if (I_Hist[p_data.period].size() == 0)
      cout << " OOPS! Where'd my data go? " << endl;

   if (I_Hist[p_data.period].size() == 1) // for lareflex, all intervals are same
   {                                      // breaks logic below. only one bin in hist
                                          // use it
      PhaseHistIter it = I_Hist[p_data.period].begin();
      int ticks = it->first;  // interval size in ticks
      p_data.i_short = ticks+1;
      p_data.i_max = ticks+1;
      it = E_Hist[p_data.period].begin();
      ticks = it->first;
      p_data.e_short = ticks+1;
      p_data.e_max = ticks+1;
      I_Hist[p_data.period].clear();
      E_Hist[p_data.period].clear();
      return;

   }
   count = I_Hist[p_data.period][1];
   maxcnt = count;

   for (size_t lo = 2; 2 * lo < I_Hist[p_data.period].size(); ++lo) 
   {
      count -= I_Hist[p_data.period][lo - 1];
      if (Debug) cout << "lo: " << lo << " count- " << count << "  maxcnt " << maxcnt << endl;
      for (size_t i = 2 * (lo - 1); i < 2 * lo; ++i) {
      count += I_Hist[p_data.period][i];
      if (Debug) cout << "  i: " << i << "  count+ " << count << "  maxcnt " << maxcnt << endl;
    }
    if (count > maxcnt) {
      maxcnt = count;
      max_at = lo;
      if (Debug) cout << "  i max_at " << max_at << endl;
    }
   }
   p_data.i_short = max_at;
   p_data.i_max = max_at * 2 - 1;

   if (Debug) cout << "I one tick : " << p_data.i_short << " to " << p_data.i_max << endl;

   count = E_Hist[p_data.period][1];
   maxcnt = count;
   max_at = 1;

   for (size_t lo = 2; 2 * lo < E_Hist[p_data.period].size(); ++lo) 
   {
      count -= E_Hist[p_data.period][lo - 1];
      if (Debug) cout << "lo: " << lo << " count- " << count << "  maxcnt " << maxcnt << endl;
      for (size_t i = 2 * (lo - 1); i < 2 * lo; ++i) {
      count += E_Hist[p_data.period][i];
      if (Debug) cout << "  e: " << i << "  count+ " << count << "  maxcnt " << maxcnt << endl;
    }
    if (count > maxcnt) {
      maxcnt = count;
      max_at = lo;
      if (Debug) cout << "  e max_at " << max_at << endl;
    }
  }

   if (Debug) cout << "E one tick : " << max_at << " to " << max_at*2-1 << endl;

   p_data.e_short = max_at;
   p_data.e_max = max_at * 2 - 1;

    // We are done with these, set up for next period/file
   I_Hist[p_data.period].clear();
   E_Hist[p_data.period].clear();
}


// index for which one tick bin (if any) the current event belongs in
int
get_bin (int j, int n_p, int n_b)
{
  double i = ((2. * j + 1) * n_b / n_p - 1) / 2;
  int bin_tick = floor ((0.5 + i) * n_p / n_b);
  if (bin_tick == j) return i;
  bin_tick = floor ((0.5 + i + 1) * n_p / n_b);
  if (bin_tick == j) return i + 1;
  return -1;
}

// If current cycle has an i or e phase that is too short or too long
// don't accumulate spikes that may be in it.
bool NeuroBins::UseCycle(PeriodInfo& p_data)
{
   bool res = true;

   if (ip_time < p_data.i_short)
   {
      ++p_data.i_min_tossed;
      res = false;
   }
   else if (ip_time > p_data.i_max)
   {
      ++p_data.i_max_tossed;
      res = false;
   }

   if (ep_time < p_data.e_short)
   {
      ++p_data.e_min_tossed;
      res = false;
   }
   else if (ep_time > p_data.e_max)
   {
      ++p_data.e_max_tossed;
      res = false;
   }
   if (res)
      ++p_data.total_accept;
   else
   {
      ++p_data.total_reject;
      if (Debug) cout << "Reject one tick for " << p_data.period << "  total:  " << p_data.total_reject << endl;
   }

   return res;
}

// I one tick bins
void NeuroBins::UpdateI_OT(DataLine& d_line)
{
   int t = d_line.GetTime();
   int j = t - i_start;
   int n_p = ip_time;
   int n_b = i_onetick.size();
   int idx;

   if (t == e_start)  // this belongs in next E phase
   {
      if (Debug) cout << "Save i_ot time " << t << " for next e_ot" << endl;
      SaveForE_OT(d_line);
   }
   else
   {
     idx = get_bin(j,n_p,n_b);
     if (idx >= 0)
     {
        ++i_onetick[idx];
        if (Debug) cout << "update i_ot time: " << t << " bin # " << idx << " so far: " << i_onetick[idx] << endl;
     }
     else
     {
       if (Debug) cout << "rejected i_ot time: " << t << " idx: " << idx << endl;
     }
   }
}

void NeuroBins::UpdateE_OT(DataLine& d_line)
{
   int t = d_line.GetTime();
   int j = t - e_start;
   int n_p = ep_time;
   int n_b = e_onetick.size();
   int idx;

   if (t == next_i)  // this belongs in next I phase
   {
      if (Debug) cout << "Save e_ot time " << t << " for next i_ot" << endl;
      SaveForI_OT(d_line);
   }
   else
   {
      idx = get_bin(j,n_p,n_b);
      if (idx >= 0)
      {
         ++e_onetick[idx];
        if (Debug) cout << "update e_ot time: " << t << " bin # " << idx << " so far: " << e_onetick[idx] << endl;
      }
     else
     {
       if (Debug) cout << "rejected e_ot idx: " << idx << endl;
     }
   }
}

// Swallow and lareflex one tick bins
void NeuroBins::UpdateS_OT(DataLine& d_line)
{
   int t = d_line.GetTime();
   int j = t - i_start;  // start of swallow or lareflex
   int n_p = ip_time+ep_time;    // entire interval
   int n_b = n_p;
   int idx;
   int e_idx;

  idx = get_bin(j,n_p,n_b);
  if (Debug) cout << "swall one tick: idx: " << idx << " i_start: " << i_start << " ip_time: " << ip_time << " size: " << i_onetick.size() << endl;
  if (idx >= 0)
  {
     if (idx < int(i_onetick.size()))
        ++i_onetick[idx];
     else 
     {
        e_idx = idx - i_onetick.size();
        if (e_idx >= 0 && e_idx < int(e_onetick.size()))
           ++e_onetick[e_idx];
        else
        {
           if (e_idx < 0)
              cout << "swallow1 or lareflex onetick bin underflow: " << idx << " (not an error)" << endl;
           else 
              cout << "swallow1 or lareflex onetick bin overflow: " << idx <<" (not an error)" <<  endl;
        }
     }
  }
  else
  {
    if (Debug) cout << "rejected s_ot time: " << t << " idx: " << idx << endl;
  }
}


/* We have reached the end of a cycle, update the stats for all
   of the neuron bins, and set up for next cycle.
*/
void NeuroBins::UpdateBinRates(PeriodInfo& p_data)
{
   BinArrayIter iter;
   for (iter = i_bins.begin(); iter != i_bins.end(); ++iter)
   {
      iter->UpdateBinStats(p_data);
      iter->ResetSpike();
   }
   for (iter = e_bins.begin(); iter != e_bins.end(); ++iter)
   {
      iter->UpdateBinStats(p_data);
      iter->ResetSpike();
   }
}

/* This does more than it has to for equal sized bins (we need just one
   set of vars to hold some values for all the bins), but it lets us have 
   different widths for each bin if we want to.
*/
void NeuroBins::UpdateBinTimes()
{
   for (BinArrayIter iter = i_bins.begin(); iter != i_bins.end(); ++iter)
      iter->SetBinTime(NeuroBins::i_time, NeuroBins::secs_per_tick);
   for (BinArrayIter iter = e_bins.begin(); iter != e_bins.end(); ++iter)
      iter->SetBinTime(NeuroBins::e_time, NeuroBins::secs_per_tick);
}

// Calculate the stats for each bin in each phase
void NeuroBins::Stats(PeriodInfo& p_data)
{
   BinArrayIter iter;

   for (iter = i_bins.begin(); iter != i_bins.end(); ++iter)
   {
      iter->CalcBinStats(p_data);
      if (Debug) cout << "I Mean rate: " << iter->GetMeanRate() << " I sd " << iter->GetSD() << endl;
   }
   for (iter = e_bins.begin(); iter != e_bins.end(); ++iter)
   {
      iter->CalcBinStats(p_data);
      if (Debug) cout << "E Mean rate: " << iter->GetMeanRate()  << " E sd " << iter->GetSD() << endl;
   }
}

// First pass through files complete.  Time for second pass.
// Create the one tick bin counters we'll need for each channel for all periods.
void SetupForTicks(AllBins& bins, PerPeriod& per_p)
{
   PerKey p_key;
   for (AllBinsIter f_bins = bins.begin(); f_bins != bins.end(); ++f_bins)
   {
      for (FBinsIter fbi = f_bins->second.begin(); fbi != f_bins->second.end(); ++fbi)
      {
         int per = fbi->second.GetPeriod();
         if (per == CtlStims[CTL_SWALL1] || per == CtlStims[CTL_LAREFLEX])
            continue;
         p_key = make_pair(f_bins->first,per);
         PerPeriodIter p_iter = per_p.find(p_key);
         if (p_iter != per_p.end())
         {
            if (Debug) cout << "exp: " << f_bins->first << " per: " <<  fbi->second.GetPeriod()  << "  One ticks i: " << p_iter->second.i_short << "  e: " << p_iter->second.e_short << endl;
            fbi->second.CreateOneTickBins(p_iter->second.i_short, p_iter->second.e_short);
         }
           // fail to find not an error, some periods do not have a cell
      }
   }
}

// We have reached end of file.  Calculate the stats for the bins in all files and 
// all periods
void UpdateStats(AllBinsIter& f_bins, PeriodInfo &p_data)
{
   int period;
   for (FBinsIter fbi = f_bins->second.begin(); fbi != f_bins->second.end(); ++fbi)
   {
      period = fbi->first.second;
      if (period == p_data.period)
         fbi->second.Stats(p_data);
   }
}

/* Find and return largest mean value in bins for point.  Used for scaling output.
*/
double NeuroBins::GetBinPeak()
{
   BinArrayIter it;
   double max = 0;
   double curr;

   for (it = i_bins.begin(); it != i_bins.end(); ++it)
   {
      curr = it->GetMeanRate();
      if (curr > max)
         max = curr;
   }
   for (it = e_bins.begin(); it != e_bins.end(); ++it)
   {
      curr = it->GetMeanRate();
      if (curr > max)
         max = curr;
   }
   return max;
}

/* Return sum of mean rates of all the bins for a point (AKA Area)
*/
double NeuroBins::GetBinMeans()
{
   BinArrayIter it;
   double mean = 0;

   if (Debug) cout << "chan: " << chan << endl;
   for (it = i_bins.begin(); it != i_bins.end(); ++it)
   {
      mean += it->GetMeanRate();
      if (Debug) cout << " mean rate:  " <<  it->GetMeanRate() << " mean: " << mean << endl;
   }
   for (it = e_bins.begin(); it != e_bins.end(); ++it)
   {
      mean += it->GetMeanRate();
      if (Debug) cout << " mean rate:  " <<  it->GetMeanRate() << " mean: " << mean << endl;
   }
   return mean;
}


// debug function to validate scale for area normalization
void NeuroBins::CheckScale(double scale, double target)
{
   BinArrayIter it;
   double mean = 0, curr = 0;

   cout << "chan: " << chan << endl;
   if (Norm == MEAN)
   { 
      for (it = i_bins.begin(); it != i_bins.end(); ++it)
      {
         mean += it->GetMeanRate() * scale;
         cout << " mean rate:  " <<  it->GetMeanRate() << " mean: " << mean << endl;
      }
      for (it = e_bins.begin(); it != e_bins.end(); ++it)
      {
         mean += it->GetMeanRate() * scale;
         cout << " mean rate:  " <<  it->GetMeanRate() << " mean: " << mean << endl;
      }
   }
   else if (Norm == UNIT)
   { 
      for (it = i_bins.begin(); it != i_bins.end(); ++it)
      {
         curr = it->GetMeanRate();
         if (mean < curr)
            mean = curr;
         cout << " mean rate:  " <<  curr << " mean: " << mean << endl;
      }
      for (it = e_bins.begin(); it != e_bins.end(); ++it)
      {
         curr = it->GetMeanRate();
         if (mean < curr)
            mean = curr;
         cout << " mean rate:  " <<  curr << " mean: " << mean << endl;
      }
      mean *= scale;
   }
   else
      cout << "Peak and raw normalizing not handled." << endl;

   cout << "checkscale: "<< scale << "   "  << mean << "  " << target << endl;
}


// Save the mean scaling value and then calculate 
// the scaled std error for each bin
void NeuroBins::SetScale(double scale)
{
   BinArrayIter it;
   for (it = i_bins.begin(); it != i_bins.end(); ++it)
   {
      it->SetScaledVars(scale);
      if (Debug) cout << "scaled std err: " << it->GetScaledStdErr() <<  " scaled mean rate: " << it->GetScaledMeanRate() << endl;
   }
   for (it = e_bins.begin(); it != e_bins.end(); ++it)
   {
      it->SetScaledVars(scale);
      if (Debug) cout << "scaled std err: " << it->GetScaledStdErr() <<  " scaled mean rate: " << it->GetScaledMeanRate() << endl;
   }
}


// count of all spikes in all bins in a cth
unsigned long NeuroBins::GetTotalSpikes()
{
   int tot = 0;
   for (BinArrayIter it = i_bins.begin(); it != i_bins.end(); ++it)
      tot += it->GetTotalSpikes();
   for (BinArrayIter it = e_bins.begin(); it != e_bins.end(); ++it)
      tot += it->GetTotalSpikes();
   if (Debug)  cout << "Total spikes: " << tot << endl;
   return tot;
}

void NeuroBins::PrintIBins()
{
  for (BinArrayIter it = i_bins.begin(); it != i_bins.end(); ++it)
    cout << "  " << "mean i rate: " << it->GetMeanRate()*secs_per_tick << "  total spikes: " << it->GetTotalSpikes() << endl;
}

void NeuroBins::PrintEBins() 
{
   for (BinArrayIter it = e_bins.begin(); it != e_bins.end(); ++it)
      cout << "  " << "mean e rate: " << it->GetMeanRate() << " avg: " << it->GetMeanRate()*secs_per_tick << "  total spikes: " << it->GetTotalSpikes() << endl;
}

void NeuroBins::PrintITickBins()
{
   cout << "I One Ticks" << " size is " << i_onetick.size() << endl;
   for (OneTickIter it = i_onetick.begin() ; it != i_onetick.end(); ++it)
      cout << *it << endl;
   cout << "I One Ticks End" << endl;
}

void NeuroBins::PrintETickBins()
{
   cout << "E One Ticks " <<" size is " << e_onetick.size() <<  endl;
   for (OneTickIter it = e_onetick.begin() ; it != e_onetick.end(); ++it)
      cout << *it << endl;
   cout << "E One Ticks End" << endl;
} 

// Write one ticks bins for this cell for Octave
void NeuroBins::WriteOneTicks(ofstream& out) 
{
   OneTickIter it;

   out << "% name: OneTicks" << endl;
   out << "% type: matrix" << endl;
   out << "% rows: 1" << endl;
   if (i_onetick.size() + e_onetick.size() == 0)
      out << "% columns: " << 2 << endl;
   else
      out << "% columns: " << i_onetick.size() + e_onetick.size() << endl;
   if (i_onetick.size())
   {
      for (it = i_onetick.begin(); it != i_onetick.end(); ++it)
         out << *it << " ";
      out << "-1 ";   // i/e boundary
   }
   else
      out << "0 " << "-1 ";

   if (e_onetick.size())
   {
      for (it = e_onetick.begin(); it != e_onetick.end()-1; ++it)
         out << *it << " ";
      out << *it << endl;
   }
   else
      out << " 0" << endl;
}

// Write one tick spike count bins for this cell to output file for R
void NeuroBins::WriteOneTickBins(ofstream& out, int seq) 
{
   OneTickIter it;

   out << seq << " ";
   if (i_onetick.size())
   {
      for (it = i_onetick.begin(); it != i_onetick.end(); ++it)
         out << *it << " ";
      out << -1 << " ";  // mark i and e boundary
   }
   else
      out << "0 -1 ";
   if (e_onetick.size())
   {
      for (it = e_onetick.begin(); it != e_onetick.end()-1; ++it)
         out << *it << " ";
      out << *it;
   }
   else
      out << "0";
}

// Tell octave how many sparse CTHs rejected.
void NeuroBins::WriteOSparse(ofstream& out)
{
   out << endl;
   out << "% name: NumSparse " << endl;
   out << "% type: scalar" << endl;
   out << NumSparse << endl;
   out << endl;

   out << "% name: SparseThresh " << endl;
   out << "% type: scalar" << endl;
   out << RejectThresh << endl;
   out << endl;
}

void NeuroBins::WriteOVersion(ofstream& out)
{
   out << "% START MARK CTH_VERSION" << endl;
   out << "% name: " << "CTH_VERSION" << endl;
   out << "% type: string" << endl;
   out << "% elements: " << 1 << endl;
   out << "% length: " << CTH_VERSION.length()  << endl;
   out << CTH_VERSION << endl;
   out << "% END MARK" << endl;
}

// Funcs to write a file octave can read
// This create a var with name format:
//   F_<filename>_<m or p>_<neuron chan>_<total spikes>_<period>_<seq #>
// The calls to write I, then E bins must immediately follow
void NeuroBins::WriteOctaveHeader(ofstream& out, string expname, string fname, int seq, int total, int period)
{
   string norm;

   if (Norm == MEAN)   // code type of scaling, mean or peak
      norm = "m";
   else if (Norm == PEAK)
      norm = "p";
   else if (Norm == UNIT)
      norm = "u";
   else
      norm = "n";      // none

    // start of block marker
   out << setfill('0') << "% START MARK F_" << fname << "_" << norm << "_" << chan << "_" <<  total << "_" << period << "_" << setw(5) << seq << endl;
     // var name
   out << setfill('0') << "% name: F_" << fname << "_" << norm << "_" << chan << "_" <<  total << "_" << period << "_" << setw(5) << seq << endl;
   out << "% type: scalar struct" << endl;  // make a struct var
   out << "% ndims: 2" << endl;
   out << " 1 1" << endl;                   // 1 row, 1 col
   out << CTH_FIELDS << endl;              // # fields in it

   out << "% name: " << "ExpFileName" << endl;  // experiment filename
   out << "% type: string" << endl;
   out << "% elements: " << 1 << endl;
   out << "% length: " << expname.length()  << endl;
   out << expname << endl;

   out << "% name: IsSparse " << endl;
   out << "% type: scalar" << endl;
   out << (is_sparse == true) << endl;
   out << endl;
}

// write octave format bin values, I then E (last one write the newline)
void NeuroBins::WriteOIBins(ofstream& out)
{
   out << "% name: NSpaceCoords" << endl;
   out << "% type: matrix" << endl;
   out << "% rows: 1" << endl;
   out << "% columns: " << num_e_bins + num_i_bins  << endl;
   for (BinArrayIter it = i_bins.begin(); it != i_bins.end(); ++it)
      out << " " << it->GetScaledMeanRate();
}

// Companion to above
void NeuroBins::WriteOEBins(ofstream& out)
{
   for (BinArrayIter it = e_bins.begin(); it != e_bins.end(); ++it)
      out << " " << it->GetScaledMeanRate();
   out << endl;
   out << endl;
}

// Write the  mean rate and scaled standard error for each bin 
// so the octave program to use for error bar plotting.
void NeuroBins::WriteOStats(ofstream& out, int seq)
{
   int bins = NumEBins+NumIBins;
   BinArrayIter it;

   out << "% name: Cthstat" << endl;
   out << "% type: matrix" << endl;
   out << "% rows: " << bins << endl;
   out << "% columns: " << 2  << endl;

   for (it = i_bins.begin(); it != i_bins.end(); ++it)
      out << " " << it->GetMeanRate() << " " << it->GetScaledStdErr() << endl;
   for (it = e_bins.begin(); it != e_bins.end(); ++it)
      out << " " << it->GetMeanRate() << " " << it->GetScaledStdErr() << endl;
   out << endl;
   out << endl;
}

// Write the stereotaxis coordinates for pts that have them.
// If no value, write (0,0,0);
void NeuroBins::WriteOCoords(ofstream& out)
{
   out << "% name: RealCoords" << endl;
   out << "% type: scalar struct" << endl;
   out << "% ndims: 2" << endl;
   out << "% 1 1" << endl;
   out << "% length 12" << endl;

   out << "% name: expname" << endl;
   out << "% type: string" << endl;
   out << "% elements: 1" << endl;
   out << "% length " << cell_coords.expname.length() << endl;
   out << cell_coords.expname << endl;

   out << "% name: expbasename" << endl;
   out << "% type: string" << endl;
   out << "% elements: 1" << endl;
   out << "% length " << cell_coords.expbasename.length() << endl;
   out << cell_coords.expbasename << endl;

   out << "% name: name" << endl;
   out << "% type: string" << endl;
   out << "% elements: 1" << endl;
   out << "% length " << cell_coords.name.length() << endl;
   out << cell_coords.name << endl;

   out << "% name: mchan" << endl;
   out << "% type: string" << endl;
   out << "% elements: 1" << endl;
   out << "% length " << cell_coords.mchan.length() << endl;
   out << cell_coords.mchan << endl;

   out << "% name: dchan" << endl;
   out << "% type: string" << endl;
   out << "% elements: 1" << endl;
   out << "% length " << cell_coords.dchan.length() << endl;
   out << cell_coords.dchan << endl;

   out << "% name: ref" << endl;
   out << "% type: string" << endl;
   out << "% elements: 1" << endl;
   out << "% length " << cell_coords.ref.length() << endl;
   out << cell_coords.ref << endl;

   out << "% name: rl" << endl;
   out << "% type: scalar" << endl;
   out << cell_coords.rl << endl;

   out << "% name: dp" << endl;
   out << "% type: scalar" << endl;
   out << cell_coords.dp << endl;

   out << "% name: ap" << endl;
   out << "% type: scalar" << endl;
   out << cell_coords.ap << endl;

   out << "% name: ap_atlas" << endl;
   out << "% type: scalar" << endl;
   out << cell_coords.ap_atlas << endl;

   out << "% name: rl_atlas" << endl;
   out << "% type: scalar" << endl;
   out << cell_coords.rl_atlas << endl;

   out << "% name: dp_atlas" << endl;
   out << "% type: scalar" << endl;
   out << cell_coords.dp_atlas << endl;
}

// Write the stats used to create the custom distance metric for current var
// so octave can use it to save subsets of clusters.  This is most correct when
// area (mean) normalization is used.  Not so much for other normalizations.
void NeuroBins::WriteDistStats(ofstream& out)
{
   BinArrayIter it;
   int unsigned binidx;
   double scale;

   out << "% name: MeanSclStdErr" << endl;
   out << "% type: matrix" << endl;
   out << "% rows: " << num_i_bins+num_e_bins << endl;
   out << "% columns: 2" << endl;

   scale = GetScale();
   for (binidx = 0; binidx < num_i_bins; ++binidx) 
   {
      out << " " << GetIMeanRate(binidx)*scale;
      out << " " << GetIScaledStdErr(binidx);
      out << endl;
   }
   for (binidx = 0; binidx < num_e_bins; ++binidx) 
   {
      out << " " << GetEMeanRate(binidx)*scale;
      out << " " << GetEScaledStdErr(binidx);
      out << endl;
   }
   out << endl;
}

void NeuroBins::WriteRVals(ofstream& out, CTHCurves &r_vals, string& key)
{
   CTHCurvesIter r_iter;
   r_iter = r_vals.find(key);
   string curve;
   bool do_fake = false;
   int fake;
   size_t num_nums;

   if (r_iter != r_vals.end())
   {
      curve = r_iter->second;
      istringstream nums(curve);
      num_nums = distance(istream_iterator<string>(nums), istream_iterator<string>()); 
      if (num_nums <= 4) // cth_curve.R will send 4 "0" chars if there is no curve data 
      {                  // (this fails if doing curves for 4 bin cths)
         do_fake = true;
         curve.clear();
      }
   }
   else
      do_fake = true;

   if (do_fake)
   {
      for (fake = 0 ; fake < NumEBins+NumIBins+2; ++fake) // bins + 3 dummy info nums
         curve += "0 ";
      curve += "0";
      num_nums = NumEBins+NumIBins+3;
   }

   out << "% name: Curve" << endl;
   out << "% type: matrix" << endl;
   out << "% rows: 1" << endl;
   out << "% columns: " << num_nums  << endl;
   out << curve << endl;
}

// We write per-point vars out in a block.  This marks them so the octave
// programs can later extract subsets of the file based on the first variable
// name, which must immediately follow the start mark.  This mostly insulates
// them from file format changes.
void NeuroBins::WriteOEndMark(ofstream& out)
{
   out << endl << "% END MARK" << endl;
}
void NeuroBins::WriteOExpName(ofstream& out, string fname, int seq)
{
   char name[1024];
   BinArrayIter it;

   out << "% START MARK EXP_NAME" << endl;
   sprintf(name,"%s_%05d","ExpName",seq);
   out << "% name: " << name << endl;
   out << "% type: string" << endl;
   out << "% elements: " << 1 << endl;
   out << "% length: " << fname.length()  << endl;
   out << fname << endl;
   out << endl << "% END MARK" << endl;
}

// Funcs to write generic csv file, assumes I first, E last
void NeuroBins::WriteCIBins(ofstream& out)
{
   BinArrayIter it;
   for (it = i_bins.begin(); it != i_bins.end()-1; ++it)
      out << it->GetScaledMeanRate() << ",";
   out << it->GetScaledMeanRate();
}
void NeuroBins::WriteCEBins(ofstream& out)
{
   for (BinArrayIter it = e_bins.begin(); it != e_bins.end(); ++it)
      out << it->GetScaledMeanRate() << ",";
   out << endl;
}

// Funcs to write csv file ggobi likes
void NeuroBins::WriteGCHeader(ofstream& out)
{
   unsigned int seq;
    // "bin0,bin1. . ."
   out << "\"" << "\"" <<  ",";
   for (seq = 1; seq < num_e_bins+num_i_bins; seq++)
      out << "\"" << "bin" << seq << "\"" <<  ",";
   out << "\"" <<  "bin" << seq << "\"" << endl;
}

void NeuroBins::WriteGCIBins(ofstream& out, int seq)
{
   BinArrayIter it = i_bins.begin(); 
   out << "\"" << seq << "\"" << "," << it->GetScaledMeanRate() << ",";
   ++it;
   for (; it != i_bins.end(); ++it)
      out << it->GetScaledMeanRate() << ",";
}

void NeuroBins::WriteGCEBins(ofstream& out)
{
   BinArrayIter it;
   for (it = e_bins.begin(); it != e_bins.end()-1; ++it)
      out << it->GetScaledMeanRate() << ",";
   out << it->GetScaledMeanRate();
   out << endl;
}


typedef struct sumry
{
   double mr;    // mean rate
   double sserr;  // scaled std err 
} SUMMARY;

typedef vector< vector<SUMMARY> > SummaryBins;

// Create a distance matrix using a non-euclidean distance metric.
// The octave program that reads this file may use this info.
void WriteDist(ofstream& out, AllBins& bins, P_LIST& periods, PerPeriod& p_data)
{
   VFilesIter   f_iter;
   AllBinsIter  f_bins;
   FBinsIter    fbi;
   BinArrayIter it;
   int size = 0;
   int row0, row1, row, col;
   int idx, binidx, summaryidx, total_bins;
   double mr0, mr1, sserr0, sserr1, d,subd, sumsqr,scale;
   total_bins = NumIBins+NumEBins;
   double per_bin = 0;
   pair <long int,long int> p;
   int curr_period;

      // figure matrix size to avoid lots of memory reallocing
   for (f_bins = bins.begin(); f_bins != bins.end(); ++f_bins)
   {
      for (fbi = f_bins->second.begin(); fbi != f_bins->second.end(); ++fbi)
      {
         curr_period = fbi->first.second;
         if (find(periods.begin(), periods.end(), curr_period) != periods.end())
         {
            if (!fbi->second.IsSparse())
               ++size;
         }
      }
   }
   ++size;  // add in slot for zero flat
   MatrixXd dist = MatrixXd::Zero(size,size);

     // since we are going to access this stuff a lot, make one
     // pass through the bins and summarize things into a set of
     // the values we'll use to calculate distance
   SummaryBins summary(size, vector<SUMMARY>(total_bins));

   for (idx = 0, f_bins = bins.begin(); f_bins != bins.end(); ++f_bins)
   {
      for (fbi = f_bins->second.begin(); fbi != f_bins->second.end(); ++fbi)
      {
         curr_period = fbi->first.second;
         if (find(periods.begin(), periods.end(), curr_period) != periods.end())
         {
            if (!fbi->second.IsSparse())
            {
               scale = fbi->second.GetScale();
               for (summaryidx = 0, binidx = 0; binidx < NumIBins; ++binidx,++summaryidx)
               {
                  summary[idx][summaryidx].mr = fbi->second.GetIMeanRate(binidx)*scale;
                  summary[idx][summaryidx].sserr = fbi->second.GetIScaledStdErr(binidx);
               }
               for (binidx = 0; binidx < NumEBins; ++binidx,++summaryidx)
               {
                  summary[idx][summaryidx].mr = fbi->second.GetEMeanRate(binidx)*scale;
                  summary[idx][summaryidx].sserr = fbi->second.GetEScaledStdErr(binidx);
               }
               ++idx;
                 // this assumes values will be the same for all included periods.
               if (!per_bin)
                  per_bin = maxCTH[CtlStims[curr_period]]/total_bins;
            }
         }
      }
   }

   for (summaryidx=0 ; summaryidx < total_bins; ++summaryidx) // add zero flat
   {
      summary[idx][summaryidx].mr = per_bin;
      summary[idx][summaryidx].sserr = 0;
   }

   // okay, do what we're here for
   for (row0 = 0; row0 < size-1; ++row0)
   {
      for (row1 = row0+1; row1 < size; ++row1)
      {
         sumsqr = 0;
         for (idx=0; idx < total_bins; idx++)
         {
            mr0   = summary[row0][idx].mr;
            sserr0 = summary[row0][idx].sserr;
            mr1   = summary[row1][idx].mr;
            sserr1 = summary[row1][idx].sserr;
            if (sserr0 == 0 && sserr1 == 0)
               continue;
                  // our distance metric
            subd  = (mr0 - mr1) / sqrt(sserr0*sserr0 + sserr1*sserr1);
// for testing, do euclidean
//            subd  = (mr0 - mr1);
            sumsqr += subd*subd;
         }
         d = sqrt(sumsqr);
         dist(row0,row1) = d;  // add to matrix
         dist(row1,row0) = d;
      }
   }
//   cout << dist << endl;   // print entire matrix (might be big)

   out << "% name: meandist" << endl;  // octave program expects this name
   out << "% type: matrix" << endl;
   out << "% rows: " << size << endl;
   out << "% columns: " << size << endl;

   for (row = 0 ; row < size; ++row)
   {
      for (col = 0; col < size; ++col)
         out << " " << dist(row,col);
      out << endl; 
   }
}


// Create a cth that has the same mean/peak/nothing as the largest cth
// (This probably does not work with anything but mean scaling)
void WriteZeroFlat(ofstream& out,int seq,double max_cth)
{
   int    bins = NumEBins+NumIBins;
   int    zf;
   double per_bin = max_cth / bins;
   string noname("Zero Flat");

   out << setfill('0') << "% START MARK ZeroFlat" << "_" << "m" << "_" << "xxx" << "_" <<  (int) max_cth << "_" << 0 << "_" << setw(5) << seq << endl;
   out << setfill('0') << "% name: ZeroFlat" << "_" << "m" << "_" << "xxx" << "_" <<  (int) max_cth << "_" << 0 << "_" << setw(5) << seq << endl;
   out << "% type: scalar struct" << endl;  // make a struct var
   out << "% ndims: 2" << endl;
   out << " 1 1" << endl;                   // 1 row, 1 col
   out << CTH_FIELDS << endl;              // # fields in it

   out << "% name: " << "ExpFileName" << endl;  // experiment filename
   out << "% type: string" << endl;
   out << "% elements: " << 1 << endl;
   out << "% length: " << noname.length()  << endl;
   out << noname << endl << endl;

   out << "% name: IsSparse " << endl;
   out << "% type: scalar" << endl;
   out << 0 << endl;
   out << endl;

   out << "% name: NSpaceCoords" << endl;
   out << "% type: matrix" << endl;
   out << "% rows: 1" << endl;
   out << "% columns: " << bins  << endl;
   for (int col=0; col < bins; col++)
      out << " " << per_bin;
   out << endl << endl;

   out << "% name: Cthstat" << endl;
   out << "% type: matrix" << endl;
   out << "% rows: " << bins << endl;
   out << "% columns: " << 2  << endl;

   for (zf=0; zf < bins; ++zf)
      out << " " << per_bin << " " << 0.0 << endl;
   out << endl << endl; 

   out << "% name: RealCoords" << endl;
   out << "% type: scalar struct" << endl;
   out << "% ndims: 2" << endl;
   out << "% 1 1" << endl;
   out << "% length 12" << endl;

   out << "% name: expname" << endl;
   out << "% type: string" << endl;
   out << "% elements: 1" << endl;
   out << "% length " << strlen("NOEXP") << endl;
   out << "NOEXP" << endl;

   out << "% name: expbasename" << endl;
   out << "% type: string" << endl;
   out << "% elements: 1" << endl;
   out << "% length " << strlen("NOBASE") << endl;
   out << "NOBASE" << endl;

   out << "% name: name" << endl;
   out << "% type: string" << endl;
   out << "% elements: 1" << endl;
   out << "% length " << strlen("zeroflat") << endl;
   out << "zeroflat" << endl << endl << endl;

   out << "% name: mchan" << endl;
   out << "% type: string" << endl;
   out << "% elements: 1" << endl;
   out << "% length " << strlen("NOMCHAN") << endl;
   out << "NOMCHAN" << endl;

   out << "% name: dchan" << endl;
   out << "% type: string" << endl;
   out << "% elements: 1" << endl;
   out << "% length " << strlen("NODCHAN") << endl;
   out << "NODCHAN" << endl;

   out << "% name: ref" << endl;
   out << "% type: string" << endl;
   out << "% elements: 1" << endl;
   out << "% length " << strlen("NOREF") << endl;
   out << "NOREF" << endl;

   out << "% name: rl" << endl;
   out << "% type: scalar" << endl;
   out << 0.0 << endl << endl << endl;

   out << "% name: dp" << endl;
   out << "% type: scalar" << endl;
   out << 0.0 << endl << endl << endl;

   out << "% name: ap" << endl;
   out << "% type: scalar" << endl;
   out << 0.0 << endl << endl << endl;

   out << "% name: ap_atlas" << endl;
   out << "% type: scalar" << endl;
   out << 0.0 << endl;

   out << "% name: rl_atlas" << endl;
   out << "% type: scalar" << endl;
   out << 0.0 << endl;

   out << "% name: dp_atlas" << endl;
   out << "% type: scalar" << endl;
   out << 0.0 << endl;

   out << "% name: MeanSclStdErr" << endl;
   out << "% type: matrix" << endl;
   out << "% rows: " << bins << endl;
   out << "% columns: " << 2  << endl;
   for (zf=0; zf < bins; ++zf)
      out << " " << per_bin << " " << 0.0 << endl;
   out << endl << endl; 

   const int oneticks = 4001;
   out << "% name: OneTicks" << endl;   // fake but rather dense zero one tick
   out << "% type: matrix" << endl;
   out << "% rows: 1" << endl;
   out << "% columns: " << oneticks << endl;
   int i;
   for (i = 0; i < oneticks/2; ++i)
      out << "1 ";
   out << "-1 ";   // fake phase boundary
   for (; i < oneticks-1; ++i)
      out << "1 ";
   out << "1" << endl;

   out << "% name: Curve" << endl;
   out << "% type: matrix" << endl;
   out << "% rows: 1" << endl;
   out << "% columns: " << bins+3  << endl;
   for (int col=0; col < bins; col++)
      out << " " << per_bin;
   out << " 1 0 1";
   out << endl << endl;

   out << endl << "% END MARK" << endl;
}

// Utility to compose the unique identifier for each cth.
// It is:  experiment file name + channel # + period
string MakeVarName(const string& exp, int chan, int period)
{
   stringstream name;
   name << exp << "_" << chan << "_" <<  period;
   return string(name.str());
}


/* Optionally write out a file that contains something the 
   R curve fitting program can read in.
   This will be a matrix, where each row is a cth with a leading seq #
   Optionally Start the R script to use this.  
   Wait for it to complete (can take a REAL long time).
   Optionally read the results back in and add them to the output file.
   Sometimes it is best to use one of the compute servers to do the
   work, so this allows you to also just write out the input file R
   wants, copy it and the cth file to another machine, do the curve
   fitting, then copy the results back and re-run this program to read
   the R output file and add it to this program's output.
   If we have no input from R, we still need to return something in the
   r_vals param.  Create vectors of 0 for each cth
*/
void GetRVals(AllBins& bins, string in_name, string out_name, CTHCurves& r_vals)
{
   unsigned long rows = 1;
   AllBinsIter  f_bins;
   FBinsIter fbi;

   r_vals.clear();  // ensure no pts

   if (DoCurves || WrCurves)  // shall we even do this?
   {
      ofstream out(in_name.c_str());
      if (!out)
      {
         cout << "Could not open " << in_name << ".  Curve file not created" << endl;
         return;
      }
      cout << "Creating cth_curve.R program input file with CTHs from all periods." << endl;
      for (f_bins = bins.begin(); f_bins != bins.end(); ++f_bins)
      {
         string exp(f_bins->first);
         for (fbi = f_bins->second.begin(); fbi != f_bins->second.end(); ++fbi, ++rows)
         {
            int chan = fbi->first.first;
            int per = fbi->first.second;
            out << MakeVarName(exp,chan,per) << endl;
            fbi->second.WriteOneTickBins(out,rows);
            out << endl;
         }
      }
      out.close();
      WrCurves = false; // we put all pts in output file, so only need one instance.
   }


   if (DoCurves)  // do calcs?
   {
      string cmd = "cth_curve.R";    // hopefully is on the path
      cmd += " " + in_name + " " + out_name;
      cout << "Starting " << cmd << " (this is a long computation, many hours, possibly days.)" << endl;
      int call_res = system(cmd.c_str());
      cout << "Back from cth_curve.R, result is " << call_res << endl;
   }

    // Since it takes so long to calculate these, if the R output file already
    // exists, assume we want to reuse it, so read it in.  No check is made for
    // valid input (e,g. # of pts, order of pts the same, etc.)
   if (DoCurves || RdCurves)
   {
      cout << "Read cth_curve.R program output file for CTH curves." << endl;
      ifstream in(out_name.c_str());
      if (in)
      {
         string key, line;
         while (!in.eof())
         {
            getline(in,key);
            key.erase(key.find_last_not_of(" ")+1);
            getline(in,line);
            if (key.empty() || line.empty())
               continue;
            r_vals.insert(make_pair(key,line));
         }
         in.close();
      }
      else
         cout << "Could not open " << out_name << endl;
   }
}


void PrintBins(AllBinsIter& f_bins)
{
   cout << " #########  BINS  ########" << endl;
   for (FBinsIter fbi = f_bins->second.begin(); fbi != f_bins->second.end(); ++fbi)
   {
      fbi->second.PrintIBins();
      fbi->second.PrintEBins();
   }
}

void PrintTickBins(AllBinsIter& f_bins)
{
   for (FBinsIter fbi = f_bins->second.begin(); fbi != f_bins->second.end(); ++fbi)
   {
      fbi->second.PrintITickBins();
      fbi->second.PrintETickBins();
   }
}

// new i phase begins
void StartNewCycle(AllBinsIter& f_bins, PeriodInfo &per_p)
{
   if (Debug) cout << "********** start new cycle ***********" << endl;

   for (FBinsIter fbi = f_bins->second.begin(); fbi != f_bins->second.end(); ++fbi)
   {
      if (fbi->first.second == per_p.period)
      {
         fbi->second.UpdateBinTimes();
           // saved E's from the last phase (which was E) belong in 
           // this (new) phase's I bins for rate bins.
         while (fbi->second.HaveI())
         {
            fbi->second.UpdateI(fbi->second.GetI());
            fbi->second.RemoveI();
            if (Debug) cout << "New Cycle: pick up saved E for new I" << endl;
         }
      }
   }
}

/* Boundary events saved from the first phase (which was i)
   belong in the next (E) phase's E bin(s) and are now put there.
*/
void StartEPhase(AllBinsIter& f_bins,PeriodInfo& p_data) 
{
   for (FBinsIter fbi = f_bins->second.begin(); fbi != f_bins->second.end(); ++fbi)
   {
      if (fbi->first.second == p_data.period)
      {
         while (fbi->second.HaveE())
         {
            fbi->second.UpdateE(fbi->second.GetE());
            fbi->second.RemoveE();
            if (Debug) cout << "picking up saved I for new E" << endl;
         }
      }
   }
}

// new one tick i phase begins
void StartNew_OT_Cycle(AllBinsIter& f_bins,PeriodInfo& p_data)
{
   for (FBinsIter fbi = f_bins->second.begin(); fbi != f_bins->second.end(); ++fbi)
   {
        // saved E's from the last phase (which was E) belong in 
        // this (new) phase's I bins for rate bins.
      if (fbi->first.second == p_data.period)
      {
         while (fbi->second.HaveI_OT())
         {
            fbi->second.UpdateI_OT(fbi->second.GetI_OT());
            fbi->second.RemoveI_OT();
            if (Debug) cout << "picking up saved E_OTfor new I_OT" << endl;
         }
      }
   }
}

/* Boundary events saved from the first phase (which was i)
   belong in the next (E) phase's E bin(s) and are now put there.
*/
void StartE_OT_Phase(AllBinsIter& f_bins, PeriodInfo& p_data)
{
   for (FBinsIter fbi = f_bins->second.begin(); fbi != f_bins->second.end(); ++fbi)
   {
      if (fbi->first.second == p_data.period)
      {
         while (fbi->second.HaveE_OT())
         {
            fbi->second.UpdateE_OT(fbi->second.GetE_OT());
            fbi->second.RemoveE_OT();
            if (Debug) cout << "picking up saved I_OT for new E_OT" << endl;
         }
      }
   }
}

void EndCurrCycle(AllBinsIter& f_bins,PeriodInfo& p_data)
{
   if (Debug) cout << "********** end of current cycle ***********" << endl;
   for (FBinsIter fbi = f_bins->second.begin(); fbi != f_bins->second.end(); ++fbi)
   {
      if (fbi->first.second == p_data.period)
         fbi->second.UpdateBinRates(p_data);
   }
}

/* If we have saved events waiting to be carried over to the next bin
   and we decide to skip the next bin, we have to toss the saved events.  This
   is because it is assumed that the next event time for the saved events is
   the same as the time of the current phase.  When we skip ahead, the saved
   event is much less than the new interval time, resulting in a negative bin
   index.  This does not happen often, so we are not affecting the results very
   much.
*/
void PurgeAllSaved(AllBinsIter& f_bins, PeriodInfo& p_data)
{
   if (Debug) cout << "Purging " << endl;

   for (FBinsIter fbi = f_bins->second.begin(); fbi != f_bins->second.end(); ++fbi)
   {
      if (fbi->first.second == p_data.period)
      {
         while (fbi->second.HaveI())
            fbi->second.RemoveI();
         while (fbi->second.HaveE())
            fbi->second.RemoveE();
         while (fbi->second.HaveI_OT())
            fbi->second.RemoveI_OT();
         while (fbi->second.HaveE_OT())
            fbi->second.RemoveE_OT();
      }
   }
}


/* We need to use a perl utility to process the info xls file and
   calculate the stereotaxic coords for us, so they will agree with the atlas
   program (which also uses the perl program) The cth_cells.pl file has been
   modified to include the chan number in the tmp file.
   This assumes the cth_cells.pl utility is on the current PATH.
*/
int GetAtlasCoords(string fname, AtlasCoords& a_coords)
{
   string tmp("cth_cells.tmp");
   int    cnum, matches;
   string line;
   ACoords coords;

   string cth_cells("cth_cells.pl " + fname + " > " + tmp);
   int ok = system(cth_cells.c_str());

   if (ok)
   {
      cout << "Error: The program " << cth_cells.c_str() << " returned an error." << endl;
      perror("Error is");
      return 1;
   }

   ifstream cell_file(tmp);
   if (cell_file)
   {
      while (!cell_file.eof())
      {
         getline(cell_file,line);
         matches = sscanf(line.c_str(),"%lf %lf %lf %d",&coords.ap, &coords.rl, &coords.dp, &cnum);
         if (matches == 4)
         {
            if (coords.ap == 0.0 && coords.rl == 0.0 && coords.dp == 0.0)
               cout << "Have cell at origin" << endl;
            a_coords.emplace(make_pair(cnum,coords));
         }
      }
      return 0;
   }
   else
   {
      cout << "Could not open cth_cells.pl output file: " << tmp << endl;
      return 1;
   }
}



// Before we start what could be lengthy processing, make sure all the files
// we need exist. From time to time, .edt files are updated and the old ones
// go away.
bool DoFilesExist(string fname)
{
   int line_count = 0;
   unsigned int worksheet_index;
   const void *handle = nullptr;
   int ret;
   unsigned int max_worksheet;
   unsigned int rows;
   unsigned short columns;
   unsigned int row;
   unsigned short col;
   const char *worksheet_name;
   string    currline;
   bool      found = false;
   bool      error = false;
   string    cell_text;
   string    tmpnum;
   string    pathfile;
   size_t    beg, end; 
   bool goodfile = false;
   size_t f_start = -1;
   size_t f_end = 0;
   size_t f_beg;
   FreeXL_CellValue cell;
   stringstream file_list;
   bool      is_good = true;

   cout << "Making sure files exist..." << endl;
   while (is_good)
   {
         /* opening the .XLS file [Workbook] */
      ret = freexl_open(fname.c_str(), &handle);
      if (ret != FREEXL_OK)
      {
         cout << "Could not open: " << fname << endl;
         is_good = false;
         break;
      }

         /* querying BIFF Worksheet entries */
      ret = freexl_get_info (handle, FREEXL_BIFF_SHEET_COUNT, &max_worksheet);
      if (ret != FREEXL_OK)
      {
         cout << "GET-INFO [FREEXL_BIFF_SHEET_COUNT] Error: " << ret << " in " << fname << endl;
         is_good = false;
         break;
      }
        // look for "all cells" or "Info File" page
      for (worksheet_index = 0; worksheet_index < max_worksheet && !error; worksheet_index++)
      {
         ret = freexl_get_worksheet_name (handle, worksheet_index, &worksheet_name);
         if (ret != FREEXL_OK)
         {
            cout << "GET-WORKSHEET-NAME Error:  " << ret << endl;
            is_good = false;
            break;;
         }
         if (worksheet_name != nullptr)
         { 
            if ((strcasestr(worksheet_name,FILE_XLS1) != nullptr)  || (strcasestr(worksheet_name,FILE_XLS2) != nullptr))
            {
               found = true;
               break;
            }
         }
      }
      if (!found)
      {
         cout << "Unknown spreadsheet contents.  Only file list spreadsheets are handled by this program." << endl;
         is_good = false;
         break;;
      }
         // select worksheet and get info
      ret = freexl_select_active_worksheet(handle, worksheet_index);
      if (ret != FREEXL_OK)
      {
         cout <<  "SELECT-ACTIVE_WORKSHEET Error: " << ret << endl;
         is_good = false;
         break;;
      }

      ret = freexl_worksheet_dimensions(handle, &rows, &columns);
      if (ret != FREEXL_OK)
      {
         cout << "WORKSHEET-DIMENSIONS Error: " << ret << endl;
         is_good = false;
         break;;
      }

      // No easy way to know # of columns, so look for header row and parse it.
      // Expect first col to be "path" and last to be "notes".  If there are
      // fewer or more columns than we expect, reject the spread
      for (row = 0; row < rows && !goodfile; ++row)
      {
         for (col = 0; col < columns && !goodfile; ++col) 
         {
            ret = freexl_get_cell_value(handle, row, col, &cell);
            if (ret != FREEXL_OK)
            {
               cout <<"CELL-VALUE-ERROR (r=" << row << "c="<< col << "err: " <<  ret << endl;
               continue;
            }
            if (cell.type == FREEXL_CELL_TEXT || cell.type == FREEXL_CELL_SST_TEXT)
            {
               cell_text = cell.value.text_value;
               f_beg = cell_text.find_first_not_of("' '\t");
               f_end = cell_text.find_last_not_of("' '\t\n");
               if (f_beg != string::npos)
               {
                  cell_text=cell_text.substr(f_beg,f_end-f_beg+1);
                  if (cell_text.compare(FILE_HEADER_START) == 0)
                     f_start = col;
                  else if (cell_text.compare(FILE_HEADER_END) == 0)
                  {
                     f_end = col;
                     goodfile = true;
                  }
               }
            }
         }
      }

      if (!goodfile || (f_end-f_start+1) != CURR_FORMAT) // we don't know what this is, toss it
      {
         cout << "The spreadsheet " << fname << " is too old or too new for this program to recognize." << endl << "If it is newer, this program will need to be modified to use the new fields ." << endl;
         break;
      }

      for (row = 0; row < rows && ! error; row++)
      {
         FreeXL_CellValue cell;
         stringstream fields;
         vector<int> ftype(columns);
         fill(ftype.begin(),ftype.end(),0);
         currline.clear();

         for (col = 0; col < columns && !error; col++) 
         {
            ret = freexl_get_cell_value(handle, row, col, &cell);
            if (ret != FREEXL_OK)
            {
               cout <<"CELL-VALUE-ERROR (r=" << row << "c="<< col << "err: " <<  ret << endl;
//               error = true;
//               break;;
            }
             // Read in the current line and keep track of the field types.
             // Later, we will apply some heuristics to see if the line
             // matches any signatures of lines we want to keep
            switch (cell.type)
            {
               case FREEXL_CELL_INT:
                   fields << cell.value.int_value << "\t";  
                   ftype[col] = FREEXL_CELL_INT;
                   break;
               case FREEXL_CELL_DOUBLE:
                   fields << setprecision(15) << cell.value.double_value << "\t";
                   ftype[col] = FREEXL_CELL_DOUBLE;
                   break;
               case FREEXL_CELL_TEXT:
               case FREEXL_CELL_SST_TEXT:
                   cell_text = cell.value.text_value;
                   beg = cell_text.find_first_not_of("' '\t");
                   end = cell_text.find_last_not_of("' '\t\n");
                   if (beg != string::npos)
                      cell_text=cell_text.substr(beg,end-beg+1);
                   else
                      cell_text.clear(); // empty or all ws
                   fields << cell_text  << "\t";
                   ftype[col] = FREEXL_CELL_TEXT;
                   break;
               case FREEXL_CELL_DATE:
               case FREEXL_CELL_DATETIME:
               case FREEXL_CELL_TIME:
                   fields << cell.value.text_value << "\t";  
                   ftype[col] = FREEXL_CELL_DATE;
                   break;
               case FREEXL_CELL_NULL:
               default:
                   ftype[col] = FREEXL_CELL_NULL;
                   fields << "\t";
                   break;
            }
         }
         fields << endl;

         if (Debug) cout << fields.str();
         if (Debug) cout << ftype[0] << " " << ftype[1] << " " << ftype[2] << " " << ftype[3] << endl;
         if (ftype[0] == FREEXL_CELL_TEXT && fields.str()[0] == '/' &&
             ftype[1] == FREEXL_CELL_TEXT && (ftype[2] == FREEXL_CELL_DOUBLE || ftype[2] == FREEXL_CELL_INT))
         {
            file_list << fields.str();
            ++line_count;
         }
      }
      break;   // done
   }

   while (!file_list.eof())
   {
      getline(file_list,currline);
      if (currline.empty())  // ignore blank lines, #comments, lines with just whitespace
         continue;
      size_t pos = currline.find_first_not_of("' '\t");
      if (pos == string::npos || currline[pos] == '#')
         continue;

      FileFacts facts;
      istringstream fields(currline);
      getline(fields,facts.path,'\t');
      getline(fields,facts.fname,'\t');
      getline(fields,tmpnum,'\t');       // skip
      getline(fields,tmpnum,'\t');
      getline(fields,tmpnum,'\t');
      getline(fields,facts.info_path,'\t');
      getline(fields,facts.info_file,'\t');
      getline(fields,facts.changes_path,'\t');
      getline(fields,facts.changes_file,'\t');
      getline(fields,facts.notes,'\n');

        // check existence, some optional
      pathfile = facts.path+"/"+facts.fname;  
      ifstream file;
      file.open(pathfile);
      if (!file)
      {
         cout << pathfile << " is missing." <<  endl;
         is_good = false;
      }
      else
         file.close();
      if (facts.info_path.length() && facts.info_file.length())
      {
         pathfile = facts.info_path+"/"+facts.info_file;
         file.open(pathfile);
         if (!file)
         {
            cout << pathfile << " is missing." <<  endl;
            is_good = false;
         }
         else
            file.close();
      }
      if (facts.changes_path.length() && facts.changes_file.length())
      {
         pathfile = facts.changes_path+"/"+facts.changes_file;
         file.open(pathfile);
         if (!file)
         {
            cout << pathfile << " is missing." <<  endl;
            is_good = false;
         }
         else
            file.close();
      }
   }

   if (handle) /* closing the .XLS file [Workbook] */
   {
      ret = freexl_close(handle);
      if (ret != FREEXL_OK)
         fprintf (stderr, "CLOSE ERROR: %d\n", ret);
   }

   if (is_good && line_count == 0)
   {
      cout << "The file: " << fname << " does not have any rows with data." << endl;
      is_good = false;
   }

   return is_good;
}


/* Read the excel spreadsheet and get info about where various files live
   and how to read them.  This uses the libfreexl-dev package.  This reads two
   types of xls files that has info we need.  The first is the file that
   contains file paths and names of .edt, etc. files and where the second type
   of .xls files are, info about format of data in these files, such as E and I
   channels, info about various control/stim segments, and so on.
   The second type of xls files contain the stereotaxic coordinates
   of the points in the .edt, etc., files.
*/
int GetRowsFromSprd(string fname, XLS_TYPE filetype, stringstream& sprd_data)
{
   int line_count = 0;
   unsigned int worksheet_index;
   const void *handle = nullptr;
   int ret;
   unsigned int max_worksheet;
   unsigned int rows;
   unsigned short columns;
   unsigned int row;
   unsigned short col;
   const char *worksheet_name;
   string    currline;
   bool      found = false;
   bool      error = false;
   string    cell_text;
   size_t    beg, end; 

   while (!error)
   {
         /* opening the .XLS file [Workbook] */
      ret = freexl_open(fname.c_str(), &handle);
      if (ret != FREEXL_OK)
      {
         cout << "Could not open: " << fname << endl;
         return line_count;
      }
      cout << "Opened spreadsheet " << fname << endl;

         /* querying BIFF Worksheet entries */
      ret = freexl_get_info (handle, FREEXL_BIFF_SHEET_COUNT, &max_worksheet);
      if (ret != FREEXL_OK)
      {
         cout << "GET-INFO [FREEXL_BIFF_SHEET_COUNT] Error: " << ret << " in " << fname << endl;
         error = true;
         break;;
      }
        // look for "all cells" or "Info File" page
      for (worksheet_index = 0; worksheet_index < max_worksheet && !error; worksheet_index++)
      {
         ret = freexl_get_worksheet_name (handle, worksheet_index, &worksheet_name);
         if (ret != FREEXL_OK)
         {
            cout << "GET-WORKSHEET-NAME Error:  " << ret << endl;
            error = true;
            break;;
         }
         if (worksheet_name != nullptr)
         { 
            if ( ((strcasestr(worksheet_name,FILE_XLS1) != nullptr && filetype == FILELIST) ||
                  (strcasestr(worksheet_name,FILE_XLS2) != nullptr && filetype == FILELIST)) ||
                 ( (strcasestr(worksheet_name,INFO_XLS1) != nullptr && 
                    strcasestr(worksheet_name,INFO_XLS2) == nullptr) && filetype == INFO))
            {
               found = true;
               break;
            }
         }
      }
      if (!found)
      {
         cout << "Unknown spreadsheet contents.  Only file list or info spreadsheets are handled by this program." << endl;
         error = true;
         break;;
      }
         // select worksheet and get info
      ret = freexl_select_active_worksheet(handle, worksheet_index);
      if (ret != FREEXL_OK)
      {
         cout <<  "SELECT-ACTIVE_WORKSHEET Error: " << ret << endl;
         error = true;
         break;;
      }

      ret = freexl_worksheet_dimensions(handle, &rows, &columns);
      if (ret != FREEXL_OK)
      {
         cout << "WORKSHEET-DIMENSIONS Error: " << ret << endl;
         error = true;
         break;;
      }

      // If a FILE type, no easy way to know # of columns, so look for header row
      // and parse it.  Expect first col to be "path" and last to be "notes".
      // If there are fewer or more columns than we expect, reject the spread
      if (filetype == FILELIST) 
      {
         bool goodfile = false;
         int start = -1;
         int end = 0;
         FreeXL_CellValue cell;
         for (row = 0; row < rows && !goodfile; ++row)
         {
            for (col = 0; col < columns && !goodfile; ++col) 
            {
               ret = freexl_get_cell_value(handle, row, col, &cell);
               if (ret != FREEXL_OK)
               {
                  cout <<"CELL-VALUE-ERROR (r=" << row << "c="<< col << "err: " <<  ret << endl;
                  continue;
               }
               if (cell.type == FREEXL_CELL_TEXT || cell.type == FREEXL_CELL_SST_TEXT)
               {
                   cell_text = cell.value.text_value;
                   beg = cell_text.find_first_not_of("' '\t");
                   end = cell_text.find_last_not_of("' '\t\n");
                   if (beg != string::npos)
                   {
                      cell_text=cell_text.substr(beg,end-beg+1);
                      if (cell_text.compare(FILE_HEADER_START) == 0)
                         start = col;
                      else if (cell_text.compare(FILE_HEADER_END) == 0)
                      {
                         end = col;
                         goodfile = true;
                      }
                   }
                }
            }
         }

         if (!goodfile || (end-start+1) != CURR_FORMAT) // we don't know what this is, toss it
         {
            cout << "The spreadsheet " << fname << " is too old or too new for this program to recognize." << endl << "If it is newer, this program will need to be modified to use the new fields ." << endl;
            return line_count;
         }
      }

      for (row = 0; row < rows && ! error; row++)
      {
         FreeXL_CellValue cell;
         stringstream fields;
         vector<int> ftype(columns);
         fill(ftype.begin(),ftype.end(),0);
         currline.clear();

         for (col = 0; col < columns && !error; col++) 
         {
            ret = freexl_get_cell_value(handle, row, col, &cell);
            if (ret != FREEXL_OK)
            {
               cout <<"CELL-VALUE-ERROR (r=" << row << "c="<< col << "err: " <<  ret << endl;
//               error = true;
//               break;;
            }
             // Read in the current line and keep track of the field types.
             // Later, we will apply some heuristics to see if the line
             // matches any signatures of lines we want to keep
            switch (cell.type)
            {
               case FREEXL_CELL_INT:
                   fields << cell.value.int_value << "\t";  
                   ftype[col] = FREEXL_CELL_INT;
                   break;
               case FREEXL_CELL_DOUBLE:
                   fields << setprecision(15) << cell.value.double_value << "\t";
                   ftype[col] = FREEXL_CELL_DOUBLE;
                   break;
               case FREEXL_CELL_TEXT:
               case FREEXL_CELL_SST_TEXT:
                   cell_text = cell.value.text_value;
                   beg = cell_text.find_first_not_of("' '\t");
                   end = cell_text.find_last_not_of("' '\t\n");
                   if (beg != string::npos)
                      cell_text=cell_text.substr(beg,end-beg+1);
                   else
                      cell_text.clear(); // empty or all ws
                   fields << cell_text  << "\t";
                   ftype[col] = FREEXL_CELL_TEXT;
                   break;
               case FREEXL_CELL_DATE:
               case FREEXL_CELL_DATETIME:
               case FREEXL_CELL_TIME:
                   fields << cell.value.text_value << "\t";  
                   ftype[col] = FREEXL_CELL_DATE;
                   break;
               case FREEXL_CELL_NULL:
               default:
                   cout << "cell case: " << cell.type << endl;
                   ftype[col] = FREEXL_CELL_NULL;
                   fields << "\t";
                   break;
            }
         }
         fields << endl;

         if (Debug) cout << fields.str();
         if (Debug) cout << ftype[0] << " " << ftype[1] << " " << ftype[2] << " " << ftype[3] << endl;
         if (filetype == FILELIST) 
         {
            if (ftype[0] == FREEXL_CELL_TEXT && fields.str()[0] == '/' &&
                ftype[1] == FREEXL_CELL_TEXT && (ftype[2] == FREEXL_CELL_DOUBLE || ftype[2] == FREEXL_CELL_INT))
            {
               sprd_data << fields.str();
               ++line_count;
            }
         }
         else if (filetype == INFO)
         {
            if (ftype[0] == FREEXL_CELL_TEXT && 
               (ftype[1] == FREEXL_CELL_INT || ftype[1] == FREEXL_CELL_DOUBLE) && 
               (ftype[2] == FREEXL_CELL_INT || ftype[2] == FREEXL_CELL_DOUBLE) && 
               (ftype[3] == FREEXL_CELL_INT || ftype[3] == FREEXL_CELL_DOUBLE) )
            {
               sprd_data << fields.str();
               ++line_count;
            }
         }
      }
      break;   // done
   }

   if (Debug) cout << sprd_data.str();

   if (handle) /* closing the .XLS file [Workbook] */
   {
      ret = freexl_close(handle);
      if (ret != FREEXL_OK)
         fprintf (stderr, "CLOSE ERROR: %d\n", ret);
   }
   if (!error)
   {
      if (!line_count)
         cout << "The file: " << fname << " does not have any rows with data." << endl;
      return line_count;
   }
   else
      return 0;
}

/*  Scan the input file list and create an array of class instances for each file with 
    interesting info.  The file list is an .xls file.  The info is put into an
    array of strings
    Also add to the AllBins container an FBin for each file.
*/
bool GetFileInfo(VFiles &file_info, AllBins& bins)
{
   string fname;
   string line;
   string tmpnum;
   string ans;
   stringstream file_list;

   if (!DoFilesExist(InName))
   {
      cout << "There are missing files." << endl;
      cout << "Do you want to continue anyway (Y/N)? ";
      getline(cin,ans);
      if (ans.length() > 0 && toupper(ans[0])=='Y')
      {
      }
      else
         return false;
   }
   else
      cout << "Found all files, looking good. . ." << endl;

   if (!GetRowsFromSprd(InName, FILELIST, file_list))
      return false;

   while (!file_list.eof())
   {
      getline(file_list,line);
      if (line.empty())  // ignore blank lines, #comments, lines with just whitespace
         continue;
      size_t pos = line.find_first_not_of("' '\t");
      if (pos == string::npos || line[pos] == '#')
         continue;

      FileFacts facts;
      istringstream fields(line);
      getline(fields,facts.path,'\t');
      getline(fields,facts.fname,'\t');
      getline(fields,tmpnum,'\t');

      try { facts.pmark = stoi(tmpnum); }
      catch(...) { }
      if(facts.pmark == -1)
      {
         cout << "Skipping line (and file): [" << line << "]" << endl;
         cout << "because period marker value is not a number" << endl;
         continue;
      }

      getline(fields,tmpnum,'\t');
      try {facts.ctl_period_len = stoi(tmpnum);}
      catch(...) { }

      getline(fields,tmpnum,'\t');
      try {facts.num_cells =stoi(tmpnum);}
      catch(...){ }

      getline(fields,facts.info_path,'\t');
      getline(fields,facts.info_file,'\t');
      getline(fields,facts.changes_path,'\t');
      getline(fields,facts.changes_file,'\t');

      getline(fields,tmpnum,'\t');
      try {facts.cco2stim_code =stoi(tmpnum);}
      catch(...){ }
      getline(fields,tmpnum,'\t');
      try {facts.cco2stim_interval =stoi(tmpnum);}
      catch(...){ }
      getline(fields,tmpnum,'\t');
      try {facts.vco2stim_code =stoi(tmpnum);}
      catch(...){ }
      getline(fields,tmpnum,'\t');
      try {facts.vco2stim_interval =stoi(tmpnum);}
      catch(...){ }

      getline(fields,tmpnum,'\t');
      try {facts.tb_cgh_stim_code =stoi(tmpnum);}
      catch(...){ }
      getline(fields,tmpnum,'\t');
      try {facts.tb_cgh_ctl_interval =stoi(tmpnum);}
      catch(...){ }
      getline(fields,tmpnum,'\t');
      try {facts.tb_cgh_stim_cycles =stoi(tmpnum);}
      catch(...){ }
      getline(fields,tmpnum,'\t');
      try {facts.lar_cgh_stim_code =stoi(tmpnum);}
      catch(...){ }
      getline(fields,tmpnum,'\t');
      try {facts.lar_cgh_ctl_interval =stoi(tmpnum);}
      catch(...){ }
      getline(fields,tmpnum,'\t');
      try {facts.lar_cgh_stim_cycles =stoi(tmpnum);}
      catch(...){ }
      getline(fields,facts.vagotomized,'\t');

      getline(fields,tmpnum,'\t');
      try {facts.swall1_start_code =stoi(tmpnum);}
      catch(...){ }
      getline(fields,tmpnum,'\t');
      try {facts.swall1_end_code =stoi(tmpnum);}
      catch(...){ }

      getline(fields,tmpnum,'\t');
      try {facts.lareflex_start_code =stoi(tmpnum);}
      catch(...){ }
      getline(fields,tmpnum,'\t');
      try {facts.lareflex_end_code =stoi(tmpnum);}
      catch(...){ }

      if (facts.cco2stim_code > 0 || facts.vco2stim_code >0 || facts.tb_cgh_stim_code > 0 || facts.lar_cgh_stim_code > 0 || facts.swall1_start_code > 0 || facts.lareflex_start_code > 0)
         facts.have_markers = true;

      getline(fields,facts.notes,'\n');  // gobble up rest of line

      if (vago_state == 0 && (facts.vagotomized.compare("Y") == 0 || facts.vagotomized.compare("y") == 0))
      {
         cout << "Skipping vagotomized experiment " << facts.fname << endl;
         continue;
      }
      else if (vago_state == 1 && (facts.vagotomized.compare("N") == 0 || facts.vagotomized.compare("N") == 0))
      {
         cout << "Skipping non-vagotomized experiment " << facts.fname << endl;
         continue;
      }

      facts.cmplt_name = facts.path+"/"+facts.fname;
      ifstream file(facts.cmplt_name);  // what kind of file?
      if (file)
      {
         getline(file,line);
         if (line.find("   33   3333333") != string::npos)
            facts.ftype =  F_TYPE::EDT;
         else if (line.find("   11 1111111") != string::npos)
            facts.ftype = F_TYPE::BDT;
         else
            facts.ftype = F_TYPE::ADT; // any way to know for sure?

         file_info.push_back(facts);   // remember the file facts

            // by the way, this ensures there is only one instance of a file
            // even if the name is in the file list more than once
         bins.insert(make_pair(facts.fname,FBins()));
      }
      else
      {
        cout << "Could not open [" << facts.cmplt_name << "]" << ", skipping. . ." <<  endl;
        continue;
      }
   }
   if (file_info.size() == 0)
   {
      cout << "Unable to open any files." << endl << "Check to see if you have access to file system where the files are kept." << endl << "Also make sure the input spreadsheet file has valid information, aborting program. . ." << endl;
      return false;
   }
   else
   {
      cout << "Using " << file_info.size() << " files." << endl; 
      return true;
   }
}

/* Peer into the stream and determine time values for:  i0----e0------i1
   Assumes:  The caller has detected the start of an i phase and the
             next line in the file stream is the i0 event.
   Returns:  if still in same period
               file times for i0, e0, i1 events and and file positioned at i0
             else (if detected end of period or eof)
                no data and
                file positioned at start of next period if period marker found
                or
                file positioned at eof if at eof.

   This also checks for and filters out loss of phase intervals, such as two
   I's in a row without an E, and also if the cycle is longer than a cutoff
   value (10 seconds at this time) or shorter than a cutoff value (1 second).
   For single tick bin counting, on the first pass through the period we 
   create some values we will use on the second pass to filter bins in or out.
*/
bool NeuroBins::FindPhaseMarks(ifstream& file,long marker, AllBinsIter& f_bins, PeriodInfo&p_data)
{
   string    line;
   streampos curr_pos;
   bool      ret = false;
   bool      found_e = false;
   bool      found_i = false;
   bool      repos = true;
   long      i_s =0;
   long      e_s =0;
   long      n_i =0;
   streampos lastline;

   curr_pos = file.tellg();     // remember i0
   while (!file.eof())
   {
      lastline = file.tellg();  // we may need to back up one line
      getline(file,line);
      if (line.empty())
         continue;
      DataLine  d_line(line);

      if (d_line.IsMarker(marker)) // partial last phase, nothing to return
      {
         file.seekg(lastline);     // caller will want to see marker, back up
         repos = false;
         break;
      }
      else if (d_line.IsIPhase() && found_i && !found_e)
      {
         cout << "Missing E phase" << endl;
         found_i = found_e = false;  // start hunting again
      }
      else if (d_line.IsEPhase() && found_e && found_i)
      {
         cout << "Missing I phase" << endl;
         found_i = found_e = false;
      }
      else if (d_line.IsIPhase() && !found_i)  // found i0 start
      {
         i_s = d_line.GetTime();
         found_i = true; 
         if (Debug) cout << line << endl << " i seg start " << i_s << endl;
      }
      else if (d_line.IsEPhase() && found_i && !found_e)  // find e0 start
      {
         e_s = d_line.GetTime();
         if (Debug) cout << line << endl << "e seg start " << e_s << endl;
         found_e = true; 
      }
      else if (d_line.IsIPhase() && found_e)  // found i1, end of cycle
      {
         n_i = d_line.GetTime();
         if (Debug) cout << line << endl << "next i seg is " << n_i  << endl;

          // now have everything we need to set up globals shared by all
          // spikes that occur in this cycle
         i_start = i_s;
         e_start = e_s;
         next_i  = n_i;
         i_time = (double) (e_s - i_s) / num_i_bins;
         e_time = (double) (n_i - e_s) / num_e_bins;

         ip_time = (e_s - i_s);   // one tick bin widths in ticks
         ep_time = (n_i - e_s);
if (Debug) cout << "FPM  i_start: " << i_start << 
                   " e_start: " << e_start << 
                   " next_i: "  << next_i << 
                   " i_time: "  << i_time << 
                   " e_time: "  << e_time << 
                   " ip_time: "  << ip_time << 
                   " ep_time: "  << ep_time <<  endl;
         if (FilePass1)    // build one tick hists on first pass
            UpdatePHist(CTL_PERIOD);  // always do this, even if we reject the cycle later

          // Reject this cycle?
         double tcheck=(next_i - i_start) * (NeuroBins::GetSecsTick()/1000.0);
         if ((tcheck > MAX_CYCLE_TIME || tcheck < MIN_CYCLE_TIME)
             || (!FilePass1 && !NeuroBins::UseCycle(p_data)))
         {
            if (Debug) cout << " reject is: " << i_s << "  reject es " << e_s << " ip_time "<< ip_time << " ep_time " << ep_time << endl; 
            file.seekg(lastline);     // backup to stay on I1 marker as new I0 
            curr_pos = file.tellg();  // remember new place
            found_i = found_e = false;
            i_start = e_start = next_i = 0;
            PurgeAllSaved(f_bins,p_data);  // toss any pending carry-overs
         }
         else
         { 
            ret = true;
            break;
         }
      }
   }
   if (file.eof()) // leave it here
      repos = false;
   if (repos)
      file.seekg(curr_pos);  // seek back to i0

   return ret;
}

/* Handle moving through a ctl/stim section of a recording.
   States:  1. Next line is marker and offset < 0
               Start of new ctl/stim section.
               Rewind back to offset.  
               Seek forward to find first i0.
               Return with file position on i0.
            2. Marker on next i0 and offset < 0
               We're working through i0 e0 i1 for ctl part of ctl/stim
               Look for next i0 e0 i1 and return line on i0
               If hit marker before finding i1, dump partial results
               and leave next line on marker
               If hit eof, dump partial results and leave at eof.
            3. Marker on ctl/stim OR marker on i0 AND offset > 0 
               Working interval stim part of ctl/stim  
               Look for next i0 e0 i1 and return line on i0
               If hit end of interval before finding i1, dump partial results
               and leave next line where it is on entry.
               If hit eof, dump partial results and leave at eof.
            4. Marker on ctl/stim OR marker on i0 AND  offset == 0 
               if cycle > 0
                  find next i0 e0 i1 and return it.
               else
                  if cycle == 0,  done with this marker, clean up and
                  leave line where we found it.
   Returns:  file times for i0, e0, i1 events and and file positioned at i0
             else if detected partial last i0e0i1 or eof, return no data 
             and file positioned at next location (marker or past stim interval/cycle)
*/
bool NeuroBins::FindCtlStimMarks(ifstream& file, int code, int interval,AllBinsIter& f_bins, PeriodInfo& p_data, int stim_cycles)
{
   string    line;
   streampos curr_pos, i0_pos;
   bool      ret = false;
   bool      found_e = false;
   bool      found_i = false;
   bool      done = false;
   long      i_s =0;
   long      e_s =0;
   long      n_i =0;
   long      backline;
   long      interval_time;
   unsigned long targ_time;
   streampos lastline, backup;
   DataLine  mark, d_line;

   interval_time = interval * ((1.0/NeuroBins::GetSecsTick() * 1000));
   curr_pos = file.tellg();     // remember current stream pos
   backup = curr_pos;
   getline(file,line);
   file.seekg(backup);
   d_line = line;

   if (interval < 0)
   {
      if (d_line.IsMarker(code)) // start of ctl period 
      {
         targ_time = d_line.GetTime() + interval_time;
         backline = line.size() + 1; // for newline
         while (!file.eof())
         {
            backup -= backline;
            file.seekg(backup);
            getline(file,line);
            d_line = line;
            if (d_line.GetTime() <= targ_time)
            {
               backup -= backline; // include this line
               file.seekg(backup);
               break;
            }
            file.seekg(backup);
         }

         if (file.eof())
         {
            PurgeAllSaved(f_bins,p_data);  // toss any pending carry-overs
            i_s = e_s = n_i = 0;
            return false;
         }
      }

      while (!done && !file.eof())
      {
         lastline = file.tellg();
         getline(file,line);
         if (line.empty())
            continue;
         d_line = line;
         if (d_line.IsMarker(code)) // found end of ctl period
         {
              // todo: need to select just cco2ctl, vco2ctl, no ctl
            PurgeAllSaved(f_bins,p_data);  // toss any pending carry-overs
            i_s = e_s = n_i = 0;
            i0_pos = lastline;    // leave on marker
            ret = false;
            done = true;
         }
         else if (d_line.IsIPhase() && found_i && !found_e)
         {
            cout << "Missing E phase" << endl;
            found_i = found_e = false;  // start hunting again
         }
         else if (d_line.IsEPhase() && found_e && found_i)
         {
            cout << "Missing I phase" << endl;
            found_i = found_e = false;
         }
         else if (d_line.IsIPhase() && !found_i)  // found i0 start
         {
            i_s = d_line.GetTime();
            i0_pos = lastline;     // remember start of i0
            found_i = true; 
            if (Debug) cout << line << endl << " i seg start " << i_s << endl;
         }
         else if (d_line.IsEPhase() && found_i && !found_e)  // find e0 start
         {
            e_s = d_line.GetTime();
            if (Debug) cout << line << endl << "e seg start " << e_s << endl;
            found_e = true; 
         }
         else if (d_line.IsIPhase() && found_e)  // found i1, end of cycle
         {
            n_i = d_line.GetTime();
            if (Debug) cout << line << endl << "next i seg is " << n_i  << endl;
             if (FilePass1)    // build one tick hists on first pass
                UpdatePHist(p_data.period);
               // In FindPhaseMarks, which this is derived from, we
               // reject cycles that are too long or too short.  We keep
               // everything in the ctl/stim periods
            ret = true;
            done = true;
         }
      }
   }
   else if (interval > 0)   // doing stim by time/interval
   {
      if (d_line.IsMarker(code)) // start of stim period, how long until end? 
         stim_end = d_line.GetTime() + interval_time;

        // now find i0
      while (!done)
      {
         lastline = file.tellg();
         getline(file,line);
         if (line.empty())
            continue;
         d_line = line;
         if (d_line.GetTime() > stim_end)
         {
              // todo: need to select just cco2ctl, vco2ctl, no ctl
            PurgeAllSaved(f_bins,p_data);  // toss any pending carry-overs
            i_s = e_s = n_i = 0;
            ret = false;
            done = true;
         }
         else if (d_line.IsIPhase() && found_i && !found_e)
         {
            cout << "Missing E phase" << endl;
            found_i = found_e = false;  // start hunting again
         }
         else if (d_line.IsEPhase() && found_e && found_i)
         {
            cout << "Missing I phase" << endl;
            found_i = found_e = false;
         }
         else if (d_line.IsIPhase() && !found_i)  // found i0 start
         {
            i_s = d_line.GetTime();
            i0_pos = lastline;     // remember start of i0
            found_i = true; 
            if (Debug) cout << line << endl << " i seg start " << i_s << endl;
         }
         else if (d_line.IsEPhase() && found_i && !found_e)  // find e0 start
         {
            e_s = d_line.GetTime();
            if (Debug) cout << line << endl << "e seg start " << e_s << endl;
            found_e = true; 
         }
         else if (d_line.IsIPhase() && found_e)  // found i1, end of cycle
         {
            n_i = d_line.GetTime();
            if (Debug) cout << line << endl << "next i seg is " << n_i  << endl;

            if (FilePass1)                  // build one tick hists on first pass
               UpdatePHist(p_data.period);
            ret = true;
            done = true;
         }
      }
   }
   else if (interval == 0)   // doing stim by cycles (just cough at this time)
   {
      if (stim_cycles == 0)
      {
            PurgeAllSaved(f_bins,p_data);  // toss any pending carry-overs
            i_s = e_s = n_i = 0;
            ret = false;
            done = true;
            i0_pos = curr_pos;   // leave stream where we found it
      }
      else
      {
         while (!done) // find i0
         {
            lastline = file.tellg();
            getline(file,line);
            if (line.empty())
               continue;
            d_line = line;
            if (d_line.IsIPhase() && found_i && !found_e)
            {
               cout << "Missing E phase" << endl;
               found_i = found_e = false;  // start hunting again
            }
            else if (d_line.IsEPhase() && found_e && found_i)
            {
               cout << "Missing I phase" << endl;
               found_i = found_e = false;
            }
            else if (d_line.IsIPhase() && !found_i)  // found i0 start
            {
               i_s = d_line.GetTime();
               i0_pos = lastline;     // remember start of i0
               found_i = true; 
               if (Debug) cout << line << endl << " i seg start " << i_s << endl;
            }
            else if (d_line.IsEPhase() && found_i && !found_e)  // find e0 start
            {
               e_s = d_line.GetTime();
               if (Debug) cout << line << endl << "e seg start " << e_s << endl;
               found_e = true; 
            }
            else if (d_line.IsIPhase() && found_e)  // found i1, end of cycle
            {
               n_i = d_line.GetTime();
               if (Debug) cout << line << endl << "next i seg is " << n_i  << endl;
               if (FilePass1)                  // build one tick hists on first pass
                  UpdatePHist(p_data.period);
               ret = true;
               done = true;
            }
         }
      }
   }

     // now have everything we need to set up globals shared by all
     // spikes that occur in this cycle
   i_start = i_s;
   e_start = e_s;
   next_i  = n_i;
   i_time = (double) (e_s - i_s) / num_i_bins;
   e_time = (double) (n_i - e_s) / num_e_bins;
   ip_time = (e_s - i_s);   // one tick bin widths in ticks
   ep_time = (n_i - e_s);

if (Debug) cout << "FCTLPM  i_start: " << i_start << 
                   "e_start: " << e_start << 
                   "next_i: "  << next_i << 
                   "i_time: "  << i_time << 
                   "e_time: "  << e_time << 
                   "ip_time: "  << ip_time << 
                   "ep_time: "  << ep_time <<  endl;

   if (!file.eof())        // leave it at eof 
      file.seekg(i0_pos);  // else seek back to i0 or marker
   return ret;
}


/* Find the end of the current swallow1 or lareflex period. This is unlike the
   other periods in that they are not bi-phasic.  There are no phases in
   swallow1 or lareflex, just a start and end marker. I and E does not figure
   into this.  So, binning is simpler. There is no phase normalization, so we
   only need to know how long the interval between the start and end are. We
   also do not have to worry about carrying the last I spike over to the first
   E bin.  Assumes the next line in the file is the appropriate start marker.
*/
bool NeuroBins::FindSwall1Marks(ifstream& file, int start, int end, PeriodInfo& p_data)
{
   string    line;
   streampos curr_pos;
   bool      ret = false;
   bool      found_start = false;
   bool      repos = true;
   long      s_s =0;
   long      s_e =0;
   streampos lastline;

   curr_pos = file.tellg();     // remember current start position
   while (!file.eof())
   {
      lastline = file.tellg();  // we may need to back up one line
      getline(file,line);
      if (line.empty())
         continue;
      DataLine  d_line(line);
      if (d_line.IsMarker(start))
      {
         if (!found_start)
         {
            found_start = true;
            i_start = s_s = d_line.GetTime();      // start time
            continue;
         }
         else  // we seem to have missed the end mark, skip this one 
         {
            s_s = d_line.GetTime();  // start over 
            curr_pos = file.tellg();
            file.seekg(lastline);    // seek back to this start code
            found_start = false;     // and look some more
         }
      }
      else if (d_line.IsMarker(end))
      {
         if (found_start) 
         {
             // now have everything we need to set up globals shared by all
             // spikes that occur in this period
             // use same values for both i & e times even though there is no i or e
             // for swallow or lareflex
            s_e = d_line.GetTime();
            i_time = e_time = (double) (s_e - s_s) / (num_i_bins + num_e_bins);
            e_start = i_time/2;  // split "period" into 2 equal intervals
            ip_time = ep_time = (s_e - s_s)/2;
   if (Debug) cout << "FPM  i_start: " << i_start << 
                      " e_start: " << e_start << 
                      " next_i: "  << next_i << 
                      " i_time: "  << i_time << 
                      " e_time: "  << e_time << 
                      " ip_time: "  << ip_time << 
                      " ep_time: "  << ep_time <<  endl;
            if (FilePass1)           // build one tick hists on first pass
            {
               UpdatePHist(p_data.period);
            }
            ret = true;
            break;
         }
         else  // end without a start, skip ahead
         { 
            continue;
         }
      }
      else if (d_line.IsMarker(start))
      {
         if (!found_start)
         {
            found_start = true;
            i_start = s_s = d_line.GetTime();      // start time
            continue;
         }
         else  // we seem to have missed the end mark, skip this one 
         {
            s_s = d_line.GetTime();  // start over 
            curr_pos = file.tellg();
            file.seekg(lastline);    // seek back to this start code
         }
      }
   }
   if (file.eof()) // leave it here
      repos = false;
   if (repos)
      file.seekg(curr_pos);  // seek back to starting position

   return ret;
}


// scan the current file from the current location until EOF and count how many
// ctl/stim markers there are.
// Return the count in the passed in f_iter struct and the file position where
// it was when we were called.
// NOTE: Adding new types of markers will required adjustments to this.
void GetNumCtlStimMarkers(ifstream &file,VFilesIter& f_iter)
{
   string        line;
   DataLine      d_line;
   streampos     curr_pos = file.tellg();  // will seek back to here
   bool start = false;
   while (!file.eof())
   {
      getline(file,line);
      if (line.empty())
         continue;
       d_line = line;
       if (d_line.IsMarker(f_iter->cco2stim_code))
          ++f_iter->cco2marker_count;
       else if (d_line.IsMarker(f_iter->vco2stim_code))
          ++f_iter->vco2marker_count;
       else if (d_line.IsMarker(f_iter->tb_cgh_stim_code))
          ++f_iter->tb_cghmarker_count;
       else if (d_line.IsMarker(f_iter->lar_cgh_stim_code))
          ++f_iter->lar_cghmarker_count;
       else if (d_line.IsMarker(f_iter->swall1_start_code))
       {
          ++f_iter->swall1marker_s_count;
          if (!start)
             start = true;
          else
          {
             cout << "missing swallow1 end marker at " << d_line.GetTime() << endl;
             start = false;
          }
       }
       else if (d_line.IsMarker(f_iter->swall1_end_code))
       {
          ++f_iter->swall1marker_e_count;
          if (start)
             start = false;
          else
          {
             cout << "missing swallow1 start marker at " << d_line.GetTime() << endl;
             start = true;
          }
       }
       else if (d_line.IsMarker(f_iter->lareflex_start_code))
       {
          ++f_iter->lareflexmarker_s_count;
          if (!start)
             start = true;
          else
          {
             cout << "missing lareflex end marker at " << d_line.GetTime() << endl;
             start = false;
          }
       }
       else if (d_line.IsMarker(f_iter->lareflex_end_code))
       {
          ++f_iter->lareflexmarker_e_count;
          if (start)
             start = false;
          else
          {
             cout << "missing lareflex start marker at " << d_line.GetTime() << endl;
             start = true;
          }
       }
   }

   file.clear();
   file.seekg(curr_pos);
   cout << "found " << f_iter->cco2marker_count << " CCO2 markers" << endl; 
   cout << "found " << f_iter->vco2marker_count << " VCO2 markers" << endl; 
   cout << "found " << f_iter->tb_cghmarker_count << " TB COUGH markers" << endl; 
   cout << "found " << f_iter->lar_cghmarker_count << " LAR COUGH markers" << endl; 
   cout << "found " << f_iter->swall1marker_s_count << " SWALLOW1 start markers" << endl; 
   cout << "found " << f_iter->swall1marker_e_count << " SWALLOW1 end markers" << endl; 
   cout << "found " << f_iter->lareflexmarker_s_count << " LAREFLEX start markers" << endl; 
   cout << "found " << f_iter->lareflexmarker_e_count << " LAREFLEX end markers" << endl; 
   if (f_iter->swall1marker_s_count || f_iter->swall1marker_e_count)
      if (f_iter->swall1marker_s_count != f_iter->swall1marker_e_count)
         cout << "ERROR:  There the number of swallow1 start and end markers are not the same!" << endl;
   if (f_iter->lareflexmarker_s_count || f_iter->lareflexmarker_e_count)
      if (f_iter->lareflexmarker_s_count != f_iter->lareflexmarker_e_count)
         cout << "ERROR:  There the number of lareflex start and end markers are not the same!" << endl;
}

// Scan through the list of CTHs to save and flag the ones that are sparse.
// Note that the threshold count for this changes depending on the # of bins.
// We still want to write out all of the CTHs so we can reuse the curve fitting
// output.
// At this time, we do not mark any of the ctl/stim CTHs as sparse, even ones
// that would be considered sparse.  It is thought that this is throwing away
// potentially useful information.
void MarkSparses(AllBins& bins)
{
   AllBinsIter   f_bins;
   FBinsIter     fbi;
   unsigned long total_spikes, rejects = 0;

   cout << "Filtering sparse. . ." << endl;
   for (f_bins = bins.begin(); f_bins != bins.end(); ++f_bins)
   {
      for (fbi = f_bins->second.begin(); fbi != f_bins->second.end(); ++fbi)
      {
         if (fbi->first.second == CtlStims[CTL_PERIOD]) // mark only control period CTHs
         {
            total_spikes = fbi->second.GetTotalSpikes();
            if (total_spikes < RejectThresh)
            {
               if (Debug)
               {
                  cout << "pair: " << fbi->first.first << "  "  << fbi->first.second << endl;
                  cout << "Total spikes: " << total_spikes << endl;
                  cout <<  "BinMeans: " << fbi->second.GetBinMeans() << endl;
                }
                fbi->second.SetSparse();
                ++rejects;
            }
         }
      }
   }
   if (rejects)
      cout << "Marked " << rejects << " as sparse Control Period CTHs because they contained less than " << RejectThresh << " spikes" << endl;
   else
      cout << "No sparse CTHs" << endl;

   NumSparse = rejects;  // save in output file later
}

// For swallow control period, we want to include the sparse cths. Unmark
// the ones we previously marked. This should be one of the last things in
// the results saving control flow. These only go into the ctl_stim pairs cth file.
void UnMarkSparses(AllBins& bins)
{
   AllBinsIter   f_bins;
   FBinsIter     fbi;
   for (f_bins = bins.begin(); f_bins != bins.end(); ++f_bins)
      for (fbi = f_bins->second.begin(); fbi != f_bins->second.end(); ++fbi)
         if (fbi->first.second == CtlStims[CTL_PERIOD])
                fbi->second.UnSetSparse();
   NumSparse = 0;
}

// Since the only way to create a NeuroBins record is if there are spikes for
// the cell, AND if we are doing ctl/stim periods, we want there to be a pair
// of CTHs, even if one of them has no spikes in it.  Find ctl/stim CTH and if
// it does not have a sibling, create a blank one for it.
void CreateBlanks(AllBins& bins)
{
   AllBinsIter   f_bins;
   FBinsIter     fbi, found_fbi;
   string        fname;
   int           chan, period;

   for (f_bins = bins.begin(); f_bins != bins.end(); ++f_bins)
   {
      for (fbi = f_bins->second.begin(); fbi != f_bins->second.end(); ++fbi)
      {
         period = fbi->first.second;
         {
             // find companion ctl/stim for this cell
            if (period == CtlStims[CTL_CCO2])
               period = CtlStims[STIM_CCO2];
            else if (period == CtlStims[STIM_CCO2])
               period = CtlStims[CTL_CCO2];
            else if (period == CtlStims[CTL_VCO2])
               period = CtlStims[STIM_VCO2];
            else if (period == CtlStims[STIM_VCO2])
               period = CtlStims[CTL_VCO2];
            else if (period == CtlStims[TB_CTL_CGH])
               period = CtlStims[TB_STIM_CGH];
            else if (period == CtlStims[TB_STIM_CGH])
               period = CtlStims[TB_CTL_CGH];
            else if (period == CtlStims[LAR_CTL_CGH])
               period = CtlStims[LAR_STIM_CGH];
            else if (period == CtlStims[LAR_STIM_CGH])
               period = CtlStims[LAR_CTL_CGH];
            else if (have_swall1 && period == CtlStims[CTL_PERIOD])
              // For swall1 and lareflex periods, the control period is the sibling.
              // There may not be a ctl cth in this case, so make one
               period = CtlStims[STIM_SWALL1];
            else if (period == CtlStims[STIM_SWALL1]) // there is no ctl_swall1 period,
               period = CtlStims[CTL_PERIOD]; // we re-use the control period for its sib
            else if (have_lareflex && period == CtlStims[CTL_PERIOD])
               period = CtlStims[STIM_LAREFLEX];
            else if (period == CtlStims[STIM_LAREFLEX])
               period = CtlStims[CTL_PERIOD];

            chan = fbi->second.GetChan();
            FBinKey key(chan,period);
            found_fbi = f_bins->second.find(key);
            if (found_fbi == f_bins->second.end())  // make the companion, a "blank" CTH
            {
               cout << "Create blank CTH for exp: " << f_bins->first << " chan: " << chan << " period: " << period << endl;
               found_fbi = f_bins->second.insert(        // create a record
                 FBins::value_type(key,
                 NeuroBins(
                 chan,
                 fbi->second.GetCoords(),
                 period,
                 0,
                 0,
                 NeuroBins::GetSecsTick()))).first;
             }
         }
      }
   }
}

// Remove all the blank CTHs for the given period
// When we save just the swall1 or lareflex CTHs, we do not want the
// the blank CTHs included, because there are no pairs in these files.
// This should be done after all other periods have been saved to files.
void RemoveBlanks(AllBins& bins, int period)
{
   AllBinsIter   f_bins;
   FBinsIter     fbi, found_fbi;
   string        fname;
   int           chan, curr_period, total_spikes;

   for (f_bins = bins.begin(); f_bins != bins.end(); ++f_bins)
   {
      for (fbi = f_bins->second.begin(); fbi != f_bins->second.end(); )
      {
         curr_period = fbi->first.second;
         if (curr_period == period)
         {
            chan = fbi->second.GetChan();
            FBinKey key(chan,period);
            found_fbi = f_bins->second.find(key);
            if (found_fbi != f_bins->second.end())
            {
               total_spikes = fbi->second.GetTotalSpikes();
               if (total_spikes == 0)
               {
                  fbi = f_bins->second.erase(found_fbi);
                  cout << "Remove blank CTH for exp: " << f_bins->first << " chan: " << chan << " period: " << period << endl;
                  continue;
               }
            }
         }
         ++fbi;
      }
   }
}




// Save the results in a format ocatave can read, and perhaps other formats.
void SaveResults(AllBins& bins, P_LIST& periods, PerPeriod &p_data)
{
   AllBinsIter   f_bins;
   FBinsIter     fbi, tmpfbi;
   unsigned long total_spikes;
   string        varname;
   string        expname;
   string        fname, r_inname,r_outname;
   string        split;
   unsigned int  pt_num;
   int           exp_num;
   int           curr_period;
   int           save_period;
   size_t        last_idx;
   pair <long int,long int> p;
   CTHCurves     Rvals;
   double        max_cth = 0;
   bool          have_cth;
   bool          ctl_swall1_period=false;
   bool          ctl_lareflex_period=false;

   if (bins.size() == 0)
   {
      cout << "There were no files read, so there are no results to save." << endl << "Check to see if you have access to file system where the files are kept." << endl << "Also make sure the input spreadsheet file has valid information." << endl;
      return;
   }

   if (AllP)
      split = PNames[ALL_PER];
   else
   {
      switch (periods[0])   // if periods[1], it will be same type of ctl/stim
      {
         case CTL_PERIOD:
            split = PNames[CTL_PER];
            break;
         case CTL_CCO2:
            split = PNames[CCO2_PER];
            break;
         case CTL_VCO2:
            split = PNames[VCO2_PER];
            break;
         case TB_CTL_CGH:
            split = PNames[TB_CGH_PER];
            break;
         case LAR_CTL_CGH:
            split = PNames[LAR_CGH_PER];
            break;
         case STIM_SWALL1:
            split = PNames[SWALL1_PER];
            break;
         case CTL_SWALL1:
            split = PNames[CTL_SWALL1_PER];
            ctl_swall1_period=true;
            break;
         case STIM_LAREFLEX:
            split = PNames[STIM_LAREFLEX_PER];
            break;
         case CTL_LAREFLEX:
            split = PNames[CTL_LAREFLEX_PER];
            ctl_lareflex_period=true;
            break;
         default:
            cout << "Unhandled type of period name in SaveResults" << endl;
            split = "unkn";
            break;
      }
   }

   if (OutName.size() == 0)
   {
      fname="cth_" + to_string(NumEBins+NumIBins) + split + ".cth";
      r_inname = "cth_" + to_string(NumEBins+NumIBins) + "_to_r.txt"; 
      r_outname = "cth_" + to_string(NumEBins+NumIBins) + "_from_r.txt";
   }
   else
   {
      fname = OutName + split + ".cth";
      r_inname  = OutName + "_to_r.txt";
      if (RdFile.length() > 0) // use over-ride filename?
         r_outname = RdFile;
      else
         r_outname = OutName + "_from_r.txt";
   }
   cout << "Saving results  to " << fname << endl;

     // always create file for octave
   ofstream oct_vars(fname.c_str());
   oct_vars << setprecision(15);
   pt_num = 1;
   exp_num = 1;

   NeuroBins::WriteOVersion(oct_vars);
   int contp,cco2,vco2,tbcgh,larcgh,ctl_swall1,stim_swall1,ctl_lareflex,stim_lareflex;
   contp=cco2=vco2=tbcgh=larcgh=ctl_swall1=stim_swall1=ctl_lareflex=stim_lareflex=0;

   GetRVals(bins,r_inname,r_outname, Rvals); // get R curve plotting points
                                             // or dummy date if not available
// cout << "bins: " << bins.size() << endl;
   for (f_bins = bins.begin(); f_bins != bins.end(); ++f_bins) // for all files
   {
// cout << "f_bins: " << f_bins->second.size() << endl;
      have_cth = false;
      for (fbi = f_bins->second.begin(); fbi != f_bins->second.end(); ++fbi)
      {                                              // all bins in file
         save_period = curr_period = fbi->first.second;
         if (curr_period > MAX_PERIODS-1)
            cout << "Unknown period: " << curr_period << endl;

         if (find(periods.begin(), periods.end(), curr_period) != periods.end())
         {                                          // all bins in period(s)
            if (curr_period==CTL_PERIOD)
               ++contp;
            else if (curr_period==CTL_CCO2 || curr_period==STIM_CCO2)
               ++cco2;
            else if (curr_period==CTL_VCO2 || curr_period==STIM_VCO2)
               ++vco2;
            else if (curr_period==TB_CTL_CGH || curr_period==TB_STIM_CGH)
               ++tbcgh;
            else if (curr_period==LAR_CTL_CGH || curr_period==LAR_STIM_CGH)
               ++larcgh;
            else if (curr_period==CTL_SWALL1)
               ++ctl_swall1;
            else if (curr_period==STIM_SWALL1)
               ++stim_swall1;
            else if (curr_period==CTL_LAREFLEX)
               ++ctl_lareflex;
            else if (curr_period==STIM_LAREFLEX)
               ++stim_lareflex;
            else
               cout << "Unknown period: " << curr_period << endl;

            if (ctl_swall1_period)
               if (curr_period == CTL_PERIOD) // use a different period# in the save
                  save_period = CTL_SWALL1;
            if (ctl_lareflex_period)
               if (curr_period == CTL_PERIOD)
                  save_period = CTL_LAREFLEX;

                 // create var names that octave likes
            have_cth = true;
            varname = f_bins->first;
            expname = varname;
            last_idx = varname.find_last_of("\\/");  // strip path
            varname.erase(0,last_idx+1);
            replace(varname.begin(), varname.end(),'.','_');
            replace(varname.begin(), varname.end(),'-','_');
            last_idx = expname.find_last_of("\\/");  // strip path
            expname.erase(0,last_idx+1);

            total_spikes = fbi->second.GetTotalSpikes();
            fbi->second.WriteOctaveHeader(oct_vars,expname,varname,pt_num,total_spikes,save_period);
            fbi->second.WriteOIBins(oct_vars);
            fbi->second.WriteOEBins(oct_vars);
            fbi->second.WriteOStats(oct_vars,pt_num);
            fbi->second.WriteOCoords(oct_vars);
            fbi->second.WriteDistStats(oct_vars);
            fbi->second.WriteOneTicks(oct_vars);
               // todo the next is probably wrong, but we don't do curve fitting
               // very much these days, so this may not matter
            string key = MakeVarName(f_bins->first,fbi->second.GetChan(),save_period);
            fbi->second.WriteRVals(oct_vars,Rvals,key);
            NeuroBins::WriteOEndMark(oct_vars);
            ++pt_num;

            tmpfbi = fbi;
            if (!max_cth)
               max_cth = maxCTH[CtlStims[curr_period]];
         }
         else
         {
            if (Debug) cout << "Skipping period " << curr_period << " fbi " << endl;
         }
      }
      if (have_cth)
      {
         tmpfbi->second.WriteOExpName(oct_vars, varname, exp_num);
         tmpfbi->second.WriteOSparse(oct_vars);
      }
      ++exp_num;
   }
   WriteDist(oct_vars,bins,periods,p_data); 
   WriteZeroFlat(oct_vars,pt_num,max_cth);
   oct_vars.close();

   cout << "found " << contp << " control cths" << endl;
   cout << "found " << cco2 << " cco2 cths" << endl;
   cout << "found " << vco2 << " vco2 cths" << endl;
   cout << "found " << tbcgh << " tbcgh cths" << endl;
   cout << "found " << larcgh << " larcgh cths" << endl;
   cout << "found " << ctl_swall1 << " ctl_swall1 cths" << endl;
   cout << "found " << stim_swall1 << " stim_swall1 cths" << endl;
   cout << "found " << ctl_lareflex << " ctl_lareflex cths" << endl;
   cout << "found " << stim_lareflex << " stim_lareflex cths" << endl;
   cout << "found 1 ZeroFlat" << endl;
   cout << "Wrote " << pt_num << " points" << endl;

   if (Csv) // generic csv?
   {
      last_idx = fname.find_last_of(".");
      fname.erase(last_idx,string::npos);
      fname += ".csv";
      cout << "csv:  " << fname << endl;
      ofstream csv_points(fname.c_str());
      csv_points << setprecision(15);

      for (f_bins = bins.begin(); f_bins != bins.end(); ++f_bins)
      {
         for (fbi = f_bins->second.begin(); fbi != f_bins->second.end(); ++fbi)
         {
            curr_period = fbi->first.second;
            if (find(periods.begin(), periods.end(), curr_period) != periods.end())
            {
               fbi->second.WriteCIBins(csv_points);
               fbi->second.WriteCEBins(csv_points);
            }
         }
      }
      csv_points.close();
   }

       // ggobi csv?
   if (Gcsv)
   {
      bool have_header = false; 
      pt_num = 1;
      last_idx = fname.find_last_of(".");
      fname.erase(last_idx,string::npos);
      fname += "_g.csv";
      cout << "gcsv: " << fname << endl;
      ofstream gcsv_points(fname.c_str());
      gcsv_points << setprecision(15);

      for (f_bins = bins.begin(); f_bins != bins.end(); ++f_bins)
      {
         for (fbi = f_bins->second.begin(); fbi != f_bins->second.end(); ++fbi)
         {
            curr_period = fbi->first.second;
            if (find(periods.begin(), periods.end(), curr_period) != periods.end())
            {
               if (!have_header)  // only need one
               {
                  fbi->second.WriteGCHeader(gcsv_points);
                  have_header = true;
               }
               fbi->second.WriteGCIBins(gcsv_points,pt_num);
               fbi->second.WriteGCEBins(gcsv_points);
               ++pt_num;
            }
         }
      }
      gcsv_points.close();
   }
}


// We have all the data, caculate and save global scale and max cth values.
// Only p_data[0] holds results, other not needed
// else
// Iterate the periods in the p_data array, such that 
// use p_data[0] sets scale,max_cth for the ctl period,
// use p_data[1] and p_data[2] for their scale,max_cth  (CTL)
// use p_data[3] and p_data[4] for their scale,max_cth  (CCO2)
// use p_data[5] and p_data[6] for their scale,max_cth  (CGH)
// use p_data[7] and p_data[8] for their scale,max_cth  (SWALL1)
// use p_data[9] and p_data[10] for their scale,max_cth (LAREFLEX)
// Update the p_data array.  [1] and [2] will be same output
// Update the p_data array.  [3] and [4] will be same output
// Update the p_data array.  [5] and [6] will be same output
// Update the p_data array.  [7] and [8] will be same output
// Update the p_data array.  [9] and [10] will be same output
void CalcsForAllFiles(AllBins& bins, PerPeriod& p_data)
{
   AllBinsIter f_bins;
   FBinsIter   fbi;
   double      curr = 0;
   double      scale;
   int         curr_period;
   int         pair_max;

   MarkSparses(bins);  // Flag sparse CTHs
   CreateBlanks(bins); // Create a CTH for a CTH pair with no spikes
   curr_period = CtlStims[CTL_PERIOD]; // assume just period 0

       // scan bins to find largest bin value.
   if (Norm == MEAN)    // AKA Area
   {
      for (f_bins = bins.begin(); f_bins != bins.end(); ++f_bins)
      {
         for (fbi = f_bins->second.begin(); fbi != f_bins->second.end(); ++fbi)
         {
            curr_period = fbi->first.second;
            if (!fbi->second.IsSparse())
            {
               curr = fbi->second.GetBinMeans();
               if (curr > maxCTH[CtlStims[curr_period]])
               {
                  maxCTH[CtlStims[curr_period]] = curr;
                  if (Debug) cout << "New max mean value for all data in period " << curr_period << " is " << maxCTH[CtlStims[curr_period]] << endl;
               }
            }
         }
      }
        // need to normalize the ctl/stim cths to the same value, pick the
        // max for both sets of CTHs
     pair_max = max(maxCTH[CTL_CCO2],maxCTH[STIM_CCO2]);
     maxCTH[CTL_CCO2] = maxCTH[STIM_CCO2] = pair_max;
     pair_max = max(maxCTH[CTL_VCO2],maxCTH[STIM_VCO2]);
     maxCTH[CTL_VCO2] = maxCTH[STIM_VCO2] = pair_max;
     pair_max = max(maxCTH[TB_CTL_CGH],maxCTH[TB_STIM_CGH]);
     maxCTH[TB_CTL_CGH] = maxCTH[TB_STIM_CGH] = pair_max;
     pair_max = max(maxCTH[LAR_CTL_CGH],maxCTH[LAR_STIM_CGH]);
     maxCTH[LAR_CTL_CGH] = maxCTH[LAR_STIM_CGH] = pair_max;
      // Since we save the control and swall1 cths as a pair
      // and we already have a CONTROL period max, use the
      // the control max value to scale the swallow ones.
     maxCTH[STIM_SWALL1] = maxCTH[CTL_SWALL1] = maxCTH[CTL_PERIOD];
     maxCTH[STIM_LAREFLEX] = maxCTH[CTL_LAREFLEX] = maxCTH[CTL_PERIOD];
   }
   else if (Norm == PEAK) 
   {
      for (f_bins = bins.begin(); f_bins != bins.end(); ++f_bins)
         for (fbi = f_bins->second.begin(); fbi != f_bins->second.end(); ++fbi)
          {
            curr_period = fbi->first.second;
            if (!fbi->second.IsSparse())
            {
               curr = fbi->second.GetBinPeak();
               if (curr > maxCTH[CtlStims[curr_period]])
               {
                  maxCTH[CtlStims[curr_period]] = curr;
                  if (Debug) cout << "New max mean value for all data in period " << curr_period << " is " << maxCTH[CtlStims[curr_period]] << endl;
               }
            }
         }
        // max for both sets of CTHs
     pair_max = max(maxCTH[CTL_CCO2],maxCTH[STIM_CCO2]);
     maxCTH[CTL_CCO2] = maxCTH[STIM_CCO2] = pair_max;
     pair_max = max(maxCTH[CTL_VCO2],maxCTH[STIM_VCO2]);
     maxCTH[CTL_VCO2] = maxCTH[STIM_VCO2] = pair_max;
     pair_max = max(maxCTH[TB_CTL_CGH],maxCTH[TB_STIM_CGH]);
     maxCTH[TB_CTL_CGH] = maxCTH[TB_STIM_CGH] = pair_max;
     pair_max = max(maxCTH[LAR_CTL_CGH],maxCTH[LAR_STIM_CGH]);
     maxCTH[LAR_CTL_CGH] = maxCTH[LAR_STIM_CGH] = pair_max;
     pair_max = max(maxCTH[CTL_PERIOD],maxCTH[STIM_SWALL1]);
     maxCTH[STIM_SWALL1] = maxCTH[CTL_SWALL1] = pair_max;
     pair_max = max(maxCTH[CTL_PERIOD],maxCTH[STIM_LAREFLEX]);
     maxCTH[STIM_LAREFLEX] = maxCTH[CTL_LAREFLEX] = pair_max;
   }
   else if (Norm == UNIT)
      maxCTH.fill(NumIBins + NumEBins); // so zeroflat == 1.0 for all bins
   else
      maxCTH.fill(1.0);  // no scaling, unscaled spikes/sec values

     // now do scaling for each set of bins
   for (f_bins = bins.begin(); f_bins != bins.end(); ++f_bins)
   {
      for (fbi = f_bins->second.begin(); fbi != f_bins->second.end(); ++fbi)
      {
         curr_period = fbi->first.second;
         if (Norm == MEAN)
            curr = fbi->second.GetBinMeans();
         else if (Norm == PEAK)
            curr = fbi->second.GetBinPeak();
         else if (Norm == UNIT)
         {
            curr = fbi->second.GetBinPeak();
         }
         else
            curr = 1.0;

         if (curr != 0)
         {
            scale = maxCTH[CtlStims[curr_period]] / curr;
            if (Debug) cout << "scale for cth  " << fbi->second.GetChan() << ":   " << scale << endl;
         }
         else
         {
            if (fbi->second.GetTotalSpikes())  // for blanks, mean is always 0
            {
               cout << "That's odd, means for all bins are zero." << endl;
               cout << "chan: " << fbi->first.first << " period: " << fbi->first.second << endl;
               cout << "means:  "  << fbi->second.GetBinMeans() << endl;
               cout << "total spikes: "  << fbi->second.GetTotalSpikes() << endl;
            }
            scale = 1.0;
            if (Debug) cout << "scale for max: " << maxCTH[CtlStims[curr_period]] << " curr: " << curr << " scale: " << scale <<  endl;
         }
         fbi->second.SetScale(scale);  // for each set of bins

         if (Debug) fbi->second.CheckScale(scale,maxCTH[CtlStims[curr_period]]);
         if (Debug) cout << "scale for max: " << maxCTH[CtlStims[curr_period]] << " curr: " << curr << " scale: " << scale <<  endl;
      }
   }
}


// If the current file has cell coordinates, put them in the passed-in
// container.  If not, ensure the container is empty.
void GetCoords(VFilesIter& f_list, AllBins& bins,CellCoords& cell_coords)
{
   ReMapIter       rmap_iter;
   string          fname, line;
   int             old_cell, new_cell;
   int             matches;
   string          name, chan, mchan, ap, rl, dp, dchan, ref, rest;
   int             cnum;
   stringstream    coords;
   AtlasCoords     atlas_coords;
   AtlasCoordsIter atlasiter;

   if (f_list->info_path.empty() || f_list->info_file.empty())
   {
      cell_coords.clear();  // nothing in this file
      return;
   }

   fname = f_list->info_path + "/" + f_list->info_file;
   if (!GetRowsFromSprd(fname, INFO, coords))
      return;
   GetAtlasCoords(fname, atlas_coords);

      // Do we have an optional cell # remapping file?  Some cell #s in older
      // ?dt files were less than 100.  Those files were edited and new cell 
      // numbers were used that were >= 100.  The info spreads were not changed,
      // so we use a remapping table to turn the spread # into the ?dt number.
      // entries are of the form cellnuminsprd-cellnumin?dt file,
      // e.g:   51-203
   if (!f_list->changes_path.empty() && !f_list->changes_file.empty())
   {
      string remap_name= f_list->changes_path + "/" + f_list->changes_file;
      ifstream remap_file(remap_name);
      if (remap_file)
      {
         cout << "Have remapping file " << remap_name << endl;
         while (!remap_file.eof())
         {
            getline(remap_file,line);
            matches = sscanf(line.c_str(),"%d-%d",&old_cell,&new_cell);
            if (matches == 2)
               f_list->remap[old_cell] = new_cell;
         }
      }
      else
         cout << "Could not open remapping file: " << remap_name << endl;
   }

     // extract the stereotaxic coords and other info
   while (!coords.eof())
   {
      getline(coords,line);
cout << line << endl;
if (coords.eof())
   cout << "EOF" << endl;
      stringstream fields(line);
      if (line.empty())
         continue;

      getline(fields,name,'\t');
      getline(fields,chan,'\t');
      getline(fields,ap,'\t');
      getline(fields,rl,'\t');
      getline(fields,dp,'\t');
      getline(fields,dchan,'\t');
      getline(fields,ref,'\t');
      getline(fields,rest,'\n');  // gobble up rest of line

      try { cnum = stoi(chan); }
      catch(...) { cout << "Conversion error for chan " << chan << " in stereotaxic data" << endl; cnum = 0; }
      mchan = chan;   // we may change this shortly, but need original, too
      atlasiter = atlas_coords.find(cnum);
      if (cnum < MIN_CELL || cnum > MAX_CELL)
      {
         rmap_iter = f_list->remap.find(cnum);
          // Not all of the rows in the info grid will have remaps because 
          // they are not neuron channels, so lookup failure is not an error
          // but it does mean this row does not belong in the output
         if (rmap_iter != f_list->remap.end())
            cnum = rmap_iter->second;
         else
            continue;
      }
      StCoords pt;

      try { pt.rl = stod(rl); }
      catch(...) { cout << "Conversion error for rl in stereotaxic data" << endl; }
      try { pt.dp = stod(dp); }
      catch(...) { cout << "Conversion error for dp in stereotaxic data" << endl; }
      try { pt.ap = stod(ap); }
      catch(...) { cout << "Conversion error for ap in stereotaxic data" << endl; }
      if (atlasiter != atlas_coords.end())
      {
         pt.rl_atlas = atlasiter->second.rl;
         pt.ap_atlas = atlasiter->second.ap;
         pt.dp_atlas = atlasiter->second.dp;
      }
      else
      {
         cout << "The atlas channel " << mchan << " was filtered out" << endl;
         pt.rl_atlas = pt.ap_atlas = pt.dp_atlas = 0.0;
      }
     
      pt.name = name;
      pt.mchan = mchan;
      pt.dchan = dchan;
      pt.ref = ref;
      pt.expname = f_list->info_file;
      vector <string> bits;
      string delim="_-";
      Tokenize(f_list->info_file,bits,delim);
      stringstream exp;
      if (bits.size())
      {
         if (isdigit(bits[0][0]))
            exp << bits[0] << "-" << bits[1] << "-" << bits[2] << "_" << bits[3];
         else
            exp << bits[0];
      }
      else
      exp << "noexpname";
      exp >> pt.expbasename;
      cell_coords[cnum]= pt;
      if (Debug) cout << "name " << pt.name << " mchan " << mchan << " chan " << cnum << "  rl: " << pt.rl << " dp: " << pt.dp << " ap: " << pt.ap << " ap atlas: " << pt.ap_atlas << "  dp atlas" << pt.ap_atlas << "  rl atlas " << pt.rl_atlas << endl;

   }
}

/* What we are here for.  Scan the files and gather info and figure average
   firing rates for each neuron, do some calcs, and save it out to file(s).

   Earlier versions of this assume that the period markers signaled the start
   of a period and that the second period was treated the same as the first
   period. This is not the case. 

   The period markers mark the start and end of the control period. This is a
   continuous stream of I and E phases. Spike counts are accumulated and stats
   are calculated, creating CTHs and other pieces of info (just "CTHs" for short).

   The data after the end of the control period is handled differently. The
   input spreadsheet may have one or more markers and intervals for a
   recording. If present, these indicate that the data of interest is in a
   region with a center on the marker and a region before and a region after
   the marker. The region before is a control period. The region after is a
   stimulation period. A CTH is created for each region. Since there can be
   more than one type of stim, is is possible to create 1, 2, 4, or more CTHs
   for a channel. The machinery that accumulates spike counts and the state
   machine that controls progress through the data streams has been adjusted
   accordingly.

   (Adding in ctl/stem periods has made this waaaay out of control, too many states...)
*/
bool ScanFiles(VFiles& f_list, AllBins& bins)
{
   PerPeriod     *per_p = new PerPeriod();
   VFilesIter    f_iter;
   AllBinsIter   f_bins;
   FBinsIter     fbi;
   PerPeriodIter p_data;
   pair <PerPeriodIter,bool> p_pair;
   PerPeriodRet  p_ret;
   PerKey        p_key;
   PerPeriodIter ctl_data;
   string        line;
   bool          done = false;
   bool          in_ctl_period; 
   bool          in_ctlstim;
   bool          in_ephase;
   bool          in_iphase;
   bool          in_ctlphase;
   bool          in_stimphase;
   bool          in_swall1;
   bool          in_swall1_phase;
   bool          in_lareflex;
   bool          in_lareflex_phase;
   bool          found_ok;
   bool          have_cco2=false;
   bool          have_vco2=false;
   bool          have_tb_cgh=false;
   bool          have_lar_cgh=false;
   int           period_markers = 0;
   long          chan;
   streampos     prevline;
   CellCoords    cell_coords;
   CellCoordsIter cc_iter;
   StCoords      coords;
   CTL_STIMS     period_idx;
   ReMapIter     rmap_iter;
   int           code = 0, interval = 0;
   int           start_code, end_code;
   int           stim_cycles=0;
   int           marker_count;
   DataLine      d_line;

   while (!done) // 2 passes, one for cths, one for one tick bin spike counts
   {
        // work the file list
      for (f_iter = f_list.begin(); f_iter != f_list.end(); ++f_iter)
      {
         ifstream file(f_iter->cmplt_name);
         if (!file)
         {
            cout << "The input file appears to have disappeared, skipping. . ." << endl;
            continue;
         }
         if (FilePass1)
            cout << "Opened: " << f_iter->cmplt_name << endl;
         else
            cout << "Re-opened: " << f_iter->cmplt_name << endl;
           // per-file init
         in_ephase = in_iphase = false;
         in_ctlphase = in_stimphase = false;
         in_ctl_period = in_ctlstim = false;
         in_swall1 = in_swall1_phase = false;
         in_lareflex = in_lareflex_phase = false;
         period_markers = 0;
         marker_count = -1;

         NeuroBins::ResetPeriodTimes();
         NeuroBins::SetSecsTick(f_iter->TickVal());
         period_idx = CTL_PERIOD;
         f_bins = bins.find(f_iter->fname);

         if (f_bins == bins.end()) // this is a bug
         {
            cout << "WARNING: You have a bug. " << endl << "Did not find bin record for file " << f_iter->cmplt_name << endl << "in " << __FILE__ << __LINE__ << endl;
            continue;
         }

         if (FilePass1)
         {
            cell_coords.clear();                 // get optional stereotaxic coords 
            GetCoords(f_iter,bins,cell_coords);  // only on 1st pass
         }
         while (!file.eof())
         {
            prevline = file.tellg();  // may need to seek back to this line
            getline(file,line);
            if (line.empty())
               continue;

            d_line = line;
            chan = d_line.GetChan();                // need to remap chan #?
            rmap_iter = f_iter->remap.find(chan);
            if (rmap_iter != f_iter->remap.end())
            {
               chan = rmap_iter->second;
               d_line.RemapChan(line,chan);
            } 
            if (!in_ctl_period)    // if set, this also means we are passed the ctl period
            {
               if (d_line.IsMarker(f_iter->pmark))  // found 1st (ctl) period start
               {
                  in_ctl_period = true;
                  ++period_markers;
                  NeuroBins::ResetPeriodTimes();
                  NeuroBins::PeriodStart(d_line.GetTime());
                    // create / fetch this period's info class
                  p_key = make_pair(f_iter->fname,period_idx);
                  p_data = per_p->emplace(make_pair(p_key,PeriodInfo(period_idx))).first;
                  if (FilePass1)
                     cout << "Found Control Period start at " << f_iter->TimeVal(d_line.GetTime())/1000 << " sec." << endl;
               }
               else
                  continue;   // keep looking for ctl period start
            }
            else if (d_line.IsMarker(f_iter->pmark) && !in_ctlstim)  // end of ctl period
            {
               if (FilePass1)
               {
                  cout << "Control ";
                  NeuroBins::PeriodEnd(d_line.GetTime());
                  NeuroBins::PrintPeriodTime();
                  NeuroBins::ResetPeriodTimes();
               }
               ++period_markers;
               in_ctlstim = true;             // second part of file
               if (FilePass1)
               {
                  PeriodCalcs(p_data->second);
                  UpdateStats(f_bins,p_data->second);
               }
               if (file.eof())  // if hit eof, prevline will be last event time 
               {
                  file.clear();  // have to do this so we can seek
                  file.seekg(prevline);
                  getline(file,line);
                  if (!line.empty())
                  {
                     DataLine d_line(line);
                     NeuroBins::PeriodEnd(d_line.GetTime());
                     NeuroBins::PrintPeriodTime();
                  }
               }
                 // if we find ctl/stim markers
               in_ephase = in_iphase = false;      // setup for ctl/stim period
               in_ctlphase = in_stimphase = false;
               if (f_iter->have_markers)
               {
                  if (FilePass1)
                     GetNumCtlStimMarkers(file,f_iter); // count # of these (if any)
               }
               else
                  cout << "No ctl/stim/swallow markers in this file." << endl;
            }
            else if (!in_lareflex && !in_swall1 && !in_ephase && !in_iphase && !in_ctlstim && d_line.IsIPhase())
            {
                 // not in any phase and is still in 1st (ctl) period
               file.seekg(prevline);  // back up so i line is next in stream
               if (NeuroBins::FindPhaseMarks(file,f_iter->pmark,f_bins,p_data->second))
               {
                  if (FilePass1)
                     StartNewCycle(f_bins,p_data->second);
                  else
                     StartNew_OT_Cycle(f_bins,p_data->second);
                  in_iphase= true;
                }
            }
            else if (f_iter->cco2stim_code && d_line.IsMarker(f_iter->cco2stim_code)) 
            {
                  // look for ctl/stim markers (start of 2nd-nth pseudo-periods)
               if (FilePass1 && Debug)
                  cout << "Found C co2 marker code " << f_iter->cco2stim_code << endl;
               code = f_iter->cco2stim_code;
               interval = f_iter->cco2stim_interval;
               if (!in_ctlphase) // process interval before/up to marker
               {
                  if (FilePass1 && !have_cco2)
                  {
                     NeuroBins::ResetPeriodTimes();
                     NeuroBins::PeriodStart(d_line.GetTime());
                     cout << "Found CCO2 period start at " << f_iter->TimeVal(d_line.GetTime())/1000 << " sec." << endl;
                  }
                  in_ctlphase = true;
                  if (marker_count == -1)
                     marker_count = f_iter->cco2marker_count;
                  period_idx = CtlStims[CTL_CCO2];
                  p_key = make_pair(f_iter->fname,period_idx); // create new one or look up existing one
                  p_data = per_p->emplace(make_pair(p_key,PeriodInfo(period_idx))).first;
                  ctl_data = p_data;
                  file.seekg(prevline);  // back up so marker line is next in stream
                  if (NeuroBins::FindCtlStimMarks(file,code, -interval,f_bins,p_data->second))
                  {
                     in_iphase=true;
                     in_ephase=false;
                     have_cco2 = true;
                  }
                  else
                     cout << "WARNING: failed to find start of Cco2 interval" << endl;
               }
               else
               {
                  in_ctlphase = false; // process interval marker/after 
                  in_stimphase = true;
                  period_idx = CtlStims[STIM_CCO2];
                  p_key = make_pair(f_iter->fname,period_idx);
                  p_data = per_p->emplace(make_pair(p_key,PeriodInfo(period_idx))).first;
                  file.seekg(prevline);  // back up so marker line is next in stream
                  if (NeuroBins::FindCtlStimMarks(file,code,interval,f_bins,p_data->second))
                  {
                     in_iphase=true;
                     in_ephase=false;
                  }
                  else
                     cout << "WARNING: failed to find end of Cco2 interval" << endl;
                  --marker_count;
                  if (marker_count == 0 && FilePass1)
                  {
                     int cycle = f_iter->cco2stim_interval * 2 * f_iter->cco2marker_count;
                     cout << "CCO2 Total Cycle Time: " << cycle << " sec.  " << cycle/60.0 << " min." << endl;
                     cout << "CCO2 Elapsed ";
                     NeuroBins::PeriodEnd(d_line.GetTime());
                     NeuroBins::PrintPeriodTime();
                     NeuroBins::ResetPeriodTimes();
                  }
               }
            }
            else if (f_iter->vco2stim_code && d_line.IsMarker(f_iter->vco2stim_code)) 
            {
               if (FilePass1 && Debug)
                  cout << "Found V co2 marker code " << f_iter->vco2stim_code << endl;
               code = f_iter->vco2stim_code;
               interval = f_iter->vco2stim_interval;
               if (!in_ctlphase)
               {
                  if (FilePass1 && !have_vco2)
                  {
                     NeuroBins::ResetPeriodTimes();
                     NeuroBins::PeriodStart(d_line.GetTime());
                     cout << "Found VCO2 period start at " << f_iter->TimeVal(d_line.GetTime())/1000 << " sec." << endl;
                  }
                  in_ctlphase = true;
                  if (marker_count == -1)
                     marker_count = f_iter->vco2marker_count;
                  period_idx = CtlStims[CTL_VCO2];
                  p_key = make_pair(f_iter->fname,period_idx);
                  p_data = per_p->emplace(make_pair(p_key,PeriodInfo(period_idx))).first;
                  ctl_data = p_data;
                  file.seekg(prevline); // back up so marker line is next in stream
                  if (NeuroBins::FindCtlStimMarks(file,code, -interval,f_bins,p_data->second))
                  {
                     in_iphase=true;
                     in_ephase=false;
                     have_vco2 = true;
                  }
                  else
                     cout << "WARNING: failed to find start of VCO2 interval" << endl;
               }
               else
               {
                  in_ctlphase = false;
                  in_stimphase = true;
                  period_idx = CtlStims[STIM_VCO2];
                  p_key = make_pair(f_iter->fname,period_idx);
                  p_data = per_p->emplace(make_pair(p_key,PeriodInfo(period_idx))).first;
                  file.seekg(prevline);
                  if (NeuroBins::FindCtlStimMarks(file,code,interval,f_bins,p_data->second))
                  {
                     in_iphase=true;
                     in_ephase=false;
                  }
                  else
                     cout << "WARNING: failed to find end of VCO2 interval" << endl;
                  --marker_count;
                  if (marker_count == 0 && FilePass1)
                  {
                     int cycle = f_iter->vco2stim_interval * 2 * f_iter->vco2marker_count;
                     cout << "VCO2 Total Cycle Time: " << cycle << " sec.  " << cycle/60.0 << " min." << endl;
                     cout << "VCO2 Elapsed ";
                     NeuroBins::PeriodEnd(d_line.GetTime());
                     NeuroBins::PrintPeriodTime();
                     NeuroBins::ResetPeriodTimes();
                  }
               }
            }
            else if (f_iter->tb_cgh_stim_code && d_line.IsMarker(f_iter->tb_cgh_stim_code)) 
            {
               if (FilePass1 && Debug)
                  cout << "Found TB COUGH marker code " << f_iter->tb_cgh_stim_code << endl;
               code = f_iter->tb_cgh_stim_code;
               if (!in_ctlphase)
               {
                  if (FilePass1 && !have_tb_cgh)
                  {
                     NeuroBins::ResetPeriodTimes();
                     NeuroBins::PeriodStart(d_line.GetTime());
                     cout << "Found TB COUGH period start at " << f_iter->TimeVal(d_line.GetTime())/1000 << " sec." << endl;
                  }
                  in_ctlphase = true;
                  interval = f_iter->tb_cgh_ctl_interval;
                  stim_cycles = 0;
                  if (marker_count == -1)
                     marker_count = f_iter->tb_cghmarker_count;
                  period_idx = CtlStims[TB_CTL_CGH];
                  p_key = make_pair(f_iter->fname,period_idx);
                  p_data = per_p->emplace(make_pair(p_key,PeriodInfo(period_idx))).first;
                  ctl_data = p_data;
                  file.seekg(prevline); // back up so marker line is next in stream
                  if (NeuroBins::FindCtlStimMarks(file,code, -interval,f_bins,p_data->second))
                  {
                     in_iphase=true;
                     in_ephase=false;
                     have_tb_cgh = true;
                  }
                  else
                     cout << "WARNING: failed to find start of TB COUGH interval" << endl;
               }
               else
               {
                  in_ctlphase = false;
                  in_stimphase = true;
                  interval = 0;
                  stim_cycles = f_iter->tb_cgh_stim_cycles;
                  ++stim_cycles;  // we always skip the 1st i0-e0-i1 cycle for cough
                  period_idx = CtlStims[TB_STIM_CGH];
                  p_key = make_pair(f_iter->fname,period_idx);
                  p_data = per_p->emplace(make_pair(p_key,PeriodInfo(period_idx))).first;
                  file.seekg(prevline);
                  if (NeuroBins::FindCtlStimMarks(file,code,interval,f_bins,p_data->second,stim_cycles))
                  {
                     --stim_cycles;  // this is the one we want to start at
                     getline(file,line); // skip this i0
                     if (NeuroBins::FindCtlStimMarks(file,code,interval,f_bins,p_data->second,stim_cycles))
                     {
                        in_iphase=true;
                        in_ephase=false;
                     }
                     else
                        cout << "WARNING: failed to find start of TB COUGH stim interval" << endl;
                  }
                  else
                     cout << "WARNING: failed to find end of TB COUGH interval" << endl;

                  --marker_count;
                  if (marker_count == 0 && FilePass1)
                  {
                     cout << "COUGH Elapsed ";
                     NeuroBins::PeriodEnd(d_line.GetTime());
                     NeuroBins::PrintPeriodTime();
                     NeuroBins::ResetPeriodTimes();
                  }
               }
            }
            else if (f_iter->lar_cgh_stim_code && d_line.IsMarker(f_iter->lar_cgh_stim_code)) 
            {
               if (FilePass1 && Debug)
                  cout << "Found LAR COUGH marker code " << f_iter->lar_cgh_stim_code << endl;
               code = f_iter->lar_cgh_stim_code;
               if (!in_ctlphase)
               {
                  if (FilePass1 && !have_lar_cgh)
                  {
                     NeuroBins::ResetPeriodTimes();
                     NeuroBins::PeriodStart(d_line.GetTime());
                     cout << "Found LAR COUGH period start at " << f_iter->TimeVal(d_line.GetTime())/1000 << " sec." << endl;
                  }
                  in_ctlphase = true;
                  interval = f_iter->lar_cgh_ctl_interval;
                  stim_cycles = 0;
                  if (marker_count == -1)
                     marker_count = f_iter->lar_cghmarker_count;
                  period_idx = CtlStims[LAR_CTL_CGH];
                  p_key = make_pair(f_iter->fname,period_idx);
                  p_data = per_p->emplace(make_pair(p_key,PeriodInfo(period_idx))).first;
                  ctl_data = p_data;
                  file.seekg(prevline); // back up so marker line is next in stream
                  if (NeuroBins::FindCtlStimMarks(file,code, -interval,f_bins,p_data->second))
                  {
                     in_iphase=true;
                     in_ephase=false;
                     have_lar_cgh = true;
                  }
                  else
                     cout << "WARNING: failed to find start of LAR COUGH interval" << endl;
               }
               else
               {
                  in_ctlphase = false;
                  in_stimphase = true;
                  interval = 0;
                  stim_cycles = f_iter->lar_cgh_stim_cycles;
                  ++stim_cycles; // skip 1st
                  period_idx = CtlStims[LAR_STIM_CGH];
                  p_key = make_pair(f_iter->fname,period_idx);
                  p_data = per_p->emplace(make_pair(p_key,PeriodInfo(period_idx))).first;
                  file.seekg(prevline);
                  if (NeuroBins::FindCtlStimMarks(file,code,interval,f_bins,p_data->second,stim_cycles))
                  {
                     --stim_cycles;  // this is the one we want to start at
                     getline(file,line); // skip this i0
                     if (NeuroBins::FindCtlStimMarks(file,code,interval,f_bins,p_data->second,stim_cycles))
                     {
                        in_iphase=true;
                        in_ephase=false;
                     }
                     else
                        cout << "WARNING: failed to find start of LAR COUGH stim interval" << endl;
                  }
                  else
                     cout << "WARNING: failed to find end of LAR COUGH ctl interval" << endl;
                  --marker_count;
                  if (marker_count == 0 && FilePass1)
                  {
                     cout << "COUGH Elapsed ";
                     NeuroBins::PeriodEnd(d_line.GetTime());
                     NeuroBins::PrintPeriodTime();
                     NeuroBins::ResetPeriodTimes();
                  }
               }
            }
            else if (f_iter->swall1_start_code && d_line.IsMarker(f_iter->swall1_start_code)) 
            {
               if (FilePass1 && Debug)
                  cout << "Found SWALLOW 1 start marker code " << f_iter->swall1_start_code << endl;
               start_code = f_iter->swall1_start_code;
               end_code = f_iter->swall1_end_code;
               if (!have_swall1) // one-time init on entering swall1 period(s)
               {
                  if (FilePass1) 
                  {
                     NeuroBins::ResetPeriodTimes();
                     NeuroBins::PeriodStart(d_line.GetTime());
                     cout << "Found SWALLOW1 period start at " << f_iter->TimeVal(d_line.GetTime())/1000 << " sec." << endl;
                     have_swall1 = true;
                     in_swall1 = true;
                  }
               }
               if (!in_swall1_phase) // starting new phase
               {
                  if (marker_count == -1)
                     marker_count = min(f_iter->swall1marker_s_count,f_iter->swall1marker_e_count);
                  in_swall1_phase = true;
                  period_idx = CtlStims[STIM_SWALL1];
                  p_key = make_pair(f_iter->fname,period_idx);
                  p_data = per_p->emplace(make_pair(p_key,PeriodInfo(period_idx))).first;
                  file.seekg(prevline); // back up so marker line is next in stream
                  if (NeuroBins::FindSwall1Marks(file,start_code,end_code,p_data->second))
                  {
                     getline(file,line); // skip marker
                  }
                  else
                     cout << "WARNING: failed to find SWALL1 interval" << endl;
                  if (FilePass1)
                     StartNewCycle(f_bins,p_data->second);
                  else
                     StartNew_OT_Cycle(f_bins,p_data->second);
               }
            }
            else if (f_iter->swall1_end_code && d_line.IsMarker(f_iter->swall1_end_code)) 
            {
               if (FilePass1 && Debug)
                  cout << "Found SWALLOW 1 end marker code " << f_iter->swall1_end_code << endl;
               if (in_swall1_phase)
               {
                  ++p_data->second.num_cycles;
                  EndCurrCycle(f_bins,p_data->second);
               }
               else
                  cout << "WARNING: Found end of SWALLOW1 interval with no start." << endl;
               in_swall1_phase = false;
               --marker_count;
               if (marker_count == 0 && FilePass1)
               {
                  cout << "SWALLOW1 Elapsed ";
                  NeuroBins::PeriodEnd(d_line.GetTime());
                  NeuroBins::PrintPeriodTime();
                  NeuroBins::ResetPeriodTimes();
                  PeriodCalcs(p_data->second);
                  UpdateStats(f_bins,p_data->second);
               }
            }
            else if (f_iter->lareflex_start_code && d_line.IsMarker(f_iter->lareflex_start_code)) 
            {
               if (FilePass1 && Debug)
                  cout << "Found LAREFLEX start marker code " << f_iter->lareflex_start_code << endl;
               start_code = f_iter->lareflex_start_code;
               end_code = f_iter->lareflex_end_code;
               if (!have_lareflex) // one-time init on entering lareflex period(s)
               {
                  if (FilePass1) 
                  {
                     NeuroBins::ResetPeriodTimes();
                     NeuroBins::PeriodStart(d_line.GetTime());
                     cout << "Found LAREFLEX period start at " << f_iter->TimeVal(d_line.GetTime())/1000 << " sec." << endl;
                     have_lareflex = true;
                     in_lareflex= true;
                  }
               }
               if (!in_lareflex_phase) // starting new phase
               {
                  if (marker_count == -1)
                     marker_count = min(f_iter->lareflexmarker_s_count,f_iter->lareflexmarker_e_count);
                  in_lareflex_phase = true;
                  period_idx = CtlStims[STIM_LAREFLEX];
                  p_key = make_pair(f_iter->fname,period_idx);
                  p_data = per_p->emplace(make_pair(p_key,PeriodInfo(period_idx))).first;
                  file.seekg(prevline); // back up so marker line is next in stream
                     // the swall1 and larflex finding is the same
                  if (NeuroBins::FindSwall1Marks(file,start_code,end_code,p_data->second))
                  {
                     getline(file,line); // skip marker
                  }
                  else
                     cout << "WARNING: failed to find LAREFLEX interval" << endl;
                  if (FilePass1)
                     StartNewCycle(f_bins,p_data->second);
                  else
                     StartNew_OT_Cycle(f_bins,p_data->second);
               }
            }
            else if (f_iter->lareflex_end_code && d_line.IsMarker(f_iter->lareflex_end_code)) 
            {
               if (FilePass1 && Debug)
                  cout << "Found LAREFLEX 1 end marker code " << f_iter->lareflex_end_code << endl;
               if (in_lareflex_phase)
               {
                  ++p_data->second.num_cycles;
                  EndCurrCycle(f_bins,p_data->second);
               }
               else
                  cout << "WARNING: Found end of LAREFLEX interval with no start." << endl;
               in_lareflex_phase = false;
               --marker_count;
               if (marker_count == 0 && FilePass1)
               {
                  cout << "LAREFLEX Elapsed ";
                  NeuroBins::PeriodEnd(d_line.GetTime());
                  NeuroBins::PrintPeriodTime();
                  NeuroBins::ResetPeriodTimes();
                  PeriodCalcs(p_data->second);
                  UpdateStats(f_bins,p_data->second);
               }
            }
            else // accumulate neuron firing counts for current period,phase
            {
               if ((!in_lareflex || !in_swall1) && in_iphase && d_line.IsEPhase())
               {
                  in_ephase= true;
                  in_iphase= false;
                  if (FilePass1)
                     StartEPhase(f_bins,p_data->second); // pick up boundary I events
                  else
                     StartE_OT_Phase(f_bins,p_data->second);
               }
               else if ((!in_lareflex || !in_swall1) && in_ephase && d_line.IsIPhase()) // end of e, set up for next cycle
               {
                  --stim_cycles;
                  in_iphase = true;
                  in_ephase = false;
                  file.seekg(prevline);  // back up so i marker is next in stream
                  if (FilePass1)
                  {
                     ++p_data->second.num_cycles;
                     EndCurrCycle(f_bins,p_data->second);
                  }
                  if (in_ctlphase)
                  {
                     found_ok = NeuroBins::FindCtlStimMarks(file,code,-interval,f_bins,p_data->second);
                  }
                  else if (in_stimphase)
                  {   // set up phase vars for next i0-e0-i1
                     found_ok = NeuroBins::FindCtlStimMarks(file,code,interval,f_bins,p_data->second,stim_cycles);
                     if (!found_ok)  // hit end of current ctl/stim period
                     {
                        if (marker_count == 0)
                        {
                           if (FilePass1)
                           {
                              PeriodCalcs(p_data->second);   // done with all periods for 
                              PeriodCalcs(ctl_data->second); // this CTH pair
                              UpdateStats(f_bins,p_data->second);
                              UpdateStats(f_bins,ctl_data->second);
                           }
                           marker_count = -1;
                        }
                        in_stimphase = false;
                        in_iphase = in_ephase = false;
                        stim_cycles = 0;
                        interval = 0;
                     }
                  }
                  else
                     found_ok = NeuroBins::FindPhaseMarks(file,f_iter->pmark,f_bins,p_data->second);
                  if (found_ok)
                  {
                     if (FilePass1)
                        StartNewCycle(f_bins,p_data->second);
                     else
                        StartNew_OT_Cycle(f_bins,p_data->second);
                 }
                 if (Debug) cout << "# of chans: " << f_bins->second.size() << endl;
               }
               else if ((in_iphase || in_ephase || in_swall1_phase || in_lareflex_phase) && d_line.IsNeuron())
               {
                  long chan = d_line.GetChan();     // a spike, count it
                  FBinKey key(chan,period_idx);
                  fbi = f_bins->second.find(key);
                  if (fbi == f_bins->second.end() && FilePass1)  // 1st spike 1st pass?
                  {
                     cc_iter = cell_coords.find(chan);
                     if (cc_iter != cell_coords.end())
                        coords = cc_iter->second;
                     else
                        coords = StCoords();
                     fbi = f_bins->second.insert(        // create a record
                               FBins::value_type(key,
                               NeuroBins(
                               chan,
                               coords,
                               period_idx,
                               NeuroBins::GetETime(),
                               NeuroBins::GetITime(),
                               NeuroBins::GetSecsTick()))).first;
                  }
                  if (Debug) cout << "[chan] " << chan << endl;

                  if (FilePass1)  // 1st pass, create cths and data for 2nd pass
                  {
                     if (!in_swall1_phase && !in_lareflex_phase)
                     {
                        if (in_iphase)
                           fbi->second.UpdateI(d_line);
                        else 
                           fbi->second.UpdateE(d_line);
                     }
                     else
                        fbi->second.UpdateS(d_line);
                  }
                  else
                  {
                     if (fbi != f_bins->second.end()) // not an error if not in map,
                     {                                // some filtered out as sparse
                        if (!in_swall1_phase && !in_lareflex_phase)
                        {
                           if (in_iphase)
                              fbi->second.UpdateI_OT(d_line);
                           else
                              fbi->second.UpdateE_OT(d_line);
                        }
                        else
                           fbi->second.UpdateS_OT(d_line);
                     }
                     else
                        cout << "Test in scanfiles, are spareses really deleted?" << endl;
                  }
               }
            } // else nothing we want
         }

         if (Debug) PrintBins(f_bins);

         if (period_markers == 0)
            cout << "WARNING:  No control period markers found" << endl;
         else if (period_markers == 1)
            cout << "WARNING:  Did not find end of control period marker" << endl;
         file.close();
      } // end of reading all files

      if (FilePass1)
      {
         CalcsForAllFiles(bins,*per_p);
         if (Debug) p_data->second.Report();
         // Now that know the shortest phases for each phase for every channel in
         // each period, calculate the one tick bin spike counts for later curve-fitting.
         SetupForTicks(bins,*per_p);
         FilePass1 = false;  // set up for second pass through files
         cout << "Rescanning files for one tick bin processing." << endl;
      }
      else
         done = true;
   }

    // put ctl/stim CTHs in separate files

   P_LIST per;
   if (AllP)
   {
      per = CtlStims;
      SaveResults(bins,per,*per_p);
   }
   else
   {
      per.push_back(CTL_PERIOD);
      SaveResults(bins,per,*per_p);
      per.clear();
      NumSparse = 0;   // kind of a hack, but we don't mark ctl/stim CTHs as sparse 
      if (have_cco2)
      {
         per.push_back(CTL_CCO2);
         per.push_back(STIM_CCO2);
         SaveResults(bins,per,*per_p);
         per.clear();
      }
      if (have_vco2)
      {
         per.push_back(CTL_VCO2);
         per.push_back(STIM_VCO2);
         SaveResults(bins,per,*per_p);
         per.clear();
      }
      if (have_tb_cgh)
      {
         per.push_back(TB_CTL_CGH);
         per.push_back(TB_STIM_CGH);
         SaveResults(bins,per,*per_p);
         per.clear();
      }
      if (have_lar_cgh)
      {
         per.push_back(LAR_CTL_CGH);
         per.push_back(LAR_STIM_CGH);
         SaveResults(bins,per,*per_p);
         per.clear();
      }
      if (have_swall1)
      {
         UnMarkSparses(bins);
         per.push_back(CTL_SWALL1);
         per.push_back(CTL_PERIOD);
         per.push_back(STIM_SWALL1);
         SaveResults(bins,per,*per_p);
         per.clear();

         RemoveBlanks(bins, STIM_SWALL1);
         per.push_back(STIM_SWALL1);
         SaveResults(bins,per,*per_p);
//         MarkSparses(bins); // leave sparseness like we found it
         per.clear();
      }
      if (have_lareflex)
      {
         UnMarkSparses(bins);
         per.push_back(CTL_LAREFLEX);
         per.push_back(CTL_PERIOD);
         per.push_back(STIM_LAREFLEX);
         SaveResults(bins,per,*per_p);
         per.clear();

         RemoveBlanks(bins, STIM_LAREFLEX);
         per.push_back(STIM_LAREFLEX);
         SaveResults(bins,per,*per_p);
//         MarkSparses(bins); // leave sparseness like we found it
         per.clear();
      }
   }
   delete per_p;
   return true;
}


/* Print out the array of info about files
*/
void dump_list(VFiles& f_list)
{
   for_each(f_list.begin(), f_list.end(), mem_fun_ref(&FileFacts::print));
}

// If a channel has fewer than this many spikes, we do not save it.
void CalcReject()
{
   double bins = NumEBins + NumIBins;
   RejectThresh = ceil (log10(0.05/bins) / log10((bins-1)/bins)); 
// testing, for larger # of bins, need a lot more spikes 
//RejectThresh = 100;
   if (Debug) cout << "Reject Threshold is " << RejectThresh << endl;
}


int main(int argc, char* argv[])
{
   string  next_file;
   VFiles    *file_list = new VFiles();
   AllBins   *bins = new AllBins();

   cout << argv[0] << " Version " << VERSION << endl;

   ParseArgs(argc,argv);
   if (!GetFileInfo(*file_list,*bins))
      exit(1);
   cout << "In Name:   " << InName << "\nOut Name:  " <<  OutName << "\nUsing: " << NumEBins << " E Bins\nUsing: " << NumIBins << " I Bins" << endl;

   // These are global and constant once we start processing
   NeuroBins::SetNumEBins(NumEBins);
   NeuroBins::SetNumIBins(NumIBins);
   CalcReject();

   if (Debug) dump_list(*file_list);

   ScanFiles(*file_list,*bins);

   delete bins;
   delete file_list;

   cout << "DONE!" << endl;
}

