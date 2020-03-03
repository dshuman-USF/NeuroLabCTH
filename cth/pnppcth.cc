//g++  -Wall -std=gnu++11 -g -O0 -o pnppcth pnppcth.cc -lboost_program_options
//g++  -Wall -std=gnu++11 -ggdb3 -o pnppcth pnppcth.cc -lboost_program_options
// phase normalized Poisson preserving CTH


#define GLIBCXX_USE_CXX11_ABI 1

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <errno.h>
#include <error.h>
#include <cmath>
#include <boost/program_options.hpp>


namespace po = boost::program_options;

using namespace std;

struct Cycle {
  int pdur[2];
  int start;
  vector<int> pst[2];
};

void printcycle(Cycle&c)
{
   cout << "start: " << c.start << " I: " << c.pdur[0] << "  E: " << c.pdur[1] << endl;
#if 0
   cout << "pdur[0]: " << c.pdur[0] << "  pdur[1]: " << c.pdur[1] << "  start: " << c.start << " pst0: ";
   for (vector<int>::iterator it = c.pst[0].begin(); it != c.pst[0].end(); ++it)
      cout << *it << "  " ;
   cout << "  pst1: ";
   for (vector<int>::iterator it = c.pst[1].begin(); it != c.pst[1].end(); ++it)
      cout << *it << "  " ;
   cout << endl;
#endif
}

vector<Cycle>
get_cycles (string filename, int cell, int ccode, int pcode)
{
  ifstream f(filename);
  if (!f)
    error (1, errno, "error opening %s", filename.c_str());
  string s;
  getline (f, s);
  getline (f, s);
  vector<int> ct;               // cell times
  vector<int> cct;              // cycle code times
  vector<int> pct;              // phase code times
  while (getline (f, s)) {
    int id = atoi (s.substr(0,5).c_str());
    int t = atoi (s.c_str() + 5);

    if (id == cell)
      ct.push_back (t);
    if (id == ccode)
      cct.push_back (t);
    if (id == pcode)
      pct.push_back (t);
  }
  f.close ();

  vector<Cycle> cycles;
  cycles.reserve (cct.size());
  for (size_t i = 1; i < cct.size(); i++) {
    int start = cct[i - 1];
    int end = cct[i];

    // find phase marker for cycle
    auto ip0 = lower_bound (pct.begin(), pct.end(), start);
    auto ipn = lower_bound (pct.begin(), pct.end(), end);
    // there can be only one
    if (ipn - ip0 != 1) {
      cerr << "cycle starting at " << cct[i - 1]
           << " has " << ipn - ip0
           << " phase markers, skipping\n";
      continue;
    }
    Cycle c;
    c.pdur[0] = *ip0 - start; // pdur[0] is length of phase0 (i)
    c.pdur[1] = end - *ip0;   // pdur[1] is length of phase1 (e)
    c.start = start;
    // find spikes for phase 0
    auto i0s0 = lower_bound (ct.begin(), ct.end(), start); // tick for 1st i spike in phase
    auto i0sn = lower_bound (ct.begin(), ct.end(), *ip0);  // tick for 1st e spike in phase
    c.pst[0].insert (c.pst[0].end(), i0s0, i0sn);
    for (auto &t : c.pst[0])
      t -= start;
    // find spikes for phase 1
    auto i1s0 = lower_bound (ct.begin(), ct.end(), *ip0);
    auto i1sn = lower_bound (ct.begin(), ct.end(), end);
    c.pst[1].insert (c.pst[1].end(), i1s0, i1sn);
    for (auto &t : c.pst[1])
      t -= *ip0;
//printcycle(c);
    cycles.push_back (move (c));
  }
  return cycles;
}

int
pick_short (vector<size_t> &h)
{
  size_t count = h[1];
  size_t maxcnt = h[1];
  int max_at = 1;

  for (size_t lo = 2; 2 * lo < h.size(); ++lo) {
    count -= h[lo - 1];
    for (size_t i = 2 * (lo - 1); i < 2 * lo; ++i) {
      count += h[i];
    }
    if (count > maxcnt) {
      maxcnt = count;
      max_at = lo;
    }
  }
  return max_at;
}

#pragma GCC diagnostic ignored "-Wformat-security"
int
dur_hist (vector<Cycle> &cycles, int phase, string &ofile,
          bool show)
{
  int max = 0;
  for (auto &c : cycles) {
    int d = c.pdur[phase];
    if (d > max){
       max = d;
    }
  }

  vector<size_t> h (max + 1);
  for (auto &c : cycles) {
    int d = c.pdur[phase];
    if (d >= 0 && d <= max)
      h[d]++;
  }

//cout << "H vector start, " << h.size() << " elements." << endl;
// int row = 0;
// for (vector<size_t>::iterator it = h.begin(); it < h.end(); ++it, ++row)
//   if (*it != 0)
//      cout << row << "  " << *it << endl;
//cout << "H vector end" << endl;

    int mindur = pick_short (h);
    cout << "phase " << phase << ": " << mindur
         << " to " << 2 * mindur - 1 << " ticks\n";
  if (!show) return mindur;

  ofstream f (ofile);
  for (int i = 0; i <= max; i++)
    f << i << " " << h[i] << endl;
  f.flush();
  string cmd = (string ("plot '") + ofile
                + "' w his; pause mouse close\n");
  FILE *g = popen ("gnuplot", "w");
  fprintf (g, cmd.c_str());
  fclose (g);
  return mindur;
}

vector<Cycle>
cull_cycles (vector<Cycle> &cycles, int* pmin)
{
  vector<Cycle> culled;
  for (auto &c : cycles)
  {
//cout << c.pdur[0] << "  " << pmin[0] << "  " << pmin[0] * 2 << endl;
//cout << c.pdur[1] << "  " << pmin[1] << "  " << pmin[1] * 2 << endl;
    if (c.pdur[0] >= pmin[0] && c.pdur[0] < pmin[0] * 2
        && c.pdur[1] >= pmin[1] && c.pdur[1] < pmin[1] * 2)
//    if (c.pdur[0] >= pmin[0] && c.pdur[1] >= pmin[1])
      culled.push_back (move (c));
    else
    {
//      cout << "skip cycle starting at: " << c.start << endl;
    }
  }
  return culled;
}

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

vector<int> 
cth (vector<Cycle> &cycles, int *pmin, int phase)
{
  vector<int> h (pmin[phase]);
  for (auto &c : cycles)
    for (int j : c.pst[phase]) {

      int bin = get_bin (j, c.pdur[phase], pmin[phase]);
//      if (phase == 0) // just i's, use 1 for just e's
//      {
//        cout << "j: " << j << " n_p: " << c.pdur[phase] << " n_b: " << pmin[phase] << " start: " << c.start;
//        cout << " bin is: " << bin << endl;
//      }
      assert (bin < int (h.size()));

// check bin value
//      int n_b = pmin[phase];
//      int i = bin;
//      int n_p = c.pdur[phase];
//
//      double c_i = (0.5 + i) / n_b * n_p;
//      double c_s = j + .5;
//
//      cout << setprecision (17);
//      if ((i >= 0 &&  abs (c_s - nextafter (c_i, c_i + 1)) > .5)
//          || (i < 0 && abs (c_s - nextafter (c_i, c_i - 1)) < .5))
//        cout << "n_p: " << n_p << endl
//             << "n_b: " << n_b << endl
//             << "i: " << i << endl
//             << "j: " << j << endl
//             << "c_i: " << c_i << endl
//             << "c_s: " << c_s << endl
//             << "abs (c_s - c_i): " << abs (c_s - c_i) << endl
//             << "nabs (c_s - nextafter (c_i, c_i + 1)): "
//             << abs (c_s - nextafter (c_i, c_i + 1)) << endl
//             << "nabs (c_s - nextafter (c_i, c_i - 1)): "
//             << abs (c_s - nextafter (c_i, c_i - 1)) << endl
//
//             << endl;

      if (bin >= 0)
        h[bin]++;
      else
         cout << "bin size is out of bounds" << endl;
    }
  return h;
}

void
write_hist (ofstream &f, vector<int> h)
{
  for (int n : h)
    f << n << endl;
  f << -1 << endl;
}

int
main (int argc, char **argv)
{
  string infile = "in.edt";
  string ofile = "pnppcth.out";
  int cell = 501;
  int ccode = 97;
  int pcode = 98;
  bool phist[2];
  po::options_description opt("\n Options");
  try {using namespace po;
    opt.add_options()
      ("help,h"  ,"show this message")
      ("in,i", value (&infile), "input file (in.edt)")
      ("out,o", value (&ofile), "output file (pnppcth.out)")
      ("cell,c", value (&cell), "cell id (501)")
      ("ccode,m", value (&ccode), "cycle marker (97)")
      ("pcode,p", value (&pcode), "phase marker (98)")
      ("p0hist,0", "phase 0 histogram")
      ("p1hist,1", "phase 1 histogram")
      ;
    variables_map vm;
    store(parse_command_line(argc, argv, opt), vm);
    notify(vm);
    phist[0] = vm.count ("p0hist");
    phist[1] = vm.count ("p1hist");
    if (vm.count("help")) {  
      cout << opt << "\n";
      return 0;
    }
  } catch(exception& e) {cerr << e.what() << "\n" << opt << endl;}
  cout << "cell: " << cell << endl;

  vector<Cycle> cycles = get_cycles (infile, cell, ccode, pcode);

  int pmin[2];                  // minimum duration of phase
  for (int p : {0,1})
    pmin[p] = dur_hist (cycles, p, ofile, phist[p]);

  size_t before = cycles.size();
  cycles = cull_cycles (cycles, pmin);
  size_t after = cycles.size();

//  vector<Cycle> cycles1  = cull_cycles (cycles, pmin);
//  size_t after = cycles1.size();

  cout << "using " << 100. * after / before << " % of cycles\n";

  vector<vector<int> > hist(2);
  ofstream f (ofile);
  for (int p : {0,1})
    write_hist (f, cth (cycles, pmin, p));

  return 0;
}
