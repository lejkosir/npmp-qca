//////////////////////////////////////////////////////////
// qca_sim: headless single-run simulator for fault     //
// sweeps. Based on main_batch_sim.c.                   //
//////////////////////////////////////////////////////////
// Runs one EXHAUSTIVE_VERIFICATION bistable simulation //
// and writes, per clock-period segment, the input      //
// values and the output polarization during the hold   //
// phase of the output cell's clock zone.               //
//                                                      //
// Usage:                                               //
//   qca_sim -f design [-n samples] [-s seed]           //
//           [-r radius] [-o result.csv]                //
//           [-T trace.csv] [-k trace_stride]           //
//////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <glib.h>
#include "fileio.h"
#include "design.h"
#include "global_consts.h"
#include "simulation.h"
#include "simulation_data.h"
#include "bistable_simulation.h"

extern bistable_OP bistable_options ;

static void usage (void)
  {
  fprintf (stderr,
    "Usage: qca_sim -f design [-n samples] [-s seed] [-r radius_of_effect]\n"
    "               [-o result.csv] [-T trace.csv] [-k trace_stride]\n") ;
  exit (1) ;
  }

int main (int argc, char **argv)
  {
  char *pszFName = NULL, *pszOut = NULL, *pszTrace = NULL ;
  unsigned int seed = 1 ;
  int trace_stride = 10 ;
  int Nix, Nix1, idxSeg, j ;
  DESIGN *design = NULL ;
  simulation_data *sim_data = NULL ;
  FILE *fp = stdout ;
  int icInputs, icOutputs, icSegments, segment_length ;
  int *out_clock = NULL ;

  // Table 1 of the reference paper; the engine defaults match except samples
  bistable_options.number_of_samples = 128000 ;
  bistable_options.animate_simulation = FALSE ;

  for (Nix = 1 ; Nix < argc ; Nix++)
    {
    if (Nix + 1 >= argc) usage () ;
    if (!strcmp (argv[Nix], "-f")) pszFName = argv[++Nix] ;
    else
    if (!strcmp (argv[Nix], "-n")) bistable_options.number_of_samples = atoi (argv[++Nix]) ;
    else
    if (!strcmp (argv[Nix], "-s")) seed = (unsigned int)strtoul (argv[++Nix], NULL, 10) ;
    else
    if (!strcmp (argv[Nix], "-r")) bistable_options.radius_of_effect = g_ascii_strtod (argv[++Nix], NULL) ;
    else
    if (!strcmp (argv[Nix], "-o")) pszOut = argv[++Nix] ;
    else
    if (!strcmp (argv[Nix], "-T")) pszTrace = argv[++Nix] ;
    else
    if (!strcmp (argv[Nix], "-k"))
      {
      trace_stride = atoi (argv[++Nix]) ;
      if (trace_stride < 1) trace_stride = 1 ;
      }
    else
      usage () ;
    }
  if (NULL == pszFName) usage () ;

  g_type_init () ;
  srand (seed) ;

  if (!open_project_file (pszFName, &design) || NULL == design)
    {
    fprintf (stderr, "Failed to open the circuit file !\n") ;
    return 3 ;
    }

  icInputs = design->bus_layout->inputs->icUsed ;
  icOutputs = design->bus_layout->outputs->icUsed ;
  if (0 == icInputs || 0 == icOutputs)
    {
    fprintf (stderr, "Design has %d inputs and %d outputs !\n", icInputs, icOutputs) ;
    return 4 ;
    }

  out_clock = g_malloc0 (sizeof (int) * icOutputs) ;
  for (Nix = 0 ; Nix < icOutputs ; Nix++)
    out_clock[Nix] = exp_array_index_1d (design->bus_layout->outputs, BUS_LAYOUT_CELL, Nix).cell->cell_options.clock ;

  bistable_options_dump (&bistable_options, stderr) ;

  if (NULL == (sim_data = run_simulation (BISTABLE, EXHAUSTIVE_VERIFICATION, design, NULL)))
    {
    fprintf (stderr, "Simulation failed !\n") ;
    return 5 ;
    }

  // The exhaustive clock runs 2^(inputs+1) periods over the whole run; inputs
  // are constant within each period, so one period == one segment.
  icSegments = 1 << (icInputs + 1) ;
  segment_length = sim_data->number_samples / icSegments ;

  if (NULL != pszOut && NULL == (fp = fopen (pszOut, "w")))
    {
    fprintf (stderr, "Failed to open %s for writing !\n", pszOut) ;
    return 6 ;
    }

  fprintf (fp, "segment") ;
  for (Nix = 0 ; Nix < icInputs ; Nix++)
    fprintf (fp, ",in_%s", sim_data->trace[Nix].data_labels) ;
  for (Nix = 0 ; Nix < icOutputs ; Nix++)
    {
    char *lbl = sim_data->trace[icInputs + Nix].data_labels ;
    fprintf (fp, ",%s_mid,%s_mean,%s_absmin,%s_hold_samples", lbl, lbl, lbl, lbl) ;
    }
  fprintf (fp, "\n") ;

  for (idxSeg = 0 ; idxSeg < icSegments ; idxSeg++)
    {
    int idxBeg = idxSeg * segment_length, idxEnd = idxBeg + segment_length ;

    fprintf (fp, "%d", idxSeg) ;
    for (Nix = 0 ; Nix < icInputs ; Nix++)
      fprintf (fp, ",%d", sim_data->trace[Nix].data[idxBeg + segment_length / 2] > 0 ? 1 : 0) ;

    for (Nix = 0 ; Nix < icOutputs ; Nix++)
      {
      double *clk = sim_data->clock_data[out_clock[Nix]].data ;
      double *pol = sim_data->trace[icInputs + Nix].data ;
      double dSum = 0.0, dAbsMin = 2.0 ;
      int icHold = 0, idxFirst = -1, idxLast = -1 ;

      // Hold phase == the clock is clamped at clock_low
      for (j = idxBeg ; j < idxEnd ; j++)
        if (clk[j] <= bistable_options.clock_low * (1.0 + 1e-9))
          {
          if (idxFirst < 0) idxFirst = j ;
          idxLast = j ;
          dSum += pol[j] ;
          dAbsMin = MIN (dAbsMin, fabs (pol[j])) ;
          icHold++ ;
          }

      if (icHold > 0)
        fprintf (fp, ",%.6f,%.6f,%.6f,%d", pol[(idxFirst + idxLast) / 2], dSum / icHold, dAbsMin, icHold) ;
      else
        fprintf (fp, ",nan,nan,nan,0") ;
      }
    fprintf (fp, "\n") ;
    }

  if (stdout != fp) fclose (fp) ;

  if (NULL != pszTrace)
    {
    if (NULL == (fp = fopen (pszTrace, "w")))
      {
      fprintf (stderr, "Failed to open %s for writing !\n", pszTrace) ;
      return 6 ;
      }
    fprintf (fp, "sample") ;
    for (Nix = 0 ; Nix < sim_data->number_of_traces ; Nix++)
      fprintf (fp, ",%s", sim_data->trace[Nix].data_labels) ;
    for (Nix = 0 ; Nix < 4 ; Nix++)
      fprintf (fp, ",clock%d", Nix) ;
    fprintf (fp, "\n") ;
    for (j = 0 ; j < sim_data->number_samples ; j += trace_stride)
      {
      fprintf (fp, "%d", j) ;
      for (Nix1 = 0 ; Nix1 < sim_data->number_of_traces ; Nix1++)
        fprintf (fp, ",%.6f", sim_data->trace[Nix1].data[j]) ;
      for (Nix1 = 0 ; Nix1 < 4 ; Nix1++)
        fprintf (fp, ",%.6e", sim_data->clock_data[Nix1].data[j]) ;
      fprintf (fp, "\n") ;
      }
    fclose (fp) ;
    }

  sim_data = simulation_data_destroy (sim_data) ;
  design = design_destroy (design) ;
  g_free (out_clock) ;

  return 0 ;
  }
