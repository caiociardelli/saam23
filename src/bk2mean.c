/*
 SphGLLTools/SphModel

 Author: Caio Ciardelli, University of São Paulo, October 2020

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along
 with this program; if not, write to the Free Software Foundation, Inc.,
 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

-----------------------------------------------------------------------------------------------

 BK2MEAN

 USAGE
   ./bin/bk2mean PARAMETER MIN_DEPTH MAX_DEPTH RADIAL_RESOLUTION HORIZONTAL_RESOLUTION

 EXAMPLE
   ./bin/bin/bk2mean vph -5 2891 0.5 1.0 10

 COMMAND LINE ARGUMENTS
   PARAMETER              - model parameter to be expanded (vph, rho, eta, vsv, etc.)
   MIN_DEPTH              - minimum depth
   MAX_DEPTH              - maximum depth
   HORIZONTAL_RESOLUTION  - spatial distance (in degrees) between grid points along the
                            horizontal direction
   RADIAL_RESOLUTION      - spatial distance (in km) between grid points along the radial
                            direction

 DESCRIPTION
   Reads the desired model parameter, the minimum and maximum depths, the horizontal resolution
   used to compute the arithmetic and geometric mean, and the radial resolution from the command
   line and creates a vertical cross-section of the model. The routine writes the output to a
   file called PARAMETER.dat.

----------------------------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>
#include <time.h>
#include <math.h>
#include "exmath.h"
#include "legendre.h"
#include "coordinates.h"
#include "structs.h"
#include "spharm.h"
#include "io.h"
#include "progress.h"
#include "constants.h"

static void computeMeans (struct MeanModel *mo,
                          unsigned np, unsigned nt,
                          double M[np][nt])
{
  /* Computes arithmetic and geometric means */
  unsigned nv = 0;

  mo->vam = 0;
  mo->vgm = 0;

  for (unsigned i = 0; i < np; i++)

    for (unsigned j = 0; j < nt; j++)
    {
      double v = M[i][j];

      if (v > WATER_LEVEL)
      {
        mo->vam += v;
        mo->vgm += log (v);

        nv++;
      }
    }

  if (nv)
  {
    mo->vam /= nv;
    mo->vgm /= nv;
    mo->vgm = exp (mo->vgm);
  }
}

static double block2Value (struct SphericalBoundaries *sb,
                           double rmin, double rmax,
                           double radius, double theta, double phi,
                           unsigned np_b, unsigned nt_b, unsigned nr_b,
                           double Bm[np_b][nt_b][nr_b])
{
  /* Computes parameter value from block model */
  int np = (int) np_b;
  int nt = (int) nt_b;
  int nr = (int) nr_b;

  if (phi    < sb->pmin) return 0;
  if (phi    > sb->pmax) return 0;
  if (theta  < sb->tmin) return 0;
  if (theta  > sb->tmax) return 0;
  if (radius < rmin)     radius = rmin;
  if (radius > rmax)     radius = rmax;

  double p = phi + 1.9198621771937625;
  double t = theta - 1.1344640137963142;
  double r = radius - rmin;

  double dp = (sb->pmax - sb->pmin) / (np - 1);
  double dt = (sb->tmax - sb->tmin) / (nt - 1);
  double dr = (rmax - rmin) / (nr - 1);

  int i0 = (int) (p / dp);
  int j0 = (int) (t / dt);
  int k0 = (int) (r / dr);

  double p0 = i0 * dp;
  double t0 = j0 * dt;
  double r0 = k0 * dr;

  int i1 = (i0 < np - 1) ? i0 + 1 : i0;
  int j1 = (j0 < nt - 1) ? j0 + 1 : j0;
  int k1 = (k0 < nr - 1) ? k0 + 1 : k0;

  double p1 = p0 + dp;
  double t1 = t0 + dt;
  double r1 = r0 + dr;

  double f000 = Bm[i0][j0][k0];
  double f100 = Bm[i1][j0][k0];
  double f010 = Bm[i0][j1][k0];
  double f110 = Bm[i1][j1][k0];
  double f001 = Bm[i0][j0][k1];
  double f101 = Bm[i1][j0][k1];
  double f011 = Bm[i0][j1][k1];
  double f111 = Bm[i1][j1][k1];

  double sum = 0; unsigned n = 0;

  if (f000 >= WATER_LEVEL) {sum += f000; n++;}
  if (f100 >= WATER_LEVEL) {sum += f100; n++;}
  if (f010 >= WATER_LEVEL) {sum += f010; n++;}
  if (f110 >= WATER_LEVEL) {sum += f110; n++;}
  if (f001 >= WATER_LEVEL) {sum += f001; n++;}
  if (f101 >= WATER_LEVEL) {sum += f101; n++;}
  if (f011 >= WATER_LEVEL) {sum += f011; n++;}
  if (f111 >= WATER_LEVEL) {sum += f111; n++;}

  if (n == 0) return 0;

  if (n < 8)
  {
    double mean = sum / n;

    if (f000 < WATER_LEVEL) f000 = mean;
    if (f100 < WATER_LEVEL) f100 = mean;
    if (f010 < WATER_LEVEL) f010 = mean;
    if (f110 < WATER_LEVEL) f110 = mean;
    if (f001 < WATER_LEVEL) f001 = mean;
    if (f101 < WATER_LEVEL) f101 = mean;
    if (f011 < WATER_LEVEL) f011 = mean;
    if (f111 < WATER_LEVEL) f111 = mean;
  }

  double k000 = f000 * (p1 - p) * (t1 - t) * (r1 - r);
  double k100 = f100 * (p - p0) * (t1 - t) * (r1 - r);
  double k010 = f010 * (p1 - p) * (t - t0) * (r1 - r);
  double k110 = f110 * (p - p0) * (t - t0) * (r1 - r);
  double k001 = f001 * (p1 - p) * (t1 - t) * (r - r0);
  double k101 = f101 * (p - p0) * (t1 - t) * (r - r0);
  double k011 = f011 * (p1 - p) * (t - t0) * (r - r0);
  double k111 = f111 * (p - p0) * (t - t0) * (r - r0);

  return (k000 + k100 + k010 + k110 +
          k001 + k101 + k011 + k111) / (dp * dt * dr);
}

static void expandValue (struct SphericalBoundaries *sb,
                         double r1, double r2,
                         unsigned nr, unsigned nt, unsigned np,
                         double R[nr], double Theta[nt], double Phi[np],
                         unsigned np_b, unsigned nt_b, unsigned nr_b,
                         double Bm[np_b][nt_b][nr_b],
                         struct MeanModel Mo[nr])
{
  /* Expands crust or mantle from block model */
  double M[np][nt];

  clock_t starttime = clock ();

  for (unsigned i = 0; i < nr; i++)
  {
    double r = R[i];

    if (r >= r1 && r <= r2)
    {
      for (unsigned j = 0; j < nt; j++)

        for (unsigned k = 0; k < np; k++)

          M[k][j] = block2Value (sb, r1, r2, r, Theta[j], Phi[k],
                                 np_b, nt_b, nr_b, Bm);

      computeMeans (&Mo[i], np, nt, M);
    }

    progressBar (i, 10, nr, starttime);
  }
}

static void helpMenu (void)
{
  char *help_menu = "\n BK2MEAN"

                    "\n\n USAGE"
                    "\n    ./bin/bk2mean PARAMETER MIN_DEPTH MAX_DEPTH RADIAL_RESOLUTION HORIZONTAL_RESOLUTION"

                    "\n\n EXAMPLE"
                    "\n    ./bin/bk2mean vph -5 2891 0.5 1.0 10"

                    "\n\n COMMAND LINE ARGUMENTS"
                    "\n    PARAMETER              - model parameter to be expanded (vph, rho, eta, vsv, etc.)"
                    "\n    MIN_DEPTH              - minimum depth"
                    "\n    MAX_DEPTH              - maximum depth"
                    "\n    HORIZONTAL_RESOLUTION  - spatial distance (in degrees) between grid points along the"
                    "\n                             horizontal direction"
                    "\n    RADIAL_RESOLUTION      - spatial distance (in km) between grid points along the radial"
                    "\n                             direction"

                    "\n\n DESCRIPTION"
                    "\n    Reads the desired model parameter, the minimum and maximum depths, the horizontal resolution"
                    "\n    used to compute the arithmetic and geometric mean, and the radial resolution from the command"
                    "\n    line and creates a vertical cross-section of the model. The routine writes the output to a"
                    "\n    file called PARAMETER.dat.\n\n";

  fprintf (stderr, "%s", help_menu);
}

int main (int argc, char *argv[])
{
  if (argc != 6)
  {
    fprintf (stderr, "\n Error: wrong number of parameters on the comand line...\n");
    helpMenu ();

    return 1;
  }

  char *prm = argv[1];

  double dmin = atof (argv[2]);
  double dmax = atof (argv[3]);
  double dh   = atof (argv[4]);
  double dr   = atof (argv[5]);

  unsigned np = (unsigned) (360.0 / dh + 1);
  unsigned nt = (unsigned) (180.0 / dh + 1);
  unsigned nr = (unsigned) (fabs (dmax - dmin) / dr + 1);

  double r1 = depth2R (dmax);
  double r2 = depth2R (dmin);

  struct SphericalBoundaries sb;

  sb.rmin =  CMB_R;
  sb.rmax =  TOP_R;
  sb.tmin =  1.1344640137963142;
  sb.tmax =  2.7401669256310974;
  sb.pmin = -1.9198621771937625;
  sb.pmax =  0.2094395102393195;

  if (r1 < sb.rmin || r1 > sb.rmax || r2 < sb.rmin || r2 > sb.rmax)
  {
    fprintf (stderr, "\n Error: depths should be between %.1lf and %.1lf km....\n",
             r2Depth (sb.rmin), r2Depth (sb.rmax)); return 1;
  }

  char *dir = "crust";

  double rmin = MOHO_R;
  double rmax = TOP_R;

  unsigned np_b, nt_b, nr_b;

  if (checkBlockModelHeaderIO (readBlockModelHeader (dir, prm,
                                                     &np_b, &nt_b, &nr_b))) return 1;

  double (*Bm1)[nt_b][nr_b] = malloc (sizeof (double[np_b][nt_b][nr_b]));

  fprintf (stderr, "\nReading block model from %s.bin file...\n", prm);

  if (checkBlockModelIO (readBlockModel (dir, prm, np_b, nt_b, nr_b, Bm1))) return 1;

  fprintf (stderr, "\nCreating output grid...\n");

  double R[nr], Theta[nt], Phi[np];

  createProfilePf (r1, r2, nr, R);
  createOutputGridLL (nt, np, Theta, Phi);

  struct MeanModel (*Mo) = malloc (sizeof (struct MeanModel[nr]));

  initializeMeanArray (nr, Mo);

  fprintf (stderr, "Zone 1 stretching from %.1lf to %.1lf km depth...\n",
           r2Depth (rmax), r2Depth (rmin));
  fprintf (stderr, "Computing mean model for zone 1...\n\n");

  expandValue (&sb, rmin, rmax, nr, nt, np,
               R, Theta, Phi, np_b, nt_b, nr_b, Bm1, Mo);

  free (Bm1);

  dir = "mantle";

  rmin = CMB_R;
  rmax = MOHO_R;

  if (checkBlockModelHeaderIO (readBlockModelHeader (dir, prm,
                                                     &np_b, &nt_b, &nr_b))) return 1;

  double (*Bm2)[nt_b][nr_b] = malloc (sizeof (double[np_b][nt_b][nr_b]));

  fprintf (stderr, "\n\nReading block model from %s.bin file...\n", prm);

  if (checkBlockModelIO (readBlockModel (dir, prm, np_b, nt_b, nr_b, Bm2))) return 1;

  fprintf (stderr, "\nZone 2 stretching from %.1lf to %.1lf km depth...\n",
           r2Depth (rmax), r2Depth (rmin));
  fprintf (stderr, "Computing mean model for zone 2...\n\n");

  expandValue (&sb, rmin, rmax, nr, nt, np,
               R, Theta, Phi, np_b, nt_b, nr_b, Bm2, Mo);

  free (Bm2);

  fprintf (stderr, "\n\nWriting expansion...\n");

  if (checkExpansionIO (writeModel1D (prm, nr, R, Mo))) return 1;

  free (Mo);

  fprintf (stderr, "Done!\n");

  return 0;
}

