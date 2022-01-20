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

 BKMNS2PF

 USAGE
   ./bin/bk2pf PARAMETER MIN_DEPTH MAX_DEPTH LAT LON RESOLUTION

 EXAMPLE
   ./bin/bk2pf rho 10 2891 -10.0 -50.0 1.0

 COMMAND LINE ARGUMENTS
   PARAMETER              - model parameter to be expanded (vph, rho, eta, vsv, etc.)
   MIN_DEPTH              - minimum depth
   MAX_DEPTH              - maximum depth
   LAT LON                - latitude and longitude of the point
   RESOLUTION             - spatial distance (in km) between grid points along the radial
                            direction

 DESCRIPTION
   Reads the desired model parameter, minimum depth, maximum depth, the coordinates of the point
   in which you want the profile and the grid spacing along the radial direction, and creates a
   one-dimensional profile of the model. In case you want the perturbations instead of the
   absolute values, just add a 'd' at beginning of the parameter code (e.g., dvs, drho, etc).
   The routine writes the output to a file called PARAMETER_PF.dat.

----------------------------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>
#include <math.h>
#include "exmath.h"
#include "legendre.h"
#include "coordinates.h"
#include "structs.h"
#include "spharm.h"
#include "io.h"
#include "constants.h"

static double gMeanModel (double r, unsigned nl, struct MeanModel Mm[nl])
{
  /* Computes the geometric mean model value using
     linear interpolation */
  int n = (int) nl;
  int i = (int) (nl * (r - Mm[0].r) / (Mm[n - 1].r - Mm[0].r) + 0.5);

  if (i <= 0) return Mm[0].vgm; else if (i > n - 1) return Mm[n - 1].vgm;

  double m = (Mm[i].vgm - Mm[i - 1].vgm) / (Mm[i].r - Mm[i - 1].r);

  return Mm[i].vgm + m * (r - Mm[i].r);
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
                         double r1, double r2, unsigned nr,
                         double R[nr], double t, double p,
                         unsigned np_b, unsigned nt_b, unsigned nr_b,
                         double Bm[np_b][nt_b][nr_b], double M[nr],
                         unsigned nl,
                         struct MeanModel Mm[nl], bool dvv)
{
  /* Expands crust from block model */
  for (unsigned i = 0; i < nr; i++)
  {
    double r = R[i];

    if (r >= r1 && r <= r2)
    {
      double v = block2Value (sb, r1, r2, r, t, p,
                              np_b, nt_b, nr_b, Bm);

      M[i] = dvv && v ? 100 * log (v / gMeanModel (r, nl, Mm)) : v;
    }
  }
}

static void helpMenu (void)
{
  char *help_menu = "\n BK2PF"

                    "\n\n USAGE"
                    "\n    ./bin/bk2pf PARAMETER MIN_DEPTH MAX_DEPTH LAT LON RESOLUTION"

                    "\n\n EXAMPLE"
                    "\n    ./bin/bk2pf rho 10 2891 -10.0 -50.0 1.0"

                    "\n\n COMMAND LINE ARGUMENTS"
                    "\n    PARAMETER              - model parameter to be expanded (vph, rho, eta, vsv, etc.)"
                    "\n    MIN_DEPTH              - minimum depth"
                    "\n    MAX_DEPTH              - maximum depth"
                    "\n    LAT LON                - latitude and longitude of the point"
                    "\n    RESOLUTION             - spatial distance (in km) between grid points along the radial"
                    "\n                             direction"

                    "\n\n DESCRIPTION"
                    "\n    Reads the desired model parameter, minimum depth, maximum depth, the coordinates of the point"
                    "\n    in which you want the profile and the grid spacing along the radial direction, and creates a"
                    "\n    one-dimensional profile of the model. In case you want the perturbations instead of the"
                    "\n    absolute values, just add a 'd' at beginning of the parameter code (e.g., dvs, drho, etc)."
                    "\n    The routine writes the output to a file called PARAMETER_PF.dat.\n\n";

  fprintf (stderr, "%s", help_menu);
}

int main (int argc, char *argv[])
{
  if (argc != 7)
  {
    fprintf (stderr, "\n Error: wrong number of parameters on the comand line...\n");
    helpMenu ();

    return 1;
  }

  char *prm = argv[1];

  bool dvv = prm[0] == 'd' ? true : false;

  if (dvv) prm = &prm[1];

  double dmin = atof (argv[2]);
  double dmax = atof (argv[3]);
  double t    = degree2Rad (90 - atof (argv[4]));
  double p    = degree2Rad (atof (argv[5]));
  double dr   = atof (argv[6]);

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

  fprintf (stderr, "\nReading mean model...\n");

  unsigned nl = 0;

  if (checkMeanModelHeaderIO (readMeanModelHeader (prm, &nl))) return 1;

  struct MeanModel Mm[nl];

  if (checkMeanModelIO (readMeanModel (prm, nl, Mm))) return 1;

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

  double R[nr];

  createProfilePf (r1, r2, nr, R);

  double (*M) = malloc (sizeof (double[nr]));

  initializeArray (nr, M);

  fprintf (stderr, "Zone 1 stretching from %.1lf to %.1lf km depth...\n",
           r2Depth (rmax), r2Depth (rmin));
  fprintf (stderr, "Expanding zone 1...");

  expandValue (&sb, rmin, rmax, nr, R, t, p,
               np_b, nt_b, nr_b, Bm1, M, nl, Mm, dvv);

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
  fprintf (stderr, "Expanding zone 2...");

  expandValue (&sb, rmin, rmax, nr, R, t, p,
               np_b, nt_b, nr_b, Bm2, M, nl, Mm, dvv);

  free (Bm2);

  fprintf (stderr, "\n\nWriting expansion...\n");

  if (checkExpansionIO (writeExpansionPf (argv, dvv,
                                          nr, R, M))) return 1;

  free (M);

  fprintf (stderr, "Done!\n");

  return 0;
}

