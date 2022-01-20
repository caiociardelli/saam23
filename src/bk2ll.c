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

 BK2LL

 USAGE
   ./bin/bk2ll PARAMETER DEPTH RESOLUTION

 EXAMPLE
   ./bin/bk2ll vsv 110 0.1

 COMMAND LINE ARGUMENTS
   PARAMETER              - model parameter to be expanded (vph, rho, eta, vsv, etc.)
   DEPTH                  - depth in which the depth slice will be created
   RESOLUTION             - spatial distance (in degrees) between both latitude and longitude
                            grid points of the depth slice

 DESCRIPTION
   Reads the desired model parameter, the depth, and the horizontal resolution from the command
   line, and creates a vertical cross-section of the model. In case you want the perturbations
   instead of the absolute values, just add a 'd' at the beginning of the parameter code (e.g.,
   dvs, drho, etc). The routine writes the output to a file called PARAMETER_DEPTH_DS.dat.

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
                         double r1, double r2, double r,
                         unsigned np, unsigned nt,
                         double Theta[nt], double Phi[np],
                         unsigned np_b, unsigned nt_b, unsigned nr_b,
                         double Bm[np_b][nt_b][nr_b], double M[np][nt],
                         unsigned nl,
                         struct MeanModel Mm[nl], bool dvv)
{
  /* Expands crust or mantle from block model */
  if (r >= r1 && r <= r2)

    for (unsigned i = 0; i < nt; i++)

      for (unsigned j = 0; j < np; j++)
      {
        double v = block2Value (sb, r1, r2, r, Theta[i], Phi[j],
                                np_b, nt_b, nr_b, Bm);

        M[j][i] = dvv && v ? 100 * log (v / gMeanModel (r, nl, Mm)) : v;
      }
}

static void helpMenu (void)
{
  char *help_menu = "\n BK2LL"

                    "\n\n USAGE"
                    "\n    ./bin/bk2ll PARAMETER DEPTH RESOLUTION"

                    "\n\n EXAMPLE"
                    "\n    ./bin/bk2ll vsv 110 0.1"

                    "\n\n COMMAND LINE ARGUMENTS"
                    "\n    PARAMETER              - model parameter to be expanded (vph, rho, eta, vsv, etc.)"
                    "\n    DEPTH                  - depth in which the depth slice will be created"
                    "\n    RESOLUTION             - spatial distance (in degrees) between both latitude and longitude"
                    "\n                             grid points of the depth slice"

                    "\n\n DESCRIPTION"
                    "\n    Reads the desired model parameter, the depth, and the horizontal resolution from the command"
                    "\n    line, and creates a vertical cross-section of the model. In case you want the perturbations"
                    "\n    instead of the absolute values, just add a 'd' at the beginning of the parameter code (e.g.,"
                    "\n    dvs, drho, etc). The routine writes the output to a file called PARAMETER_DEPTH_DS.dat.\n\n";

  fprintf (stderr, "%s", help_menu);
}

int main (int argc, char *argv[])
{
  if (argc != 4)
  {
    fprintf (stderr, "\n Error: wrong number of parameters on the comand line...\n");
    helpMenu ();

    return 1;
  }

  char *prm = argv[1];

  bool dvv = prm[0] == 'd' ? true : false;

  if (dvv) prm = &prm[1];

  double r  = depth2R (atof (argv[2]));
  double dh = atof (argv[3]);

  unsigned np = (unsigned) (122.0 / dh + 1);
  unsigned nt = (unsigned) (92.0 / dh + 1);

  struct SphericalBoundaries sb;

  sb.rmin =  CMB_R;
  sb.rmax =  TOP_R;
  sb.tmin =  1.1344640137963142;
  sb.tmax =  2.7401669256310974;
  sb.pmin = -1.9198621771937625;
  sb.pmax =  0.2094395102393195;

  if (r < sb.rmin || r > sb.rmax)
  {
    fprintf (stderr, "\n Error: depth should be between %.1lf and %.1lf km....\n",
             r2Depth (sb.rmin), r2Depth (sb.rmax)); return 1;
  }

  fprintf (stderr, "\nReading mean model...\n");

  unsigned nl = 0;

  if (checkMeanModelHeaderIO (readMeanModelHeader (prm, &nl))) return 1;

  struct MeanModel Mm[nl];

  if (checkMeanModelIO (readMeanModel (prm, nl, Mm))) return 1;

  char *dir;
  unsigned zone;
  double rmin, rmax;

  if (r >= MOHO_R)
  {
    dir = "crust";
    zone = 1;

    rmin = MOHO_R; rmax = TOP_R;
  }

  else
  {
    dir = "mantle";
    zone = 2;

    rmin = CMB_R; rmax = MOHO_R;
  }

  unsigned np_b, nt_b, nr_b;

  if (checkBlockModelHeaderIO (readBlockModelHeader (dir, prm,
                                                     &np_b, &nt_b, &nr_b))) return 1;

  double (*Bm)[nt_b][nr_b] = malloc (sizeof (double[np_b][nt_b][nr_b]));

  fprintf (stderr, "\nReading block model from %s.bin file...\n", prm);

  if (checkBlockModelIO (readBlockModel (dir, prm, np_b, nt_b, nr_b, Bm))) return 1;

  fprintf (stderr, "\nCreating output grid...\n");

  double Theta[nt], Phi[np];

  createOutputGridLL (nt, np, Theta, Phi);

  double (*M)[nt] = malloc (sizeof (double[np][nt]));

  initialize2DArray (np, nt, M);

  fprintf (stderr, "Zone %u stretching from %.1lf to %.1lf km depth...\n",
           zone, r2Depth (rmax), r2Depth (rmin));
  fprintf (stderr, "Expanding zone %u...", zone);

  expandValue (&sb, rmin, rmax, r, np, nt, Theta, Phi, np_b, nt_b, nr_b, Bm, M, nl, Mm, dvv);

  free (Bm);

  fprintf (stderr, "\n\nWriting expansion...\n");

  if (checkExpansionIO (writeExpansionLL (prm, dvv, &sb, r, np, nt, M))) return 1;

  free (M);

  fprintf (stderr, "Done!\n");

  return 0;
}

