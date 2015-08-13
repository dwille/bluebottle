/*******************************************************************************
 ********************************* BLUEBOTTLE **********************************
 *******************************************************************************
 *
 *  Copyright 2012 - 2015 Adam Sierakowski, The Johns Hopkins University
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *  Please contact the Johns Hopkins University to use Bluebottle for
 *  commercial and/or for-profit applications.
 ******************************************************************************/

#include <time.h>

#include "bluebottle.h"
#include "domain.h"
#include "particle.h"

void seeder_read_input(void)
{

  int N;               // number of parts
  real loa;            // interaction length

  real a;              // particle radius
  real xs,ys,zs;       // particle positions; used for start of arrays
  real aFx, aFy, aFz;  // particle linear forcing
  real aLx, aLy, aLz;  // particle angular forcing
  real rho;            // density
  real E;              // youngs modulus
  real sigma;          // poisson ratio
  real e_dry;          // dry coefficient of restitution
  real l_rough;        // particle surface roughness
  int order;           // lamb truncation order
  real rs_r;           // cage ratio extents
  real spring_k;       // particle spring constant
  real spring_x;       // spring attachment locations
  real spring_y;
  real spring_z;
  real spring_l;       // spring length
  int trans;           // particle is allowed to translate
  int rot;             // particle is allowed to rotate

  int fret = 0;
  fret = fret;         // prevent compiler warning

  // open configuration file for reading
  char fname[FILE_NAME_SIZE] = "";
  sprintf(fname, "%s/input/part.config", ROOT_DIR);
  FILE *infile = fopen(fname, "r");
  if(infile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  // read particle list
  fret = fscanf(infile, "n %d\n", &N);

#ifdef DOUBLE
  fret = fscanf(infile, "(l/a) %lf\n", &loa);
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "r %lf\n", &a);
  fret = fscanf(infile, "(x, y, z) %lf %lf %lf\n", &xs, &ys, &zs);
  fret = fscanf(infile, "(aFx, aFy, aFz) %lf %lf %lf\n",
    &aFx, &aFy, &aFz);
  fret = fscanf(infile, "(aLx, aLy, aLz) %lf %lf %lf\n",
    &aLx, &aLy, &aLz);
  fret = fscanf(infile, "rho %lf\n", &rho);
  fret = fscanf(infile, "E %lf\n", &E);
  fret = fscanf(infile, "sigma %lf\n", &sigma);
  fret = fscanf(infile, "e_dry %lf\n", &e_dry);
  fret = fscanf(infile, "l_rough %lf\n", &l_rough);
#else
  fret = fscanf(infile, "loa %f\n", &loa);
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "r %f\n", &a);
  fret = fscanf(infile, "(x, y, z) %f %f %f\n", &xs, &ys, &zs);
  fret = fscanf(infile, "(aFx, aFy, aFz) %f %f %f\n",
    &aFx, &aFy, &aFz);
  fret = fscanf(infile, "(aLx, aLy, aLz) %f %f %f\n",
    &aLx, &aLy, &aLz);
  fret = fscanf(infile, "rho %f\n", &rho);
  fret = fscanf(infile, "E %f\n", &E);
  fret = fscanf(infile, "sigma %f\n", &sigma);
  fret = fscanf(infile, "e_dry %f\n", &e_dry);
  fret = fscanf(infile, "l_rough %f\n", &l_rough);
#endif
  fret = fscanf(infile, "order %d\n", &order);
#ifdef DOUBLE
  fret = fscanf(infile, "rs/r %lf\n", &rs_r);
  fret = fscanf(infile, "spring_k %lf\n", &spring_k);
  fret = fscanf(infile, "spring (x, y, z) %lf %lf %lf\n",
    &spring_x, &spring_y, &spring_z);
  fret = fscanf(infile, "spring_l %lf\n", &spring_l);
#else // single precision
  fret = fscanf(infile, "rs/r %f\n", &rs_r);
  fret = fscanf(infile, "spring_k %f\n", &spring_k);
  fret = fscanf(infile, "spring (x, y, z) %f %f %f\n",
    &spring_x, &spring_y, &spring_z);
  fret = fscanf(infile, "spring_l %f\n", &spring_l);
#endif
  fret = fscanf(infile, "translating %d\n", &trans);
  fret = fscanf(infile, "rotating %d\n", &rot);

  // check parameters

  if (N < 1) {
    printf("Error: N must be > 1\n");
    exit(EXIT_FAILURE);
  } else if (a < 0) {
    printf("Error: a must be > 0\n");
    exit(EXIT_FAILURE);
  } else if (rho < 0) {
    printf("Error: rho must be > 0\n");
    exit(EXIT_FAILURE);
  } else if (E < 0) {
    printf("Error: E must be > 0\n");
    exit(EXIT_FAILURE);
  } else if (sigma > 0.5 || sigma <= -1) {
    printf("Error: sigma must be between -1 < sigma <= 0.5\n");
    exit(EXIT_FAILURE);
  } else if (order < 0) {
    printf("Error: order must be > 0\n");
    exit(EXIT_FAILURE);
  } else if (trans > 1) {
    printf("Error: translating must be 0 or 1\n");
    exit(EXIT_FAILURE);
  } else if (rot > 1) {
    printf("Error: rotating must be 0 or 1\n");
    exit(EXIT_FAILURE);
  }

  // SEEDER

  printf("Requested Parameters:\n");
  printf("       N = %d\n", N);
#ifdef DOUBLE
  printf("       (l/a) = %lf\n", loa);
  printf("       r = %.32f\n", a);
  printf("       (x, y, z) = (%lf, %lf, %lf)\n", xs, ys, zs);
  printf("       (aFx, aFy, aFz) = (%lf, %lf, %lf)\n", aFx, aFy, aFz);
  printf("       (aLx, aLy, aLz) = (%lf, %lf, %lf)\n", aLx, aLy, aLz);
  printf("       rho = %lf\n", rho);
  printf("       E = %lf\n", E);
  printf("       sigma = %lf\n", sigma);
  printf("       e_dry = %lf\n", e_dry);
  printf("       l_rough = %lf\n", l_rough);
  printf("       order = %d\n", order);
#else
  printf("       (l/a) = %f\n", loa);
  printf("       r = %f\n", a);
  printf("       (x, y, z) = (%f, %f, %f)\n", xs, ys, zs);
  printf("       (aFx, aFy, aFz) = (%f, %f, %f)\n", aFx, aFy, aFz);
  printf("       (aLx, aLy, aLz) = (%f, %f, %f)\n", aLx, aLy, aLz);
  printf("       rho = %f\n", rho);
  printf("       E = %f\n", E);
  printf("       sigma = %f\n", sigma);
  printf("       e_dry = %f\n", e_dry);
  printf("       l_rough = %f\n", l_rough);
  printf("       order = %d\n", order);
#endif
#ifdef DOUBLE
  printf("       rs_r = %lf\n", rs_r);
  printf("       spring_k = %lf\n", spring_k);
  printf("       spring (x, y, z) = (%lf, %lf, %lf)\n", 
    spring_x, spring_y, spring_z);
  printf("       spring_l = %lf\n", spring_l);
#else
  printf("       rs_r = %f\n", rs_r);
  printf("       spring_k = %f\n", spring_k);
  printf("       spring (x, y, z) = (%f, %f, %f)\n", 
    spring_x, spring_y, spring_z);
  printf("       spring_l = %f\n", spring_l);
#endif
  printf("       translating = %d\n", trans);
  printf("       rotating = %d\n\n", rot);
  fflush(stdout);

  printf("Seed particles according to parameters specified in");
  printf(" parts.config? (y/N)\n");
  fflush(stdout);
  int c = getchar();

  if (c == 'Y' || c == 'y') {
    printf("Seed particles for which kind of array?\n");
    printf("\t(r)andom / (a)rray / (h)ex / (p)erturbed?\n");
    fflush(stdout);
    //TODO: ask adam about this
    int tmp = getchar();
    tmp = tmp;
    int type = getchar();

    int Nx = 0; int Ny = 0; int Nz = 0;   // number of particles in each dir
    real dx = 0.; real dy = 0.; real dz = 0.; // interparticle spacing
    // TODO: maybe input Nx, Ny, Nz for hex? hard to tell b/c of spacing...
    real xExtent = 0.; real yExtent = 0.; real zExtent = 0; // array extents
    real alpha = 0.;      // volume fraction
    real ratio = 0.;      // max percent of radius that particle can beperturbed
    int nperturb = 0;     // n times to perturb

    switch(type) {
      case 'r':   // random
        seeder(N, loa, a, aFx, aFy, aFz, aLx, aLy, aLz, rho, E, sigma, e_dry, 
          l_rough, order, rs_r, spring_k, spring_x, spring_y, spring_z, 
          spring_l, trans, rot);   
        break;

      case 'a':   // array
#ifdef DOUBLE
        printf("Starting position taken from part.config: (%lf, %lf, %lf)\n",
          xs, ys, zs);
#else        
        printf("Starting position taken from part.config: (%f, %f, %f)\n",
          xs, ys, zs);
#endif
        // particle numbers
        printf("Input number of particles in X direction:\n");
        fflush(stdout);
        fret = scanf("%d", &Nx);
        printf("Input number of particles in Y direction:\n");
        fflush(stdout);
        fret = scanf("%d", &Ny);
        printf("Input number of particles in Z direction:\n");
        fflush(stdout);
        fret = scanf("%d", &Nz);
        // particle spacing (between surfaces of particle)
        printf("Input interparticle spacing in X direction:\n");
        fflush(stdout);
        fret = scanf("%lf", &dx);
        printf("Input interparticle spacing in Y direction:\n");
        fflush(stdout);
        fret = scanf("%lf", &dy);
        printf("Input interparticle spacing in Z direction:\n");
        fflush(stdout);
        fret = scanf("%lf", &dz);

        seeder_array(Nx, Ny, Nz, dx, dy, dx, xs, ys, zs, 
          loa, a, aFx, aFy, aFz, aLx, aLy, aLz, rho, E, sigma, e_dry, l_rough, 
          order, rs_r, spring_k, spring_x, spring_y, spring_z, spring_l,trans, 
          rot);      
        break;

      case 'h':   // hexagonal packed
#ifdef DOUBLE
        printf("Starting position taken from part.config: (%lf, %lf, %lf)\n",
          xs, ys, zs);
#else        
        printf("Starting position taken from part.config: (%f, %f, %f)\n",
          xs, ys, zs);
#endif
        // extents of array
        printf("Please input X extent\n");
        fflush(stdout);
        fret = scanf("%lf", &xExtent);
        printf("Please input Y extent\n");
        fflush(stdout);
        fret = scanf("%lf", &yExtent);
        printf("Please input Z extent\n");
        fflush(stdout);
        fret = scanf("%lf", &zExtent);
        // desired volume fraction of array
        printf("Please input the volume fraction (-1 for HCP)\n");
        fflush(stdout);
        fret = scanf("%lf", &alpha);

        seeder_hex(xExtent, yExtent, zExtent, alpha, xs, ys, zs,
          loa, a, aFx, aFy, aFz, aLx, aLy, aLz, rho, E, sigma, e_dry, l_rough, 
          order, rs_r, spring_k, spring_x, spring_y, spring_z, spring_l, trans, 
          rot);
        break;

      case 'p':   // perturbed
        printf("Please input the particle number in x direction\n");
        fflush(stdout);
        fret = scanf("%d", &Nx);
        printf("Please input the particle number in y direction\n");
        fflush(stdout);
        fret = scanf("%d", &Ny);
        printf("Please input the particle number in z direction\n");
        fflush(stdout);
        fret = scanf("%d", &Nz);
        printf("Please input the perturbation magnitude");
        printf(" (between 0 and 1)\n");
        fflush(stdout);
        fret = scanf("%lf", &ratio);
        printf("Please input the number of perturbations");
        printf(" (should be larger than 100,000)\n");
        fflush(stdout);
        fret = scanf("%d", &nperturb);

        seeder_high_vol_random(Nx, Ny, Nz, ratio, nperturb, 
          loa, a, aFx, aFy, aFz, aLx, aLy, aLz, rho, E, sigma, e_dry, l_rough, 
          order, rs_r, spring_k, spring_x, spring_y, spring_z, spring_l, trans, 
          rot);    
        break;

      default:    // any other options
        printf("Option unrecognized.\n");
        fflush(stdout);
        exit(EXIT_FAILURE);
    }

  } else {
    printf("Please specify the desired parameters in parts.config\n\n");
    fflush(stdout);
    exit(EXIT_FAILURE);
  }
}

void seeder(int nparts, real loa, real a, real aFx, real aFy, real aFz, 
  real aLx, real aLy, real aLz, real rho, real E, real sigma, real e_dry,
  real l_rough, int o, real rs, real spring_k, real spring_x, real spring_y,
  real spring_z, real spring_l, int t, int r)
{

  printf("Running random seeder for %d particles...\n\n", nparts);
  fflush(stdout);
  real xx, yy, zz;
  int fits = 1;
  int attempts = 1;
  int fail = 0;
  int redo = 1;

  // seed the random number generator
  srand(time(NULL));

  // read domain input
  domain_read_input();
  domain_init();

  // domain size accounting for screen
  real xs = Dom.xs + bc.dsW;
  real xe = Dom.xe - bc.dsE;
  real xl = Dom.xl - bc.dsE - bc.dsW;
  real ys = Dom.ys + bc.dsS;
  real ye = Dom.ye - bc.dsN;
  real yl = Dom.yl - bc.dsN - bc.dsS;
  real zs = Dom.zs + bc.dsB;
  real ze = Dom.ze - bc.dsT;
  real zl = Dom.zl - bc.dsT - bc.dsB;
  
  // allocate particle list
  parts = (part_struct*) malloc(nparts * sizeof(part_struct));
  cpumem += nparts * sizeof(part_struct);

  real gap = 1.00;

  // place the first particle
  parts[0].r = a;
  redo = 1;
  while(redo == 1) {
    redo = 0;
    parts[0].x = rand() / (real)RAND_MAX;
    parts[0].x *= xl;
    parts[0].x += xs;
    if((bc.uW != PERIODIC) && (parts[0].x < (xs + gap*parts[0].r)))
      redo = 1;
    if((bc.uE != PERIODIC) && (parts[0].x > (xe - gap*parts[0].r)))
      redo = 1;
  }
  redo = 1;
  while(redo == 1) {
    redo = 0;
    parts[0].y = rand() / (real)RAND_MAX;
    parts[0].y *= yl;
    parts[0].y += ys;
    if((bc.vS != PERIODIC) && (parts[0].y < (ys + gap*parts[0].r)))
      redo = 1;
    if((bc.vN != PERIODIC) && (parts[0].y > (ye - gap*parts[0].r)))
      redo = 1;
  }
  redo = 1;
  while(redo == 1) {
    redo = 0;
    parts[0].z = rand() / (real)RAND_MAX;
    //parts[0].z = acos(2.*(parts[0].z-0.5))/PI;
    parts[0].z *= zl;
    parts[0].z += zs;
    if((bc.wB != PERIODIC) && (parts[0].z < (zs + gap*parts[0].r)))
      redo = 1;
    if((bc.wT != PERIODIC) && (parts[0].z > (ze - gap*parts[0].r)))
      redo = 1;
  }

  parts[0].u = 0;
  parts[0].v = 0;
  parts[0].w = 0;
  parts[0].aFx = aFx;
  parts[0].aFy = aFy;
  parts[0].aFz = aFz;
  parts[0].aLx = aLx;
  parts[0].aLy = aLy;
  parts[0].aLz = aLz;
  parts[0].rho = rho;
  parts[0].E = E;
  parts[0].sigma = sigma;
  parts[0].e_dry = e_dry;
  parts[0].l_rough = l_rough;
  parts[0].order = o;
  parts[0].rs = rs;
  parts[0].spring_k = spring_k;
  parts[0].spring_x = spring_x;
  parts[0].spring_y = spring_y;
  parts[0].spring_z = spring_z;
  parts[0].spring_l = spring_l;
  parts[0].ncoeff = 0;
  parts[0].translating = t;
  parts[0].rotating = r;

  // place the rest of the particles
  int i = 0;
  for(i = 1; i < nparts; i++) {
    fits = !fits;
    if(fail) break;
    while(!fits) {
      attempts++;
      // place particle
      parts[i].r = a;
      redo = 1;
      while(redo == 1) {
        redo = 0;
        parts[i].x = rand() / (real)RAND_MAX;
        parts[i].x *= xl;
        parts[i].x += xs;
        if((bc.uW != PERIODIC) && (parts[i].x < (xs + gap*parts[i].r)))
          redo = 1;
        if((bc.uE != PERIODIC) && (parts[i].x > (xe - gap*parts[i].r)))
          redo = 1;
      }
      redo = 1;
      while(redo == 1) {
        redo = 0;
        parts[i].y = rand() / (real)RAND_MAX;
        parts[i].y *= yl;
        parts[i].y += ys;
        if((bc.vS != PERIODIC) && (parts[i].y < (ys + gap*parts[i].r)))
          redo = 1;
        if((bc.vN != PERIODIC) && (parts[i].y > (ye - gap*parts[i].r)))
          redo = 1;
      }
      redo = 1;
      while(redo == 1) {
        redo = 0;
        parts[i].z = rand() / (real)RAND_MAX;
        parts[i].z *= zl;
        parts[i].z += zs;
        if((bc.wB != PERIODIC) && (parts[i].z < (zs + gap*parts[i].r)))
          redo = 1;
        if((bc.wT != PERIODIC) && (parts[i].z > (ze - gap*parts[i].r)))
          redo = 1;
      }

      parts[i].u = 0;
      parts[i].v = 0;
      parts[i].w = 0;
      parts[i].aFx = aFx;
      parts[i].aFy = aFy;
      parts[i].aFz = aFz;
      parts[i].aLx = aLx;
      parts[i].aLy = aLy;
      parts[i].aLz = aLz;
      parts[i].rho = rho;
      parts[i].E = E;
      parts[i].sigma = sigma;
      parts[i].e_dry = e_dry;
      parts[i].l_rough = l_rough;
      parts[i].order = o;
      parts[i].rs = rs;
      parts[i].spring_k = spring_k;
      parts[i].spring_x = spring_x;
      parts[i].spring_y = spring_y;
      parts[i].spring_z = spring_z;
      parts[i].spring_l = spring_l;
      parts[i].ncoeff = 0;
      parts[i].translating = t;
      parts[i].rotating = r;

      // check that this particle does not intersect any other particle
      fits = !fits;
      for(int j = 0; j < i; j++) {
        xx = parts[i].x - parts[j].x;
        xx = xx * xx;
        yy = parts[i].y - parts[j].y;
        yy = yy * yy;
        zz = parts[i].z - parts[j].z;
        zz = zz * zz;
        if(sqrt(xx + yy + zz) < (gap*parts[i].r + gap*parts[j].r)) {
          fits = !fits;
          break;
        }

        // also use virtual particle to check if a particle is too close in
        // a periodic direction
        // x only
        if(bc.uW == PERIODIC) {
          xx = parts[i].x - parts[j].x;
          xx = xx * xx;
          yy = parts[i].y - parts[j].y;
          yy = yy * yy;
          zz = parts[i].z - parts[j].z;
          zz = zz * zz;
          if(parts[i].x < (xs + parts[i].r))
            xx = parts[i].x + xl - parts[j].x;
          if(parts[i].x > (xe - parts[i].r))
            xx = parts[i].x - xl - parts[j].x;
          xx = xx * xx;
          yy = parts[i].y - parts[j].y;
          yy = yy * yy;
          zz = parts[i].z - parts[j].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (gap*parts[i].r + gap*parts[j].r)) {
            fits = !fits;
            break;
          }
        }

        // y only
        if(bc.vS == PERIODIC) {
          xx = parts[i].x - parts[j].x;
          xx = xx * xx;
          yy = parts[i].y - parts[j].y;
          yy = yy * yy;
          zz = parts[i].z - parts[j].z;
          zz = zz * zz;
          xx = parts[i].x - parts[j].x;
          xx = xx * xx;
          if(parts[i].y < (ys + parts[i].r))
            yy = parts[i].y + yl - parts[j].y;
          if(parts[i].y > (ye - parts[i].r))
            yy = parts[i].y - yl - parts[j].y;
          yy = yy * yy;
          zz = parts[i].z - parts[j].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (gap*parts[i].r + gap*parts[j].r)) {
            fits = !fits;
            break;
          }
        }

        // z only
        if(bc.wB == PERIODIC) {
          xx = parts[i].x - parts[j].x;
          xx = xx * xx;
          yy = parts[i].y - parts[j].y;
          yy = yy * yy;
          zz = parts[i].z - parts[j].z;
          zz = zz * zz;
          xx = parts[i].x - parts[j].x;
          xx = xx * xx;
          yy = parts[i].y - parts[j].y;
          yy = yy * yy;
          if(parts[i].z < (zs + parts[i].r))
            zz = parts[i].z + zl - parts[j].z;
          if(parts[i].z > (ze - parts[i].r))
            zz = parts[i].z - zl - parts[j].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (gap*parts[i].r + gap*parts[j].r)) {
            fits = !fits;
            break;
          }
        }

        // x and y
        if(bc.uW == PERIODIC && bc.vS == PERIODIC) {
          xx = parts[i].x - parts[j].x;
          xx = xx * xx;
          yy = parts[i].y - parts[j].y;
          yy = yy * yy;
          zz = parts[i].z - parts[j].z;
          zz = zz * zz;
          if(parts[i].x < (xs + parts[i].r))
            xx = parts[i].x + xl - parts[j].x;
          if(parts[i].x > (xe - parts[i].r))
            xx = parts[i].x - xl - parts[j].x;
          xx = xx * xx;
          if(parts[i].y < (ys + parts[i].r))
            yy = parts[i].y + yl - parts[j].y;
          if(parts[i].y > (ye - parts[i].r))
            yy = parts[i].y - yl - parts[j].y;
          yy = yy * yy;
          zz = parts[i].z - parts[j].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (gap*parts[i].r + gap*parts[j].r)) {
            fits = !fits;
            break;
          }
        }

        // y and z
        if(bc.vS == PERIODIC && bc.wB == PERIODIC) {
          xx = parts[i].x - parts[j].x;
          xx = xx * xx;
          yy = parts[i].y - parts[j].y;
          yy = yy * yy;
          zz = parts[i].z - parts[j].z;
          zz = zz * zz;
          xx = parts[i].x - parts[j].x;
          xx = xx * xx;
          if(parts[i].y < (ys + parts[i].r))
            yy = parts[i].y + yl - parts[j].y;
          if(parts[i].y > (ye - parts[i].r))
            yy = parts[i].y - yl - parts[j].y;
          yy = yy * yy;
          if(parts[i].z < (zs + parts[i].r))
            zz = parts[i].z + zl - parts[j].z;
          if(parts[i].z > (ze - parts[i].r))
            zz = parts[i].z - zl - parts[j].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (gap*parts[i].r + gap*parts[j].r)) {
            fits = !fits;
            break;
          }
        }

        // z and x
        if(bc.uW == PERIODIC && bc.wB == PERIODIC) {
          xx = parts[i].x - parts[j].x;
          xx = xx * xx;
          yy = parts[i].y - parts[j].y;
          yy = yy * yy;
          zz = parts[i].z - parts[j].z;
          zz = zz * zz;
          if(parts[i].x < (xs + parts[i].r))
            xx = parts[i].x + xl - parts[j].x;
          if(parts[i].x > (xe - parts[i].r))
            xx = parts[i].x - xl - parts[j].x;
          xx = xx * xx;
          yy = parts[i].y - parts[j].y;
          yy = yy * yy;
          if(parts[i].z < (zs + parts[i].r))
            zz = parts[i].z + zl - parts[j].z;
          if(parts[i].z > (ze - parts[i].r))
            zz = parts[i].z - zl - parts[j].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (gap*parts[i].r + gap*parts[j].r)) {
            fits = !fits;
            break;
          }
        }

        // x, y, and z
        if(bc.uW == PERIODIC && bc.vS == PERIODIC && bc.wB == PERIODIC) {
          xx = parts[i].x - parts[j].x;
          xx = xx * xx;
          yy = parts[i].y - parts[j].y;
          yy = yy * yy;
          zz = parts[i].z - parts[j].z;
          zz = zz * zz;
          if(parts[i].x < (xs + parts[i].r))
            xx = parts[i].x + xl - parts[j].x;
          if(parts[i].x > (xe - parts[i].r))
            xx = parts[i].x - xl - parts[j].x;
          xx = xx * xx;
          if(parts[i].y < (ys + parts[i].r))
            yy = parts[i].y + yl - parts[j].y;
          if(parts[i].y > (ye - parts[i].r))
            yy = parts[i].y - yl - parts[j].y;
          yy = yy * yy;
          if(parts[i].z < (zs + parts[i].r))
            zz = parts[i].z + zl - parts[j].z;
          if(parts[i].z > (ze - parts[i].r))
            zz = parts[i].z - zl - parts[j].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (gap*parts[i].r + gap*parts[j].r)) {
            fits = !fits;
            break;
          }
        }

        // check both ways
        // x only
        if(bc.uW == PERIODIC) {
          xx = parts[j].x - parts[i].x;
          xx = xx * xx;
          yy = parts[j].y - parts[i].y;
          yy = yy * yy;
          zz = parts[j].z - parts[i].z;
          zz = zz * zz;
          if(parts[j].x < (xs + parts[j].r))
            xx = parts[j].x + xl - parts[i].x;
          if(parts[j].x > (xe - parts[j].r))
            xx = parts[j].x - xl - parts[i].x;
          xx = xx * xx;
          yy = parts[j].y - parts[i].y;
          yy = yy * yy;
          zz = parts[j].z - parts[i].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (gap*parts[j].r + gap*parts[i].r)) {
            fits = !fits;
            break;
          }
        }

        // y only
        if(bc.vS == PERIODIC) {
          xx = parts[j].x - parts[i].x;
          xx = xx * xx;
          yy = parts[j].y - parts[i].y;
          yy = yy * yy;
          zz = parts[j].z - parts[i].z;
          zz = zz * zz;
          xx = parts[j].x - parts[i].x;
          xx = xx * xx;
          if(parts[j].y < (ys + parts[j].r))
            yy = parts[j].y + yl - parts[i].y;
          if(parts[j].y > (ye - parts[j].r))
            yy = parts[j].y - yl - parts[i].y;
          yy = yy * yy;
          zz = parts[j].z - parts[i].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (gap*parts[j].r + gap*parts[i].r)) {
            fits = !fits;
            break;
          }
        }

        // z only
        if(bc.wB == PERIODIC) {
          xx = parts[j].x - parts[i].x;
          xx = xx * xx;
          yy = parts[j].y - parts[i].y;
          yy = yy * yy;
          zz = parts[j].z - parts[i].z;
          zz = zz * zz;
          xx = parts[j].x - parts[i].x;
          xx = xx * xx;
          yy = parts[j].y - parts[i].y;
          yy = yy * yy;
          if(parts[j].z < (zs + parts[j].r))
            zz = parts[j].z + zl - parts[i].z;
          if(parts[j].z > (ze - parts[j].r))
            zz = parts[j].z - zl - parts[i].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (gap*parts[j].r + gap*parts[i].r)) {
            fits = !fits;
            break;
          }
        }

        // x and y
        if(bc.uW == PERIODIC && bc.vS == PERIODIC) {
          xx = parts[j].x - parts[i].x;
          xx = xx * xx;
          yy = parts[j].y - parts[i].y;
          yy = yy * yy;
          zz = parts[j].z - parts[i].z;
          zz = zz * zz;
          if(parts[j].x < (xs + parts[j].r))
            xx = parts[j].x + xl - parts[i].x;
          if(parts[j].x > (xe - parts[j].r))
            xx = parts[j].x - xl - parts[i].x;
          xx = xx * xx;
          if(parts[j].y < (ys + parts[j].r))
            yy = parts[j].y + yl - parts[i].y;
          if(parts[j].y > (ye - parts[j].r))
            yy = parts[j].y - yl - parts[i].y;
          yy = yy * yy;
          zz = parts[j].z - parts[i].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (gap*parts[j].r + gap*parts[i].r)) {
            fits = !fits;
            break;
          }
        }

        // y and z
        if(bc.vS == PERIODIC && bc.wB == PERIODIC) {
          xx = parts[j].x - parts[i].x;
          xx = xx * xx;
          yy = parts[j].y - parts[i].y;
          yy = yy * yy;
          zz = parts[j].z - parts[i].z;
          zz = zz * zz;
          xx = parts[j].x - parts[i].x;
          xx = xx * xx;
          if(parts[j].y < (ys + parts[j].r))
            yy = parts[j].y + yl - parts[i].y;
          if(parts[j].y > (ye - parts[j].r))
            yy = parts[j].y - yl - parts[i].y;
          yy = yy * yy;
          if(parts[j].z < (zs + parts[j].r))
            zz = parts[j].z + zl - parts[i].z;
          if(parts[j].z > (ze - parts[j].r))
            zz = parts[j].z - zl - parts[i].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (gap*parts[j].r + gap*parts[i].r)) {
            fits = !fits;
            break;
          }
        }

        // z and x
        if(bc.uW == PERIODIC && bc.wB == PERIODIC) {
          xx = parts[j].x - parts[i].x;
          xx = xx * xx;
          yy = parts[j].y - parts[i].y;
          yy = yy * yy;
          zz = parts[j].z - parts[i].z;
          zz = zz * zz;
          if(parts[j].x < (xs + parts[j].r))
            xx = parts[j].x + xl - parts[i].x;
          if(parts[j].x > (xe - parts[j].r))
            xx = parts[j].x - xl - parts[i].x;
          xx = xx * xx;
          yy = parts[j].y - parts[i].y;
          yy = yy * yy;
          if(parts[j].z < (zs + parts[j].r))
            zz = parts[j].z + zl - parts[i].z;
          if(parts[j].z > (ze - parts[j].r))
            zz = parts[j].z - zl - parts[i].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (gap*parts[j].r + gap*parts[i].r)) {
            fits = !fits;
            break;
          }
        }

        // x, y, and z
        if(bc.uW == PERIODIC && bc.vS == PERIODIC && bc.wB == PERIODIC) {
          xx = parts[j].x - parts[i].x;
          xx = xx * xx;
          yy = parts[j].y - parts[i].y;
          yy = yy * yy;
          zz = parts[j].z - parts[i].z;
          zz = zz * zz;
          if(parts[j].x < (xs + parts[j].r))
            xx = parts[j].x + xl - parts[i].x;
          if(parts[j].x > (xe - parts[j].r))
            xx = parts[j].x - xl - parts[i].x;
          xx = xx * xx;
          if(parts[j].y < (ys + parts[j].r))
            yy = parts[j].y + yl - parts[i].y;
          if(parts[j].y > (ye - parts[j].r))
            yy = parts[j].y - yl - parts[i].y;
          yy = yy * yy;
          if(parts[j].z < (zs + parts[j].r))
            zz = parts[j].z + zl - parts[i].z;
          if(parts[j].z > (ze - parts[j].r))
            zz = parts[j].z - zl - parts[i].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (gap*parts[j].r + gap*parts[i].r)) {
            fits = !fits;
            break;
          }
        }

      }
      if(attempts == 1e5*nparts) {
        fail = !fail;
        break;
      }
    }
  }

  if(fail) {
    printf("After %d attempts, the seeder has placed", attempts);
    printf(" %d of %d particles (a = %f).\n\n", i-1, nparts, a);
    printf("...bluebottle seeder done.\n\n");
    exit(EXIT_FAILURE);
  }

  printf("It took %d attempts to place %d", attempts, nparts);
  printf(" particles (a = %f) with no intersections.\n\n", a);
  fflush(stdout);

  printf("Writing part_seeder.config...");
  fflush(stdout);
  // write particle configuration to file
  char fname[FILE_NAME_SIZE] = "";
  // open file for writing
  sprintf(fname, "%spart_seeder.config", INPUT_DIR);
  FILE *ofile = fopen(fname, "w");
  if(ofile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  // write the number of particles and compact support length
  fprintf(ofile, "n %d\n", nparts);
  fprintf(ofile, "(l/a) %f\n", loa);

  // write each particle configuration
  for(int i = 0; i < nparts; i++) {
    fprintf(ofile, "\n");
    fprintf(ofile, "r %f\n", parts[i].r);
    fprintf(ofile, "(x, y, z) %f %f %f\n", parts[i].x, parts[i].y, parts[i].z);
    fprintf(ofile, "(aFx, aFy, aFz) %f %f %f\n", parts[i].aFx, parts[i].aFy,
      parts[i].aFz);
    fprintf(ofile, "(aLx, aLy, aLz) %f %f %f\n", parts[i].aLx, parts[i].aLy,
      parts[i].aLz);
    fprintf(ofile, "rho %f\n", parts[i].rho);
    fprintf(ofile, "E %f\n", parts[i].E);
    fprintf(ofile, "sigma %f\n", parts[i].sigma);
    fprintf(ofile, "e_dry %f\n", parts[i].e_dry);
    fprintf(ofile, "l_rough %f\n", parts[i].l_rough);
    fprintf(ofile, "order %d\n", parts[i].order);
    fprintf(ofile, "rs/r %f\n", parts[i].rs);
    fprintf(ofile, "spring_k %f\n", parts[i].spring_k);
    fprintf(ofile, "spring (x, y, z) %f %f %f\n",
      parts[i].spring_x, parts[i].spring_y, parts[i].spring_z);
    fprintf(ofile, "spring_l %f\n", parts[i].spring_l);
    fprintf(ofile, "translating %d\n", parts[i].translating);
    fprintf(ofile, "rotating %d\n", parts[i].rotating);
  }

  // close the file
  fclose(ofile);
  printf("done.\n");
  printf("\n...bluebottle seeder done.\n\n");
  fflush(stdout);

  // clean up
  domain_clean();
  parts_clean();
}

void seeder_array(int Nx, int Ny, int Nz, real dx, real dy, real dz,
  real xs, real ys, real zs, real loa, real a, real aFx, real aFy, real aFz, 
  real aLx, real aLy, real aLz, real rho, real E, real sigma, real e_dry, 
  real l_rough, int o, real rs, real spring_k, real spring_x, real spring_y, 
  real spring_z, real spring_l, int t, int r)
{

  int nparts = Nx*Ny*Nz;
  printf("Running array seeder for %d particles...\n\n", nparts);
  fflush(stdout);

  // read domain input
  domain_read_input();
  domain_init();

  // domain size accounting for screen
  real dom_xs = Dom.xs + bc.dsW;
  real dom_xl = Dom.xl - bc.dsE - bc.dsW;
  real dom_ys = Dom.ys + bc.dsS;
  real dom_yl = Dom.yl - bc.dsN - bc.dsS;
  real dom_zs = Dom.zs + bc.dsB;
  real dom_zl = Dom.zl - bc.dsT - bc.dsB;

  // allocate particle list
  parts = (part_struct*) malloc(nparts * sizeof(part_struct));
  cpumem += nparts * sizeof(part_struct);

  // check for interactions
  if (dx < 0 || dy < 0 || dz < 0){
    // check particle spacing
    printf("Particles intersect! (dx, dy, dz) must be positive\n");
    fflush(stdout);
    printf("...bluebottle seeder failed.\n\n");
    exit(EXIT_FAILURE);
  } else if ((xs - a) < dom_xs || (ys - a) < dom_ys || (zs - a) < dom_zs) {
    // check domain boundaries
    printf("Starting position creates particle outside of domain\n");
    fflush(stdout);
    printf("...bluebottle seeder failed.\n\n");
    exit(EXIT_FAILURE);
  } else if ((Nx*a + (Nx - 1)*dx) > dom_xl ||
             (Ny*a + (Ny - 1)*dy) > dom_yl ||
             (Nz*a + (Nz - 1)*dz) > dom_zl) {
    // check if particles fit in domain
    printf("Too many particles in the domain for given numbers and spacing\n");
    fflush(stdout);
    printf("...bluebottle seeder failed.\n\n");
    exit(EXIT_FAILURE);
  }

  // distance between particle centers
  real Rx = dx + 2*a;
  real Ry = dy + 2*a;
  real Rz = dz + 2*a;

  for (int k = 0; k < Nz; k++) {
    for (int j = 0; j < Ny; j++) {
      for (int i = 0; i < Nx; i++) {
        parts[i + j*Nx + k*(Nx*Ny)].x = xs + i*Rx;
        parts[i + j*Nx + k*(Nx*Ny)].y = ys + j*Ry;
        parts[i + j*Nx + k*(Nx*Ny)].z = zs + k*Rz;
        parts[i + j*Nx + k*(Nx*Ny)].r = a;
        parts[i + j*Nx + k*(Nx*Ny)].u = 0.;
        parts[i + j*Nx + k*(Nx*Ny)].v = 0.;
        parts[i + j*Nx + k*(Nx*Ny)].w = 0.;
        parts[i + j*Nx + k*(Nx*Ny)].aFx = aFx;
        parts[i + j*Nx + k*(Nx*Ny)].aFy = aFy;
        parts[i + j*Nx + k*(Nx*Ny)].aFz = aFz;
        parts[i + j*Nx + k*(Nx*Ny)].aLx = aLx;
        parts[i + j*Nx + k*(Nx*Ny)].aLy = aLy;
        parts[i + j*Nx + k*(Nx*Ny)].aLz = aLz;
        parts[i + j*Nx + k*(Nx*Ny)].rho = rho;
        parts[i + j*Nx + k*(Nx*Ny)].E = E;
        parts[i + j*Nx + k*(Nx*Ny)].sigma = sigma;
        parts[i + j*Nx + k*(Nx*Ny)].e_dry = e_dry;
        parts[i + j*Nx + k*(Nx*Ny)].l_rough = l_rough; 
        parts[i + j*Nx + k*(Nx*Ny)].order = o;
        parts[i + j*Nx + k*(Nx*Ny)].rs = rs;
        parts[i + j*Nx + k*(Nx*Ny)].spring_k = spring_k;
        parts[i + j*Nx + k*(Nx*Ny)].spring_x = spring_x;
        parts[i + j*Nx + k*(Nx*Ny)].spring_y = spring_y;
        parts[i + j*Nx + k*(Nx*Ny)].spring_z = spring_z;
        parts[i + j*Nx + k*(Nx*Ny)].spring_l = spring_l;
        parts[i + j*Nx + k*(Nx*Ny)].ncoeff = 0.;
        parts[i + j*Nx + k*(Nx*Ny)].translating = t;
        parts[i + j*Nx + k*(Nx*Ny)].rotating = r;       
      }
    }
  }
  printf("Writing part_seeder_array.config...");
  fflush(stdout);
  // write particle configuration to file
  char fname[FILE_NAME_SIZE] = "";
  // open file for writing
  sprintf(fname, "%spart_seeder_array.config", INPUT_DIR);
  FILE *ofile = fopen(fname, "w");
  if(ofile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  // write the number of particles and compact support length
  fprintf(ofile, "n %d\n", nparts);
  fprintf(ofile, "(l/a) %f\n", loa);

  // write each particle configuration
  for(int i = 0; i < nparts; i++) {
    fprintf(ofile, "\n");
    fprintf(ofile, "r %f\n", parts[i].r);
    fprintf(ofile, "(x, y, z) %f %f %f\n", parts[i].x, parts[i].y, parts[i].z);
    fprintf(ofile, "(aFx, aFy, aFz) %f %f %f\n", parts[i].aFx, parts[i].aFy,
      parts[i].aFz);
    fprintf(ofile, "(aLx, aLy, aLz) %f %f %f\n", parts[i].aLx, parts[i].aLy,
      parts[i].aLz);
    fprintf(ofile, "rho %f\n", parts[i].rho);
    fprintf(ofile, "E %f\n", parts[i].E);
    fprintf(ofile, "sigma %f\n", parts[i].sigma);
    fprintf(ofile, "e_dry %f\n", parts[i].e_dry);
    fprintf(ofile, "l_rough %f\n", parts[i].l_rough);
    fprintf(ofile, "order %d\n", parts[i].order);
    fprintf(ofile, "rs/r %f\n", parts[i].rs);
    fprintf(ofile, "spring_k %f\n", parts[i].spring_k);
    fprintf(ofile, "spring (x, y, z) %f %f %f\n",
      parts[i].spring_x, parts[i].spring_y, parts[i].spring_z);
    fprintf(ofile, "spring_l %f\n", parts[i].spring_l);
    fprintf(ofile, "translating %d\n", parts[i].translating);
    fprintf(ofile, "rotating %d\n", parts[i].rotating);
  }

  // close the file
  fclose(ofile);
  printf("done.\n");
  printf("\n...bluebottle seeder done.\n\n");
  fflush(stdout);

  // clean up
  domain_clean();
  parts_clean(); 
}

void seeder_hex(real xExtent, real yExtent, real zExtent, real alpha, real xs, 
  real ys, real zs,
  real loa, real a, real aFx, real aFy, real aFz, real aLx, real aLy, real aLz, 
  real rho, real E, real sigma, real e_dry, real l_rough, int o, real rs, 
  real spring_k, real spring_x, real spring_y, real spring_z, real spring_l, 
  int t, int r)
{

  // read domain input
  domain_read_input();
  domain_init();

  // find separation distance (from wolframalpha, hex close pack proof)
  // Vsp = 8*pi*a^3
  // Acell = 3.2 * R^2 * sqrt(2)
  // hcell = 2*R*sqrt(2/3)
  // alpha = 8*pi*a^3/(3*R^3*sqrt(2))
  // --> R = (8*pi*a^3/(3*eta*sqrt(2)))^(1/3)

  real R = 0.; // horizontal separation between centers
  real h = 0.; // vertical separation between centers

  if (alpha > 0) {
    R = pow(8.*PI*pow(a,3.)/(3.*alpha*sqrt(2.)), 1./3.);
    h = R*sqrt(3.)/2.;

    // check overlap
    if (R < 2.*a) {
      printf("Overlap is wrong. Volume fraction range is 0 to 0.74...\n");
      exit(EXIT_FAILURE);
    }
  } else if (alpha == -1) {
    R = 2.*a;
    h = R*sqrt(3.)/2.;
  } else {
    printf("Unrecognized option.\n");
    exit(EXIT_FAILURE);
  }

  // clear nonsig figs of inputs:
  // -- If domain size is exact multiple of radius and HCP is chosen,
  // -- floating point arithmetic will be slighly off and so 'float' used below
  // -- will round to the wrong number. This should fix that issue.
  real a16 = a*1e16;
  real R16 = R*1e16;
  real h16 = h*1e16;
  real xExtent16 = xExtent*1e16;
  real yExtent16 = yExtent*1e16;
  real zExtent16 = zExtent*1e16;
  
  // find total number of particles in large and small layers
  int Nx_lg = floor((xExtent16 - 2.*a16)/R16) + 1;
  int Ny_lg = floor((yExtent16 - 2.*a16)/R16) + 1;
  int Nx_sm = Nx_lg - 1;
  int Ny_sm = Ny_lg - 1;
  int N_lg = Nx_lg*Ny_lg;
  int N_sm = Nx_sm*Ny_sm;

  // find number of layers of each type
  int nlayers = floor((zExtent16 - 2.*a16)/h16);         // total layers
  int nlay_lg = floor(0.5*nlayers) + (nlayers % 2); // large layers
  int nlay_sm = nlayers - nlay_lg;                  // small layers

  printf("%d large layers: %d x %d\n", nlay_lg, Nx_lg, Ny_lg);
  printf("%d small layers: %d x %d\n", nlay_sm, Nx_sm, Ny_sm);

  nparts = N_lg*nlay_lg + N_sm*nlay_sm;

  printf("The total number of particles is %d \n", nparts);
  printf("Running bluebottle seeder for %d particles...\n", nparts);
  fflush(stdout);

  // allocate particle list
  parts = (part_struct*) malloc(nparts * sizeof(part_struct));
  cpumem += nparts * sizeof(part_struct);

  // calculate starting positions
  real dom_xs = Dom.xs + bc.dsW;
  real dom_xl = Dom.xl - bc.dsE - bc.dsW;
  real dom_ys = Dom.ys + bc.dsS;
  real dom_yl = Dom.yl - bc.dsN - bc.dsS;
  real dom_zs = Dom.zs + bc.dsB;
  real dom_zl = Dom.zl - bc.dsT - bc.dsB;

  // check for interactions
  if ((xs - a) < dom_xs || (ys - a) < dom_ys || (zs - a) < dom_zs) {
    // check domain boundaries
    printf("Starting position creates particle outside of domain\n");
    fflush(stdout);
    printf("...bluebottle seeder failed.\n\n");
    exit(EXIT_FAILURE);
  } else if (((Nx_lg - 1)*R + 2*a) > dom_xl ||
             ((Ny_lg - 1)*R + 2*a) > dom_yl ||
             ((nlayers - 1)*h + 2*a) > dom_zl) {
    // check if particles fit in domain
    printf("Too many particles in the domain for given numbers and spacing\n");
    fflush(stdout);
    printf("...bluebottle seeder failed.\n\n");
    exit(EXIT_FAILURE);
  }

  int index = 0;
  real X = 0;
  real Y = 0;
  real Z = 0;
  // large layers
  for (int k = 0; k < nlay_lg; k++) {
    for (int j = 0; j < Ny_lg; j++) {
      for (int i = 0; i < Nx_lg; i++ ) {
        X = xs + i*R;
        Y = ys + j*R;
        Z = zs + 2.*h*k;

        parts[index].x = X;
        parts[index].y = Y;
        parts[index].z = Z;
        index++;
      }
    }
  }

  // small layers
  for (int k = 0; k < nlay_sm; k++) {
    for (int j = 0; j < Ny_sm; j++) {
      for (int i = 0; i < Nx_sm; i++ ) {

        X = xs + 0.5*R + i*R;
        Y = ys + 0.5*R + j*R;
        Z = zs + h + 2.*h*k;

        parts[index].x = X;
        parts[index].y = Y;
        parts[index].z = Z;
        index++;
      }
    }
  }

  for (int nn = 0; nn < nparts; nn++) {
        parts[nn].u = 0.;
        parts[nn].v = 0.;
        parts[nn].w = 0.;
        parts[nn].aFx = aFx;
        parts[nn].aFy = aFy;
        parts[nn].aFz = aFz;
        parts[nn].aLx = aLx;
        parts[nn].aLy = aLy;
        parts[nn].aLz = aLz;
        parts[nn].rho = rho;
        parts[nn].E = E;
        parts[nn].sigma = sigma;
        parts[nn].e_dry = e_dry;
        parts[nn].l_rough = l_rough;         
        parts[nn].order = o;
        parts[nn].rs = rs;
        parts[nn].spring_k = spring_k;
        parts[nn].spring_x = spring_x;
        parts[nn].spring_y = spring_y;
        parts[nn].spring_z = spring_z;
        parts[nn].spring_l = spring_l;        
        parts[nn].ncoeff = 0.;
        parts[nn].translating = t;
        parts[nn].rotating = r;
  }

  printf("Writing part_seeder_hex.config... ");
  fflush(stdout);
  // write particle configuration to file
  char fname[FILE_NAME_SIZE] = "";
  // open file for writing
  sprintf(fname, "%spart_seeder_hex.config", INPUT_DIR);
  FILE *ofile = fopen(fname, "w");
  if(ofile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  // write the number of particles and compact support length
  fprintf(ofile, "n %d\n", nparts);
  fprintf(ofile, "(l/a) %f\n", loa);

  // write each particle configuration
  for(int i = 0; i < nparts; i++) {
    fprintf(ofile, "\n");
    fprintf(ofile, "r %f\n", a);
    fprintf(ofile, "(x, y, z) %f %f %f\n", parts[i].x, parts[i].y, parts[i].z);
    fprintf(ofile, "(aFx, aFy, aFz) %f %f %f\n", parts[i].aFx, parts[i].aFy,
      parts[i].aFz);
    fprintf(ofile, "(aLx, aLy, aLz) %f %f %f\n", parts[i].aLx, parts[i].aLy,
      parts[i].aLz);
    fprintf(ofile, "rho %f\n", parts[i].rho);
    fprintf(ofile, "E %f\n", parts[i].E);
    fprintf(ofile, "sigma %f\n", parts[i].sigma);
    fprintf(ofile, "e_dry %f\n", parts[i].e_dry);
    fprintf(ofile, "l_rough %f\n", parts[i].l_rough);
    fprintf(ofile, "order %d\n", parts[i].order);
    fprintf(ofile, "rs/r %f\n", parts[i].rs);
    fprintf(ofile, "spring_k %f\n", parts[i].spring_k);
    fprintf(ofile, "spring (x, y, z) %f %f %f\n",
      parts[i].spring_x, parts[i].spring_y, parts[i].spring_z);
    fprintf(ofile, "spring_l %f\n", parts[i].spring_l);
    fprintf(ofile, "translating %d\n", parts[i].translating);
    fprintf(ofile, "rotating %d\n", parts[i].rotating);
  }

  // close the file
  fclose(ofile);
  printf("done.\n");
  printf("\n...bluebottle seeder done.\n\n");
  fflush(stdout);

  // clean up
  domain_clean();
  parts_clean();    
}

void seeder_high_vol_random(int Nx, int Ny, int Nz, real ratio, int nperturb, 
  real loa, real a, real aFx, real aFy, real aFz, real aLx, real aLy, real aLz,
  real rho, real E, real sigma, real e_dry, real l_rough, int o, real rs, 
  real spring_k, real spring_x, real spring_y, real spring_z, real spring_l, 
  int t, int r)
{

  // Generate a random field by perturbing a regular arrary nperturb times and
  // checking for interactions
  // ratio is the max a particle can be perturbed, given as percentage of radius

  nparts = Nx*Ny*Nz;
  printf("Running bluebottle seeder for %d particles...\n\n", nparts);
  fflush(stdout);
  //int fits = 1;
  //int attempts = 1;
  int fail = 0;
  //int redo = 1;
  
  // read domain input
  domain_read_input();
  domain_init();

  // allocate particle list
  parts = (part_struct*) malloc(nparts * sizeof(part_struct));
  cpumem += nparts * sizeof(part_struct);

  // dx in the distance between centers of two nearby particles in x direction
  real dx = Dom.xl/Nx;
  real dy = Dom.yl/Ny;
  real dz = Dom.zl/Nz;

  if(dx < 2.*a) {
    printf(" Too many particles in x direction\n");
    fail = 1;
  }
  if(dy < 2.*a) {
    printf("Too many particles in y direction\n");
    fail = 1; 
  }
  if(dz < 2.*a) {
    printf("Too many particles in z direction\n");
    fail = 1; 
  } 
  if(fail == 1) {
    printf("...bluebottle seeder failed.\n\n");
    exit(EXIT_FAILURE);
  }
  // Set the initial regular domain
  for(int k = 0; k < Nz; k++) {
    for (int j = 0; j < Ny; j++) {
      for (int i = 0; i < Nx; i++) {
        parts[i + j*Nx + k*(Nx*Ny)].x = Dom.xs + (2.*i + 1)*dx*0.5; 
        parts[i + j*Nx + k*(Nx*Ny)].y = Dom.ys + (2.*j + 1)*dy*0.5;
        parts[i + j*Nx + k*(Nx*Ny)].z = Dom.zs + (2.*k + 1)*dz*0.5;  
      }
    }
  }

  real x_new = 0.;
  real y_new = 0.;
  real z_new = 0.;
  real d_min = 100.*a;
  real d_pair= 0.;

  real x_pert = 0.;
  real y_pert = 0.;
  real z_pert = 0.;

  for (int t = 0; t < nperturb; t++) {
    for(int k = 0; k < Nz; k++) {
      for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
          d_min = 100.*a;
          
          x_pert = -1. + 2.*rand() / (real)RAND_MAX;
          y_pert = -1. + 2.*rand() / (real)RAND_MAX;
          z_pert = -1. + 2.*rand() / (real)RAND_MAX;

          x_new = parts[i + j*Nx + k*(Nx*Ny)].x + ratio*a*x_pert;
          y_new = parts[i + j*Nx + k*(Nx*Ny)].y + ratio*a*y_pert;
          z_new = parts[i + j*Nx + k*(Nx*Ny)].z + ratio*a*z_pert;
          if (x_new > Dom.xs && x_new < Dom.xe && 
              y_new > Dom.ys && y_new < Dom.ye && 
              z_new > Dom.zs && z_new < Dom.ze) {
            for (int n = 0; n < Nz; n++) {
              for (int m = 0; m < Ny; m++) {
                for (int l = 0; l < Nx; l++) {
                  // if it calculates the distance to itself
                  if(i == l && j == m && k == n) { 
                    d_pair = 100.*a;
                  } else {
                    d_pair = (x_new - parts[l + m*Nx + n*(Nx*Ny)].x)*
                             (x_new - parts[l + m*Nx + n*(Nx*Ny)].x) +
                             (y_new - parts[l + m*Nx + n*(Nx*Ny)].y)*
                             (y_new - parts[l + m*Nx + n*(Nx*Ny)].y) +
                             (z_new - parts[l + m*Nx + n*(Nx*Ny)].z)*
                             (z_new - parts[l + m*Nx + n*(Nx*Ny)].z);
                    d_pair = sqrt(d_pair);
                  } 
                  if (d_pair < d_min) { 
                    //find the minimum distance between particle pairs
                    d_min = d_pair;
                  } 
                }
              }
            }
            if (d_min > 2.*a) {
            // accept the perturbation if no particle interactions
              parts[i + j*Nx + k*(Nx*Ny)].x = x_new;
              parts[i + j*Nx + k*(Nx*Ny)].y = y_new;
              parts[i + j*Nx + k*(Nx*Ny)].z = z_new;
            } 
          }
        }
      }
    }
  }

 for(int k = 0; k < Nz; k++) {
    for (int j = 0; j < Ny; j++) {
      for (int i = 0; i < Nx; i++) {
        parts[i + j*Nx + k*(Nx*Ny)].r = a;
        parts[i + j*Nx + k*(Nx*Ny)].u = 0.;
        parts[i + j*Nx + k*(Nx*Ny)].v = 0.;
        parts[i + j*Nx + k*(Nx*Ny)].w = 0.;
        parts[i + j*Nx + k*(Nx*Ny)].aFx = aFx;
        parts[i + j*Nx + k*(Nx*Ny)].aFy = aFy;
        parts[i + j*Nx + k*(Nx*Ny)].aFz = aFz;
        parts[i + j*Nx + k*(Nx*Ny)].aLx = aLx;
        parts[i + j*Nx + k*(Nx*Ny)].aLy = aLy;
        parts[i + j*Nx + k*(Nx*Ny)].aLz = aLz;
        parts[i + j*Nx + k*(Nx*Ny)].rho = rho;
        parts[i + j*Nx + k*(Nx*Ny)].E = E;
        parts[i + j*Nx + k*(Nx*Ny)].sigma = sigma;
        parts[i + j*Nx + k*(Nx*Ny)].e_dry = e_dry;
        parts[i + j*Nx + k*(Nx*Ny)].l_rough = l_rough; 
        parts[i + j*Nx + k*(Nx*Ny)].order = o;
        parts[i + j*Nx + k*(Nx*Ny)].rs = rs;
        parts[i + j*Nx + k*(Nx*Ny)].spring_k = spring_k;
        parts[i + j*Nx + k*(Nx*Ny)].spring_x = spring_x;
        parts[i + j*Nx + k*(Nx*Ny)].spring_y = spring_y;
        parts[i + j*Nx + k*(Nx*Ny)].spring_z = spring_z;
        parts[i + j*Nx + k*(Nx*Ny)].spring_l = spring_l;
        parts[i + j*Nx + k*(Nx*Ny)].ncoeff = 0.;
        parts[i + j*Nx + k*(Nx*Ny)].translating = t;
        parts[i + j*Nx + k*(Nx*Ny)].rotating = r;       
      }
    }
  }
  printf("Writing part_seeder_perturbed.config...");
  fflush(stdout);
  // write particle configuration to file
  char fname[FILE_NAME_SIZE] = "";
  // open file for writing
  sprintf(fname, "%spart_seeder_perturbed.config", INPUT_DIR);
  FILE *ofile = fopen(fname, "w");
  if(ofile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  // write the number of particles and compact support length
  fprintf(ofile, "n %d\n", nparts);
  fprintf(ofile, "(l/a) %f\n", loa);

  // write each particle configuration
  for(int i = 0; i < nparts; i++) {
    fprintf(ofile, "\n");
    fprintf(ofile, "r %f\n", parts[i].r);
    fprintf(ofile, "(x, y, z) %f %f %f\n", parts[i].x, parts[i].y, parts[i].z);
    fprintf(ofile, "(aFx, aFy, aFz) %f %f %f\n", parts[i].aFx, parts[i].aFy,
      parts[i].aFz);
    fprintf(ofile, "(aLx, aLy, aLz) %f %f %f\n", parts[i].aLx, parts[i].aLy,
      parts[i].aLz);
    fprintf(ofile, "rho %f\n", parts[i].rho);
    fprintf(ofile, "E %f\n", parts[i].E);
    fprintf(ofile, "sigma %f\n", parts[i].sigma);
    fprintf(ofile, "e_dry %f\n", parts[i].e_dry);
    fprintf(ofile, "l_rough %f\n", parts[i].l_rough);
    fprintf(ofile, "order %d\n", parts[i].order);
    fprintf(ofile, "rs/r %f\n", parts[i].rs);
    fprintf(ofile, "spring_k %f\n", parts[i].spring_k);
    fprintf(ofile, "spring (x, y, z) %f %f %f\n",
      parts[i].spring_x, parts[i].spring_y, parts[i].spring_z);
    fprintf(ofile, "spring_l %f\n", parts[i].spring_l);
    fprintf(ofile, "translating %d\n", parts[i].translating);
    fprintf(ofile, "rotating %d\n", parts[i].rotating);
  }

  // close the file
  fclose(ofile);
  printf("done.\n");
  printf("\n...bluebottle seeder done.\n\n");
  fflush(stdout);

  // clean up
  domain_clean();
  parts_clean(); 
}  
