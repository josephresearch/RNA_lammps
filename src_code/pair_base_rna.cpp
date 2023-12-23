/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// Contributing author: Dilimulati Aierken, Princeton University, d.aierken@princeton.edu

#include "pair_base_rna.h"

#include "atom.h"

#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"
#include <cmath>
#include <cstring>
#include <math.h>

using namespace LAMMPS_NS;
static constexpr double SMALL = 0.00001;

/* ---------------------------------------------------------------------- */

PairBaseRNA::PairBaseRNA(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
  comm_forward = 1;
}

PairBaseRNA::~PairBaseRNA()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(ubp0);
    memory->destroy(k_theta);
    memory->destroy(k_phi);
    memory->destroy(theta1);
    memory->destroy(theta2);
    memory->destroy(phi1);
    memory->destroy(phi2);
    memory->destroy(offset);
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairBaseRNA::allocate()
{
  allocated = 1;
  int np1 = atom->ntypes + 1;

  memory->create(setflag, np1, np1, "pair:setflag");
  for (int i = 1; i < np1; i++)
    for (int j = i; j < np1; j++) setflag[i][j] = 0;

  memory->create(cutsq, np1, np1, "pair:cutsq");
  memory->create(cut, np1, np1, "pair:cut");
  memory->create(ubp0, np1, np1, "pair:ubp0");
  memory->create(kr, np1, np1, "pair:kr");
  memory->create(r0, np1, np1, "pair:r0");
  memory->create(k_theta, np1, np1, "pair:k_theta");
  memory->create(k_phi, np1, np1, "pair:k_phi");
  memory->create(theta1, np1, np1, "pair:theta1");
  memory->create(theta2, np1, np1, "pair:theta2");
  memory->create(phi1, np1, np1, "pair:phi1");
  memory->create(phi2, np1, np1, "pair:phi2");
  memory->create(offset, np1, np1, "pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairBaseRNA::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR, "Pair style base/angle must have exactly one argument");
  cut_global = utils::numeric(FLERR, arg[0], false, lmp);

  // reset per-type pair cutoffs that have been explicitly set previously

  if (allocated) {
    for (int i = 1; i <= atom->ntypes; i++)
      for (int j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairBaseRNA::coeff(int narg, char **arg)
{
  if (narg < 11 || narg > 12) error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double ubp0_one = utils::numeric(FLERR, arg[2], false, lmp);
  double kr_one = utils::numeric(FLERR, arg[3], false, lmp);
  double r0_one = utils::numeric(FLERR, arg[4], false, lmp);
  double ktheta_one = utils::numeric(FLERR, arg[5], false, lmp);
  double kphi_one = utils::numeric(FLERR, arg[6], false, lmp);
  double theta1_one = utils::numeric(FLERR, arg[7], false, lmp);
  double theta2_one = utils::numeric(FLERR, arg[8], false, lmp);
  double phi1_one = utils::numeric(FLERR, arg[9], false, lmp);
  double phi2_one = utils::numeric(FLERR, arg[10], false, lmp);
  double cut_one = cut_global;

  if (narg == 12) cut_one = utils::numeric(FLERR, arg[11], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      ubp0[i][j] = ubp0_one;
      kr[i][j] = kr_one;
      r0[i][j] = r0_one;
      cut[i][j] = cut_one;
      k_theta[i][j] = ktheta_one;
      k_phi[i][j] = kphi_one;
      theta1[i][j] = theta1_one;
      theta2[i][j] = theta2_one;
      phi1[i][j] = phi1_one;
      phi2[i][j] = phi2_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairBaseRNA::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");

  if (offset_flag) {
    offset[i][j] = 0.00;
  } else
    offset[i][j] = 0.0;

  ubp0[j][i] = ubp0[i][j];
  kr[j][i] = kr[i][j];
  r0[j][i] = r0[i][j];
  k_theta[j][i] = k_theta[i][j];
  k_phi[j][i] = k_phi[i][j];
  theta1[j][i] = theta1[i][j];
  theta2[j][i] = theta2[i][j];
  phi1[j][i] = phi1[i][j];
  phi2[j][i] = phi2[i][j];
  offset[j][i] = offset[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

void PairBaseRNA::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;
  double rsq1, rsq2, r, dr, uexp, factor_lj, ubp1, ubp2, ubp3;
  int *ilist, *jlist, *numneigh, **firstneigh;

  double thetaIJJM, thetaIMIJ, thetaIJJP, thetaIPIJ;
  double phiJMJIIM, phiJPJIIP;
  double costheta, sintheta, cosPhi, sinPhi;
  double fiIJJM[3], fiIMIJ[3], fiIJJP[3], fiIPIJ[3];
  double fjIJJM[3], fjIMIJ[3], fjIJJP[3], fjIPIJ[3];
  double fkIJJM[3], fkIMIJ[3], fkIJJP[3], fkIPIJ[3];
  double xjiIJJM[3], xjkIJJM[3];
  double xjiIMIJ[3], xjkIMIJ[3];
  double xjiIJJP[3], xjkIJJP[3];
  double xjiIPIJ[3], xjkIPIJ[3];
  double xjiJMJIIM[3], xjkJMJIIM[3], xjlJMJIIM[3];
  double xjiJPJIIP[3], xjkJPJIIP[3], xjlJPJIIP[3];
  double fiJMJIIM[3], fiJPJIIP[3];
  double fjJMJIIM[3], fjJPJIIP[3];
  double fkJMJIIM[3], fkJPJIIP[3];
  double flJMJIIM[3], flJPJIIP[3];
  double fir[3], fjr[3], t[3], u[3];
  double xij[3], xkj[3], xlk[3];
  double rij, rkj, rlk, abst, tsq, absu, usq, tcrossu[3];
  double dtheta1, dtheta2, dtheta3, dtheta4;
  double factorTheta, factorTheta11, factorTheta12, factorTheta22;
  double factorPhi, factorPhi11, factorPhi12, factorPhi22;
  int localIP, localIM, localJP, localJM;
  int iGlobal, jGlobal, iChain, jChain;
  double angleSign;
  evdwl = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;    //force
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  //the number of neighbor lists that contains atoms owned by this domain
  inum = list->inum;

  //the (local) indices of the atoms for which neighbor lists have been created.
  ilist = list->ilist;

  // inum sized array with the number of entries of each list of neighbors
  numneigh = list->numneigh;

  //a list of pointers to those lists
  firstneigh = list->firstneigh;

  // *-----------------
  comm->forward_comm(this);    //communicate number of bonds


  // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];     // list[ii]th bead
    xtmp = x[i][0];    //ith lists
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    
    if (atom->nspecial[i][0] == 2) {    //should have two bonds for torsion calcularion
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;
        if (atom->nspecial[j][0] == 2) {    //should have two bonds for torsion calcularion

          iGlobal = atom->tag[i];
          jGlobal = atom->tag[j];
          iChain = atom->molecule[i];
          jChain = atom->molecule[j];
          localIM = domain->closest_image(i, atom->map(iGlobal - 1));      // local id of i-1
          localIP = domain->closest_image(i, atom->map(iGlobal + 1));      // local id of i+1
          localJM = domain->closest_image(i, atom->map(jGlobal - 1));      // local id of j-1
          localJP = domain->closest_image(i, atom->map(jGlobal + 1));      // local id of j+1
          if (abs(iChain - jChain) > 0 || abs(iGlobal - jGlobal) > 4) {    // at least 5 beads apart
            delx = xtmp - x[j][0];
            dely = ytmp - x[j][1];
            delz = ztmp - x[j][2];
            rsq1 = delx * delx + dely * dely + delz * delz;
            jtype = type[j];
            r = sqrt(rsq1);    // |rij| or |rij|
            if (rsq1 < cutsq[itype][jtype] && ubp0[itype][jtype] > 0.1) {
              if (localIM < 0 || localIP < 0 || localJM < 0 || localJP < 0) {

                utils::logmesg(lmp, "{:d}, {:d}, {:d}, {:d}, {:d}, {:d},\t", nlocal + nghost,
                               iGlobal, iChain, jGlobal, jChain, i, j);
                utils::logmesg(lmp, "{:d}, {:d}, {:d}, {:d}, \t", localIM, localIP, localJM,
                               localJP);
                utils::logmesg(lmp, "{:d}, {:d}\t", atom->special[j][0], atom->special[j][1]);
              }

              // *--------harmonic bond part----------------
              //! bonded term
              dr = r - r0[itype][jtype];
              fpair = 2.0 * kr[itype][jtype] * dr;
              fpair *= 1.0 / r;
              fir[0] = fpair * delx;
              fir[1] = fpair * dely;
              fir[2] = fpair * delz;
              fjr[0] = -fir[0];
              fjr[1] = -fir[1];
              fjr[2] = -fir[2];

              // *--------harmonic angle part----------------
              // ! angle between i,j,j-1
              xij[0] = delx;
              xij[1] = dely;
              xij[2] = delz;
              rij = r;

              xkj[0] = x[localJM][0] - x[j][0];
              xkj[1] = x[localJM][1] - x[j][1];
              xkj[2] = x[localJM][2] - x[j][2];

              xjiIJJM[0] = xij[0];
              xjiIJJM[1] = xij[1];
              xjiIJJM[2] = xij[2];

              xjkIJJM[0] = xkj[0];
              xjkIJJM[1] = xkj[1];
              xjkIJJM[2] = xkj[2];

              rsq2 = xkj[0] * xkj[0] + xkj[1] * xkj[1] + xkj[2] * xkj[2];
              rkj = sqrt(rsq2);

              costheta = (xij[0] * xkj[0] + xij[1] * xkj[1] + xij[2] * xkj[2]) / (rij * rkj);

              if (costheta > 1.0) costheta = 1.0;
              if (costheta < -1.0) costheta = -1.0;
              sintheta = sqrt(1 - costheta * costheta);
              if (sintheta < SMALL) sintheta = SMALL;
              thetaIJJM = acos(costheta);

              dtheta1 = thetaIJJM - theta1[itype][jtype];
              factorTheta = -2 * k_theta[itype][jtype] * dtheta1 / sintheta;
              factorTheta11 = factorTheta * costheta / (rsq1);
              factorTheta12 = factorTheta / (rij * rkj);
              factorTheta22 = factorTheta * costheta / (rsq2);

              fiIJJM[0] = factorTheta12 * xkj[0] - factorTheta11 * xij[0];
              fiIJJM[1] = factorTheta12 * xkj[1] - factorTheta11 * xij[1];
              fiIJJM[2] = factorTheta12 * xkj[2] - factorTheta11 * xij[2];

              fkIJJM[0] = factorTheta12 * xij[0] - factorTheta22 * xkj[0];
              fkIJJM[1] = factorTheta12 * xij[1] - factorTheta22 * xkj[1];
              fkIJJM[2] = factorTheta12 * xij[2] - factorTheta22 * xkj[2];

              fjIJJM[0] = 0 - fiIJJM[0] - fkIJJM[0];
              fjIJJM[1] = 0 - fiIJJM[1] - fkIJJM[1];
              fjIJJM[2] = 0 - fiIJJM[2] - fkIJJM[2];

              // ! angle between i-1,i,j
              xij[0] = x[localIM][0] - xtmp;
              xij[1] = x[localIM][1] - ytmp;
              xij[2] = x[localIM][2] - ztmp;
              rsq2 = xij[0] * xij[0] + xij[1] * xij[1] + xij[2] * xij[2];
              rij = sqrt(rsq2);

              xkj[0] = -delx;
              xkj[1] = -dely;
              xkj[2] = -delz;
              rkj = r;

              xjiIMIJ[0] = xij[0];
              xjiIMIJ[1] = xij[1];
              xjiIMIJ[2] = xij[2];

              xjkIMIJ[0] = xkj[0];
              xjkIMIJ[1] = xkj[1];
              xjkIMIJ[2] = xkj[2];

              costheta = (xij[0] * xkj[0] + xij[1] * xkj[1] + xij[2] * xkj[2]) / (rij * rkj);
              if (costheta > 1.0) costheta = 1.0;
              if (costheta < -1.0) costheta = -1.0;
              sintheta = sqrt(1 - costheta * costheta);
              if (sintheta < SMALL) sintheta = SMALL;
              thetaIMIJ = acos(costheta);

              dtheta2 = thetaIMIJ - theta1[itype][jtype];
              factorTheta = -2 * k_theta[itype][jtype] * dtheta2 / sintheta;
              factorTheta11 = factorTheta * costheta / (rsq2);
              factorTheta12 = factorTheta / (rij * rkj);
              factorTheta22 = factorTheta * costheta / (rsq1);

              fiIMIJ[0] = factorTheta12 * xkj[0] - factorTheta11 * xij[0];
              fiIMIJ[1] = factorTheta12 * xkj[1] - factorTheta11 * xij[1];
              fiIMIJ[2] = factorTheta12 * xkj[2] - factorTheta11 * xij[2];

              fkIMIJ[0] = factorTheta12 * xij[0] - factorTheta22 * xkj[0];
              fkIMIJ[1] = factorTheta12 * xij[1] - factorTheta22 * xkj[1];
              fkIMIJ[2] = factorTheta12 * xij[2] - factorTheta22 * xkj[2];

              fjIMIJ[0] = 0 - fiIMIJ[0] - fkIMIJ[0];
              fjIMIJ[1] = 0 - fiIMIJ[1] - fkIMIJ[1];
              fjIMIJ[2] = 0 - fiIMIJ[2] - fkIMIJ[2];

              // ! angle between i,j,j+1
              xij[0] = delx;
              xij[1] = dely;
              xij[2] = delz;
              rij = r;

              xkj[0] = x[localJP][0] - x[j][0];
              xkj[1] = x[localJP][1] - x[j][1];
              xkj[2] = x[localJP][2] - x[j][2];

              xjiIJJP[0] = xij[0];
              xjiIJJP[1] = xij[1];
              xjiIJJP[2] = xij[2];

              xjkIJJP[0] = xkj[0];
              xjkIJJP[1] = xkj[1];
              xjkIJJP[2] = xkj[2];

              rsq2 = xkj[0] * xkj[0] + xkj[1] * xkj[1] + xkj[2] * xkj[2];
              rkj = sqrt(rsq2);

              costheta = (xij[0] * xkj[0] + xij[1] * xkj[1] + xij[2] * xkj[2]) / (rij * rkj);
              if (costheta > 1.0) costheta = 1.0;
              if (costheta < -1.0) costheta = -1.0;
              sintheta = sqrt(1 - costheta * costheta);
              if (sintheta < SMALL) sintheta = SMALL;
              thetaIJJP = acos(costheta);

              dtheta3 = thetaIJJP - theta2[itype][jtype];
              factorTheta = -2 * k_theta[itype][jtype] * dtheta3 / sintheta;
              factorTheta11 = factorTheta * costheta / (rsq1);
              factorTheta12 = factorTheta / (rij * rkj);
              factorTheta22 = factorTheta * costheta / (rsq2);

              fiIJJP[0] = factorTheta12 * xkj[0] - factorTheta11 * xij[0];
              fiIJJP[1] = factorTheta12 * xkj[1] - factorTheta11 * xij[1];
              fiIJJP[2] = factorTheta12 * xkj[2] - factorTheta11 * xij[2];

              fkIJJP[0] = factorTheta12 * xij[0] - factorTheta22 * xkj[0];
              fkIJJP[1] = factorTheta12 * xij[1] - factorTheta22 * xkj[1];
              fkIJJP[2] = factorTheta12 * xij[2] - factorTheta22 * xkj[2];

              fjIJJP[0] = 0 - fiIJJP[0] - fkIJJP[0];
              fjIJJP[1] = 0 - fiIJJP[1] - fkIJJP[1];
              fjIJJP[2] = 0 - fiIJJP[2] - fkIJJP[2];

              // ! angle between i+1,i,j
              xij[0] = x[localIP][0] - xtmp;
              xij[1] = x[localIP][1] - ytmp;
              xij[2] = x[localIP][2] - ztmp;
              rsq2 = xij[0] * xij[0] + xij[1] * xij[1] + xij[2] * xij[2];
              rij = sqrt(rsq2);

              xkj[0] = -delx;
              xkj[1] = -dely;
              xkj[2] = -delz;

              xjiIPIJ[0] = xij[0];
              xjiIPIJ[1] = xij[1];
              xjiIPIJ[2] = xij[2];

              xjkIPIJ[0] = xkj[0];
              xjkIPIJ[1] = xkj[1];
              xjkIPIJ[2] = xkj[2];

              rkj = r;

              costheta = (xij[0] * xkj[0] + xij[1] * xkj[1] + xij[2] * xkj[2]) / (rij * rkj);
              if (costheta > 1.0) costheta = 1.0;
              if (costheta < -1.0) costheta = -1.0;
              sintheta = sqrt(1 - costheta * costheta);
              if (sintheta < SMALL) sintheta = SMALL;
              thetaIPIJ = acos(costheta);

              dtheta4 = thetaIPIJ - theta2[itype][jtype];
              factorTheta = -2 * k_theta[itype][jtype] * dtheta4 / sintheta;
              factorTheta11 = factorTheta * costheta / (rsq2);
              factorTheta12 = factorTheta / (rij * rkj);
              factorTheta22 = factorTheta * costheta / (rsq1);

              fiIPIJ[0] = factorTheta12 * xkj[0] - factorTheta11 * xij[0];
              fiIPIJ[1] = factorTheta12 * xkj[1] - factorTheta11 * xij[1];
              fiIPIJ[2] = factorTheta12 * xkj[2] - factorTheta11 * xij[2];

              fkIPIJ[0] = factorTheta12 * xij[0] - factorTheta22 * xkj[0];
              fkIPIJ[1] = factorTheta12 * xij[1] - factorTheta22 * xkj[1];
              fkIPIJ[2] = factorTheta12 * xij[2] - factorTheta22 * xkj[2];

              fjIPIJ[0] = 0 - fiIPIJ[0] - fkIPIJ[0];
              fjIPIJ[1] = 0 - fiIPIJ[1] - fkIPIJ[1];
              fjIPIJ[2] = 0 - fiIPIJ[2] - fkIPIJ[2];

              // *-----------------dihedral angle part---------------------------

              // ! dihedral between j-1,j,i,i-1
              // r_i - r_j
              xij[0] = x[localJM][0] - x[j][0];
              xij[1] = x[localJM][1] - x[j][1];
              xij[2] = x[localJM][2] - x[j][2];
              // r_k - r_j
              xkj[0] = delx;
              xkj[1] = dely;
              xkj[2] = delz;
              // r_l - r_k
              xlk[0] = x[localIM][0] - xtmp;
              xlk[1] = x[localIM][1] - ytmp;
              xlk[2] = x[localIM][2] - ztmp;

              xjiJMJIIM[0] = xij[0];
              xjiJMJIIM[1] = xij[1];
              xjiJMJIIM[2] = xij[2];

              xjkJMJIIM[0] = xkj[0];
              xjkJMJIIM[1] = xkj[1];
              xjkJMJIIM[2] = xkj[2];

              xjlJMJIIM[0] = x[localIM][0] - x[j][0];
              xjlJMJIIM[1] = x[localIM][1] - x[j][1];
              xjlJMJIIM[2] = x[localIM][2] - x[j][2];

              t[0] = xij[1] * xkj[2] - xij[2] * xkj[1];
              t[1] = xij[2] * xkj[0] - xij[0] * xkj[2];
              t[2] = xij[0] * xkj[1] - xij[1] * xkj[0];
              tsq = t[0] * t[0] + t[1] * t[1] + t[2] * t[2];
              abst = sqrt(tsq);

              u[0] = xlk[1] * xkj[2] - xlk[2] * xkj[1];
              u[1] = xlk[2] * xkj[0] - xlk[0] * xkj[2];
              u[2] = xlk[0] * xkj[1] - xlk[1] * xkj[0];
              usq = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
              absu = sqrt(usq);

              if (tsq > SMALL && usq > SMALL) {
                cosPhi = (t[0] * u[0] + t[1] * u[1] + t[2] * u[2]) / (abst * absu);
                if (cosPhi > 1.0) cosPhi = 1.0;
                if (cosPhi < -1.0) cosPhi = -1.0;
                sinPhi = sqrt(1 - cosPhi * cosPhi);
                if (sinPhi < SMALL) sinPhi = SMALL;
                tcrossu[0] = t[1] * u[2] - t[2] * u[1];
                tcrossu[1] = t[2] * u[0] - t[0] * u[2];
                tcrossu[2] = t[0] * u[1] - t[1] * u[0];
                angleSign = xkj[0] * tcrossu[0] + xkj[1] * tcrossu[1] + xkj[2] * tcrossu[2];
                phiJMJIIM = acos(cosPhi);
                // angleSign = 1;
                if (angleSign < 0) {
                  sinPhi *= -1.0;
                  phiJMJIIM *= -1.0;
                }

                factorPhi = k_phi[itype][jtype] * sin(phiJMJIIM + phi1[itype][jtype]) / sinPhi;
                factorPhi11 = factorPhi * cosPhi / (tsq);
                factorPhi12 = factorPhi / (abst * absu);
                factorPhi22 = factorPhi * cosPhi / (usq);

                fiJMJIIM[0] = factorPhi12 * (xkj[1] * u[2] - xkj[2] * u[1]) +
                    factorPhi11 * (xkj[2] * t[1] - xkj[1] * t[2]);
                fiJMJIIM[1] = factorPhi12 * (xkj[2] * u[0] - xkj[0] * u[2]) +
                    factorPhi11 * (xkj[0] * t[2] - xkj[2] * t[0]);
                fiJMJIIM[2] = factorPhi12 * (xkj[0] * u[1] - xkj[1] * u[0]) +
                    factorPhi11 * (xkj[1] * t[0] - xkj[0] * t[1]);

                flJMJIIM[0] = factorPhi12 * (xkj[1] * t[2] - xkj[2] * t[1]) +
                    factorPhi22 * (xkj[2] * u[1] - xkj[1] * u[2]);
                flJMJIIM[1] = factorPhi12 * (xkj[2] * t[0] - xkj[0] * t[2]) +
                    factorPhi22 * (xkj[0] * u[2] - xkj[2] * u[0]);
                flJMJIIM[2] = factorPhi12 * (xkj[0] * t[1] - xkj[1] * t[0]) +
                    factorPhi22 * (xkj[1] * u[0] - xkj[0] * u[1]);

                fjJMJIIM[0] = factorPhi12 *
                        ((xkj[2] - xij[2]) * u[1] + (xij[1] - xkj[1]) * u[2] - xlk[2] * t[1] +
                         xlk[1] * t[2]) +
                    factorPhi11 * (-(xkj[2] - xij[2]) * t[1] - (xij[1] - xkj[1]) * t[2]) +
                    factorPhi22 * (xlk[2] * u[1] - xlk[1] * u[2]);

                fjJMJIIM[1] = factorPhi12 *
                        ((xij[2] - xkj[2]) * u[0] + (xkj[0] - xij[0]) * u[2] + xlk[2] * t[0] -
                         xlk[0] * t[2]) +
                    factorPhi11 * (-(xij[2] - xkj[2]) * t[0] - (xkj[0] - xij[0]) * t[2]) +
                    factorPhi22 * (-xlk[2] * u[0] + xlk[0] * u[2]);

                fjJMJIIM[2] = factorPhi12 *
                        ((xkj[1] - xij[1]) * u[0] + (xij[0] - xkj[0]) * u[1] - xlk[1] * t[0] +
                         xlk[0] * t[1]) +
                    factorPhi11 * (-(xkj[1] - xij[1]) * t[0] - (xij[0] - xkj[0]) * t[1]) +
                    factorPhi22 * (xlk[1] * u[0] - xlk[0] * u[1]);

                fkJMJIIM[0] = 0 - fiJMJIIM[0] - fjJMJIIM[0] - flJMJIIM[0];
                fkJMJIIM[1] = 0 - fiJMJIIM[1] - fjJMJIIM[1] - flJMJIIM[1];
                fkJMJIIM[2] = 0 - fiJMJIIM[2] - fjJMJIIM[2] - flJMJIIM[2];

                // ! dihedral between j+1,j,i,i+1
                xij[0] = x[localJP][0] - x[j][0];
                xij[1] = x[localJP][1] - x[j][1];
                xij[2] = x[localJP][2] - x[j][2];

                xkj[0] = delx;
                xkj[1] = dely;
                xkj[2] = delz;

                xlk[0] = x[localIP][0] - xtmp;
                xlk[1] = x[localIP][1] - ytmp;
                xlk[2] = x[localIP][2] - ztmp;

                xjiJPJIIP[0] = xij[0];
                xjiJPJIIP[1] = xij[1];
                xjiJPJIIP[2] = xij[2];

                xjkJPJIIP[0] = xkj[0];
                xjkJPJIIP[1] = xkj[1];
                xjkJPJIIP[2] = xkj[2];

                xjlJPJIIP[0] = x[localIP][0] - x[j][0];
                xjlJPJIIP[1] = x[localIP][1] - x[j][1];
                xjlJPJIIP[2] = x[localIP][2] - x[j][2];

                t[0] = xij[1] * xkj[2] - xij[2] * xkj[1];
                t[1] = xij[2] * xkj[0] - xij[0] * xkj[2];
                t[2] = xij[0] * xkj[1] - xij[1] * xkj[0];
                tsq = t[0] * t[0] + t[1] * t[1] + t[2] * t[2];
                abst = sqrt(tsq);

                u[0] = xlk[1] * xkj[2] - xlk[2] * xkj[1];
                u[1] = xlk[2] * xkj[0] - xlk[0] * xkj[2];
                u[2] = xlk[0] * xkj[1] - xlk[1] * xkj[0];
                usq = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
                absu = sqrt(usq);

                if (tsq > SMALL && usq > SMALL) {
                  cosPhi = (t[0] * u[0] + t[1] * u[1] + t[2] * u[2]) / (abst * absu);
                  if (cosPhi > 1.0) cosPhi = 1.0;
                  if (cosPhi < -1.0) cosPhi = -1.0;
                  sinPhi = sqrt(1 - cosPhi * cosPhi);
                  if (sinPhi < SMALL) sinPhi = SMALL;
                  tcrossu[0] = t[1] * u[2] - t[2] * u[1];
                  tcrossu[1] = t[2] * u[0] - t[0] * u[2];
                  tcrossu[2] = t[0] * u[1] - t[1] * u[0];
                  angleSign = xkj[0] * tcrossu[0] + xkj[1] * tcrossu[1] + xkj[2] * tcrossu[2];
                  phiJPJIIP = acos(cosPhi);
                  // angleSign = 1;
                  if (angleSign < 0) {
                    sinPhi *= -1.0;
                    phiJPJIIP *= -1.0;
                  }

                  factorPhi = k_phi[itype][jtype] * sin(phiJPJIIP + phi2[itype][jtype]) / sinPhi;
                  factorPhi11 = factorPhi * cosPhi / (tsq);
                  factorPhi12 = factorPhi / (abst * absu);
                  factorPhi22 = factorPhi * cosPhi / (usq);

                  fiJPJIIP[0] = factorPhi12 * (xkj[1] * u[2] - xkj[2] * u[1]) +
                      factorPhi11 * (xkj[2] * t[1] - xkj[1] * t[2]);
                  fiJPJIIP[1] = factorPhi12 * (xkj[2] * u[0] - xkj[0] * u[2]) +
                      factorPhi11 * (xkj[0] * t[2] - xkj[2] * t[0]);
                  fiJPJIIP[2] = factorPhi12 * (xkj[0] * u[1] - xkj[1] * u[0]) +
                      factorPhi11 * (xkj[1] * t[0] - xkj[0] * t[1]);

                  flJPJIIP[0] = factorPhi12 * (xkj[1] * t[2] - xkj[2] * t[1]) +
                      factorPhi22 * (xkj[2] * u[1] - xkj[1] * u[2]);
                  flJPJIIP[1] = factorPhi12 * (xkj[2] * t[0] - xkj[0] * t[2]) +
                      factorPhi22 * (xkj[0] * u[2] - xkj[2] * u[0]);
                  flJPJIIP[2] = factorPhi12 * (xkj[0] * t[1] - xkj[1] * t[0]) +
                      factorPhi22 * (xkj[1] * u[0] - xkj[0] * u[1]);

                  fjJPJIIP[0] = factorPhi12 *
                          ((xkj[2] - xij[2]) * u[1] + (xij[1] - xkj[1]) * u[2] - xlk[2] * t[1] +
                           xlk[1] * t[2]) +
                      factorPhi11 * (-(xkj[2] - xij[2]) * t[1] - (xij[1] - xkj[1]) * t[2]) +
                      factorPhi22 * (xlk[2] * u[1] - xlk[1] * u[2]);

                  fjJPJIIP[1] = factorPhi12 *
                          ((xij[2] - xkj[2]) * u[0] + (xkj[0] - xij[0]) * u[2] + xlk[2] * t[0] -
                           xlk[0] * t[2]) +
                      factorPhi11 * (-(xij[2] - xkj[2]) * t[0] - (xkj[0] - xij[0]) * t[2]) +
                      factorPhi22 * (-xlk[2] * u[0] + xlk[0] * u[2]);

                  fjJPJIIP[2] = factorPhi12 *
                          ((xkj[1] - xij[1]) * u[0] + (xij[0] - xkj[0]) * u[1] - xlk[1] * t[0] +
                           xlk[0] * t[1]) +
                      factorPhi11 * (-(xkj[1] - xij[1]) * t[0] - (xij[0] - xkj[0]) * t[1]) +
                      factorPhi22 * (xlk[1] * u[0] - xlk[0] * u[1]);

                  fkJPJIIP[0] = 0 - fiJPJIIP[0] - fjJPJIIP[0] - flJPJIIP[0];
                  fkJPJIIP[1] = 0 - fiJPJIIP[1] - fjJPJIIP[1] - flJPJIIP[1];
                  fkJPJIIP[2] = 0 - fiJPJIIP[2] - fjJPJIIP[2] - flJPJIIP[2];

                  //! Energy and Forces
                  //*-------------------------------------------------------------------

                  ubp1 = -kr[itype][jtype] * dr * dr;
                  ubp2 =
                      dtheta1 * dtheta1 + dtheta2 * dtheta2 + dtheta3 * dtheta3 + dtheta4 * dtheta4;
                  ubp2 *= -k_theta[itype][jtype];
                  ubp3 =
                      cos(phiJMJIIM + phi1[itype][jtype]) + cos(phiJPJIIP + phi2[itype][jtype]) + 2;
                  ubp3 *= -k_phi[itype][jtype];

                  uexp = -ubp0[itype][jtype] * exp(ubp1 + ubp2 + ubp3);

                  // ! force on i-1
                  f[localIM][0] += uexp * (fiIMIJ[0] + flJMJIIM[0]);
                  f[localIM][1] += uexp * (fiIMIJ[1] + flJMJIIM[1]);
                  f[localIM][2] += uexp * (fiIMIJ[2] + flJMJIIM[2]);

                  // ! force on i
                  f[i][0] += uexp *
                      (fir[0] + fiIJJM[0] + fjIMIJ[0] + fiIJJP[0] + fjIPIJ[0] + fkJMJIIM[0] +
                       fkJPJIIP[0]);
                  f[i][1] += uexp *
                      (fir[1] + fiIJJM[1] + fjIMIJ[1] + fiIJJP[1] + fjIPIJ[1] + fkJMJIIM[1] +
                       fkJPJIIP[1]);
                  f[i][2] += uexp *
                      (fir[2] + fiIJJM[2] + fjIMIJ[2] + fiIJJP[2] + fjIPIJ[2] + fkJMJIIM[2] +
                       fkJPJIIP[2]);

                  // ! force on i+1
                  f[localIP][0] += uexp * (fiIPIJ[0] + flJPJIIP[0]);
                  f[localIP][1] += uexp * (fiIPIJ[1] + flJPJIIP[1]);
                  f[localIP][2] += uexp * (fiIPIJ[2] + flJPJIIP[2]);

                  // ! force on j-1
                  f[localJM][0] += uexp * (fkIJJM[0] + fiJMJIIM[0]);
                  f[localJM][1] += uexp * (fkIJJM[1] + fiJMJIIM[1]);
                  f[localJM][2] += uexp * (fkIJJM[2] + fiJMJIIM[2]);
                  // }

                  // ! force on J
                  f[j][0] += uexp *
                      (fjr[0] + fjIJJM[0] + fkIMIJ[0] + fjIJJP[0] + fkIPIJ[0] + fjJMJIIM[0] +
                       fjJPJIIP[0]);
                  f[j][1] += uexp *
                      (fjr[1] + fjIJJM[1] + fkIMIJ[1] + fjIJJP[1] + fkIPIJ[1] + fjJMJIIM[1] +
                       fjJPJIIP[1]);
                  f[j][2] += uexp *
                      (fjr[2] + fjIJJM[2] + fkIMIJ[2] + fjIJJP[2] + fkIPIJ[2] + fjJMJIIM[2] +
                       fjJPJIIP[2]);
                  // }

                  // ! force on j+1
                  f[localJP][0] += uexp * (fkIJJP[0] + fiJPJIIP[0]);
                  f[localJP][1] += uexp * (fkIJJP[1] + fiJPJIIP[1]);
                  f[localJP][2] += uexp * (fkIJJP[2] + fiJPJIIP[2]);

                  if (eflag) evdwl = uexp;
                  if (evflag) {
                    ev_tally6(i, localIM, localIP, j, localJM, localJP, evdwl, delx, dely, delz,
                              fir, xjiIJJM, xjkIJJM, fiIJJM, fkIJJM, xjiIMIJ, xjkIMIJ, fiIMIJ,
                              fkIMIJ, xjiIJJP, xjkIJJP, fiIJJP, fkIJJP, xjiIPIJ, xjkIPIJ, fiIPIJ,
                              fkIPIJ, xjiJMJIIM, xjkJMJIIM, xjlJMJIIM, fiJMJIIM, fkJMJIIM, flJMJIIM,
                              xjiJPJIIP, xjkJPJIIP, xjlJPJIIP, fiJPJIIP, fkJPJIIP, flJPJIIP);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

int PairBaseRNA::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int i, j, m;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = atom->nspecial[j][0];
  }
  // utils::logmesg(lmp, "Packed at {:d}, {:d}\n", comm->me, m);
  // utils::logmesg(lmp, "Packed at {:d}, {:d}, {:d}, {:d}, {:d}, {:d}\t", comm->me, m, atom->tag[j], atom->molecule[j], atom->nspecial[j][0]);

  return m;
}

/* ---------------------------------------------------------------------- */

void PairBaseRNA::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) atom->nspecial[i][0] = buf[m++];
  // utils::logmesg(lmp, "unPacked at {:d}, {:d}\n", comm->me, m);

  // utils::logmesg(lmp, "unPacked at {:d}, {:d}, {:d}, {:d}, {:d}, {:d}\n", comm->me, m, atom->tag[i], atom->molecule[i], atom->nspecial[i][0]);
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBaseRNA::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&ubp0[i][j], sizeof(double),1,fp); 
        fwrite(&kr[i][j], sizeof(double),1,fp); 
        fwrite(&r0[i][j], sizeof(double),1,fp);  
        fwrite(&k_theta[i][j], sizeof(double),1,fp); 
        fwrite(&k_phi[i][j], sizeof(double),1,fp); 
        fwrite(&theta1[i][j], sizeof(double),1,fp); 
        fwrite(&theta2[i][j], sizeof(double),1,fp); 
        fwrite(&phi1[i][j], sizeof(double),1,fp); 
        fwrite(&phi2[i][j], sizeof(double),1,fp);
        fwrite(&cut[i][j], sizeof(double),1,fp); 
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */
void PairBaseRNA::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR, &setflag[i][j],sizeof(int),1,fp, nullptr, error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR, &ubp0[i][j], sizeof(double),1,fp, nullptr, error); 
          utils::sfread(FLERR, &kr[i][j], sizeof(double),1,fp, nullptr, error); 
          utils::sfread(FLERR, &r0[i][j], sizeof(double),1,fp, nullptr, error);  
          utils::sfread(FLERR, &k_theta[i][j], sizeof(double),1,fp, nullptr, error); 
          utils::sfread(FLERR, &k_phi[i][j], sizeof(double),1,fp, nullptr, error); 
          utils::sfread(FLERR, &theta1[i][j], sizeof(double),1,fp, nullptr, error); 
          utils::sfread(FLERR, &theta2[i][j], sizeof(double),1,fp, nullptr, error); 
          utils::sfread(FLERR, &phi1[i][j], sizeof(double),1,fp, nullptr, error); 
          utils::sfread(FLERR, &phi2[i][j], sizeof(double),1,fp, nullptr, error);
          utils::sfread(FLERR, &cut[i][j], sizeof(double),1,fp, nullptr, error); 
        }
        MPI_Bcast(&ubp0[i][j], 1,MPI_DOUBLE, 0, world); 
        MPI_Bcast(&kr[i][j], 1,MPI_DOUBLE, 0, world); 
        MPI_Bcast(&r0[i][j], 1,MPI_DOUBLE, 0, world);  
        MPI_Bcast(&k_theta[i][j], 1,MPI_DOUBLE, 0, world); 
        MPI_Bcast(&k_phi[i][j], 1,MPI_DOUBLE, 0, world); 
        MPI_Bcast(&theta1[i][j], 1,MPI_DOUBLE, 0, world); 
        MPI_Bcast(&theta2[i][j], 1,MPI_DOUBLE, 0, world); 
        MPI_Bcast(&phi1[i][j], 1,MPI_DOUBLE, 0, world); 
        MPI_Bcast(&phi2[i][j], 1,MPI_DOUBLE, 0, world);
        MPI_Bcast(&cut[i][j], 1,MPI_DOUBLE, 0, world); 
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBaseRNA::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBaseRNA::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR, &cut_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR, &offset_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR, &mix_flag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairBaseRNA::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g %g %g %g %g %g %g\n",i, ubp0[i][i], kr[i][i], r0[i][i], 
              k_theta[i][i], theta1[i][i], theta2[i][i], phi1[i][i], phi2[i][i], cut[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairBaseRNA::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g %g %g %g %g %g\n",i, j, ubp0[i][j], kr[i][j], r0[i][j], 
              k_theta[i][j], theta1[i][j], theta2[i][j], phi1[i][j], phi2[i][j], cut[i][j]);
}




/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into global or per-atom accumulators
   called by BaseRNA potential, newton_pair is always on
 ------------------------------------------------------------------------- */

void PairBaseRNA::ev_tally6(int i, int im, int ip, int j, int jm, int jp, double evdwl,
double delx, double dely, double delz, double *fpair,
double *xjiIJJM, double *xjkIJJM, double *fiIJJM, double *fkIJJM, 
double *xjiIMIJ, double *xjkIMIJ, double *fiIMIJ, double *fkIMIJ, 
double *xjiIJJP, double *xjkIJJP, double *fiIJJP, double *fkIJJP, 
double *xjiIPIJ, double *xjkIPIJ, double *fiIPIJ, double *fkIPIJ, 
double *xjiJMJIIM, double *xjkJMJIIM, double *xjlJMJIIM,
double *fiJMJIIM, double *fkJMJIIM, double *flJMJIIM,
double *xjiJPJIIP, double *xjkJPJIIP, double *xjlJPJIIP,
double *fiJPJIIP, double *fkJPJIIP, double *flJPJIIP
)
{
  double epairsixth,v[6];
  if (eflag_either) {
    if (eflag_global) eng_vdwl += evdwl;
    if (eflag_atom) {
      epairsixth = 0.16666666666 * evdwl;
      eatom[i] += epairsixth;
      eatom[im] += epairsixth;
      eatom[ip] += epairsixth;
      eatom[j] += epairsixth;
      eatom[jm] += epairsixth;
      eatom[jp] += epairsixth;
    }
  }

  if (vflag_either) {
    //! from the attraction part between i and j
    v[0] = delx*fpair[0];
    v[1] = dely*fpair[1];
    v[2] = delz*fpair[2];
    v[3] = delx*fpair[1];
    v[4] = delx*fpair[2];
    v[5] = dely*fpair[2];

    // ! angle between i,j,j-1
    v[0] += xjiIJJM[0]*fiIJJM[0] + xjkIJJM[0]*fkIJJM[0];
    v[1] += xjiIJJM[1]*fiIJJM[1] + xjkIJJM[1]*fkIJJM[1];
    v[2] += xjiIJJM[2]*fiIJJM[2] + xjkIJJM[2]*fkIJJM[2];
    v[3] += xjiIJJM[0]*fiIJJM[1] + xjkIJJM[0]*fkIJJM[1];
    v[4] += xjiIJJM[0]*fiIJJM[2] + xjkIJJM[0]*fkIJJM[2];
    v[5] += xjiIJJM[1]*fiIJJM[2] + xjkIJJM[1]*fkIJJM[2];

    // ! angle between i-1,i,j
    v[0] += xjiIMIJ[0]*fiIMIJ[0] + xjkIMIJ[0]*fkIMIJ[0];
    v[1] += xjiIMIJ[1]*fiIMIJ[1] + xjkIMIJ[1]*fkIMIJ[1];
    v[2] += xjiIMIJ[2]*fiIMIJ[2] + xjkIMIJ[2]*fkIMIJ[2];
    v[3] += xjiIMIJ[0]*fiIMIJ[1] + xjkIMIJ[0]*fkIMIJ[1];
    v[4] += xjiIMIJ[0]*fiIMIJ[2] + xjkIMIJ[0]*fkIMIJ[2];
    v[5] += xjiIMIJ[1]*fiIMIJ[2] + xjkIMIJ[1]*fkIMIJ[2];

    // ! angle between i,j,j+1
    v[0] += xjiIJJP[0]*fiIJJP[0] + xjkIJJP[0]*fkIJJP[0];
    v[1] += xjiIJJP[1]*fiIJJP[1] + xjkIJJP[1]*fkIJJP[1];
    v[2] += xjiIJJP[2]*fiIJJP[2] + xjkIJJP[2]*fkIJJP[2];
    v[3] += xjiIJJP[0]*fiIJJP[1] + xjkIJJP[0]*fkIJJP[1];
    v[4] += xjiIJJP[0]*fiIJJP[2] + xjkIJJP[0]*fkIJJP[2];
    v[5] += xjiIJJP[1]*fiIJJP[2] + xjkIJJP[1]*fkIJJP[2];

    // ! angle between i+1,i,j
    v[0] += xjiIPIJ[0]*fiIPIJ[0] + xjkIPIJ[0]*fkIPIJ[0];
    v[1] += xjiIPIJ[1]*fiIPIJ[1] + xjkIPIJ[1]*fkIPIJ[1];
    v[2] += xjiIPIJ[2]*fiIPIJ[2] + xjkIPIJ[2]*fkIPIJ[2];
    v[3] += xjiIPIJ[0]*fiIPIJ[1] + xjkIPIJ[0]*fkIPIJ[1];
    v[4] += xjiIPIJ[0]*fiIPIJ[2] + xjkIPIJ[0]*fkIPIJ[2];
    v[5] += xjiIPIJ[1]*fiIPIJ[2] + xjkIPIJ[1]*fkIPIJ[2];


    // ! dihedral between j-1,j,i,i-1
    v[0] += xjiJMJIIM[0]*fiJMJIIM[0] + xjkJMJIIM[0]*fkJMJIIM[0] + xjlJMJIIM[0]*flJMJIIM[0];
    v[1] += xjiJMJIIM[1]*fiJMJIIM[1] + xjkJMJIIM[1]*fkJMJIIM[1] + xjlJMJIIM[1]*flJMJIIM[1];
    v[2] += xjiJMJIIM[2]*fiJMJIIM[2] + xjkJMJIIM[2]*fkJMJIIM[2] + xjlJMJIIM[2]*flJMJIIM[2];
    v[3] += xjiJMJIIM[0]*fiJMJIIM[1] + xjkJMJIIM[0]*fkJMJIIM[1] + xjlJMJIIM[0]*flJMJIIM[1];
    v[4] += xjiJMJIIM[0]*fiJMJIIM[2] + xjkJMJIIM[0]*fkJMJIIM[2] + xjlJMJIIM[0]*flJMJIIM[2];
    v[5] += xjiJMJIIM[1]*fiJMJIIM[2] + xjkJMJIIM[1]*fkJMJIIM[2] + xjlJMJIIM[1]*flJMJIIM[2];
    

    // ! dihedral between j+1,j,i,i+1
    v[0] += xjiJPJIIP[0]*fiJPJIIP[0] + xjkJPJIIP[0]*fkJPJIIP[0] + xjlJPJIIP[0]*flJPJIIP[0];
    v[1] += xjiJPJIIP[1]*fiJPJIIP[1] + xjkJPJIIP[1]*fkJPJIIP[1] + xjlJPJIIP[1]*flJPJIIP[1];
    v[2] += xjiJPJIIP[2]*fiJPJIIP[2] + xjkJPJIIP[2]*fkJPJIIP[2] + xjlJPJIIP[2]*flJPJIIP[2];
    v[3] += xjiJPJIIP[0]*fiJPJIIP[1] + xjkJPJIIP[0]*fkJPJIIP[1] + xjlJPJIIP[0]*flJPJIIP[1];
    v[4] += xjiJPJIIP[0]*fiJPJIIP[2] + xjkJPJIIP[0]*fkJPJIIP[2] + xjlJPJIIP[0]*flJPJIIP[2];
    v[5] += xjiJPJIIP[1]*fiJPJIIP[2] + xjkJPJIIP[1]*fkJPJIIP[2] + xjlJPJIIP[1]*flJPJIIP[2];

    //! coeffecient
    v[0] *=  evdwl;
    v[1] *=  evdwl;
    v[2] *=  evdwl;
    v[3] *=  evdwl;
    v[4] *=  evdwl;
    v[5] *=  evdwl;
    

    if (vflag_global) {
      virial[0] += v[0];
      virial[1] += v[1];
      virial[2] += v[2];
      virial[3] += v[3];
      virial[4] += v[4];
      virial[5] += v[5];
    }

    if (vflag_atom) {
      v[0] *=  0.16666666666;
      v[1] *=  0.16666666666;
      v[2] *=  0.16666666666;
      v[3] *=  0.16666666666;
      v[4] *=  0.16666666666;
      v[5] *=  0.16666666666;

      vatom[i][0] += v[0]; vatom[i][1] += v[1]; vatom[i][2] += v[2];
      vatom[i][3] += v[3]; vatom[i][4] += v[4]; vatom[i][5] += v[5];
      
      vatom[j][0] += v[0]; vatom[j][1] += v[1]; vatom[j][2] += v[2];
      vatom[j][3] += v[3]; vatom[j][4] += v[4]; vatom[j][5] += v[5];

      vatom[im][0] += v[0]; vatom[im][1] += v[1]; vatom[im][2] += v[2];
      vatom[im][3] += v[3]; vatom[im][4] += v[4]; vatom[im][5] += v[5];
      
      vatom[jm][0] += v[0]; vatom[jm][1] += v[1]; vatom[jm][2] += v[2];
      vatom[jm][3] += v[3]; vatom[jm][4] += v[4]; vatom[jm][5] += v[5];

      vatom[ip][0] += v[0]; vatom[ip][1] += v[1]; vatom[ip][2] += v[2];
      vatom[ip][3] += v[3]; vatom[ip][4] += v[4]; vatom[ip][5] += v[5];
      
      vatom[jp][0] += v[0]; vatom[jp][1] += v[1]; vatom[jp][2] += v[2];
      vatom[jp][3] += v[3]; vatom[jp][4] += v[4]; vatom[jp][5] += v[5];
    }
  }
}