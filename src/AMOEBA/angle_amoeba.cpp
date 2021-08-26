// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "angle_amoeba.h"

#include <cmath>
#include <cstring>
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define SMALL 0.001

/* ---------------------------------------------------------------------- */

AngleAmoeba::AngleAmoeba(LAMMPS *lmp) : Angle(lmp) 
{
  pflag = nullptr;
  theta0 = nullptr;
  k2 = nullptr;
  k3 = nullptr;
  k4 = nullptr;
  k5 = nullptr;
  k6 = nullptr;
}

/* ---------------------------------------------------------------------- */

AngleAmoeba::~AngleAmoeba()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);

    memory->destroy(pflag);
    memory->destroy(theta0);
    memory->destroy(k2);
    memory->destroy(k3);
    memory->destroy(k4);
    memory->destroy(k5);
    memory->destroy(k6);
  }
}

/* ---------------------------------------------------------------------- */

void AngleAmoeba::compute(int eflag, int vflag)
{
  int i1,i2,i3,n,type,tflag;
  double delx1,dely1,delz1,delx2,dely2,delz2;
  double eangle,f1[3],f3[3];
  double dtheta,dtheta2,dtheta3,dtheta4,dtheta5,dtheta6,de_angle;
  double dr1,dr2,tk1,tk2,aa1,aa2,aa11,aa12,aa21,aa22;
  double rsq1,rsq2,r1,r2,c,s,a,a11,a12,a22,b1,b2;

  eangle = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nanglelist; n++) {
    i1 = anglelist[n][0];
    i2 = anglelist[n][1];
    i3 = anglelist[n][2];
    type = anglelist[n][3];

    // tflag = 0 for "angle", 1 for "anglep" in Tinker PRM file
    // atom 2 must have 3 bond partners to invoke "anglep" option

    tflag = pflag[type];

    if (tflag && atom->num_bond[i2] == 3) {
      anglep(i1,i2,i3,type,eflag);
      continue;
    }

    // 1st bond

    delx1 = x[i1][0] - x[i2][0];
    dely1 = x[i1][1] - x[i2][1];
    delz1 = x[i1][2] - x[i2][2];

    rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    r1 = sqrt(rsq1);

    // 2nd bond

    delx2 = x[i3][0] - x[i2][0];
    dely2 = x[i3][1] - x[i2][1];
    delz2 = x[i3][2] - x[i2][2];

    rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
    r2 = sqrt(rsq2);

    // angle (cos and sin)

    c = delx1*delx2 + dely1*dely2 + delz1*delz2;
    c /= r1*r2;

    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    s = sqrt(1.0 - c*c);
    if (s < SMALL) s = SMALL;
    s = 1.0/s;

    // force & energy for angle term

    dtheta = acos(c) - theta0[type];
    dtheta2 = dtheta*dtheta;
    dtheta3 = dtheta2*dtheta;
    dtheta4 = dtheta3*dtheta;
    dtheta5 = dtheta4*dtheta;
    dtheta6 = dtheta5*dtheta;

    de_angle = 2.0*k2[type]*dtheta + 3.0*k3[type]*dtheta2 +
      4.0*k4[type]*dtheta3 + 5.0*k5[type]*dtheta4 + 6.0*k6[type]*dtheta5;

    a = -de_angle*s;
    a11 = a*c / rsq1;
    a12 = -a / (r1*r2);
    a22 = a*c / rsq2;

    f1[0] = a11*delx1 + a12*delx2;
    f1[1] = a11*dely1 + a12*dely2;
    f1[2] = a11*delz1 + a12*delz2;

    f3[0] = a22*delx2 + a12*delx1;
    f3[1] = a22*dely2 + a12*dely1;
    f3[2] = a22*delz2 + a12*delz1;

    if (eflag) eangle = k2[type]*dtheta2 + k3[type]*dtheta3 + 
                 k4[type]*dtheta4 + k5[type]*dtheta5 + k6[type]*dtheta6;

    // apply force to each of 3 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += f1[0];
      f[i1][1] += f1[1];
      f[i1][2] += f1[2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= f1[0] + f3[0];
      f[i2][1] -= f1[1] + f3[1];
      f[i2][2] -= f1[2] + f3[2];
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += f3[0];
      f[i3][1] += f3[1];
      f[i3][2] += f3[2];
    }

    if (evflag) ev_tally(i1,i2,i3,nlocal,newton_bond,eangle,f1,f3,
                         delx1,dely1,delz1,delx2,dely2,delz2);
  }
}

/* ---------------------------------------------------------------------- */

void AngleAmoeba::anglep(int i1, int i2, int i3, int type, int eflag)
{
  int i4;
  tagint i1tag,i3tag,i4tag;
  double xia,yia,zia,xib,yib,zib,xic,yic,zic,xid,yid,zid;
  double xad,yad,zad,xbd,ybd,zbd,xcd,ycd,zcd;
  double xt,yt,zt,rt2;
  double xip,yip,zip,xap,yap,zap,xcp,ycp,zcp;
  double rap2,rcp2;
  double dtheta,dtheta2,dtheta3,dtheta4,dtheta5,dtheta6;
  double xm,ym,zm,rm,dot;
  double cosine,angle;
  double eangle,deddt;
  double dedxip,dedyip,dedzip,dpdxia,dpdyia,dpdzia,dpdxic,dpdyic,dpdzic;
  double delta,delta2,ptrt2,term,terma,termc;
  double f1[3],f2[3],f3[3],f4[3];

  double **x = atom->x;
  double **f = atom->f;
  tagint **bond_atom = atom->bond_atom;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  // i4 = index of third atom that i2 is bonded to
  
  i1tag = atom->tag[i1];
  i3tag = atom->tag[i3];

  for (int ibond = 0; ibond < 3; ibond++) {
    i4tag = bond_atom[i2][ibond];
    if (i4tag != i1tag && i4tag != i3tag) break;
  }

  i4 = atom->map(i4tag);
  i4 = domain->closest_image(i2,i4);

  // anglep out-of-plane calculation from Tinker

  xia = x[i1][0];
  yia = x[i1][1];
  zia = x[i1][2];
  xib = x[i2][0];
  yib = x[i2][1];
  zib = x[i2][2];
  xic = x[i3][0];
  yic = x[i3][1];
  zic = x[i3][2];
  xid = x[i4][0];
  yid = x[i4][1];
  zid = x[i4][2];

  xad = xia - xid;
  yad = yia - yid;
  zad = zia - zid;
  xbd = xib - xid;
  ybd = yib - yid;
  zbd = zib - zid;
  xcd = xic - xid;
  ycd = yic - yid;
  zcd = zic - zid;

  xt = yad*zcd - zad*ycd;
  yt = zad*xcd - xad*zcd;
  zt = xad*ycd - yad*xcd;
  rt2 = xt*xt + yt*yt + zt*zt;
  delta = -(xt*xbd + yt*ybd + zt*zbd) / rt2;
  xip = xib + xt*delta;
  yip = yib + yt*delta;
  zip = zib + zt*delta;
  xap = xia - xip;
  yap = yia - yip;
  zap = zia - zip;
  xcp = xic - xip;
  ycp = yic - yip;
  zcp = zic - zip;
  rap2 = xap*xap + yap*yap + zap*zap;
  rcp2 = xcp*xcp + ycp*ycp + zcp*zcp;

  // NOTE: can these be 0.0 ?  what to do?

  if (rap2 == 0.0 || rcp2 == 0.0) return;

  xm = ycp*zap - zcp*yap;
  ym = zcp*xap - xcp*zap;
  zm = xcp*yap - ycp*xap;
  rm = sqrt(xm*xm + ym*ym + zm*zm);
  rm = MAX(rm,0.0001);
  dot = xap*xcp + yap*ycp + zap*zcp;
  cosine = dot / sqrt(rap2*rcp2);
  cosine = MIN(1.0,MAX(-1.0,cosine));

  // force & energy for angle term
  
  dtheta = acos(cosine) - theta0[type];
  dtheta2 = dtheta*dtheta;
  dtheta3 = dtheta2*dtheta;
  dtheta4 = dtheta3*dtheta;
  dtheta5 = dtheta4*dtheta;
  dtheta6 = dtheta5*dtheta;

  deddt = 2.0*k2[type]*dtheta + 3.0*k3[type]*dtheta2 +
    4.0*k4[type]*dtheta3 + 5.0*k5[type]*dtheta4 + 6.0*k6[type]*dtheta5;

  if (eflag) eangle = k2[type]*dtheta2 + k3[type]*dtheta3 + 
               k4[type]*dtheta4 + k5[type]*dtheta5 + k6[type]*dtheta6;

  printf("ANGLEP: %d %d %d %d: eng %g\n",
         atom->tag[i1],
         atom->tag[i2],
         atom->tag[i3],
         atom->tag[i4],
         eangle);

  // chain rule terms for first derivative components

  terma = -deddt / (rap2*rm);
  termc = deddt / (rcp2*rm);
  f1[0] = terma * (yap*zm-zap*ym);
  f1[1] = terma * (zap*xm-xap*zm);
  f1[2] = terma * (xap*ym-yap*xm);
  f3[0] = termc * (ycp*zm-zcp*ym);
  f3[1] = termc * (zcp*xm-xcp*zm);
  f3[2] = termc * (xcp*ym-ycp*xm);
  dedxip = -f1[0] - f3[0];
  dedyip = -f1[1] - f3[1];
  dedzip = -f1[2] - f3[2];

  // chain rule components for the projection of the central atom

  delta2 = 2.0 * delta;
  ptrt2 = (dedxip*xt + dedyip*yt + dedzip*zt) / rt2;
  term = (zcd*ybd-ycd*zbd) + delta2*(yt*zcd-zt*ycd);
  dpdxia = delta*(ycd*dedzip-zcd*dedyip) + term*ptrt2;
  term = (xcd*zbd-zcd*xbd) + delta2*(zt*xcd-xt*zcd);
  dpdyia = delta*(zcd*dedxip-xcd*dedzip) + term*ptrt2;
  term = (ycd*xbd-xcd*ybd) + delta2*(xt*ycd-yt*xcd);
  dpdzia = delta*(xcd*dedyip-ycd*dedxip) + term*ptrt2;
  term = (yad*zbd-zad*ybd) + delta2*(zt*yad-yt*zad);
  dpdxic = delta*(zad*dedyip-yad*dedzip) + term*ptrt2;
  term = (zad*xbd-xad*zbd) + delta2*(xt*zad-zt*xad);
  dpdyic = delta*(xad*dedzip-zad*dedxip) + term*ptrt2;
  term = (xad*ybd-yad*xbd) + delta2*(yt*xad-xt*yad);
  dpdzic = delta*(yad*dedxip-xad*dedyip) + term*ptrt2;

  // compute derivative components for this interaction

  f1[0] += dpdxia;
  f1[1] += dpdyia;
  f1[2] += dpdzia;
  f2[0] = dedxip;
  f2[1] = dedyip;
  f2[2] = dedzip;
  f3[0] += dpdxic;
  f3[1] += dpdyic;
  f3[2] += dpdzic;
  f4[0] = -f1[0] - f2[0] - f3[0];
  f4[1] = -f1[1] - f2[1] - f3[1];
  f4[2] = -f1[2] - f2[2] - f3[2];

  // apply force to each of 4 atoms

  if (newton_bond || i1 < nlocal) {
    f[i1][0] += f1[0];
    f[i1][1] += f1[1];
    f[i1][2] += f1[2];
  }

  if (newton_bond || i2 < nlocal) {
    f[i2][0] += f2[0];
    f[i2][1] += f2[1];
    f[i2][2] += f2[2];
  }

  if (newton_bond || i3 < nlocal) {
    f[i3][0] += f3[0];
    f[i3][1] += f3[1];
    f[i3][2] += f3[2];
  }

  if (newton_bond || i4 < nlocal) {
    f[i4][0] += f4[0];
    f[i4][1] += f4[1];
    f[i4][2] += f4[2];
  }

  if (evflag) ev_tally4(i1,i2,i3,14,nlocal,newton_bond,eangle,f1,f2,f3,f4);
}

/* ---------------------------------------------------------------------- */

void AngleAmoeba::allocate()
{
  allocated = 1;
  int n = atom->nangletypes;

  memory->create(pflag,n+1,"angle:pflag");
  memory->create(theta0,n+1,"angle:theta0");
  memory->create(k2,n+1,"angle:k2");
  memory->create(k3,n+1,"angle:k3");
  memory->create(k4,n+1,"angle:k4");
  memory->create(k5,n+1,"angle:k5");
  memory->create(k6,n+1,"angle:k6");

  memory->create(setflag,n+1,"angle:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void AngleAmoeba::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  int ilo,ihi;
  utils::bounds(FLERR,arg[0],1,atom->nangletypes,ilo,ihi,error);

  int count = 0;
  
  if (narg != 8) error->all(FLERR,"Incorrect args for angle coefficients");

  int pflag_one = utils::inumeric(FLERR,arg[1],false,lmp);
  double theta0_one = utils::numeric(FLERR,arg[2],false,lmp);
  double k2_one = utils::numeric(FLERR,arg[3],false,lmp);
  double k3_one = utils::numeric(FLERR,arg[4],false,lmp);
  double k4_one = utils::numeric(FLERR,arg[5],false,lmp);
  double k5_one = utils::numeric(FLERR,arg[6],false,lmp);
  double k6_one = utils::numeric(FLERR,arg[7],false,lmp);

  // convert theta0 from degrees to radians
  
  for (int i = ilo; i <= ihi; i++) {
    pflag[i] = pflag_one;
    theta0[i] = theta0_one/180.0 * MY_PI;
    k2[i] = k2_one;
    k3[i] = k3_one;
    k4[i] = k4_one;
    k5[i] = k5_one;
    k6[i] = k6_one;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for angle coefficients");

  for (int i = ilo; i <= ihi; i++) setflag[i] = 1;
}

/* ---------------------------------------------------------------------- */

double AngleAmoeba::equilibrium_angle(int i)
{
  return theta0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void AngleAmoeba::write_restart(FILE *fp)
{
  fwrite(&pflag[1],sizeof(int),atom->nangletypes,fp);
  fwrite(&theta0[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&k2[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&k3[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&k4[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&k5[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&k6[1],sizeof(double),atom->nangletypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void AngleAmoeba::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR,&pflag[1],sizeof(int),atom->nangletypes,fp,nullptr,error);
    utils::sfread(FLERR,&theta0[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
    utils::sfread(FLERR,&k2[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
    utils::sfread(FLERR,&k3[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
    utils::sfread(FLERR,&k4[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
    utils::sfread(FLERR,&k5[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
    utils::sfread(FLERR,&k6[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
  }

  MPI_Bcast(&pflag[1],atom->nangletypes,MPI_INT,0,world);
  MPI_Bcast(&theta0[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&k2[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&k3[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&k4[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&k5[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&k6[1],atom->nangletypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nangletypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void AngleAmoeba::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nangletypes; i++)
    fprintf(fp,"%d %d %g %g %g %g %g %g\n",
            i,pflag[i],theta0[i]/MY_PI*180.0,k2[i],k3[i],k4[i],k5[i],k6[i]);
}

/* ---------------------------------------------------------------------- */

double AngleAmoeba::single(int type, int i1, int i2, int i3)
{
  double **x = atom->x;

  double delx1 = x[i1][0] - x[i2][0];
  double dely1 = x[i1][1] - x[i2][1];
  double delz1 = x[i1][2] - x[i2][2];
  domain->minimum_image(delx1,dely1,delz1);
  double r1 = sqrt(delx1*delx1 + dely1*dely1 + delz1*delz1);

  double delx2 = x[i3][0] - x[i2][0];
  double dely2 = x[i3][1] - x[i2][1];
  double delz2 = x[i3][2] - x[i2][2];
  domain->minimum_image(delx2,dely2,delz2);
  double r2 = sqrt(delx2*delx2 + dely2*dely2 + delz2*delz2);

  double c = delx1*delx2 + dely1*dely2 + delz1*delz2;
  c /= r1*r2;
  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

  double s = sqrt(1.0 - c*c);
  if (s < SMALL) s = SMALL;
  s = 1.0/s;

  double dtheta = acos(c) - theta0[type];
  double dtheta2 = dtheta*dtheta;
  double dtheta3 = dtheta2*dtheta;
  double dtheta4 = dtheta3*dtheta;
  double dtheta5 = dtheta4*dtheta;
  double dtheta6 = dtheta5*dtheta;

  double energy = k2[type]*dtheta2 + k3[type]*dtheta3 + k4[type]*dtheta4
           + k5[type]*dtheta5 + k6[type]*dtheta6;

  return energy;
}