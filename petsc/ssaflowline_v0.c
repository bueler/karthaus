static const char help[] = 
"\nComputes the velocity from the SSA in the flow-line case.\n\
See mfiles/{flowline.m,ssaflowline.m,testshelf.m} for the problem being solved.\n\n\
This minimal version shows an early-but-working form of the code.  It has\n\
no added user options.  It is representative of what you might get when you\n\
write a first draft based on tutorial PETSc examples.\n\n\
It only works with non-obvious options controlling nonlinear and linear tools,\n\
like these:\n\n\
un-preconditioned matrix-free Newton-Krylov (gmres):\n\
  ./ssaflowline_v0 -snes_mf\n\
un-preconditioned matrix-free Newton-Krylov (gmres) on finer grid, and with\n\
information on convergence; STRUGGLES without these options:\n\
  ./ssaflowline_v0 -da_grid_x 200 -snes_mf -ksp_gmres_restart 300 -mat_mffd_type ds -ksp_converged_reason -snes_monitor \n\
slow assembly of finite-difference matrix, but then ILU-preconditioned\n\
Newton-Krylov; scales better than above but won't scale up much more:\n\
  ./ssaflowline_v0 -snes_fd -mat_fd_type ds -ksp_converged_reason -snes_monitor -da_grid_x 1000 \n\n";

/*
   Include "petscda.h" so that we can use distributed arrays (DAs).
   Include "petscsnes.h" so that we can use SNES solvers.  Note that this
   file automatically includes:
     petsc.h       - base PETSc routines   petscvec.h - vectors
     petscsys.h    - system routines       petscmat.h - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners
     petscksp.h   - linear solvers
*/
#include "petscda.h"
#include "petscsnes.h"

/* User-defined application context - contains data needed by the 
   application-provided call-back routines, esp. FormFunctionLocal().  */
typedef struct {
  DA          da;       /* one-dimensional distributed array data structure for soln and residual */
  PetscReal   secpera;
  PetscReal   n;
  PetscReal   rho;
  PetscReal   rhow;
  PetscReal   g;
  PetscReal   A;
  PetscReal   ug;
  PetscReal   Hg;
  PetscReal   accum;
  PetscReal   L;
  PetscReal   Hcalv;
  PetscReal   uexactcalv;
  PetscReal   gamma;
  Vec         H;
  Vec         uexact;
} AppCtx;

/* Declare the user-defined routines  */
static PetscErrorCode FillThicknessAndExactSoln(AppCtx*);
static PetscErrorCode FormFunctionLocal(DALocalInfo*,PetscScalar*,PetscScalar*,AppCtx*);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  PetscErrorCode         ierr;

  SNES                   snes;                 /* nonlinear solver */
  Vec                    u,r;                  /* solution, residual vectors */
  AppCtx                 user;                 /* user-defined work context */
  PetscInt               its;                  /* iteration count */
  SNESConvergedReason    reason;               /* Check convergence */
  Mat                    B;                    /* preconditioning matrix (only needed with -snes_fd) */

  PetscInitialize(&argc,&argv,(char *)0,help);

  user.n       = 3.0;          /* Glen flow law exponent */
  user.secpera = 31556926.0;
  user.rho     = 900.0;        /* kg m^-3 */
  user.rhow    = 1000.0;       /* kg m^-3 */
  user.g       = 9.8;          /* m s^-2 */
  user.A       = 1.4579e-25;   /* s^-1 Pa^-3; used by MacAyeal et al (1996) */
  user.ug      = 50.0 / user.secpera; /* m s^-1 */
  user.Hg      = 500.0;        /* m */
  user.accum   = 0.3 / user.secpera; /* m s^-1 */
  user.L       = 600.0e3;      /* m */

  /* these are set in FillThicknessAndExactSoln(): */
  user.Hcalv      = -1.0;
  user.uexactcalv = -1.0;
  user.gamma      = -1.0;

  /* Create machinery for parallel grid management (DA), nonlinear solver (SNES), 
     and Vecs for fields (thickness, velocity, RHS).  Note default M=20 grid
     points, degrees of freedom = 1, and stencil width = 1.                                    */
  ierr = DACreate1d(PETSC_COMM_WORLD,DA_NONPERIODIC,-20,1,1,PETSC_NULL,&user.da);CHKERRQ(ierr);
  ierr = DASetUniformCoordinates(user.da,0.0,user.L,
                                 PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);

  /* Extract global vectors from DA and duplicate for remaining same types */
  ierr = DACreateGlobalVector(user.da,&u);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&r);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&user.H);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&user.uexact);CHKERRQ(ierr);

  ierr = DASetLocalFunction(user.da,(DALocalFunction1)FormFunctionLocal);CHKERRQ(ierr);

  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);

  ierr = SNESSetFunction(snes,r,SNESDAFormFunction,&user);CHKERRQ(ierr);

  /* setting up a matrix is only actually needed for -snes_fd case */
  ierr = DAGetMatrix(user.da,MATAIJ,&B);CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,B,B,SNESDAComputeJacobian,PETSC_NULL);CHKERRQ(ierr);

  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  /* We use a particular formula for the thickness H(x).  This formula has the
     property that our numerical velocity solution (u) should converge to a
     known exact solution, which is also computed here.  Note that the 
     thickness is computed on the staggered grid and the exact solution on the
     regular grid.  For more, see notes lecture.pdf.   */
  ierr = FillThicknessAndExactSoln(&user);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,  /* minimal checks on thickness, exact */
           "ssaflowline checks:\n"
           "  user.Hcalv      = %.4f (?= 279.2801 m)\n"
           "  user.uexactcalv = %.4f (?= 734.0302 m/a)\n",
           user.Hcalv, user.uexactcalv * user.secpera);CHKERRQ(ierr);

  /* use exact solution as initial guess; this is cheating but a reasonable
     first step in a scientific-computing project: you *choose* to solve a
     problem to which you already know the exact solution! */
  ierr = VecCopy(user.uexact,u);CHKERRQ(ierr);

  /* Solve nonlinear system  */
  ierr = SNESSolve(snes,PETSC_NULL,u);CHKERRQ(ierr);

  ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr);
  ierr = SNESGetConvergedReason(snes,&reason);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%s Number of Newton iterations = %D\n",
                     SNESConvergedReasons[reason],its);CHKERRQ(ierr);

  ierr = VecDestroy(u);CHKERRQ(ierr);
  ierr = VecDestroy(r);CHKERRQ(ierr);
  ierr = VecDestroy(user.H);CHKERRQ(ierr);
  ierr = VecDestroy(user.uexact);CHKERRQ(ierr);
  ierr = MatDestroy(B);CHKERRQ(ierr);
  ierr = SNESDestroy(snes);CHKERRQ(ierr);
  ierr = DADestroy(user.da);CHKERRQ(ierr);

  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "GetUEx"
static inline PetscScalar GetUEx(PetscScalar ug, PetscScalar qg, PetscScalar a, PetscScalar n, 
                                 PetscScalar Cs, PetscScalar flux) {
  /* MATLAB:
        flux   = a * x + qg;
        uexact = ( ug^(n+1) + (Cs/a) * (flux.^(n+1) - qg^(n+1)) ).^(1/(n+1));    */
  return  PetscPowScalar( PetscPowScalar(ug,n+1.0) + (Cs / a)
                * (PetscPowScalar(flux,n+1.0) - PetscPowScalar(qg,n+1.0)), 1.0 / (n+1.0));
}


#undef __FUNCT__
#define __FUNCT__ "FillThicknessAndExactSoln"
/*  Compute the right thickness H = H(x).  See lecture notes for analysis
    leading to exact ice shelf shape. */
static PetscErrorCode FillThicknessAndExactSoln(AppCtx *user)
{
  PetscErrorCode ierr;
  PetscInt       i,Mx,xs,xm;
  PetscReal      hx, n, r, Cs, xx, flux, qg;
  PetscScalar    *H, *uex;

  PetscFunctionBegin;
  ierr = DAGetInfo(user->da,PETSC_IGNORE,&Mx,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                   PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);
  ierr = DAGetCorners(user->da,&xs,PETSC_NULL,PETSC_NULL,&xm,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);

  /* constants, independent of x */
  hx = user->L / (PetscReal)(Mx-1);
  n  = user->n;
  r  = user->rho / user->rhow;
  Cs = user->A * PetscPowScalar( ( 0.25 * user->rho * user->g * (1.0 - r) ), n);
  qg = user->ug * user->Hg;

  /* Compute regular grid exact soln and staggered-grid thickness over the
     locally-owned part of the grid */
  ierr = DAVecGetArray(user->da,user->uexact,&uex);CHKERRQ(ierr);
  ierr = DAVecGetArray(user->da,user->H,&H);CHKERRQ(ierr);
  for (i=xs; i<xs+xm; i++) {
    xx = hx * (PetscReal)i;  /* = x_i = distance from grounding line */
    flux = user->accum * xx + qg;
    uex[i] = GetUEx(user->ug, qg, user->accum, n, Cs, flux);

    flux += user->accum * hx * 0.5;
    H[i] = flux / GetUEx(user->ug, qg, user->accum, n, Cs, flux);
  }
  ierr = DAVecRestoreArray(user->da,user->uexact,&uex);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(user->da,user->H,&H);CHKERRQ(ierr);

  /* separately compute and store calving-front values */
  flux = user->accum * user->L + qg;
  user->uexactcalv = GetUEx(user->ug, qg, user->accum, n, Cs, flux);
  user->Hcalv = flux / user->uexactcalv;

  /* MATLAB: gamma = ( 0.25 * p.A^(1/n) * (1 - r) * p.rho * p.g * H(end) )^n; */
  user->gamma = 0.25 * PetscPowScalar(user->A,1.0/user->n) * (1.0 - r)
                  * user->rho * user->g * user->Hcalv;
  user->gamma = PetscPowScalar(user->gamma,user->n);
  
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FormFunctionLocal"
/* Compute pointwise residual f(x) over the locally-owned part of the grid
   This is a finite difference method. In TeX, the formula we compute is

f_i = \left|\frac{u_{i+1}-u_i}{dx}\right|^{p-2} H_{i+1/2} (u_{i+1}-u_i)
      - \left|\frac{u_i-u_{i-1}{dx}\right|^{p-2} H_{i-1/2} (u_i-u_{i-1})
      - dx K (H_{i+1/2}^2 - H_{i-1/2}^2)

where

  dx = L / Mx, p = 1+1/n, K = rho * g * (1-r) / (4 * B), r = rho / rhow, B = A^{1/n}
*/
static PetscErrorCode FormFunctionLocal(DALocalInfo *info,PetscScalar *u,PetscScalar *f,AppCtx *user)
{
  PetscErrorCode ierr;
  PetscReal      hx, p, K, B, r, epssqr, duL, duR, sL, sR;
  PetscScalar    *H;
  PetscInt       i, Mx;
  Vec            localH;

  PetscFunctionBegin;

  ierr = DAGetLocalVector(user->da,&localH);CHKERRQ(ierr);
  ierr = DAGlobalToLocalBegin(user->da,user->H,INSERT_VALUES,localH); CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd(user->da,user->H,INSERT_VALUES,localH); CHKERRQ(ierr);

  r  = user->rho / user->rhow;
  p  = 1.0 + 1.0 / user->n;
  B  = PetscPowScalar(user->A,-1.0/user->n);
  K  = user->rho * user->g * (1.0 - r) / (4.0 * B);
  Mx = info->mx;
  hx = user->L / ((PetscReal)Mx - 1.0);

  /* regularize using strain rate of 1.0 / L  per year */
  epssqr = (1.0 / user->secpera) / user->L; 
  epssqr *= epssqr;

  ierr = DAVecGetArray(user->da,localH,&H);CHKERRQ(ierr);
  for (i=info->xs; i<info->xs+info->xm; i++) {
    if (i == 0) {
      f[0] = u[0] - user->ug;  /* Dirichlet condition */
    } else {
      duL = u[i] - u[i-1];
      sL = PetscPowScalar((duL/hx) * (duL/hx) + epssqr,(p-2.0)/2.0) * H[i-1] * duL;
      if (i == Mx-1) {  /* Neumann: calving front stress boundary condition */
        duR = u[Mx-2] + 2.0 * hx * user->gamma - u[Mx-1];
      } else {
        duR = u[i+1] - u[i];
      }
      sR = PetscPowScalar((duR/hx) * (duR/hx) + epssqr,(p-2.0)/2.0) * H[i] * duR;
      f[i] = sR - sL - hx * K * (H[i]*H[i] - H[i-1]*H[i-1]);
    }
  }
  ierr = DAVecRestoreArray(user->da,localH,&H);CHKERRQ(ierr);

  ierr = DARestoreLocalVector(user->da,&localH);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

