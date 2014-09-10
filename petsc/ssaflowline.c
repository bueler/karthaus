static const char help[] = 
"\nComputes the velocity from the SSA model for ice, in a flow-line case.\n\
See ../notes/notes.pdf and ../slides.pdf and\n\
../mfiles/{flowline.m,ssaflowline.m,testshelf.m} for the Matlab/Octave version\n\
of the same problem.\n\
\n\
It works in PETSc 3.5.2.\n\
Compare examples src/snes/examples/tutorials/{ex5.c,ex33.c}.\n\
\n\
Program usage:  mpiexec -n <procs> ./ssaflowline [-help] [all PETSc options]\n\
\n\
Examples:\n\
*  see all SSA-related options, and default values of constants:\n\
      ./ssaflowline -help |grep ssa_\n\
*  use defaults (Newton-Krylov w analytical Jacobian and M=21 grid points):\n\
      ./ssaflowline\n\
*  see info on PETSc objects, and run in parallel:\n\
      mpiexec -n 2 ./ssaflowline -snes_view -snes_monitor -ksp_converged_reason\n\
*  much finer grid\n\
      ./ssaflowline -snes_monitor -da_grid_x 20001\n\
*  finer grid and Picard iterations: slower than Newton!\n\
      ./ssaflowline -snes_monitor -da_grid_x 20001 -ssa_picard\n\
*  finite difference Jacobian:\n\
      ./ssaflowline -snes_monitor -snes_fd\n\
*  use exact solution as initial guess:\n\
      ./ssaflowline -snes_monitor -ssa_guess 2\n\
*  ask PETSc to check our analytical Jacobian in M=6 case (ignore all but the\n\
          first displayed case):\n\
      ./ssaflowline -da_grid_x 6 -mat_fd_type ds -snes_type test\n\
*  use X graphics to show velocity during iteration:\n\
      ./ssaflowline -snes_monitor_solution -draw_pause 1\n\
*  use X graphics to show thickness and exact solution:\n\
      ./ssaflowline -ssa_show -draw_pause 2\n\
*  extreme accuracy on 10 cm grid:\n\
      ./ssaflowline -snes_monitor -da_grid_x 2000001 -ssa_epsilon 1.0e-8\n\
*  major idea: compare performance using analytical Jacobian, straight Picard\n\
          iteration, and matrix-free with Picard-as-preconditioner\n\
          (all on 1 m grid):\n\
      time ./ssaflowline -snes_monitor -da_grid_x 200001\n\
      time ./ssaflowline -snes_monitor -da_grid_x 200001 -ssa_picard \n\
      time ./ssaflowline -snes_monitor -da_grid_x 200001 -ssa_picard -snes_mf_operator \n\
\n";

#include "petscdmda.h"
#include "petscsnes.h"
#include "petscdraw.h"
#include "petscviewer.h"


/* User-defined application context - contains data needed by the 
   application-provided call-back routines, esp. FormFunctionLocal().  */
typedef struct {
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
  PetscReal   epsilon;
  PetscReal   Hcalv;
  PetscReal   uexactcalv;
  PetscReal   gamma;
  Vec         H;
  Vec         uexact;
  Vec         viscosity;
  Mat         J;
} AppCtx;


/* declare the user-written routines ... */
extern PetscErrorCode FillThicknessAndExactSoln(DM,AppCtx*);
extern PetscErrorCode FormInitialGuess(DM,AppCtx*,Vec);

//static PetscErrorCode FormFunctionLocal(DMDALocalInfo*,PetscScalar*,PetscScalar*,AppCtx*);
extern PetscErrorCode FormFunctionLocal(DMDALocalInfo*,PetscScalar*,PetscScalar*,AppCtx*);

/* ... including two versions of "analytical" Jacobian; these functions are only
   evaluated if *not* using one of the methods -fd, -snes_fd, -snes_mf */
extern PetscErrorCode FormTrueJacobianMatrixLocal(DMDALocalInfo*,PetscScalar*,Mat,Mat,AppCtx*);
extern PetscErrorCode FormPicardMatrixLocal(DMDALocalInfo*,PetscScalar*,Mat,Mat,AppCtx*);


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  PetscErrorCode         ierr;

  DM                     da;
  SNES                   snes;                 /* nonlinear solver */
  Vec                    u,r;                  /* solution, residual vectors */
  AppCtx                 user;                 /* user-defined work context */
  PetscInt               guess=1,              /* initial guess type */
                         its,Mx;               /* iteration count, num of pts */
  PetscReal              err1,errinf;          /* average and max norm of error */
  PetscBool              checks = PETSC_FALSE, /* display values at calving front */
                         show = PETSC_FALSE,   /* show some fields in X viewers;
                                                  use -draw_pause N for N sec delay */
                         matview = PETSC_FALSE;/* dump preconditioner matrix */
  SNESConvergedReason    reason;               /* Check convergence */
  PetscBool              picard = PETSC_FALSE,
                         eps_set = PETSC_FALSE;

  PetscInitialize(&argc,&argv,(char *)0,help);

  ierr = PetscPrintf(PETSC_COMM_WORLD,
    "SSAFLOWLINE solves for velocity in 1D ice shelf with steady-state thickness\n"
    "  [run with -help for info and options]\n");CHKERRQ(ierr);

  user.n       = 3.0;          /* Glen flow law exponent */
  user.secpera = 31556926.0;
  user.rho     = 900.0;        /* kg m^-3 */
  user.rhow    = 1000.0;       /* kg m^-3 */
  user.g       = 9.8;          /* m s^-2 */
  user.A       = 1.4579e-25;   /* s^-1 Pa^-3; used by MacAyeal et al (1996) */
  user.ug      = 50.0 / user.secpera; /* m s^-1 */
  user.Hg      = 500.0;        /* m */
  user.accum   = 0.3 / user.secpera; /* m s^-1 */
  user.L       = 200.0e3;      /* m */
  user.epsilon = (1.0 / user.secpera) / user.L; /* regularize using strain rate
                                                   of 1/L = 5e-6 per year */
  /* these are set in FillThicknessAndExactSoln(): */
  user.Hcalv      = -1.0;
  user.uexactcalv = -1.0;
  user.gamma      = -1.0;
  /* this matrix points to the preconditioner if user does -ssa_mat_view */
  user.J = PETSC_NULL;

  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,NULL,
                           "SSA flow line options",__FILE__);CHKERRQ(ierr);
  {
    ierr = PetscOptionsReal("-ssa_glen","Glen flow law exponent for SSA","",
                            user.n,&user.n,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-ssa_rho","ice density (kg m^-3) for SSA","",
                            user.rho,&user.rho,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-ssa_rhow","sea water density (kg m^-3) for SSA","",
                            user.rhow,&user.rhow,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-ssa_softness","A = ice softness (s^-1 Pa^-3) for SSA","",
                            user.A,&user.A,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-ssa_ug","velocity across grounding line (m s^-1) for SSA","",
                            user.ug,&user.ug,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-ssa_Hg","thickness at grounding line (m) for SSA","",
                            user.Hg,&user.Hg,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-ssa_accum","constant accumulation rate over shelf for SSA","",
                            user.accum,&user.accum,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-ssa_L","ice shelf length for SSA","",
                            user.L,&user.L,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-ssa_epsilon","regularization (a strain rate in units of 1/a)","",
                            user.epsilon * user.secpera,&user.epsilon,&eps_set);CHKERRQ(ierr);
    if (eps_set) {  user.epsilon *= 1.0 / user.secpera;  }
    ierr = PetscOptionsBool("-ssa_picard","compute Picard matrix instead of analytical Jacobian","",
                             picard,&picard,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsBool("-ssa_checks","print values at calving front as minimal check","",
                             checks,&checks,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsBool("-ssa_show","show thickness, exact viscosity, and exact solution","",
                             show,&show,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsBool("-ssa_mat_view","put preconditioner matrix in Matlab form to stdout","",
                             matview,&matview,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsInt( "-ssa_guess","initial guess: 1=(linear; default), 2=(exact soln)","",
                             guess,&guess,NULL);CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);

  /* check some user parameters for reasonableness ... there could be more ... */
  if (user.n < 1.0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "SSA WARNING: glen exponent n=%g is out of range (require n>= 1)\n",
                       user.n);CHKERRQ(ierr);
  }
  if ((guess != 1) && (guess != 2)) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "SSA ERROR: -ssa_guess must be 1 or 2 ... ending\n");CHKERRQ(ierr);
    PetscEnd();
  }

  /* Create machinery for parallel grid management (DMDA), nonlinear solver (SNES), 
     and Vecs for fields (thickness, velocity, RHS).  Note default Mx=21 is
     number of grid points.  Also degrees of freedom = 1 (scalar problem) and
     stencil radius = ghost width = 1.                                    */
  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,
                      -21,1,1,PETSC_NULL,&da);CHKERRQ(ierr);
  ierr = DMDASetUniformCoordinates(da,0.0,user.L, 0.0, 1.0, 0.0, 1.0);CHKERRQ(ierr);
  ierr = DMSetApplicationContext(da,&user);CHKERRQ(ierr);

  /* Extract global vectors from DA and duplicate (allocate) for remaining same types */
  ierr = DMCreateGlobalVector(da,&u);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&r);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&user.H);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&user.uexact);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&user.viscosity);CHKERRQ(ierr);

  ierr = DMDASNESSetFunctionLocal(da,INSERT_VALUES,(DMDASNESFunction)FormFunctionLocal,&user);CHKERRQ(ierr);
  if (picard) {
    ierr = DMDASNESSetJacobianLocal(da,(DMDASNESJacobian)FormPicardMatrixLocal,&user);CHKERRQ(ierr);
  } else {
    ierr = DMDASNESSetJacobianLocal(da,(DMDASNESJacobian)FormTrueJacobianMatrixLocal,&user);CHKERRQ(ierr);
  }

  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
  ierr = SNESSetDM(snes,da);CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  ierr = DMDAGetInfo(da,PETSC_IGNORE,&Mx,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE);
                     CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"  grid of Mx=%d points and dx=%.3f m\n",
                     Mx,user.L/(Mx-1));CHKERRQ(ierr);

  /* We use a particular formula for the thickness (user.H).  This formula has the
     property that our numerical velocity solution (u) should converge to a
     known exact solution (user.uexact), which is also computed here.  Note that the 
     thickness is computed on the staggered grid and the exact solution on the
     regular grid.  For more, see notes lecture.pdf.   */
  ierr = FillThicknessAndExactSoln(da,&user);CHKERRQ(ierr);

  if (checks || show) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,  /* minimal checks on thickness, exact */
           "  ssaflowline checks:\n"
           "    user.Hcalv      = %.4f (?= 279.7397 m, if L=200 km [default])\n"
           "    user.uexactcalv = %.4f (?= 303.8539 m/a, if L=200 km [default])\n",
           user.Hcalv, user.uexactcalv * user.secpera);CHKERRQ(ierr);
  }  
  if (show) {
    PetscDraw draw;
    ierr = PetscViewerDrawGetDraw(PETSC_VIEWER_DRAW_WORLD,0,&draw); CHKERRQ(ierr);

    ierr = PetscDrawSetTitle(draw,"ice thickness (m)"); CHKERRQ(ierr);
    ierr = VecView(user.H,PETSC_VIEWER_DRAW_WORLD); CHKERRQ(ierr);

    ierr = PetscDrawSetTitle(draw,"exact viscosity (10^15 Pa s)"); CHKERRQ(ierr);
    ierr = VecScale(user.viscosity,1.0e-15); CHKERRQ(ierr);
    ierr = VecView(user.viscosity,PETSC_VIEWER_DRAW_WORLD); CHKERRQ(ierr);
    ierr = VecScale(user.viscosity,1.0e-15); CHKERRQ(ierr);

    ierr = PetscDrawSetTitle(draw,"exact ice velocity (m/a)"); CHKERRQ(ierr);
    ierr = VecScale(user.uexact,user.secpera); CHKERRQ(ierr);
    ierr = VecView(user.uexact,PETSC_VIEWER_DRAW_WORLD); CHKERRQ(ierr);
    ierr = VecScale(user.uexact,1.0 / user.secpera); CHKERRQ(ierr);

    ierr = PetscDrawSetTitle(draw,""); CHKERRQ(ierr);  // clear title
  }

  /* Evaluate initial guess.  The user needs to initialize the solution vector
     (x) with the initial guess prior to calling SNESSolve().   */
  if (guess == 1) { /* linear guess: see notes */
    ierr = PetscPrintf(PETSC_COMM_WORLD,"  using linear initial guess\n");
             CHKERRQ(ierr);
    ierr = FormInitialGuess(da,&user,u);CHKERRQ(ierr);
  } else if (guess == 2) { /* use exact soln as initial guess: */
    ierr = PetscPrintf(PETSC_COMM_WORLD,"  using exact solution as initial guess\n");
             CHKERRQ(ierr);
    ierr = VecCopy(user.uexact,u);CHKERRQ(ierr);
  } else {
    SETERRQ(PETSC_COMM_SELF,1,"invalid integer for guess = (initial solution guess case)\n");
  }

  /************ SOLVE NONLINEAR SYSTEM!  ************/
  ierr = SNESSolve(snes,PETSC_NULL,u);CHKERRQ(ierr);

  ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr);
  ierr = SNESGetConvergedReason(snes,&reason);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
           "  %s number of Newton iterations = %D\n",
           SNESConvergedReasons[reason],its);CHKERRQ(ierr);

  /* evaluate error relative to exact solution */
  ierr = VecAXPY(u,-1.0,user.uexact);CHKERRQ(ierr); /* "y:=ax+y"  so   u := u - uexact */
  ierr = VecNorm(u,NORM_1,&err1);CHKERRQ(ierr);
  ierr = VecNorm(u,NORM_INFINITY,&errinf);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
           "  numerical errors in velocity: %.4e m/a average\n"
           "                                %.4e m/a maximum\n"
           "                                %.4e relative maximum\n",
           err1*user.secpera/Mx,errinf*user.secpera,errinf/user.uexactcalv);CHKERRQ(ierr);
  
  ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&r);CHKERRQ(ierr);
  ierr = VecDestroy(&(user.H));CHKERRQ(ierr);
  ierr = VecDestroy(&(user.uexact));CHKERRQ(ierr);
  ierr = VecDestroy(&(user.viscosity));CHKERRQ(ierr);

  ierr = SNESDestroy(&snes);CHKERRQ(ierr);
  ierr = DMDestroy(&da);CHKERRQ(ierr);

  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "GetUEx"
static inline void GetUEx(PetscScalar ug, PetscScalar qg, PetscScalar a, PetscScalar n, 
                          PetscScalar Cs, PetscScalar flux,
                          PetscScalar *UEx, PetscScalar *dUdxEx) {
  const PetscScalar inside = PetscPowScalar(ug,n+1.0) + (Cs / a)
                * (PetscPowScalar(flux,n+1.0) - PetscPowScalar(qg,n+1.0));
  *UEx    = PetscPowScalar(inside, 1.0 / (n+1.0));
  *dUdxEx = Cs * PetscPowScalar(flux, n) * PetscPowScalar(inside, (1.0 / (n+1.0)) - 1.0);
}


#undef __FUNCT__
#define __FUNCT__ "GetViscosityFromStrainRate"
static inline PetscScalar GetViscosityFromStrainRate(PetscScalar dudx, PetscScalar p, 
                                                     PetscScalar B) {
  return  B * PetscPowScalar(dudx * dudx, (p - 2.0) / 2.0 );
}


#undef __FUNCT__
#define __FUNCT__ "FillThicknessAndExactSoln"
/*  Compute the right thickness H = H(x).  See lecture notes for exact solution. */
PetscErrorCode FillThicknessAndExactSoln(DM da, AppCtx *user)
{
  PetscErrorCode ierr;
  PetscInt       i,Mx,xs,xm;
  PetscReal      hx, n, r, Cs, xx, flux, qg, p, B, dudx, ustag;
  PetscScalar    *H, *uex, *visc;

  PetscFunctionBegin;
  ierr = DMDAGetInfo(da,PETSC_IGNORE,&Mx,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);
  ierr = DMDAGetCorners(da,&xs,PETSC_NULL,PETSC_NULL,&xm,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);

  /* constants, independent of x */
  hx = user->L / (PetscReal)(Mx-1);
  n  = user->n;
  r  = user->rho / user->rhow;
  Cs = user->A * PetscPowScalar( ( 0.25 * user->rho * user->g * (1.0 - r) ), n);
  qg = user->ug * user->Hg;
  p  = 1.0 + 1.0 / user->n;
  B  = PetscPowScalar(user->A,-1.0/user->n);

  /* Compute regular grid exact soln and staggered-grid thickness over the
     locally-owned part of the grid */
  ierr = DMDAVecGetArray(da,user->uexact,&uex);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,user->H,&H);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,user->viscosity,&visc);CHKERRQ(ierr);
  for (i=xs; i<xs+xm; i++) {
    /* get exact velocity and strain rate on regular grid */
    xx = hx * (PetscReal)i;  /* = x_i = distance from grounding line */
    flux = user->accum * xx + qg; /* flux at x_i */
    GetUEx(user->ug, qg, user->accum, n, Cs, flux, &(uex[i]), &dudx);

    /* exact viscosity on regular grid */
    visc[i] = GetViscosityFromStrainRate(dudx, p, B);

    /* exact thickness on staggered grid */
    flux += user->accum * hx * 0.5; /* flux at x_{i+1/2} */
    GetUEx(user->ug, qg, user->accum, n, Cs, flux, &ustag, &dudx);
    H[i] = flux / ustag;
  }
  ierr = DMDAVecRestoreArray(da,user->uexact,&uex);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da,user->H,&H);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da,user->viscosity,&visc);CHKERRQ(ierr);

  /* separately compute and store calving-front values */
  flux = user->accum * user->L + qg;
  GetUEx(user->ug, qg, user->accum, n, Cs, flux, &(user->uexactcalv), &dudx);
  user->Hcalv = flux / user->uexactcalv;

  /* strain rate at calving front */
  /* MATLAB: gamma = ( 0.25 * p.A^(1/n) * (1 - r) * p.rho * p.g * H(end) )^n; */
  user->gamma = 0.25 * PetscPowScalar(user->A,1.0/user->n) * (1.0 - r)
                  * user->rho * user->g * user->Hcalv;
  user->gamma = PetscPowScalar(user->gamma,user->n);

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FormInitialGuess"
PetscErrorCode FormInitialGuess(DM da, AppCtx *user,Vec X)
{
  PetscErrorCode ierr;
  PetscInt       i,Mx,xs,xm;
  PetscReal      hx, xx;
  PetscScalar    *u;

  PetscFunctionBegin;
  ierr = DMDAGetInfo(da,PETSC_IGNORE,&Mx,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);
  ierr = DMDAGetCorners(da,&xs,PETSC_NULL,PETSC_NULL,&xm,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);

  hx = user->L / (PetscReal)(Mx-1);

  /* Compute initial guess over the locally owned part of the grid.  See lecture
     notes for comments on why this linear function is a reasonable initial guess. */
  ierr = DMDAVecGetArray(da,X,&u);CHKERRQ(ierr);
  for (i=xs; i<xs+xm; i++) {
    xx = hx * (PetscReal)i;  /* distance from grounding line */
    u[i] = user->ug + user->gamma * xx;
  }
  ierr = DMDAVecRestoreArray(da,X,&u);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "GetEta"
static inline PetscScalar GetEta(PetscScalar Z, PetscScalar dx, PetscScalar p, PetscScalar eps) {
  return  PetscPowScalar((Z/dx) * (Z/dx) + eps * eps, (p - 2.0) / 2.0 );
}


#undef __FUNCT__
#define __FUNCT__ "FormFunctionLocal"
/* Compute pointwise residual f(x) over the locally-owned part of the grid
   This is a finite difference method. In TeX, the formula we compute is

f_i = \eta_{i+1/2} H_{i+1/2} (u_{i+1}-u_i) - \eta_{i-1/2} H_{i-1/2} (u_i-u_{i-1})
      - dx K (H_{i+1/2}^2 - H_{i-1/2}^2)

where

\eta_{i+1/2} = \left|\frac{u_{i+1}-u_i}{dx}\right|^{p-2}

with some regularization using user.epsilon, and

  dx = L / Mx
  p = 1+1/n
  K = rho * g * (1-r) / (4 * B)
  r = rho / rhow
  B = A^{1/n}

*/
PetscErrorCode FormFunctionLocal(DMDALocalInfo *info,PetscScalar *u,PetscScalar *f,AppCtx *user)
{
  PetscErrorCode ierr;
  PetscReal      hx, p, K, B, duL, duR, sL, sR;
  PetscScalar    *H;
  PetscInt       i, Mx;
  Vec            localH;

  PetscFunctionBegin;

  ierr = DMGetLocalVector(info->da,&localH);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(info->da,user->H,INSERT_VALUES,localH); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(info->da,user->H,INSERT_VALUES,localH); CHKERRQ(ierr);

  p  = 1.0 + 1.0 / user->n;
  B  = PetscPowScalar(user->A,-1.0/user->n);
  K  = user->rho * user->g * (1.0 - user->rho/user->rhow) / (4.0 * B);
  Mx = info->mx;
  hx = user->L / ((PetscReal)Mx - 1.0);

  ierr = DMDAVecGetArray(info->da,localH,&H);CHKERRQ(ierr);
  for (i=info->xs; i<info->xs+info->xm; i++) {
    if (i == 0) {
      f[0] = u[0] - user->ug;  /* Dirichlet condition */
    } else {
      if (i == 1) { /* use Dirichlet condition as value for neighbor, which symmetrizes */
        duL = u[i] - user->ug;
      } else {
        duL = u[i] - u[i-1];
      }
      sL = GetEta(duL, hx, p, user->epsilon) * H[i-1] * duL;
      if (i == Mx-1) {  /* Neumann: calving front stress boundary condition */
        duR = u[Mx-2] + 2.0 * hx * user->gamma - u[Mx-1];
      } else {
        duR = u[i+1] - u[i];
      }
      sR = GetEta(duR, hx, p, user->epsilon) * H[i] * duR;
      f[i] = sR - sL - hx * K * (H[i]*H[i] - H[i-1]*H[i-1]);
    }
  }
  ierr = DMDAVecRestoreArray(info->da,localH,&H);CHKERRQ(ierr);

  ierr = DMRestoreLocalVector(info->da,&localH);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FormPicardMatrixLocal"
/* FormPicardMatrixLocal - Evaluates an approximation to the Jacobian matrix. 
The result of this routine can be used as a preconditioner to a finite-
difference matrix-free approach *or* directly as an approximation to the Jacobian.
The former happens with options
  -ssa_picard -snes_mf_operator */
PetscErrorCode FormPicardMatrixLocal(DMDALocalInfo *info,PetscScalar *u,Mat picpre,Mat pic,AppCtx *user)
{
  PetscErrorCode ierr;
  PetscInt       i, Mx;
  PetscInt       col[3],row[1];
  PetscScalar    v[3];
  PetscReal      hx, p, duL, duR, sLfact, sRfact;
  PetscScalar    *H;
  Vec            localH;

  PetscFunctionBegin;

  ierr = DMGetLocalVector(info->da,&localH);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(info->da,user->H,INSERT_VALUES,localH); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(info->da,user->H,INSERT_VALUES,localH); CHKERRQ(ierr);

  p  = 1.0 + 1.0 / user->n;
  Mx = info->mx;
  hx = user->L / ((PetscReal)Mx - 1.0);

  ierr = DMDAVecGetArray(info->da,localH,&H);CHKERRQ(ierr);
  for (i=info->xs; i<info->xs+info->xm; i++) {
    row[0] = i;
    if (i == 0) {
      col[0] = 0;
      v[0] = 1.0;
      ierr = MatSetValues(pic,1,row,1,col,v,INSERT_VALUES);CHKERRQ(ierr);
    } else {
      if (i == 1) {
        duL = u[i] - user->ug;
      } else {
        duL = u[i] - u[i-1];
      }
      sLfact = GetEta(duL, hx, p, user->epsilon) * H[i-1];
      if (i == Mx-1) {  /* Neumann: calving front stress boundary condition */
        duR = u[Mx-2] + 2.0 * hx * user->gamma - u[Mx-1];
      } else {
        duR = u[i+1] - u[i];
      }
      sRfact = GetEta(duR, hx, p, user->epsilon) * H[i];
      
      if (i == 1) {
        col[0] = 1;  col[1] = 2;
        v[0] = -sLfact - sRfact;  v[1] = sRfact;
        ierr = MatSetValues(pic,1,row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
      } else if (i == Mx-1) {
        col[0] = i-1;  col[1] = i;
        v[0] = sLfact + sRfact; v[1] = -sLfact - sRfact;
        ierr = MatSetValues(pic,1,row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
      } else {
        col[0] = i-1;  col[1] = i; col[2] = i+1;
        v[0] = sLfact; v[1] = -sLfact - sRfact;  v[2] = sRfact;
        ierr = MatSetValues(pic,1,row,3,col,v,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
  }
  ierr = DMDAVecRestoreArray(info->da,localH,&H);CHKERRQ(ierr);

  ierr = DMRestoreLocalVector(info->da,&localH);CHKERRQ(ierr);

  /* Assemble matrix, using the 2-step process */
  ierr = MatAssemblyBegin(pic,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(pic,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  /* Tell the matrix we will never add a new nonzero location to the
     matrix. If we do, it will generate an error.                    */
  ierr = MatSetOption(pic,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);

  if (user->J != PETSC_NULL) {
    ierr = PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) user->J,"Picard_matrix"); CHKERRQ(ierr);
    ierr = MatView(user->J, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "GetDEta"
static inline PetscScalar GetDEta(PetscScalar Z, PetscScalar dx, PetscScalar p, PetscScalar eps) {
  return  ((p - 2.0) / (dx * dx)) * Z * 
              PetscPowScalar((Z/dx) * (Z/dx) + eps * eps, (p - 4.0) / 2.0 );
}


#undef __FUNCT__
#define __FUNCT__ "GetOmega"
static inline PetscScalar GetOmega(PetscScalar Z, PetscScalar dx, PetscScalar p, PetscScalar eps) {
  return  Z * GetDEta(Z,dx,p,eps) + GetEta(Z,dx,p,eps);
}


#undef __FUNCT__
#define __FUNCT__ "FormTrueJacobianMatrixLocal"
/* FormTrueJacobianMatrixLocal - Evaluates analytical Jacobian matrix. */
PetscErrorCode FormTrueJacobianMatrixLocal(DMDALocalInfo *info,PetscScalar *u,Mat jacpre,Mat jac,AppCtx *user)
{
  PetscErrorCode ierr;
  PetscInt       i, Mx;
  PetscInt       col[3],row[1];
  PetscScalar    v[3];
  PetscReal      hx, p, duL, duR, omHL, omHR;
  PetscScalar    *H;
  Vec            localH;

  PetscFunctionBegin;

  ierr = DMGetLocalVector(info->da,&localH);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(info->da,user->H,INSERT_VALUES,localH); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(info->da,user->H,INSERT_VALUES,localH); CHKERRQ(ierr);

  p  = 1.0 + 1.0 / user->n;
  Mx = info->mx;
  hx = user->L / ((PetscReal)Mx - 1.0);

  ierr = DMDAVecGetArray(info->da,localH,&H);CHKERRQ(ierr);
  for (i=info->xs; i<info->xs+info->xm; i++) {
    row[0] = i;
    if (i == 0) {
      col[0] = 0;
      v[0] = 1.0;
      ierr = MatSetValues(jac,1,row,1,col,v,INSERT_VALUES);CHKERRQ(ierr);
    } else {

      if (i == 1) {
        duL = u[i] - user->ug;
      } else {
        duL = u[i] - u[i-1];
      }
      omHL = GetOmega(duL, hx, p, user->epsilon) * H[i-1];
      if (i == Mx-1) {  /* Neumann: calving front stress boundary condition */
        duR = u[Mx-2] + 2.0 * hx * user->gamma - u[Mx-1];
      } else {
        duR = u[i+1] - u[i];
      }
      omHR = GetOmega(duR, hx, p, user->epsilon) * H[i];
      
      if (i == 1) {
        col[0] = 1;  col[1] = 2;
        v[0] = - omHL - omHR;  v[1] = omHR;
        ierr = MatSetValues(jac,1,row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
      } else if (i == Mx-1) {
        col[0] = i-1;  col[1] = i;
        v[0] = omHL + omHR; v[1] = - omHL - omHR;
        ierr = MatSetValues(jac,1,row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
      } else {
        col[0] = i-1;  col[1] = i; col[2] = i+1;
        v[0] = omHL; v[1] = - omHL - omHR;  v[2] = omHR;
        ierr = MatSetValues(jac,1,row,3,col,v,INSERT_VALUES);CHKERRQ(ierr);
      }
    }

  }
  ierr = DMDAVecRestoreArray(info->da,localH,&H);CHKERRQ(ierr);

  ierr = DMRestoreLocalVector(info->da,&localH);CHKERRQ(ierr);

  /* Assemble matrix, using the 2-step process */
  ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  /* Tell the matrix we will never add a new nonzero location to the
     matrix. If we do, it will generate an error.                    */
  ierr = MatSetOption(jac,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);

  if (user->J != PETSC_NULL) {
    ierr = PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) user->J,"Jacobian_matrix"); CHKERRQ(ierr);
    ierr = MatView(user->J, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}
