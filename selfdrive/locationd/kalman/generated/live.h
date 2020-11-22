/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_3158214232202474891);
void inv_err_fun(double *nom_x, double *true_x, double *out_5725068012840555629);
void H_mod_fun(double *state, double *out_8657057670133681652);
void f_fun(double *state, double dt, double *out_8736522924417776153);
void F_fun(double *state, double dt, double *out_3328850354945798330);
void h_3(double *state, double *unused, double *out_7720750652569981479);
void H_3(double *state, double *unused, double *out_341792269985784737);
void h_4(double *state, double *unused, double *out_8142258943073980802);
void H_4(double *state, double *unused, double *out_3240307902281911554);
void h_9(double *state, double *unused, double *out_3917925883842028233);
void H_9(double *state, double *unused, double *out_3817800836669103682);
void h_10(double *state, double *unused, double *out_7971293568959411420);
void H_10(double *state, double *unused, double *out_8404913797888176650);
void h_12(double *state, double *unused, double *out_5135172896049446395);
void H_12(double *state, double *unused, double *out_6832469082739381262);
void h_13(double *state, double *unused, double *out_1126848983439859532);
void H_13(double *state, double *unused, double *out_7847509844363063938);
void h_14(double *state, double *unused, double *out_3917925883842028233);
void H_14(double *state, double *unused, double *out_3817800836669103682);
void h_19(double *state, double *unused, double *out_7576122994831089039);
void H_19(double *state, double *unused, double *out_954613562813734758);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);