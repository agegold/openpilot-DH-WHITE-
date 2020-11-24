/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_6527096871983429504);
void inv_err_fun(double *nom_x, double *true_x, double *out_7406771233427917374);
void H_mod_fun(double *state, double *out_8208304452669981300);
void f_fun(double *state, double dt, double *out_5789550752280723550);
void F_fun(double *state, double dt, double *out_7618743868083865277);
void h_3(double *state, double *unused, double *out_4798363570708426436);
void H_3(double *state, double *unused, double *out_8457878772352757720);
void h_4(double *state, double *unused, double *out_8310912896270388338);
void H_4(double *state, double *unused, double *out_5315640738519338261);
void h_9(double *state, double *unused, double *out_4314120975370170928);
void H_9(double *state, double *unused, double *out_3762568502684500560);
void h_10(double *state, double *unused, double *out_2412124971989410055);
void H_10(double *state, double *unused, double *out_6410014017175275927);
void h_12(double *state, double *unused, double *out_390967357260881350);
void H_12(double *state, double *unused, double *out_7422196250258547830);
void h_31(double *state, double *unused, double *out_846408729819177188);
void H_31(double *state, double *unused, double *out_4747106475115128009);
void h_32(double *state, double *unused, double *out_2343051480456404082);
void H_32(double *state, double *unused, double *out_1553259269145273177);
void h_13(double *state, double *unused, double *out_4182319111389486617);
void H_13(double *state, double *unused, double *out_653447131660764693);
void h_14(double *state, double *unused, double *out_4314120975370170928);
void H_14(double *state, double *unused, double *out_3762568502684500560);
void h_19(double *state, double *unused, double *out_7879957134962037784);
void H_19(double *state, double *unused, double *out_670111348417930621);
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
const static double MAHA_THRESH_31 = 7.814728;
void update_31(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_32 = 9.487729;
void update_32(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);