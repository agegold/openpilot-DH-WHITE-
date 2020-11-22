/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_6182305149975481084);
void inv_err_fun(double *nom_x, double *true_x, double *out_5772716229829584760);
void H_mod_fun(double *state, double *out_6573258343469333187);
void f_fun(double *state, double dt, double *out_9174724217193728746);
void F_fun(double *state, double dt, double *out_7673946190717308446);
void h_25(double *state, double *unused, double *out_628975419662729201);
void H_25(double *state, double *unused, double *out_7448964055686931099);
void h_24(double *state, double *unused, double *out_8549256257243174563);
void H_24(double *state, double *unused, double *out_7925276905366183369);
void h_26(double *state, double *unused, double *out_1612116933575715082);
void H_26(double *state, double *unused, double *out_2926951939376508717);
void h_27(double *state, double *unused, double *out_356100148548021822);
void H_27(double *state, double *unused, double *out_6872708239176032209);
void h_29(double *state, double *unused, double *out_6898344091781638272);
void H_29(double *state, double *unused, double *out_5593307003910287637);
void h_28(double *state, double *unused, double *out_6095167101378246893);
void H_28(double *state, double *unused, double *out_3375372154212119283);
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
