/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_5456919468856549453);
void inv_err_fun(double *nom_x, double *true_x, double *out_5406950546569109440);
void H_mod_fun(double *state, double *out_8055124209405881978);
void f_fun(double *state, double dt, double *out_7061105761666017332);
void F_fun(double *state, double dt, double *out_5005698098331862813);
void h_25(double *state, double *unused, double *out_1051757273230348996);
void H_25(double *state, double *unused, double *out_6043417528456999617);
void h_24(double *state, double *unused, double *out_5022686355337100841);
void H_24(double *state, double *unused, double *out_5973588777829285175);
void h_30(double *state, double *unused, double *out_5719077953120001940);
void H_30(double *state, double *unused, double *out_3570481382492183663);
void h_26(double *state, double *unused, double *out_7725992852928704031);
void H_26(double *state, double *unused, double *out_6698733390313342726);
void h_27(double *state, double *unused, double *out_108247420328686561);
void H_27(double *state, double *unused, double *out_4858063370328808975);
void h_29(double *state, double *unused, double *out_7078705484129371768);
void H_29(double *state, double *unused, double *out_2995922768636647488);
void h_28(double *state, double *unused, double *out_3445327860317372799);
void H_28(double *state, double *unused, double *out_5455322817643771575);
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
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
