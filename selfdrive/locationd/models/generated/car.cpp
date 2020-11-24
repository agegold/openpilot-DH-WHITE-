
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_5456919468856549453) {
   out_5456919468856549453[0] = delta_x[0] + nom_x[0];
   out_5456919468856549453[1] = delta_x[1] + nom_x[1];
   out_5456919468856549453[2] = delta_x[2] + nom_x[2];
   out_5456919468856549453[3] = delta_x[3] + nom_x[3];
   out_5456919468856549453[4] = delta_x[4] + nom_x[4];
   out_5456919468856549453[5] = delta_x[5] + nom_x[5];
   out_5456919468856549453[6] = delta_x[6] + nom_x[6];
   out_5456919468856549453[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_5406950546569109440) {
   out_5406950546569109440[0] = -nom_x[0] + true_x[0];
   out_5406950546569109440[1] = -nom_x[1] + true_x[1];
   out_5406950546569109440[2] = -nom_x[2] + true_x[2];
   out_5406950546569109440[3] = -nom_x[3] + true_x[3];
   out_5406950546569109440[4] = -nom_x[4] + true_x[4];
   out_5406950546569109440[5] = -nom_x[5] + true_x[5];
   out_5406950546569109440[6] = -nom_x[6] + true_x[6];
   out_5406950546569109440[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_8055124209405881978) {
   out_8055124209405881978[0] = 1.0;
   out_8055124209405881978[1] = 0.0;
   out_8055124209405881978[2] = 0.0;
   out_8055124209405881978[3] = 0.0;
   out_8055124209405881978[4] = 0.0;
   out_8055124209405881978[5] = 0.0;
   out_8055124209405881978[6] = 0.0;
   out_8055124209405881978[7] = 0.0;
   out_8055124209405881978[8] = 0.0;
   out_8055124209405881978[9] = 1.0;
   out_8055124209405881978[10] = 0.0;
   out_8055124209405881978[11] = 0.0;
   out_8055124209405881978[12] = 0.0;
   out_8055124209405881978[13] = 0.0;
   out_8055124209405881978[14] = 0.0;
   out_8055124209405881978[15] = 0.0;
   out_8055124209405881978[16] = 0.0;
   out_8055124209405881978[17] = 0.0;
   out_8055124209405881978[18] = 1.0;
   out_8055124209405881978[19] = 0.0;
   out_8055124209405881978[20] = 0.0;
   out_8055124209405881978[21] = 0.0;
   out_8055124209405881978[22] = 0.0;
   out_8055124209405881978[23] = 0.0;
   out_8055124209405881978[24] = 0.0;
   out_8055124209405881978[25] = 0.0;
   out_8055124209405881978[26] = 0.0;
   out_8055124209405881978[27] = 1.0;
   out_8055124209405881978[28] = 0.0;
   out_8055124209405881978[29] = 0.0;
   out_8055124209405881978[30] = 0.0;
   out_8055124209405881978[31] = 0.0;
   out_8055124209405881978[32] = 0.0;
   out_8055124209405881978[33] = 0.0;
   out_8055124209405881978[34] = 0.0;
   out_8055124209405881978[35] = 0.0;
   out_8055124209405881978[36] = 1.0;
   out_8055124209405881978[37] = 0.0;
   out_8055124209405881978[38] = 0.0;
   out_8055124209405881978[39] = 0.0;
   out_8055124209405881978[40] = 0.0;
   out_8055124209405881978[41] = 0.0;
   out_8055124209405881978[42] = 0.0;
   out_8055124209405881978[43] = 0.0;
   out_8055124209405881978[44] = 0.0;
   out_8055124209405881978[45] = 1.0;
   out_8055124209405881978[46] = 0.0;
   out_8055124209405881978[47] = 0.0;
   out_8055124209405881978[48] = 0.0;
   out_8055124209405881978[49] = 0.0;
   out_8055124209405881978[50] = 0.0;
   out_8055124209405881978[51] = 0.0;
   out_8055124209405881978[52] = 0.0;
   out_8055124209405881978[53] = 0.0;
   out_8055124209405881978[54] = 1.0;
   out_8055124209405881978[55] = 0.0;
   out_8055124209405881978[56] = 0.0;
   out_8055124209405881978[57] = 0.0;
   out_8055124209405881978[58] = 0.0;
   out_8055124209405881978[59] = 0.0;
   out_8055124209405881978[60] = 0.0;
   out_8055124209405881978[61] = 0.0;
   out_8055124209405881978[62] = 0.0;
   out_8055124209405881978[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_7061105761666017332) {
   out_7061105761666017332[0] = state[0];
   out_7061105761666017332[1] = state[1];
   out_7061105761666017332[2] = state[2];
   out_7061105761666017332[3] = state[3];
   out_7061105761666017332[4] = state[4];
   out_7061105761666017332[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_7061105761666017332[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_7061105761666017332[7] = state[7];
}
void F_fun(double *state, double dt, double *out_5005698098331862813) {
   out_5005698098331862813[0] = 1;
   out_5005698098331862813[1] = 0;
   out_5005698098331862813[2] = 0;
   out_5005698098331862813[3] = 0;
   out_5005698098331862813[4] = 0;
   out_5005698098331862813[5] = 0;
   out_5005698098331862813[6] = 0;
   out_5005698098331862813[7] = 0;
   out_5005698098331862813[8] = 0;
   out_5005698098331862813[9] = 1;
   out_5005698098331862813[10] = 0;
   out_5005698098331862813[11] = 0;
   out_5005698098331862813[12] = 0;
   out_5005698098331862813[13] = 0;
   out_5005698098331862813[14] = 0;
   out_5005698098331862813[15] = 0;
   out_5005698098331862813[16] = 0;
   out_5005698098331862813[17] = 0;
   out_5005698098331862813[18] = 1;
   out_5005698098331862813[19] = 0;
   out_5005698098331862813[20] = 0;
   out_5005698098331862813[21] = 0;
   out_5005698098331862813[22] = 0;
   out_5005698098331862813[23] = 0;
   out_5005698098331862813[24] = 0;
   out_5005698098331862813[25] = 0;
   out_5005698098331862813[26] = 0;
   out_5005698098331862813[27] = 1;
   out_5005698098331862813[28] = 0;
   out_5005698098331862813[29] = 0;
   out_5005698098331862813[30] = 0;
   out_5005698098331862813[31] = 0;
   out_5005698098331862813[32] = 0;
   out_5005698098331862813[33] = 0;
   out_5005698098331862813[34] = 0;
   out_5005698098331862813[35] = 0;
   out_5005698098331862813[36] = 1;
   out_5005698098331862813[37] = 0;
   out_5005698098331862813[38] = 0;
   out_5005698098331862813[39] = 0;
   out_5005698098331862813[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_5005698098331862813[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_5005698098331862813[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_5005698098331862813[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_5005698098331862813[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_5005698098331862813[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_5005698098331862813[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_5005698098331862813[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_5005698098331862813[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_5005698098331862813[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_5005698098331862813[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5005698098331862813[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5005698098331862813[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_5005698098331862813[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_5005698098331862813[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_5005698098331862813[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5005698098331862813[56] = 0;
   out_5005698098331862813[57] = 0;
   out_5005698098331862813[58] = 0;
   out_5005698098331862813[59] = 0;
   out_5005698098331862813[60] = 0;
   out_5005698098331862813[61] = 0;
   out_5005698098331862813[62] = 0;
   out_5005698098331862813[63] = 1;
}
void h_25(double *state, double *unused, double *out_1051757273230348996) {
   out_1051757273230348996[0] = state[6];
}
void H_25(double *state, double *unused, double *out_6043417528456999617) {
   out_6043417528456999617[0] = 0;
   out_6043417528456999617[1] = 0;
   out_6043417528456999617[2] = 0;
   out_6043417528456999617[3] = 0;
   out_6043417528456999617[4] = 0;
   out_6043417528456999617[5] = 0;
   out_6043417528456999617[6] = 1;
   out_6043417528456999617[7] = 0;
}
void h_24(double *state, double *unused, double *out_5022686355337100841) {
   out_5022686355337100841[0] = state[4];
   out_5022686355337100841[1] = state[5];
}
void H_24(double *state, double *unused, double *out_5973588777829285175) {
   out_5973588777829285175[0] = 0;
   out_5973588777829285175[1] = 0;
   out_5973588777829285175[2] = 0;
   out_5973588777829285175[3] = 0;
   out_5973588777829285175[4] = 1;
   out_5973588777829285175[5] = 0;
   out_5973588777829285175[6] = 0;
   out_5973588777829285175[7] = 0;
   out_5973588777829285175[8] = 0;
   out_5973588777829285175[9] = 0;
   out_5973588777829285175[10] = 0;
   out_5973588777829285175[11] = 0;
   out_5973588777829285175[12] = 0;
   out_5973588777829285175[13] = 1;
   out_5973588777829285175[14] = 0;
   out_5973588777829285175[15] = 0;
}
void h_30(double *state, double *unused, double *out_5719077953120001940) {
   out_5719077953120001940[0] = state[4];
}
void H_30(double *state, double *unused, double *out_3570481382492183663) {
   out_3570481382492183663[0] = 0;
   out_3570481382492183663[1] = 0;
   out_3570481382492183663[2] = 0;
   out_3570481382492183663[3] = 0;
   out_3570481382492183663[4] = 1;
   out_3570481382492183663[5] = 0;
   out_3570481382492183663[6] = 0;
   out_3570481382492183663[7] = 0;
}
void h_26(double *state, double *unused, double *out_7725992852928704031) {
   out_7725992852928704031[0] = state[7];
}
void H_26(double *state, double *unused, double *out_6698733390313342726) {
   out_6698733390313342726[0] = 0;
   out_6698733390313342726[1] = 0;
   out_6698733390313342726[2] = 0;
   out_6698733390313342726[3] = 0;
   out_6698733390313342726[4] = 0;
   out_6698733390313342726[5] = 0;
   out_6698733390313342726[6] = 0;
   out_6698733390313342726[7] = 1;
}
void h_27(double *state, double *unused, double *out_108247420328686561) {
   out_108247420328686561[0] = state[3];
}
void H_27(double *state, double *unused, double *out_4858063370328808975) {
   out_4858063370328808975[0] = 0;
   out_4858063370328808975[1] = 0;
   out_4858063370328808975[2] = 0;
   out_4858063370328808975[3] = 1;
   out_4858063370328808975[4] = 0;
   out_4858063370328808975[5] = 0;
   out_4858063370328808975[6] = 0;
   out_4858063370328808975[7] = 0;
}
void h_29(double *state, double *unused, double *out_7078705484129371768) {
   out_7078705484129371768[0] = state[1];
}
void H_29(double *state, double *unused, double *out_2995922768636647488) {
   out_2995922768636647488[0] = 0;
   out_2995922768636647488[1] = 1;
   out_2995922768636647488[2] = 0;
   out_2995922768636647488[3] = 0;
   out_2995922768636647488[4] = 0;
   out_2995922768636647488[5] = 0;
   out_2995922768636647488[6] = 0;
   out_2995922768636647488[7] = 0;
}
void h_28(double *state, double *unused, double *out_3445327860317372799) {
   out_3445327860317372799[0] = state[5];
   out_3445327860317372799[1] = state[6];
}
void H_28(double *state, double *unused, double *out_5455322817643771575) {
   out_5455322817643771575[0] = 0;
   out_5455322817643771575[1] = 0;
   out_5455322817643771575[2] = 0;
   out_5455322817643771575[3] = 0;
   out_5455322817643771575[4] = 0;
   out_5455322817643771575[5] = 1;
   out_5455322817643771575[6] = 0;
   out_5455322817643771575[7] = 0;
   out_5455322817643771575[8] = 0;
   out_5455322817643771575[9] = 0;
   out_5455322817643771575[10] = 0;
   out_5455322817643771575[11] = 0;
   out_5455322817643771575[12] = 0;
   out_5455322817643771575[13] = 0;
   out_5455322817643771575[14] = 1;
   out_5455322817643771575[15] = 0;
}
}

extern "C"{
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
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
