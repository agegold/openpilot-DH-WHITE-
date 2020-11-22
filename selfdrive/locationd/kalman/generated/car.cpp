
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
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_6182305149975481084) {
   out_6182305149975481084[0] = delta_x[0] + nom_x[0];
   out_6182305149975481084[1] = delta_x[1] + nom_x[1];
   out_6182305149975481084[2] = delta_x[2] + nom_x[2];
   out_6182305149975481084[3] = delta_x[3] + nom_x[3];
   out_6182305149975481084[4] = delta_x[4] + nom_x[4];
   out_6182305149975481084[5] = delta_x[5] + nom_x[5];
   out_6182305149975481084[6] = delta_x[6] + nom_x[6];
   out_6182305149975481084[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_5772716229829584760) {
   out_5772716229829584760[0] = -nom_x[0] + true_x[0];
   out_5772716229829584760[1] = -nom_x[1] + true_x[1];
   out_5772716229829584760[2] = -nom_x[2] + true_x[2];
   out_5772716229829584760[3] = -nom_x[3] + true_x[3];
   out_5772716229829584760[4] = -nom_x[4] + true_x[4];
   out_5772716229829584760[5] = -nom_x[5] + true_x[5];
   out_5772716229829584760[6] = -nom_x[6] + true_x[6];
   out_5772716229829584760[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_6573258343469333187) {
   out_6573258343469333187[0] = 1.0;
   out_6573258343469333187[1] = 0.0;
   out_6573258343469333187[2] = 0.0;
   out_6573258343469333187[3] = 0.0;
   out_6573258343469333187[4] = 0.0;
   out_6573258343469333187[5] = 0.0;
   out_6573258343469333187[6] = 0.0;
   out_6573258343469333187[7] = 0.0;
   out_6573258343469333187[8] = 0.0;
   out_6573258343469333187[9] = 1.0;
   out_6573258343469333187[10] = 0.0;
   out_6573258343469333187[11] = 0.0;
   out_6573258343469333187[12] = 0.0;
   out_6573258343469333187[13] = 0.0;
   out_6573258343469333187[14] = 0.0;
   out_6573258343469333187[15] = 0.0;
   out_6573258343469333187[16] = 0.0;
   out_6573258343469333187[17] = 0.0;
   out_6573258343469333187[18] = 1.0;
   out_6573258343469333187[19] = 0.0;
   out_6573258343469333187[20] = 0.0;
   out_6573258343469333187[21] = 0.0;
   out_6573258343469333187[22] = 0.0;
   out_6573258343469333187[23] = 0.0;
   out_6573258343469333187[24] = 0.0;
   out_6573258343469333187[25] = 0.0;
   out_6573258343469333187[26] = 0.0;
   out_6573258343469333187[27] = 1.0;
   out_6573258343469333187[28] = 0.0;
   out_6573258343469333187[29] = 0.0;
   out_6573258343469333187[30] = 0.0;
   out_6573258343469333187[31] = 0.0;
   out_6573258343469333187[32] = 0.0;
   out_6573258343469333187[33] = 0.0;
   out_6573258343469333187[34] = 0.0;
   out_6573258343469333187[35] = 0.0;
   out_6573258343469333187[36] = 1.0;
   out_6573258343469333187[37] = 0.0;
   out_6573258343469333187[38] = 0.0;
   out_6573258343469333187[39] = 0.0;
   out_6573258343469333187[40] = 0.0;
   out_6573258343469333187[41] = 0.0;
   out_6573258343469333187[42] = 0.0;
   out_6573258343469333187[43] = 0.0;
   out_6573258343469333187[44] = 0.0;
   out_6573258343469333187[45] = 1.0;
   out_6573258343469333187[46] = 0.0;
   out_6573258343469333187[47] = 0.0;
   out_6573258343469333187[48] = 0.0;
   out_6573258343469333187[49] = 0.0;
   out_6573258343469333187[50] = 0.0;
   out_6573258343469333187[51] = 0.0;
   out_6573258343469333187[52] = 0.0;
   out_6573258343469333187[53] = 0.0;
   out_6573258343469333187[54] = 1.0;
   out_6573258343469333187[55] = 0.0;
   out_6573258343469333187[56] = 0.0;
   out_6573258343469333187[57] = 0.0;
   out_6573258343469333187[58] = 0.0;
   out_6573258343469333187[59] = 0.0;
   out_6573258343469333187[60] = 0.0;
   out_6573258343469333187[61] = 0.0;
   out_6573258343469333187[62] = 0.0;
   out_6573258343469333187[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_9174724217193728746) {
   out_9174724217193728746[0] = state[0];
   out_9174724217193728746[1] = state[1];
   out_9174724217193728746[2] = state[2];
   out_9174724217193728746[3] = state[3];
   out_9174724217193728746[4] = state[4];
   out_9174724217193728746[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_9174724217193728746[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_9174724217193728746[7] = state[7];
}
void F_fun(double *state, double dt, double *out_7673946190717308446) {
   out_7673946190717308446[0] = 1;
   out_7673946190717308446[1] = 0;
   out_7673946190717308446[2] = 0;
   out_7673946190717308446[3] = 0;
   out_7673946190717308446[4] = 0;
   out_7673946190717308446[5] = 0;
   out_7673946190717308446[6] = 0;
   out_7673946190717308446[7] = 0;
   out_7673946190717308446[8] = 0;
   out_7673946190717308446[9] = 1;
   out_7673946190717308446[10] = 0;
   out_7673946190717308446[11] = 0;
   out_7673946190717308446[12] = 0;
   out_7673946190717308446[13] = 0;
   out_7673946190717308446[14] = 0;
   out_7673946190717308446[15] = 0;
   out_7673946190717308446[16] = 0;
   out_7673946190717308446[17] = 0;
   out_7673946190717308446[18] = 1;
   out_7673946190717308446[19] = 0;
   out_7673946190717308446[20] = 0;
   out_7673946190717308446[21] = 0;
   out_7673946190717308446[22] = 0;
   out_7673946190717308446[23] = 0;
   out_7673946190717308446[24] = 0;
   out_7673946190717308446[25] = 0;
   out_7673946190717308446[26] = 0;
   out_7673946190717308446[27] = 1;
   out_7673946190717308446[28] = 0;
   out_7673946190717308446[29] = 0;
   out_7673946190717308446[30] = 0;
   out_7673946190717308446[31] = 0;
   out_7673946190717308446[32] = 0;
   out_7673946190717308446[33] = 0;
   out_7673946190717308446[34] = 0;
   out_7673946190717308446[35] = 0;
   out_7673946190717308446[36] = 1;
   out_7673946190717308446[37] = 0;
   out_7673946190717308446[38] = 0;
   out_7673946190717308446[39] = 0;
   out_7673946190717308446[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_7673946190717308446[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_7673946190717308446[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_7673946190717308446[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_7673946190717308446[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_7673946190717308446[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_7673946190717308446[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_7673946190717308446[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_7673946190717308446[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_7673946190717308446[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_7673946190717308446[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7673946190717308446[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7673946190717308446[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_7673946190717308446[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_7673946190717308446[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_7673946190717308446[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7673946190717308446[56] = 0;
   out_7673946190717308446[57] = 0;
   out_7673946190717308446[58] = 0;
   out_7673946190717308446[59] = 0;
   out_7673946190717308446[60] = 0;
   out_7673946190717308446[61] = 0;
   out_7673946190717308446[62] = 0;
   out_7673946190717308446[63] = 1;
}
void h_25(double *state, double *unused, double *out_628975419662729201) {
   out_628975419662729201[0] = state[6];
}
void H_25(double *state, double *unused, double *out_7448964055686931099) {
   out_7448964055686931099[0] = 0;
   out_7448964055686931099[1] = 0;
   out_7448964055686931099[2] = 0;
   out_7448964055686931099[3] = 0;
   out_7448964055686931099[4] = 0;
   out_7448964055686931099[5] = 0;
   out_7448964055686931099[6] = 1;
   out_7448964055686931099[7] = 0;
}
void h_24(double *state, double *unused, double *out_8549256257243174563) {
   out_8549256257243174563[0] = state[4];
   out_8549256257243174563[1] = state[5];
}
void H_24(double *state, double *unused, double *out_7925276905366183369) {
   out_7925276905366183369[0] = 0;
   out_7925276905366183369[1] = 0;
   out_7925276905366183369[2] = 0;
   out_7925276905366183369[3] = 0;
   out_7925276905366183369[4] = 1;
   out_7925276905366183369[5] = 0;
   out_7925276905366183369[6] = 0;
   out_7925276905366183369[7] = 0;
   out_7925276905366183369[8] = 0;
   out_7925276905366183369[9] = 0;
   out_7925276905366183369[10] = 0;
   out_7925276905366183369[11] = 0;
   out_7925276905366183369[12] = 0;
   out_7925276905366183369[13] = 1;
   out_7925276905366183369[14] = 0;
   out_7925276905366183369[15] = 0;
}
void h_26(double *state, double *unused, double *out_1612116933575715082) {
   out_1612116933575715082[0] = state[7];
}
void H_26(double *state, double *unused, double *out_2926951939376508717) {
   out_2926951939376508717[0] = 0;
   out_2926951939376508717[1] = 0;
   out_2926951939376508717[2] = 0;
   out_2926951939376508717[3] = 0;
   out_2926951939376508717[4] = 0;
   out_2926951939376508717[5] = 0;
   out_2926951939376508717[6] = 0;
   out_2926951939376508717[7] = 1;
}
void h_27(double *state, double *unused, double *out_356100148548021822) {
   out_356100148548021822[0] = state[3];
}
void H_27(double *state, double *unused, double *out_6872708239176032209) {
   out_6872708239176032209[0] = 0;
   out_6872708239176032209[1] = 0;
   out_6872708239176032209[2] = 0;
   out_6872708239176032209[3] = 1;
   out_6872708239176032209[4] = 0;
   out_6872708239176032209[5] = 0;
   out_6872708239176032209[6] = 0;
   out_6872708239176032209[7] = 0;
}
void h_29(double *state, double *unused, double *out_6898344091781638272) {
   out_6898344091781638272[0] = state[1];
}
void H_29(double *state, double *unused, double *out_5593307003910287637) {
   out_5593307003910287637[0] = 0;
   out_5593307003910287637[1] = 1;
   out_5593307003910287637[2] = 0;
   out_5593307003910287637[3] = 0;
   out_5593307003910287637[4] = 0;
   out_5593307003910287637[5] = 0;
   out_5593307003910287637[6] = 0;
   out_5593307003910287637[7] = 0;
}
void h_28(double *state, double *unused, double *out_6095167101378246893) {
   out_6095167101378246893[0] = state[5];
   out_6095167101378246893[1] = state[6];
}
void H_28(double *state, double *unused, double *out_3375372154212119283) {
   out_3375372154212119283[0] = 0;
   out_3375372154212119283[1] = 0;
   out_3375372154212119283[2] = 0;
   out_3375372154212119283[3] = 0;
   out_3375372154212119283[4] = 0;
   out_3375372154212119283[5] = 1;
   out_3375372154212119283[6] = 0;
   out_3375372154212119283[7] = 0;
   out_3375372154212119283[8] = 0;
   out_3375372154212119283[9] = 0;
   out_3375372154212119283[10] = 0;
   out_3375372154212119283[11] = 0;
   out_3375372154212119283[12] = 0;
   out_3375372154212119283[13] = 0;
   out_3375372154212119283[14] = 1;
   out_3375372154212119283[15] = 0;
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
