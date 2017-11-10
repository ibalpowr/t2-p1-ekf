#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

//predict the state
void KalmanFilter::Predict() {

  //noise w/ zero mean
  x_ = F_ * x_; //sect8 lesson5
  MatrixXd Ft = F_.transpose(); 
  P_ = F_ * P_ * Ft + Q_; //sect9 lesson5
}

//update laser
void KalmanFilter::Update(const VectorXd &z) {

	//sect7 lesson5
  VectorXd y = z - H_ * x_;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * Ht * Si;

	//new state lesson12
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);

  P_ = (I - K * H_) * P_;
}

//update radar
void KalmanFilter::UpdateEKF(const VectorXd &z) {

  //sect 14 lesson 5

  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);

  float rho = sqrt(px*px+py*py);
  float phi = atan2(py,px);
  float rho_dot = (px*vx+py*vy)/rho;

  VectorXd z_pred = VectorXd(3);
  z_pred << rho, phi, rho_dot;
  VectorXd y = z - z_pred;

  //angle normalization
  y[1] = atan2(sin(y[1]),cos(y[1]));
	
	//sect7 lesson 5
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * Ht * Si;

  //new state
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
