#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  //initializing matrices
  R_laser_ = MatrixXd(2, 2); //covariance matrix for laser
  R_radar_ = MatrixXd(3, 3); //covariance matrix for radar

  H_laser_ = MatrixXd(2, 4); //measurement matrix
  Hj_ = MatrixXd(3, 4); //jacobian matrix for radar

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  //measurement matrix for laser
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  //state transition 4x4 matrix, lesson 8
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;

  //state covariance 4x4 matrix for process noise
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1000, 0,
             0, 0, 0, 1000;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack){

  /*************************************************************
   *  Initialization
  *************************************************************/
  if (!is_initialized_) {

    //first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1; //important to the RMSE

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
 
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      float rho_dot = measurement_pack.raw_measurements_[2];
			
      //convert radar from polar to cartesian coordinates and initialize state
      ekf_.x_ << rho*cos(phi), rho*sin(phi), 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
  
      //initialize state ekf_.x_
      ekf_.x_<< measurement_pack.raw_measurements_[0], 
                measurement_pack.raw_measurements_[1], 0, 0;
    }

    previous_timestamp_ = measurement_pack.timestamp_;

    //done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**********************************************************
   *  Prediction
   **********************************************************/

  //set the acceleration noise components
  float noise_ax = 9; //sect 13 lesson 5
  float noise_ay = 9;

  //convert dt to seconds
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; 

  previous_timestamp_ = measurement_pack.timestamp_;

  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  //update state transition F matrix w/ elaspsed time
  //sect 8 lesson 5
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  //update process noise covariance matrix Q 
  //sect 9 lesson 5
  ekf_.Q_ = MatrixXd(4,4);
  ekf_.Q_ << dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
             0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
             dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
             0, dt_3/2*noise_ay, 0, dt_2*noise_ay;


  //call kalman filter Predict function
  ekf_.Predict();

  /**********************************************************
   *  Update
   **********************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    //radar updates
    ekf_.R_ = R_radar_;

    //call jacobian function
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);

    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    //laser updates
    ekf_.R_ = R_laser_;
    ekf_.H_ = H_laser_;

    ekf_.Update(measurement_pack.raw_measurements_);  
  }

  //print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
