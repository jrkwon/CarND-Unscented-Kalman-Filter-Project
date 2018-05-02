#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2; //30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3; //30;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  ///* State dimension
  n_x_ = 5;

  ///* Augmented state dimension
  n_aug_ = 7;

  ///* Sigma point spreading parameter
  lambda_ = 3 - n_x_;

  ///* Weights of sigma points
  weights_ = VectorXd(2*n_aug_ + 1);
  weights_(0) = lambda_/(lambda_ + n_aug_);
  double weight = 0.5/(lambda_ + n_aug_);

  for (int i = 1; i < 2*n_aug_ + 1; i++) {
    weights_(i) = weight;
  }
 
  ///* initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  ////////////////////
  ///* Initialization

  if (is_initialized_ == false) {
    ///* set initial values

    // initialize states
    x_.fill(1.0);
    
    // initialize prediction covariance matrix
    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

    double px = 0.0, py = 0.0;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = meas_package.raw_measurements_[0];
      double psi = meas_package.raw_measurements_[1];
      //double rho_dot = meas_package.raw_measurements_[2];

      px = rho*cos(psi);
      py = rho*sin(psi);
      //double vx = rho_dot*cos(psi);
      //double vy = rho_dot*sin(psi);

      //double psi_dot = atan2(vy, vx);

    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      px = meas_package.raw_measurements_[0];
      py = meas_package.raw_measurements_[1];

    }

    // px and py shouldn't be 0.
    if (px == 0 && py == 0)
    {
      px = py = 0.0001;
    }

    x_ << px, py, 0, 0, 0;

    time_us_ = meas_package.timestamp_;

    // done initializing
    is_initialized_ = true;
    return;
  }

  ////////////////
  ///* Prediction
  
  //compute the time elapsed between the current and previous measurements
  double dt_in_sec = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
  time_us_ = meas_package.timestamp_;

  // make prediction of sigma points
  Prediction(dt_in_sec);

  ////////////
  ///* Update

  if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  }
  else if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  }
  
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double dt_in_sec) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  ///////////////////////////
  ///* Generate sigma points

  // Augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;

  // Augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_*std_a_;
  P_aug(n_aug_ - 1, n_aug_ - 1) = std_yawdd_*std_yawdd_;
  
  // Square root matrix
  MatrixXd A_aug = P_aug.llt().matrixL();
  // Sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_ + 1);
  // set first column of sigma point matrix
  Xsig_aug.col(0) = x_aug;

  // set remaining sigma points
  for (int i = 0; i < n_aug_; i++){
    Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_)*A_aug.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_)*A_aug.col(i);
  }
  
  ////////////////////////////
  ///* Sigma point predcition
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_ + 1);

  // Predict sigma points
  for (int i = 0; i < 2*n_aug_ + 1; i++){
    // Just for readability
    double px = Xsig_aug(0, i);
    double py = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    // Predicted states
    double px_p, py_p;

    // To avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = px + v/yawd*(sin(yaw + yawd*dt_in_sec) - sin(yaw));
      py_p = py + v/yawd*(cos(yaw) - cos(yaw + yawd*dt_in_sec));
    }
    else {
      px_p = px + v*dt_in_sec*cos(yaw);
      py_p = py + v*dt_in_sec*sin(yaw);
    }

    // Predicted v, yaw, yawd
    double v_p = v;
    double yaw_p = yaw + yawd*dt_in_sec;
    double yawd_p = yawd;

    // Add noise
    px_p += 0.5*nu_a*dt_in_sec*dt_in_sec*cos(yaw);
    py_p += 0.5*nu_a*dt_in_sec*dt_in_sec*sin(yaw);
    v_p  += nu_a*dt_in_sec;
    yaw_p += 0.5*nu_yawdd*dt_in_sec*dt_in_sec;
    yawd_p += nu_yawdd*dt_in_sec;

    // Predited sigma points goes into right column
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }

  // Predict state mean
  x_.fill(0.0);
  for (int i = 0; i < 2*n_aug_ + 1; i++) {  
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2*n_aug_ + 1; i++) {  
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    NormalizeAngle(x_diff(3));
    P_ += weights_(i)*x_diff*x_diff.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  // Sensor state dimension
  int n_z = 2;

  ///*

  // Sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_ + 1);

  // Transform sigma points into measurement space 
  for (int i = 0; i < 2*n_aug_ + 1; i++) {  
    // Measurement model
    Zsig(0, i) = Xsig_pred_(0, i);
    Zsig(1, i) = Xsig_pred_(1, i);
  }

  // Mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < 2*n_aug_ + 1; i++) {
    z_pred += weights_(i)*Zsig.col(i);
  }

  // Measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);

  for (int i = 0; i < 2*n_aug_ + 1; i++) { 
    VectorXd z_diff = Zsig.col(i) - z_pred;
    NormalizeAngle(z_diff(1));
    S += weights_(i)*z_diff*z_diff.transpose();
  }

  // Measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_laspx_*std_laspx_, 0, 0, std_laspy_*std_laspy_;
  S += R;

  ///* 

  // Matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  //calculate cross correlation matrix
  for (int i = 0; i < 2*n_aug_ + 1; i++) { 
    VectorXd z_diff = Zsig.col(i) - z_pred;
    NormalizeAngle(z_diff(1));
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    NormalizeAngle(x_diff(3));
    
    //Cross-correlation matrix
    Tc += weights_(i)*x_diff*z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc*S.inverse();

  // actual measurement
  VectorXd z = VectorXd(n_z);
  z = meas_package.raw_measurements_;
  
  // residual
  VectorXd err = z - z_pred;
  NormalizeAngle(err(1));

  //update state mean and covariance matrix
  x_ += K*err;
  P_ -= K*S*K.transpose();

  //Calculate the lidar NIS.
  laser_NIS_ = err.transpose()*S.inverse()*err;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

   // sensor state dimension
  int n_z = 3;

  ///*

  // Matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_ + 1);

  for (int i = 0; i < 2*n_aug_ + 1; i++) {
    // Just for readability
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    // measurement model
    Zsig(0, i) = sqrt(px*px + py*py);    //rho 
    Zsig(1, i) = atan2(py, px);          //phi
    
    if (fabs(Zsig(0, i)) > 0.0001)
      Zsig(2, i) = (px*cos(yaw)*v + py*sin(yaw)*v) / Zsig(0, i);   //rho_dot
    else
      Zsig(2, i) = 0;
  }

  // Mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < 2*n_aug_ + 1; i++) {
    z_pred += weights_(i)*Zsig.col(i);
  }

  // Measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < 2*n_aug_ + 1; i++) {  //2n+1 simga points
    VectorXd z_diff = Zsig.col(i) - z_pred;
    NormalizeAngle(z_diff(1));
    S += weights_(i)*z_diff*z_diff.transpose();
  }

  // Add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_radr_*std_radr_, 0, 0, 0, std_radphi_*std_radphi_, 0,
        0, 0, std_radrd_*std_radrd_;
  S += R;

  ///* 

  // Create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  // Calculate cross correlation matrix
  for (int i = 0; i < 2*n_aug_ + 1; i++) {  //2n+1 simga points
    VectorXd z_diff = Zsig.col(i) - z_pred;
    NormalizeAngle(z_diff(1));        

    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    NormalizeAngle(x_diff(3));    
    
    Tc += weights_(i)*x_diff*z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc*S.inverse();

  // Actual measurement
  VectorXd z = VectorXd(n_z);
  z = meas_package.raw_measurements_;

  VectorXd err = z - z_pred;

  NormalizeAngle(err(1));

  //update state mean and covariance matrix
  x_ += K*err;
  P_ -= K*S*K.transpose();

  radar_NIS_ = err.transpose()*S.inverse()*err;
}

double UKF::NormalizeAngle(double& angle) 
{
  while (angle >  M_PI) angle -= 2.*M_PI;
  while (angle < -M_PI) angle += 2.*M_PI;

  return angle;
}