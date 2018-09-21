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
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0, 
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;

  Q_ = MatrixXd(5, 5);
  Q_ << 0.1, 0, 0, 0, 0,
        0, 0.1, 0, 0, 0,
        0, 0, 0.1, 0, 0, 
        0, 0, 0, 0.1, 0,
        0, 0, 0, 0, 0.1;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
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
  
  R_lidar_ = MatrixXd(2, 2);
  R_lidar_ << std_laspx_ * std_laspx_, 0,
              0, std_laspy_ * std_laspy_;

  R_radar_ = MatrixXd(3,3);
  R_radar_ << std_radr_ * std_radr_, 0, 0,
              0, std_radphi_ * std_radphi_, 0,
              0, 0, std_radrd_ * std_radrd_;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  n_x_ = 5;
  n_aug_ = 2 * n_x_ + 1;
  is_initialized_ = false;
  Xsig_pred_ = MatrixXd(n_aug_, n_x_);
  weights_ = VectorXd(n_aug_);
  float alpha = 1.2;
  float kappa = 0.0;
  float beta = 2;
  lambda_ = alpha * alpha * (n_x_ + kappa) - n_x_;

  weights_(0) = lambda_ / (n_x_ + lambda_);

  for (int i = 1; i < n_aug_; ++i) {
    weights_(i) = 0.5/(n_x_ + lambda_);
  }
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

  if (!is_initialized_)
  {
    //initialize the state

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      float ro = meas_package.raw_measurements_(0);
      float theta = meas_package.raw_measurements_(1);
      float ro_dot = meas_package.raw_measurements_(2);

      float x = ro * cos(theta);
      float y = ro * sin(theta);
      float vx = ro_dot * cos(theta);
      float vy = ro_dot * sin(theta);

      float v = sqrt(vx*vx + vy*vy);
      x_ << x, y, v, 0, 0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      float x = meas_package.raw_measurements_(0);
      float y = meas_package.raw_measurements_(1);

      x_ << x, y, 0, 0, 0;
    }
  
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;

    cout << "UKF \n";
    cout << "State X : \n";
    cout << x_ << endl;
    cout << "covariance P : \n";
    cout << P_ << endl;

    return;
  }

  double dt = (meas_package.timestamp_ - time_us_)/1000000.0;
  time_us_ = meas_package.timestamp_;

  Prediction(dt);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    UpdateRadar(meas_package);
  } 
  else 
  {
    UpdateLidar(meas_package);
  }

  cout << "UKF \n";
  cout << "State X : \n";
  cout << x_ << endl;
  cout << "covariance P : \n";
  cout << P_ << endl;

}

VectorXd UKF::ProcessModel(VectorXd X, double dt)
{
  VectorXd out_X = VectorXd(n_x_);
  
  out_X(0) = X(0) + X(2) * cos(X(3)) * dt;
  out_X(1) = X(1) + X(2) * sin(X(3)) * dt;
  out_X(2) = X(2);
  out_X(3) = X(3) + X(4) * dt;
  out_X(4) = X(4);

  return out_X;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  ComputeSigmaPoints();
  
  // pass the sigma points through the motion model
  for (int i = 0; i < n_aug_; ++i) {
    Xsig_pred_.row(i) = ProcessModel(Xsig_pred_.row(i), delta_t);
  }

  x_.fill(0.0);
  // compute the mean
  for (int i = 0; i < n_aug_; ++i) {
    VectorXd X_sigma = Xsig_pred_.row(i);
    x_ = x_ + weights_(i) * X_sigma;
  }

  P_.fill(0.0);
  // compute the covariance
  for (int i = 0; i < n_aug_; ++i) {
    VectorXd X_sigma = Xsig_pred_.row(i);
    VectorXd x_diff = X_sigma - x_;
    while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

    P_ = P_ + weights_(i) * (x_diff) * (x_diff).transpose();
  }
  P_ += Q_;
}

void UKF::ComputeSigmaPoints() {
  // calculate the square root for the addition term
  MatrixXd matrix_sqrt = (n_x_ + lambda_) * P_;
  matrix_sqrt = matrix_sqrt.llt().matrixL();
  MatrixXd sigma_add = matrix_sqrt;

  // Set the first sigma point
  Xsig_pred_.row(0) = x_;

  // Set the other sigma points
  for (int i = 1; i <= n_x_; ++i){
    VectorXd add = sigma_add.col(i-1);
    Xsig_pred_.row(i) = x_ + add;
    Xsig_pred_.row(n_x_+i) = x_ - add;
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

  // Compute sigma points for measurement
  ComputeSigmaPoints();
  MatrixXd Sigma_z = MatrixXd(n_aug_,2);
  MatrixXd C_lidar = MatrixXd(2, 5);
  C_lidar << 1, 0, 0, 0, 0,
             0, 1, 0, 0, 0;

  // Pass the sigma points through measurement model
  for (int i = 0; i < n_aug_; ++i) {
    VectorXd X_sigma = Xsig_pred_.row(i);
    Sigma_z.row(i) = C_lidar * X_sigma;
  }

  // Compute the mean of the measurement sigma points
  VectorXd Z = VectorXd(2);
  Z << 0, 0;
  for (int i = 0; i < n_aug_; i++) {
    VectorXd Z_sigma = Sigma_z.row(i);
    Z = Z + weights_(i) * Z_sigma;
  }

  // Compute covariance of the measurement sigma points
  MatrixXd S_lidar = MatrixXd(2,2);
  S_lidar << 0, 0,
             0, 0;
  for (int i = 0; i < n_aug_; ++i) {
    VectorXd Z_sigma = Sigma_z.row(i);
    S_lidar = S_lidar + weights_(i) * (Z_sigma - Z) * (Z_sigma - Z).transpose();
  }
  S_lidar += R_lidar_;

  MatrixXd S = MatrixXd(5,2);
  S.fill(0.0);
  
  for (int i = 0; i < n_aug_; ++i) {
    VectorXd X_sigma = Xsig_pred_.row(i);
    VectorXd Z_sigma = Sigma_z.row(i);
    
    VectorXd x_diff = X_sigma - x_;
    while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

    S = S + weights_(i) * (x_diff) * (Z_sigma - Z).transpose();
  }
  MatrixXd Kalman_gain = S * S_lidar.inverse();
  x_ = x_ + Kalman_gain * (meas_package.raw_measurements_ - Z);
  P_ = P_ - Kalman_gain * S_lidar * Kalman_gain.transpose();
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

  ComputeSigmaPoints();
  MatrixXd Sigma_z = MatrixXd(n_aug_,3);

  // Pass the sigma points through measurement model
  for (int i = 0; i < n_aug_; ++i) {
    VectorXd X_sigma = Xsig_pred_.row(i);
    Sigma_z(i,0) = sqrt(X_sigma(0)*X_sigma(0) + X_sigma(1)*X_sigma(1));
    Sigma_z(i,1) = atan2(X_sigma(1), X_sigma(0));
    float vx = X_sigma(2) * cos(X_sigma(3));
    float vy = X_sigma(2) * sin(X_sigma(3));
    Sigma_z(i,2) = (X_sigma(0) * vx + X_sigma(1) * vy) / Sigma_z(i,0);
  }

  // Compute the mean of the measurement sigma points
  VectorXd Z = VectorXd(3);
  Z << 0, 0, 0;
  for (int i = 0; i < n_aug_; i++) {
    VectorXd Z_sigma = Sigma_z.row(i);
    Z = Z + weights_(i) * Z_sigma;
  }

  // Compute covariance of the measurement sigma points
  MatrixXd S_radar = MatrixXd(3,3);
  S_radar << 0, 0, 0,
             0, 0, 0,
             0, 0, 0;

  for (int i = 0; i < n_aug_; ++i) {
    VectorXd Z_sigma = Sigma_z.row(i);
    VectorXd z_diff = Z_sigma - Z;
    while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;
    S_radar = S_radar + weights_(i) * (z_diff) * (z_diff).transpose();
  }
  S_radar += R_radar_;
  MatrixXd S = MatrixXd(5,3);
  S << 0, 0, 0,   
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0;
  
  for (int i = 0; i < n_aug_; ++i) {
    VectorXd X_sigma = Xsig_pred_.row(i);
    VectorXd Z_sigma = Sigma_z.row(i);

    VectorXd z_diff = Z_sigma - Z;
    while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

    VectorXd x_diff = X_sigma - x_;
    while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;
    S = S + weights_(i) * (x_diff) * (z_diff).transpose();
  }

  MatrixXd Kalman_gain = S * S_radar.inverse();
  VectorXd z_diff = meas_package.raw_measurements_ - Z; 
  while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
  while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;   
  x_ = x_ + Kalman_gain * (z_diff);
  P_ = P_ - Kalman_gain * S_radar * Kalman_gain.transpose();
}
