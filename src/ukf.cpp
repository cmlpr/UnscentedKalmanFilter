#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
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
  // std_a_ = 30;
  // Change to a reasonable value
  std_a_ = 0.8;

  // Process noise standard deviation yaw acceleration in rad/s^2
  // std_yawdd_ = 30;
  // Change to reasonable value
  std_yawdd_ = 0.6;

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  std_a_ and std_yawdd_
  */

  is_initialized_ = false;

  ///* State dimension
  n_x_ = 5;
  
  ///* Augmented state dimension
  n_aug_ = n_x_ + 2;

  ///* Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  ///* predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  ///* time when the state is true, in us
  time_us_ = 0.0;

  ///* Weights of sigma points
  weights_ = VectorXd(2 * n_aug_ + 1);

  // Laser
  H_laser_ = MatrixXd(2, 5);
  R_laser_ = MatrixXd(2, 2);

  // NIS
  NIS_radar_ = 0.0;
  NIS_laser_ = 0.0;

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

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_ ) {

    // Initialize state matrix
    x_ << 0.1, 0.1, 0.1, 0.1, 0.01;

    //initial state covariance matrix
    P_ << 0.2, 0,   0,   0,   0,
          0,   0.2, 0,   0,   0,
          0,   0,   0.2, 0,   0,
          0,   0,   0,   0.2, 0,
          0,   0,   0,   0,   0.3;

    time_us_ = meas_package.timestamp_;

    if(meas_package.raw_measurements_[0]==0 || meas_package.raw_measurements_[1]==0){
      cout << "Px or Py is equal to zero - skipping measurement point" << endl;
      return;
    }

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho       = meas_package.raw_measurements_[0];
      float theta     = meas_package.raw_measurements_[1];
      float rho_dot   = meas_package.raw_measurements_[2];
      
      // px
      x_(0) = rho * cos(theta);
      // py
      x_(1) = rho * sin(theta);
      // v
      // x_(2) = rho_dot * rho / (x_(0) * cos(theta) + x_(1) * sin(theta));
      x_(2) = sqrt(pow(rho_dot * cos(theta), 2.0) +  
                   pow(rho_dot * sin(theta), 2.0));
      // // yaw
      // x_(3) = 0.1;
      // if(rho_dot * cos(theta) != 0){
      //   x_(3) = (rho_dot * sin(theta)) / (rho_dot * cos(theta));
      // }
      // // yaw rate
      // x_(4) = 0.01;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
      /**
      Initialize state.
      */
      x_(0) = meas_package.raw_measurements_[0];
      x_(1) = meas_package.raw_measurements_[1];
      // x_(2) = 0.1;
      // x_(3) = 0.1;
      // x_(4) = 0.1;
    }

    is_initialized_ = true;

    return;
  }

  /*****************************************************************************
  *  Prediction
  ****************************************************************************/

  //compute the time elapsed between the current and previous measurements
  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
  
  Prediction(dt);

  if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    UpdateLidar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    UpdateRadar(meas_package);
  }

  time_us_ = meas_package.timestamp_;

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

  /*****************************************************************************
   *  Create Augmented Sigma Points
   ****************************************************************************/

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  
  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  
  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  
  //create augmented mean state - added two noise values to the mean of the augmented state. 
  // The mean of the noise values is zero, so I set these numbers to zero.
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5, 5) = P_;
  P_aug.bottomRightCorner(2,2) << std_a_ * std_a_, 0.0, 0.0, std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd A_aug = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;

  for (int i = 0; i < n_aug_; i++)
  {
    Xsig_aug.col(i+1)     = x_aug + sqrt(lambda_ + n_aug_) * A_aug.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A_aug.col(i);
  }

  // std::cout << "x_aug = " << std::endl << x_aug << std::endl;
  // std::cout << "A_aug = " << std::endl << A_aug << std::endl;
  // std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;

  /*****************************************************************************
  *  Predict Sigma Points
  ****************************************************************************/

  for (int i = 0; i < 2 * n_aug_ + 1; i++) 
  {
    double px = Xsig_aug(0, i);
    double py = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yaw_rate = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yaw_rate = Xsig_aug(6, i);

    double a_noise = 0.5 * pow(delta_t, 2.0) * cos(yaw) * nu_a;
    double b_noise = 0.5 * pow(delta_t, 2.0) * sin(yaw) * nu_a;
    double c_noise = delta_t * nu_a;
    double d_noise = 0.5 * pow(delta_t, 2.0) * nu_yaw_rate;
    double e_noise = delta_t * nu_yaw_rate;
    
    // Avoid division by zero
    if (fabs(yaw_rate) > 0.001) {
      Xsig_pred_(0,i) = px + (v/yaw_rate) * (sin(yaw + yaw_rate * delta_t) - sin(yaw)) + a_noise;
      Xsig_pred_(1,i) = py + (v/yaw_rate) * (-cos(yaw + yaw_rate * delta_t) + cos(yaw)) + b_noise;
      Xsig_pred_(2,i) = v + c_noise;
      Xsig_pred_(3,i) = yaw + yaw_rate * delta_t + d_noise;
      Xsig_pred_(4,i) = yaw_rate + e_noise;
    } else {
      Xsig_pred_(0,i) = px + v * cos(yaw) * delta_t + a_noise;
      Xsig_pred_(1,i) = py + v * sin(yaw) * delta_t + b_noise;
      Xsig_pred_(2,i) = v + c_noise;
      Xsig_pred_(3,i) = yaw + d_noise;
      Xsig_pred_(4,i) = yaw_rate + e_noise;
    } 
  }

  // std::cout << "Xsig_pred_ = " << std::endl << Xsig_pred_ << std::endl;

  /*****************************************************************************
  *  Predict Mean and Covariance
  ****************************************************************************/

  //set weights
  for (int i = 0; i < 2 * n_aug_ + 1; i++) 
  {
      if (i == 0) {
        weights_(i) = lambda_ / (lambda_ + n_aug_);
      }
      else {
        weights_(i) = 1.0 / 2.0 / (lambda_ + n_aug_);
      }
  }
  //predict state mean
  x_.fill(0.0);
  for (int j = 0; j < 2 * n_aug_ + 1; j++) {
      x_ += Xsig_pred_.col(j) * weights_(j);
  }
  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

      // state difference
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      //angle normalization
      // In the equation we always need the difference between the mean predicted state and a sigma points. 
      // The problem here is that the state contains an angle. As you have learned before, subtracting angles is 
      // a problem for Kalman filters, because the result might be 2π plus a small angle, instead of just a small 
      // angle. That’s why I normalize the angle here.
      while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
      while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

      P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }

  // std::cout << "x_ = " << std::endl << x_ << std::endl;
  // std::cout << "P_ = " << std::endl << P_ << std::endl;
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

  /*****************************************************************************
  *  Update Lidar Measurements 
  ****************************************************************************/

  // The mapping from state space to Lidar is linear. Go back and look at the equations
  H_laser_ << 1, 0, 0, 0, 0,
              0, 1, 0, 0, 0;

  R_laser_ << std_laspx_ * std_laspx_, 0,
              0, std_laspx_ * std_laspy_;

  VectorXd z_pred = H_laser_ * x_;
	VectorXd y = meas_package.raw_measurements_ - z_pred;
	MatrixXd Ht = H_laser_.transpose();
	MatrixXd S = H_laser_ * P_ * Ht + R_laser_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_laser_) * P_;

  // NIS - normalized innovation square
  NIS_laser_ = y.transpose() * S.inverse() * y;

  std::cout << "NIS_laser_ = " << NIS_laser_ << std::endl;
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

  /*****************************************************************************
  *  Predict Radar Sigma Points
  ****************************************************************************/

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
          
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  //transform sigma points into measurement space
  Zsig.fill(0.0);
  for (int j = 0; j < 2 * n_aug_ + 1; j++) {
      Zsig(0,j) = sqrt(Xsig_pred_(0,j) * Xsig_pred_(0,j) + Xsig_pred_(1,j) * Xsig_pred_(1,j));
      Zsig(1,j) = atan2(Xsig_pred_(1,j), Xsig_pred_(0,j));
      Zsig(2,j) = (Xsig_pred_(0,j) * cos(Xsig_pred_(3,j)) * Xsig_pred_(2,j) + Xsig_pred_(1,j) * sin(Xsig_pred_(3,j)) * Xsig_pred_(2,j));
      Zsig(2,j) = Zsig(2,j) / Zsig(0,j);
  }
  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int j = 0; j < 2 * n_aug_ + 1; j++) {
      z_pred += Zsig.col(j) * weights_(j);
  }
  //calculate measurement covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

      // state difference
      VectorXd x_diff = Zsig.col(i) - z_pred;
      //angle normalization
      while (x_diff(1)> M_PI) x_diff(1)-=2.*M_PI;
      while (x_diff(1)<-M_PI) x_diff(1)+=2.*M_PI;

      S = S + weights_(i) * x_diff * x_diff.transpose() ;
  }

  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_ * std_radr_, 0, 0,
          0, std_radphi_ * std_radphi_, 0,
          0, 0, std_radrd_ * std_radrd_;
  S = S + R;

  /*****************************************************************************
  *  Update Radar
  ****************************************************************************/
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

      // state difference
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      VectorXd z_diff = Zsig.col(i) - z_pred;
      //angle normalization
      while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
      while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
      while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
      while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

      Tc += weights_(i) * x_diff * z_diff.transpose() ;
  }
  //calculate Kalman gain K;
  MatrixXd Kg = Tc * S.inverse();

  //update state mean and covariance matrix
  //residual
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;
  
  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  
  //update state mean and covariance matrix
  x_ = x_ + Kg * z_diff;
  P_ = P_ - Kg * S * Kg.transpose();

  // NIS - normalized innovation square
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;

  std::cout << "NIS_radar_ = " << NIS_radar_ << std::endl;
}