#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  MatrixXd P_ = MatrixXd(4,4);
  P_ << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1000, 0,
        0, 0, 0, 5000;
  MatrixXd F_ = MatrixXd(4,4);
  F_ << 1, 0, 1, 0,
	    0, 1, 0, 1,
    	0, 0, 1, 0,
	    0, 0, 0, 1;
  H_laser_ <<  1, 0, 0, 0,
               0, 1, 0, 0;
  VectorXd x_in = VectorXd(4);
  MatrixXd q_in = MatrixXd(4,4); 
  ekf_ = KalmanFilter();
  ekf_.Init(x_in, P_, F_, H_laser_, R_laser_, q_in);
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

static void convertPolarToCartesian(const VectorXd &pack, VectorXd &x) {
   float rho = pack(0);
   float bearing = pack(1);
   float rho_dot = pack(2);
   x << rho * cos(bearing), rho * sin(bearing), rho_dot * cos(bearing), rho_dot * sin(bearing);  	
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      convertPolarToCartesian(measurement_pack.raw_measurements_, ekf_.x_);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_ << measurement_pack.raw_measurements_(0),  measurement_pack.raw_measurements_(1), 
    		     2, 2; 
    
   }

    previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  // convert timestamp to seconds
  float noise_ax =  7; //1000.0;
  float noise_ay = 5; //200.0;
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  ekf_.F_ << 1, 0, dt, 0,
             0, 1, 0, dt,
             0, 0, 1, 0,
             0, 0, 0, 1;
  ekf_.Q_ << pow(dt, 4) / 4 * noise_ax, 0, pow(dt,3) / 2 * noise_ax, 0,
             0, pow(dt,4)/4 * noise_ay, 0, pow(dt,3)/2 * noise_ay,
             pow(dt, 3)/2 *noise_ax, 0, pow(dt, 2) * noise_ax, 0,
             0, pow(dt, 3)/2 * noise_ay, 0, pow(dt, 2)*noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.R_ = R_radar_;
    ekf_.H_ = Hj_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    Hj_ = ekf_.H_;
    R_radar_ = ekf_.R_;
  } else {
    // Laser updates
    ekf_.R_ = R_laser_;
    ekf_.H_ = H_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
    R_laser_ = ekf_.R_;
    H_laser_ = ekf_.H_;
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
  previous_timestamp_ = measurement_pack.timestamp_;
}
