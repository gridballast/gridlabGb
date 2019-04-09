/*
 * gridballastcontroller.h
 *
 *  Created on: May 8, 2017
 *      Author: jingkungao
 */

#ifndef RESIDENTIAL_GRIDBALLASTCONTROLLER_H_
#define RESIDENTIAL_GRIDBALLASTCONTROLLER_H_

#include "residential.h"


namespace gridballastcontroller {

class gridballastcontroller {
private:
	// once initially set, we don't want to change them
	double freq_lowlimit;			// lower tripping limit of the frequency
	double freq_uplimit;			// upper tripping limit of the frequency

	double volt_lowlimit;
	double volt_uplimit;

	double T_setpoint;
	double T_deadband;

public:
	/* we don't seem to need initialize these variables as they are inputs to the controllers
	 * However, we do need to define them inside the class where the controller is implemented.

	// frequency control variables
	double freq_t; 					// grid frequency at time t

	bool circuit_status;			// True - ON; False - OFF, the returned variable to decide ON/OFF status
	bool enable_lock;				// force the circuit to be ON/OFF
	bool lock_STATUS;				// True - force ON; Flase - force OFF
	bool lock_OVERRIDE_TS;			// whether to override the thermostat setpoint for TCLs

	// thermostat control variables
	double T_t;						// temperature at time t

	// if True, we let the circuit ON when T_t < lower thermal band,
	// and let the circuit OFF when T_t > upper thermal band (for water heater)
	// if False, it is the other way around (for refrigerator, e.g.)
	bool reverse_ON_OFF;

	// we use this variable to toggle with/without frequency control
	bool enable_freq_control;
	*/
public:
	gridballastcontroller();

	gridballastcontroller(double freq_lowlimit,
			double freq_uplimit,double volt_lowlimit,
			double volt_uplimit,double T_setpoint,double T_deadband);

	virtual ~gridballastcontroller();

	void set_parameters(double freq_low,
			double freq_up,double volt_low,
			double volt_up,double T_sp,double T_db);

	bool check_freq_violation(double freq_t);
	bool check_volt_violation(double volt_t);
	bool check_thermal_violation(double T_t);

	bool lock_mode_controller(bool circuit_status, bool enable_lock, bool lock_STATUS);
	bool frequency_controller(bool circuit_status, double freq_t);
	bool voltage_controller(bool circuit_status, double volt_t);
	bool thermostat_controller(bool circuit_status,double T_t,bool reverse_ON_OFF);
};

} /* namespace gridballastcontroller */

#endif /* RESIDENTIAL_GRIDBALLASTCONTROLLER_H_ */
