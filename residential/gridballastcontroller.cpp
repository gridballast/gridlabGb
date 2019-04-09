/*
 * gridballastcontroller.cpp
 *
 *  Created on: May 8, 2017
 *      Author: jingkungao
 *  The GridBallast Controller to return ON/OFF status bases on frequency
 *  the controller will only affect the m_t, the idea is that,
 *  given the frequency information,
 *  decide whether to bring up the load or shut down the load.
 *  the controller can also be called directly to override the thermostat_controller
 *  by letting enable_lock to be True and providing values for lock_STATUS
 */

#include "gridballastcontroller.h"

#define TSTAT_PRECISION 0.01

namespace gridballastcontroller {

gridballastcontroller::gridballastcontroller() {
	// set default values for the private variables
	freq_lowlimit = 59.9;
	freq_uplimit = 60.1;
	volt_lowlimit = 119.1;
	volt_uplimit = 121.1;
	T_setpoint = 120;
	T_deadband = 2;
}

//overload constructor
gridballastcontroller::gridballastcontroller(double freq_low,
		double freq_up,double volt_low,
		double volt_up,double T_sp,double T_db) {
	// set private variables with chosen values during constructions
	freq_lowlimit = freq_low;
	freq_uplimit = freq_up;
	volt_lowlimit = volt_low;
	volt_uplimit = volt_up;
	T_setpoint = T_sp;
	T_deadband = T_db;
}

void gridballastcontroller::set_parameters(double freq_low,
		double freq_up,double volt_low,
		double volt_up,double T_sp,double T_db) {
	// methods which can be called to set private variables
	freq_lowlimit = freq_low;
	freq_uplimit = freq_up;
	volt_lowlimit = volt_low;
	volt_uplimit = volt_up;
	T_setpoint = T_sp;
	T_deadband = T_db;
}

gridballastcontroller::~gridballastcontroller() {
	// TODO Auto-generated destructor stub
}

bool gridballastcontroller::check_freq_violation(double freq_t) {
	return (freq_t < freq_lowlimit) || (freq_t > freq_uplimit);
}

bool gridballastcontroller::check_volt_violation(double volt_t) {
	return (volt_t < volt_lowlimit) || (volt_t > volt_uplimit);
}

bool gridballastcontroller::check_thermal_violation(double T_t) {
	return ((T_t - TSTAT_PRECISION) < (T_setpoint - T_deadband/2)) || ((T_t  + TSTAT_PRECISION) > (T_setpoint + T_deadband/2));
}

bool gridballastcontroller::lock_mode_controller(bool circuit_status, bool enable_lock, bool lock_STATUS){
	if (enable_lock) {
		circuit_status = lock_STATUS;
		return circuit_status;
	}
	else {
		// if lock is not enabled, using the previous circuit_stauts
		return circuit_status;
	}
}

bool gridballastcontroller::frequency_controller(bool circuit_status, double freq_t){
	// compare the frequency with lower and upper limit
	if (freq_t < freq_lowlimit){
		// turn OFF load if freq_t is too low
		circuit_status = false;
	}
	else if (freq_t > freq_uplimit){
		// turn ON the load if too high
		circuit_status = true;
	}
	// if none of the condition triggered, using the previous circuit_stauts
	return circuit_status;
}

bool gridballastcontroller::voltage_controller(bool circuit_status, double volt_t){
	// compare the voltage with lower and upper limit
	if (volt_t < volt_lowlimit){
		// turn OFF load if freq_t is too low
		circuit_status = false;
	}
	else if (volt_t > volt_uplimit){
		// turn ON the load if too high
		circuit_status = true;
	}
	// if none of the condition triggered, using the previous circuit_stauts
	return circuit_status;
}

bool gridballastcontroller::thermostat_controller(bool circuit_status,double T_t,bool reverse_ON_OFF){

	// compare the thermostat setpoint with the deadband
	if ((T_t - TSTAT_PRECISION) < (T_setpoint - T_deadband/2)) {
		// turn on the load if reverse_ON_OFF false otherwise turn off
		circuit_status = !reverse_ON_OFF;
	}
	else if ((T_t  + TSTAT_PRECISION) > (T_setpoint + T_deadband/2)) {
		// turn off the load if reverse_ON_OFF false otheweise turn on
		circuit_status = reverse_ON_OFF;
	}
	// if none of the conditions triggered, using the previous circuit_status
	return circuit_status;
}
} /* namespace gridballastcontroller */
