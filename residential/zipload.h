/** $Id: zipload.h 4738 2014-07-03 00:55:39Z dchassin $
	Copyright (C) 2009 Battelle Memorial Institute
	@file zipload.h
	@addtogroup ZIPload
	@ingroup residential

 @{
 **/

#ifndef _ZIPLOAD_H
#define _ZIPLOAD_H
#include "residential.h"
#include "residential_enduse.h"
#include "gridballastcontroller.h"
#include <vector>

using namespace std;
using std::vector;

class ZIPload : public residential_enduse
{
public:
	double base_power;			///< Base real power of the system
	double power_pf;			///< power factor of constant power load
	double current_pf;			///< power factor of constant current load
	double impedance_pf;		///< power factor of constant impedance load
	bool is_240;				///< load connected at 220 
	double breaker_val;			///< Amperage limit for connected breaker
	complex actual_power;		///< Actual load after adjusted for voltage factors

	bool demand_response_mode;	///< Activates equilibrium dynamic representation of demand response 
	int64 N;					///< Number of devices to model - base power is per device 
	int16 L;					///< Range of the thermostat's control operation 
	double N_off;				///< Number of devices that are off 
	double N_on;					///< Number of devices that are on 
	double noff;				///< Density of devices that are off per unit of temperature 
	double non;					///< Density of devices that are on per unit of temperature 
	double roff;				///< rate at which devices cool down 
	double ron;					///< rate at which devices heat up 
	double t; 					///< total cycle time of a thermostatic device 
	double toff;				///< total off time of device 
	double ton;					///< total on time of device 
	int16 x;					///< temperature of the device's controlled media (eg air temp or water temp) 
	double phi;					///< duty cycle of the device 
	double PHI;					///< diversity of a population of devices 
	double eta;					///< consumer demand rate that prematurely turns on a device or population 
	double rho;					///< effect rate at which devices heats up or cools down under consumer demand 
	double nominal_power;
	int64 next_time, last_time; ///< used to keep track of time in "special" modes - DR, duty-cycle
	double duty_cycle;			///< effective duty cycle of device
	double last_duty_cycle;
	double period;				///< period at which duty cycle is applied
	double phase;				///< phase of the duty cycle in terms of 0-1
	double multiplier;			///< static multiplier to modify base power ( load = base_power * multiplier )
	double recovery_duty_cycle; ///< duty cycle during recovery interval
	bool heatgain_only;			///< Activates a heat only mode - no load electric load is assigned to the load

	// define the controller
	gridballastcontroller::gridballastcontroller gbcontroller;

	// frequency control variables
	double main_frequency;			// grid frequency accessed from the parent triplex node
	double measured_frequency; 		// grid frequency from measurement at each time t
	double freq_lowlimit;			// lower tripping limit of the frequency
	double freq_uplimit;			// upper tripping limit of the frequency
	// voltage control variables
	double measured_voltage;		// grid voltage from measurement at each time t
	double volt_lowlimit;			// lower tripping limit of the voltage
	double volt_uplimit;			// upper tripping limit of the voltage

	// we use this variable to toggle with/without frequency/voltage control
	bool enable_freq_control;
	bool enable_volt_control;

	bool prev_status;
	bool circuit_status;			// True - ON; False - OFF, the returned variable to decide ON/OFF status
	bool temp_status;

	// jitter function,  jitter is enabled by default once freq/volt controller is enabled, can be set to 0
	// we give the same parameter to freq/volt jitters
	double average_delay_time;   			// in seconds, parameter for the uniform distribution
	// freq jitter variabels
	int freq_jitter_counter;			  	// a jitter counter generated based on Uniform Distribution each time the frequency violation happened
	bool freq_circuit_status_after_delay;  	// boolen to keep track of the circuit status after certain delay
	bool freq_jitter_toggler;				// indicate whether the freq jitter is activated (freq_jitter_counter>0)
	// volt jitter variables
	int volt_jitter_counter;			  	// a jitter counter generated based on Uniform Distribution each time the voltage violation happened
	bool volt_circuit_status_after_delay;  	// boolen to keep track of the circuit status after certain delay
	bool volt_jitter_toggler;				// indicate whether the volt jitter is activated (volt_jitter_counter>0)
	int temp_cnt;

	// force the circuit to be ON/OFF, we don't need it here
	int enable_lock;
	int lock_STATUS;

	// controller_priority
	int controller_priority;
	bool status_confirmed;
	vector<pair<int,int> > controller_array;

	typedef struct {
		double *on;
		double *off;
		int16 nbins;
	} DRMODEL;

	DRMODEL drm;
	DRMODEL previous_drm;			///< structures to save drm population and previous population

private:
	int first_pass;

public:
	static CLASS *oclass, *pclass;
	static ZIPload *defaults;

	ZIPload(MODULE *module);
	~ZIPload();
	int create();
	int init(OBJECT *parent);
	int isa(char *classname);
	bool get_status(int controller_number);
	TIMESTAMP sync(TIMESTAMP t0, TIMESTAMP t1);

};

#endif // _ZIPLOAD_H

/**@}**/
