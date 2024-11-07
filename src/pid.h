#pragma once

typedef struct
{
	double lower;
	double upper;
} pid_limits_t;

typedef struct
{
	double kP;
	double kI;
	double kD;
} pid_gain_t;

typedef struct
{
	double proportional;
	double integral;
	double derivative;
} pid_term_t;

typedef struct
{
	// PID clamping limits
	pid_limits_t limits;

	double lastInput;
	double lastOutput;
	double dt;

	double setpoint;

	pid_gain_t gain;

	pid_term_t terms;

	double error;
} pid_t;

int pidInit(pid_t *pid);
int pidGain(pid_t *pid, double kP, double kI, double kD);
int pidLimits(pid_t *pid, double upper, double lower);
int pidSetpoint(pid_t *pid, double setpoint, double dt);
double pidUpdate(pid_t *pid, double input);
