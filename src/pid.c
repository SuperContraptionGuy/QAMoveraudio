#include <stddef.h>
#include <math.h>

#include "pid.h"

static double pidClamp(pid_limits_t *limits, double input)
{
	if (limits == NULL)
	{
		return NAN;
	}

	// Limit is disabled if set to NAN

	if ((isnan(limits->lower) == 0) && (limits->lower > input))
	{
		// Clamp to lower limit
		return limits->lower;
	}
	else if ((isnan(limits->upper) == 0) && (limits->upper < input))
	{
		// Clamp to upper limit
		return limits->upper;
	}

	return input;
}

int pidInit(pid_t *pid)
{
	if (pid == NULL)
	{
		return -1;
	}

	pid->limits.lower = NAN;
	pid->limits.upper = NAN;

	pid->lastInput = 0.0;
	pid->lastOutput = 0.0;
	pid->dt = 0.0;

	pid->setpoint = 0.0;

	pid->gain.kP = 0.0;
	pid->gain.kI = 0.0;
	pid->gain.kD = 0.0;

	pid->terms.proportional = 0.0;
	pid->terms.integral = 0.0;
	pid->terms.derivative = 0.0;

	pid->error = 0.0;

	return 0;
}

int pidGain(pid_t *pid, double kP, double kI, double kD)
{
	if (pid == NULL)
	{
		return -1;
	}

	if (isnan(kP) == 0)
	{
		pid->gain.kP = kP;
	}

	if (isnan(kI) == 0)
	{
		pid->gain.kI = kI;
	}

	if (isnan(kD) == 0)
	{
		pid->gain.kD = kD;
	}

	return 0;
}

int pidLimits(pid_t *pid, double upper, double lower)
{
	if (pid == NULL)
	{
		return -1;
	}

	if (isnan(upper) == 0)
	{
		pid->limits.upper = upper;
	}

	if (isnan(lower) == 0)
	{
		pid->limits.lower = lower;
	}

	return 0;
}

int pidSetpoint(pid_t *pid, double setpoint, double dt)
{
	if (pid == NULL)
	{
		return -1;
	}

	if (isnan(setpoint) == 0)
	{
		pid->setpoint = setpoint;
	}

	if (isnan(dt) == 0)
	{
		pid->dt = dt;
	}

	return 0;
}

double pidUpdate(pid_t *pid, double input)
{
	if (pid == NULL)
	{
		return NAN;
	}

	if (isnan(input) != 0)
	{
		return NAN;
	}

	pid->error = pid->setpoint - input;

	// Compute input delta
	double dInput = input - pid->lastInput;

	// Compute proportional
	pid->terms.proportional = pid->gain.kP * pid->error;
	// Compute and clamp integral
	pid->terms.integral = pidClamp(&pid->limits, pid->gain.kI * pid->error * pid->dt);
	// Compute derivative
	pid->terms.derivative = (pid->gain.kD * -1.0) * (dInput / pid->dt);

	// Save the output value
	pid->lastOutput = pidClamp(
		&pid->limits, 
		(
			pid->terms.proportional +
			pid->terms.integral +
			pid->terms.derivative
		)
	);

	// Save the input value
	pid->lastInput = input;

	return pid->lastOutput;
}
