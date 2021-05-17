// To use this library create a wave of parameters with the following form:
//	FWHM	// Gaussian FWHM
//	t0		// Offset on t, where the heaviside funtion steps
//	Other parameters required by the fit function

// EXAMPLE of parameter wave for the fit function Rise1expFall2exp:
//	FWHM
//	t0
//	tau_r1
//	I_f1
//	tau_f1
//	I_f2
//	tau_f2
//	y0_f
//	A
//		See the function Rise1expFall2exp to see what each parameter represents

// By changing the fit function in QuickConv it's possible to quickly test different
// function shapes, like 1 rising exp and 2 falling or 2 rising and 2 falling




#pragma rtGlobals=1		// Use modern global access method


/// <summary>Gaussian function g(t) with peak normalized to 1.
/// "param" is the fwhm, "t" the input wave on which the function is calculated</summary>
Function Gaussian(param,t)
	wave param;		// Parameter vector
	variable t;		// Time vector

	// Parameters of the gaussian
	variable fwhm = abs(param[0]);		// Gaussian FWHM

	//g(t) = exp[-2.773 (t/fwhm)^2], where -2.773 normalized the peak to 1
	variable g = exp(-2.773*(t/fwhm)^2);

	return(g)
End




/// <summary>Function f(t) = 1 rising exp * sum of 2 falling exp.
/// "param" is the wave of parameters, "t" the input wave on which the function is calculated</summary>
Function Rise1expFall2exp(param,t)
	wave param;		// Parameter vector
	variable t;		// Time vector

	// Parameters of f(t)
	variable t0 = param[1];			// Offset t
	variable tau_r1 = abs(param[2]);	// Rise time constant
	variable I_f1 = abs(param[3]);		// Fall intensity 1
	variable tau_f1 = abs(param[4]);	// Fall time constant 1
	variable I_f2 = abs(param[5]);		// Fall intensity 2
	variable tau_f2 = abs(param[6]);	// Fall time constant 2
	variable y0_f = param[7];		// Fall offset y
	variable A = param[8];			// Constant the whole expression	
	
	// f(t) =        (______rise exp 1_____) *          ( _______fall exp 1______  + _______fall exp 2______  )    
	variable f = A * (1-exp(-(t-t0)/tau_r1)) * ( y0_f + ( I_f1*exp(-(t-t0)/tau_f1) + I_f2*exp(-(t-t0)/tau_f2) ) );

	return (f);
End




/// <summary>Function f(t) = 1 rising exp * 1 falling exp.
/// "param" is the wave of parameters, "t" the input wave on which the function is calculated</summary>
Function Rise1expFall1exp(param,t)
	wave param;		// Parameter vector
	variable t;		// Time vector

	// Parameters of f(t)
	variable t0 = param[1];			// Offset t
	variable tau_r1 = abs(param[2]);	// Rise time constant
	variable I_f1 = abs(param[3]);		// Fall intensity 1
	variable tau_f1 = abs(param[4]);	// Fall time constant 1
	variable y0_f = param[5];		// Fall offset y
	variable A = param[6];			// Constant the whole expression
	
	// f(t) =        (______rise exp 1_____) *          ( _______fall exp 1_______ )    
	//variable f = A * (1-exp(-(t-t0)/tau_r1)) * ( y0_f + ( I_f1*exp(-(t-t0)/tau_f1) ) );
	
	variable f = A * (1-exp(-(t-t0)/tau_r1)) * ( y0_f + I_f1*exp(-(t-t0)/tau_f1)  );

	return (f);
End




/// <summary>Function f(t) = sum of 2 rising exp * 1 falling exp.
/// "param" is the wave of parameters, "t" the input wave on which the function is calculated</summary>
Function Rise2expFall1exp(param,t)
	wave param;		// Parameter vector
	variable t;		// Time vector

	// Parameters of f(t)
	variable t0 = param[1];			// Offset t
	variable I_r1 = abs(param[2]);		// Rise intensity 1
	variable tau_r1 = abs(param[3]);	// Rise time constant 1
	variable I_r2 = abs(param[4]);		// Rise intensity 2
	variable tau_r2 = abs(param[5]);	// Rise time constant 2
	variable I_f1 = abs(param[6]);		// Fall intensity 1
	variable tau_f1 = abs(param[7]);	// Fall time constant 1
	variable y0_f = param[8];		// Fall offset y
	variable A = param[9];			// Constant the whole expression	
	
	// f(t) =        (______rise exp 1_____) *          ( _______fall exp 1______  + _______fall exp 2______  )    
	variable f = A * ( (1-I_r1*exp(-(t-t0)/tau_r1)) +  (1-I_r2*exp(-(t-t0)/tau_r2))) * ( y0_f + I_f1*exp(-(t-t0)/tau_f1)  );

	return (f);
End





/// <summary>Convolutes the function f(t) with g(t).
/// Computes the integral of f(t')g(t-t')H(t-t0)dt', where H(t) is the Heaviside function.
/// "param" is the vector with all the parameters of both functions.
/// "t" the input wave on which the functions f(t) and g(t) are defined.</summary>
Function QuickConv(param, t)//, interval_start, interval_end, step)
	wave param;		// Parameter vector
	variable t;		// Time vector

	variable t0 = param[1];			// Offset on t

	variable interval_start = -400;	// Convolution interval start
	variable interval_end = 3600;	// Convolution interval end
	variable step = 2;				// Numeric integration step
	variable y; 		// Convolution variable
	variable conv = 0;	// Final result

	for(y = interval_start; y < interval_end; y += step)
		if (y < t0)		// This if statement is a Heaviside function
			continue;
		endif
		conv += Gaussian(param,t-y) * Rise1expFall2exp(param,y) * step;
	endfor

	return (conv);
End




/// <summary>Convolutes the function f(t) with the g(t) with IGOR built in function.
/// Computes the integral of f(t')g(t-t')dt'.
/// "param" is the vector with all the parameters of both functions</summary>
Function QuickIgorConv(param, t)
	wave param;		// Parameter vector
	wave t;			// Time wave

	// Duplicate t to make the wave for g and f the right size
	// /D creates a double precision vector, /O overwrites
	Duplicate/D/O t g
	Duplicate/D/O t f

	// Assign to f and g the functions you prefer
	f = Rise1expFall2exp(param,t)
	g = Gaussian(param,t)

	// Convolve with causality (/A) f and g. The result of the convolution is stored in g
	// Note: can't explain result
	Convolve/A f g
End

EndMacro
