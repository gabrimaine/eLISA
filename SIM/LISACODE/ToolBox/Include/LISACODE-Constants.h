// $Id:  $
/*
 *  LISACODE-Constants.h
 *  V 2.0
 *
 *  Created on 19/09/05 by A.Petiteau (APC)
 *  Last modification on 26/01/10 by A.Petiteau (AEI) 
 */

/**\defgroup ToolBox ToolBox  
 * This group (directory ToolBox) include the tools used 
 * by the rest of LISACode
 *\{
 */

/** \defgroup Constants Constants  
 *  Physical constants, reference values and unit conversions.
 *  \author A. Petiteau
 *  \version 2.0
 *  \date 14/02/11
 *	\{
 *
 */

#ifndef __LCCONSTANTS_H
#define __LCCONSTANTS_H

#include <math.h>
#include <float.h>
#include <cstring>


namespace LC{
	
	// ***********
	// * Version *
	// ***********
	
	/**\brief Version of simulator */
	const char LCVersion[] = "2.0.beta6";
	
	/**\brief Date of last modification */
	const char DateOfLastUpdate[] = "30/08/2013"; 
	
	// ************************
	// *  Physical constants  *
	// ************************
	
	/**\brief  Light's speed in \f$m\cdot s^{-1}\f$.*/
	const double c_SI = 299792458; // Light's speed in m.s-1
	
	/**\brief Gravitational constant in \f$m^3\cdot Kg^{-1} \cdot s^{-2}\f$.*/
	const double G_SI = 6.67259e-11; // Gravitationnal constant in m3.kg-1.s-2
	
	/**\brief Boltzmann's constant in \f$joule\cdot kelvin^{-1}\f$.*/
	const double k_SI = 1.381e-23; // Boltzman constant in J.K-1
	
	/**\brief Planck's constant in \f$joule\cdot second\f$.*/
	const double h_SI = 6.62620e-34; // Planck constant in J.s
	
	/**\brief Planck's constant divided by \f$2 \cdot \pi\f$ in \f$ joule \cdot second\f$.*/
	const double hb_SI = h_SI/(2*M_PI); // Quantum de moment cinetique en J.s
	
	/**\brief Stefan's constant in \f$Watt\cdot meter^{-2} \cdot kelvin^{-1}\f$.*/
	const double S_SI = 5.671e-8; // Stefan constant in W.m-2.K-4
	
	/**\brief Permeability of free space or magnetic constant \f$\mu_0\f$ 
	 * in \f$Henry \cdot meter^{-1}\f$.*/
	const double mu0_SI = 4.e-7*M_PI; // Permeabilite magnetique du vide en H.m-1
	
	/**\brief Permittivity of free space or permittivity \f$\epsilon_0 \f$ 
	 * in \f$Farad \cdot meter^{-1}\f$.*/
	const double eps0_SI = 1.0/(4.e-7*M_PI*c_SI*c_SI); // Permittivite du vide en F.m-1
	
	/**\brief Avogadro's number in \f$mol^{-1}\f$.*/
	const double Na_SI = 6.02252e-27; // Avogadro number in mol-1
	
	/**\brief Hubble's constant in CGS. */
	const double H0_cgs = 0.7*3.24e-18; // Hubble constant in cgs 
	
	/**\brief Euler's constant.*/
	const double CE_RG = 0.57721566490153286060651209008240243104215933593994; // Constante d'Euler
	
	
	// ***********************
	// *  References values  *
	// ***********************
	
	/**\brief Electron rest mass in Kg.*/
	const double me_SI = 9.1091e-31; 
	
	/**\brief Proton rest mass in Kg.*/
	const double mp_SI = 1.6726e-27; 
	
	/**\brief Neutron rest mass in Kg.*/
	const double mn_SI = 1.6748e-27; 
	
	/**\brief Sun's mass in Kg.*/
	const double MS_SI = 1.9889e30; 
	
	/**\brief Sun's mass in Seconds.*/
	const double TSUN = 4.92549232189886339689643862e-6;
	
	/**\brief Sun's radius in meter.*/
	const double RS_SI = 6.95e8;
	
	/**\brief Sun's energy flux in Watt.*/
	const double LS_SI = 3.83e26; 
	
	/**\brief Sidereal year in seconds. */
	const double YrSi_SI = 3.15581498e7; 
	
	/**\brief Astronomical year in seconds. */
	const double Yr_SI = 3.1556926e7; 
	
	/**\brief Day duration in seconds.*/
	const double Dy_SI = 24.0*3600.0; 
	
	/**\brief Half of the Schwarzhild radius in \f$\frac{GM}{c^2}\f$.
	 *
	 * In Schwarzhild radius \f$G\f$ is the gravitational constant, \f$m\f$ 
	 * is the mass of the black hole, and \f$c\f$ is the speed of light.*/
	const double RSchw = 1.47664e3; // Half of the Schwartzchild radius GM/c^2
	
	/*!\brief Post-Newtonian constant. */
	const double gamma_u = 1. ;
	
	/*!\brief gamma : value used in PN ??? */
	const double PN_GAMMA =     0.5772156649015328606065120900824024L; 
	
	// ***********************
	// *  Units conversions  *
	// ***********************
	
	
	/**\brief Astronomical unit in meters.*/
	const double au_m = 1.49597870660e11; 
	
	/**\brief Astronomical unit in seconds.*/
	const double au_s = au_m/c_SI;
	
	/**\brief Light year in meters (\f$ 9.460730472580800 \cdot 10^{15}\f$).*/
	const double ly_m = c_SI*365.25*24.0*3600.0; // Light year in m (9.460730472580800e15)
	
	/**\brief Light year in astronomical units (\f$63240.17695575401\f$).*/ 
	const double ly_au = ly_m/au_m; // Light year in astronomic unit (63240.17695575401)
	
	/**\brief Parsec in astronomical unit (\f$206265\f$).*/
	const double pc_au = M_PI/(3600.0*180.0); // Parsec in astronomic unit (206265) 
	
	/**\brief Parsec in meters (\f$3.086\cdot 10^{16}\f$).*/
	const double pc_m = pc_au*au_m; //Parsec in meter (3.086e16) #### TO BE CHECK IT'S STRANGE
	
	/**\brief Parsec in light year  (\f$3.262\f$).*/
	const double pc_ly = pc_au/ly_au; //Parsec in light year (3.262)
	
	/**\brief Number of meters in a kiloparsec (kpc).*/
	//const double kpc_m = 3.0856675807e19; // Number of meters in a kpc
	const double kpc_m = 3.08568025e19; // Number of meters in a kpc
	
	/*!\brief Number of second in a kiloparsec (kpc) : kpc in m / c .*/
	const double kpc_s = 1.02927214e11; // Number of seconds in a kpc
	
	/** \brief Angular velocity of Earth */
	const double omegaYr = 2.*M_PI/Yr_SI; 
	
	
	
	// ************************
	// *  Computer constants  *
	// ************************
	
	/** \brief Acceptable precision error on doubles.*/
	const double PRECISION = 100.0*DBL_EPSILON;
	
	/** \brief Maximal limit of double : a value inferior is an error ... or something specific */
	const double DBLMAXLIM = 1.1e300; 
	/** \brief Minimal limit of double : an inferior value is an error ... or something specific */
	const double DBLMINLIM = -1.1e300;
	/** \brief Maximal allowed value for a double  */
	const double DBLMAXALLOW = 1.0e300;
	/** \brief Minimal allowed value for a double  */
	const double DBLMINALLOW = -1.0e300;
    
    /**\brief Maximal order for Lagrange interpolation (use for storage of coefficients)*/
    const int ORDERMAXLAG = 30;
	
    // ************
	// *  Others  *
	// ************

    /**\brief Date of last modification */
	const char xmlind[] = "    "; 
}




#endif // __LCCONSTANTS_H

/**\}*/ //end of group constants
/**\}*/ //end of group toolbox 

//end of LISACODE-Constants.h
