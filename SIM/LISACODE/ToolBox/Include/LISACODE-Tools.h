// $Id:  $
/*
 *  LISACODE-Tools.h
 *  V 2.0
 *
 *  Created on 28/05/09 by  Antoine Petiteau (AEI)
 *  Last modification on 15/02/10 by Antoine Petiteau (AEI)
 *
 */

/** \ingroup ToolBox 
 * \defgroup LCTools_  Class LCTools
 * (See class #LCTools for a detailed description)
 * \{
 */

#ifndef __LCTOOLS_H
#define __LCTOOLS_H

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <unistd.h>
#include <complex>
#include <time.h>
#include "LISACODE-Constants.h"
#include "ezxml.h"
#include "randlib.h"


#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define SWAP2(a,b) {temp=(a);(a)=(b);(b)=temp;}
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define Cout (*MT->o)
#define Coutm (*MT.o)
#define Endl std::endl

typedef std::complex<double> dcomplex;

const dcomplex I(0.0,1.0);


/*! \brief Type of endianess */
enum TypeEncoding{
	BIGENDIAN,
	LITTLEENDIAN,
	NOENDIAN
};


/*! \brief Type of file */
enum TypeFile{
	ASCII,
	BINARY,
	XML,
	UNDEFINED
};


/** \brief Class for managing memory and providing a number of tools.
 * \author A. Petiteau
 * \version 2.0
 * \date 14/02/2011
 *
 * This class is a toolbox which is known and used by all the others classes.
 * The memory is managed using the fftw allocator.
 * This class manages the screen display.
 * 
 *
 */
class LCTools
	{
	private:
		
		/** \brief Total memories allocated using LCTools::AllocMemory function */
		size_t _MemAllocSize_;
		
		/** \brief Maximal total memories which has been allocated */
		size_t _MemAllocSizeMax_;
		
		
		/** \brief Type of log display :
		 *	\arg 0  NULL
		 *	\arg 1	cout
		 *	\arg 2  file
		 */
		int DispLog;
		
		/** \brief Pointer display output file */
		std::ofstream * fdout;
		
		/** Verbose : Display full details */
		bool DispDetails;
		
		/** Seed of random generator */
		long RandSeed;
		
		/** Local encoding */
		TypeEncoding LocalEndianness;
        
        
        /**********  For Lagrange interpolation  **********/
        
        /*! \brief Table containing Lagrange coefficient (filled at the first used). */
        double * LagPolDen[LC::ORDERMAXLAG];
		
		
		/**********  For FFTW  **********/
		/******  Forward  ******/
		
		/*! \brief fftw forward plan */ 
		fftw_plan FTFwdPlan;
		
		/*! \brief Input data of forward fftw (typically time data) */ 
		double *FTFwdtimeData;
		
		/*! \brief Number of input data of forward fftw */
		int FTFwdNbtDat;
		
		/*! \brief Output data of forward fftw (typically frequency data) */ 
		dcomplex *FTFwdfreqData;
		
		/*! \brief Number of output data of forward fftw */
		int FTFwdNbfDat;
		
		/*! \brief Preision of forward fftw : if true FFTW_MEASURE else FFTW_ESTIMATE */
		bool FTFwdBestFFTW;
		
		/******  Backward  ******/
		
		/*! \brief fftw backward plan */ 
		fftw_plan FTBckPlan;
		
		/*! \brief Input data of backward fftw (typically frequency data) */ 
		dcomplex *FTBckfreqData;
		
		/*! \brief Number of input data of backward fftw */
		int FTBckNbtDat;
		
		/*! \brief Output data of backward fftw (typically time data) */
		double *FTBcktimeData;
		
		/*! \brief Number of output data of backward fftw */
		int FTBckNbfDat;
		
		/*! \brief Preision of backward fftw : if true FFTW_MEASURE else FFTW_ESTIMATE */
		bool FTBckBestFFTW; 
		
		/**********   MPI  ********** /
		
		/** \brief Process identification */
		int mpi_id;
		
		/** \brief Number of process */
		int mpi_np;
		
	public:
		
		/** \brief Pointer on output log */
		std::ostream * o; 
		
		/** \brief True if it's a tool which has been localy created */
		bool LocTools;
		
		/**********  Constructor  **********/
		
		/** \brief Default constructor. */
		LCTools();
		
		/** \brief Destructor. */
		~LCTools(); 
		
		
		/**********  Memory  **********/
		
		/** \brief Allocate memory. 
		 * @param[in] size  Size of memory to be allocatted */
		void* AllocMemory(size_t size);
		
		/** \brief Allocate memory. 
		 * @param[in] size  Old pointer 
		 * @param[in] size  Old size of memory allocatted 
		 * @param[in] size  New size of memory to be allocatted */
		void* ReAllocMemory(void * oldptr, size_t oldsize, size_t newsize);
		
		/** \brief Free memory. 
		 * @param[in] ptr   Pointer on the variable to be removed
		 * @param[in] size  Size of the variable to be removed */
		void Free(void *ptr, size_t size);
		
		/** \brief Display the current use of memory. */
		void MemDisplay() const;
		
		
		/**********  Display  **********/
		
		/** \brief Return true if log output should be displayed. */
		bool Disp() const {return((DispLog==0)?(false):(true));};
		
		/** \brief Return true if log output should be displayed. */
		int getDispLog() const {return(DispLog);};
		
		/** \brief Set the fact that log output should not be displayed. */
		void unsetDisp();
		
		/** \brief Set the fact that log output should be displayed in the standard log output (std::cout). */
		void setDispLog();
		
		/** \brief Set output display in file OutFile 
		 * @param[in] OutFile  Name of file which will contain standard log output
		 * @param[in] NewFile  If true a new file is created overwritting the old one */
		void setDispInFile(char * OutFile, bool NewFile);
		
		/** \brief Set output display in file OutFile+cpuId 
		 * @param[in] OutFile Base name of file which will contain standard log output
		 * @param[in] NewFile  If true a new file is created overwritting the old one */
		void setDispInFileCoreid(char * OutFile, bool NewFile);
		
		/** \brief Return true if details can be displayed */
		bool DispDet();
		
		/** \brief Set (swicth on) the display of full details (verbose) */
		void setDispDetails() {DispDetails = true; };
		
		/** \brief Unset (swicth off) the display of full details (verbose) */
		void unsetDispDetails() {DispDetails = false; };
		
		/**********  Random generator  **********/
		
		/** \brief Set seed of random generator 
		 * @param[in] RandSeed_n  New seed of random generator */
		void setRandSeed(long RandSeed_n);
		
		/** Return the seed used for the random generator */
		long getRandSeed() const {return(RandSeed);};
		
		/** \brief Draw randomly a real value from a uniform distribution in a range \f$ [ x_{min} , x_{max} ] \f$  
		 * @param[in] xmin  Minimal value, \f$ x_{min} \f$
		 * @param[in] xmax  Maximal value, \f$ x_{max} \f$ */
		double RandUniform(double xmin, double xmax);
		
		/** \brief Draw randomly a real value from a normal distribution  \f$ \mathcal{N} ( x_{0} , \sigma_x ) \f$  
		 * @param[in] center  Mean value (center of the gaussian), \f$ x_{0} \f$
		 * @param[in] sigma   Standard deviation (sigma of guassian), \f$ \sigma_x \f$ */
		double RandGaussian(double center, double sigma);
		
		
		/**********  Double/integer/complex tools  **********/
		
		/** \brief Return lower integer corresponding to the truncated value */
		inline int ifloor(double f) {return((int)(f));}; 
		
		/** \brief Return higher integer corresponding to the truncated value */
		inline int iceil(double f) {return((int)(f+1.0-LC::PRECISION));};
		
		/** \brief Return true if the two real are identic relatively to the numerical precision */
		inline bool deq(double a, double b) { return( (fabs(a-b) < LC::PRECISION) ); };
		
		
		/**********  Mathematical tools  **********/
		/*! \brief Linear equation solution by Gauss-Jordan elimination
		 * a[1..n][1..n] is the input matrix.
		 * b[1..m][1..m] is input containing the m right-hand side std::vectors.
		 * On output, a is replaced by its matrix inverse, and b is replaced by the corresponding set of solution std::vectors.
		 */
		int Gaussj(double **a, int n, double **b, int m);
		
		
		/**********  Physic tools  **********/
		/*! \brief  Compute Hubble expansion parameter for the given parameters :
		 *	@param[in]	z	redshift
		 *	@param[in]	H0	Current value of Hubble expansion parameter
		 *	@param[in]	Omm	Omega matter
		 *	@param[in]	Oml	Omega lambda
		 *	@param[in]	w	(- 1 - omDE) : parameter of dark energy equation of state  
		 *	return	Hubble expansion parameter
		 */
		double Hz(double z, double H0, double Omm, double Oml, double w);
		
		/*! \brief  Compute luminosity distance (in kpc) for the given parameters :
		 *	@param[in]	z	redshift
		 *	@param[in]	H0	Current value of Hubble expansion parameter
		 *	@param[in]	Omm	Omega matter
		 *	@param[in]	Oml	Omega lambda
		 *	@param[in]	w	(- 1 - omDE) : parameter of dark energy equation of state 
		 *	return luminosity distance (in kpc)
		 */
		double DL(double z, double H0, double Omm, double Oml, double w);
		
		/*! \brief  Compute redshift for the given parameters :
		 *	@param[in]	DL	luminosity distance (in kpc)
		 *	@param[in]	H0	Current value of Hubble expansion parameter
		 *	@param[in]	Omm	Omega matter
		 *	@param[in]	Oml	Omega lambda
		 *	@param[in]	w	(- 1 - omDE) : parameter of dark energy equation of state 
		 *	return redshift
		 */
		double redshift(double DL, double H0, double Omm, double Oml, double w);
		
		
		/**********  Word tools  **********/
		
		/** \brief Convert word in upper case */
		char * uppercase(const char * motmin) const;
		
		/** \brief Convert word in lower case */
		char * lowercase(const char * MOTMAJ) const;
		
		/** \brief Convert letter in majuscule (upper case letter, capital) */
		void majuscule(char &lettre) const;
		
		/** \brief Convert letter in minuscule (lower case letter) */
		void minuscule(char &lettre) const;
		
		/** \brief Compare two word : Return true if the two words are identical */
		bool wcmp(const char * mot1, const char * mot2) const ;
		
		/** \brief Extract words in a string of chararcters 
		 * @param[in]  InString String to decompose,
		 * @param[out] Words    result, i.e. list of the words */
		void wextract(char * InString, std::vector<char*> &Words);
		
		
		/** \brief Extract words in a string of chararcters where the words are separated by coma ','
		 * @param[in]  InAll	String to decompose,
		 * @param[in]  ipIA		Position in the string to decompose,
		 * @param[out] OutName  List of extracted words,
		 * @param[out] NOutName Number of output
		 */
		void wextractcoma(const char * InString, int & ipIS, char ** & Words, int & NWords, int NCharWords);
		
		/** \brief Makes a copy of the string "orig", stripping all whitespace-type characters
		 * @param[in]  orig initial string
		 * @param[out] ret  result string : it is a private copy that must be deallocated with Free(xxx , strlen(xxx)+1 )
		 */
		 void stripcopy(const char * orig, char * &ret);
		
		/** \brief Test if a string char is numerical value or not */
		bool isValue(const char * Str);
		
		
		/** \brief Extract name from path */
		void pathExtractName(char * Path, char * & Name);
		
		
		/**********  Files and System encoding  **********/
		
		/*! \brief Return the type of encoding of the local system */
		TypeEncoding testbyteorder();
		
		/*! \brief Return true if ... */
		bool IsLocalBigEndian() {return(LocalEndianness==BIGENDIAN);};
		
		/*! \brief Convert the encoding of the value */
		void convertendianness(double *val);
		
		/*! \brief  */
		void DispEncoding(TypeEncoding Encoding, std::ostream * OutDisp);
		
		/*! \brief Return the type of encoding read in string char */
		TypeEncoding ReadEncoding(char * Str);
		
		/*! \brief Return the type of encoding read in constant string char */
		TypeEncoding ReadEncodingConst(const char * Str);
		
		/*! \brief Display the type of file */
		void DispFileType(TypeFile FType, std::ostream * OutDisp);
		
		/*! \brief Check if file with name build from directory path + file name exist */
		void CheckFile(char * fDir, char * fName);
		
		
		/**********  XML tools  **********/
		
		/*! \brief Return unit 
		 *	@param[in]	Input constant char string (read from parameter)
		 *	@param[out] Value of unit
		 *	return Description of unit
		 */
		char * gXMLUnit(const char In[], double & Fact);
		
		/*! \brief Return a time value read in XML file. Unit is Second. */
		double gXMLTime(ezxml_t param);
		
		/*! \brief Return an angle value read in XML file. Unit is radians. */
		double gXMLAngle(ezxml_t param);
		
		/*! \brief Return a length value read in XML file. Unit is meter. */
		double gXMLLength(ezxml_t param);
		
		/*! \brief Return an frequency value read in XML file. Unit is hertz. */
		double gXMLFrequency(ezxml_t param);
		
		/*! \brief Return an astronomic mass value read in XML file. Unit is SolarMass. */
		double gXMLAstroMass(ezxml_t param);
		
		/*! \brief Return an astronomic distance value read in XML file. Unit is KiloParsec. */
		double gXMLAstroDistance(ezxml_t param);
		
		/*!\brief Return time in convention ISO-8601 */
		char * TimeISO8601();
		
		
        /**********  For Lagrange interpolation  **********/
        
        /**  \brief Initialize the lagrange coefficient 
         * @param[in] order Order of Lagrange interpolation
         */
        void InitLagPolDen(int order);
        
        /**  \brief Return data value corresponding to reference \f$x\f$ using Lagrange interpolation at specified order 
         * @param[in] x     Reference value \f$x\f$ (function found the corresponding data)
         * @param[in] order Order of Lagrange interpolation
         */
        double InterLagrange(double * data, int Nmax, double xr, int order);
        
        
		
		/**********  For FFTW  **********/
		
		/******  Forward  ******/
		
		/*! \brief For preparing fftw forward 
		 *	@param[in]	NtDat_n Number of input data
		 */
		void FTPrepFwdPlan(int NtDat_n);
		
		/*! \brief Compute one forward fftw
		 *	@param[in]	tData	Input data (tipically time data)
		 *	@param[out]	fData	Output data (tipically frequency data)
		 *	@param[in]	NtDat_n Number of input data
		 */
		void FTMakeFwd(double * tData, dcomplex * &fData, int NtDat_n);
		
		/*! \brief Compute multiple forward fftw
		 *	@param[in]	tData	List of input data (tipically time data)
		 *	@param[out]	fData	List of output data (tipically frequency data)
		 *	@param[in]	NtDat_n Number of input data
		 *	@param[in]	NSig	Number of lists (number of signals)
		 */
		void FTMakeFwdMulti(double ** tData, dcomplex ** &fData, int NtDat_n, int NSig);
		
		void setBestFTFwd() {FTFwdBestFFTW = true;};
		void unsetBestFTFwd() {FTFwdBestFFTW = false;};
		
		/*****  Backward  *****/
		
		/*! \brief For preparing fftw forward 
		 *	@param[in]	NfDat_n Number of input data
		 */
		void FTPrepBckPlan(int NfDat_n);
		
		/*! \brief Compute one backward fftw
		 *	@param[in]	fData	Input data (tipically frequency data)
		 *	@param[out]	tData	Output data (tipically time data)
		 *	@param[in]	NfDat_n Number of input data
		 */
		void FTMakeBck(dcomplex * fData, double * &tData, int NfDat_n);
		
		void setBestFTBck() {FTBckBestFFTW = true;};
		void unsetBestFTBck() {FTBckBestFFTW = false;};
		
		/***** Others *****/
		int  getNfFTreal(int NtDat) {return(NtDat/2+1);};
		int  getNtFTreal(int NfDat) {return(2*(NfDat-1));};
		
		
		/**********  For Power Spectral Density  **********/
		/*! \brief Compute power spectral density :
		 *	@param[in]	tDat		Input data : time data
		 *	@param[in]	NtDat		Number of input 
		 *	@param[in]	dt			Time step
		 *	@param[in]	TypeConvWin	Type of window : 1 : Bartlett ,  2 : Square ,  3 : Welch , 4 : Hann , 5 : Sinus
		 *	@param[in]	NSeg		Number segments : number of meaning  (typically 10)
		 *	@param[out] PSDDat		Output data : power spectral density
		 *	@param[out] NDatPSD		Number of output data : power spectral density
		 *	@param[out] df			Frequency step
		 *	@param[in]	AllocPSD	True if need to allocate the psd inside
		 */
		void PSD(double * tDat, int NtDat, double dt, int TypeConvWin, int NSeg, double * & PSDDat, int & NDatPSD, double & df, bool AllocPSD);
		
		
		/*! \brief Apply a convolution window  **********/
		double AppConvWin(int j, double a, double b, int TypeConvWin);
		
		
		/**********  MPI tools : Basic tools  **********/
		
		/** \brief MPI : Set id of the process and the number of process
		 * @param[in] mpi_id_n  Id of the process (core),
		 * @param[in] mpi_np_n  Number of processes (cores). */
		void MPIset(int mpi_id_n, int mpi_np_n);
		
		/** \brief  MPI : Return true if it's the master process  */
		bool MPImaster();
		
		/** \brief  MPI : Return id of local process */
		int MPIid() {return(mpi_id);};
		
		/** \brief  MPI : Return number of processes */
		int MPInbproc() {return(mpi_np);};
		
	};

#endif //__LCTOOLS_H

/**\}*/

// end of LISACODE-Tools.h
