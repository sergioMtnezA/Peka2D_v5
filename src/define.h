#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#define LINUX_COMPILATION 1

#if LINUX_COMPILATION==1

	/** Selection of the solver
	 *  1 SWE-Roe 
	*/
	#ifndef SOLVER
		#define SOLVER 1
	#endif

	/** Enables GPU compilation 
	 * 0 Disabled
	 * 1 Enabled
	*/
	#ifndef SET_SIMGPU
		#define SET_SIMGPU 1
	#endif	

	#if SET_SIMGPU
		#define SET_OPENMP 0 //always disableb for GPU compilation
	#else
		/** Enables OMP for CPU compilation.
		 * 0 Sequential compilation 
		 * NCORE parallelized compilation
		 */
		#define SET_OPENMP 4		
	#endif

	/** Enables solutes computation 
	 * 0 Disabled
	 * 1 Enabled
	*/
	#ifndef SET_SOLUTE
		#define SET_SOLUTE 1
		#define SET_MULTILAYER 0
	#endif		

#endif

/** Active elemens arrangement
 * 0 Disabled
 * 1 Enabled
 */
#ifndef RECONSTRUC_ACTIVE
	#define RECONSTRUC_ACTIVE 1
	#define nIterArrangeActElem 1000
#endif

/** Enables the generation of VTK files in binary codification
 * 0 Disabled
 * 1 Enabled
 */
#define BINARY_VTK 1

/**  Indicates the compilation mode. It affects, for instance, to the way
 * the messages are displayed and stored as well as the generation of
 * some outfiles that are not useful when using the GUI
 * APP_MODE == 1 indicates that the code will be genrated using a GUI
 * APP_MODE == 0 indicates that the code is compiled for console
 */
#define APP_MODE 0

/** Expor DLL mode
 * 
*/
#define EXPORT_DLL_MODE 0
#if EXPORT_DLL_MODE==1
	#ifndef EXPORT_DLL
		#define EXPORT_DLL extern "C" __declspec(dllexport)
	#endif
#else
	#define EXPORT_DLL
#endif


#define FT2MFACT (1/3.2808)
#define M2FTFACT (3.2808)
#define LB2KGFACT (1/2.2046)
#define KG2LBFACT (2.2046)
#define FT32M3FACT ((1/3.2808)*(1/3.2808)*(1/3.2808))
#define M32FT3FACT (3.2808*3.2808*3.2808)
#define PSI2PASCFACT (6894.75729)
#define FT22M2FACT ((1/3.2808)*(1/3.2808))
#define BTU2JOULE (1055.056)



//Program's config
#define WRITE_STATE 1
#define WRITE_MASS  1
#define threadsPerBlock 256
#define threadsPerOBC 256


//Other definitions
#define IBC_WEIR  1
#define IBC_GATE  2
#define IBC_RATING_TABLE 3
#define IBC_DAMBREACH 4

//Constants and flags
#ifndef SIGN
	#define SIGN(x) ((x > 0.) - (x < 0.))
#endif

#ifndef MIN
	#define MIN(x,y) (((x)<(y))?(x):(y))
#endif

#ifndef MAX
	#define MAX(x,y) (((x)>(y))?(x):(y))
#endif

#ifndef MINMAX
	#define MINMAX(x,y,z) (((z) < (x)) ? (x) : (((z) > (y)) ? (y) : (z)))
#endif

#if defined(_WIN32) || defined(_WIN64)
	#ifndef WIN32
		#define WIN32
	#endif
	#ifndef CBRT
		#define CBRT(x) (pow(x,1./3.))
	#endif
	#pragma warning (disable:4996)
	#define snprintf sprintf_s
#else
	#ifndef CBRT
		#define CBRT(x) (cbrt(x))
	#endif
#endif

#define MAXHMIN 0.01
#define MINHMIN 0.0005
#define MAX_ANGLE_HYDRONIA_MESH 15.0

#define _X_ 0
#define _Y_ 1

#define TOL3 1e-3
#define TOL4 1e-4
#define TOL5 1e-5
#define TOL6 1e-6
#define TOL7 1e-7
#define TOL9 1e-9
#define TOL12 1e-12
#define TOL14 1e-14
#define TOL16 1e-16

#ifndef _PI_
	#define _PI_ 3.1415926535897
#endif
#ifndef _g_
	#define _g_ 9.80665
#endif
#ifndef _e_
	#define _e_ 2.71828
#endif
#define _raizg_ 3.1315571206669694 //sqrt(_g_)
#define _atan1_ 0.785398163397448 //atan(1.0)
#define _dospi_ 6.283185307179584 //2*pi

#define _rhow_ 1000.0	// Water density (kg/m3)
#define _nu_	0.00000151	// Water kinematic viscosity (m2/s)
#define _cp_ 4182.0

// NOTIFICATIONS
#define MSG_ERROR -10
#define MSG_WARN -1
#define MSG_L0 0
#define MSG_L1 1
#define MSG_L2 2
#define MSG_L3 3

//DEBUGS
#define DISPLAY 1

// colors
#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"
// bold
#define BOLD  "\x1B[1m"
#define BRED  "\x1B[1;31m"
#define BGRN  "\x1B[1;32m"
#define BYEL  "\x1B[1;33m"
#define BBLU  "\x1B[1;34m"
#define BMAG  "\x1B[1;35m"
#define BCYN  "\x1B[1;36m"
#define BWHT  "\x1B[1;37m"

#ifndef MSGERROR
	#if defined(_WIN32)
	#define MSGERROR " [FAIL] "
	#else
	#define MSGERROR " \033[1;31m[FAIL]\033[0m "
	#endif
#endif

#ifndef MSGOK
	#if defined(_WIN32)
	#define MSGOK " [_OK_] "
	#else
	#define MSGOK " \033[1;32m[_OK_]\033[0m "
	#endif
#endif

#ifndef MSGINFO
	#if defined(_WIN32)
	#define MSGINFO " [INFO] "
	#else
	#define MSGINFO " \033[1;35m__ii__\033[0m "
	#endif
#endif

#ifndef MSGLOAD
	#if defined(_WIN32)
	#define MSGLOAD " [LOAD] "
	#else
	#define MSGLOAD " \033[1;36m[LOAD]\033[0m "
	#endif
#endif

#ifndef MSGMEMO
	#if defined(_WIN32)
	#define MSGMEMO " [MEMO] "
	#else
	#define MSGMEMO " \033[1;34m[MEMO]\033[0m "
	#endif
#endif

#ifndef MSGREAD
	#if defined(_WIN32)
	#define MSGREAD " [READ] "
	#else
	#define MSGREAD " \033[1;34m[READ]\033[0m "
	#endif
#endif

#ifndef MSGGPU
	#if defined(_WIN32)
	#define MSGGPU " [GPU_] "
	#else
	#define MSGGPU " \033[1;36m[CUDA]\033[0m "
	#endif
#endif

#ifndef MSGEXC
	#if defined(_WIN32)
	#define MSGEXC " [EXC_] "
	#else
	#define MSGEXC " \033[1;34m[EXEC]\033[0m "
	#endif
#endif

#ifndef MSGWARN
	#if defined(_WIN32)
	#define MSGWARN " [WARN] "
	#else
	#define MSGWARN " \033[1;33m[WARN]\033[0m "
	#endif
#endif


#define DEBUGMSG printf("file '%s' in line %i\n",\
	__FILE__, __LINE__ );

#define CUT_CHECK_ERROR(errorMessage) do {\
	cudaThreadSynchronize();\
	cudaError_t err = cudaGetLastError();\
		if( cudaSuccess != err) {\
			fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",\
					errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) );\
			exit(EXIT_FAILURE);\
		} } while (0)


#define DUMPDELTAT 0
#define DEBUG 0
#define DEBUG_READ_BOUNDS 0
#define DEBUG_READ_SECTIONS 0
#define DEBUG_READ_SOURCES 0
#define WRITE_FROUDE 0


//flux limiters
#define FLUX_LIMITER 0 //0--> minMod    1--> Van Albada    2-->sweby beta=1.5   MinMod recommended

#define functionPhi0(y) MAX(0.0,MIN(1.0,y)) //minMod
#define functionPhi1(y) (y*y + y)/(y*y+1) //van Albada
#define functionPhi2(y) MAX(0.0,(MIN(1.5*y,1.0),MIN(y,1.5))) //sweby

#if FLUX_LIMITER==0
	#define functionPhi(y) functionPhi0(y)
#endif

#if FLUX_LIMITER==1
	#define functionPhi(y) functionPhi1(y)
#endif

#if FLUX_LIMITER==2
	#define functionPhi(y) functionPhi2(y)
#endif


// Hydrodynamic friction models
#define HFM_MANNING 0
#define HFM_CHEZY 1
#define HFM_DW 2
// By default the Manning model is used. You should choose the friction model at compilation time. You can do this when calling the makefile by specifying FRICTION=your_chosen_model_index
#ifndef HYDRO_FRICTION
    #define HYDRO_FRICTION HGM_MANNING
#endif
