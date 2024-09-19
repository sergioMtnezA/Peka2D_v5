/**
 * @file peka2dwindows.c
 * @author S. Martinez-Aranda
 * @date July 2024
 * @brief Execution file for Peka2D-v5.0.
 */

#ifdef _WIN32
	#if APP_MODE == 1
	#using <mscorlib.dll>
	#using <System.dll>
	#using <System.Windows.Forms.dll>
	#using <System.Drawing.dll>
	#if SIMULATEGPU
	#include "dev_CUDA_app\p2dGUI.h"
	#else
	#include "dev\p2dGUI.h"
	#endif

	using namespace System;
	using namespace System::Collections::Generic;
	using namespace System::IO;
	using namespace System::Windows::Forms;
	using namespace System::Drawing;
	using namespace dev;
	#endif
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "lib/define.h"
#include "lib/mainKernel.h"
//#include "lib/read_io.h"
//#include "lib/hydronia.h"

Peka2D_Setup *pksetup;

t_parameters spar;
t_message *msg;
t_mesh *mesh;

int main (int argc, char * argv[]) {
	
    char arguments[2][1024];
	sprintf(arguments[0],"C:\\svnPeka2D\\Casos\\YubaOBCP\\");
	sprintf(arguments[1],"YubaP");

	pksetup=(Peka2D_Setup*)malloc(sizeof(Peka2D_Setup));

	msg=(t_message*)malloc(sizeof(t_message));
	mesh=(t_mesh*)malloc(sizeof(t_mesh));

	#ifdef _WIN32
		#if APP_MODE == 1
			Application::Run(gcnew dev::MyForm(arguments[0],arguments[1])); //,msg
		#else
			runMainKernel(argc,argv);
		#endif
	#else
		runMainKernel(argc,argv);
	#endif

	return 0;

}

