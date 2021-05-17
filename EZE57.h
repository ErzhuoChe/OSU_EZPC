#pragma once
#include "E57_io.h"
#include "Las_io.h"
#include "EZPC_util.h"

namespace EZE57
{
	void E57_to_LAS(char** input_list, int i_input_first, int i_input_last, bool LAS_or_LAZ, bool report);
	void Extract_E57_Info(char** input_list, int i_input_first, int i_input_last);
	void Split_E57_Scan(char** input_list, int i_input_first, int i_input_last);
};