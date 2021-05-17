#pragma once
#include "LAS_io.h"
#include "EZPC_util.h"

namespace EZLAS
{
	void Remove_Invalid_RGB(char** input_list, int i_input_first, int i_input_last);
	void Merge_Split_SourceID(char** input_list, int i_input_first, int i_input_last);
	void Merge_Split_Classification(char** input_list, int i_input_first, int i_input_last);
	void Meters_to_Feet(char** input_list, int i_input_first, int i_input_last);
	void Feet_to_Meters(char** input_list, int i_input_first, int i_input_last);
	void LAS14_to_LAS13(char** input_list, int i_input_first, int i_input_last);
	void LAS_to_LAZ(char** input_list, int i_input_first, int i_input_last);
	void LAZ_to_LAS(char** input_list, int i_input_first, int i_input_last);
	void Merge(char** input_list, int i_input_first, int i_input_last);

	void Extra_Bytes_Writing_Example(char** input_list, int i_input_first, int i_input_last);
	void Extra_Bytes_Reading_Example(char** input_list, int i_input_first, int i_input_last);

	void Shift_Classification(char** input_list, int i_input_first, int i_input_last);
	void Remove_ExtraBytes(char** input_list, int i_input_first, int i_input_last);

	void Remove_Offset(char** input_list, int i_input_first, int i_input_last);
	void Show_Header(char** input_list, int i_input_first, int i_input_last);
}