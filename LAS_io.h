#pragma once
/*
LAS_io.h   - Ezra Che

06/18/2020
Support Extra Bytes Reading and Writing. 


11/08/2019
it does not support recognizing the waveform data packet.
it does not support coordinate reference system representation.
need to work on pre-defined VLR

*/
#ifdef _MSC_VER 
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS 
#endif // !_CRT_SECURE_NO_WARNINGS
#endif

#if defined(__linux__)
#define _fseeki64 fseeko64
#endif

#include "laszip_api.h"
#include <stdlib.h>
#include <cfloat>
#include <climits>
#include <cmath>
#include "EZPC_util.h"

#define LASIO_POINT_BLOCK_SIZE_R 100000  // for reading
#define LASIO_POINT_BLOCK_SIZE_W 500000  // for writing

#define LASIO_OFFSET_TO_VERSION_MAJOR 24

#define LASIO_VLR_Header_Size 54
#define LASIO_EVLR_Header_Size 54

#define LASIO_VLR_Extra_Bytes_Size 192

#define INVALID_FILE_FLAG 0
#define LAS_FILE_FLAG 1
#define LAZ_FILE_FLAG 2

#define LASIO_DEFAULT_VERSION_MODE true
#define LASIO_DEFAULT_FORMAT_MODE true
#define LASIO_DEFAULT_SCALE 0.000001

#define LASIO_N_POINT_DATA_FORMAT 11
#define LASIO_N_PDF_OPT_FIELD 6

typedef unsigned char LAS_BYTE;
typedef unsigned long long LAS_BYTE_8;

const unsigned short LAS_Header_Size_1_X[5]
{
	227,   // 1.0
	227,   // 1.1
	227,   // 1.2
	235,   // 1.3
	375	   // 1.4
};

const bool LAS_Optional_Field[LASIO_N_POINT_DATA_FORMAT][LASIO_N_PDF_OPT_FIELD]
{
	// t-2/f-1	ang(2/1) time	color	nir		wave   
		{false,	false,	false,	false,	false,	false	},	// 00	
		{false,	false,	true,	false,	false,	false	},	// 01
		{false,	false,	false,	true,	false,	false	},	// 02
		{false,	false,	true,	true,	false,	false	},	// 03
		{false,	false,	true,	false,	false,	true	},	// 04
		{false,	false,	true,	true,	false,	false	},	// 05
		{true,	true ,	true,	false,	false,	false	},	// 06
		{true,	true ,	true,	true,	false,	false	},	// 07
		{true,	true ,	true,	true,	true,	false	},	// 08
		{true,	true ,	true,	true,	false,	true	},	// 09
		{true,	true ,	true,	true,	true,	true	},	// 10
};

enum LASIO_PDF_OPT_FIELD : unsigned char
{
	LASIO_PDF_OPT_Extended_Flight_Info = 0,
	LASIO_PDF_OPT_Extended_Scan_Angle = 1,
	LASIO_PDF_OPT_GPS_Time = 2,
	LASIO_PDF_OPT_Color = 3,
	LASIO_PDF_OPT_NIR = 4,
	LASIO_PDF_OPT_Waveform = 5
};

const unsigned short LAS_Point_Size[LASIO_N_POINT_DATA_FORMAT]
{
	20, // 00
	28, // 01
	26, // 02
	34, // 03
	57, // 04
	63, // 05
	30, // 06
	36, // 07
	38, // 08
	59, // 09
	67  // 10
};

enum PreDefined_VLR_Record_ID
{
	VLR_Record_ID_Classification_Lookup = 0,
	VLR_Record_ID_Text_Area_Description = 3,
	VLR_Record_ID_Extra_Bytes = 4,
	VLR_Record_ID_Superseded = 7,
	//	Waveform_Packet_Descriptor = 99 // 99 - 255 waveform packet descriptor

	VLR_Record_ID_OGC_Math_Transform_WKT = 2111,
	VLR_Record_ID_OGC_Cooridnate_System_WKT = 2112,
	VLR_Record_ID_GeoKeyDirectoryTag = 34735,
	VLR_Record_ID_GeoDoubleParamsTag = 34736,
	VLR_Record_ID_GeoAsciiParamsTag = 34737
};

struct LAS_Header
{
	char File_Signature[4] = { 'L', 'A', 'S', 'F' };
	unsigned short File_Source_ID = 0;
	unsigned short Global_Encoding = 0;
	unsigned int Project_ID_GUID_data_1 = 0;
	unsigned short Project_ID_GUID_data_2 = 0;
	unsigned short Project_ID_GUID_data_3 = 0;
	unsigned char Project_ID_GUID_data_4[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
	unsigned char Version_Major = 1;
	unsigned char Version_Minor = 0;
	char System_Identifier[32];
	char Generating_Software[32] = { "EZLAS by Ezra Che" };
	unsigned short File_Creation_Day_of_Year = 0;
	unsigned short File_Creation_Year = 0;
	unsigned short Header_Size = LAS_Header_Size_1_X[0];
	unsigned int Offset_to_point_data = LAS_Header_Size_1_X[0];
	unsigned int Number_of_Variable_Length_Records = 0;
	unsigned char Point_Data_Format_ID = 0;
	unsigned short Point_Data_Record_Length = LAS_Point_Size[0];

	unsigned int Legacy_Number_of_point_records = 0;
	unsigned int Legacy_Number_of_points_by_return[5] = { 0, 0, 0, 0, 0 };

	double X_scale_factor = 0.;
	double Y_scale_factor = 0.;
	double Z_scale_factor = 0.;
	double X_offset = 0.;
	double Y_offset = 0.;
	double Z_offset = 0.;
	double Max_X = -DBL_MAX;
	double Min_X = DBL_MAX;
	double Max_Y = -DBL_MAX;
	double Min_Y = DBL_MAX;
	double Max_Z = -DBL_MAX;
	double Min_Z = DBL_MAX;
	unsigned long long Start_of_Waveform_Data_Packet_Record = 0;
	unsigned long long Start_of_First_Extended_Variable_Length_Record = 0;
	unsigned int Number_of_Extended_Variable_Length_Records = 0;

	unsigned long long Number_of_point_records = 0;
	unsigned long long Number_of_points_by_return[15] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 , 0, 0, 0, 0, 0 };
};

typedef struct LAS_VLR_Header
{
	unsigned short Reserved = 0;
	char User_ID[16] = "\0";
	unsigned short Record_ID = 0;
	unsigned short Record_Length_After_Header = 0;
	char Description[32] = "\0";
} LAS_EVLR_Header;

typedef struct LAS_VLR
{
	LAS_VLR_Header vlr_header;
	LAS_BYTE* vlr_buffer = NULL;
} LAS_EVLR;

enum LAS_Classification : unsigned char
{
	LAS_Class_Never_Classified = 0,
	LAS_Class_Unclassified = 1,
	LAS_Class_Ground = 2,
	LAS_Class_Low_Vegetation = 3,
	LAS_Class_Medium_Vegetation = 4,
	LAS_Class_High_Vegetation = 5,
	LAS_Class_Building = 6,
	LAS_Class_Low_Point = 7,
	LAS_Class_Reserved_1 = 8,
	LAS_Class_Water = 9,
	LAS_Class_Rail = 10,
	LAS_Class_Road_Surface = 11,
	LAS_Class_Reserved_2 = 12,
	LAS_Class_Wire_Guard = 13,
	LAS_Class_Wire_Conductor = 14,
	LAS_Class_Transmission_Tower = 15,
	LAS_Class_Wire_Structure_Connector = 16,
	LAS_Class_Bridge_Deck = 17,
	LAS_Class_High_Noise = 18,
	LAS_Class_Overhead_Structure = 19,
	LAS_Class_Ignored_Ground = 20,
	LAS_Class_Snow = 21,
	LAS_Class_Temporal_Exclusion = 22,
};

struct LAS_Point
{
	int X = 0;
	int Y = 0;
	int Z = 0;
	unsigned short Intensity = 0;

	unsigned char Return_Number = 0;
	unsigned char Number_of_Returns = 0;
	unsigned char Scan_Direction_Flag = 0;
	unsigned char Edge_of_Flight_Line = 0;
	unsigned char Classification_Flags = 0;
	unsigned char Scanner_Channel = 0;

	unsigned char Classification = LAS_Class_Never_Classified;
	short Scan_Angle_Rank = 0;

	unsigned char User_Data = 0;
	unsigned short Point_Source_ID = 0;
	double GPS_Time = 0;
	unsigned short Red = 0;
	unsigned short Green = 0;
	unsigned short Blue = 0;

	unsigned short NIR = 0;

	unsigned char Wave_Packet_Discriptor_Index = 0;
	unsigned long long Byte_offset_to_waveform_data = 0;
	unsigned int Waveform_packet_size_in_bytes = 0;
	float Return_Point_Waveform_Location = 0;
	float Xt = 0;
	float Yt = 0;
	float Zt = 0;

	LAS_BYTE *Extra_Bytes = NULL;
};

struct sGeoKeys
{
	unsigned short wKeyDirectoryVersion = 1;
	unsigned short wKeyRevision = 1;
	unsigned short wMinorRevision = 0;
	unsigned short wNumberOfKeys;
	struct sKeyEntry
	{
		unsigned short wKeyID;
		unsigned short wTIFFTagLocation;
		unsigned short wCount;
		unsigned short wValue_Offset;
	} pKey[1];
};

struct Waveform_Packet_Descriptor
{
	unsigned char Bits_per_Sample;
	unsigned char Waveform_Compression_Type = 0;
	unsigned int Number_of_Samples;
	unsigned int Temporal_Sample_Spacing;
	double Digitizer_Gain;
	double Digitizer_Offset;
};

enum LAS_Extra_Bytes_Data_Type : unsigned char
{
	EB_Type_Undocumented = 0,
	EB_Type_Unsigned_Char = 1,
	EB_Type_Char = 2,
	EB_Type_Unsigned_Short = 3,
	EB_Type_Short = 4,
	EB_Type_Unsigned_Long = 5,
	EB_Type_Long = 6,
	EB_Type_Unsigned_Long_Long = 7,
	EB_Type_Long_Long = 8,
	EB_Type_Float = 9,
	EB_Type_Double = 10,
	EB_Type_Unsigned_Char_2 = 11,
	EB_Type_Char_2 = 12,
	EB_Type_Unsigned_Short_2 = 13,
	EB_Type_Short_2 = 14,
	EB_Type_Unsigned_Long_2 = 15,
	EB_Type_Long_2 = 16,
	EB_Type_Unsigned_Long_Long_2 = 17,
	EB_Type_Long_Long_2 = 18,
	EB_Type_Float_2 = 19,
	EB_Type_Double_2 = 20,
	EB_Type_Unsigned_Char_3 = 21,
	EB_Type_Char_3 = 22,
	EB_Type_Unsigned_Short_3 = 23,
	EB_Type_Short_3 = 24,
	EB_Type_Unsigned_Long_3 = 25,
	EB_Type_Long_3 = 26,
	EB_Type_Unsigned_Long_Long_3 = 27,
	EB_Type_Long_Long_3 = 28,
	EB_Type_Float_3 = 29,
	EB_Type_Double_3 = 30,
};

const unsigned short LAS_Extra_Bytes_Type_Info[][4]
{
	//each num floating singed 
	{ 0, 1, 1, 1 },	//EB_Undocumented
	{ 1, 1, 0, 0 },	//EB_Unsigned_Char
	{ 1, 1, 0, 1 },	//EB_Char
	{ 2, 1, 0, 0 },	//EB_Unsigned_Short
	{ 2, 1, 0, 1 },	//EB_Short
	{ 4, 1, 0, 0 },	//EB_Unsigned_Long
	{ 4, 1, 0, 1 },	//EB_Long
	{ 8, 1, 0, 0 },	//EB_Unsigned_Long_Long
	{ 8, 1, 0, 1 },	//EB_Long_Long
	{ 4, 1, 1, 1 },	//EB_Float
	{ 8, 1, 1, 1 },	//EB_Double
	{ 1, 2, 0, 0 },	//EB_Unsigned_Char_2
	{ 1, 2, 0, 1 },	//EB_Char_2
	{ 2, 2, 0, 0 },	//EB_Unsigned_Short_2
	{ 2, 2, 0, 1 },	//EB_Short_2
	{ 4, 2, 0, 0 },	//EB_Unsigned_Long_2
	{ 4, 2, 0, 1 },	//EB_Long_2
	{ 8, 2, 0, 0 },	//EB_Unsigned_Long_Long_2
	{ 8, 2, 0, 1 },	//EB_Long_Long_2
	{ 4, 2, 1, 1 },	//EB_Float_2
	{ 8, 2, 1, 1 },	//EB_Double_2
	{ 1, 3, 0, 0 },	//EB_Unsigned_Char_3
	{ 1, 3, 0, 1 },	//EB_Char_3
	{ 2, 3, 0, 0 },	//EB_Unsigned_Short_3
	{ 2, 3, 0, 1 },	//EB_Short_3
	{ 4, 3, 0, 0 },	//EB_Unsigned_Long_3
	{ 4, 3, 0, 1 },	//EB_Long_3
	{ 8, 3, 0, 0 },	//EB_Unsigned_Long_Long_3
	{ 8, 3, 0, 1 },	//EB_Long_Long_3
	{ 4, 3, 1, 1 },	//EB_Float_3
	{ 8, 3, 1, 1 },	//EB_Double_3
};

struct LAS_Extra_Bytes_Info
{
	unsigned char reserved[2] = { 0, 0 };
	unsigned char data_type = EB_Type_Undocumented;

	bool options_no_data_bit = false;
	bool options_min_bit = false;
	bool options_max_bit = false;
	bool options_scale_bit = false;
	bool options_offset_bit = false;

	char name[32] = "\0";
	unsigned char unused[4] = { 0, 0, 0, 0 };
	LAS_BYTE_8 no_data[3] = { 0, 0, 0 };
	LAS_BYTE_8 min[3] = { 0, 0, 0 };
	LAS_BYTE_8 max[3] = { 0, 0, 0 };
	double scale[3] = { 0, 0, 0 };
	double offset[3] = { 0, 0, 0 };
	char description[32] = "\0";

	bool floating_or_integer = true;
	bool signed_or_unsigned = true;

	unsigned short size_value = 0;
	unsigned short n_values = 1;
	unsigned short byte_pos = 0;
};

class LAS_io
{
public:
	int n_thread;
	LAS_Header las_header;
	LAS_VLR *las_vlr = NULL;
	LAS_Point *las_points = NULL;
	LAS_EVLR *las_evlr = NULL;

	void Initialize(LAS_Header header);
	void Initialize(LAS_Header header, bool point, bool vlr, bool evlr);
	void Initialize(bool point, bool vlr, bool evlr);
	void Initialize();
	void Initialize(double X_scale_factor,
		double Y_scale_factor,
		double Z_scale_factor,
		double X_offset,
		double Y_offset,
		double Z_offset,
		unsigned long long Number_of_point_records,
		unsigned int Number_of_Variable_Length_Records,
		unsigned int Number_of_Extended_Variable_Length_Records);

	void Autofill_Header(unsigned char Version_Major,
		unsigned char Version_Minor,
		unsigned char Point_Data_Format_ID,
		unsigned long long Number_of_point_records,
		unsigned int Number_of_Variable_Length_Records,
		unsigned int Number_of_Extended_Variable_Length_Records);

	void Autofill_Header();

	LAS_Header* Read_LAS_Headers(char** file_list, int i_file_start, int i_file_end, bool get_extent);
	LAS_Header Read_LAS_Header(char* file_name, bool get_extent);

	LAS_Header Merge_LAS_Headers(LAS_Header* headers, int n_headers, bool version_mode, bool format_mode); 
		// meta info takes from latest file created, version mode: merge to true(the latest version)/false(the oldest verision)  format mode: true(most fields) false(least fields)  This will ignore extrabytes and vlr
	LAS_Header Merge_LAS_Headers(char** file_list, int i_file_start, int i_file_end, LAS_VLR*& vlr, LAS_EVLR*& evlr);
	LAS_Header Merge_LAS_Headers(char** file_list, int i_file_start, int i_file_end, LAS_VLR*& vlr);
	LAS_Header Merge_LAS_Headers(char** file_list, int i_file_start, int i_file_end);

	LAS_VLR* Read_LAS_VLRs(char* file_name, unsigned int&n_vlrs);
	LAS_VLR* Read_LAS_EVLRs(char* file_name, unsigned int&n_evlrs);
	void Merge_LAS_VLRs(LAS_VLR*& vlr_des, unsigned int& n_vlr_des, LAS_VLR* vlr2, unsigned int n_vlr2);
	void Merge_LAS_EVLRs(LAS_VLR*& evlr_des, unsigned int& n_evlr_des, LAS_VLR* evlr2, unsigned int n_evlr2);

	void Copy_from_LAS_io(LAS_io las_src, bool points, bool vlrs, bool evlrs);
	void Copy_from_LAS_io(LAS_io las_src);

	void Read_LAS_Files(char** file_list, int i_file_start, int i_file_end, bool version_mode, bool format_mode, bool load_points, bool get_extent, bool load_vlrs, bool load_evlrs);
	void Read_LAS_Files(char** file_list, int i_file_start, int i_file_end);
	void Read_LAS_File(char* file_name, bool load_points, bool get_extent, bool load_vlrs, bool load_evlrs);
	void Read_LAS_File(char* file_name);

	void Sample_LAS_Points(bool* flag_sample);
	void Sample_LAS_File(char* file_name, unsigned long long n_sample);
	void Sample_LAS_Files(char** file_list, int i_file_start, int i_file_end, unsigned long long n_sample);

	void Add_LAS_Points(char* file_name, unsigned long long in_i_start, unsigned long long in_i_end, unsigned long long out_i_satrt, bool load_points, bool get_extent); // add points from file.
	void Add_LAS_Points(char* file_name, unsigned long long out_i_start, bool load_points, bool get_extent); // add points from file.
	void Add_LAS_VLRs(char* file_name, unsigned int in_i_start, unsigned int in_i_end, unsigned int out_i_start); // add vlrs from file.
	void Add_LAS_VLRs(char* file_name, unsigned int out_i_start); // add vlrs from file.
	void Add_LAS_EVLRs(char* file_name, unsigned int in_i_start, unsigned int in_i_end, unsigned int out_i_start); // add evlrs from file.
	void Add_LAS_EVLRs(char* file_name, unsigned int out_i_start); // add evlrs from file.

	LAS_Header Load_LAS_Header(FILE* fp_in);
	LAS_Header Load_LAS_Header_1_X(FILE* fp_in);

	LAS_Point Load_A_Point(FILE* fp_in, unsigned short Point_Format_ID);
	LAS_VLR Load_A_VLR(FILE* fp_in, unsigned int &offset);
	LAS_EVLR Load_A_EVLR(FILE* fp_in, unsigned long long &offset);

	void Sort_by_Time();
	void Sort_by_Source_ID();
	void Sort_by_Classification();

	void Write_LAS(char* file_name);
	void Write_LAS_Header(FILE* fp_out);
	void Write_LAS_Header_1_X(FILE* fp_out);

	void Write_LAS_Points(FILE* fp_out);
	void Write_LAS_VLRs(FILE* fp_out);
	void Write_LAS_EVLRs(FILE* fp_out);

	void Point_to_Block(LAS_BYTE *point_bytes, unsigned long long &offset, LAS_Point p); // write to bytes
	void Write_A_Point(FILE* fp_out, LAS_Point p); // write to file
	void Write_A_VLR(FILE* fp_out, unsigned int &offset, LAS_VLR vlr);
	void Write_A_EVLR(FILE* fp_out, unsigned long long &offset, LAS_EVLR evlr);

	//vlr
	LAS_VLR create_a_vlr(char* user_id, unsigned short record_id, unsigned short record_length_after_header);
	LAS_VLR create_a_vlr(unsigned short record_length_after_header);
	void Set_VLR_User_ID(LAS_VLR &vlr, char* user_id);
	void Set_VLR_Record_ID(LAS_VLR &vlr, unsigned short record_id);
	void Set_VLR_Description(LAS_VLR &vlr, char* description);
	void Set_VLR_Bytes(LAS_VLR &vlr, LAS_BYTE* bytes, unsigned short size_bytes);

	void Add_New_VLRs(LAS_VLR* vlr, unsigned int n_vlr);
	void Add_New_VLR(LAS_VLR vlr);

	bool is_vlr_matched(LAS_VLR vlr_1, LAS_VLR vlr_2, bool header_only);

	LAS_VLR OGC_Math_Transform_WKT_Record_VLR();
	LAS_VLR OGC_Coordinate_System_WKT_Record_VLR();
	LAS_VLR GeoKeyDirectoryTag_Record_VLR();
	LAS_VLR GeoDoubleParamsTag_Record_VLR();
	LAS_VLR GeoAsciiParamsTag_Record_VLR();
	LAS_VLR Classification_Lookup_VLR();
	LAS_VLR Text_Area_Description_VLR();
	LAS_VLR Extra_Bytes_VLR();
	LAS_VLR Superseded_VLR();
	LAS_VLR Waveform_Packet_Descriptor_VLR(unsigned short idx);
	LAS_EVLR Waveform_Data_Packets_EVLR();
	
	//extra bytes
	unsigned short size_extra_bytes(bool header_or_vlr); // how many bytes that are extra
	unsigned short size_extra_bytes(unsigned short Point_Format_ID, unsigned short Point_Record_Length);  // how many bytes that are extra	
	unsigned short size_extra_bytes(LAS_Extra_Bytes_Info *eb_info, unsigned short n_eb);  // how many bytes that are extra

	long long idx_vlr_extra_bytes(); // vlr index for extra bytes
	unsigned short n_extra_bytes(); // number of extra bytes/attributes
	long long Get_All_Extra_Bytes_Info(LAS_Extra_Bytes_Info* &eb_info, unsigned short &n_eb); // get all the extra bytes info
	long long Get_One_Extra_Bytes_Info(LAS_Extra_Bytes_Info &eb_info, unsigned short idx_eb); // get an extra bytes info
	int Get_Extra_Bytes_Info(LAS_Extra_Bytes_Info& eb_info, char* name, LAS_Extra_Bytes_Data_Type type_id);

	void Bytes_to_Extra_Bytes_Info(LAS_BYTE* bytes, LAS_Extra_Bytes_Info &eb_info);
	void Bytes_from_Extra_Bytes_Info(LAS_BYTE* bytes, LAS_Extra_Bytes_Info eb_info);

	unsigned short Add_Extra_Bytes_VLR(LAS_Extra_Bytes_Info *extra_bytes_info, unsigned short n_extra_bytes); // add fields to vlr and allocate extra bytes for points, return starting idx of extra bytes
	
	void Initialize_Point_Extra_Bytes(); // allocate extra bytes from header
	void Clear_Extra_Bytes(); // clear all extra bytes in both vlrs and points
	void Remove_Extra_Bytes(char* name, LAS_Extra_Bytes_Data_Type type_id);


	LAS_Extra_Bytes_Info create_extra_bytes_info(char* name, char* description, LAS_Extra_Bytes_Data_Type data_type); // create an extra bytes info structure
	LAS_Extra_Bytes_Info create_extra_bytes_info(LAS_Extra_Bytes_Data_Type data_type); // create an extra bytes info structure
	void Set_Extra_Bytes_Type_Info(LAS_Extra_Bytes_Info &eb_info); // this part is auto
	void Set_Extra_Bytes_Name(LAS_Extra_Bytes_Info &eb_info, char* name);
	void Set_Extra_Bytes_Description(LAS_Extra_Bytes_Info &eb_info, char* description);

	void Set_Extra_Bytes_Options_Invalid(LAS_Extra_Bytes_Info &eb_info);
	void Set_Extra_Bytes_No_Data_Invalid(LAS_Extra_Bytes_Info &eb_info);
	void Set_Extra_Bytes_Min_Invalid(LAS_Extra_Bytes_Info &eb_info);
	void Set_Extra_Bytes_Max_Invalid(LAS_Extra_Bytes_Info &eb_info);
	void Set_Extra_Bytes_Scale_Invalid(LAS_Extra_Bytes_Info &eb_info);
	void Set_Extra_Bytes_Offset_Invalid(LAS_Extra_Bytes_Info &eb_info);

	void Set_Extra_Bytes_No_Data(LAS_Extra_Bytes_Info &eb_info, double no_data);
	void Set_Extra_Bytes_No_Data(LAS_Extra_Bytes_Info &eb_info, long long no_data);
	void Set_Extra_Bytes_No_Data(LAS_Extra_Bytes_Info &eb_info, unsigned long long no_data);

	void Set_Extra_Bytes_Min(LAS_Extra_Bytes_Info &eb_info, double no_data);
	void Set_Extra_Bytes_Min(LAS_Extra_Bytes_Info &eb_info, long long no_data);
	void Set_Extra_Bytes_Min(LAS_Extra_Bytes_Info &eb_info, unsigned long long no_data);

	void Set_Extra_Bytes_Max(LAS_Extra_Bytes_Info &eb_info, double no_data);
	void Set_Extra_Bytes_Max(LAS_Extra_Bytes_Info &eb_info, long long no_data);
	void Set_Extra_Bytes_Max(LAS_Extra_Bytes_Info &eb_info, unsigned long long no_data);

	void Set_Extra_Bytes_No_Data(LAS_Extra_Bytes_Info &eb_info, void* no_data_x, void* no_data_y, void* no_data_z);
	void Set_Extra_Bytes_No_Data(LAS_Extra_Bytes_Info &eb_info, void* no_data_x, void* no_data_y);
	void Set_Extra_Bytes_No_Data(LAS_Extra_Bytes_Info &eb_info, void* no_data);

	void Set_Extra_Bytes_Min(LAS_Extra_Bytes_Info &eb_info, void* min_x, void* min_y, void* min_z);
	void Set_Extra_Bytes_Min(LAS_Extra_Bytes_Info &eb_info, void* min_x, void* min_y);
	void Set_Extra_Bytes_Min(LAS_Extra_Bytes_Info &eb_info, void* min);

	void Set_Extra_Bytes_Max(LAS_Extra_Bytes_Info &eb_info, void* max);
	void Set_Extra_Bytes_Max(LAS_Extra_Bytes_Info &eb_info, void* max_x, void* max_y);
	void Set_Extra_Bytes_Max(LAS_Extra_Bytes_Info &eb_info, void* max_x, void* max_y, void* max_z);

	void Set_Extra_Bytes_Scale(LAS_Extra_Bytes_Info &eb_info, double scale_x, double scale_y, double scale_z);
	void Set_Extra_Bytes_Scale(LAS_Extra_Bytes_Info &eb_info, double scale_x, double scale_y);
	void Set_Extra_Bytes_Scale(LAS_Extra_Bytes_Info &eb_info, double scale);

	void Set_Extra_Bytes_Offset(LAS_Extra_Bytes_Info &eb_info, double offset_x, double offset_y, double offset_z);
	void Set_Extra_Bytes_Offset(LAS_Extra_Bytes_Info &eb_info, double offset_x, double offset_y);
	void Set_Extra_Bytes_Offset(LAS_Extra_Bytes_Info &eb_info, double offset);

	void Set_Extra_Bytes_Value_No_Data(unsigned long long pid, LAS_Extra_Bytes_Info eb_info);

	void Set_Extra_Bytes_Value(unsigned long long pid, LAS_Extra_Bytes_Info eb_info, double in_value); // add a field value
	void Set_Extra_Bytes_Value(unsigned long long pid, LAS_Extra_Bytes_Info eb_info, long long in_value); // add a field value
	void Set_Extra_Bytes_Value(unsigned long long pid, LAS_Extra_Bytes_Info eb_info, unsigned long long in_value); // add a field value

	void Set_Extra_Bytes_Value(unsigned long long pid, LAS_Extra_Bytes_Info eb_info, void* in_value); // add a field value
	void Set_Extra_Bytes_Value(unsigned long long pid, LAS_Extra_Bytes_Info eb_info, void* in_value_x, void* in_value_y); // add a field value
	void Set_Extra_Bytes_Value(unsigned long long pid, LAS_Extra_Bytes_Info eb_info, void* in_value_x, void* in_value_y, void* in_value_z); // add a field value

	bool is_extra_bytes_value_valid(unsigned long long pid, LAS_Extra_Bytes_Info eb_info);

	void Get_Extra_Bytes_Value(unsigned long long pid, LAS_Extra_Bytes_Info eb_info, double &val);
	void Get_Extra_Bytes_Value(unsigned long long pid, LAS_Extra_Bytes_Info eb_info, long long &val);
	void Get_Extra_Bytes_Value(unsigned long long pid, LAS_Extra_Bytes_Info eb_info, unsigned long long &val);

	void Get_Extra_Bytes_Value_Byte_8(unsigned long long pid, LAS_Extra_Bytes_Info eb_info, LAS_BYTE_8 &val);
	void Get_Extra_Bytes_Value_Byte_8(unsigned long long pid, LAS_Extra_Bytes_Info eb_info, LAS_BYTE_8 &val_x, LAS_BYTE_8 &val_y);
	void Get_Extra_Bytes_Value_Byte_8(unsigned long long pid, LAS_Extra_Bytes_Info eb_info, LAS_BYTE_8 &val_x, LAS_BYTE_8 &val_y, LAS_BYTE_8 &val_z);

	void Upcasting_to_Byte_8(LAS_Extra_Bytes_Info eb_info, LAS_BYTE_8 &bytes, void* value);
	void Downcasting_to_Bytes(LAS_Extra_Bytes_Info eb_info, LAS_BYTE* bytes, void* value);

	////////////////

	void Show_LAS_Header();
	void Show_Point(unsigned long long i);
	void Show_VLR(unsigned int i);
	void Show_EVLR(unsigned int i);
	void Show_Extra_Bytes_VLR();

	int x_coord2record(double x);
	int y_coord2record(double y);
	int z_coord2record(double z);
	double x_record2coord(int x);
	double y_record2coord(int y);
	double z_record2coord(int z);
	void Delete_Data();

	//laszip

	laszip_POINTER Open_LAZ_Reader(char* filename);
	LAS_Header Load_LAZ_Header(laszip_POINTER laszip_reader);
	LAS_VLR Load_A_VLR(laszip_POINTER laszip_reader, unsigned int i_vlr);
	LAS_Point Load_A_Point(laszip_POINTER laszip_reader, unsigned short Point_Format_ID); //load a laz point
	void Close_LAZ_Reader(laszip_POINTER laszip_reader);
	LAS_Header Convert_Header(laszip_header laz_header); // from laszip to ezlas
	LAS_VLR Convert_VLR(laszip_vlr laz_vlr); // from laszip to ezlas
	LAS_Point Convert_Point(laszip_point *laz_point, unsigned short Point_Format_ID);

	void Write_LAZ_File(char* filename);

	void Set_LAZ_Header(laszip_POINTER laszip_writer);
	void Set_LAZ_VLRs(laszip_POINTER laszip_writer);
	void Set_LAZ_Points(laszip_POINTER laszip_writer);

	laszip_header Convert_Header(LAS_Header las_header); // from ezlas to laszip
	void LAS2LAZ_Point(LAS_Point las_point, laszip_point *laz_point, unsigned short Point_Format_ID);
	void LAS2LAS_Point(LAS_Point las_point_src, LAS_Point &las_point_des, unsigned short size_extra_bytes);
	void LAS2LAS_VLR(LAS_VLR las_vlr_src, LAS_VLR& las_vlr_des);
	void LAS2LAS_EVLR(LAS_EVLR las_evlr_src, LAS_EVLR& las_evlr_des);

	int File_Flag(char* filename); // LAZ

	int coord2record(double offset, double scale_factor, double coord);
	double record2coord(double offset, double scale_factor, int record);
	void Sort_LAS_by_GPS_Time(LAS_Point* pts, unsigned long long npts);
	void Sort_LAS_by_Point_Source_ID(LAS_Point* pts, unsigned long long npts);
	void Sort_LAS_by_Classification(LAS_Point* pts, unsigned long long npts);

	unsigned char Union_Point_Data_Format(unsigned char Format_ID_1, unsigned char Format_ID_2);
	unsigned char Intersect_Point_Data_Format(unsigned char Format_ID_1, unsigned char Format_ID_2);
	unsigned char Match_Point_Data_Format(bool opt_field[LASIO_N_PDF_OPT_FIELD]);
	unsigned char Match_Point_Data_Format(unsigned short point_size);
	unsigned char add_opt_field(unsigned char Format_ID, LASIO_PDF_OPT_FIELD opt_field);
	unsigned char remove_opt_field(unsigned char Format_ID, LASIO_PDF_OPT_FIELD opt_field);
	void Set_Point_Data_Format(LAS_Header &header, unsigned char nu_format_id);

	void Print_LAS_Header_1_X(LAS_Header h);
	void Print_LAS_Point(LAS_Point p, unsigned short Point_Format_ID);
	void Print_LAS_VLR(LAS_VLR vlr);
	void Print_LAS_EVLR(LAS_EVLR evlr);

	bool Is_LAZ_File(char* filename); // if the file is a LAZ.
	bool Is_LAS_File(char* filename); // if the file is a LAS.

	void Data2Bytes(LAS_BYTE* bytes, void* data, unsigned short data_size, unsigned short n_data, unsigned long long &pos);
	void Bytes2Data(LAS_BYTE* bytes, void* data, unsigned short data_size, unsigned short n_data, unsigned long long &pos);

	double double_from_byte_8(LAS_BYTE_8 bytes);
	long long int64_from_byte_8(LAS_BYTE_8 bytes);
	unsigned long long uint64_from_byte_8(LAS_BYTE_8 bytes);
};

bool lasio_compare_gps_time(LAS_Point a, LAS_Point b);
bool lasio_compare_point_source_id(LAS_Point a, LAS_Point b);
bool lasio_compare_classification(LAS_Point a, LAS_Point b);

