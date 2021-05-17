#pragma once
#ifdef _MSC_VER 
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS 
#endif // !_CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <algorithm>
#include <math.h>
#include "EZPC_util.h"
#include "E57Foundation.h"
#include "E57Simple.h"

#define E57_IO_MAX_STRING_SIZE 256
#define E57_IO_BLOCK_SIZE 10000000 
#define E57_IO_COORDINATE_METADATA "EZE57 by Ezra"

struct E57_Header
{
	// point size
	__int64 point_count;

	// data info
	char name[E57_IO_MAX_STRING_SIZE];
	char guid[E57_IO_MAX_STRING_SIZE];
	char original_guid[E57_IO_MAX_STRING_SIZE][E57_IO_MAX_STRING_SIZE];
	int original_guid_count = 0;
	char description[E57_IO_MAX_STRING_SIZE];

	// sensor info
	char vendor[E57_IO_MAX_STRING_SIZE];
	char model[E57_IO_MAX_STRING_SIZE];
	char serial_number[E57_IO_MAX_STRING_SIZE];
	char hardware_version[E57_IO_MAX_STRING_SIZE];
	char software_version[E57_IO_MAX_STRING_SIZE];
	char firmware_version[E57_IO_MAX_STRING_SIZE];

	// environment
	float temperature;
	float relative_humidity;
	float atmospheric_pressure;

	// data acquisition time
	double date_time_start;
	double date_time_end;
	bool data_time_start_isAtomicClockReferenced;
	bool data_time_end_isAtomicClockReferenced;

	// pose
	double translation_x;
	double translation_y;
	double translation_z;

	double quaternion_x;
	double quaternion_y;
	double quaternion_z;
	double quaternion_w;

	// xyz
	bool x_field;
	bool y_field;
	bool z_field;
	bool xyz_invalid_field;
	double x_minimum;
	double x_maximum;
	double y_minimum;
	double y_maximum;
	double z_minimum;
	double z_maximum;

	// spherical
	bool range_field;
	bool azimuth_field;
	bool elevation_field;
	bool spherical_invalid_field;
	double range_minimum;
	double range_maximum;
	double azimuth_start;
	double azimuth_end;
	double elevation_minimum;
	double elevation_maximum;

	// index
	bool row_index_field;
	bool column_index_field;
	bool return_index_field;
	bool return_count_field;
	__int64 row_maximum_limit;
	__int64 column_maximum_limit;
	__int64 return_maximum_limit;
	__int64 row_minimum;
	__int64 row_maximum;
	__int64 column_minimum;
	__int64 column_maximum;
	__int64 return_minimum;
	__int64 return_maximum;

	// intensity
	bool intensity_field;
	bool intensity_invalid_field;
	double intensity_scale_factor;  // -1 = integer; 0 = float; others = scale factor.	
	double intensity_minimum_limit;
	double intensity_maximum_limit;

	// color 
	bool red_field;
	bool green_field;
	bool blue_field;
	bool rgb_invalid_field;
	double red_minimum_limit;
	double red_maximum_limit;
	double green_minimum_limit;
	double green_maximum_limit;
	double blue_minimum_limit;
	double blue_maximum_limit;
	
	// scale
	double distance_minimum_limit;			
	double distance_maximum_limit;			
	double distance_scale_factor; // -1 = integer; 0 = float; others = scale factor.	

	double angle_minimum_limit;
	double angle_maximum_limit;
	double angle_scale_factor; // -1 = integer; 0 = float; others = scale factor.	

	// time
	bool time_stamp_field;
	bool time_stamp_invalid_field;
	double time_stamp_maximum_limit;

};

struct E57_Block
{
	double* x = NULL;
	double* y = NULL;
	double* z = NULL;
	int8_t* invalid_xyz = NULL;

	double* intensity = NULL;
	int8_t* invalid_intensity = NULL;

	uint16_t* red = NULL;
	uint16_t* green = NULL;
	uint16_t* blue = NULL;
	int8_t* invalid_rgb = NULL;

	double* range = NULL;
	double* azimuth = NULL;
	double* elevation = NULL;
	int8_t* invalid_spherical = NULL;

	double* time_stamp = NULL;
	int8_t* invalid_time_stamp = NULL;

	int32_t* row_index = NULL;
	int32_t* column_index = NULL;

	int8_t* return_index = NULL;
	int8_t* return_count = NULL;
};

struct E57_Point
{
	double x;
	double y;
	double z;
	
	double intensity;

	unsigned short red;
	unsigned short green;
	unsigned short blue;

	double range;
	double azimuth;
	double elevation;

	double time_stamp;

	__int64 row_index;
	__int64 column_index;

	char return_index;
	char return_count;
	bool invalid_xyz;
	bool invalid_intensity;
	bool invalid_rgb;
	bool invalid_spherical;
	bool invalid_time_stamp;
};

class E57_Scan
{
public:
	char source_file[E57_IO_MAX_STRING_SIZE];
	__int64 scan_index = -1;
	E57_Header e57_header;
	E57_Point* e57_point = NULL;

	void Read_E57_Scan(char* filename, __int64 scn_idx);
	void Read_E57_Scan();

	void Read_E57_Header(char* filename, __int64 scn_idx);
	void Read_E57_Header();
	void Read_E57_Point();

	void Copy_To_E57_Scan(E57_Scan& dst_scn);
	void Copy_From_E57_Scan(E57_Scan src_scn);
	void Auto_Fill_E57_Point_Field();
	void Set_XYZ_From_Spherical();

	E57_Point global_e57_point(E57_Point local_point);
	E57_Point local_e57_point(E57_Point global_point);
	E57_Point e57_scan_point(__int64 point_index, bool local_or_global);

	void Write_E57_Scan(char* filename);
	void Write_E57_Scan(e57::Writer eWriter); // this will keep the eWriter open.

	e57::Data3D Data3D_from_E57_Header();
	void E57_Header_from_Data3D(e57::Data3D header);

	E57_Block create_e57_block(__int64 block_size);
	void Delete_E57_Block(E57_Block &e57_block);

	E57_Point e57_point_from_block(E57_Block block, __int64  block_idx);
	void Block_from_E57_Point(E57_Block &block, __int64 block_idx, E57_Point point);

	void Initialize();
	void Initialize(E57_Header header);
	void Initialize(__int64 n_pts);

	void Delete();

	void Show_E57_Scan_Info();
	void Show_E57_Point(__int64 i_p);
};

class E57_io
{
public:
	__int64 scan_count = 0; 
	E57_Scan* e57_scan = NULL;

	void Add_E57_Scan(char** file_list, int i_file_start, int i_file_end);
	void Add_E57_Scan(char* file_name);
	void Add_E57_Scan(char* file_name, __int64 scan_index);
	void Add_E57_Scan(E57_Scan add_scans);
	void Add_E57_Scan(E57_Scan* add_scans, __int64 n_add_scn);
	
	void Initialize();
	void Initialize(int n_scn);

	E57_Header e57_scan_header(__int64 scan_index);
	E57_Point e57_scan_point(__int64 scan_index, __int64 point_index, bool local_or_global);

	void Write_E57_Scan(char* file_name, bool merge_or_separate);

	void Show_E57_Scan_Info(__int64 i_scan);
	void Delete();	
};

__int64 e57_scan_count(char* file_name);
E57_Header* Get_E57_Header(char** file_list, int i_file_start, int i_file_end, __int64 &header_count);
E57_Header* Get_E57_Header(char* file_name, __int64 &header_count);
void Print_E57_Header(E57_Header header);