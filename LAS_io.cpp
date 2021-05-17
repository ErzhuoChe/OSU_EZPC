#include "LAS_io.h"

void LAS_io::Initialize(LAS_Header header)
{
	Initialize(header, true, true, true);
}

void LAS_io::Initialize(LAS_Header header, bool point, bool vlr, bool evlr)
{
	las_header = header;
	if (las_points && point)
	{
		printf("\n\nThere have been las points initilized!");
		getchar();
	}
	if (las_vlr && vlr)
	{
		printf("\n\nThere have been las vlr initilized!");
		getchar();
	}
	if (las_evlr && evlr)
	{
		printf("\n\nThere have been las evlr initilized!");
		getchar();
	}

	if (header.Number_of_point_records > 0 && point)
	{
		las_points = new LAS_Point[header.Number_of_point_records];
		Initialize_Point_Extra_Bytes();
	}

	if (header.Number_of_Variable_Length_Records > 0 && vlr)
	{
		las_vlr = new LAS_VLR[header.Number_of_Variable_Length_Records];
	}

	if (header.Number_of_Extended_Variable_Length_Records > 0 && evlr)
	{
		las_evlr = new LAS_EVLR[header.Number_of_Extended_Variable_Length_Records];
	}
}

void LAS_io::Initialize(bool point, bool vlr, bool evlr)
{
	Initialize(las_header, point, vlr, evlr);
}

void LAS_io::Initialize()
{
	Initialize(las_header);
}

void LAS_io::Initialize(double X_scale_factor,
	double Y_scale_factor,
	double Z_scale_factor,
	double X_offset,
	double Y_offset,
	double Z_offset,
	unsigned long long Number_of_point_records,
	unsigned int Number_of_Variable_Length_Records,
	unsigned int Number_of_Extended_Variable_Length_Records)
{
	las_header.X_scale_factor = X_scale_factor;
	las_header.Y_scale_factor = Y_scale_factor;
	las_header.Z_scale_factor = Z_scale_factor;
	las_header.X_offset = X_offset;
	las_header.Y_offset = Y_offset;
	las_header.Z_offset = Z_offset;
	las_header.Number_of_point_records = Number_of_point_records;
	las_header.Number_of_Variable_Length_Records = Number_of_Variable_Length_Records;
	las_header.Number_of_Extended_Variable_Length_Records = Number_of_Extended_Variable_Length_Records;

	Initialize();
}

void LAS_io::Autofill_Header(unsigned char Version_Major,
	unsigned char Version_Minor,
	unsigned char Point_Data_Format_ID,
	unsigned long long Number_of_point_records,
	unsigned int Number_of_Variable_Length_Records,
	unsigned int Number_of_Extended_Variable_Length_Records)
{
	time_t rawtime;
	struct tm *timeinfo;

	time(&rawtime);
	timeinfo = gmtime(&rawtime);

	las_header.File_Signature[0] = 'L';
	las_header.File_Signature[1] = 'A';
	las_header.File_Signature[2] = 'S';
	las_header.File_Signature[3] = 'F';

	//las_header.File_Source_ID;
	//las_header.Global_Encoding;

	//las_header.Project_ID_GUID_data_1;
	//las_header.Project_ID_GUID_data_2;
	//las_header.Project_ID_GUID_data_3;
	//las_header.Project_ID_GUID_data_4[8];

	las_header.Version_Major = Version_Major;
	las_header.Version_Minor = Version_Minor;

	//las_header.System_Identifier;
	strncpy(las_header.System_Identifier, "EZLAS Header Autofilling", 32);
	strncpy(las_header.Generating_Software, "EZLAS by Ezra Che", 32);

	las_header.File_Creation_Day_of_Year = timeinfo->tm_yday + 1;
	las_header.File_Creation_Year = timeinfo->tm_year + 1900;

	switch (las_header.Version_Major)
	{
		case 1:
			las_header.Header_Size = LAS_Header_Size_1_X[Version_Minor];
			break;
	}

	size_t tot_vlr_size = 0;	

	if (las_vlr)
	{
		for (unsigned long i = 0; i < Number_of_Variable_Length_Records; i++)
		{
			tot_vlr_size += LASIO_VLR_Header_Size + las_vlr[i].vlr_header.Record_Length_After_Header;
		}
	}

	las_header.Offset_to_point_data = (unsigned long)(las_header.Header_Size + tot_vlr_size);
	las_header.Number_of_Variable_Length_Records = Number_of_Variable_Length_Records;
	las_header.Point_Data_Format_ID = Point_Data_Format_ID;

	las_header.Point_Data_Record_Length = LAS_Point_Size[Point_Data_Format_ID] + size_extra_bytes(false);

	////////////////////////////////// Find Min and Max, count the points by return
	if (las_points)
	{
		int max_x = INT_MIN;
		int max_y = INT_MIN;
		int max_z = INT_MIN;

		int min_x = INT_MAX;
		int min_y = INT_MAX;
		int min_z = INT_MAX;

		for (int i = 0; i < 15; i++)
		{
			las_header.Number_of_points_by_return[i] = 0;
		}

		for (unsigned long long i = 0; i < Number_of_point_records; i++)
		{
			max_x = las_points[i].X > max_x ? las_points[i].X : max_x;
			max_y = las_points[i].Y > max_y ? las_points[i].Y : max_y;
			max_z = las_points[i].Z > max_z ? las_points[i].Z : max_z;

			min_x = las_points[i].X < min_x ? las_points[i].X : min_x;
			min_y = las_points[i].Y < min_y ? las_points[i].Y : min_y;
			min_z = las_points[i].Z < min_z ? las_points[i].Z : min_z;

			las_header.Number_of_points_by_return[las_points[i].Return_Number - 1]++;
		}

		las_header.Max_X = max_x * las_header.X_scale_factor + las_header.X_offset;
		las_header.Min_X = min_x * las_header.X_scale_factor + las_header.X_offset;

		las_header.Max_Y = max_y * las_header.Y_scale_factor + las_header.Y_offset;
		las_header.Min_Y = min_y * las_header.Y_scale_factor + las_header.Y_offset;

		las_header.Max_Z = max_z * las_header.Z_scale_factor + las_header.Z_offset;
		las_header.Min_Z = min_z * las_header.Z_scale_factor + las_header.Z_offset;
	}

	las_header.Number_of_point_records = Number_of_point_records;
	las_header.Legacy_Number_of_point_records = (unsigned long)Number_of_point_records;

	for (int i = 0; i < 5; i++)
	{
		las_header.Legacy_Number_of_points_by_return[i] = (unsigned long)las_header.Number_of_points_by_return[i];
	}

	///////////////////////////////////////////////////////

//	las_header.Start_of_Waveform_Data_Packet_Record = las_header.Offset_to_point_data + las_header.Point_Data_Record_Length * Number_of_point_records;

	las_header.Start_of_First_Extended_Variable_Length_Record = las_header.Offset_to_point_data + las_header.Point_Data_Record_Length * Number_of_point_records;
	las_header.Number_of_Extended_Variable_Length_Records = Number_of_Extended_Variable_Length_Records;
}

void LAS_io::Autofill_Header()
{
	Autofill_Header(las_header.Version_Major,
		las_header.Version_Minor,
		las_header.Point_Data_Format_ID,
		las_header.Number_of_point_records,
		las_header.Number_of_Variable_Length_Records,
		las_header.Number_of_Extended_Variable_Length_Records);
}

LAS_Header* LAS_io::Read_LAS_Headers(char** file_list, int i_file_start, int i_file_end, bool get_extent)
{
	int n_file = i_file_end - i_file_start + 1;
	LAS_Header* headers = new LAS_Header[n_file];

	for (int i_file = i_file_start; i_file <= i_file_end; i_file++)
	{
		headers[i_file - i_file_start] = Read_LAS_Header(file_list[i_file], get_extent);
	}

	return headers;
}

LAS_Header LAS_io::Read_LAS_Header(char* file_name, bool get_extent)
{
	LAS_io las_tmp;

	switch (File_Flag(file_name))
	{
		case LAS_FILE_FLAG:
		{
			FILE *fp = fopen(file_name, "rb");
			las_tmp.las_header = las_tmp.Load_LAS_Header(fp);
			fclose(fp);
			break;
		}
		case LAZ_FILE_FLAG:
		{
			laszip_POINTER laszip_reader = Open_LAZ_Reader(file_name);
			las_tmp.las_header = las_tmp.Load_LAZ_Header(laszip_reader);
			Close_LAZ_Reader(laszip_reader);
			break;
		}
		default:
			break;
	}

	if (get_extent)
	{
		las_tmp.las_header.Min_X = DBL_MAX;
		las_tmp.las_header.Min_Y = DBL_MAX;
		las_tmp.las_header.Min_Z = DBL_MAX;
	
		las_tmp.las_header.Max_X = -DBL_MAX;
		las_tmp.las_header.Max_Y = -DBL_MAX;
		las_tmp.las_header.Max_Z = -DBL_MAX;

		las_tmp.Add_LAS_Points(file_name, 0, false, get_extent);
	}
	return las_tmp.las_header;
}

LAS_Header LAS_io::Merge_LAS_Headers(LAS_Header* headers, int n_headers, bool version_mode, bool format_mode) // version mode: merge to true(the latest version)/false(the oldest verision)  format mode: true(most fields) false(least fields) 
{
	LAS_Header h_tmp = headers[0];
	
	if (n_headers > 1)
	{
		for (int i = 1; i < n_headers; i++)
		{
			if (h_tmp.File_Creation_Year < headers[i].File_Creation_Year ||
				(h_tmp.File_Creation_Year < headers[i].File_Creation_Year && h_tmp.File_Creation_Day_of_Year < headers[i].File_Creation_Day_of_Year))
			{
				h_tmp.File_Creation_Year = headers[i].File_Creation_Year;
				h_tmp.File_Creation_Day_of_Year = headers[i].File_Creation_Day_of_Year;

				memcpy(h_tmp.System_Identifier, headers[i].System_Identifier, 32);
				memcpy(h_tmp.Generating_Software, headers[i].Generating_Software, 32);

				h_tmp.Global_Encoding = h_tmp.Global_Encoding;
				h_tmp.File_Source_ID = headers[i].File_Source_ID;
				h_tmp.Project_ID_GUID_data_1 = headers[i].Project_ID_GUID_data_1;
				h_tmp.Project_ID_GUID_data_2 = headers[i].Project_ID_GUID_data_2;
				h_tmp.Project_ID_GUID_data_3 = headers[i].Project_ID_GUID_data_3;

				for (int j = 0; j < 8; j++)
				{
					h_tmp.Project_ID_GUID_data_4[i] = headers[i].Project_ID_GUID_data_4[j];
				}
			}

			h_tmp.Max_X = h_tmp.Max_X > headers[i].Max_X ? h_tmp.Max_X : headers[i].Max_X;
			h_tmp.Max_Y = h_tmp.Max_Y > headers[i].Max_Y ? h_tmp.Max_Y : headers[i].Max_Y;
			h_tmp.Max_Z = h_tmp.Max_Z > headers[i].Max_Z ? h_tmp.Max_Z : headers[i].Max_Z;

			h_tmp.Min_X = h_tmp.Min_X < headers[i].Min_X ? h_tmp.Min_X : headers[i].Min_X;
			h_tmp.Min_Y = h_tmp.Min_Y < headers[i].Min_Y ? h_tmp.Min_Y : headers[i].Min_Y;
			h_tmp.Min_Z = h_tmp.Min_Z < headers[i].Min_Z ? h_tmp.Min_Z : headers[i].Min_Z;

			h_tmp.X_scale_factor = h_tmp.X_scale_factor > headers[i].X_scale_factor ? h_tmp.X_scale_factor : headers[i].X_scale_factor;
			h_tmp.Y_scale_factor = h_tmp.Y_scale_factor > headers[i].Y_scale_factor ? h_tmp.Y_scale_factor : headers[i].Y_scale_factor;
			h_tmp.Z_scale_factor = h_tmp.Z_scale_factor > headers[i].Z_scale_factor ? h_tmp.Z_scale_factor : headers[i].Z_scale_factor;

			h_tmp.X_offset = h_tmp.X_offset < headers[i].X_offset ? h_tmp.X_offset : headers[i].X_offset;
			h_tmp.Y_offset = h_tmp.Y_offset < headers[i].Y_offset ? h_tmp.Y_offset : headers[i].Y_offset;
			h_tmp.Z_offset = h_tmp.Z_offset < headers[i].Z_offset ? h_tmp.Z_offset : headers[i].Z_offset;

			h_tmp.Legacy_Number_of_point_records += headers[i].Legacy_Number_of_point_records;
			for (int j = 0; j < 5; j++)
			{
				h_tmp.Legacy_Number_of_points_by_return[j] += headers[i].Legacy_Number_of_points_by_return[j];
			}

			h_tmp.Number_of_point_records += headers[i].Number_of_point_records;
			for (int j = 0; j < 15; j++)
			{
				h_tmp.Number_of_points_by_return[j] += headers[i].Number_of_points_by_return[j];
			}

			h_tmp.Number_of_Variable_Length_Records = 0;
			h_tmp.Number_of_Extended_Variable_Length_Records = 0;

			if (version_mode)
			{
				h_tmp.Version_Major = h_tmp.Version_Major > headers[i].Version_Major ? h_tmp.Version_Major : headers[i].Version_Major;
				h_tmp.Version_Minor = h_tmp.Version_Minor > headers[i].Version_Minor ? h_tmp.Version_Minor : headers[i].Version_Minor;
			}
			else
			{
				h_tmp.Version_Major = h_tmp.Version_Major < headers[i].Version_Major ? h_tmp.Version_Major : headers[i].Version_Major;
				h_tmp.Version_Minor = h_tmp.Version_Minor < headers[i].Version_Minor ? h_tmp.Version_Minor : headers[i].Version_Minor;
			}

			unsigned short size_eb_hdr = size_extra_bytes(headers[i].Point_Data_Format_ID, headers[i].Point_Data_Record_Length);
			unsigned short size_eb_tmp = size_extra_bytes(h_tmp.Point_Data_Format_ID, h_tmp.Point_Data_Record_Length);
			if (size_eb_hdr != size_eb_tmp || h_tmp.Number_of_Variable_Length_Records == 0)
			{
				size_eb_tmp = 0;
			}

			if (format_mode)
			{
				h_tmp.Point_Data_Format_ID = Union_Point_Data_Format(h_tmp.Point_Data_Format_ID, headers[i].Point_Data_Format_ID);
			}
			else
			{
				h_tmp.Point_Data_Format_ID = Intersect_Point_Data_Format(h_tmp.Point_Data_Format_ID, headers[i].Point_Data_Format_ID);
			}

			h_tmp.Point_Data_Record_Length = LAS_Point_Size[h_tmp.Point_Data_Format_ID] + size_eb_tmp;

		}

		//offset
		unsigned long long size_pts = h_tmp.Number_of_point_records * h_tmp.Point_Data_Record_Length;
		h_tmp.Offset_to_point_data = h_tmp.Header_Size;
		h_tmp.Start_of_First_Extended_Variable_Length_Record = h_tmp.Offset_to_point_data + size_pts;
	}
	else
	{
		h_tmp = headers[0];
	}

	return h_tmp;
}

LAS_Header LAS_io::Merge_LAS_Headers(char** file_list, int i_file_start, int i_file_end, LAS_VLR* &vlr, LAS_EVLR*& evlr) // the vlrs need to be exactly the same
{
	int n_file = i_file_end - i_file_start + 1;
	LAS_Header* h = Read_LAS_Headers(file_list, i_file_start, i_file_end, false);

	LAS_Header h_tmp = Merge_LAS_Headers(h, n_file, false, false);

	delete[] h;
	h = NULL;

	unsigned int n_vlr = 0;
	unsigned int n_evlr = 0;

	vlr = Read_LAS_VLRs(file_list[i_file_start], n_vlr);
	evlr = Read_LAS_EVLRs(file_list[i_file_start], n_evlr);

	for (int i_file = i_file_start + 1; i_file <= i_file_end; i_file++)
	{
		unsigned int n_vlr_tmp = 0;
		LAS_VLR* vlr_tmp = Read_LAS_VLRs(file_list[i_file], n_vlr_tmp);
		Merge_LAS_VLRs(vlr, n_vlr, vlr_tmp, n_vlr_tmp);
		delete[] vlr_tmp;

		unsigned int n_evlr_tmp = 0;
		LAS_VLR* evlr_tmp = Read_LAS_EVLRs(file_list[i_file], n_evlr_tmp);
		Merge_LAS_VLRs(evlr, n_evlr, evlr_tmp, n_evlr_tmp);
		delete[] evlr_tmp;
	}
	
	h_tmp.Number_of_Variable_Length_Records = n_vlr;

	LAS_io las_tmp;
	las_tmp.las_header = h_tmp;
	las_tmp.las_vlr = vlr;
	las_tmp.las_evlr = evlr;
	las_tmp.Autofill_Header();

	return las_tmp.las_header;
}

LAS_Header LAS_io::Merge_LAS_Headers(char** file_list, int i_file_start, int i_file_end, LAS_VLR*& vlr)
{
	LAS_EVLR* evlr;
	LAS_Header h = Merge_LAS_Headers(file_list, i_file_start, i_file_end, vlr, evlr);
	delete[] evlr;
	return h;
}

LAS_Header LAS_io::Merge_LAS_Headers(char** file_list, int i_file_start, int i_file_end)
{
	LAS_VLR* vlr = NULL;
	LAS_EVLR* evlr = NULL;
	LAS_Header h = Merge_LAS_Headers(file_list, i_file_start, i_file_end, vlr, evlr);
	if (vlr)
	{
		delete[] vlr;
	}
	if (evlr)
	{
		delete[] evlr;
	}
	return h;
}

LAS_VLR* LAS_io::Read_LAS_VLRs(char* file_name, unsigned int &n_vlrs)
{
	LAS_io las_tmp; 
	las_tmp.las_header = Read_LAS_Header(file_name, false);
	las_tmp.Initialize(false, true, false);
	las_tmp.Add_LAS_VLRs(file_name, 0);
	n_vlrs = las_tmp.las_header.Number_of_Variable_Length_Records;
	return las_tmp.las_vlr;
}

LAS_VLR* LAS_io::Read_LAS_EVLRs(char* file_name, unsigned int &n_evlrs)
{
	LAS_io las_tmp;
	las_tmp.las_header = Read_LAS_Header(file_name, false);
	las_tmp.Initialize(false, false, true);
	las_tmp.Add_LAS_EVLRs(file_name, 0);
	n_evlrs = las_tmp.las_header.Number_of_Extended_Variable_Length_Records;
	return las_tmp.las_evlr;
}

void LAS_io::Merge_LAS_VLRs(LAS_VLR* &vlr1, unsigned int &n_vlr1, LAS_VLR* vlr2, unsigned int n_vlr2)
{
	LAS_VLR* vlr = NULL;
	unsigned int n_vlr = 0;
	unsigned int *idx_vlr1 = new unsigned int [min(n_vlr1, n_vlr2)];
	for (unsigned int i_vlr1 = 0; i_vlr1 < n_vlr1; i_vlr1++)
	{
		for (unsigned int i_vlr2 = 0; i_vlr2 < n_vlr2; i_vlr2++)
		{
			if (is_vlr_matched(vlr1[i_vlr1], vlr2[i_vlr2], false))
			{
				idx_vlr1[n_vlr] = i_vlr1;
				n_vlr++;
			}
		}
	}
	
	if (n_vlr > 0)
	{
		vlr = new LAS_VLR[n_vlr];
		for (unsigned int i_vlr = 0; i_vlr < n_vlr; i_vlr++)
		{
			vlr[i_vlr] = vlr1[idx_vlr1[i_vlr]];
		}
	}

	delete[] idx_vlr1;

	if (vlr1)
	{
		delete[] vlr1;
		vlr1 = NULL;
	}
	vlr1 = vlr;
	n_vlr1 = n_vlr;

}

void LAS_io::Merge_LAS_EVLRs(LAS_VLR*& evlr_des, unsigned int& n_evlr_des, LAS_VLR* evlr2, unsigned int n_evlr2)
{
	Merge_LAS_VLRs(evlr_des, n_evlr_des, evlr2, n_evlr2);
}

void LAS_io::Copy_from_LAS_io(LAS_io las_src, bool points, bool vlrs, bool evlrs)
{
	las_header = las_src.las_header;
	las_header.Number_of_point_records = points ? las_header.Number_of_point_records : 0;
	las_header.Number_of_Variable_Length_Records = vlrs ? las_header.Number_of_Variable_Length_Records : 0;
	las_header.Number_of_Extended_Variable_Length_Records = vlrs ? las_header.Number_of_Extended_Variable_Length_Records : 0;

	Initialize(points && las_src.las_points, vlrs && las_src.las_vlr, evlrs && las_src.las_evlr);

	if (points && las_src.las_points)
	{
#pragma omp parallel for schedule (static)
		for (long long i = 0; i < (long long)las_header.Number_of_point_records; i++)
		{
			LAS2LAS_Point(las_src.las_points[i], las_points[i], size_extra_bytes(true));
		}
	}

	if (vlrs && las_src.las_vlr)
	{
#pragma omp parallel for schedule (static)
		for (long long i = 0; i < las_header.Number_of_Variable_Length_Records; i++)
		{
			LAS2LAS_VLR(las_src.las_vlr[i], las_vlr[i]);
		}
	}

	if (evlrs && las_src.las_evlr)
	{
#pragma omp parallel for schedule (static)
		for (long long i = 0; i < las_header.Number_of_Extended_Variable_Length_Records; i++)
		{
			LAS2LAS_EVLR(las_src.las_evlr[i], las_evlr[i]);
		}
	}
}

void LAS_io::Copy_from_LAS_io(LAS_io las_src)
{
	Copy_from_LAS_io(las_src, true, true, true);
}

void LAS_io::Read_LAS_Files(char** file_list, int i_file_start, int i_file_end, bool version_mode, bool format_mode, bool load_points, bool get_extent, bool load_vlrs, bool load_evlrs)
{
	int n_file = i_file_end - i_file_start + 1;
	
	LAS_VLR* vlr_tmp;
	LAS_EVLR* evlr_tmp;
	las_header = Merge_LAS_Headers(file_list, i_file_start, i_file_end, vlr_tmp, evlr_tmp);

	if (load_vlrs)
	{
		las_vlr = vlr_tmp;
	}
	else
	{
		delete[] vlr_tmp;
	}

	if (load_evlrs)
	{
		las_evlr = evlr_tmp;
	}
	else
	{
		delete[] evlr_tmp;
	}

	Initialize(load_points, false, false);

	unsigned long long out_i_start = 0;

	for (int i_file = 0; i_file < n_file; i_file++)
	{
		LAS_Header h = Read_LAS_Header(file_list[i_file + i_file_start], false);

		if (load_points)
		{
			Add_LAS_Points(file_list[i_file + i_file_start], out_i_start, load_points, get_extent);
		}

		out_i_start += h.Number_of_point_records;
	}
}

void LAS_io::Read_LAS_Files(char** file_list, int i_file_start, int i_file_end)
{
	Read_LAS_Files(file_list, i_file_start, i_file_end, LASIO_DEFAULT_VERSION_MODE, LASIO_DEFAULT_FORMAT_MODE, true, true, true, true);
}

void LAS_io::Read_LAS_File(char* file_name, bool load_points, bool get_extent, bool load_vlrs, bool load_evlrs)
{
	Read_LAS_Files(&file_name, 0, 0, LASIO_DEFAULT_VERSION_MODE, LASIO_DEFAULT_FORMAT_MODE, load_points, get_extent, load_vlrs, load_evlrs);
}

void LAS_io::Read_LAS_File(char* file_name)
{
	Read_LAS_Files(&file_name, 0, 0, LASIO_DEFAULT_VERSION_MODE, LASIO_DEFAULT_FORMAT_MODE, true, true, true, true);
}

void LAS_io::Sample_LAS_File(char* file_name, unsigned long long n_sample)
{
	Sample_LAS_Files(&file_name, 0, 0, n_sample);
}

void LAS_io::Sample_LAS_Files(char** file_list, int i_file_start, int i_file_end, unsigned long long n_sample)
{
	int n_file = i_file_end - i_file_start + 1;
	LAS_Header *headers = Read_LAS_Headers(file_list, i_file_start, i_file_end, false);
	LAS_Header header_m = Merge_LAS_Headers(headers, n_file, LASIO_DEFAULT_VERSION_MODE, LASIO_DEFAULT_FORMAT_MODE);

	unsigned long long sample_intvl = (unsigned long long)ceil((double)header_m.Number_of_point_records / n_sample);

	las_header = header_m;
	las_header.Number_of_point_records = n_sample + n_file;
	las_header.Max_X = -DBL_MAX;
	las_header.Max_Y = -DBL_MAX;
	las_header.Max_Z = -DBL_MAX;
	las_header.Min_X = DBL_MAX;
	las_header.Min_Y = DBL_MAX;
	las_header.Min_Z = DBL_MAX;
	Initialize();

	long long i_sample = 0;
	long long ip = 0;
	
	for (int i_file = 0; i_file < n_file; i_file++)
	{
		long long n_spts = (long long)ceil((double)headers[i_file].Number_of_point_records / sample_intvl);

#pragma omp parallel for schedule (static)
		for (long long spid = 0; spid < n_spts; spid++)
		{
			Add_LAS_Points(file_list[i_file + i_file_start], spid * sample_intvl, spid * sample_intvl, i_sample + spid, true, true);
		}

		i_sample += n_spts;
	}

	las_header.Number_of_point_records = i_sample;
	
	for (long long pid = 0; pid < i_sample; pid++)
	{
		double x = x_record2coord(las_points[pid].X);
		double y = y_record2coord(las_points[pid].Y);
		double z = z_record2coord(las_points[pid].Z);

		las_header.Max_X = x > las_header.Max_X ? x : las_header.Max_X;
		las_header.Max_Y = y > las_header.Max_Y ? y : las_header.Max_Y;
		las_header.Max_Z = z > las_header.Max_Z ? z : las_header.Max_Z;
		las_header.Min_X = x < las_header.Min_X ? x : las_header.Min_X;
		las_header.Min_Y = y < las_header.Min_Y ? y : las_header.Min_Y;
		las_header.Min_Z = z < las_header.Min_Z ? z : las_header.Min_Z;
	}
	
	delete[] headers;
	headers = NULL;
}

void LAS_io::Sample_LAS_Points(bool* flag_sample)
{
	unsigned long long npts_sample = 0;
	unsigned long long* pid_sample = new unsigned long long[las_header.Number_of_point_records];
#pragma omp parallel for schedule (static)
	for (long long pid = 0; pid < (long long)las_header.Number_of_point_records; pid++)
	{
		if (flag_sample[pid])
		{
#pragma omp critical
			{
				pid_sample[npts_sample] = pid;
				npts_sample++;
			}
		}
	}

	LAS_Point* point_sample = new LAS_Point[npts_sample];
#pragma omp parallel for schedule (static)
	for (long long i_sample = 0; i_sample < (long long)npts_sample; i_sample++)
	{
		point_sample[i_sample] = las_points[pid_sample[i_sample]];
	}
	delete[] las_points;
	las_points = point_sample;
	las_header.Number_of_point_records = npts_sample;
}

void LAS_io::Add_LAS_Points(char* file_name, unsigned long long out_i_start, bool load_points, bool get_extent)
{
	LAS_Header header = Read_LAS_Header(file_name, false);
	Add_LAS_Points(file_name, 0, header.Number_of_point_records - 1, out_i_start, load_points, get_extent);
}

void LAS_io::Add_LAS_Points(char* file_name, unsigned long long in_i_start, unsigned long long in_i_end, unsigned long long out_i_start, bool load_points, bool get_extent)
{
	if (in_i_start <= in_i_end && (load_points || get_extent))
	{
		LAS_Header header = Read_LAS_Header(file_name, false);

		if (header.Number_of_point_records > 0)
		{
			unsigned long long in_npts = in_i_end - in_i_start + 1;
			unsigned short point_length = header.Point_Data_Record_Length;

			long long n_blocks = (long long)ceil(in_npts / (double)(LASIO_POINT_BLOCK_SIZE_R));

			int* max;
			int* min;

			min = new int[n_blocks * 3];
			max = new int[n_blocks * 3];
			if (get_extent)
			{
#pragma omp parallel for schedule (static)		
				for (long long i_block = 0; i_block < n_blocks; i_block++)
				{
					min[i_block * 3 + 0] = INT_MAX;
					min[i_block * 3 + 1] = INT_MAX;
					min[i_block * 3 + 2] = INT_MAX;

					max[i_block * 3 + 0] = INT_MIN;
					max[i_block * 3 + 1] = INT_MIN;
					max[i_block * 3 + 2] = INT_MIN;
				}
			}

#pragma omp parallel for schedule (static)
			for (long long i_block = 0; i_block < n_blocks; i_block++)
			{
				LAS_Point point;
				FILE *fp = NULL;
				laszip_POINTER laszip_reader = NULL;

				switch (File_Flag(file_name))
				{
					case LAS_FILE_FLAG:
					{
						fp = fopen(file_name, "rb");
						unsigned long long file_offset = (unsigned long long)header.Offset_to_point_data + point_length * (in_i_start + i_block * LASIO_POINT_BLOCK_SIZE_R);
						_fseeki64(fp, file_offset, SEEK_CUR);
						break;
					}
					case LAZ_FILE_FLAG:
					{
						laszip_reader = Open_LAZ_Reader(file_name);
						laszip_seek_point(laszip_reader, in_i_start + i_block * LASIO_POINT_BLOCK_SIZE_R);
						break;
					}
					default:
					break;
				}
				unsigned long long i_p_start = i_block * LASIO_POINT_BLOCK_SIZE_R;
				unsigned long long i_p_end = i_p_start + LASIO_POINT_BLOCK_SIZE_R;
				i_p_end = i_p_end < in_npts ? i_p_end : in_npts - 1;

				for (unsigned long long i_p = i_p_start; i_p <= i_p_end; i_p++)
				{
					switch (File_Flag(file_name))
					{
						case LAS_FILE_FLAG:
						{
							point = Load_A_Point(fp, header.Point_Data_Format_ID);
							break;
						}
						case LAZ_FILE_FLAG:
						{
							point = Load_A_Point(laszip_reader, header.Point_Data_Format_ID);
							break;
						}
						default:
							break;
					}
					
					if (get_extent)
					{
						min[i_block * 3 + 0] = point.X < min[i_block * 3 + 0] ? point.X : min[i_block * 3 + 0];
						min[i_block * 3 + 1] = point.Y < min[i_block * 3 + 1] ? point.Y : min[i_block * 3 + 1];
						min[i_block * 3 + 2] = point.Z < min[i_block * 3 + 2] ? point.Z : min[i_block * 3 + 2];

						max[i_block * 3 + 0] = point.X > max[i_block * 3 + 0] ? point.X : max[i_block * 3 + 0];
						max[i_block * 3 + 1] = point.Y > max[i_block * 3 + 1] ? point.Y : max[i_block * 3 + 1];
						max[i_block * 3 + 2] = point.Z > max[i_block * 3 + 2] ? point.Z : max[i_block * 3 + 2];
					}

					if (load_points)
					{
						LAS2LAS_Point(point, las_points[i_p + out_i_start], size_extra_bytes(true));

						double x = record2coord(header.X_offset, header.X_scale_factor, point.X);
						double y = record2coord(header.Y_offset, header.Y_scale_factor, point.Y);
						double z = record2coord(header.Z_offset, header.Z_scale_factor, point.Z);

						las_points[i_p + out_i_start].X = x_coord2record(x);
						las_points[i_p + out_i_start].Y = y_coord2record(y);
						las_points[i_p + out_i_start].Z = z_coord2record(z);
					}

					if (point.Extra_Bytes)
					{
						delete[] point.Extra_Bytes;
						point.Extra_Bytes = NULL;
					}
				}

				switch (File_Flag(file_name))
				{
					case LAS_FILE_FLAG:
					{
						fclose(fp);
						break;
					}
					case LAZ_FILE_FLAG:
					{
						Close_LAZ_Reader(laszip_reader);
						break;
					}
					default:
						break;
				}
			}
	
			if (get_extent)
			{
				int min_x = INT_MAX;
				int min_y = INT_MAX;
				int min_z = INT_MAX;

				int max_x = INT_MIN;
				int max_y = INT_MIN;
				int max_z = INT_MIN;

				for (long long i_block = 0; i_block < n_blocks; i_block++)
				{
					min_x = min[i_block * 3 + 0] < min_x ? min[i_block * 3 + 0] : min_x;
					min_y = min[i_block * 3 + 1] < min_y ? min[i_block * 3 + 1] : min_y;
					min_z = min[i_block * 3 + 2] < min_z ? min[i_block * 3 + 2] : min_z;

					max_x = max[i_block * 3 + 0] > max_x ? max[i_block * 3 + 0] : max_x;
					max_y = max[i_block * 3 + 1] > max_y ? max[i_block * 3 + 1] : max_y;
					max_z = max[i_block * 3 + 2] > max_z ? max[i_block * 3 + 2] : max_z;
				}

				double Min_X = record2coord(header.X_offset, header.X_scale_factor, min_x);
				double Min_Y = record2coord(header.Y_offset, header.Y_scale_factor, min_y);
				double Min_Z = record2coord(header.Z_offset, header.Z_scale_factor, min_z);

				double Max_X = record2coord(header.X_offset, header.X_scale_factor, max_x);
				double Max_Y = record2coord(header.Y_offset, header.Y_scale_factor, max_y);
				double Max_Z = record2coord(header.Z_offset, header.Z_scale_factor, max_z);


				las_header.Min_X = Min_X < las_header.Min_X ? Min_X : las_header.Min_X;
				las_header.Min_Y = Min_Y < las_header.Min_Y ? Min_Y : las_header.Min_Y;
				las_header.Min_Z = Min_Z < las_header.Min_Z ? Min_Z : las_header.Min_Z;

				las_header.Max_X = Max_X > las_header.Max_X ? Max_X : las_header.Max_X;
				las_header.Max_Y = Max_Y > las_header.Max_Y ? Max_Y : las_header.Max_Y;
				las_header.Max_Z = Max_Z > las_header.Max_Z ? Max_Z : las_header.Max_Z;
			}
			delete[] min;
			delete[] max;
			min = NULL;
			max = NULL;
		}
	}
}

void LAS_io::Add_LAS_VLRs(char* file_name, unsigned int in_i_start, unsigned int in_i_end, unsigned int out_i_start)
{
	if (in_i_start <= in_i_end)
	{
		LAS_Header header = Read_LAS_Header(file_name, false);
		if (header.Number_of_Variable_Length_Records > 0)
		{
			switch (File_Flag(file_name))
			{
				case LAS_FILE_FLAG:
				{
					FILE* fp_in = fopen(file_name, "rb");
					unsigned int i_vlr = out_i_start;
					unsigned int offset = las_header.Header_Size;
					for (unsigned int i = 0; i <= in_i_end; i++)
					{
						LAS_VLR vlr = Load_A_VLR(fp_in, offset);

						if (i >= in_i_start)
						{
							las_vlr[i_vlr] = vlr;
							i_vlr++;
						}
					}
					fclose(fp_in);
					break;
				}
				case LAZ_FILE_FLAG:
				{
					laszip_POINTER laszip_reader = Open_LAZ_Reader(file_name);
					unsigned int i_vlr = out_i_start;
					for (unsigned int i = 0; i <= in_i_end; i++)
					{
						LAS_VLR vlr = Load_A_VLR(laszip_reader, i);

						if (i >= in_i_start)
						{
							las_vlr[i_vlr] = vlr;
							i_vlr++;
						}
					}

					Close_LAZ_Reader(laszip_reader);
					break;
				}
				default:
					break;
			}
		}
	}
}

void LAS_io::Add_LAS_VLRs(char* file_name, unsigned int out_i_start)
{
	LAS_Header header = Read_LAS_Header(file_name, false);		
	Add_LAS_VLRs(file_name, 0, header.Number_of_Variable_Length_Records - 1, out_i_start);
}

void LAS_io::Add_LAS_EVLRs(char* file_name, unsigned int in_i_start, unsigned int in_i_end, unsigned int out_i_start)
{
	LAS_Header header = Read_LAS_Header(file_name, false);
	if (header.Number_of_Extended_Variable_Length_Records > 0)
	{
		switch (File_Flag(file_name))
		{
			case LAS_FILE_FLAG:
			{
				FILE* fp_in = fopen(file_name, "rb");
				unsigned long long offset = las_header.Start_of_First_Extended_Variable_Length_Record;

				unsigned int i_evlr = out_i_start;
				for (unsigned int i = 0; i <= in_i_end; i++)
				{
					LAS_EVLR evlr = Load_A_EVLR(fp_in, offset);
					if (i >= in_i_start)
					{
						las_evlr[i_evlr] = evlr;
						i_evlr++;
					}
				}
				fclose(fp_in);
				break;
			}

			default:
				break;
		}
	}
}

void LAS_io::Add_LAS_EVLRs(char* file_name, unsigned int out_i_start)
{
	LAS_Header header = Read_LAS_Header(file_name, false);
	Add_LAS_EVLRs(file_name, 0, header.Number_of_Extended_Variable_Length_Records - 1, out_i_start);
}

LAS_Header LAS_io::Load_LAS_Header(FILE* fp_in)
{
	unsigned char Version_Major;
	rewind(fp_in);
	_fseeki64(fp_in, LASIO_OFFSET_TO_VERSION_MAJOR, SEEK_CUR);

	fread(&Version_Major, sizeof(unsigned char), 1, fp_in);

	switch (Version_Major)
	{
		case 1:
			return Load_LAS_Header_1_X(fp_in);
			break;

		default:
			printf("\ninvalid version!");
			getchar();
			LAS_Header h_empty;
			return h_empty;
			break;
	}
}

LAS_Header LAS_io::Load_LAS_Header_1_X(FILE* fp_in)
{
	rewind(fp_in);
	LAS_Header header;

	fread(&header.File_Signature, sizeof(header.File_Signature), 1, fp_in);
	fread(&header.File_Source_ID, sizeof(header.File_Source_ID), 1, fp_in);
	fread(&header.Global_Encoding, sizeof(header.Global_Encoding), 1, fp_in);
	fread(&header.Project_ID_GUID_data_1, sizeof(header.Project_ID_GUID_data_1), 1, fp_in);
	fread(&header.Project_ID_GUID_data_2, sizeof(header.Project_ID_GUID_data_2), 1, fp_in);
	fread(&header.Project_ID_GUID_data_3, sizeof(header.Project_ID_GUID_data_3), 1, fp_in);
	fread(header.Project_ID_GUID_data_4, sizeof(header.Project_ID_GUID_data_4), 1, fp_in);
	fread(&header.Version_Major, sizeof(header.Version_Major), 1, fp_in);
	fread(&header.Version_Minor, sizeof(header.Version_Minor), 1, fp_in);
	fread(header.System_Identifier, sizeof(header.System_Identifier), 1, fp_in);
	fread(header.Generating_Software, sizeof(header.Generating_Software), 1, fp_in);
	fread(&header.File_Creation_Day_of_Year, sizeof(header.File_Creation_Day_of_Year), 1, fp_in);
	fread(&header.File_Creation_Year, sizeof(header.File_Creation_Year), 1, fp_in);
	fread(&header.Header_Size, sizeof(header.Header_Size), 1, fp_in);
	fread(&header.Offset_to_point_data, sizeof(header.Offset_to_point_data), 1, fp_in);
	fread(&header.Number_of_Variable_Length_Records, sizeof(header.Number_of_Variable_Length_Records), 1, fp_in);

	fread(&header.Point_Data_Format_ID, sizeof(header.Point_Data_Format_ID), 1, fp_in);
	fread(&header.Point_Data_Record_Length, sizeof(header.Point_Data_Record_Length), 1, fp_in);

	//unsigned char id_from_length = Match_Point_Data_Format(header.Point_Data_Record_Length);
	//if (header.Point_Data_Format_ID != id_from_length)
	//{
	//	header.Point_Data_Format_ID = id_from_length;
	//}

	fread(&header.Legacy_Number_of_point_records, sizeof(header.Legacy_Number_of_point_records), 1, fp_in);	
	fread(header.Legacy_Number_of_points_by_return, sizeof(header.Legacy_Number_of_points_by_return), 1, fp_in);

	header.Number_of_point_records = header.Legacy_Number_of_point_records;
	header.Legacy_Number_of_point_records = 0;

	for (int i = 0; i < 5; i++)
	{
		header.Number_of_points_by_return[i] = header.Legacy_Number_of_points_by_return[i];
		header.Legacy_Number_of_points_by_return[i] = 0;
	}

	fread(&header.X_scale_factor, sizeof(header.X_scale_factor), 1, fp_in);
	fread(&header.Y_scale_factor, sizeof(header.Y_scale_factor), 1, fp_in);
	fread(&header.Z_scale_factor, sizeof(header.Z_scale_factor), 1, fp_in);

	fread(&header.X_offset, sizeof(header.X_offset), 1, fp_in);
	fread(&header.Y_offset, sizeof(header.Y_offset), 1, fp_in);
	fread(&header.Z_offset, sizeof(header.Z_offset), 1, fp_in);

	fread(&header.Max_X, sizeof(header.Max_X), 1, fp_in);
	fread(&header.Min_X, sizeof(header.Min_X), 1, fp_in);
	fread(&header.Max_Y, sizeof(header.Max_Y), 1, fp_in);
	fread(&header.Min_Y, sizeof(header.Min_Y), 1, fp_in);
	fread(&header.Max_Z, sizeof(header.Max_Z), 1, fp_in);
	fread(&header.Min_Z, sizeof(header.Min_Z), 1, fp_in);

	if (header.Version_Minor >= 3)
	{
		fread(&header.Start_of_Waveform_Data_Packet_Record, sizeof(header.Start_of_Waveform_Data_Packet_Record), 1, fp_in);
	}

	if (header.Version_Minor >= 4)
	{
		fread(&header.Start_of_First_Extended_Variable_Length_Record, sizeof(header.Start_of_First_Extended_Variable_Length_Record), 1, fp_in);
		fread(&header.Number_of_Extended_Variable_Length_Records, sizeof(header.Number_of_Extended_Variable_Length_Records), 1, fp_in);

		fread(&header.Number_of_point_records, sizeof(header.Number_of_point_records), 1, fp_in);
		fread(header.Number_of_points_by_return, sizeof(header.Number_of_points_by_return), 1, fp_in);
	}

	return header;
}

LAS_Point LAS_io::Load_A_Point(FILE* fp_in, unsigned short Point_Format_ID)
{
	LAS_Point p;
	const bool* opt_field = LAS_Optional_Field[Point_Format_ID];

	fread(&p.X, sizeof(p.X), 1, fp_in);
	fread(&p.Y, sizeof(p.Y), 1, fp_in);
	fread(&p.Z, sizeof(p.Z), 1, fp_in);
	fread(&p.Intensity, sizeof(p.Intensity), 1, fp_in);

	if (opt_field[0])
	{
		unsigned char Bit8_1, Bit8_2;
		fread(&Bit8_1, sizeof(Bit8_1), 1, fp_in);
		fread(&Bit8_2, sizeof(Bit8_2), 1, fp_in);

		p.Return_Number = (Bit8_1 >> 0) & 15;
		p.Number_of_Returns = (Bit8_1 >> 4) & 15;

		p.Classification_Flags = (Bit8_2 >> 0) & 15;
		p.Scanner_Channel = (Bit8_2 >> 4) & 3;
		p.Scan_Direction_Flag = (Bit8_2 >> 6) & 1;
		p.Edge_of_Flight_Line = (Bit8_2 >> 7) & 1;
	}
	else
	{
		unsigned char Bit8_1;
		fread(&Bit8_1, sizeof(Bit8_1), 1, fp_in);

		p.Return_Number = (Bit8_1 >> 0) & 7;
		p.Number_of_Returns = (Bit8_1 >> 3) & 7;
		p.Scan_Direction_Flag = (Bit8_1 >> 6) & 1;
		p.Edge_of_Flight_Line = (Bit8_1 >> 7) & 1;
	}

	fread(&p.Classification, sizeof(p.Classification), 1, fp_in);

	if (opt_field[1])
	{
		fread(&p.Scan_Angle_Rank, sizeof(p.Scan_Angle_Rank), 1, fp_in);
	}
	else
	{
		char scan_angle_rank = (char)p.Scan_Angle_Rank;
		fread(&scan_angle_rank, sizeof((char)scan_angle_rank), 1, fp_in);
	}

	fread(&p.User_Data, sizeof(p.User_Data), 1, fp_in);
	fread(&p.Point_Source_ID, sizeof(p.Point_Source_ID), 1, fp_in);

	if (opt_field[2])
	{
		fread(&p.GPS_Time, sizeof(p.GPS_Time), 1, fp_in);
	}

	if (opt_field[3])
	{
		fread(&p.Red, sizeof(p.Red), 1, fp_in);
		fread(&p.Green, sizeof(p.Green), 1, fp_in);
		fread(&p.Blue, sizeof(p.Blue), 1, fp_in);
	}

	if (opt_field[4])
	{
		fread(&p.NIR, sizeof(p.NIR), 1, fp_in);
	}

	if (opt_field[5])
	{
		fread(&p.Wave_Packet_Discriptor_Index, sizeof(p.Wave_Packet_Discriptor_Index), 1, fp_in);
		fread(&p.Byte_offset_to_waveform_data, sizeof(p.Byte_offset_to_waveform_data), 1, fp_in);
		fread(&p.Waveform_packet_size_in_bytes, sizeof(p.Waveform_packet_size_in_bytes), 1, fp_in);
		fread(&p.Return_Point_Waveform_Location, sizeof(p.Return_Point_Waveform_Location), 1, fp_in);
		fread(&p.Xt, sizeof(p.Xt), 1, fp_in);
		fread(&p.Yt, sizeof(p.Yt), 1, fp_in);
		fread(&p.Zt, sizeof(p.Zt), 1, fp_in);
	}

	if (size_extra_bytes(Point_Format_ID, las_header.Point_Data_Record_Length) > 0)
	{
		if (!p.Extra_Bytes)
		{
			p.Extra_Bytes = new LAS_BYTE[size_extra_bytes(Point_Format_ID, las_header.Point_Data_Record_Length)];
		}
		fread(p.Extra_Bytes, sizeof(LAS_BYTE), size_extra_bytes(Point_Format_ID, las_header.Point_Data_Record_Length), fp_in);
	}

	return p;
}

LAS_VLR LAS_io::Load_A_VLR(FILE* fp_in, unsigned int &offset)
{
	rewind(fp_in);
	_fseeki64(fp_in, offset, SEEK_CUR);

	LAS_VLR vlr;
	fread(&vlr.vlr_header.Reserved, sizeof(vlr.vlr_header.Reserved), 1, fp_in);
	fread(vlr.vlr_header.User_ID, sizeof(vlr.vlr_header.User_ID), 1, fp_in);
	fread(&vlr.vlr_header.Record_ID, sizeof(vlr.vlr_header.Record_ID), 1, fp_in);
	fread(&vlr.vlr_header.Record_Length_After_Header, sizeof(vlr.vlr_header.Record_Length_After_Header), 1, fp_in);
	fread(vlr.vlr_header.Description, sizeof(vlr.vlr_header.Description), 1, fp_in);
	if (vlr.vlr_header.Record_Length_After_Header > 0)
	{
		vlr.vlr_buffer = new LAS_BYTE[vlr.vlr_header.Record_Length_After_Header];
		fread(vlr.vlr_buffer, sizeof(LAS_BYTE), vlr.vlr_header.Record_Length_After_Header, fp_in);
	}

	offset = ftell(fp_in);
	return vlr;
}

LAS_EVLR LAS_io::Load_A_EVLR(FILE* fp_in, unsigned long long &offset)
{
	rewind(fp_in);
	_fseeki64(fp_in, offset, SEEK_CUR);

	LAS_EVLR evlr;
	fread(&evlr.vlr_header.Reserved, sizeof(evlr.vlr_header.Reserved), 1, fp_in);
	fread(evlr.vlr_header.User_ID, sizeof(evlr.vlr_header.User_ID), 1, fp_in);
	fread(&evlr.vlr_header.Record_ID, sizeof(evlr.vlr_header.Record_ID), 1, fp_in);
	fread(&evlr.vlr_header.Record_Length_After_Header, sizeof(evlr.vlr_header.Record_Length_After_Header), 1, fp_in);
	fread(evlr.vlr_header.Description, sizeof(evlr.vlr_header.Description), 1, fp_in);
	if (evlr.vlr_header.Record_Length_After_Header > 0)
	{
		evlr.vlr_buffer = new LAS_BYTE[evlr.vlr_header.Record_Length_After_Header];
		fread(evlr.vlr_buffer, sizeof(LAS_BYTE), evlr.vlr_header.Record_Length_After_Header, fp_in);
	}

	offset = ftell(fp_in);
	return evlr;
}

void LAS_io::Sort_by_Time()
{
	Sort_LAS_by_GPS_Time(las_points, las_header.Number_of_point_records);
}

void LAS_io::Sort_by_Source_ID()
{
	Sort_LAS_by_Point_Source_ID(las_points, las_header.Number_of_point_records);
}

void LAS_io::Sort_by_Classification()
{
	Sort_LAS_by_Classification(las_points, las_header.Number_of_point_records);
}

void LAS_io::Write_LAS(char* file_name)
{
	Autofill_Header();		
	switch (File_Flag(file_name))
	{
		case LAS_FILE_FLAG:
		{			
			FILE* fp_out = fopen(file_name, "wb");
			Write_LAS_Header(fp_out);
			Write_LAS_VLRs(fp_out);	
			Write_LAS_Points(fp_out);
			Write_LAS_EVLRs(fp_out);

			fclose(fp_out);
			break;
		}
		case LAZ_FILE_FLAG:
		{
			// create the writer
			Write_LAZ_File(file_name);
			break;
		}
		default:
			break;
	}
}

void LAS_io::Write_LAS_Header(FILE* fp_out)
{
	switch (las_header.Version_Major)
	{
		case 1:
			Write_LAS_Header_1_X(fp_out);
			break;
	}
}

void LAS_io::Write_LAS_Header_1_X(FILE* fp_out)
{
	rewind(fp_out);
	
	fwrite(&las_header.File_Signature, sizeof(las_header.File_Signature), 1, fp_out);
	fwrite(&las_header.File_Source_ID, sizeof(las_header.File_Source_ID), 1, fp_out);
	fwrite(&las_header.Global_Encoding, sizeof(las_header.Global_Encoding), 1, fp_out);
	fwrite(&las_header.Project_ID_GUID_data_1, sizeof(las_header.Project_ID_GUID_data_1), 1, fp_out);
	fwrite(&las_header.Project_ID_GUID_data_2, sizeof(las_header.Project_ID_GUID_data_2), 1, fp_out);
	fwrite(&las_header.Project_ID_GUID_data_3, sizeof(las_header.Project_ID_GUID_data_3), 1, fp_out);
	fwrite(las_header.Project_ID_GUID_data_4, sizeof(las_header.Project_ID_GUID_data_4), 1, fp_out);
	fwrite(&las_header.Version_Major, sizeof(las_header.Version_Major), 1, fp_out);
	fwrite(&las_header.Version_Minor, sizeof(las_header.Version_Minor), 1, fp_out);
	fwrite(las_header.System_Identifier, sizeof(las_header.System_Identifier), 1, fp_out);
	fwrite(las_header.Generating_Software, sizeof(las_header.Generating_Software), 1, fp_out);
	fwrite(&las_header.File_Creation_Day_of_Year, sizeof(las_header.File_Creation_Day_of_Year), 1, fp_out);
	fwrite(&las_header.File_Creation_Year, sizeof(las_header.File_Creation_Year), 1, fp_out);
	fwrite(&las_header.Header_Size, sizeof(las_header.Header_Size), 1, fp_out);
	fwrite(&las_header.Offset_to_point_data, sizeof(las_header.Offset_to_point_data), 1, fp_out);
	fwrite(&las_header.Number_of_Variable_Length_Records, sizeof(las_header.Number_of_Variable_Length_Records), 1, fp_out);

	fwrite(&las_header.Point_Data_Format_ID, sizeof(las_header.Point_Data_Format_ID), 1, fp_out);
	fwrite(&las_header.Point_Data_Record_Length, sizeof(las_header.Point_Data_Record_Length), 1, fp_out);

	if (las_header.Version_Minor < 4)
	{
		fwrite(&las_header.Number_of_point_records, sizeof(las_header.Legacy_Number_of_point_records), 1, fp_out);
		fwrite(las_header.Number_of_points_by_return, sizeof(las_header.Legacy_Number_of_points_by_return), 1, fp_out);
	}
	else
	{
		unsigned int Legacy_Number_of_point_records = 0;
		unsigned int Legacy_Number_of_points_by_return[5] = { 0, 0, 0, 0, 0 };

		fwrite(&las_header.Legacy_Number_of_point_records, sizeof(las_header.Legacy_Number_of_point_records), 1, fp_out);
		fwrite(las_header.Legacy_Number_of_points_by_return, sizeof(las_header.Legacy_Number_of_points_by_return), 1, fp_out);
	}

	fwrite(&las_header.X_scale_factor, sizeof(las_header.X_scale_factor), 1, fp_out);
	fwrite(&las_header.Y_scale_factor, sizeof(las_header.Y_scale_factor), 1, fp_out);
	fwrite(&las_header.Z_scale_factor, sizeof(las_header.Z_scale_factor), 1, fp_out);
	fwrite(&las_header.X_offset, sizeof(las_header.X_offset), 1, fp_out);
	fwrite(&las_header.Y_offset, sizeof(las_header.Y_offset), 1, fp_out);
	fwrite(&las_header.Z_offset, sizeof(las_header.Z_offset), 1, fp_out);

	fwrite(&las_header.Max_X, sizeof(las_header.Max_X), 1, fp_out);
	fwrite(&las_header.Min_X, sizeof(las_header.Min_X), 1, fp_out);
	fwrite(&las_header.Max_Y, sizeof(las_header.Max_Y), 1, fp_out);
	fwrite(&las_header.Min_Y, sizeof(las_header.Min_Y), 1, fp_out);
	fwrite(&las_header.Max_Z, sizeof(las_header.Max_Z), 1, fp_out);
	fwrite(&las_header.Min_Z, sizeof(las_header.Min_Z), 1, fp_out);

	if (las_header.Version_Minor >= 3)
	{
		fwrite(&las_header.Start_of_Waveform_Data_Packet_Record, sizeof(las_header.Start_of_Waveform_Data_Packet_Record), 1, fp_out);
	}

	if (las_header.Version_Minor >= 4)
	{
		fwrite(&las_header.Start_of_First_Extended_Variable_Length_Record, sizeof(las_header.Start_of_First_Extended_Variable_Length_Record), 1, fp_out);
		fwrite(&las_header.Number_of_Extended_Variable_Length_Records, sizeof(las_header.Number_of_Extended_Variable_Length_Records), 1, fp_out);

		fwrite(&las_header.Number_of_point_records, sizeof(las_header.Number_of_point_records), 1, fp_out);
		fwrite(las_header.Number_of_points_by_return, sizeof(las_header.Number_of_points_by_return), 1, fp_out);
	}
}

void LAS_io::Write_LAS_Points(FILE* fp_out)
{
	long long n_blocks = (long long)ceil(las_header.Number_of_point_records / (double)LASIO_POINT_BLOCK_SIZE_W);
	unsigned long long block_bytes = LASIO_POINT_BLOCK_SIZE_W * las_header.Point_Data_Record_Length;

#pragma omp parallel for schedule (static)
	for (long long i_block = 0; i_block < n_blocks; i_block++)
	{
		LAS_BYTE* block = new LAS_BYTE[block_bytes];
		unsigned long long offset_block = 0;
		unsigned long long i_p_start = i_block * LASIO_POINT_BLOCK_SIZE_W;
		unsigned long long i_p_end = i_p_start + LASIO_POINT_BLOCK_SIZE_W - 1;
		i_p_end = i_p_end < las_header.Number_of_point_records ? i_p_end : las_header.Number_of_point_records - 1;

		for (unsigned long long i_p = i_p_start; i_p <= i_p_end; i_p++)
		{
			Point_to_Block(block, offset_block, las_points[i_p]);
		}
		
#pragma omp critical
		{
			rewind(fp_out);
			unsigned long long file_offset = (unsigned long long)las_header.Offset_to_point_data + i_block * block_bytes;
			_fseeki64(fp_out, file_offset, SEEK_CUR);
			fwrite(block, las_header.Point_Data_Record_Length, i_p_end - i_p_start + 1, fp_out);
		}
		delete[] block;
		block = NULL;
	}
}

void LAS_io::Write_LAS_VLRs(FILE* fp_out)
{
	rewind(fp_out);
	unsigned int offset = las_header.Header_Size;
	for (unsigned int i = 0; i < las_header.Number_of_Variable_Length_Records; i++)
	{
		Write_A_VLR(fp_out, offset, las_vlr[i]);
	}
}

void LAS_io::Write_LAS_EVLRs(FILE* fp_out)
{
	rewind(fp_out);
	unsigned long long offset = las_header.Start_of_First_Extended_Variable_Length_Record;
	for (unsigned int i = 0; i < las_header.Number_of_Extended_Variable_Length_Records; i++)
	{
		Write_A_EVLR(fp_out, offset, las_evlr[i]);
	}
}

void LAS_io::Point_to_Block(LAS_BYTE *point_bytes, unsigned long long &offset, LAS_Point p)
{
	const bool* opt_field = LAS_Optional_Field[las_header.Point_Data_Format_ID];

	Data2Bytes(point_bytes, &p.X, sizeof(p.X), 1, offset);
	Data2Bytes(point_bytes, &p.Y, sizeof(p.Y), 1, offset);
	Data2Bytes(point_bytes, &p.Z, sizeof(p.Z), 1, offset);
	Data2Bytes(point_bytes, &p.Intensity, sizeof(p.Intensity), 1, offset);

	if (opt_field[0])
	{
		LAS_BYTE Bit8_1 = (p.Return_Number >> 0)
			| (p.Number_of_Returns >> 4);
		LAS_BYTE Bit8_2 = (p.Classification_Flags >> 0)
			| (p.Scanner_Channel >> 4)
			| (p.Scan_Direction_Flag >> 6)
			| (p.Edge_of_Flight_Line >> 7);

		Data2Bytes(point_bytes, &Bit8_1, sizeof(Bit8_1), 1, offset);
		Data2Bytes(point_bytes, &Bit8_2, sizeof(Bit8_2), 1, offset);
	}
	else
	{
		LAS_BYTE Bit8_1 = (p.Return_Number >> 0)
			| (p.Number_of_Returns >> 3)
			| (p.Scan_Direction_Flag >> 6)
			| (p.Edge_of_Flight_Line >> 7);

		Data2Bytes(point_bytes, &Bit8_1, sizeof(Bit8_1), 1, offset);
	}

	Data2Bytes(point_bytes, &p.Classification, sizeof(p.Classification), 1, offset);

	if (opt_field[1])
	{
		Data2Bytes(point_bytes, &p.Scan_Angle_Rank, sizeof(p.Scan_Angle_Rank), 1, offset);
	}
	else
	{
		char scan_angle_rank = (char)p.Scan_Angle_Rank;
		Data2Bytes(point_bytes, &scan_angle_rank, sizeof(scan_angle_rank), 1, offset); 
	}
	Data2Bytes(point_bytes, &p.User_Data, sizeof(p.User_Data), 1, offset);
	Data2Bytes(point_bytes, &p.Point_Source_ID, sizeof(p.Point_Source_ID), 1, offset);

	if (opt_field[2])
	{
		Data2Bytes(point_bytes, &p.GPS_Time, sizeof(p.GPS_Time), 1, offset);
	}

	if (opt_field[3])
	{
		Data2Bytes(point_bytes, &p.Red, sizeof(p.Red), 1, offset);
		Data2Bytes(point_bytes, &p.Green, sizeof(p.Green), 1, offset);
		Data2Bytes(point_bytes, &p.Blue, sizeof(p.Blue), 1, offset);
	}

	if (opt_field[4])
	{
		Data2Bytes(point_bytes, &p.NIR, sizeof(p.NIR), 1, offset);
	}

	if (opt_field[5])
	{
		Data2Bytes(point_bytes, &p.Wave_Packet_Discriptor_Index, sizeof(p.Wave_Packet_Discriptor_Index), 1, offset);
		Data2Bytes(point_bytes, &p.Byte_offset_to_waveform_data, sizeof(p.Byte_offset_to_waveform_data), 1, offset);
		Data2Bytes(point_bytes, &p.Waveform_packet_size_in_bytes, sizeof(p.Waveform_packet_size_in_bytes), 1, offset);
		Data2Bytes(point_bytes, &p.Return_Point_Waveform_Location, sizeof(p.Return_Point_Waveform_Location), 1, offset);
		Data2Bytes(point_bytes, &p.Xt, sizeof(p.Xt), 1, offset);
		Data2Bytes(point_bytes, &p.Yt, sizeof(p.Yt), 1, offset);
		Data2Bytes(point_bytes, &p.Zt, sizeof(p.Zt), 1, offset);
	}

	if (size_extra_bytes(true) > 0)
	{
		Data2Bytes(point_bytes, p.Extra_Bytes, sizeof(LAS_BYTE), size_extra_bytes(true), offset);
	}
}

void LAS_io::Write_A_Point(FILE* fp_out, LAS_Point p)
{
	LAS_BYTE* point_bytes = new LAS_BYTE[las_header.Point_Data_Record_Length];
	unsigned long long pos = 0;
	Point_to_Block(point_bytes, pos, p);
	fwrite(point_bytes, 1, las_header.Point_Data_Record_Length, fp_out);
	delete[] point_bytes;
	point_bytes = NULL;
}

void LAS_io::Write_A_VLR(FILE* fp_out, unsigned int &offset, LAS_VLR vlr)
{
	rewind(fp_out);
	_fseeki64(fp_out, offset, SEEK_CUR);
	fwrite(&vlr.vlr_header.Reserved, sizeof(vlr.vlr_header.Reserved), 1, fp_out);
	fwrite(vlr.vlr_header.User_ID, sizeof(vlr.vlr_header.User_ID), 1, fp_out);
	fwrite(&vlr.vlr_header.Record_ID, sizeof(vlr.vlr_header.Record_ID), 1, fp_out);
	fwrite(&vlr.vlr_header.Record_Length_After_Header, sizeof(vlr.vlr_header.Record_Length_After_Header), 1, fp_out);
	fwrite(vlr.vlr_header.Description, sizeof(vlr.vlr_header.Description), 1, fp_out);

	fwrite(vlr.vlr_buffer, 1, vlr.vlr_header.Record_Length_After_Header, fp_out);
	
	offset = ftell(fp_out);
}

void LAS_io::Write_A_EVLR(FILE* fp_out, unsigned long long &offset, LAS_EVLR evlr)
{
	rewind(fp_out);
	_fseeki64(fp_out, offset, SEEK_CUR);

	fwrite(&evlr.vlr_header.Reserved, sizeof(evlr.vlr_header.Reserved), 1, fp_out);
	fwrite(evlr.vlr_header.User_ID, sizeof(evlr.vlr_header.User_ID), 1, fp_out);
	fwrite(&evlr.vlr_header.Record_ID, sizeof(evlr.vlr_header.Record_ID), 1, fp_out);
	fwrite(&evlr.vlr_header.Record_Length_After_Header, sizeof(evlr.vlr_header.Record_Length_After_Header), 1, fp_out);
	fwrite(evlr.vlr_header.Description, sizeof(evlr.vlr_header.Description), 1, fp_out);

	for (unsigned int i = 0; i < evlr.vlr_header.Record_Length_After_Header; i++)
	{
		fwrite(&evlr.vlr_buffer[i], 1, 1, fp_out);
	}

	offset = ftell(fp_out);
}


///////////  vlr
LAS_VLR LAS_io::create_a_vlr(unsigned short record_length_after_header)
{
	LAS_VLR vlr;
	vlr.vlr_header.Record_Length_After_Header = record_length_after_header;
	if (record_length_after_header > 0)
	{
		vlr.vlr_buffer = new LAS_BYTE[record_length_after_header];
	}
	return vlr;
}

LAS_VLR LAS_io::create_a_vlr(char* user_id, unsigned short record_id, unsigned short record_length_after_header)
{
	LAS_VLR vlr = create_a_vlr(record_length_after_header);
	Set_VLR_User_ID(vlr, user_id);
	Set_VLR_Record_ID(vlr, record_id);
	return vlr;
}

void LAS_io::Set_VLR_User_ID(LAS_VLR &vlr, char* user_id)
{
	size_t len = strlen(user_id);
	memcpy(vlr.vlr_header.User_ID, user_id, min(len, (size_t)16));
	if (len < 16)
	{
		vlr.vlr_header.User_ID[len] = '\0';
	}
}

void LAS_io::Set_VLR_Record_ID(LAS_VLR &vlr, unsigned short record_id)
{
	vlr.vlr_header.Record_ID = record_id;
}

void LAS_io::Set_VLR_Description(LAS_VLR &vlr, char* description)
{
	size_t len = strlen(description);
	memcpy(vlr.vlr_header.Description, description, min(len, (size_t)32));
	if (len < 32)
	{
		vlr.vlr_header.Description[len] = '\0';
	}
}

void LAS_io::Set_VLR_Bytes(LAS_VLR &vlr, LAS_BYTE* bytes, unsigned short size_bytes)
{
	if (vlr.vlr_header.Record_Length_After_Header < size_bytes)
	{
		if (!vlr.vlr_buffer)
		{
			delete[] vlr.vlr_buffer;
			vlr.vlr_buffer = NULL;
		}
		vlr.vlr_header.Record_Length_After_Header = size_bytes;
		vlr.vlr_buffer = new LAS_BYTE[size_bytes];
	}
	memcpy(vlr.vlr_buffer, bytes, size_bytes);
}

void LAS_io::Add_New_VLRs(LAS_VLR* vlr, unsigned int n_vlr)
{
	if (las_header.Number_of_Variable_Length_Records == 0)
	{
		las_vlr = new LAS_VLR[n_vlr];
		for (unsigned int idx = 0; idx < n_vlr; idx++)
		{
			las_vlr[idx].vlr_header = vlr[idx].vlr_header;
			if (vlr[idx].vlr_header.Record_Length_After_Header > 0)
			{
				las_vlr[idx].vlr_buffer = new LAS_BYTE[vlr[idx].vlr_header.Record_Length_After_Header];
			}
		}
	}
	else
	{
		LAS_VLR *vlr_cpy = las_vlr;
		las_vlr = new LAS_VLR[las_header.Number_of_Variable_Length_Records + n_vlr];
		
#pragma omp parallel for schedule(static)
		for (long long idx = 0; idx < las_header.Number_of_Variable_Length_Records; idx++)
		{
			las_vlr[idx].vlr_header = vlr_cpy[idx].vlr_header;
			las_vlr[idx].vlr_buffer = vlr_cpy[idx].vlr_buffer;
		}

#pragma omp parallel for schedule(static)
		for (long long idx = 0; idx < n_vlr; idx++)
		{
			unsigned int i_vlr = (unsigned int)idx + las_header.Number_of_Variable_Length_Records;
			las_vlr[i_vlr].vlr_header = vlr[idx].vlr_header;
			if (vlr[i_vlr].vlr_header.Record_Length_After_Header > 0)
			{
				las_vlr[i_vlr].vlr_buffer = new LAS_BYTE[vlr[i_vlr].vlr_header.Record_Length_After_Header];
			}
		}
		delete[] vlr_cpy;
	}
	las_header.Number_of_Variable_Length_Records += n_vlr;
}

void LAS_io::Add_New_VLR(LAS_VLR vlr)
{
	Add_New_VLRs(&vlr, 1);
}

bool LAS_io::is_vlr_matched(LAS_VLR vlr_1, LAS_VLR vlr_2, bool header_only)
{
	bool flag = false;
	if (EZPC::c_str_cmpi(vlr_1.vlr_header.User_ID, vlr_2.vlr_header.User_ID) == 0
		&& vlr_1.vlr_header.Record_ID == vlr_2.vlr_header.Record_ID
		&& vlr_1.vlr_header.Record_Length_After_Header == vlr_2.vlr_header.Record_Length_After_Header)
	{
		if (header_only || memcmp(vlr_1.vlr_buffer, vlr_2.vlr_buffer, vlr_1.vlr_header.Record_Length_After_Header) == 0)
		{
			flag = true;
		}
	}
	return flag;
}

LAS_VLR LAS_io::OGC_Math_Transform_WKT_Record_VLR()
{
	char lasf_prj[] = "LASF_Projection";
	return create_a_vlr(lasf_prj, 2111, 0);
}

LAS_VLR LAS_io::OGC_Coordinate_System_WKT_Record_VLR()
{
	char lasf_prj[] = "LASF_Projection";
	return create_a_vlr(lasf_prj, 2112, 0);
}

LAS_VLR LAS_io::GeoKeyDirectoryTag_Record_VLR()
{
	char lasf_prj[] = "LASF_Projection";
	return create_a_vlr(lasf_prj, 34735, sizeof(sGeoKeys));
}

LAS_VLR LAS_io::GeoDoubleParamsTag_Record_VLR()
{
	char lasf_prj[] = "LASF_Projection";
	return create_a_vlr(lasf_prj, 34736, 0);
}

LAS_VLR LAS_io::GeoAsciiParamsTag_Record_VLR()
{
	char lasf_prj[] = "LASF_Projection";
	return create_a_vlr(lasf_prj, 34737, 0);
}

LAS_VLR LAS_io::Classification_Lookup_VLR()
{
	char lasf_spec[] = "LASF_Spec";
	return create_a_vlr(lasf_spec, 0, 0);
}

LAS_VLR LAS_io::Text_Area_Description_VLR()
{
	char lasf_spec[] = "LASF_Spec";
	return create_a_vlr(lasf_spec, 3, 0);
}

LAS_VLR LAS_io::Extra_Bytes_VLR()
{
	char lasf_spec[] = "LASF_Spec";
	return create_a_vlr(lasf_spec, 4, 0);
}

LAS_VLR LAS_io::Superseded_VLR()
{
	char lasf_spec[] = "LASF_Spec";
	return create_a_vlr(lasf_spec, 7, 0);
}

LAS_VLR LAS_io::Waveform_Packet_Descriptor_VLR(unsigned short idx)
{
	char lasf_spec[] = "LASF_Spec";
	return create_a_vlr(lasf_spec, 99 + idx, sizeof(Waveform_Packet_Descriptor));
}

LAS_EVLR LAS_io::Waveform_Data_Packets_EVLR()
{
	char lasf_spec[] = "LASF_Spec";
	return create_a_vlr(lasf_spec, 65535, 0);
}

////////////////////////// extra bytes

unsigned short LAS_io::size_extra_bytes(bool header_or_vlr) // how many bytes that are extra) // how many bytes that are extra from header
{
	unsigned short size = 0;
	if (header_or_vlr)
	{
		size = size_extra_bytes(las_header.Point_Data_Format_ID, las_header.Point_Data_Record_Length);
	}
	else
	{
		unsigned short n_eb;
		LAS_Extra_Bytes_Info* eb_info;
		Get_All_Extra_Bytes_Info(eb_info, n_eb);
		if (n_eb > 0)
		{
			size = size_extra_bytes(eb_info, n_eb);
			delete[] eb_info;
		}
	}
	return size;
}

unsigned short LAS_io::size_extra_bytes(unsigned short Point_Format_ID, unsigned short Point_Record_Length) // how many bytes that are extra from header
{
	return (Point_Record_Length - LAS_Point_Size[Point_Format_ID]);
}
	
unsigned short LAS_io::size_extra_bytes(LAS_Extra_Bytes_Info *eb_info, unsigned short n_eb)  // how many bytes that are extra
{
	unsigned short len_eb = 0;
	for (unsigned int i_eb = 0; i_eb < n_eb; i_eb++)
	{
		len_eb += eb_info[i_eb].size_value * eb_info[i_eb].n_values;
	}

	return len_eb;
}

long long LAS_io::idx_vlr_extra_bytes() // vlr index for extra bytes
{
	long long idx_vlr = -1;
#pragma omp parallel for schedule(static)
	for (long long idx = 0; idx < las_header.Number_of_Variable_Length_Records; idx++)
	{
		if (las_vlr[idx].vlr_header.Record_ID == VLR_Record_ID_Extra_Bytes)
		{
			idx_vlr = idx;
			break;
		}
	}

	return idx_vlr;
}

unsigned short LAS_io::n_extra_bytes()  // number of extra bytes/attributes
{
	long long idx_vlr = idx_vlr_extra_bytes();
	unsigned short n_eb = 0;

	if (idx_vlr > -1)
	{
		n_eb = las_vlr[idx_vlr].vlr_header.Record_Length_After_Header / LASIO_VLR_Extra_Bytes_Size;
	} 
	return n_eb;
}

long long LAS_io::Get_All_Extra_Bytes_Info(LAS_Extra_Bytes_Info* &eb_info, unsigned short &n_eb) // get an extra bytes info
{
	long long idx_vlr = idx_vlr_extra_bytes();	
	n_eb = n_extra_bytes();

	if (n_eb > 0)
	{
		eb_info = new LAS_Extra_Bytes_Info[n_eb];
		unsigned short pos = 0;
		for (unsigned short i_eb = 0; i_eb < n_eb; i_eb++)
		{
			Get_One_Extra_Bytes_Info(eb_info[i_eb], i_eb);
			eb_info[i_eb].byte_pos = pos;
			pos += eb_info[i_eb].size_value * eb_info[i_eb].n_values;
		}
	}

	return idx_vlr;
}

long long LAS_io::Get_One_Extra_Bytes_Info(LAS_Extra_Bytes_Info &eb_info, unsigned short idx_eb) // get an extra bytes info
{
	long long idx_vlr = idx_vlr_extra_bytes();
	unsigned long long pos = idx_eb * LASIO_VLR_Extra_Bytes_Size;
	Bytes_to_Extra_Bytes_Info(las_vlr[idx_vlr].vlr_buffer + pos, eb_info);
	Set_Extra_Bytes_Type_Info(eb_info);
	return idx_vlr;
}

int LAS_io::Get_Extra_Bytes_Info(LAS_Extra_Bytes_Info &eb_info, char* name, LAS_Extra_Bytes_Data_Type type_id)
{
	int idx_eb = -1;
	LAS_Extra_Bytes_Info *eb;
	unsigned short n_eb = 0;
	Get_All_Extra_Bytes_Info(eb, n_eb);
	
	if (n_eb > 0)
	{
		for (unsigned short i_eb = 0; i_eb < n_eb; i_eb++)
		{
			if (EZPC::c_str_cmpi(eb[i_eb].name, name) == 0 && eb[i_eb].data_type == type_id)
			{
				idx_eb = i_eb;
				eb_info = eb[i_eb];
			}
		}
	}
	if (n_eb > 0)
	{
		delete[] eb;
		eb = NULL;
	}

	return idx_eb;
}

void LAS_io::Bytes_to_Extra_Bytes_Info(LAS_BYTE* bytes, LAS_Extra_Bytes_Info &eb_info)
{
	unsigned long long pos = 0;
	Bytes2Data(bytes, eb_info.reserved, sizeof(eb_info.reserved), 1, pos);
	Bytes2Data(bytes, &eb_info.data_type, sizeof(eb_info.data_type), 1, pos);

	LAS_BYTE options;
	Bytes2Data(bytes, &options, sizeof(options), 1, pos);

	eb_info.options_no_data_bit = (options >> 0 & 1) > 0 ? true : false;
	eb_info.options_min_bit = (options >> 1 & 1) > 0 ? true : false;
	eb_info.options_max_bit = (options >> 2 & 1) > 0 ? true : false;
	eb_info.options_scale_bit = (options >> 3 & 1) > 0 ? true : false;
	eb_info.options_offset_bit = (options >> 4 & 1) > 0 ? true : false;

	Bytes2Data(bytes, eb_info.name, sizeof(eb_info.name), 1, pos);
	Bytes2Data(bytes, eb_info.unused, sizeof(eb_info.unused), 1, pos);
	Bytes2Data(bytes, eb_info.no_data, sizeof(eb_info.no_data), 1, pos);
	Bytes2Data(bytes, eb_info.min, sizeof(eb_info.min), 1, pos);
	Bytes2Data(bytes, eb_info.max, sizeof(eb_info.max), 1, pos);
	Bytes2Data(bytes, eb_info.scale, sizeof(eb_info.scale), 1, pos);
	Bytes2Data(bytes, eb_info.offset, sizeof(eb_info.offset), 1, pos);
	Bytes2Data(bytes, eb_info.description, sizeof(eb_info.description), 1, pos);
}

void LAS_io::Bytes_from_Extra_Bytes_Info(LAS_BYTE* bytes, LAS_Extra_Bytes_Info eb_info)
{
	unsigned long long pos = 0;
	Data2Bytes(bytes, eb_info.reserved, sizeof(eb_info.reserved), 1, pos);
	Data2Bytes(bytes, &eb_info.data_type, sizeof(eb_info.data_type), 1, pos);

	LAS_BYTE options = 0;
	options |= (eb_info.options_no_data_bit ? 1 : 0) << 0;
	options |= (eb_info.options_min_bit ? 1 : 0) << 1;
	options |= (eb_info.options_max_bit ? 1 : 0) << 2;
	options |= (eb_info.options_scale_bit ? 1 : 0) << 3;
	options |= (eb_info.options_offset_bit ? 1 : 0) << 4;

	Data2Bytes(bytes, &options, sizeof(options), 1, pos);
	Data2Bytes(bytes, eb_info.name, sizeof(eb_info.name), 1, pos);
	Data2Bytes(bytes, eb_info.unused, sizeof(eb_info.unused), 1, pos);
	Data2Bytes(bytes, eb_info.no_data, sizeof(eb_info.no_data), 1, pos);
	Data2Bytes(bytes, eb_info.min, sizeof(eb_info.min), 1, pos);
	Data2Bytes(bytes, eb_info.max, sizeof(eb_info.max), 1, pos);
	Data2Bytes(bytes, eb_info.scale, sizeof(eb_info.scale), 1, pos);
	Data2Bytes(bytes, eb_info.offset, sizeof(eb_info.offset), 1, pos);
	Data2Bytes(bytes, eb_info.description, sizeof(eb_info.description), 1, pos);
}

unsigned short LAS_io::Add_Extra_Bytes_VLR(LAS_Extra_Bytes_Info* eb_info, unsigned short n_eb) // add fields to vlr and allocate extra bytes for points, return the starting idndex of the new extra bytes 
{
	long long idx_vlr = idx_vlr_extra_bytes();

	if (idx_vlr < 0)
	{
		idx_vlr = las_header.Number_of_Variable_Length_Records;
		LAS_VLR eb_vlr = create_a_vlr(0);
		char usr_id[] = "LASF_Spec";
		char desc[] = "Extra_Bytes";
		Set_VLR_User_ID(eb_vlr, usr_id);
		Set_VLR_Record_ID(eb_vlr, VLR_Record_ID_Extra_Bytes);
		Set_VLR_Description(eb_vlr, desc);
		Add_New_VLR(eb_vlr);
	}

	LAS_BYTE* buffer_cpy = las_vlr[idx_vlr].vlr_buffer;

	las_vlr[idx_vlr].vlr_buffer = new LAS_BYTE[las_vlr[idx_vlr].vlr_header.Record_Length_After_Header + LASIO_VLR_Extra_Bytes_Size * n_eb];
	if (las_vlr[idx_vlr].vlr_header.Record_Length_After_Header > 0)
	{
		memcpy(las_vlr[idx_vlr].vlr_buffer, buffer_cpy, las_vlr[idx_vlr].vlr_header.Record_Length_After_Header);
		delete[] buffer_cpy;
	}

	unsigned long long offset = las_vlr[idx_vlr].vlr_header.Record_Length_After_Header;
	unsigned short len_eb = size_extra_bytes(true);

	LAS_BYTE* bytes_tmp = new LAS_BYTE[LASIO_VLR_Extra_Bytes_Size];
	for (unsigned int i_eb = 0; i_eb < n_eb; i_eb++)
	{
		Bytes_from_Extra_Bytes_Info(bytes_tmp, eb_info[i_eb]);
		Data2Bytes(las_vlr[idx_vlr].vlr_buffer, bytes_tmp, LASIO_VLR_Extra_Bytes_Size, 1, offset);
		eb_info[i_eb].byte_pos = len_eb;
		len_eb += eb_info[i_eb].size_value * eb_info[i_eb].n_values;
	}
	delete[] bytes_tmp;

	LAS_BYTE** eb_bytes = new LAS_BYTE*[las_header.Number_of_point_records];
#pragma omp parallel for schedule (static)
	for (long long pid = 0; pid < (long long)las_header.Number_of_point_records; pid++)
	{
		eb_bytes[pid] = new LAS_BYTE[len_eb];
	}

#pragma omp parallel for schedule (static)
	for (long long pid = 0; pid < (long long)las_header.Number_of_point_records; pid++)
	{
		if (size_extra_bytes(true) > 0 && las_points[pid].Extra_Bytes)
		{
			memcpy(eb_bytes[pid], las_points[pid].Extra_Bytes, size_extra_bytes(true));
			delete[] las_points[pid].Extra_Bytes;
		}
		las_points[pid].Extra_Bytes = eb_bytes[pid];
	}
	delete[] eb_bytes;

	las_vlr[idx_vlr].vlr_header.Record_Length_After_Header += n_eb * LASIO_VLR_Extra_Bytes_Size;
	las_header.Point_Data_Record_Length = LAS_Point_Size[las_header.Point_Data_Format_ID] + len_eb;

	return (unsigned short)(n_extra_bytes() - n_eb);
}

void LAS_io::Initialize_Point_Extra_Bytes() // allocate extra bytes from header
{
	if (size_extra_bytes(true) > 0)
	{
#pragma omp parallel for schedule (static)
		for (long long pid = 0; pid < (long long)las_header.Number_of_point_records; pid++)
		{
			if (las_points[pid].Extra_Bytes)
			{
				delete[] las_points[pid].Extra_Bytes;
			}
			las_points[pid].Extra_Bytes = new LAS_BYTE[size_extra_bytes(true)];
		}
	}
}

void LAS_io::Clear_Extra_Bytes() // clear all extra bytes in both vlrs and points
{
	long long idx_eb = idx_vlr_extra_bytes();

	if (idx_eb > -1)
	{
		LAS_VLR *vlr_cpy = las_vlr;
		las_header.Number_of_Variable_Length_Records--;
		las_vlr = new LAS_VLR[las_header.Number_of_Variable_Length_Records];

#pragma omp parallel for schedule (static)
		for (long long vid = 0; vid < las_header.Number_of_Variable_Length_Records; vid++)
		{
			unsigned int vid_cpy = vid < idx_eb ? (unsigned int)vid : (unsigned int)(vid + 1);

			las_vlr[vid].vlr_header = vlr_cpy[vid_cpy].vlr_header;
			las_vlr[vid].vlr_buffer = vlr_cpy[vid_cpy].vlr_buffer;
		}

		delete[] vlr_cpy;
	}

	if(size_extra_bytes(true) > 0)
	{
#pragma omp parallel for schedule (static)
		for (long long pid = 0; pid < (long long)las_header.Number_of_point_records; pid++)
		{
			if (las_points[pid].Extra_Bytes)
			{
				delete[] las_points[pid].Extra_Bytes;
				las_points[pid].Extra_Bytes = NULL;
			}
		}

		las_header.Point_Data_Record_Length = LAS_Point_Size[las_header.Point_Data_Format_ID];
	}
}

void LAS_io::Remove_Extra_Bytes(char* name, LAS_Extra_Bytes_Data_Type type_id)
{
	LAS_Extra_Bytes_Info eb_info;
	int idx_eb = Get_Extra_Bytes_Info(eb_info, name, type_id);

	if (idx_eb > -1)
	{
		if (n_extra_bytes() == 1)
		{
			Clear_Extra_Bytes();
		}
		else
		{
			long long idx_vlr_eb = idx_vlr_extra_bytes();
			if (idx_vlr_eb > -1)
			{
				LAS_BYTE* nu_buff = NULL;
				if (las_vlr[idx_vlr_eb].vlr_header.Record_Length_After_Header < LASIO_VLR_Extra_Bytes_Size)
				{
					nu_buff = new LAS_BYTE[las_vlr[idx_vlr_eb].vlr_header.Record_Length_After_Header - LASIO_VLR_Extra_Bytes_Size];
					memcpy(nu_buff, las_vlr[idx_vlr_eb].vlr_buffer, idx_eb * LASIO_VLR_Extra_Bytes_Size);
					memcpy(nu_buff + idx_eb * LASIO_VLR_Extra_Bytes_Size,
						las_vlr[idx_vlr_eb].vlr_buffer + (idx_eb + 1) * LASIO_VLR_Extra_Bytes_Size,
						las_vlr[idx_vlr_eb].vlr_header.Record_Length_After_Header - (idx_eb + 1) * LASIO_VLR_Extra_Bytes_Size);

					las_vlr[idx_vlr_eb].vlr_header.Record_Length_After_Header -= LASIO_VLR_Extra_Bytes_Size;
				}

				delete[] las_vlr[idx_vlr_eb].vlr_buffer;
				las_vlr[idx_vlr_eb].vlr_buffer = nu_buff;
			}

#pragma omp parallel for schedule (static)
			for (long long pid = 0; pid < (long long)las_header.Number_of_point_records; pid++)
			{
				if (las_points[pid].Extra_Bytes)
				{
					LAS_BYTE* nu_eb = NULL;
					if (size_extra_bytes(true) > eb_info.size_value)
					{
						nu_eb = new LAS_BYTE[size_extra_bytes(true) - eb_info.size_value * eb_info.n_values];
						memcpy(nu_eb, las_points[pid].Extra_Bytes, eb_info.byte_pos);
						memcpy(nu_eb + eb_info.byte_pos,
							las_points[pid].Extra_Bytes + eb_info.byte_pos + eb_info.size_value * eb_info.n_values,
							size_extra_bytes(true) - eb_info.size_value * eb_info.n_values - eb_info.byte_pos);

						las_header.Point_Data_Record_Length -= eb_info.size_value * eb_info.n_values;
					}

					delete[] las_points[pid].Extra_Bytes;
					las_points[pid].Extra_Bytes = nu_eb;
				}
			}
		}
	}
}

LAS_Extra_Bytes_Info LAS_io::create_extra_bytes_info(char* name, char* description, LAS_Extra_Bytes_Data_Type data_type) // create an extra bytes info structure
{
	LAS_Extra_Bytes_Info eb_info;
	eb_info.data_type = data_type;
	Set_Extra_Bytes_Type_Info(eb_info);
	Set_Extra_Bytes_Name(eb_info, name);
	Set_Extra_Bytes_Description(eb_info, description);
	return eb_info;
}

LAS_Extra_Bytes_Info LAS_io::create_extra_bytes_info(LAS_Extra_Bytes_Data_Type data_type) // create an extra bytes info structure
{
	LAS_Extra_Bytes_Info eb_info;
	eb_info.data_type = data_type;
	Set_Extra_Bytes_Type_Info(eb_info);
	return eb_info;
}

void LAS_io::Set_Extra_Bytes_Type_Info(LAS_Extra_Bytes_Info &eb_info)
{
	eb_info.size_value = LAS_Extra_Bytes_Type_Info[eb_info.data_type][0];
	eb_info.n_values = LAS_Extra_Bytes_Type_Info[eb_info.data_type][1];
	eb_info.floating_or_integer = LAS_Extra_Bytes_Type_Info[eb_info.data_type][2] > 0 ? true : false;
	eb_info.signed_or_unsigned = LAS_Extra_Bytes_Type_Info[eb_info.data_type][3] > 0 ? true : false;
}

void LAS_io::Set_Extra_Bytes_Options_Invalid(LAS_Extra_Bytes_Info &eb_info)
{
	Set_Extra_Bytes_No_Data_Invalid(eb_info);
	Set_Extra_Bytes_Min_Invalid(eb_info);
	Set_Extra_Bytes_Max_Invalid(eb_info);
	Set_Extra_Bytes_Scale_Invalid(eb_info);
	Set_Extra_Bytes_Offset_Invalid(eb_info);
}

void LAS_io::Set_Extra_Bytes_No_Data_Invalid(LAS_Extra_Bytes_Info &eb_info)
{
	eb_info.options_no_data_bit = false;
	eb_info.no_data[0] = 0;
	eb_info.no_data[1] = 0;
	eb_info.no_data[2] = 0;
}

void LAS_io::Set_Extra_Bytes_Min_Invalid(LAS_Extra_Bytes_Info &eb_info)
{
	eb_info.options_min_bit = false;
	eb_info.min[0] = 0;
	eb_info.min[1] = 0;
	eb_info.min[2] = 0;
}

void LAS_io::Set_Extra_Bytes_Max_Invalid(LAS_Extra_Bytes_Info &eb_info)
{
	eb_info.options_max_bit = false;
	eb_info.max[0] = 0;
	eb_info.max[1] = 0;
	eb_info.max[2] = 0;
}

void LAS_io::Set_Extra_Bytes_Scale_Invalid(LAS_Extra_Bytes_Info &eb_info)
{
	eb_info.options_scale_bit = false;
	eb_info.scale[0] = 0;
	eb_info.scale[1] = 0;
	eb_info.scale[2] = 0;
}

void LAS_io::Set_Extra_Bytes_Offset_Invalid(LAS_Extra_Bytes_Info &eb_info)
{
	eb_info.options_offset_bit = false;
	eb_info.offset[0] = 0;
	eb_info.offset[1] = 0;
	eb_info.offset[2] = 0;
}

void LAS_io::Set_Extra_Bytes_Name(LAS_Extra_Bytes_Info &eb_info, char* name)
{
	size_t len = strlen(name);
	memcpy(eb_info.name, name, min(len, (size_t)32));
	if (len < 32)
	{
		eb_info.name[len] = '\0';
	}
}

void LAS_io::Set_Extra_Bytes_No_Data(LAS_Extra_Bytes_Info &eb_info, double no_data)
{
	if (!eb_info.floating_or_integer)
	{
		printf("\n\nSet extra bytes no data: input cannot be a floating point");
		getchar();
	}
	else
	{
		Set_Extra_Bytes_No_Data(eb_info, &no_data);
	}
}

void LAS_io::Set_Extra_Bytes_No_Data(LAS_Extra_Bytes_Info &eb_info, long long no_data)
{
	if (eb_info.floating_or_integer || !eb_info.signed_or_unsigned)
	{
		printf("\n\nSet extra bytes no data: input cannot be a long long");
		getchar();
	}
	else
	{
		Set_Extra_Bytes_No_Data(eb_info, &no_data);
	}
}

void LAS_io::Set_Extra_Bytes_No_Data(LAS_Extra_Bytes_Info &eb_info, unsigned long long no_data)
{
	if (eb_info.floating_or_integer || eb_info.signed_or_unsigned)
	{
		printf("\n\nSet extra bytes no data: input cannot be a unsigned long long");
		getchar();
	}
	else
	{
		Set_Extra_Bytes_No_Data(eb_info, &no_data);
	}
}

void LAS_io::Set_Extra_Bytes_No_Data(LAS_Extra_Bytes_Info &eb_info, void* no_data_x, void* no_data_y, void* no_data_z)
{
	Set_Extra_Bytes_No_Data_Invalid(eb_info);
	eb_info.options_no_data_bit = true;

	Upcasting_to_Byte_8(eb_info, eb_info.no_data[0], no_data_x);
	if (eb_info.n_values > 1)
	{
		Upcasting_to_Byte_8(eb_info, eb_info.no_data[1], no_data_y);
		if (eb_info.n_values > 2)
		{
			Upcasting_to_Byte_8(eb_info, eb_info.no_data[2], no_data_z);
		}
	}
}

void LAS_io::Set_Extra_Bytes_No_Data(LAS_Extra_Bytes_Info &eb_info, void* no_data_x, void* no_data_y)
{
	void* tmp = NULL;
	Set_Extra_Bytes_No_Data(eb_info, no_data_x, no_data_y, tmp);
}

void LAS_io::Set_Extra_Bytes_No_Data(LAS_Extra_Bytes_Info &eb_info, void* no_data)
{
	void* tmp = NULL;
	Set_Extra_Bytes_No_Data(eb_info, no_data, tmp, tmp);
}

void LAS_io::Set_Extra_Bytes_Min(LAS_Extra_Bytes_Info &eb_info, double min)
{
	if (!eb_info.floating_or_integer)
	{
		printf("\n\nSet extra bytes min: input cannot be a floating point");
		getchar();
	}
	else
	{
		Set_Extra_Bytes_Min(eb_info, &min);
	}
}

void LAS_io::Set_Extra_Bytes_Min(LAS_Extra_Bytes_Info &eb_info, long long min)
{
	if (eb_info.floating_or_integer || !eb_info.signed_or_unsigned)
	{
		printf("\n\nSet extra bytes min: input cannot be a long long");
		getchar();
	}
	else
	{
		Set_Extra_Bytes_Min(eb_info, &min);
	}
}

void LAS_io::Set_Extra_Bytes_Min(LAS_Extra_Bytes_Info &eb_info, unsigned long long min)
{
	if (eb_info.floating_or_integer || eb_info.signed_or_unsigned)
	{
		printf("\n\nSet extra bytes min: input cannot be a unsigned long long");
		getchar();
	}
	else
	{
		Set_Extra_Bytes_Min(eb_info, &min);
	}
}

void LAS_io::Set_Extra_Bytes_Min(LAS_Extra_Bytes_Info &eb_info, void* min_x, void* min_y, void* min_z)
{
	Set_Extra_Bytes_Min_Invalid(eb_info);
	eb_info.options_min_bit = true;
	Upcasting_to_Byte_8(eb_info, eb_info.min[0], min_x);
	if (eb_info.n_values > 1)
	{
		Upcasting_to_Byte_8(eb_info, eb_info.min[1], min_y);
		if (eb_info.n_values > 2)
		{
			Upcasting_to_Byte_8(eb_info, eb_info.min[2], min_z);
		}
	}
}

void LAS_io::Set_Extra_Bytes_Min(LAS_Extra_Bytes_Info &eb_info, void* min_x, void* min_y)
{
	void* tmp = NULL;
	Set_Extra_Bytes_Min(eb_info, min_x, min_y, tmp);
}

void LAS_io::Set_Extra_Bytes_Min(LAS_Extra_Bytes_Info &eb_info, void* min)
{
	void* tmp = NULL;
	Set_Extra_Bytes_Min(eb_info, min, tmp, tmp);
}

void LAS_io::Set_Extra_Bytes_Max(LAS_Extra_Bytes_Info &eb_info, double max)
{
	if (!eb_info.floating_or_integer)
	{
		printf("\n\nSet extra bytes max: input cannot be a floating point");
		getchar();
	}
	else
	{
		Set_Extra_Bytes_Max(eb_info, &max);
	}
}

void LAS_io::Set_Extra_Bytes_Max(LAS_Extra_Bytes_Info &eb_info, long long max)
{
	if (eb_info.floating_or_integer || !eb_info.signed_or_unsigned)
	{
		printf("\n\nSet extra bytes max: input cannot be a long long");
		getchar();
	}
	else
	{
		Set_Extra_Bytes_Max(eb_info, &max);
	}
}

void LAS_io::Set_Extra_Bytes_Max(LAS_Extra_Bytes_Info &eb_info, unsigned long long max)
{
	if (eb_info.floating_or_integer || eb_info.signed_or_unsigned)
	{
		printf("\n\nSet extra bytes max: input cannot be a unsigned long long");
		getchar();
	}
	else
	{
		Set_Extra_Bytes_Max(eb_info, &max);
	}
}

void LAS_io::Set_Extra_Bytes_Max(LAS_Extra_Bytes_Info &eb_info, void* max_x, void* max_y, void* max_z)
{
	Set_Extra_Bytes_Max_Invalid(eb_info);
	eb_info.options_max_bit = true;
	Upcasting_to_Byte_8(eb_info, eb_info.max[0], max_x);
	if (eb_info.n_values > 1)
	{
		Upcasting_to_Byte_8(eb_info, eb_info.max[1], max_y);
		if (eb_info.n_values > 2)
		{
			Upcasting_to_Byte_8(eb_info, eb_info.max[2], max_z);
		}
	}
}

void LAS_io::Set_Extra_Bytes_Max(LAS_Extra_Bytes_Info &eb_info, void* max_x, void* max_y)
{
	void* tmp = NULL;
	Set_Extra_Bytes_Max(eb_info, max_x, max_y, tmp);
}

void LAS_io::Set_Extra_Bytes_Max(LAS_Extra_Bytes_Info &eb_info, void* max)
{
	void* tmp = NULL;
	Set_Extra_Bytes_Max(eb_info, max, tmp, tmp);
}

void LAS_io::Set_Extra_Bytes_Scale(LAS_Extra_Bytes_Info &eb_info, double scale_x, double scale_y, double scale_z)
{
	Set_Extra_Bytes_Scale_Invalid(eb_info);
	eb_info.options_scale_bit = true;
	eb_info.scale[0] = scale_x;
	if (eb_info.n_values > 1)
	{
		eb_info.scale[1] = scale_y;
		if (eb_info.n_values > 2)
		{
			eb_info.scale[2] = scale_z;
		}
	}
}

void LAS_io::Set_Extra_Bytes_Scale(LAS_Extra_Bytes_Info &eb_info, double scale_x, double scale_y)
{
	Set_Extra_Bytes_Scale(eb_info, scale_x, scale_y, 0);
}

void LAS_io::Set_Extra_Bytes_Scale(LAS_Extra_Bytes_Info &eb_info, double scale)
{
	Set_Extra_Bytes_Scale(eb_info, scale, 0, 0);
}

void LAS_io::Set_Extra_Bytes_Offset(LAS_Extra_Bytes_Info &eb_info, double offset_x, double offset_y, double offset_z)
{
	Set_Extra_Bytes_Offset_Invalid(eb_info);
	eb_info.options_offset_bit = true;
	eb_info.offset[0] = offset_x;
	if (eb_info.n_values > 1)
	{
		eb_info.offset[1] = offset_y;
		if (eb_info.n_values > 2)
		{
			eb_info.offset[2] = offset_z;
		}
	}
}

void LAS_io::Set_Extra_Bytes_Offset(LAS_Extra_Bytes_Info &eb_info, double offset_x, double offset_y)
{
	Set_Extra_Bytes_Offset(eb_info, offset_x, offset_y, 0);
}

void LAS_io::Set_Extra_Bytes_Offset(LAS_Extra_Bytes_Info &eb_info, double offset)
{
	Set_Extra_Bytes_Offset(eb_info, offset, 0, 0);
}

void LAS_io::Set_Extra_Bytes_Description(LAS_Extra_Bytes_Info &eb_info, char* description)
{
	size_t len = strlen(description);
	memcpy(eb_info.description, description, min(len, (size_t)32));
	if (len < 32)
	{
		eb_info.description[len] = '\0';
	}
}

void LAS_io::Set_Extra_Bytes_Value_No_Data(unsigned long long pid, LAS_Extra_Bytes_Info eb_info)
{
	if (eb_info.floating_or_integer)
	{	
		Set_Extra_Bytes_Value(pid, eb_info, double_from_byte_8(eb_info.no_data[0]));
	}
	else
	{
		if (eb_info.signed_or_unsigned)
		{
			Set_Extra_Bytes_Value(pid, eb_info, int64_from_byte_8(eb_info.no_data[0]));
		}
		else
		{
			Set_Extra_Bytes_Value(pid, eb_info, uint64_from_byte_8(eb_info.no_data[0]));
		}
	}
}

void LAS_io::Set_Extra_Bytes_Value(unsigned long long pid, LAS_Extra_Bytes_Info eb_info, double in_value) // add a field value
{
	if (!eb_info.floating_or_integer)
	{
		printf("\n\nSet extra bytes value: input cannot be a floating point. data type = %i", eb_info.data_type);
		getchar();
	}
	else
	{
		Set_Extra_Bytes_Value(pid, eb_info, &in_value);
	}
}

void LAS_io::Set_Extra_Bytes_Value(unsigned long long pid, LAS_Extra_Bytes_Info eb_info, long long in_value) // add a field value
{
	if (eb_info.floating_or_integer || !eb_info.signed_or_unsigned)
	{
		printf("\n\nSet extra bytes value: input cannot be a long long. data type = %i", eb_info.data_type);		
		getchar();
	}
	else
	{
		Set_Extra_Bytes_Value(pid, eb_info, &in_value);
	}
}

void LAS_io::Set_Extra_Bytes_Value(unsigned long long pid, LAS_Extra_Bytes_Info eb_info, unsigned long long in_value) // add a field value
{
	if (eb_info.floating_or_integer || eb_info.signed_or_unsigned)
	{
		printf("\n\nSet extra bytes value: input cannot be a unsigned long long. data type = %i", eb_info.data_type);	
		getchar();
	}
	else
	{
		Set_Extra_Bytes_Value(pid, eb_info, &in_value);
	}
}

void LAS_io::Set_Extra_Bytes_Value(unsigned long long pid, LAS_Extra_Bytes_Info eb_info, void* in_value) // add a field value
{
	void* tmp = NULL;
	Set_Extra_Bytes_Value(pid, eb_info, in_value, tmp, tmp);
}

void LAS_io::Set_Extra_Bytes_Value(unsigned long long pid, LAS_Extra_Bytes_Info eb_info, void* in_value_x, void* in_value_y)// add a field value
{
	void* tmp = NULL;
	Set_Extra_Bytes_Value(pid, eb_info, in_value_x, in_value_y, tmp);
}

void LAS_io::Set_Extra_Bytes_Value(unsigned long long pid, LAS_Extra_Bytes_Info eb_info, void* in_value_x, void* in_value_y, void* in_value_z)// add a field value
{
	unsigned long long pos = (unsigned long long)eb_info.byte_pos;

	LAS_BYTE byte_x[8];
	LAS_BYTE byte_y[8];
	LAS_BYTE byte_z[8];

	Downcasting_to_Bytes(eb_info, byte_x, in_value_x);
	Data2Bytes(las_points[pid].Extra_Bytes, byte_x, eb_info.size_value, 1, pos);

	if (eb_info.n_values > 1)
	{
		Downcasting_to_Bytes(eb_info, byte_y, in_value_y);
		Data2Bytes(las_points[pid].Extra_Bytes, byte_y, eb_info.size_value, 1, pos);
		if (eb_info.n_values > 2)
		{
			Downcasting_to_Bytes(eb_info, byte_z, in_value_z);
			Data2Bytes(las_points[pid].Extra_Bytes, byte_z, eb_info.size_value, 1, pos);
		}
	}
}

bool LAS_io::is_extra_bytes_value_valid(unsigned long long pid, LAS_Extra_Bytes_Info eb_info)
{
	LAS_BYTE_8 val_x, val_y, val_z; 
	Get_Extra_Bytes_Value_Byte_8(pid, eb_info, val_x, val_y, val_z);
	return (val_x == eb_info.no_data[0] && val_y == eb_info.no_data[1] && val_z == eb_info.no_data[2]);
}

void LAS_io::Get_Extra_Bytes_Value(unsigned long long pid, LAS_Extra_Bytes_Info eb_info, double &val)
{
	if (!eb_info.floating_or_integer)
	{
		printf("\n\nGet extra bytes value: input cannot be a floating point");
		getchar();
	}
	else
	{
		LAS_BYTE_8 tmp;
		Get_Extra_Bytes_Value_Byte_8(pid, eb_info, tmp);
		val = double_from_byte_8(tmp);
	}
}

void LAS_io::Get_Extra_Bytes_Value(unsigned long long pid, LAS_Extra_Bytes_Info eb_info, long long &val)
{
	if (eb_info.floating_or_integer || !eb_info.signed_or_unsigned)
	{
		printf("\n\nGet extra bytes value: input cannot be a long long");
		getchar();
	}
	else
	{
		LAS_BYTE_8 tmp;
		Get_Extra_Bytes_Value_Byte_8(pid, eb_info, tmp);
		val = int64_from_byte_8(tmp);
	}
}

void LAS_io::Get_Extra_Bytes_Value(unsigned long long pid, LAS_Extra_Bytes_Info eb_info, unsigned long long &val)
{
	if (eb_info.floating_or_integer || eb_info.signed_or_unsigned)
	{
		printf("\n\nGet extra bytes value: input cannot be a unsigned long long");
		getchar();
	}
	else
	{
		LAS_BYTE_8 tmp;
		Get_Extra_Bytes_Value_Byte_8(pid, eb_info, tmp);
		val = uint64_from_byte_8(tmp);
	}
}

void LAS_io::Get_Extra_Bytes_Value_Byte_8(unsigned long long pid, LAS_Extra_Bytes_Info eb_info, LAS_BYTE_8 &val)
{
	LAS_BYTE_8 tmp;
	Get_Extra_Bytes_Value_Byte_8(pid, eb_info, val, tmp, tmp);
}

void LAS_io::Get_Extra_Bytes_Value_Byte_8(unsigned long long pid, LAS_Extra_Bytes_Info eb_info, LAS_BYTE_8 &val_x, LAS_BYTE_8 &val_y)
{
	LAS_BYTE_8 tmp;
	Get_Extra_Bytes_Value_Byte_8(pid, eb_info, val_x, val_y, tmp);
}

void LAS_io::Get_Extra_Bytes_Value_Byte_8(unsigned long long pid, LAS_Extra_Bytes_Info eb_info, LAS_BYTE_8 &val_x, LAS_BYTE_8 &val_y, LAS_BYTE_8 &val_z)
{
	LAS_BYTE *bytes = new LAS_BYTE[eb_info.size_value];
	unsigned long long pos = (unsigned long long)eb_info.byte_pos;
	LAS_BYTE_8 val[3];
	val[0] = eb_info.no_data[0];
	val[1] = eb_info.no_data[1];
	val[2] = eb_info.no_data[2];

	for (int i_val = 0; i_val < eb_info.n_values; i_val++)
	{
		Bytes2Data(las_points[pid].Extra_Bytes, bytes, eb_info.size_value, 1, pos);
		if (eb_info.floating_or_integer)
		{
			double result;
			switch (eb_info.data_type)
			{
				case EB_Type_Float: case EB_Type_Float_2: case EB_Type_Float_3:
				{
					float tmp;
					memcpy(&tmp, bytes, eb_info.size_value);
					result = tmp;
					break;
				}
				case EB_Type_Double: case EB_Type_Double_2: case EB_Type_Double_3:
				{
					double tmp;
					memcpy(&tmp, bytes, eb_info.size_value);
					result = tmp;
					break;
				}
			}
			memcpy(&val[i_val], &result, sizeof(LAS_BYTE_8));
		}
		else
		{
			if (eb_info.signed_or_unsigned)
			{
				long long result;
				switch (eb_info.data_type)
				{
					case EB_Type_Char: case EB_Type_Char_2: case EB_Type_Char_3:
					{
						char tmp;
						memcpy(&tmp, bytes, eb_info.size_value);
						result = tmp;
						break;
					}
					case EB_Type_Short: case EB_Type_Short_2: case EB_Type_Short_3:
					{
						short tmp;
						memcpy(&tmp, bytes, eb_info.size_value);
						result = tmp;
						break;
					}
					case EB_Type_Long: case EB_Type_Long_2: case EB_Type_Long_3:
					{
						long tmp;
						memcpy(&tmp, bytes, eb_info.size_value);
						result = tmp;
						break;
					}
					case EB_Type_Long_Long: case EB_Type_Long_Long_2: case EB_Type_Long_Long_3:
					{
						long long tmp;
						memcpy(&tmp, bytes, eb_info.size_value);
						result = tmp;
						break;
					}
				}
				memcpy(&val[i_val], &result, sizeof(LAS_BYTE_8));
			}
			else
			{
				unsigned long long result;
				switch (eb_info.data_type)
				{
					case EB_Type_Unsigned_Char: case EB_Type_Unsigned_Char_2: case EB_Type_Unsigned_Char_3:
					{
						unsigned char tmp;
						memcpy(&tmp, bytes, eb_info.size_value);
						result = tmp;
						break;
					}
					case EB_Type_Unsigned_Short: case EB_Type_Unsigned_Short_2:	case EB_Type_Unsigned_Short_3:
					{
						unsigned short tmp;
						memcpy(&tmp, bytes, eb_info.size_value);
						result = tmp;
						break;
					}
					case EB_Type_Unsigned_Long: case EB_Type_Unsigned_Long_2: case EB_Type_Unsigned_Long_3:
					{
						unsigned int tmp;
						memcpy(&tmp, bytes, eb_info.size_value);
						result = tmp;
						break;
					}
					case EB_Type_Unsigned_Long_Long: case EB_Type_Unsigned_Long_Long_2: case EB_Type_Unsigned_Long_Long_3:
					{
						unsigned long long tmp;
						memcpy(&tmp, bytes, eb_info.size_value);
						result = tmp;
						break;
					}
				}
				memcpy(&val[i_val], &result, sizeof(LAS_BYTE_8));
			}
		}
	}

	val_x = val[0];
	val_y = val[1];
	val_z = val[2];
	delete[] bytes;
}

void LAS_io::Upcasting_to_Byte_8(LAS_Extra_Bytes_Info eb_info, LAS_BYTE_8 &bytes, void* value)
{
	if (eb_info.floating_or_integer)
	{
		double tmp = *(double*)value;
		memcpy(&bytes, &tmp, sizeof(LAS_BYTE_8));
	}
	else
	{
		if (eb_info.signed_or_unsigned)
		{
			long long tmp = *(long long*)value;

			printf("\ntmp %lli", tmp);
			getchar();

			memcpy(&bytes, &tmp, sizeof(LAS_BYTE_8));
		}
		else
		{
			unsigned long long tmp = *(unsigned long long*)value;
			memcpy(&bytes, &tmp, sizeof(LAS_BYTE_8));
		}
	}
}

void LAS_io::Downcasting_to_Bytes(LAS_Extra_Bytes_Info eb_info, LAS_BYTE* bytes, void* value)
{
	if (eb_info.data_type == EB_Type_Undocumented)
	{
		printf("\n\nUndocumented Extra Bytes Type Unsupported Yet!");
		getchar();
	}

	if (eb_info.floating_or_integer)
	{
		switch (eb_info.data_type)
		{
			case EB_Type_Float: case EB_Type_Float_2: case EB_Type_Float_3:
			{
				float tmp = *(float*)value;
				memcpy(bytes, &tmp, eb_info.size_value);
				break;
			}
			case EB_Type_Double: case EB_Type_Double_2: case EB_Type_Double_3:
			{
				double tmp = *(double*)value;
				memcpy(bytes, &tmp, eb_info.size_value);
				break;
			}
		}
	}
	else
	{
		if (eb_info.signed_or_unsigned)
		{
			switch (eb_info.data_type)
			{
				case EB_Type_Char: case EB_Type_Char_2: case EB_Type_Char_3:
				{
					char tmp = *(char*)value;
					memcpy(bytes, &tmp, eb_info.size_value);
					break;
				}
				case EB_Type_Short: case EB_Type_Short_2: case EB_Type_Short_3:
				{
					short tmp = *(short*)value;
					memcpy(bytes, &tmp, eb_info.size_value);
					break;
				}
				case EB_Type_Long: case EB_Type_Long_2: case EB_Type_Long_3:
				{
					int tmp = *(int*)value;
					memcpy(bytes, &tmp, eb_info.size_value);
					break;
				}
				case EB_Type_Long_Long: case EB_Type_Long_Long_2: case EB_Type_Long_Long_3:
				{
					long long tmp = *(long long*)value;
					memcpy(bytes, &tmp, eb_info.size_value);
					break;
				}
			}
		}
		else
		{
			switch (eb_info.data_type)
			{
				case EB_Type_Unsigned_Char: case EB_Type_Unsigned_Char_2: case EB_Type_Unsigned_Char_3:
				{
					char tmp = *(char*)value;
					memcpy(bytes, &tmp, eb_info.size_value);
					break;
					break;
				}
				case EB_Type_Unsigned_Short: case EB_Type_Unsigned_Short_2:	case EB_Type_Unsigned_Short_3:
				{
					short tmp = *(short*)value;
					memcpy(bytes, &tmp, eb_info.size_value);
					break;
					break;
				}
				case EB_Type_Unsigned_Long: case EB_Type_Unsigned_Long_2: case EB_Type_Unsigned_Long_3:
				{
					int tmp = *(int*)value;
					memcpy(bytes, &tmp, eb_info.size_value);
					break;
					break;
				}
				case EB_Type_Unsigned_Long_Long: case EB_Type_Unsigned_Long_Long_2: case EB_Type_Unsigned_Long_Long_3:
				{		
					long long tmp = *(long long*)value;
					memcpy(bytes, &tmp, eb_info.size_value);
					break;
					break;
				}
			}
		}
	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////



///////laszip function

int LAS_io::File_Flag(char* filename)
{
	if (Is_LAS_File(filename))
	{
		return LAS_FILE_FLAG;
	}
	else if (Is_LAZ_File(filename))
	{
		return LAZ_FILE_FLAG;
	}
	else
	{
		printf("\nInvalid File Format! %s", filename);
		getchar();
		return INVALID_FILE_FLAG;
	}
}

laszip_POINTER LAS_io::Open_LAZ_Reader(char* filename)
{
	laszip_POINTER laszip_reader;
	if (laszip_create(&laszip_reader))
	{
		fprintf(stderr, "\nDLL ERROR: creating laszip reader\n");
		getchar();
	}

	laszip_BOOL is_compressed = 0;
	if (laszip_open_reader(laszip_reader, filename, &is_compressed))
	{
		fprintf(stderr, "\nDLL ERROR: opening laszip reader for '%s'\n", filename);
		getchar();
	}

	return laszip_reader;
}

LAS_Header LAS_io::Load_LAZ_Header(laszip_POINTER laszip_reader)
{
	laszip_header* laz_header;

	if (laszip_get_header_pointer(laszip_reader, &laz_header))
	{
		fprintf(stderr, "\nDLL ERROR: getting header pointer from laszip reader\n");
		getchar();
	}

	if (laz_header->number_of_point_records > 0)
	{
		laz_header->extended_number_of_point_records = laz_header->number_of_point_records;
		for (int i = 0; i < 5; i++)
		{
			laz_header->extended_number_of_points_by_return[i] = laz_header->number_of_points_by_return[i];
		}
	}

	return Convert_Header(*laz_header);
}

LAS_Point LAS_io::Load_A_Point(laszip_POINTER laszip_reader, unsigned short Point_Format_ID)
{
	laszip_point* point;
	if (laszip_get_point_pointer(laszip_reader, &point))
	{
		fprintf(stderr, "DLL ERROR: getting point pointer from laszip reader\n");
		getchar();
	}

	if (laszip_read_point(laszip_reader))
	{
		fprintf(stderr, "DLL ERROR: reading point\n");
		getchar();
	}
	
	return Convert_Point(point, Point_Format_ID);
}

LAS_VLR LAS_io::Load_A_VLR(laszip_POINTER laszip_reader, unsigned int i_vlr)
{
	laszip_header* header;

	if (laszip_get_header_pointer(laszip_reader, &header))
	{
		fprintf(stderr, "\nDLL ERROR: getting header pointer from laszip reader\n");
		getchar();
	}

	return Convert_VLR(header->vlrs[i_vlr]);
}

void LAS_io::Close_LAZ_Reader(laszip_POINTER laszip_reader)
{
	// close the reader

	if (laszip_close_reader(laszip_reader))
	{
		fprintf(stderr, "\nDLL ERROR: closing laszip reader\n");
		getchar();
	}

	// destroy the reader

	if (laszip_destroy(laszip_reader))
	{
		fprintf(stderr, "\nDLL ERROR: destroying laszip reader\n");
		getchar();
	}
}

LAS_Header LAS_io::Convert_Header(laszip_header laz_header)
{
	LAS_Header las_header;
	las_header.File_Signature[0] = 'L';
	las_header.File_Signature[1] = 'A';
	las_header.File_Signature[2] = 'S';
	las_header.File_Signature[3] = 'F';

	las_header.File_Source_ID = laz_header.file_source_ID;
	las_header.Global_Encoding = laz_header.global_encoding;
	las_header.Project_ID_GUID_data_1 = laz_header.project_ID_GUID_data_1;
	las_header.Project_ID_GUID_data_2 = laz_header.project_ID_GUID_data_2;
	las_header.Project_ID_GUID_data_3 = laz_header.project_ID_GUID_data_3;

	for (int i = 0; i < 8; i++)
	{
		las_header.Project_ID_GUID_data_4[i] = (unsigned char)laz_header.project_ID_GUID_data_4[i];
	}

	las_header.Version_Major = laz_header.version_major;
	las_header.Version_Minor = laz_header.version_minor;

	memcpy(las_header.System_Identifier, laz_header.system_identifier, 32);
	memcpy(las_header.Generating_Software, laz_header.generating_software, 32);

	las_header.File_Creation_Year = laz_header.file_creation_year;
	las_header.File_Creation_Day_of_Year = laz_header.file_creation_day;

	las_header.Header_Size = laz_header.header_size;
	las_header.Offset_to_point_data = laz_header.offset_to_point_data;

	las_header.Number_of_Variable_Length_Records = laz_header.number_of_variable_length_records;
	las_header.Point_Data_Format_ID = laz_header.point_data_format;

	las_header.Point_Data_Record_Length = laz_header.point_data_record_length;
	las_header.Legacy_Number_of_point_records = laz_header.number_of_point_records;
	for (int i = 0; i < 5; i++)
	{
		las_header.Legacy_Number_of_points_by_return[i] = laz_header.number_of_points_by_return[i];
	}

	las_header.Number_of_point_records = las_header.Legacy_Number_of_point_records;
	las_header.Legacy_Number_of_point_records = 0;

	for (int i = 0; i < 5; i++)
	{
		las_header.Number_of_points_by_return[i] = las_header.Legacy_Number_of_points_by_return[i];
		las_header.Legacy_Number_of_points_by_return[i] = 0;
	}

	las_header.X_scale_factor = laz_header.x_scale_factor;
	las_header.Y_scale_factor = laz_header.y_scale_factor;
	las_header.Z_scale_factor = laz_header.z_scale_factor;

	las_header.X_offset = laz_header.x_offset;
	las_header.Y_offset = laz_header.y_offset;
	las_header.Z_offset = laz_header.z_offset;

	las_header.Max_X = laz_header.max_x;
	las_header.Min_X = laz_header.min_x;
	las_header.Max_Y = laz_header.max_y;
	las_header.Min_Y = laz_header.min_y;
	las_header.Max_Z = laz_header.max_z;
	las_header.Min_Z = laz_header.min_z;

	if (las_header.Version_Major == 1 && las_header.Version_Minor >= 4)
	{
		las_header.Start_of_Waveform_Data_Packet_Record = laz_header.start_of_waveform_data_packet_record;
		las_header.Start_of_First_Extended_Variable_Length_Record = laz_header.start_of_first_extended_variable_length_record;
		las_header.Number_of_Extended_Variable_Length_Records = laz_header.number_of_extended_variable_length_records;

		las_header.Number_of_point_records = laz_header.extended_number_of_point_records;
		for (int i = 0; i < 15; i++)
		{
			las_header.Number_of_points_by_return[i] = laz_header.extended_number_of_points_by_return[i];
		}
	}

	return las_header;
}

LAS_Point LAS_io::Convert_Point(laszip_point *laz_point, unsigned short Point_Format_ID)
{
	LAS_Point las_point;
	const bool* opt_field = LAS_Optional_Field[Point_Format_ID];

	las_point.X = laz_point->X;
	las_point.Y = laz_point->Y;
	las_point.Z = laz_point->Z;
	las_point.Intensity = laz_point->intensity;


	if (opt_field[0])
	{
		las_point.Return_Number = laz_point->extended_return_number;
		las_point.Number_of_Returns = laz_point->extended_number_of_returns;

		las_point.Classification_Flags = laz_point->extended_classification_flags;
		las_point.Scanner_Channel = laz_point->extended_scanner_channel;
		las_point.Scan_Direction_Flag = laz_point->scan_direction_flag;
		las_point.Edge_of_Flight_Line = laz_point->edge_of_flight_line;
	}
	else
	{
		las_point.Return_Number = laz_point->return_number;
		las_point.Number_of_Returns = laz_point->number_of_returns;
		las_point.Scan_Direction_Flag = laz_point->scan_direction_flag;
		las_point.Edge_of_Flight_Line = laz_point->edge_of_flight_line;
	}

	las_point.Classification = laz_point->classification;
	if (opt_field[1])
	{
		las_point.Scan_Angle_Rank = laz_point->extended_scan_angle;
	}
	else
	{
		las_point.Scan_Angle_Rank = laz_point->scan_angle_rank;
	}
	las_point.User_Data = laz_point->user_data;
	las_point.Point_Source_ID = laz_point->point_source_ID;

	if (opt_field[2])
	{
		las_point.GPS_Time = laz_point->gps_time;
	}

	if (opt_field[3])
	{
		las_point.Red = laz_point->rgb[0];
		las_point.Green = laz_point->rgb[1];
		las_point.Blue = laz_point->rgb[2];
	}

	if (opt_field[4])
	{
		las_point.NIR = laz_point->rgb[3];
	}

	if (opt_field[5])
	{
		int offset = 0;
		memcpy(&las_point.Wave_Packet_Discriptor_Index, laz_point->wave_packet + offset, sizeof(las_point.Wave_Packet_Discriptor_Index));
		offset += sizeof(las_point.Wave_Packet_Discriptor_Index);
		memcpy(&las_point.Byte_offset_to_waveform_data, laz_point->wave_packet + offset, sizeof(las_point.Byte_offset_to_waveform_data));
		offset += sizeof(las_point.Byte_offset_to_waveform_data);
		memcpy(&las_point.Waveform_packet_size_in_bytes, laz_point->wave_packet + offset, sizeof(las_point.Waveform_packet_size_in_bytes));
		offset += sizeof(las_point.Waveform_packet_size_in_bytes);
		memcpy(&las_point.Return_Point_Waveform_Location, laz_point->wave_packet + offset, sizeof(las_point.Return_Point_Waveform_Location));
		offset += sizeof(las_point.Return_Point_Waveform_Location);

		memcpy(&las_point.Xt, laz_point->wave_packet + offset, sizeof(las_point.Xt));
		offset += sizeof(las_point.Xt);
		memcpy(&las_point.Yt, laz_point->wave_packet + offset, sizeof(las_point.Yt));
		offset += sizeof(las_point.Yt);
		memcpy(&las_point.Zt, laz_point->wave_packet + offset, sizeof(las_point.Zt));
		offset += sizeof(las_point.Zt);
	}

	if (laz_point->num_extra_bytes > 0)
	{
		las_point.Extra_Bytes = new LAS_BYTE[laz_point->num_extra_bytes];
		memcpy(las_point.Extra_Bytes, laz_point->extra_bytes, laz_point->num_extra_bytes);
	}

	return las_point;
}

LAS_VLR LAS_io::Convert_VLR(laszip_vlr laz_vlr)
{
	LAS_VLR las_vlr;

	las_vlr.vlr_header.Reserved = laz_vlr.reserved;
	for (int i = 0; i < 16; i++)
	{
		las_vlr.vlr_header.User_ID[i] = laz_vlr.user_id[i];
	}

	las_vlr.vlr_header.Record_ID = laz_vlr.record_id;
	las_vlr.vlr_header.Record_Length_After_Header = laz_vlr.record_length_after_header;
	for (int i = 0; i < 32; i++)
	{
		las_vlr.vlr_header.Description[i] = laz_vlr.description[i];
	}

	las_vlr.vlr_buffer = new LAS_BYTE[laz_vlr.record_length_after_header];
	memcpy(las_vlr.vlr_buffer, laz_vlr.data, laz_vlr.record_length_after_header);

	return las_vlr;
}

//void LAS_io::Write_LAZ_File(char* filename)
//{
//	// create the writer
//	laszip_POINTER laszip_writer;
//
//	if (laszip_create(&laszip_writer))
//	{
//		fprintf(stderr, "DLL ERROR: creating laszip writer\n");
//		getchar();
//	}
//
//	laszip_header* header;
//	if (laszip_get_header_pointer(laszip_writer, &header))
//	{
//		fprintf(stderr, "DLL ERROR: getting header pointer from laszip writer\n");
//		getchar();
//	}
//
//	*header = Convert_Header(las_header);
//	// add VLRs
//	//for (unsigned long i = 0; i < las_header.Number_of_Variable_Length_Records; i++)
//	//{
//	//	laszip_vlr laz_vlr = Convert_VLR(las_vlr[i]); 
//	//	if (laszip_add_vlr(laszip_writer, 
//	//		laz_vlr.user_id,
//	//		laz_vlr.record_id,
//	//		laz_vlr.record_length_after_header,
//	//		laz_vlr.description,
//	//		laz_vlr.data))
//	//	{
//	//		fprintf(stderr, "DLL ERROR: adding VLR to the header\n");
//	//		getchar();
//	//	}
//	//}
//	// open the writer
//	printf("!");  
//	fprintf(stderr, "offset_to_point_data before adding 'height above ground' is : %d\n", (laszip_I32)header->offset_to_point_data);
//	getchar();
//	header->point_data_record_length += 2;
//	printf("\n%i", header->point_data_record_length);
//	laszip_add_attribute(laszip_writer, 3, "test", "test", 1.0, 0.0);
//	fprintf(stderr, "offset_to_point_data before adding 'height above ground' is : %d\n", (laszip_I32)header->offset_to_point_data);
//
//	printf("?");
//
//	laszip_BOOL compress = (strstr(filename, ".laz") != 0);
//	printf("\na");
//	getchar();
//	if (laszip_open_writer(laszip_writer, filename, compress))
//	{
//		fprintf(stderr, "DLL ERROR: opening laszip writer for '%s'\n", filename);
//		getchar();
//	}
//
//	//// get a pointer to the point of the writer that we will populate and write
//	printf("\nb");
//	getchar();
//	laszip_point* point;
//
//	if (laszip_get_point_pointer(laszip_writer, &point))
//	{
//		fprintf(stderr, "DLL ERROR: getting point pointer from laszip writer\n");
//		getchar();
//	}
//	printf("\nc");
//	getchar();
//
//	for (unsigned long long i = 0; i < las_header.Number_of_point_records; i++)
//	{
//		*point = Convert_Point(las_points[i], las_header.Point_Data_Format_ID);
//	
//		// write the point
//		if (laszip_write_point(laszip_writer))
//		{
//			fprintf(stderr, "DLL ERROR: writing point\n");
//			getchar();
//		}
//		
//	}
//
////	*point = Convert_Point(las_points[0], las_header.Point_Data_Format_ID);
//
//	// write the point
//	if (laszip_write_point(laszip_writer))
//	{
//		fprintf(stderr, "DLL ERROR: writing point\n");
//		getchar();
//	}
//
//
//
//
//
//
//	// close the writer
//	printf("\nd");
//	getchar();
//	if (laszip_close_writer(laszip_writer))
//	{
//		fprintf(stderr, "DLL ERROR: closing laszip writer\n");
//		getchar();
//	}
//	printf("\ne");
//	getchar();
//	// destroy the writer
//	if (laszip_destroy(laszip_writer))
//	{
//		fprintf(stderr, "DLL ERROR: destroying laszip writer\n");
//		getchar();
//	}
//	printf("\nf");
//	getchar();
//}

//void LAS_io::Write_LAZ_File(char* filename)
//{
//	// create the writer
//	laszip_POINTER laszip_writer;
//
//	if (laszip_create(&laszip_writer))
//	{
//		fprintf(stderr, "DLL ERROR: creating laszip writer\n");
//		getchar();
//	}
//
//	laszip_header* header;
//	if (laszip_get_header_pointer(laszip_writer, &header))
//	{
//		fprintf(stderr, "DLL ERROR: getting header pointer from laszip writer\n");
//		getchar();
//	}
//
//	header->file_source_ID = las_header.File_Source_ID;
//	header->global_encoding = las_header.Global_Encoding;     // see LAS specification for details
//	header->version_major = las_header.Version_Major;
//	header->version_minor = las_header.Version_Minor;
//	strcpy(header->system_identifier, las_header.System_Identifier);
//	header->file_creation_day = las_header.File_Creation_Day_of_Year;
//	header->file_creation_year = las_header.File_Creation_Year;
//	header->header_size = las_header.Header_Size;
//	header->offset_to_point_data = las_header.Header_Size;
//	header->point_data_format = las_header.Point_Data_Format_ID;
//	header->point_data_record_length = 34 + 1; // three "extra bytes" per point store two additional attributes
//	header->number_of_point_records = 5;           // legacy 32-bit counters should be zero for new point types > 5
//	for (int i = 0; i < 5; i++)
//	{
//		header->number_of_points_by_return[i] = 0;
//	}
//	header->extended_number_of_point_records = 0;
//	
//	header->max_x = las_header.Max_X;
//	header->min_x = las_header.Min_X;
//
//	header->max_y = las_header.Max_Y;
//	header->min_y = las_header.Min_Y;
//
//	header->max_z = las_header.Max_Z;
//	header->min_z = las_header.Min_Z;
//
//	header->x_offset = las_header.X_offset;
//	header->y_offset = las_header.Y_offset;
//	header->z_offset = las_header.Z_offset;
//
//	header->x_scale_factor = las_header.X_scale_factor;
//	header->y_scale_factor = las_header.Y_scale_factor;
//	header->z_scale_factor = las_header.Z_scale_factor;
//
//	LAS_Extra_Bytes_Info* eb_info;
//	unsigned short n_eb;
//	Get_All_Extra_Bytes_Info(eb_info, n_eb);
//	
//	printf("\ndata length: %i\n", header->point_data_record_length);
//
//	laszip_add_attribute(laszip_writer, 0, "name", "description", 1.0, 1);
//	//	laszip_add_attribute(laszip_writer, eb_info[1].data_type, eb_info[1].name, eb_info[1].description, 1.0, 0);
//	if (laszip_add_vlr(laszip_writer, "LASF_Projection", 2112, 0, "intentionally empty OGC WKT", 0))
//	{
//		fprintf(stderr, "DLL ERROR: adding intentionally empty OGC WKT VLR to the header\n");
//	}
//
//	fprintf(stderr, "offset_to_point_data before adding funny VLR is             : %d\n", (laszip_I32)header->offset_to_point_data);
//
//	laszip_BOOL request = 1;
//	if (laszip_request_compatibility_mode(laszip_writer, request))
//	{
//		fprintf(stderr, "DLL ERROR: enabling laszip LAS 1.4 compatibility mode\n");
//	}
//
//
//	laszip_BOOL compress = (strstr(filename, ".laz") != 0);
//
//	if (laszip_open_writer(laszip_writer, filename, compress))
//	{
//		fprintf(stderr, "DLL ERROR: opening laszip writer for '%s'\n", filename);
//		getchar();
//	}
//
//	laszip_point* point;
//	if (laszip_get_point_pointer(laszip_writer, &point))
//	{
//		fprintf(stderr, "DLL ERROR: getting point pointer from laszip writer\n");
//	}
//
//	for (int i = 0; i < 5; i++)
//	{
//		laszip_I64 p_count = 0;
//		laszip_F64 coordinates[3];
//
//		// populate the first point
//
//		coordinates[0] = 630499.95;
//		coordinates[1] = 4834749.17;
//		coordinates[2] = 62.15;
//
//		if (laszip_set_coordinates(laszip_writer, coordinates))
//		{
//			fprintf(stderr, "DLL ERROR: setting coordinates for point %I64d\n", p_count);
//		}
//
//		point->intensity = 60;
//		point->extended_return_number = 2;
//		point->extended_number_of_returns = 2;
//		point->classification = 2;                // it must be set because it "fits" in 5 bits
//		point->extended_classification = 2;
//		point->extended_scan_angle = 3500;
//		point->extended_scanner_channel = 1;
//		point->extended_classification_flags = 8; // overflag flag is set
//		point->gps_time = 53413162.560400;
//		printf("\nc");
//
//		// set attribute 'height above ground' quantized to 0.05 m
//		printf("\n%i", point->num_extra_bytes);
//		getchar();
//		*((unsigned char*)(point->extra_bytes + 0)) = (unsigned char)(12.50 / 0.05);
//		printf("\n%i", point->num_extra_bytes);
//		getchar();
//		printf("\nd");
//		getchar();
//		// write the first point
//
//		if (laszip_write_point(laszip_writer))
//		{
//			fprintf(stderr, "DLL ERROR: writing point %I64d\n", p_count);
//		}	
//		printf("\ne");
//		getchar();
//		p_count++;
//	}
//
//	fprintf(stderr, "offset_to_point_data before adding 'height above ground' is : %d\n", (laszip_I32)header->offset_to_point_data);
//
//
//	//// get a pointer to the point of the writer that we will populate and write
//	//	*point = Convert_Point(las_points[0], las_header.Point_Data_Format_ID);
//
//	// write the point
//	if (laszip_write_point(laszip_writer))
//	{
//		fprintf(stderr, "DLL ERROR: writing point\n");
//		getchar();
//	}
//
//	// close the writer
//	printf("\nd");
//	getchar();
//	if (laszip_close_writer(laszip_writer))
//	{
//		fprintf(stderr, "DLL ERROR: closing laszip writer\n");
//		getchar();
//	}
//	printf("\ne");
//	getchar();
//	// destroy the writer
//	if (laszip_destroy(laszip_writer))
//	{
//		fprintf(stderr, "DLL ERROR: destroying laszip writer\n");
//		getchar();
//	}
//	printf("\nf");
//	getchar();
//}

void LAS_io::Write_LAZ_File(char* filename)
{
	// create the writer
	laszip_POINTER laszip_writer;
	
	if (laszip_create(&laszip_writer))
	{
		fprintf(stderr, "DLL ERROR: creating laszip writer\n");
		getchar();
	}

	Set_LAZ_Header(laszip_writer);
	Set_LAZ_VLRs(laszip_writer);

	laszip_BOOL compress = true;
	if (laszip_open_writer(laszip_writer, filename, compress))
	{
		fprintf(stderr, "DLL ERROR: opening laszip writer for '%s'\n", filename);
		getchar();
	}

	Set_LAZ_Points(laszip_writer);

	// close the writer

	if (laszip_close_writer(laszip_writer))
	{
		fprintf(stderr, "DLL ERROR: closing laszip writer\n");
		getchar();
	}

	// destroy the writer
	if (laszip_destroy(laszip_writer))
	{
		fprintf(stderr, "DLL ERROR: destroying laszip writer\n");
		getchar();
	}
}

void LAS_io::Set_LAZ_Header(laszip_POINTER laszip_writer)
{
	laszip_header* header;
	if (laszip_get_header_pointer(laszip_writer, &header))
	{
		fprintf(stderr, "DLL ERROR: getting header pointer from laszip writer\n");
		getchar();
	}
	*header = Convert_Header(las_header);
	header->offset_to_point_data = header->header_size;
	header->number_of_variable_length_records = 0;
}

void LAS_io::Set_LAZ_VLRs(laszip_POINTER laszip_writer)
{
	LAS_Extra_Bytes_Info* eb_info;
	unsigned short n_eb;
	long long idx_vlr_eb = Get_All_Extra_Bytes_Info(eb_info, n_eb);

	if (n_eb > 0)
	{
		for (unsigned short i_eb = 0; i_eb < n_eb; i_eb++)
		{
			for (int i_val = 0; i_val < eb_info[i_eb].n_values; i_val++)
			{
				unsigned short type = 0;
				char name[32];
				char description[32];
				double scale = eb_info[i_eb].options_scale_bit ? eb_info[i_eb].scale[i_val] : 1.0;
				double offset = eb_info[i_eb].options_offset_bit ? eb_info[i_eb].offset[i_val] : 0.0;

				memcpy(name, eb_info[i_eb].name, 32);
				memcpy(description, eb_info[i_eb].description, 32);

				if (eb_info[i_eb].n_values > 1)
				{
					type = eb_info[i_eb].data_type % 10 - 1;
					snprintf(name, 32, "%s_%i%c", eb_info[i_eb].name, i_val + 1, '\0');
					snprintf(description, 32, "%s_%i%c", eb_info[i_eb].description, i_val + 1, '\0');	
				}
				else
				{
					type = eb_info[i_eb].data_type - 1; // this is the way laszip writes.
				}

				if (laszip_add_attribute(laszip_writer, type, name, description, scale, offset))
				{
					fprintf(stderr, "DLL ERROR: adding %s attribute\n", eb_info[i_eb].name);
					getchar();
				}
			}
		}
		delete[] eb_info;
	}

	for (unsigned int i_vlr = 0; i_vlr < las_header.Number_of_Variable_Length_Records; i_vlr++)
	{
		if (i_vlr != idx_vlr_eb)
		{
			if (laszip_add_vlr(laszip_writer,
				las_vlr[i_vlr].vlr_header.User_ID,
				las_vlr[i_vlr].vlr_header.Record_ID,
				las_vlr[i_vlr].vlr_header.Record_Length_After_Header,
				las_vlr[i_vlr].vlr_header.Description,
				las_vlr[i_vlr].vlr_buffer))
			{
				fprintf(stderr, "DLL ERROR: adding VLR to the header\n");
			}
		}
	}
}

void LAS_io::Set_LAZ_Points(laszip_POINTER laszip_writer)
{
	laszip_point* point;
	if (laszip_get_point_pointer(laszip_writer, &point))
	{
		fprintf(stderr, "DLL ERROR: getting point pointer from laszip writer\n");
		getchar();
	}

	for (unsigned long long pid = 0; pid < las_header.Number_of_point_records; pid++)
	{
		LAS2LAZ_Point(las_points[pid], point, las_header.Point_Data_Format_ID);
		if (laszip_write_point(laszip_writer))
		{
			fprintf(stderr, "DLL ERROR: writing point %llu\n", pid);
			getchar();
		}
	}
}

laszip_header LAS_io::Convert_Header(LAS_Header header)
{
	laszip_header laz_header;
	////////////
		laz_header.user_data_after_header_size = 0;
		laz_header.user_data_after_header = 0;
		laz_header.user_data_in_header_size = 0;
		laz_header.user_data_in_header = 0;
		laz_header.vlrs = 0;
	///////////////
	laz_header.file_source_ID = header.File_Source_ID;
	laz_header.global_encoding = header.Global_Encoding;
	laz_header.project_ID_GUID_data_1 = header.Project_ID_GUID_data_1;
	laz_header.project_ID_GUID_data_3 = header.Project_ID_GUID_data_2;
	laz_header.project_ID_GUID_data_2 = header.Project_ID_GUID_data_3;

	for (int i = 0; i < 8; i++)
	{
		laz_header.project_ID_GUID_data_4[i] = (char)header.Project_ID_GUID_data_4[i];
	}

	laz_header.version_major = header.Version_Major;
	laz_header.version_minor = header.Version_Minor;

	memcpy(laz_header.system_identifier, header.System_Identifier, sizeof(header.System_Identifier));
	memcpy(laz_header.generating_software, header.Generating_Software, sizeof(header.Generating_Software));

	laz_header.file_creation_year = header.File_Creation_Year;
	laz_header.file_creation_day = header.File_Creation_Day_of_Year;

	laz_header.header_size = header.Header_Size;
	laz_header.offset_to_point_data = header.Offset_to_point_data;

	laz_header.number_of_variable_length_records = header.Number_of_Variable_Length_Records;
	laz_header.point_data_format = header.Point_Data_Format_ID;
	laz_header.point_data_record_length = header.Point_Data_Record_Length;

	if (header.Version_Major == 1 && header.Version_Minor < 4)
	{
		laz_header.number_of_point_records = (unsigned int)header.Number_of_point_records;

		for (int i = 0; i < 5; i++)
		{
			laz_header.number_of_points_by_return[i] = (unsigned int)header.Number_of_points_by_return[i];
		}
	}
	else
	{
		laz_header.number_of_point_records = header.Legacy_Number_of_point_records;
		for (int i = 0; i < 5; i++)
		{
			laz_header.number_of_points_by_return[i] = header.Legacy_Number_of_points_by_return[i];
		}
	}

	laz_header.x_scale_factor = header.X_scale_factor;
	laz_header.y_scale_factor = header.Y_scale_factor;
	laz_header.z_scale_factor = header.Z_scale_factor;

	laz_header.x_offset = header.X_offset;
	laz_header.y_offset = header.Y_offset;
	laz_header.z_offset = header.Z_offset;

	laz_header.max_x = header.Max_X;
	laz_header.min_x = header.Min_X;

	laz_header.max_y = header.Max_Y;
	laz_header.min_y = header.Min_Y;

	laz_header.max_z = header.Max_Z;
	laz_header.min_z = header.Min_Z;

	laz_header.start_of_waveform_data_packet_record = header.Start_of_Waveform_Data_Packet_Record;
	laz_header.start_of_first_extended_variable_length_record = header.Start_of_First_Extended_Variable_Length_Record;
	laz_header.number_of_extended_variable_length_records = header.Number_of_Extended_Variable_Length_Records;

	if (header.Version_Major == 1 && header.Version_Minor >= 4)
	{
		laz_header.extended_number_of_point_records = header.Number_of_point_records;
		for (int i = 0; i < 15; i++)
		{
			laz_header.extended_number_of_points_by_return[i] = header.Number_of_points_by_return[i];
		}
	}
	else
	{
		laz_header.extended_number_of_point_records = 0;
		for (int i = 0; i < 15; i++)
		{
			laz_header.extended_number_of_points_by_return[i] = 0;
		}
	}

	return laz_header;
}

void LAS_io::LAS2LAZ_Point(LAS_Point las_point, laszip_point *laz_point, unsigned short Point_Format_ID)
{
	const bool* opt_field = LAS_Optional_Field[Point_Format_ID];

	laz_point->X = las_point.X;
	laz_point->Y = las_point.Y;
	laz_point->Z = las_point.Z;
	laz_point->intensity = las_point.Intensity;

	if (opt_field[LASIO_PDF_OPT_Extended_Flight_Info])
	{
		laz_point->extended_return_number = las_point.Return_Number;
		laz_point->extended_number_of_returns = las_point.Number_of_Returns;

		laz_point->extended_classification_flags = las_point.Classification_Flags;
		laz_point->extended_scanner_channel = las_point.Scanner_Channel;
		laz_point->scan_direction_flag = las_point.Scan_Direction_Flag;
		laz_point->edge_of_flight_line = las_point.Edge_of_Flight_Line;
	}
	else
	{
		laz_point->return_number = las_point.Return_Number;
		laz_point->number_of_returns = las_point.Number_of_Returns;
		laz_point->scan_direction_flag = las_point.Scan_Direction_Flag;
		laz_point->edge_of_flight_line = las_point.Edge_of_Flight_Line;
	}

	laz_point->classification = las_point.Classification;

	if (opt_field[LASIO_PDF_OPT_Extended_Scan_Angle])
	{
		laz_point->extended_scan_angle = las_point.Scan_Angle_Rank;
	}
	else
	{
		laz_point->scan_angle_rank = (char)las_point.Scan_Angle_Rank;
	}

	laz_point->user_data = las_point.User_Data;
	laz_point->point_source_ID = las_point.Point_Source_ID;

	if (opt_field[LASIO_PDF_OPT_GPS_Time])
	{
		laz_point->gps_time = las_point.GPS_Time;
	}

	if (opt_field[LASIO_PDF_OPT_Color])
	{
		laz_point->rgb[0] = las_point.Red;
		laz_point->rgb[1] = las_point.Green;
		laz_point->rgb[2] = las_point.Blue;
	}

	if (opt_field[LASIO_PDF_OPT_NIR])
	{
		laz_point->rgb[3] = las_point.NIR;
	}

	if (opt_field[LASIO_PDF_OPT_Waveform])
	{
		int offset = 0;
		memcpy(laz_point->wave_packet + offset, &las_point.Wave_Packet_Discriptor_Index, sizeof(las_point.Wave_Packet_Discriptor_Index));
		offset += sizeof(las_point.Wave_Packet_Discriptor_Index);
		memcpy(laz_point->wave_packet + offset, &las_point.Byte_offset_to_waveform_data, sizeof(las_point.Byte_offset_to_waveform_data));
		offset += sizeof(las_point.Byte_offset_to_waveform_data);
		memcpy(laz_point->wave_packet + offset, &las_point.Waveform_packet_size_in_bytes, sizeof(las_point.Waveform_packet_size_in_bytes));
		offset += sizeof(las_point.Waveform_packet_size_in_bytes);
		memcpy(laz_point->wave_packet + offset, &las_point.Return_Point_Waveform_Location, sizeof(las_point.Return_Point_Waveform_Location));
		offset += sizeof(las_point.Return_Point_Waveform_Location);

		memcpy(laz_point->wave_packet + offset, &las_point.Xt, sizeof(las_point.Xt));
		offset += sizeof(las_point.Xt);
		memcpy(laz_point->wave_packet + offset, &las_point.Yt, sizeof(las_point.Yt));
		offset += sizeof(las_point.Yt);
		memcpy(laz_point->wave_packet + offset, &las_point.Zt, sizeof(las_point.Zt));
		offset += sizeof(las_point.Zt);
	}

	if (laz_point->num_extra_bytes > 0)
	{
		memcpy(laz_point->extra_bytes, las_point.Extra_Bytes, laz_point->num_extra_bytes);
	}
}

void LAS_io::LAS2LAS_Point(LAS_Point las_point_src, LAS_Point &las_point_des, unsigned short size_extra_bytes)
{
	LAS_BYTE* extra_bytes = las_point_des.Extra_Bytes;
	memcpy(extra_bytes, las_point_src.Extra_Bytes, size_extra_bytes);
	las_point_des = las_point_src;
	las_point_des.Extra_Bytes = extra_bytes;
}

void LAS_io::LAS2LAS_VLR(LAS_VLR las_vlr_src, LAS_VLR& las_vlr_des)
{
	las_vlr_des.vlr_header = las_vlr_src.vlr_header;
	las_vlr_des.vlr_buffer = new LAS_BYTE[las_vlr_src.vlr_header.Record_Length_After_Header];
	memcpy(las_vlr_des.vlr_buffer, las_vlr_src.vlr_buffer, las_vlr_src.vlr_header.Record_Length_After_Header);
}

void LAS_io::LAS2LAS_EVLR(LAS_EVLR las_evlr_src, LAS_EVLR& las_evlr_des)
{
	LAS2LAS_VLR(las_evlr_src, las_evlr_des);
}
/////////////////////////

void LAS_io::Show_LAS_Header()
{
	printf("\nHeader Info:");

	switch (las_header.Version_Major)
	{
		case 1:
			Print_LAS_Header_1_X(las_header);
			break;
	}
}

void LAS_io::Show_Point(unsigned long long i)
{
	LAS_Point p = las_points[i];
	double x = x_record2coord(p.X);
	double y = x_record2coord(p.Y);
	double z = x_record2coord(p.Z);

	printf("\nPoint #%llu: ", i);
	printf("\n\tGX: %lf", x);
	printf("\n\tGY: %lf", y);
	printf("\n\tGZ: %lf", z);

	Print_LAS_Point(p, las_header.Point_Data_Format_ID);
}

void LAS_io::Show_VLR(unsigned int i)
{
	LAS_VLR vlr = las_vlr[i];
	printf("\nVariable Length Record #%i: ", i);
	Print_LAS_VLR(vlr);
}

void LAS_io::Show_EVLR(unsigned int i)
{
	LAS_EVLR evlr = las_evlr[i];
	printf("\nExtended Variable Length Record #%i: ", i);
	Print_LAS_EVLR(evlr);
}

void LAS_io::Show_Extra_Bytes_VLR()
{
	printf("\nExtra Bytes Info: ");
	long long idx = idx_vlr_extra_bytes();

	if (idx > -1)
	{
		Show_VLR((unsigned int)idx);

		unsigned short n_eb = 0;
		LAS_Extra_Bytes_Info* eb_info = NULL;
		Get_All_Extra_Bytes_Info(eb_info, n_eb);
		
		for (unsigned short i_eb = 0; i_eb < n_eb; i_eb++)
		{
			printf("\nExtra Bytes #%i: ", i_eb);
			printf("\n\treserved:    %i %i", eb_info[i_eb].reserved[0], eb_info[i_eb].reserved[1]);
			printf("\n\tdata_type:   %i ", eb_info[i_eb].data_type);

			switch (eb_info[i_eb].data_type)
			{
			case EB_Type_Undocumented:
			{
				printf("[ Undocumented ]");
				break;
			}
			case EB_Type_Unsigned_Char:
			{
				printf("[ unsigned char ]");
				break;
			}
			case EB_Type_Unsigned_Char_2:
			{
				printf("[ unsigned char[2] ]");
				break;
			}
			case EB_Type_Unsigned_Char_3:
			{
				printf("[ unsigned char[3] ]");
				break;
			}
			case EB_Type_Unsigned_Short:
			{
				printf("[ unsigned short ]");
				break;
			}
			case EB_Type_Unsigned_Short_2:
			{
				printf("[ unsigned short[2] ]");
				break;
			}
			case EB_Type_Unsigned_Short_3:
			{
				printf("[ unsigned short[3] ]");
				break;
			}
			case EB_Type_Unsigned_Long:
			{
				printf("[ unsigned long ]");
				break;
			}
			case EB_Type_Unsigned_Long_2:
			{
				printf("[ unsigned long[2] ]");
				break;
			}
			case EB_Type_Unsigned_Long_3:
			{
				printf("[ unsigned long[3]]");
				break;
			}
			case EB_Type_Unsigned_Long_Long:
			{
				printf("[ unsigned long long ]");
				break;
			}
			case EB_Type_Unsigned_Long_Long_2:
			{
				printf("[ unsigned long long[2] ]");
				break;
			}
			case EB_Type_Unsigned_Long_Long_3:
			{
				printf("[ unsigned long long[3] ]");
				break;
			}
			case EB_Type_Char:
			{
				printf("[ char ]");
				break;
			}
			case EB_Type_Char_2:
			{
				printf("[ char[2] ]");
				break;
			}
			case EB_Type_Char_3:
			{
				printf("[ char[3] ]");
				break;
			}
			case EB_Type_Short:
			{
				printf("[ short ]");
				break;
			}
			case EB_Type_Short_2:
			{
				printf("[ short[2] ]");
				break;
			}
			case EB_Type_Short_3:
			{
				printf("[ short[3] ]");
				break;
			}
			case EB_Type_Long:
			{
				printf("[ long ]");
				break;
			}
			case EB_Type_Long_2:
			{
				printf("[ long[2] ]");
				break;
			}
			case EB_Type_Long_3:
			{
				printf("[ long[3] ]");
				break;
			}
			case EB_Type_Long_Long:
			{
				printf("[ long long ]");
				break;
			}
			case EB_Type_Long_Long_2:
			{
				printf("[ Ulong long[2] ]");
				break;
			}
			case EB_Type_Long_Long_3:
			{
				printf("[ long long[3] ]");
				break;
			}
			case EB_Type_Float:
			{
				printf("[ float ]");
				break;
			}
			case EB_Type_Float_2:
			{
				printf("[ float[2] ]");
				break;
			}
			case EB_Type_Float_3:
			{
				printf("[ float[3] ]");
				break;
			}
			case EB_Type_Double:
			{
				printf("[ double ]");
				break;
			}
			case EB_Type_Double_2:
			{
				printf("[ double[2] ]");
				break;
			}
			case EB_Type_Double_3:
			{
				printf("[ Double[3] ]");
				break;
			}
			};

			printf("\n\tname:        %s", eb_info[i_eb].name);
			printf("\n\tunused:      %i %i %i %i", eb_info[i_eb].unused[0], eb_info[i_eb].unused[1], eb_info[i_eb].unused[2], eb_info[i_eb].unused[3]);

			printf("\n\tno_data:     ");

			if (eb_info[i_eb].options_no_data_bit)
			{
				if (eb_info[i_eb].floating_or_integer)
				{
					printf("%lf ", double_from_byte_8(eb_info[i_eb].no_data[0]));
					if (eb_info[i_eb].n_values > 1)
					{
						printf("%lf ", double_from_byte_8(eb_info[i_eb].no_data[1]));
						if (eb_info[i_eb].n_values > 2)
						{
							printf("%lf", double_from_byte_8(eb_info[i_eb].no_data[2]));
						}
					}
				}
				else
				{
					if (eb_info[i_eb].signed_or_unsigned)
					{
						printf("%lli ", int64_from_byte_8(eb_info[i_eb].no_data[0]));
						if (eb_info[i_eb].n_values > 1)
						{
							printf("%lli ", int64_from_byte_8(eb_info[i_eb].no_data[1]));
							if (eb_info[i_eb].n_values > 2)
							{
								printf("%lli", int64_from_byte_8(eb_info[i_eb].no_data[2]));
							}
						}
					}
					else
					{
						printf("%lli ", uint64_from_byte_8(eb_info[i_eb].no_data[0]));
						if (eb_info[i_eb].n_values > 1)
						{
							printf("%lli ", uint64_from_byte_8(eb_info[i_eb].no_data[1]));
							if (eb_info[i_eb].n_values > 2)
							{
								printf("%lli", uint64_from_byte_8(eb_info[i_eb].no_data[2]));
							}
						}
					}
				}
			}
			else
			{
				printf("Invalid");
			}

			printf("\n\tmin:         ");
			if (eb_info[i_eb].options_min_bit)
			{
				if (eb_info[i_eb].floating_or_integer)
				{
					printf("%lf ", double_from_byte_8(eb_info[i_eb].min[0]));
					if (eb_info[i_eb].n_values > 1)
					{
						printf("%lf ", double_from_byte_8(eb_info[i_eb].min[1]));
						if (eb_info[i_eb].n_values > 2)
						{
							printf("%lf", double_from_byte_8(eb_info[i_eb].min[2]));
						}
					}
				}
				else
				{
					if (eb_info[i_eb].signed_or_unsigned)
					{
						printf("%lli ", int64_from_byte_8(eb_info[i_eb].min[0]));
						if (eb_info[i_eb].n_values > 1)
						{
							printf("%lli ", int64_from_byte_8(eb_info[i_eb].min[1]));
							if (eb_info[i_eb].n_values > 2)
							{
								printf("%lli", int64_from_byte_8(eb_info[i_eb].min[2]));
							}
						}
					}
					else
					{
						printf("%lli ", uint64_from_byte_8(eb_info[i_eb].min[0]));
						if (eb_info[i_eb].n_values > 1)
						{
							printf("%lli ", uint64_from_byte_8(eb_info[i_eb].min[1]));
							if (eb_info[i_eb].n_values > 2)
							{
								printf("%lli", uint64_from_byte_8(eb_info[i_eb].min[2]));
							}
						}
					}
				}
			}
			else
			{
				printf("Invalid");
			}

			printf("\n\tmax:         ");
			if (eb_info[i_eb].options_max_bit)
			{
				if (eb_info[i_eb].floating_or_integer)
				{
					printf("%lf ", double_from_byte_8(eb_info[i_eb].no_data[0]));
					if (eb_info[i_eb].n_values > 1)
					{
						printf("%lf ", double_from_byte_8(eb_info[i_eb].no_data[1]));
						if (eb_info[i_eb].n_values > 2)
						{
							printf("%lf", double_from_byte_8(eb_info[i_eb].no_data[2]));
						}
					}
				}
				else
				{
					if (eb_info[i_eb].signed_or_unsigned)
					{
						printf("%lli ", int64_from_byte_8(eb_info[i_eb].no_data[0]));
						if (eb_info[i_eb].n_values > 1)
						{
							printf("%lli ", int64_from_byte_8(eb_info[i_eb].no_data[1]));
							if (eb_info[i_eb].n_values > 2)
							{
								printf("%lli", int64_from_byte_8(eb_info[i_eb].no_data[2]));
							}
						}
					}
					else
					{
						printf("%lli ", uint64_from_byte_8(eb_info[i_eb].no_data[0]));
						if (eb_info[i_eb].n_values > 1)
						{
							printf("%lli ", uint64_from_byte_8(eb_info[i_eb].no_data[1]));
							if (eb_info[i_eb].n_values > 2)
							{
								printf("%lli", uint64_from_byte_8(eb_info[i_eb].no_data[2]));
							}
						}
					}
				}
			}
			else
			{
				printf("Invalid");
			}

			printf("\n\tscale:       ");
			if (eb_info[i_eb].options_scale_bit)
			{
				printf("%lf ", eb_info[i_eb].scale[0]);
				if (eb_info[i_eb].n_values > 1)
				{
					printf("%lf ", eb_info[i_eb].scale[1]);
					if (eb_info[i_eb].n_values > 2)
					{
						printf("%lf", eb_info[i_eb].scale[2]);
					}
				}
			}
			else
			{
				printf("Invalid");
			}

			printf("\n\toffset:      ");
			if (eb_info[i_eb].options_offset_bit)
			{
				printf("%lf ", eb_info[i_eb].offset[0]);
				if (eb_info[i_eb].n_values > 1)
				{
					printf("%lf ", eb_info[i_eb].offset[1]);
					if (eb_info[i_eb].n_values > 2)
					{
						printf("%lf", eb_info[i_eb].offset[2]);
					}
				}
			}
			else
			{
				printf("Invalid");
			}

			printf("\n\tdescription: %s", eb_info[i_eb].description);

			printf("\n\tbyte offset: %i", eb_info[i_eb].byte_pos);
		}

		if (eb_info)
		{
			delete[] eb_info;
			eb_info = NULL;
		}
	}
	else
	{
		printf("\nNo Extra Bytes VLR!");
	}

	printf("\n");
}

int LAS_io::x_coord2record(double x)
{
	return coord2record(las_header.X_offset, las_header.X_scale_factor, x);
}

int LAS_io::y_coord2record(double y)
{
	return coord2record(las_header.Y_offset, las_header.Y_scale_factor, y);
}

int LAS_io::z_coord2record(double z)
{
	return coord2record(las_header.Z_offset, las_header.Z_scale_factor, z);
}

double LAS_io::x_record2coord(int x)
{
	return record2coord(las_header.X_offset, las_header.X_scale_factor, x);
}

double LAS_io::y_record2coord(int y)
{
	return record2coord(las_header.Y_offset, las_header.Y_scale_factor, y);
}

double LAS_io::z_record2coord(int z)
{
	return record2coord(las_header.Z_offset, las_header.Z_scale_factor, z);
}

void LAS_io::Delete_Data()
{
	if (las_points)
	{
		if (size_extra_bytes(true) > 0)
		{
			Clear_Extra_Bytes();
		}
		delete[] las_points;
		las_points = NULL;
	}

	if (las_vlr)
	{
		for (unsigned int i = 0; i < las_header.Number_of_Variable_Length_Records; i++)
		{
			if (las_vlr[i].vlr_buffer)
			{
				delete[] las_vlr[i].vlr_buffer;
				las_vlr[i].vlr_buffer = NULL;
			}
		}
		delete[] las_vlr;
		las_vlr = NULL;
	}

	if (las_evlr)
	{
		delete[] las_evlr;
		las_evlr = NULL;
	}
}

int LAS_io::coord2record(double offset, double scale_factor, double coord)
{
	return (int)round((coord - offset) / scale_factor);
}

double LAS_io::record2coord(double offset, double scale_factor, int record)
{
	return (record * scale_factor + offset);
}

void LAS_io::Sort_LAS_by_GPS_Time(LAS_Point *pts, unsigned long long npts)
{
	sort(pts, pts + npts, lasio_compare_gps_time);
}

void LAS_io::Sort_LAS_by_Point_Source_ID(LAS_Point *pts, unsigned long long npts)
{
	sort(pts, pts + npts, lasio_compare_point_source_id);
}

void LAS_io::Sort_LAS_by_Classification(LAS_Point* pts, unsigned long long npts)
{
	sort(pts, pts + npts, lasio_compare_classification);
}

bool lasio_compare_point_source_id(LAS_Point a, LAS_Point b)
{
	return (a.Point_Source_ID < b.Point_Source_ID);
}

bool lasio_compare_gps_time(LAS_Point a, LAS_Point b)
{
	return (a.GPS_Time < b.GPS_Time);
}

bool lasio_compare_classification(LAS_Point a, LAS_Point b)
{
	return (a.Classification < b.Classification);
}

unsigned char LAS_io::Union_Point_Data_Format(unsigned char Format_ID_1, unsigned char Format_ID_2)
{
	bool field[LASIO_N_PDF_OPT_FIELD];
	for (unsigned short i = 0; i < LASIO_N_PDF_OPT_FIELD; i++)
	{
		field[i] = LAS_Optional_Field[Format_ID_1][i] || LAS_Optional_Field[Format_ID_2][i];
	}
	return Match_Point_Data_Format(field);
}

unsigned char LAS_io::Intersect_Point_Data_Format(unsigned char Format_ID_1, unsigned char Format_ID_2)
{
	bool field[LASIO_N_PDF_OPT_FIELD];
	for (unsigned short i = 0; i < LASIO_N_PDF_OPT_FIELD; i++)
	{
		field[i] = LAS_Optional_Field[Format_ID_1][i] && LAS_Optional_Field[Format_ID_2][i];
	}
	return Match_Point_Data_Format(field);
}

unsigned char LAS_io::Match_Point_Data_Format(bool opt_field[LASIO_N_PDF_OPT_FIELD])
{
	for (unsigned char i = 0; i < LASIO_N_POINT_DATA_FORMAT; i++)
	{
		bool match = true;
		for (unsigned short j = 0; j < LASIO_N_PDF_OPT_FIELD; j++)
		{
			if (LAS_Optional_Field[i][j] != opt_field[j])
			{
				match = false;
				break;
			}
		}
		if (match)
		{
			return i;
		}
	}
	return LASIO_N_POINT_DATA_FORMAT;
}

unsigned char LAS_io::Match_Point_Data_Format(unsigned short point_size)
{
	for (unsigned char i = 0; i < LASIO_N_POINT_DATA_FORMAT; i++)
	{
		if (point_size == LAS_Point_Size[i])
		{
			return i;
		}
	}
	return LASIO_N_POINT_DATA_FORMAT;
}

unsigned char LAS_io::add_opt_field(unsigned char Format_ID, LASIO_PDF_OPT_FIELD las_opt_field)
{
	bool opt_field[LASIO_N_PDF_OPT_FIELD];
	memcpy(opt_field, LAS_Optional_Field[Format_ID], LASIO_N_PDF_OPT_FIELD);
	opt_field[las_opt_field] = true;
	return Match_Point_Data_Format(opt_field);
}

unsigned char LAS_io::remove_opt_field(unsigned char Format_ID, LASIO_PDF_OPT_FIELD las_opt_field)
{
	bool opt_field[LASIO_N_PDF_OPT_FIELD];
	memcpy(opt_field, LAS_Optional_Field[Format_ID], LASIO_N_PDF_OPT_FIELD);
	opt_field[las_opt_field] = false;
	return Match_Point_Data_Format(opt_field);
}

void LAS_io::Set_Point_Data_Format(LAS_Header &header, unsigned char nu_format_id)
{
	unsigned short size = LAS_Point_Size[header.Point_Data_Format_ID];
	unsigned short nu_size = LAS_Point_Size[nu_format_id];

	header.Point_Data_Format_ID = nu_format_id;
	header.Point_Data_Record_Length += nu_size - size;
}

void LAS_io::Print_LAS_Header_1_X(LAS_Header h)
{
	printf("\n\tFile_Signature: %c%c%c%c", h.File_Signature[0], h.File_Signature[1], h.File_Signature[2], h.File_Signature[3]);
	printf("\n\tFile_Source_ID: %i", h.File_Source_ID);
	printf("\n\tGlobal_Encoding: %i", h.Global_Encoding);
	printf("\n\tProject_ID_GUID_data: {%08X-%04hX-%04hX-%02X %02X-%02X%02X%02X%02X%02X%02X}",
		h.Project_ID_GUID_data_1, h.Project_ID_GUID_data_2, h.Project_ID_GUID_data_3,
		h.Project_ID_GUID_data_4[0], h.Project_ID_GUID_data_4[1], h.Project_ID_GUID_data_4[2], h.Project_ID_GUID_data_4[3],
		h.Project_ID_GUID_data_4[4], h.Project_ID_GUID_data_4[5], h.Project_ID_GUID_data_4[6], h.Project_ID_GUID_data_4[7]);

	printf("\n\tSystem_Identifier: %s", h.System_Identifier);
	printf("\n\tGenerating_Software: %s", h.Generating_Software);
	printf("\n\tFile_Creation_Time_Year (Day): %i (%i)", h.File_Creation_Year, h.File_Creation_Day_of_Year);

	printf("\n\tVersion (Header Size): %i.%i (%i)", h.Version_Major, h.Version_Minor, h.Header_Size);
	printf("\n\tPoint_Data_Format_ID (Record_Length): %i (%i)", h.Point_Data_Format_ID, h.Point_Data_Record_Length);

	printf("\n\tOffset_to_point_data: %i", h.Offset_to_point_data);
	printf("\n\tNumber_of_Variable_Length_Records: %i", h.Number_of_Variable_Length_Records);

	if (h.Version_Minor < 4)
	{
		printf("\n\tNumber_of_point_records: %llu", h.Number_of_point_records);
		for (int i = 0; i < 5; i++)
		{
			printf("\n\tNumber_of_points_by_return[%i]: %llu", i, h.Number_of_points_by_return[i]);
		}
	}
	else
	{
		printf("\n\tLegacy_Number_of_point_records: %u", h.Legacy_Number_of_point_records);
		for (int i = 0; i < 5; i++)
		{
			printf("\n\tLegacy_Number_of_points_by_return[%i]: %u", i, h.Legacy_Number_of_points_by_return[i]);
		}
	}

	printf("\n\tScale_factor: [%.6lf, %.6lf, %.6lf]", h.X_scale_factor, h.Y_scale_factor, h.Z_scale_factor);
	printf("\n\tOffset:       [%.6lf, %.6lf, %.6lf]", h.X_offset, h.Y_offset, h.Z_offset);
	printf("\n\tMin:          [%.6lf, %.6lf, %.6lf]", h.Min_X, h.Min_Y, h.Min_Z);
	printf("\n\tMax:          [%.6lf, %.6lf, %.6lf]", h.Max_X, h.Max_Y, h.Max_Z);

	if (h.Version_Minor >= 4)
	{
		printf("\n\tStart_of_Waveform_Data_Packet_Record: %llu", h.Start_of_Waveform_Data_Packet_Record);
		printf("\n\tStart_of_First_Extended_Variable_Length_Record: %llu", h.Start_of_First_Extended_Variable_Length_Record);
		printf("\n\tNumber_of_Extended_Variable_Length_Records: %i", h.Number_of_Extended_Variable_Length_Records);
		printf("\n\tNumber_of_point_records: %llu", h.Number_of_point_records);
		for (int i = 0; i < 15; i++)
		{
			printf("\n\tNumber_of_points_by_return[%i]: %llu", i, h.Number_of_points_by_return[i]);
		}
	}

	printf("\n");
}

void LAS_io::Print_LAS_Point(LAS_Point p, unsigned short Point_Format_ID)
{
	const bool* field = LAS_Optional_Field[Point_Format_ID];

	printf("\n\tX: %i", p.X);
	printf("\n\tY: %i", p.Y);
	printf("\n\tZ: %i", p.Z);

	printf("\n\tIntensity: %i", p.Intensity);

	printf("\n\tReturn_Number: %i", p.Return_Number);
	printf("\n\tNumber_of_Returns: %i", p.Number_of_Returns);

	if (field[0])
	{
		printf("\n\tClassification_Flags: %i", p.Classification_Flags);
		printf("\n\tScanner_Channel: %i", p.Scanner_Channel);
	}

	printf("\n\tScan_Direction_Flag: %i", p.Scan_Direction_Flag);
	printf("\n\tEdge_of_Flight_Line: %i", p.Edge_of_Flight_Line);

	printf("\n\tClassification: %i", p.Classification);
	printf("\n\tScan_Angle_Rank: %i", p.Scan_Angle_Rank);
	printf("\n\tUser_Data: %i", p.User_Data);
	printf("\n\tPoint_Source_ID: %i", p.X);

	if (field[2])
	{
		printf("\n\tGPS_Time: %lf", p.GPS_Time);
	}

	if (field[3])
	{
		printf("\n\tRed: %i", p.Red);
		printf("\n\tGreen: %i", p.Green);
		printf("\n\tBlue: %i", p.Blue);
	}

	if (field[4])
	{
		printf("\n\tNIR: %i", p.NIR);
	}

	if (field[5])
	{
		printf("\n\tWave_Packet_Discriptor_Index: %i", p.Wave_Packet_Discriptor_Index);
		printf("\n\tByte_offset_to_waveform_data: %llu", p.Byte_offset_to_waveform_data);
		printf("\n\tWaveform_packet_size_in_bytes: %i", p.Waveform_packet_size_in_bytes);
		printf("\n\tReturn_Point_Waveform_Location: %f", p.Return_Point_Waveform_Location);

		printf("\n\tXt: %f", p.Xt);
		printf("\n\tYt: %f", p.Yt);
		printf("\n\tZt: %f", p.Zt);
	}

	printf("\n");
}

void LAS_io::Print_LAS_VLR(LAS_VLR vlr)
{
	printf("\n\tReserved: %i", vlr.vlr_header.Reserved);
	printf("\n\tUser_ID: %s", vlr.vlr_header.User_ID);
	printf("\n\tRecord_ID: %i", vlr.vlr_header.Record_ID);
	printf("\n\tRecord_Length_After_Header: %i", vlr.vlr_header.Record_Length_After_Header);
	printf("\n\tDescription: %s", vlr.vlr_header.Description);

//	printf("\n\tVariable Length Record (String Buffer): %s", vlr.vlr_buffer);
}

void LAS_io::Print_LAS_EVLR(LAS_EVLR evlr)
{
	Print_LAS_VLR(evlr);
}

bool LAS_io::Is_LAZ_File(char* filename)
{
	int len_name = (int)strlen(filename);
	char d = filename[len_name - 4];
	char l = filename[len_name - 3];
	char a = filename[len_name - 2];
	char z = filename[len_name - 1];

	if (d == '.'
		&& (l == 'l' || l == 'L')
		&& (a == 'a' || a == 'A')
		&& (z == 'z') || z == 'Z')
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool LAS_io::Is_LAS_File(char* filename)
{
	int len_name = (int)strlen(filename);
	char d = filename[len_name - 4];
	char l = filename[len_name - 3];
	char a = filename[len_name - 2];
	char s = filename[len_name - 1];

	if (d == '.'
		&& (l == 'l' || l == 'L')
		&& (a == 'a' || a == 'A')
		&& (s == 's') || s == 'S')
	{
		return true;
	}
	else
	{
		return false;
	}
}

void LAS_io::Data2Bytes(LAS_BYTE* bytes, void* data, unsigned short data_size, unsigned short n_data, unsigned long long &pos)
{
	memcpy(bytes + pos, data, data_size * n_data);
	pos += data_size * n_data;
}

void LAS_io::Bytes2Data(LAS_BYTE* bytes, void* data, unsigned short data_size, unsigned short n_data, unsigned long long &pos)
{
	memcpy(data, bytes + pos, data_size * n_data);
	pos += data_size * n_data;
}

double LAS_io::double_from_byte_8(LAS_BYTE_8 bytes)
{
	double val;
	memcpy(&val, &bytes, sizeof(LAS_BYTE_8));
	return val;
}

long long LAS_io::int64_from_byte_8(LAS_BYTE_8 bytes)
{
	long long val;
	memcpy(&val, &bytes, sizeof(LAS_BYTE_8));
	return val;
}

unsigned long long LAS_io::uint64_from_byte_8(LAS_BYTE_8 bytes)
{
	unsigned long long val;
	memcpy(&val, &bytes, sizeof(LAS_BYTE_8));
	return val;
}

