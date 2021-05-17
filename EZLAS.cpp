#include "EZLAS.h"

using namespace EZPC;

namespace EZLAS
{
	void Remove_Invalid_RGB(char** input_list, int i_input_first, int i_input_last)
	{
		Activate_Multi_Threading();
		printf("\nRemove Points without Color (i.e., RGB = [ 0, 0, 0 ]):\n");
		double t0 = omp_get_wtime();
		double t1, t2;

		int n_file = 0;
		char** file_list = file_list_from_file(input_list, i_input_first, i_input_last, n_file);

		for (int i_file = 0; i_file < n_file; i_file++)
		{
			printf("\n\tInput File #%i: %s", i_file, file_list[i_file]);		

			LAS_io las_in;

			printf("\n\t\tReading Points . . . ");
			t1 = omp_get_wtime();
			las_in.Read_LAS_File(file_list[i_file]);
			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);

			printf("\n\t\tFiltering Points . . . ");
			t1 = omp_get_wtime();
			bool* sample_flag = new bool[las_in.las_header.Number_of_point_records];
			
#pragma omp parallel for schedule (static)
			for (long long pid = 0; pid < (long long)las_in.las_header.Number_of_point_records; pid++)
			{
				if (las_in.las_points[pid].Red == 0
					&& las_in.las_points[pid].Green == 0
					&& las_in.las_points[pid].Blue == 0)
				{
					sample_flag[pid] = false;
				}
				else
				{
					sample_flag[pid] = true;
				}
			}
			las_in.Sample_LAS_Points(sample_flag);
			delete[] sample_flag;
			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);

			char* ext_out = match_extension("_NoBlkPts", file_list[i_file]);

			char* file_out = replace_extension(file_list[i_file], ext_out, true);
			printf("\n\t\tWriting Points . . . ");		
			t1 = omp_get_wtime();
			las_in.Write_LAS(file_out);
			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);

			printf("\n\t\tFlushing . . . ");
			t1 = omp_get_wtime();
			las_in.Delete_Data();
			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);
			printf("\n\tOutput File #%i: %s", i_file, file_out);
			delete[] ext_out;
			delete[] file_out;
			printf("\n");
		}
		t2 = omp_get_wtime();
		printf("\nAll Done! (%.2lf seconds)", t2 - t0);
		getchar();
	}

	void Merge_Split_SourceID(char** input_list, int i_input_first, int i_input_last)
	{
		Activate_Multi_Threading();
		printf("\nMerge Las File and Split by Point Source ID:\n");
		double t0 = omp_get_wtime();
		double t1, t2;

		int n_file = 0;
		char** file_list = file_list_from_file(input_list, i_input_first, i_input_last, n_file);
		c_sort_list(file_list, 0, n_file - 1);

		c_print_list("\n\t", file_list, n_file);
		printf("\n");

		char* file_name;
		if (n_file > 0)
		{
			file_name = add_extension(file_list[0], "_merge", true);
		}
		else
		{
			file_name = add_extension(file_list[0], "", true);
		}

		LAS_io las_in;

		printf("\n\tReading Points . . . ");
		t1 = omp_get_wtime();
		las_in.Read_LAS_Files(file_list, 0, n_file - 1);
		t2 = omp_get_wtime();
		printf("Done! (%.2lf seconds)", t2 - t1);

		printf("\n\tSorting Points . . . ");
		t1 = omp_get_wtime();
		las_in.Sort_by_Source_ID();		
		t2 = omp_get_wtime();
		printf("Done! (%.2lf seconds)", t2 - t1);
	
		printf("\n\tExtracting Points . . . ");
		t1 = omp_get_wtime();
		unsigned short min_id = las_in.las_points[0].Point_Source_ID;
		unsigned short max_id = las_in.las_points[las_in.las_header.Number_of_point_records - 1].Point_Source_ID;

		long long* id_start = new long long [max_id - min_id + 1];
		long long* id_end = new long long [max_id - min_id + 1];

#pragma omp parallel for schedule (static)
		for (int sid = min_id; sid <= max_id; sid++)
		{
			id_start[sid - min_id] = -1;
			id_end[sid - min_id] = -2;
		}
		id_start[0] = 0;
		id_end[max_id - min_id] = las_in.las_header.Number_of_point_records - 1;
#pragma omp parallel for schedule (guided)
		for (long long pid = 1; pid < (long long)las_in.las_header.Number_of_point_records; pid++)
		{
			LAS_Point p_cur = las_in.las_points[pid];
			LAS_Point p_pre = las_in.las_points[pid - 1];

			if (p_cur.Point_Source_ID > p_pre.Point_Source_ID)
			{
				id_end[p_pre.Point_Source_ID - min_id] = pid - 1;
				id_start[p_cur.Point_Source_ID - min_id] = pid;
			}
		}
		t2 = omp_get_wtime();
		printf("Done! (%.2lf seconds)", t2 - t1);


		printf("\n\tWriting Points . . . ");
		t1 = omp_get_wtime();

		for (int sid = min_id; sid <= max_id; sid++)
		{
			if (id_start[sid - min_id] <= id_end[sid - min_id])
			{
				LAS_io las_out;
				las_out.Copy_from_LAS_io(las_in, false, true, true);
				las_out.las_header.Number_of_point_records = id_end[sid - min_id] - id_start[sid - min_id] + 1;
				las_out.Initialize(true, false, false);

#pragma omp parallel for schedule (static)
				for (long long pid = id_start[sid - min_id]; pid <= id_end[sid - min_id]; pid++)
				{
					las_out.LAS2LAS_Point(las_in.las_points[pid], las_out.las_points[pid - id_start[sid - min_id]], las_out.size_extra_bytes(true));
				}

				char ext[20];
				sprintf(ext, "_srcid[%03i]%c", sid, '\0');

				char* f_out = add_extension(file_name, ext, true);
				las_out.Write_LAS(f_out);

				delete[] f_out;
				f_out = NULL;
				las_out.Delete_Data();
			}
		}
		t2 = omp_get_wtime();
		printf("Done! (%.2lf seconds)", t2 - t1);
		
		printf("\n\t\tFlushing . . . ");
		t1 = omp_get_wtime();
		las_in.Delete_Data();
		delete[] file_name;
		file_name = NULL;
		delete[] id_start;
		delete[] id_end;
		id_start = NULL;
		id_end = NULL;
		t2 = omp_get_wtime();
		printf("Done! (%.2lf seconds)", t2 - t1);

		t2 = omp_get_wtime();
		printf("\nAll Done! (%.2lf seconds)", t2 - t0);
		getchar();
	}

	void Merge_Split_Classification(char** input_list, int i_input_first, int i_input_last)
	{
		Activate_Multi_Threading();
		printf("\nMerge Las File and Split by Class:\n");
		double t0 = omp_get_wtime();
		double t1, t2;

		int n_file = 0;
		char** file_list = file_list_from_file(input_list, i_input_first, i_input_last, n_file);
		c_sort_list(file_list, 0, n_file - 1);

		c_print_list("\n\t", file_list, n_file);
		printf("\n");

		char* file_name;
		if (n_file > 0)
		{
			file_name = add_extension(file_list[0], "_merge", true);
		}
		else
		{
			file_name = add_extension(file_list[0], "", true);
		}

		LAS_io las_in;

		printf("\n\tReading Points . . . ");
		t1 = omp_get_wtime();
		las_in.Read_LAS_Files(file_list, 0, n_file - 1);
		t2 = omp_get_wtime();
		printf("Done! (%.2lf seconds)", t2 - t1);

		printf("\n\tSorting Points . . . ");
		t1 = omp_get_wtime();
		las_in.Sort_by_Classification();
		t2 = omp_get_wtime();
		printf("Done! (%.2lf seconds)", t2 - t1);

		printf("\n\tExtracting Points . . . ");
		t1 = omp_get_wtime();
		unsigned short min_id = las_in.las_points[0].Classification;
		unsigned short max_id = las_in.las_points[las_in.las_header.Number_of_point_records - 1].Classification;

		long long* id_start = new long long[max_id - min_id + 1];
		long long* id_end = new long long[max_id - min_id + 1];

#pragma omp parallel for schedule (static)
		for (int sid = min_id; sid <= max_id; sid++)
		{
			id_start[sid - min_id] = -1;
			id_end[sid - min_id] = -2;
		}
		id_start[0] = 0;
		id_end[max_id - min_id] = las_in.las_header.Number_of_point_records - 1;
#pragma omp parallel for schedule (guided)
		for (long long pid = 1; pid < (long long)las_in.las_header.Number_of_point_records; pid++)
		{
			LAS_Point p_cur = las_in.las_points[pid];
			LAS_Point p_pre = las_in.las_points[pid - 1];

			if (p_cur.Classification > p_pre.Classification)
			{
				id_end[p_pre.Classification - min_id] = pid - 1;
				id_start[p_cur.Classification - min_id] = pid;
			}
		}
		t2 = omp_get_wtime();
		printf("Done! (%.2lf seconds)", t2 - t1);


		printf("\n\tWriting Points . . . ");
		t1 = omp_get_wtime();

		for (int sid = min_id; sid <= max_id; sid++)
		{
			if (id_start[sid - min_id] <= id_end[sid - min_id])
			{
				LAS_io las_out;
				las_out.Copy_from_LAS_io(las_in, false, true, true);
				las_out.las_header.Number_of_point_records = id_end[sid - min_id] - id_start[sid - min_id] + 1;
				las_out.Initialize(true, false, false);

#pragma omp parallel for schedule (static)
				for (long long pid = id_start[sid - min_id]; pid <= id_end[sid - min_id]; pid++)
				{
					las_out.LAS2LAS_Point(las_in.las_points[pid], las_out.las_points[pid - id_start[sid - min_id]], las_out.size_extra_bytes(true));
				}

				char ext[20];
				sprintf(ext, "_class[%03i]%c", sid, '\0');

				char* f_out = add_extension(file_name, ext, true);
				las_out.Write_LAS(f_out);

				delete[] f_out;
				f_out = NULL;
				las_out.Delete_Data();
			}
		}
		t2 = omp_get_wtime();
		printf("Done! (%.2lf seconds)", t2 - t1);

		printf("\n\t\tFlushing . . . ");
		t1 = omp_get_wtime();
		las_in.Delete_Data();
		delete[] file_name;
		file_name = NULL;
		delete[] id_start;
		delete[] id_end;
		id_start = NULL;
		id_end = NULL;
		t2 = omp_get_wtime();
		printf("Done! (%.2lf seconds)", t2 - t1);

		t2 = omp_get_wtime();
		printf("\nAll Done! (%.2lf seconds)", t2 - t0);
		getchar();
	}

	void Meters_to_Feet(char** input_list, int i_input_first, int i_input_last)
	{
		Activate_Multi_Threading();
		printf("\nMeters to Feet:\n");
		double t0 = omp_get_wtime();
		double t1, t2;

		int n_file = 0;
		char** file_list = file_list_from_file(input_list, i_input_first, i_input_last, n_file);
		c_sort_list(file_list, 0, n_file - 1);

		for (int i_file = 0; i_file < n_file; i_file++)
		{
			printf("\n\tInput File #%i: %s", i_file, file_list[i_file]);

			LAS_io las_in;

			printf("\n\t\tReading Points . . . ");
			t1 = omp_get_wtime();
			las_in.Read_LAS_File(file_list[i_file]);
			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);

			printf("\n\t\tConverting . . . ");
			t1 = omp_get_wtime();

			las_in.las_header.Min_X = las_in.las_header.Min_X / 0.3048;
			las_in.las_header.Min_Y = las_in.las_header.Min_Y / 0.3048;
			las_in.las_header.Min_Z = las_in.las_header.Min_Z / 0.3048;
			las_in.las_header.Max_X = las_in.las_header.Max_X / 0.3048;
			las_in.las_header.Max_Y = las_in.las_header.Max_Y / 0.3048;
			las_in.las_header.Max_Z = las_in.las_header.Max_Z / 0.3048;

			las_in.las_header.X_offset = las_in.las_header.X_offset / 0.3048;
			las_in.las_header.Y_offset = las_in.las_header.Y_offset / 0.3048;
			las_in.las_header.Z_offset = las_in.las_header.Z_offset / 0.3048;

			las_in.las_header.X_scale_factor = las_in.las_header.X_scale_factor / 0.3048;
			las_in.las_header.Y_scale_factor = las_in.las_header.Y_scale_factor / 0.3048;
			las_in.las_header.Z_scale_factor = las_in.las_header.Z_scale_factor / 0.3048;

			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);

			char* file_out = replace_extension(file_list[i_file], "_feet.las", true);
			printf("\n\t\tWriting Points . . . ");
			t1 = omp_get_wtime();
			las_in.Write_LAS(file_out);
			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);

			printf("\n\t\tFlushing . . . ");
			t1 = omp_get_wtime();
			las_in.Delete_Data();
			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);
			printf("\n\tOutput File #%i: %s", i_file, file_out);
			delete[] file_out;
			printf("\n");
		}
		t2 = omp_get_wtime();
		printf("\nAll Done! (%.2lf seconds)", t2 - t0);
		getchar();
	}

	void Feet_to_Meters(char** input_list, int i_input_first, int i_input_last)
	{
		Activate_Multi_Threading();
		printf("\nMeters to Feet:\n");
		double t0 = omp_get_wtime();
		double t1, t2;

		int n_file = 0;
		char** file_list = file_list_from_file(input_list, i_input_first, i_input_last, n_file);
		c_sort_list(file_list, 0, n_file - 1);

		for (int i_file = 0; i_file < n_file; i_file++)
		{
			printf("\n\tInput File #%i: %s", i_file, file_list[i_file]);

			LAS_io las_in;

			printf("\n\t\tReading Points . . . ");
			t1 = omp_get_wtime();
			las_in.Read_LAS_File(file_list[i_file]);
			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);

			printf("\n\t\tConverting . . . ");
			t1 = omp_get_wtime();

			las_in.las_header.Min_X = las_in.las_header.Min_X * 0.3048;
			las_in.las_header.Min_Y = las_in.las_header.Min_Y * 0.3048;
			las_in.las_header.Min_Z = las_in.las_header.Min_Z * 0.3048;
			las_in.las_header.Max_X = las_in.las_header.Max_X * 0.3048;
			las_in.las_header.Max_Y = las_in.las_header.Max_Y * 0.3048;
			las_in.las_header.Max_Z = las_in.las_header.Max_Z * 0.3048;

			las_in.las_header.X_offset = las_in.las_header.X_offset * 0.3048;
			las_in.las_header.Y_offset = las_in.las_header.Y_offset * 0.3048;	
			las_in.las_header.Z_offset = las_in.las_header.Z_offset * 0.3048;

			las_in.las_header.X_scale_factor = las_in.las_header.X_scale_factor * 0.3048;
			las_in.las_header.Y_scale_factor = las_in.las_header.Y_scale_factor * 0.3048;
			las_in.las_header.Z_scale_factor = las_in.las_header.Z_scale_factor * 0.3048;

			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);

			char* file_out = replace_extension(file_list[i_file], "_meters.las", true);
			printf("\n\t\tWriting Points . . . ");
			t1 = omp_get_wtime();
			las_in.Write_LAS(file_out);
			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);

			printf("\n\t\tFlushing . . . ");
			t1 = omp_get_wtime();
			las_in.Delete_Data();
			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);
			printf("\n\tOutput File #%i: %s", i_file, file_out);
			delete[] file_out;
			printf("\n");
		}
		t2 = omp_get_wtime();
		printf("\nAll Done! (%.2lf seconds)", t2 - t0);
		getchar();
	}

	void LAS14_to_LAS13(char** input_list, int i_input_first, int i_input_last)
	{
		Activate_Multi_Threading();
		printf("\nMeters to Feet:\n");
		double t0 = omp_get_wtime();
		double t1, t2;

		int n_file = 0;
		char** file_list = file_list_from_file(input_list, i_input_first, i_input_last, n_file);
		c_sort_list(file_list, 0, n_file - 1);

		for (int i_file = 0; i_file < n_file; i_file++)
		{
			printf("\n\tInput File #%i: %s", i_file, file_list[i_file]);

			LAS_io las_in;

			printf("\n\t\tReading Points . . . ");
			t1 = omp_get_wtime();
			las_in.Read_LAS_File(file_list[i_file]);
			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);

			printf("\n\t\tConverting . . . ");
			t1 = omp_get_wtime();
			las_in.Show_LAS_Header();
			las_in.Show_VLR(0);
			las_in.Show_VLR(1);
			getchar();
			las_in.las_header.Version_Minor = 3;
		//	las_in.Autofill_Header();

			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);

			char* file_out = replace_extension(file_list[i_file], "_las13.las", true);
			printf("\n\t\tWriting Points . . . ");
			t1 = omp_get_wtime();
			las_in.Write_LAS(file_out);
			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);

			printf("\n\t\tFlushing . . . ");
			t1 = omp_get_wtime();
			las_in.Delete_Data();
			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);
			printf("\n\tOutput File #%i: %s", i_file, file_out);
			delete[] file_out;
			printf("\n");
		}
		t2 = omp_get_wtime();
		printf("\nAll Done! (%.2lf seconds)", t2 - t0);
		getchar();
	}

	void LAS_to_LAZ(char** input_list, int i_input_first, int i_input_last)
	{
		Activate_Multi_Threading();
		printf("\nLAS to LAZ:\n");
		double t0 = omp_get_wtime();
		double t1, t2;

		int n_file = 0;
		char** file_list = file_list_from_file(input_list, i_input_first, i_input_last, n_file);
		c_sort_list(file_list, 0, n_file - 1);

		for (int i_file = 0; i_file < n_file; i_file++)
		{
			printf("\n\tInput File #%i: %s", i_file, file_list[i_file]);

			if (is_extension_matched(file_list[i_file], ".laz"))
			{
				printf("\n\tSkipped!");
				continue;
			}

			LAS_io las_in;

			printf("\n\t\tReading Points . . . ");
			t1 = omp_get_wtime();
			las_in.Read_LAS_File(file_list[i_file]);
			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);

			char* file_out = replace_extension(file_list[i_file], ".laz", true);
			printf("\n\t\tWriting Points . . . ");
			t1 = omp_get_wtime();
			las_in.Write_LAS(file_out);
			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);

			printf("\n\t\tFlushing . . . ");
			t1 = omp_get_wtime();
			las_in.Delete_Data();
			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);
			printf("\n\tOutput File #%i: %s", i_file, file_out);
			delete[] file_out;
			printf("\n");
		}
		t2 = omp_get_wtime();
		printf("\nAll Done! (%.2lf seconds)", t2 - t0);
		getchar();
	}

	void LAZ_to_LAS(char** input_list, int i_input_first, int i_input_last)
	{
		Activate_Multi_Threading();
		printf("\nLAZ to LAS:\n");
		double t0 = omp_get_wtime();
		double t1, t2;

		int n_file = 0;
		char** file_list = file_list_from_file(input_list, i_input_first, i_input_last, n_file);
		c_sort_list(file_list, 0, n_file - 1);

		for (int i_file = 0; i_file < n_file; i_file++)
		{
			printf("\n\tInput File #%i: %s", i_file, file_list[i_file]);

			if (is_extension_matched(file_list[i_file], ".las"))
			{
				printf("\n\tSkipped!");
				continue;
			}

			LAS_io las_in;

			printf("\n\t\tReading Points . . . ");
			t1 = omp_get_wtime();
			las_in.Read_LAS_File(file_list[i_file]);
			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);

			char* file_out = replace_extension(file_list[i_file], ".las", true);
			printf("\n\t\tWriting Points . . . ");
			t1 = omp_get_wtime();
			las_in.Write_LAS(file_out);
			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);

			printf("\n\t\tFlushing . . . ");
			t1 = omp_get_wtime();
			las_in.Delete_Data();
			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);
			printf("\n\tOutput File #%i: %s", i_file, file_out);
			delete[] file_out;
			printf("\n");
		}
		t2 = omp_get_wtime();
		printf("\nAll Done! (%.2lf seconds)", t2 - t0);
		getchar();
	}

	void Merge(char** input_list, int i_input_first, int i_input_last)
	{
		Activate_Multi_Threading();
		printf("\nMerge Las/Laz Files:\n");
		double t0 = omp_get_wtime();
		double t1, t2;

		int n_file = 0;
		char** file_list = file_list_from_file(input_list, i_input_first, i_input_last, n_file);
		c_sort_list(file_list, 0, n_file - 1);

		c_print_list("\n\t", file_list, n_file);
		printf("\n");

		char* file_name = add_extension(file_list[0], "_merge", true);

		LAS_io las_in;

		printf("\n\tReading Points . . . ");
		t1 = omp_get_wtime();
		las_in.Read_LAS_Files(file_list, 0, n_file - 1);
		t2 = omp_get_wtime();
		printf("Done! (%.2lf seconds)", t2 - t1);

		printf("\n\t\tWriting Points . . . ");
		t1 = omp_get_wtime();
		las_in.Write_LAS(file_name);
		t2 = omp_get_wtime();
		printf("Done! (%.2lf seconds)", t2 - t1);

		printf("\n\t\tFlushing . . . ");
		t1 = omp_get_wtime();
		las_in.Delete_Data();
		t2 = omp_get_wtime();
		printf("Done! (%.2lf seconds)", t2 - t1);
		printf("\n\tOutput File: %s", file_name);
		delete[] file_name;
		printf("\n");
		t2 = omp_get_wtime();
		printf("\nAll Done! (%.2lf seconds)", t2 - t0);
		getchar();
	}

	void Shift_Classification(char** input_list, int i_input_first, int i_input_last)
	{
		Activate_Multi_Threading();
		printf("\nShift Classification:\n");
		printf("\nOffset: ");
		int offset = 0;
		scanf("%i", &offset);
		getchar();
		double t0 = omp_get_wtime();
		double t1, t2;

		int n_file = 0;
		char** file_list = file_list_from_file(input_list, i_input_first, i_input_last, n_file);	
		c_sort_list(file_list, 0, n_file - 1);


		char ext[EZPC_PATH_BUFFER];
		sprintf(ext, "_class%c%i", offset < 0 ? '-' : '+', offset);

		for (int i_file = 0; i_file < n_file; i_file++)
		{
			printf("\n\tInput File #%i: %s", i_file, file_list[i_file]);

			LAS_io las_in;

			printf("\n\t\tReading Points . . . ");
			t1 = omp_get_wtime();
			las_in.Read_LAS_File(file_list[i_file]);
			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);

			printf("\n\t\tApplying Offset . . . ");
			t1 = omp_get_wtime();

#pragma omp parallel for schedule (static)
			for (long long i = 0; i < (long long)las_in.las_header.Number_of_point_records; i++)
			{
				las_in.las_points[i].Classification += offset;
			}

			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);

			char* file_out = add_extension(file_list[i_file], ext, true);
			printf("\n\t\tWriting Points . . . ");
			t1 = omp_get_wtime();
			las_in.Write_LAS(file_out);
			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);
	
			printf("\n\t\tFlushing . . . ");
			t1 = omp_get_wtime();
			las_in.Delete_Data();
			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);
			printf("\n\tOutput File #%i: %s", i_file, file_out);
			delete[] file_out;
			printf("\n");
		}

		t2 = omp_get_wtime();
		printf("\nAll Done! (%.2lf seconds)", t2 - t0);
		getchar();
	}

	void Extra_Bytes_Writing_Example(char** input_list, int i_input_first, int i_input_last)
	{
		Activate_Multi_Threading();
		printf("\nExtra Bytes Writing:\n");
		double t0 = omp_get_wtime();
		double t1, t2;

		int n_file = 0;
		char** file_list = file_list_from_file(input_list, i_input_first, i_input_last, n_file);
		c_sort_list(file_list, 0, n_file - 1);

		for (int i_file = 0; i_file < n_file; i_file++)
		{
			printf("\n\tInput File #%i: %s", i_file, file_list[i_file]);

			LAS_io las_in;

			printf("\n\t\tReading Points . . . ");
			t1 = omp_get_wtime();
			las_in.Read_LAS_File(file_list[i_file]);
			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);

			printf("\n\t\tAdding Extra Bytes . . . ");
			t1 = omp_get_wtime();

			unsigned short n_eb = 2;
			LAS_Extra_Bytes_Info* eb_info = new LAS_Extra_Bytes_Info[n_eb];

			// extra bytes #0
			eb_info[0] = las_in.create_extra_bytes_info(EB_Type_Double);
			// in thoery, all the following should be optional. CloudCompare may want you to have a name and description. 
			las_in.Set_Extra_Bytes_No_Data(eb_info[0], (double)-1);
			las_in.Set_Extra_Bytes_Min(eb_info[0], (double)-1);
			las_in.Set_Extra_Bytes_Max(eb_info[0], (double)1);

			las_in.Set_Extra_Bytes_Name(eb_info[0], "Random_Number");
			las_in.Set_Extra_Bytes_Description(eb_info[0], "Random_Double_Extra_Bytes");

			////////////////////////////////////////////////////
			// extra bytes #1
			eb_info[1] = las_in.create_extra_bytes_info(EB_Type_Long_3);

			// in thoery, all the following should be optional. CloudCompare may want you to have a name and description. 
			las_in.Set_Extra_Bytes_Scale(eb_info[1], las_in.las_header.X_scale_factor, las_in.las_header.Y_scale_factor, las_in.las_header.Z_scale_factor);
			las_in.Set_Extra_Bytes_Offset(eb_info[1], las_in.las_header.X_offset, las_in.las_header.Y_offset, las_in.las_header.Z_offset);

			las_in.Set_Extra_Bytes_Name(eb_info[1], "Coordinates");
			las_in.Set_Extra_Bytes_Description(eb_info[1], "Coordinates_Extra_Bytes");

			/////////////////////////////////////////////////////
			las_in.Add_Extra_Bytes_VLR(eb_info, n_eb); // know the start of the extra bytes index. 

#pragma omp parallel for schedule (static)
			for (long long pid = 0; pid < (long long)las_in.las_header.Number_of_point_records; pid++)
			{
				// extra bytes #0
				if (pid % 100 == 0)
				{
					las_in.Set_Extra_Bytes_Value(pid, eb_info[0], las_in.double_from_byte_8(eb_info[0].no_data[0]));  // set no data every 1000 point
				}
				else
				{
					double random = rand() % 100 / 100.;
					las_in.Set_Extra_Bytes_Value(pid, eb_info[0], &random);
				}

				// extra bytes #1
				int record_x = las_in.las_points[pid].X;
				int record_y = las_in.las_points[pid].Y;
				int record_z = las_in.las_points[pid].Z;
				las_in.Set_Extra_Bytes_Value(pid, eb_info[1], &record_x, &record_y, &record_z); // warning, it is extremely hard to handle the casting automatically. So at least make sure not to mix integer and floating points.
			}

			delete[] eb_info;

			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);

			
			las_in.Show_LAS_Header();
			las_in.Show_Extra_Bytes_VLR();
			
			
			char* file_out = replace_extension(file_list[i_file], "_eb.laz", true);
			printf("\n\t\tWriting Points . . . ");
			t1 = omp_get_wtime();
			las_in.Write_LAS(file_out);
			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);

			printf("\n\t\tFlushing . . . ");
			t1 = omp_get_wtime();
			las_in.Delete_Data();
			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);
			printf("\n\tOutput File #%i: %s", i_file, file_out);
			delete[] file_out;
			printf("\n");
		}
		t2 = omp_get_wtime();
		printf("\nAll Done! (%.2lf seconds)", t2 - t0);
		getchar();
	}

//	void Extra_Bytes_Reading_Example(char** file_list, int i_file_start, int i_file_end)  // reading extra bytes and remove points with no_data on any of fields.
//	{
//		Multi_Threading();
//		printf("\nExtra Bytes Reading:\n");
//		double t0 = omp_get_wtime();
//		double t1, t2;
//
//		for (int i_file = i_file_start; i_file <= i_file_end; i_file++)
//		{
//			printf("\n\tInput File #%i: %s", i_file, file_list[i_file]);
//
//			LAS_io las_in;
//
//			printf("\n\t\tReading Points . . . ");
//			t1 = omp_get_wtime();
//			las_in.Read_LAS_File(file_list[i_file]);
//			t2 = omp_get_wtime();
//			printf("Done! (%.2lf seconds)", t2 - t1);
//
//
//			las_in.Show_LAS_Header();
//			las_in.Show_Extra_Bytes_VLR();
//
//
//			printf("\n\t\tReading Extra Bytes . . . ");
//			t1 = omp_get_wtime();
//	
//			unsigned short n_eb = 0;
//			Extra_Bytes_Info* eb_info;
//			las_in.Get_All_Extra_Bytes_Info(eb_info, n_eb);
//
//			if (n_eb > 0)
//			{
//				las_in.Show_Extra_Bytes_VLR();
//
//				bool *flag_sample = new bool[las_in.las_header.Number_of_point_records];
//#pragma omp parallel for schedule (static)
//				for (long long pid = 0; pid < (long long)las_in.las_header.Number_of_point_records; pid++)
//				{
//					flag_sample[pid] = true;
//				}
//
//#pragma omp parallel for schedule (static)
//				for (long long pid = 0; pid < (long long)las_in.las_header.Number_of_point_records; pid++)
//				{
//					if (pid % (las_in.las_header.Number_of_point_records / 10) == 0)
//					{
//#pragma omp critical
//						printf("\npid: %lli", pid);
//						for (int i_eb = 0; i_eb < n_eb; i_eb++)
//						{
//							LAS_BYTE_8 eb_x;
//							LAS_BYTE_8 eb_y;
//							LAS_BYTE_8 eb_z;
//							las_in.Get_Extra_Bytes_Value(pid, eb_info[i_eb], eb_x, eb_y, eb_z);
//
//							if (eb_x == eb_info[i_eb].no_data[0] && eb_y == eb_info[i_eb].no_data[1] && eb_z == eb_info[i_eb].no_data[2])
//							{
//								flag_sample[pid] = false;
//							}
//							else
//							{
//
//								{
//									printf("\nextra bytes: ");
//									if (eb_info[i_eb].floating_or_integer)
//									{
//										printf("%lf ", double_from_byte_8(eb_x));
//										if (eb_info[i_eb].n_values > 1)
//										{
//											printf("%lf ", double_from_byte_8(eb_y));
//											if (eb_info[i_eb].n_values > 2)
//											{
//												printf("%lf ", double_from_byte_8(eb_z));
//											}
//										}
//									}
//									else
//									{
//										if (eb_info[i_eb].signed_or_unsigned)
//										{
//											printf("%lli ", int64_from_byte_8(eb_x));
//											if (eb_info[i_eb].n_values > 1)
//											{
//												printf("%lli ", int64_from_byte_8(eb_y));
//												if (eb_info[i_eb].n_values > 2)
//												{
//													printf("%lli ", int64_from_byte_8(eb_z));
//												}
//											}
//										}
//										else
//										{
//											printf("%lli ", uint64_from_byte_8(eb_x));
//											if (eb_info[i_eb].n_values > 1)
//											{
//												printf("%lli ", uint64_from_byte_8(eb_y));
//												if (eb_info[i_eb].n_values > 2)
//												{
//													printf("%lli ", uint64_from_byte_8(eb_z));
//												}
//											}
//										}
//									}
//
//									double x, y, z;
//
//									if (eb_info[i_eb].floating_or_integer)
//									{
//										if (eb_info[i_eb].options_scale_bit)
//										{
//											x = double_from_byte_8(eb_x) * eb_info[i_eb].scale[0];
//											y = double_from_byte_8(eb_y) * eb_info[i_eb].scale[1];
//											z = double_from_byte_8(eb_z) * eb_info[i_eb].scale[2];
//										}
//										else
//										{
//											x = double_from_byte_8(eb_x);
//											y = double_from_byte_8(eb_y);
//											z = double_from_byte_8(eb_z);
//										}
//									}
//									else
//									{
//										if (eb_info[i_eb].options_scale_bit)
//										{
//											x = int64_from_byte_8(eb_x) * eb_info[i_eb].scale[0];
//											y = int64_from_byte_8(eb_y) * eb_info[i_eb].scale[1];
//											z = int64_from_byte_8(eb_z) * eb_info[i_eb].scale[2];
//										}
//										else
//										{
//											x = (double)int64_from_byte_8(eb_x);
//											y = (double)int64_from_byte_8(eb_y);
//											z = (double)int64_from_byte_8(eb_z);
//										}
//									}
//
//									if (eb_info[i_eb].options_offset_bit)
//									{
//										x += eb_info[i_eb].offset[0];
//										y += eb_info[i_eb].offset[1];
//										z += eb_info[i_eb].offset[2];
//									}
//
//									if (eb_info[i_eb].options_scale_bit || eb_info[i_eb].options_offset_bit)
//									{
//										printf("\t-> %lf ", x);
//										if (eb_info[i_eb].n_values > 1)
//										{
//											printf("%lf ", y);
//											if (eb_info[i_eb].n_values > 2)
//											{
//												printf("%lf ", z);
//											}
//										}
//									}
//								}
//							}
//						}
//						printf("\n");
//					}
//				}
//
//				las_in.Sample_LAS_Points(flag_sample);  //sample the points
//
//				delete[] eb_info;
//				delete[] flag_sample;
//			}
//
//			t2 = omp_get_wtime();
//			printf("Done! (%.2lf seconds)", t2 - t1);
//
//			char* file_out = Replace_Extension(file_list[i_file], "_fltr.las", true, false);
//			printf("\n\t\tWriting Points . . . ");
//			t1 = omp_get_wtime();
//			las_in.Write_LAS(file_out);
//			t2 = omp_get_wtime();
//			printf("Done! (%.2lf seconds)", t2 - t1);
//
//			printf("\n\t\tFlushing . . . ");
//			t1 = omp_get_wtime();
//			las_in.Delete_Data();
//			t2 = omp_get_wtime();
//			printf("Done! (%.2lf seconds)", t2 - t1);
//			printf("\n\tOutput File #%i: %s", i_file, file_out);
//			delete[] file_out;
//			printf("\n");
//		}
//		t2 = omp_get_wtime();
//		printf("\nAll Done! (%.2lf seconds)", t2 - t0);
//		getchar();
//	}

	void Remove_ExtraBytes(char** input_list, int i_input_first, int i_input_last)
	{
		Activate_Multi_Threading();
		printf("\nRemove ExtraBytes:\n");

		double t0 = omp_get_wtime();
		double t1, t2;

		int n_file = 0;
		char** file_list = file_list_from_file(input_list, i_input_first, i_input_last, n_file);
		c_sort_list(file_list, 0, n_file - 1);


		for (int i_file = 0; i_file < n_file; i_file++)
		{
			printf("\n\tInput File #%i: %s", i_file, file_list[i_file]);

			LAS_io las_in;

			printf("\n\t\tReading Points . . . ");
			t1 = omp_get_wtime();
			las_in.Read_LAS_File(file_list[i_file]);
			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);

			printf("\n\t\tRemoving ExtraBytes . . . ");
			t1 = omp_get_wtime();

			las_in.Clear_Extra_Bytes();

			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);

			char* file_out = add_extension(file_list[i_file], "_no_extrabytes", true);
			printf("\n\t\tWriting Points . . . ");
			t1 = omp_get_wtime();
			las_in.Write_LAS(file_out);
			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);

			printf("\n\t\tFlushing . . . ");
			t1 = omp_get_wtime();
			las_in.Delete_Data();
			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);
			printf("\n\tOutput File #%i: %s", i_file, file_out);
			delete[] file_out;
			printf("\n");
		}

		t2 = omp_get_wtime();
		printf("\nAll Done! (%.2lf seconds)", t2 - t0);
		getchar();
	}

	void Remove_Offset(char** input_list, int i_input_first, int i_input_last)
	{
		Activate_Multi_Threading();
		printf("\nRemove Offset:\n");
		double t0 = omp_get_wtime();
		double t1, t2;

		int n_file = 0;
		char** file_list = file_list_from_file(input_list, i_input_first, i_input_last, n_file);
		c_sort_list(file_list, 0, n_file - 1);

		char ext[EZPC_PATH_BUFFER];
		sprintf(ext, "_no_offset");

		LAS_io las_tmp;
		las_tmp.las_header = las_tmp.Merge_LAS_Headers(file_list, 0, n_file - 1);

		double offset_x = las_tmp.las_header.X_offset;
		double offset_y = las_tmp.las_header.Y_offset;
		double offset_z = las_tmp.las_header.Z_offset;

		printf("\nOffset_X: %lf", offset_x);
		printf("\nOffset_Y: %lf", offset_y);
		printf("\nOffset_Z: %lf", offset_z);
		printf("\n");

		las_tmp.Delete_Data();

		for (int i_file = 0; i_file < n_file; i_file++)
		{
			printf("\n\tInput File #%i: %s", i_file, file_list[i_file]);

			LAS_io las_in;

			printf("\n\t\tReading Points . . . ");
			t1 = omp_get_wtime();
			las_in.Read_LAS_File(file_list[i_file]);
			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);

			printf("\n\t\tRemoving Offset . . . ");
			t1 = omp_get_wtime();
			las_in.las_header.X_offset -= offset_x;
			las_in.las_header.Y_offset -= offset_y;
			las_in.las_header.Z_offset -= offset_z;

			las_in.las_header.Min_X -= offset_x;
			las_in.las_header.Min_Y -= offset_y;
			las_in.las_header.Min_Z -= offset_z;

			las_in.las_header.Max_X -= offset_x;
			las_in.las_header.Max_Y -= offset_y;
			las_in.las_header.Max_Z -= offset_z;

			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);

			char* file_out = add_extension(file_list[i_file], ext, true);
			printf("\n\t\tWriting Points . . . ");
			t1 = omp_get_wtime();
			las_in.Write_LAS(file_out);
			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);

			printf("\n\t\tFlushing . . . ");
			t1 = omp_get_wtime();
			las_in.Delete_Data();
			t2 = omp_get_wtime();
			printf("Done! (%.2lf seconds)", t2 - t1);
			printf("\n\tOutput File #%i: %s", i_file, file_out);
			delete[] file_out;
			printf("\n");
		}

		t2 = omp_get_wtime();
		printf("\nAll Done! (%.2lf seconds)", t2 - t0);
		getchar();
	}

	void Show_Header(char** input_list, int i_input_first, int i_input_last)
	{
		Activate_Multi_Threading();
		printf("\nShow Headers:\n");

		int n_file = 0;
		char** file_list = file_list_from_file(input_list, i_input_first, i_input_last, n_file);
		c_sort_list(file_list, 0, n_file - 1);

		for (int i_file = 0; i_file < n_file; i_file++)
		{
			printf("\n\tInput File #%i: %s", i_file, file_list[i_file]);

			LAS_io las_in;
			las_in.las_header = las_in.Read_LAS_Header(file_list[i_file], false);
			las_in.Show_LAS_Header();
			getchar();
		}

		printf("\nAll Done!");
		getchar();
	}
};