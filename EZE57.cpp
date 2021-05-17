#include "EZE57.h"

using namespace EZPC;

namespace EZE57
{
	void E57_to_LAS(char** input_list, int i_input_first, int i_input_last, bool LAS_or_LAZ, bool report)
	{
		printf("\n===== Convert E57 to LAS/LAZ! =====\n");

		int n_file = 0;
		char** file_list = file_list_from_file(input_list, i_input_first, i_input_last, n_file);
		c_sort_list(file_list, 0, n_file - 1);

		printf("\nInput Files: ");
		for (int i_file = 0; i_file < n_file; i_file++)
		{
			printf("\n\t%s", file_list[i_file]);
		}
		printf("\n");

		printf("\nProcess Started!\n");

		Activate_Multi_Threading();

		double time0, time1;
		time0 = omp_get_wtime();

		//read each file
		FILE* fp_report = NULL;
		if (report)
		{
			fp_report = fopen("[EZE57]_e57_to_las_report.csv", "a+");
			if (!fp_report)
			{
				printf("\n\nCannot open/create \"[EZE57]_e57_to_las_report.csv\"!!!");
				getchar();
			}
			fseek(fp_report, 0, SEEK_END);
			if (ftell(fp_report) == 0)
			{
				fprintf(fp_report, "Processing_Date, Processing_Time, ");
				fprintf(fp_report, "E57_Path, E57_File, LAS_File, ");
				fprintf(fp_report, "Scan_Name, Scan_Description, Scan_Vendor, Scan_Model, ");
				fprintf(fp_report, "Scan_Date, Scan_Time, ");
				fprintf(fp_report, "Translation_X, Translation_Y, Translation_Z, ");
				fprintf(fp_report, "Quaternion_x, Quaternion_y, Quaternion_z, Quaternion_w, ");
				fprintf(fp_report, "Row_Count, Column_Count, Point_Count, ");
				fprintf(fp_report, "Intensity_Min, Intensity_Max, RGB_Field, Time_Field");
				fprintf(fp_report, "\n");
			}
		}

		for (int i_file = 0; i_file < n_file; i_file++)
		{
			__int64 n_scn = e57_scan_count(file_list[i_file]);
			printf("\nReading File (%i/%i) %s", i_file - i_file + 1, n_file, file_list[i_file]);
			for (__int64 i_s = 0; i_s < n_scn; i_s++)
			{
				E57_io e57_io;
				e57_io.Add_E57_Scan(file_list[i_file], i_s);
				printf("\n\tScan (%lli/%lli) %s", i_s + 1, n_scn, e57_io.e57_scan[0].e57_header.name);
				E57_Header e57_h = e57_io.e57_scan[0].e57_header;

				LAS_io las_io;

				double scale_factor = e57_h.distance_scale_factor == 0 ? LASIO_DEFAULT_SCALE : e57_h.distance_scale_factor;
				las_io.Initialize(scale_factor, scale_factor, scale_factor, e57_h.translation_x, e57_h.translation_y, e57_h.translation_z, e57_h.point_count, 0, 0);

				//copy points over.
				__int64 pid_las = 0;
				bool rgb_valid = e57_h.red_field && e57_h.green_field && e57_h.blue_field;
				bool flag_grey = true;
				double intensity_range = e57_h.intensity_maximum_limit - e57_h.intensity_minimum_limit;

				bool index_valid = e57_h.row_index_field || e57_h.column_index_field;
				
				for (__int64 i_p = 0; i_p < e57_h.point_count; i_p++)
				{
					E57_Point local_point = e57_io.e57_scan_point(0, i_p, true);

					if ((!e57_h.xyz_invalid_field || !local_point.invalid_xyz)
						&& (!index_valid || (local_point.x != 0 || local_point.y != 0 || local_point.z != 0)))
					{
						las_io.las_points[pid_las].X = las_io.x_coord2record(e57_io.e57_scan_point(0, i_p, false).x);
						las_io.las_points[pid_las].Y = las_io.y_coord2record(e57_io.e57_scan_point(0, i_p, false).y);
						las_io.las_points[pid_las].Z = las_io.z_coord2record(e57_io.e57_scan_point(0, i_p, false).z);


						if (!e57_h.intensity_invalid_field || !local_point.invalid_intensity)
						{
							las_io.las_points[pid_las].Intensity = unsigned short(local_point.intensity / intensity_range * USHRT_MAX);
						}
						else
						{
							las_io.las_points[pid_las].Intensity = 0;
						}

						if (!e57_h.rgb_invalid_field || !local_point.invalid_rgb)
						{
							las_io.las_points[pid_las].Red = local_point.red;
							las_io.las_points[pid_las].Green = local_point.green;
							las_io.las_points[pid_las].Blue = local_point.blue;
							if (!(local_point.red == local_point.green &&  local_point.green == local_point.blue))
							{
								flag_grey = false;
							}
						}
						else
						{
							las_io.las_points[pid_las].Red = 0;
							las_io.las_points[pid_las].Green = 0;
							las_io.las_points[pid_las].Blue = 0;
						}

						if (!e57_h.time_stamp_invalid_field || !local_point.invalid_time_stamp)
						{
							las_io.las_points[pid_las].GPS_Time = local_point.time_stamp;
						}
						else
						{
							las_io.las_points[pid_las].GPS_Time = 0;
						}
						pid_las++;
					}
				}
				char* e57_dir = get_directory(file_list[i_file]);
				char las_ext[E57_IO_MAX_STRING_SIZE];
				sprintf(las_ext, "_%03lli%s%c\n", e57_io.e57_scan[0].scan_index + 1, LAS_or_LAZ ? ".las" : ".laz", '\0');
				char* las_file =  replace_extension(file_list[i_file], las_ext, false);
				char* e57_file = remove_directory(file_list[i_file]);

				if (report)
				{
					time_t rawtime; 
					struct tm * ptm;
					time(&rawtime);
					ptm = gmtime(&rawtime);
					
					fprintf(fp_report, "%04i%c%02i%c%02i, %02i%c%02i%c%02i, ", ptm->tm_year + 1900, '-', ptm->tm_mon + 1, '-', ptm->tm_mday, ptm->tm_hour, ':', ptm->tm_min, ':', ptm->tm_sec);
					fprintf(fp_report, "\"%s\", \"%s\", \"%s\",", e57_dir, e57_file, las_file);
					fprintf(fp_report, "\"%s\", \"%s\", \"%s\", \"%s\",", e57_h.name, e57_h.description, e57_h.vendor, e57_h.model);

					int year, month, day, hour, minute; 
					float seconds;
					e57::DateTime datetime;
					datetime.dateTimeValue = e57_h.date_time_start;
					datetime.GetUTCDateTime(year, month, day, hour, minute, seconds);
					fprintf(fp_report, "%04i%c%02i%c%02i, %02i%c%02i%c%02i, ", year, '-', month, '-', day, hour, ':', minute, ':', (int)seconds);
					fprintf(fp_report, "%lf, %lf, %lf, ", e57_h.translation_x, e57_h.translation_y, e57_h.translation_z);
					fprintf(fp_report, "%lf, %lf, %lf, %lf, ", e57_h.quaternion_x, e57_h.quaternion_y, e57_h.quaternion_z, e57_h.quaternion_w);
					fprintf(fp_report, "%lli, %lli, %lli, ", e57_h.row_maximum, e57_h.column_maximum - e57_h.column_minimum, pid_las);
					fprintf(fp_report, "%lf, %lf, ", e57_h.intensity_minimum_limit, e57_h.intensity_maximum_limit);
					fprintf(fp_report, "%s, %s", (rgb_valid && !flag_grey) ? "Yes" : "No", e57_h.time_stamp_field ? "Yes" : "No");
					fprintf(fp_report, "\n");
				}

				unsigned char format_id = e57_h.time_stamp_field ? 3 : 2;

				las_io.Autofill_Header(1, 3, format_id, pid_las, 0, 0);
				las_io.Write_LAS(las_file);

				las_io.Delete_Data();
				delete[] e57_dir;
				delete[] las_file;
				delete[] e57_file;

				e57_io.Delete();
			}
			printf("\n");
		}

		if (report)
		{
			fclose(fp_report);
		}

		time1 = omp_get_wtime();
		printf("\n\nAll Done! (Total Time: %.2lf seconds)", time1 - time0);
		getchar();
	}
	
	void Extract_E57_Info(char** input_list, int i_input_first, int i_input_last)
	{
		printf("\n===== Extract E57 Info! =====\n");

		int n_file = 0;
		char** file_list = file_list_from_file(input_list, i_input_first, i_input_last, n_file);
		c_sort_list(file_list, 0, n_file - 1);

		printf("\nInput Files: ");
		for (int i_file = 0; i_file < n_file; i_file++)
		{
			printf("\n\t%s", file_list[i_file]);
		}
		printf("\n");

		printf("\nProcess Started!\n");

		Activate_Multi_Threading();

		double time0, time1;
		time0 = omp_get_wtime();

		for (int i_file = 0; i_file < n_file; i_file++)
		{
			char* fname_report = EZPC::replace_extension(file_list[i_file], ".txt", true);
			FILE* fp = fopen(fname_report, "wt");
			fprintf(fp, "\"%s\"\n", file_list[i_file]);
			__int64 n_scn = 0;

			printf("\nReading File (%i/%i) %s", i_file + 1, n_file, file_list[i_file]);
			E57_Header* hdr = Get_E57_Header(file_list[i_file], n_scn);
	
			for (__int64 i_s = 0; i_s < n_scn; i_s++)
			{
				fprintf(fp, "\nScan #%lli: ", i_s);
				fprintf(fp, "\n\tScan Info: ");
				fprintf(fp, "\n\t\tName: %s", hdr[i_s].name);
				fprintf(fp, "\n\t\tGUID: %s", hdr[i_s].guid);

				fprintf(fp, "\n\t\tOriginal GUIDs: ");
				for (int i = 0; i < hdr[i_s].original_guid_count; i++)
				{
					fprintf(fp, "\n\t\t\t#%i: %s", i, hdr[i_s].original_guid[i]);
				}
				fprintf(fp, "\n\t\tDescription: %s", hdr[i_s].description);

				fprintf(fp, "\n\tSensor Info: ");
				fprintf(fp, "\n\t\tVendor: %s", hdr[i_s].vendor);
				fprintf(fp, "\n\t\tModel: %s", hdr[i_s].model);
				fprintf(fp, "\n\t\tSerial Number: %s", hdr[i_s].serial_number);
				fprintf(fp, "\n\t\tHardware Version: %s", hdr[i_s].hardware_version);
				fprintf(fp, "\n\t\tSoftware Version: %s", hdr[i_s].software_version);
				fprintf(fp, "\n\t\tFirmware Version: %s", hdr[i_s].firmware_version);

				fprintf(fp, "\n\tAcquisition Info: ");
				fprintf(fp, "\n\t\tTemperautre: %f", hdr[i_s].temperature);
				fprintf(fp, "\n\t\tRelative Humidity: %f", hdr[i_s].relative_humidity);
				fprintf(fp, "\n\t\tAtmospheric Pressure: %f", hdr[i_s].atmospheric_pressure);

				e57::DateTime datetime;
				int year, month, day, hour, minute;
				float seconds;

				datetime.dateTimeValue = hdr[i_s].date_time_start;
				datetime.GetUTCDateTime(year, month, day, hour, minute, seconds);
				fprintf(fp, "\n\t\tDate Time Start: %04i%c%02i%c%02i, %02i%c%02i%c%02i", year, '-', month, '-', day, hour, ':', minute, ':', (int)seconds);
				fprintf(fp, "\n\t\tAtomic Clock Referenced: %s", hdr[i_s].data_time_start_isAtomicClockReferenced ? "Yes" : "No");

				datetime.dateTimeValue = hdr[i_s].date_time_end;
				datetime.GetUTCDateTime(year, month, day, hour, minute, seconds);
				fprintf(fp, "\n\t\tDate Time End: %04i%c%02i%c%02i, %02i%c%02i%c%02i", year, '-', month, '-', day, hour, ':', minute, ':', (int)seconds);
				fprintf(fp, "\n\t\tAtomic Clock Referenced: %s", hdr[i_s].data_time_end_isAtomicClockReferenced ? "Yes" : "No");

				fprintf(fp, "\n\tScan Pose: ");
				fprintf(fp, "\n\t\tTranslation X: %lf", hdr[i_s].translation_x);
				fprintf(fp, "\n\t\tTranslation Y: %lf", hdr[i_s].translation_y);
				fprintf(fp, "\n\t\tTranslation Z: %lf", hdr[i_s].translation_z);

				fprintf(fp, "\n\t\tQuaternion x: %.12lf", hdr[i_s].quaternion_x);
				fprintf(fp, "\n\t\tQuaternion y: %.12lf", hdr[i_s].quaternion_y);
				fprintf(fp, "\n\t\tQuaternion z: %.12lf", hdr[i_s].quaternion_z);
				fprintf(fp, "\n\t\tQuaternion w: %.12lf", hdr[i_s].quaternion_w);

				fprintf(fp, "\n\tScale Factor: ");
				fprintf(fp, "\n\t\tDistance Minimum Limit: %lf", hdr[i_s].distance_minimum_limit);
				fprintf(fp, "\n\t\tDistance Maximum Limit: %lf", hdr[i_s].distance_maximum_limit);
				fprintf(fp, "\n\t\tDistance Scale: ");
				if (hdr[i_s].distance_scale_factor == -1)
				{
					fprintf(fp, "Integer");
				}
				else if (hdr[i_s].distance_scale_factor == 0)
				{
					fprintf(fp, "Float");
				}
				else
				{
					fprintf(fp, "\n%lf", hdr[i_s].distance_scale_factor);
				}

				fprintf(fp, "\n\t\tAngle Minimum Limit: %lf", hdr[i_s].angle_minimum_limit);
				fprintf(fp, "\n\t\tAngle Maximum Limit: %lf", hdr[i_s].angle_maximum_limit);
				fprintf(fp, "\n\t\tAngle Scale: ");
				if (hdr[i_s].angle_scale_factor == -1)
				{
					fprintf(fp, "Integer");
				}
				else if (hdr[i_s].angle_scale_factor == 0)
				{
					fprintf(fp, "Float");
				}
				else
				{
					fprintf(fp, "\n%lf", hdr[i_s].angle_scale_factor);
				}

				fprintf(fp, "\n\t\tIntsenity Minimum Limit: %lf", hdr[i_s].intensity_minimum_limit);
				fprintf(fp, "\n\t\tIntensity Maximum Limit: %lf", hdr[i_s].intensity_maximum_limit);
				fprintf(fp, "\n\t\tIntensity Scale: ");
				if (hdr[i_s].intensity_scale_factor == -1)
				{
					fprintf(fp, "Integer");
				}
				else if (hdr[i_s].intensity_scale_factor == 0)
				{
					fprintf(fp, "Float");
				}
				else
				{
					fprintf(fp, "\n%lf", hdr[i_s].intensity_scale_factor);
				}

				fprintf(fp, "\n\tXYZ Coorinate Info: ");
				fprintf(fp, "\n\t\tX: %s", hdr[i_s].x_field ? "Yes" : "No");
				fprintf(fp, "\n\t\tY: %s", hdr[i_s].y_field ? "Yes" : "No");
				fprintf(fp, "\n\t\tZ: %s", hdr[i_s].z_field ? "Yes" : "No");
				fprintf(fp, "\n\t\tInvalid Flag: %s", hdr[i_s].xyz_invalid_field ? "Yes" : "No");

				fprintf(fp, "\n\t\tX Minimum: %lf", hdr[i_s].x_minimum);
				fprintf(fp, "\n\t\tX Maximum: %lf", hdr[i_s].x_maximum);
				fprintf(fp, "\n\t\tY Minimum: %lf", hdr[i_s].y_minimum);
				fprintf(fp, "\n\t\tY Maximum: %lf", hdr[i_s].y_maximum);
				fprintf(fp, "\n\t\tZ Minimum: %lf", hdr[i_s].z_minimum);
				fprintf(fp, "\n\t\tZ Maximum: %lf", hdr[i_s].z_maximum);

				fprintf(fp, "\n\tSpherical Coordinate Info: ");
				fprintf(fp, "\n\t\tRange: %s", hdr[i_s].range_field ? "Yes" : "No");
				fprintf(fp, "\n\t\tAzimuth: %s", hdr[i_s].azimuth_field ? "Yes" : "No");
				fprintf(fp, "\n\t\tElevation: %s", hdr[i_s].elevation_field ? "Yes" : "No");
				fprintf(fp, "\n\t\tInvalid Flag: %s", hdr[i_s].spherical_invalid_field ? "Yes" : "No");

				fprintf(fp, "\n\t\tRange Minimum: %lf", hdr[i_s].range_minimum);
				fprintf(fp, "\n\t\tRange Maximum: %lf", hdr[i_s].range_maximum);
				fprintf(fp, "\n\t\tAzimuth Start: %lf", hdr[i_s].azimuth_start);
				fprintf(fp, "\n\t\tAzimuth End: %lf", hdr[i_s].azimuth_end);
				fprintf(fp, "\n\t\tElevation Minimum: %lf", hdr[i_s].elevation_minimum);
				fprintf(fp, "\n\t\tElevation Maximum: %lf", hdr[i_s].elevation_maximum);

				fprintf(fp, "\n\tRGB Info: ");
				fprintf(fp, "\n\t\tRed: %s", hdr[i_s].red_field ? "Yes" : "No");
				fprintf(fp, "\n\t\tGreen: %s", hdr[i_s].green_field ? "Yes" : "No");
				fprintf(fp, "\n\t\tBlue: %s", hdr[i_s].blue_field ? "Yes" : "No");
				fprintf(fp, "\n\t\tInvalid Flag: %s", hdr[i_s].rgb_invalid_field ? "Yes" : "No");

				fprintf(fp, "\n\t\tRed Minimum Limit: %lf", hdr[i_s].red_minimum_limit);
				fprintf(fp, "\n\t\tRed Maximum Limit: %lf", hdr[i_s].red_maximum_limit);
				fprintf(fp, "\n\t\tGreen Minimum Limit: %lf", hdr[i_s].green_minimum_limit);
				fprintf(fp, "\n\t\tGreen Maximum Limit: %lf", hdr[i_s].green_maximum_limit);
				fprintf(fp, "\n\t\tBlue Minimum Limit: %lf", hdr[i_s].blue_minimum_limit);
				fprintf(fp, "\n\t\tBlue Maximum Limit: %lf", hdr[i_s].blue_maximum_limit);

				fprintf(fp, "\n\tIntensity Info: ");
				fprintf(fp, "\n\t\tIntensity: %s", hdr[i_s].intensity_field ? "Yes" : "No");
				fprintf(fp, "\n\t\tInvalid Flag: %s", hdr[i_s].intensity_invalid_field ? "Yes" : "No");

				fprintf(fp, "\n\tIndex Info: ");
				fprintf(fp, "\n\t\tRow Index: %s", hdr[i_s].row_index_field ? "Yes" : "No");
				fprintf(fp, "\n\t\tColumn Index: %s", hdr[i_s].column_index_field ? "Yes" : "No");
				fprintf(fp, "\n\t\tReturn Index: %s", hdr[i_s].return_index_field ? "Yes" : "No");

				fprintf(fp, "\n\t\tRow Maximum Limit: %lli", hdr[i_s].row_maximum);
				fprintf(fp, "\n\t\tColumn Maximum Limit: %lli", hdr[i_s].column_maximum);
				fprintf(fp, "\n\t\tReturn Maximum Limit: %lli", hdr[i_s].return_maximum);
				fprintf(fp, "\n\t\tRow Minimum: %lli", hdr[i_s].row_minimum);
				fprintf(fp, "\n\t\tRow Maximum: %lli", hdr[i_s].row_maximum);
				fprintf(fp, "\n\t\tColumn Minimum: %lli", hdr[i_s].column_minimum);
				fprintf(fp, "\n\t\tColumn Maximum: %lli", hdr[i_s].column_maximum);
				fprintf(fp, "\n\t\tReturn Minimum: %lli", hdr[i_s].return_minimum);
				fprintf(fp, "\n\t\tReturn Maximum: %lli", hdr[i_s].return_maximum);

				fprintf(fp, "\n\tTime Info: ");
				fprintf(fp, "\n\t\tTime Stamp: %s", hdr[i_s].time_stamp_field ? "Yes" : "No");
				fprintf(fp, "\n\t\tInvalid Flag: %s", hdr[i_s].time_stamp_invalid_field ? "Yes" : "No");
				fprintf(fp, "\n\t\tTime Maximum Limit: %.12lf", hdr[i_s].time_stamp_maximum_limit);

				fprintf(fp, "\n\tPoint Count: %lli", hdr[i_s].point_count);
				fprintf(fp, "\n");
			}

			delete[] hdr;
			fclose(fp);
			delete[] fname_report;
		}

		time1 = omp_get_wtime();
		printf("\n\nAll Done! (Total Time: %.2lf seconds)", time1 - time0);
		getchar();
	}

	void Split_E57_Scan(char** input_list, int i_input_first, int i_input_last)
	{
		printf("\n===== Split E57 to Individual Scans! =====\n");

		int n_file = 0;
		char** file_list = file_list_from_file(input_list, i_input_first, i_input_last, n_file);
		c_sort_list(file_list, 0, n_file - 1);

		printf("\nInput Files: ");
		for (int i_file = 0; i_file < n_file; i_file++)
		{
			printf("\n\t%s", file_list[i_file]);
		}
		printf("\n");

		printf("\nProcess Started!\n");

		Activate_Multi_Threading();

		double time0, time1;
		time0 = omp_get_wtime();

		for (int i_file = 0; i_file < n_file; i_file++)
		{
			__int64 n_scn = e57_scan_count(file_list[i_file]);
			printf("\nReading File (%i/%i) %s", i_file + 1, n_file, file_list[i_file]);
			for (__int64 i_s = 0; i_s < n_scn; i_s++)
			{
				E57_io e57_io;
				e57_io.Add_E57_Scan(file_list[i_file], i_s);
				printf("\n\tScan (%lli/%lli) %s", i_s + 1, n_scn, e57_io.e57_scan[0].e57_header.name);

				char fname[E57_IO_MAX_STRING_SIZE];
				char* no_extension = remove_extension(file_list[i_file]);
				if (n_scn < 1000)
				{
					sprintf(fname, "%s_%03lli.e57", no_extension, i_s + 1);
				}
				else if (n_scn < 1000000)
				{
					sprintf(fname, "%s_%06lli.e57", no_extension, i_s + 1);
				}
				else
				{
					sprintf(fname, "%s_%lli.e57", no_extension, i_s + 1);
				}

				e57_io.e57_scan[0].Write_E57_Scan(fname);

				delete[] no_extension;
				e57_io.Delete();
			}
		}

		time1 = omp_get_wtime();
		printf("\n\nAll Done! (Total Time: %.2lf seconds)", time1 - time0);
		getchar();
	}
}