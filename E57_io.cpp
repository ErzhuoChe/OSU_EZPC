#include "E57_io.h"

using namespace EZPC;

// E57_Scan

void E57_Scan::Read_E57_Scan(char* filename, __int64 scn_idx)
{
	sprintf(source_file, "%s%c", filename, '\0');
	scan_index = scn_idx;
	Read_E57_Scan();
}

void E57_Scan::Read_E57_Scan()
{
	Read_E57_Header();
	Initialize();
	Read_E57_Point();
}

void E57_Scan::Read_E57_Header(char* filename, __int64 scn_idx)
{
	sprintf(source_file, "%s%c", filename, '\0');
	scan_index = scn_idx;
	Read_E57_Header();
}

void E57_Scan::Read_E57_Header()
{
	e57::Reader eReader(source_file);
	e57::Data3D header;
	eReader.ReadData3D((int32_t)scan_index, header);
	E57_Header_from_Data3D(header);
	eReader.Close();
}

void E57_Scan::Read_E57_Point()
{
	__int64 block_size = e57_header.point_count < E57_IO_BLOCK_SIZE ? e57_header.point_count : E57_IO_BLOCK_SIZE;
	__int64 block_count = (__int64)ceil((double)e57_header.point_count / (double)block_size);

	E57_Block block = create_e57_block(block_size);

	e57::Reader eReader(source_file);
	e57::CompressedVectorReader dataReader = eReader.SetUpData3DPointsData(
		(int32_t)scan_index,
		(int64_t)block_size,
		block.x,
		block.y,
		block.z,
		block.invalid_xyz,
		block.intensity,
		block.invalid_intensity,
		block.red,
		block.green,
		block.blue,
		block.invalid_rgb,
		block.range,
		block.azimuth,
		block.elevation,
		block.invalid_spherical,
		block.row_index,
		block.column_index,
		block.return_index,
		block.return_count,
		block.time_stamp,
		block.invalid_time_stamp	
	);

	for (__int64 i_block = 0; i_block < block_count; i_block++)
	{
		dataReader.read();
		__int64 i_p_start = i_block * block_size;
		__int64 npts = min(e57_header.point_count - i_p_start, block_size);

#pragma omp parallel for schedule(static)
		for (__int64 i_p = 0; i_p < npts; i_p++)
		{
			e57_point[i_p + i_p_start] = e57_point_from_block(block, i_p);
		}	
	}

	if ((!e57_header.x_field || !e57_header.y_field || !e57_header.z_field)
		&& (e57_header.range_field && e57_header.azimuth_field && e57_header.elevation_field))
	{
		Set_XYZ_From_Spherical();
	}

	dataReader.close();
	eReader.Close();
	Delete_E57_Block(block);
}

void E57_Scan::Copy_To_E57_Scan(E57_Scan& dst_scn)
{
	dst_scn.Delete();
	sprintf(dst_scn.source_file, "%s%c", source_file, '\0');
	dst_scn.scan_index = scan_index;
	dst_scn.e57_header = e57_header;
	dst_scn.Initialize();
#pragma omp parallel for schedule(static)
	for (__int64 i_p = 0; i_p < e57_header.point_count; i_p++)
	{
		dst_scn.e57_point[i_p] = e57_point[i_p];
	}
}

void E57_Scan::Copy_From_E57_Scan(E57_Scan src_scn)
{
	Delete();
	sprintf(source_file, "%s%c", src_scn.source_file, '\0');
	scan_index = src_scn.scan_index;
	e57_header = src_scn.e57_header;
	Initialize(); 
#pragma omp parallel for schedule(static)
	for (__int64 i_p = 0; i_p < e57_header.point_count; i_p++)
	{
		e57_point[i_p] = src_scn.e57_point[i_p];
	}
}

void E57_Scan::Set_XYZ_From_Spherical()
{
#pragma omp parallel for schedule (static)
	for (__int64 i_p = 0; i_p < e57_header.point_count; i_p++)
	{
		E57_Point p = e57_scan_point(i_p, true);
		e57_point[i_p].x = p.range * cos(p.elevation) * cos(p.azimuth);
		e57_point[i_p].y = p.range * cos(p.elevation) * sin(p.azimuth);
		e57_point[i_p].z = p.range * sin(p.elevation);
	}
	e57_header.x_field = true;
	e57_header.y_field = true;
	e57_header.z_field = true;
}

E57_Point E57_Scan::global_e57_point(E57_Point local_point)
{
	E57_Point global_point = local_point;

	double tx = e57_header.translation_x;
	double ty = e57_header.translation_y;
	double tz = e57_header.translation_z;

	double qx = e57_header.quaternion_x;
	double qy = e57_header.quaternion_y;
	double qz = e57_header.quaternion_z;
	double qw = e57_header.quaternion_w;

	double rot[9];
		
	rot[0] = 1. - 2. * qy * qy - 2. * qz * qz;
	rot[1] = 2. * qx * qy - 2. * qz * qw;
	rot[2] = 2. * qx * qz + 2. * qy * qw;

	rot[3] = 2. * qx * qy + 2. * qz * qw;
	rot[4] = 1. - 2. * qx * qx - 2. * qz * qz;
	rot[5] = 2. * qy * qz - 2. * qx * qw;

	rot[6] = 2. * qx * qz - 2. * qy * qw;
	rot[7] = 2. * qy * qz + 2. * qx * qw;
	rot[8] = 1. - 2. * qx * qx - 2. * qy * qy;

	global_point.x = local_point.x * rot[0] + local_point.y * rot[1] + local_point.z * rot[2] + tx;
	global_point.y = local_point.x * rot[3] + local_point.y * rot[4] + local_point.z * rot[5] + ty;
	global_point.z = local_point.x * rot[6] + local_point.y * rot[7] + local_point.z * rot[8] + tz;

	return global_point;
}

E57_Point E57_Scan::local_e57_point(E57_Point global_point)
{
	E57_Point local_point = global_point;

	double tx = e57_header.translation_x;
	double ty = e57_header.translation_y;
	double tz = e57_header.translation_z;

	double qx = e57_header.quaternion_x;
	double qy = e57_header.quaternion_y;
	double qz = e57_header.quaternion_z;
	double qw = e57_header.quaternion_w;

	double rot[9];

	rot[0] = 1. - 2. * qy * qy - 2. * qz * qz;
	rot[1] = 2. * qx * qy + 2. * qz * qw;
	rot[2] = 2. * qx * qz - 2. * qy * qw;

	rot[3] = 2. * qx * qy - 2. * qz * qw;
	rot[4] = 1. - 2. * qx * qx - 2. * qz * qz;
	rot[5] = 2. * qy * qz + 2. * qx * qw;

	rot[6] = 2. * qx * qz + 2. * qy * qw;
	rot[7] = 2. * qy * qz - 2. * qx * qw;
	rot[8] = 1. - 2. * qx * qx - 2. * qy * qy;

	local_point.x = (global_point.x - tx) * rot[0] + (global_point.y - ty) * rot[1] + (global_point.z - tz) * rot[2];
	local_point.y = (global_point.x - tx) * rot[3] + (global_point.y - ty) * rot[4] + (global_point.z - tz) * rot[5];
	local_point.z = (global_point.x - tx) * rot[6] + (global_point.y - ty) * rot[7] + (global_point.z - tz) * rot[8];

	return local_point;
}

E57_Point E57_Scan::e57_scan_point(__int64 point_index, bool local_or_global)
{
	if (local_or_global)
	{
		return e57_point[point_index];
	}
	else
	{
		return global_e57_point(e57_point[point_index]);
	}
}

void E57_Scan::Write_E57_Scan(char* filename)
{
	e57::Writer eWriter(filename, E57_IO_COORDINATE_METADATA); // it was supposed to be coordinate metadata
	Write_E57_Scan(eWriter);
	eWriter.Close();
}

void E57_Scan::Write_E57_Scan(e57::Writer eWriter) // this will keep the eWriter open.
{
	__int64 block_size = e57_header.point_count < E57_IO_BLOCK_SIZE ? e57_header.point_count : E57_IO_BLOCK_SIZE;
	__int64 block_count = (__int64)ceil((double)e57_header.point_count / (double)block_size);
	E57_Block block = create_e57_block(block_size);

	e57::Data3D header = Data3D_from_E57_Header();
	int32_t scan_index = eWriter.NewData3D(header);

	e57::CompressedVectorWriter dataWriter = eWriter.SetUpData3DPointsData(
		(int32_t)scan_index,
		(int64_t)block_size,
		block.x,
		block.y,
		block.z,
		block.invalid_xyz,
		block.intensity,
		block.invalid_intensity,
		block.red,
		block.green,
		block.blue,
		block.invalid_rgb,
		block.range,
		block.azimuth,
		block.elevation,
		block.invalid_spherical,
		block.row_index,
		block.column_index,
		block.return_index,
		block.return_count,
		block.time_stamp,
		block.invalid_time_stamp
	);

	for (__int64 i_block = 0; i_block < block_count; i_block++)
	{
		__int64 i_p_start = i_block * block_size;
		__int64 npts = min(e57_header.point_count - i_p_start, block_size);

#pragma omp parallel for schedule(static)
		for (__int64 i_p = 0; i_p < npts; i_p++)
		{
			Block_from_E57_Point(block, i_p, e57_point[i_p + i_p_start]);
		}

		dataWriter.write(block_size);
	}
	dataWriter.close();
	Delete_E57_Block(block);
}

e57::Data3D E57_Scan::Data3D_from_E57_Header()
{
	e57::Data3D header;
	header.pointsSize = (int64_t)e57_header.point_count;
	header.name = e57_header.name;
	header.guid = e57_header.guid;

	header.originalGuids.resize((size_t)e57_header.original_guid_count);
	for (int i = 0; i < header.originalGuids.size(); i++)
	{
		header.originalGuids[i] = e57_header.original_guid[i];
	}
	header.description = e57_header.description; 

	header.sensorVendor = e57_header.vendor; 
	header.sensorModel = e57_header.model; 

	header.sensorSerialNumber = e57_header.serial_number;
	header.sensorHardwareVersion = e57_header.hardware_version; 
	header.sensorSoftwareVersion = e57_header.software_version; 
	header.sensorFirmwareVersion = e57_header.firmware_version; 

	header.temperature = e57_header.temperature;
	header.relativeHumidity = e57_header.relative_humidity;
	header.atmosphericPressure = e57_header.atmospheric_pressure;

	header.acquisitionStart.dateTimeValue = e57_header.date_time_start;
	header.acquisitionEnd.dateTimeValue = e57_header.date_time_end;
	header.acquisitionStart.isAtomicClockReferenced = e57_header.data_time_start_isAtomicClockReferenced ? 1 : 0;
	header.acquisitionEnd.isAtomicClockReferenced = e57_header.data_time_end_isAtomicClockReferenced ? 1 : 0;

	header.pose.translation.x = e57_header.translation_x;
	header.pose.translation.y = e57_header.translation_y;
	header.pose.translation.z = e57_header.translation_z;

	header.pose.rotation.x = e57_header.quaternion_x;
	header.pose.rotation.y = e57_header.quaternion_y;
	header.pose.rotation.z = e57_header.quaternion_z;
	header.pose.rotation.w = e57_header.quaternion_w;

	header.cartesianBounds.xMinimum = e57_header.x_minimum;
	header.cartesianBounds.xMaximum = e57_header.x_maximum;
	header.cartesianBounds.yMinimum = e57_header.y_minimum;
	header.cartesianBounds.yMaximum = e57_header.y_maximum;
	header.cartesianBounds.zMinimum = e57_header.z_minimum;
	header.cartesianBounds.zMaximum = e57_header.z_maximum;

	header.sphericalBounds.rangeMinimum = e57_header.range_minimum;
	header.sphericalBounds.rangeMaximum = e57_header.range_maximum;
	header.sphericalBounds.azimuthStart = e57_header.azimuth_start;
	header.sphericalBounds.azimuthEnd = e57_header.azimuth_end;
	header.sphericalBounds.elevationMinimum = e57_header.elevation_minimum;
	header.sphericalBounds.elevationMaximum = e57_header.elevation_maximum;

	header.indexBounds.rowMinimum = (int64_t)e57_header.row_minimum;
	header.indexBounds.rowMaximum = (int64_t)e57_header.row_maximum;
	header.indexBounds.columnMinimum = (int64_t)e57_header.column_minimum;
	header.indexBounds.columnMaximum = (int64_t)e57_header.column_maximum;
	header.indexBounds.returnMinimum = (int64_t)e57_header.return_minimum;
	header.indexBounds.returnMaximum = (int64_t)e57_header.return_maximum;

	header.intensityLimits.intensityMinimum = e57_header.intensity_minimum_limit;
	header.intensityLimits.intensityMaximum = e57_header.intensity_maximum_limit;

	header.colorLimits.colorRedMinimum = e57_header.red_minimum_limit;
	header.colorLimits.colorRedMaximum = e57_header.red_maximum_limit;
	header.colorLimits.colorGreenMinimum = e57_header.green_minimum_limit;
	header.colorLimits.colorGreenMaximum = e57_header.green_maximum_limit;
	header.colorLimits.colorBlueMinimum = e57_header.blue_minimum_limit;
	header.colorLimits.colorBlueMaximum = e57_header.blue_maximum_limit;

	header.pointFields.cartesianXField = e57_header.x_field;
	header.pointFields.cartesianYField = e57_header.y_field;
	header.pointFields.cartesianZField = e57_header.z_field;
	header.pointFields.cartesianInvalidStateField = e57_header.xyz_invalid_field;

	header.pointFields.sphericalRangeField = e57_header.range_field;
	header.pointFields.sphericalAzimuthField = e57_header.azimuth_field;
	header.pointFields.sphericalElevationField = e57_header.elevation_field;
	header.pointFields.sphericalInvalidStateField = e57_header.spherical_invalid_field;

	header.pointFields.pointRangeMinimum = e57_header.distance_minimum_limit;
	header.pointFields.pointRangeMaximum = e57_header.distance_maximum_limit;
	header.pointFields.pointRangeScaledInteger = e57_header.distance_scale_factor; // -1 = integer; 0 = float; others = scale factor.	

	header.pointFields.angleMinimum = e57_header.angle_minimum_limit;
	header.pointFields.angleMaximum = e57_header.angle_maximum_limit;
	header.pointFields.angleScaledInteger = e57_header.angle_scale_factor; // -1 = integer; 0 = float; others = scale factor.	

	header.pointFields.rowIndexField = e57_header.row_index_field;
	header.pointFields.columnIndexField = e57_header.column_index_field;
	header.pointFields.rowIndexMaximum = (uint32_t)e57_header.row_maximum_limit;
	header.pointFields.columnIndexMaximum = (uint32_t)e57_header.column_maximum_limit;

	header.pointFields.returnIndexField = e57_header.return_index_field;
	header.pointFields.returnCountField = e57_header.return_count_field;
	header.pointFields.returnMaximum = (uint8_t)e57_header.return_maximum;

	header.pointFields.timeStampField = e57_header.time_stamp_field;
	header.pointFields.isTimeStampInvalidField = e57_header.time_stamp_invalid_field;
	header.pointFields.timeMaximum = e57_header.time_stamp_maximum_limit;

	header.pointFields.intensityField = e57_header.intensity_field;
	header.pointFields.isIntensityInvalidField = e57_header.intensity_invalid_field;
	header.pointFields.isIntensityInvalidField = e57_header.intensity_scale_factor;  // -1 = integer; 0 = float; others = scale factor.		

	header.pointFields.colorRedField = e57_header.red_field;
	header.pointFields.colorGreenField = e57_header.green_field;
	header.pointFields.colorBlueField = e57_header.blue_field;
	header.pointFields.isColorInvalidField = e57_header.rgb_invalid_field;

	return header;
}

void E57_Scan::E57_Header_from_Data3D(e57::Data3D header)
{
	e57_header.point_count = header.pointsSize;

	sprintf(e57_header.name, "%s%c", header.name.c_str(), '\0');
	sprintf(e57_header.guid, "%s%c", header.guid.c_str(), '\0');
	e57_header.original_guid_count = (int)header.originalGuids.size();
	for (int i = 0; i < header.originalGuids.size(); i++)
	{
		sprintf(e57_header.original_guid[i], "%s%c", header.originalGuids[i].c_str(), '\0');
	}
	sprintf(e57_header.description, "%s%c", header.description.c_str(), '\0');

	sprintf(e57_header.vendor, "%s%c", header.sensorVendor.c_str(), '\0');
	sprintf(e57_header.model, "%s%c", header.sensorModel.c_str(), '\0');
	sprintf(e57_header.serial_number, "%s%c", header.sensorSerialNumber.c_str(), '\0');
	sprintf(e57_header.hardware_version, "%s%c", header.sensorHardwareVersion.c_str(), '\0');
	sprintf(e57_header.software_version, "%s%c", header.sensorSoftwareVersion.c_str(), '\0');
	sprintf(e57_header.firmware_version, "%s%c", header.sensorFirmwareVersion.c_str(), '\0');

	e57_header.temperature = header.temperature;
	e57_header.relative_humidity = header.relativeHumidity;
	e57_header.atmospheric_pressure = header.atmosphericPressure;

	e57_header.date_time_start = header.acquisitionStart.dateTimeValue;
	e57_header.date_time_end = header.acquisitionEnd.dateTimeValue;
	e57_header.data_time_start_isAtomicClockReferenced = (header.acquisitionStart.isAtomicClockReferenced == 1);
	e57_header.data_time_end_isAtomicClockReferenced = (header.acquisitionEnd.isAtomicClockReferenced == 1);

	e57_header.translation_x = header.pose.translation.x;
	e57_header.translation_y = header.pose.translation.y;
	e57_header.translation_z = header.pose.translation.z;

	e57_header.quaternion_x = header.pose.rotation.x;
	e57_header.quaternion_y = header.pose.rotation.y;
	e57_header.quaternion_z = header.pose.rotation.z;
	e57_header.quaternion_w = header.pose.rotation.w;

	e57_header.x_minimum = header.cartesianBounds.xMinimum;
	e57_header.x_maximum = header.cartesianBounds.xMaximum;
	e57_header.y_minimum = header.cartesianBounds.yMinimum;
	e57_header.y_maximum = header.cartesianBounds.yMaximum;
	e57_header.z_minimum = header.cartesianBounds.zMinimum;
	e57_header.z_maximum = header.cartesianBounds.zMaximum;

	e57_header.range_minimum = header.sphericalBounds.rangeMinimum;
	e57_header.range_maximum = header.sphericalBounds.rangeMaximum;
	e57_header.azimuth_start = header.sphericalBounds.azimuthStart;
	e57_header.azimuth_end = header.sphericalBounds.azimuthEnd;
	e57_header.elevation_minimum = header.sphericalBounds.elevationMinimum;
	e57_header.elevation_maximum = header.sphericalBounds.elevationMaximum;

	e57_header.row_minimum = header.indexBounds.rowMinimum;
	e57_header.row_maximum = header.indexBounds.rowMaximum;
	e57_header.column_minimum = header.indexBounds.columnMinimum;
	e57_header.column_maximum = header.indexBounds.columnMaximum;
	e57_header.return_minimum = header.indexBounds.returnMinimum;
	e57_header.return_maximum = header.indexBounds.returnMaximum;

	e57_header.intensity_minimum_limit = header.intensityLimits.intensityMinimum;
	e57_header.intensity_maximum_limit = header.intensityLimits.intensityMaximum;

	e57_header.red_minimum_limit = header.colorLimits.colorRedMinimum;
	e57_header.red_maximum_limit = header.colorLimits.colorRedMaximum;
	e57_header.green_minimum_limit = header.colorLimits.colorGreenMinimum;
	e57_header.green_maximum_limit = header.colorLimits.colorGreenMaximum;
	e57_header.blue_minimum_limit = header.colorLimits.colorBlueMinimum;
	e57_header.blue_maximum_limit = header.colorLimits.colorBlueMaximum;

	e57_header.x_field = header.pointFields.cartesianXField;
	e57_header.y_field = header.pointFields.cartesianYField;
	e57_header.z_field = header.pointFields.cartesianZField;
	e57_header.xyz_invalid_field = header.pointFields.cartesianInvalidStateField;

	e57_header.range_field = header.pointFields.sphericalRangeField;
	e57_header.azimuth_field = header.pointFields.sphericalAzimuthField;
	e57_header.elevation_field = header.pointFields.sphericalElevationField;
	e57_header.spherical_invalid_field = header.pointFields.sphericalInvalidStateField;

	e57_header.distance_minimum_limit = header.pointFields.pointRangeMinimum;
	e57_header.distance_maximum_limit = header.pointFields.pointRangeMaximum;
	e57_header.distance_scale_factor = header.pointFields.pointRangeScaledInteger; // -1 = integer; 0 = float; others = scale factor.	

	e57_header.angle_minimum_limit = header.pointFields.angleMinimum;
	e57_header.angle_maximum_limit = header.pointFields.angleMaximum;
	e57_header.angle_scale_factor = header.pointFields.angleScaledInteger; // -1 = integer; 0 = float; others = scale factor.	

	e57_header.row_index_field = header.pointFields.rowIndexField;
	e57_header.column_index_field = header.pointFields.columnIndexField;
	e57_header.row_maximum_limit = header.pointFields.rowIndexMaximum;
	e57_header.column_maximum_limit = header.pointFields.columnIndexMaximum;

	e57_header.return_index_field = header.pointFields.returnIndexField;
	e57_header.return_count_field = header.pointFields.returnCountField;
	e57_header.return_maximum = header.pointFields.returnMaximum;

	e57_header.time_stamp_field = header.pointFields.timeStampField;
	e57_header.time_stamp_invalid_field = header.pointFields.isTimeStampInvalidField;
	e57_header.time_stamp_maximum_limit = header.pointFields.timeMaximum;

	e57_header.intensity_field = header.pointFields.intensityField;
	e57_header.intensity_invalid_field = header.pointFields.isIntensityInvalidField;
	e57_header.intensity_scale_factor = header.pointFields.isIntensityInvalidField;  // -1 = integer; 0 = float; others = scale factor.		

	e57_header.red_field = header.pointFields.colorRedField;
	e57_header.green_field = header.pointFields.colorGreenField;
	e57_header.blue_field = header.pointFields.colorBlueField;
	e57_header.rgb_invalid_field = header.pointFields.isColorInvalidField;
}

E57_Block E57_Scan::create_e57_block(__int64 block_size)
{
	E57_Block block;

	block.x = e57_header.x_field ? new double[block_size] : NULL;
	block.y = e57_header.y_field ? new double[block_size] : NULL;
	block.z = e57_header.z_field ? new double[block_size] : NULL;
	block.invalid_xyz = e57_header.xyz_invalid_field ? new int8_t[block_size] : NULL;

	block.intensity = e57_header.intensity_field ? new double[block_size] : NULL;
	block.invalid_intensity = e57_header.intensity_invalid_field ? new int8_t[block_size] : NULL;

	block.red = e57_header.red_field ? new uint16_t[block_size] : NULL;
	block.green = e57_header.green_field ? new uint16_t[block_size] : NULL;
	block.blue = e57_header.blue_field ? new uint16_t[block_size] : NULL;
	block.invalid_rgb = e57_header.rgb_invalid_field ? new int8_t[block_size] : NULL;

	block.range = e57_header.range_field ? new double[block_size] : NULL;
	block.azimuth = e57_header.azimuth_field ? new double[block_size] : NULL;
	block.elevation = e57_header.elevation_field ? new double[block_size] : NULL;
	block.invalid_spherical = e57_header.spherical_invalid_field ? new int8_t[block_size] : NULL;

	block.time_stamp = e57_header.time_stamp_field ? new double[block_size] : NULL;
	block.invalid_time_stamp = e57_header.time_stamp_invalid_field ? new int8_t[block_size] : NULL;
	
	block.row_index = e57_header.row_index_field ? new int32_t[block_size] : NULL;
	block.column_index = e57_header.column_index_field ? new int32_t[block_size] : NULL;

	block.return_index = e57_header.return_index_field ? new int8_t[block_size] : NULL;
	block.return_count = e57_header.return_count_field ? new int8_t[block_size] : NULL;

	return block;
}

void E57_Scan::Delete_E57_Block(E57_Block &e57_block)
{
	if (e57_block.x)
	{
		delete[] e57_block.x;
		e57_block.x = NULL;
	}
	if (e57_block.y)
	{
		delete[] e57_block.y;
		e57_block.y = NULL;
	}
	if (e57_block.z)
	{
		delete[] e57_block.z;
		e57_block.z = NULL;
	}
	if (e57_block.invalid_xyz)
	{
		delete[] e57_block.invalid_xyz;
		e57_block.invalid_xyz = NULL;
	}
	if (e57_block.intensity)
	{
		delete[] e57_block.intensity;
		e57_block.intensity = NULL;
	}
	if (e57_block.invalid_intensity)
	{
		delete[] e57_block.invalid_intensity;
		e57_block.invalid_intensity = NULL;
	}
	if (e57_block.red)
	{
		delete[] e57_block.red;
		e57_block.red = NULL;
	}
	if (e57_block.green)
	{
		delete[] e57_block.green;
		e57_block.green = NULL;
	}
	if (e57_block.blue)
	{
		delete[] e57_block.blue;
		e57_block.blue = NULL;
	}
	if (e57_block.invalid_rgb)
	{
		delete[] e57_block.invalid_rgb;
		e57_block.invalid_rgb = NULL;
	}

	if (e57_block.range)
	{
		delete[] e57_block.range;
		e57_block.range = NULL;
	}
	if (e57_block.azimuth)
	{
		delete[] e57_block.azimuth;
		e57_block.azimuth = NULL;
	}
	if (e57_block.elevation)	
	{
		delete[] e57_block.elevation;
		e57_block.elevation = NULL;
	}
	if (e57_block.invalid_spherical)
	{
		delete[] e57_block.invalid_spherical;
		e57_block.invalid_spherical = NULL;
	}

	if (e57_block.time_stamp)
	{
		delete[] e57_block.time_stamp;
		e57_block.time_stamp = NULL;
	}
	if (e57_block.invalid_time_stamp)
	{
		delete[] e57_block.invalid_time_stamp;
		e57_block.invalid_time_stamp = NULL;
	}

	if (e57_block.row_index)
	{
		delete[] e57_block.row_index;
		e57_block.row_index = NULL;
	}
	if (e57_block.column_index)
	{
		delete[] e57_block.column_index;
		e57_block.column_index = NULL;
	}
	if (e57_block.return_index)
	{
		delete[] e57_block.return_index;
		e57_block.return_index = NULL;
	}
	if (e57_block.return_count)
	{
		delete[] e57_block.return_count;
		e57_block.return_count = NULL;
	}
}

E57_Point E57_Scan::e57_point_from_block(E57_Block block, __int64 block_idx)
{
	E57_Point point; 
	if (block.x)
		point.x = block.x[block_idx];
	if (block.y)
		point.y = block.y[block_idx];
	if(block.z)
		point.z = block.z[block_idx]; 
	if(block.invalid_xyz)
		point.invalid_xyz = block.invalid_xyz[block_idx] == 0 ? false : true;
	if(block.intensity)
		point.intensity = block.intensity[block_idx];
	if(block.invalid_intensity)
		point.invalid_intensity = block.invalid_intensity[block_idx] == 0 ? false : true;
	if (block.red)
		point.red = block.red[block_idx];
	if (block.green)
		point.green = block.green[block_idx];
	if (block.blue)
		point.blue = block.blue[block_idx];
	if (block.invalid_rgb)
		point.invalid_rgb = block.invalid_rgb[block_idx] == 0 ? false : true;
	if (block.range)
		point.range = block.range[block_idx];
	if (block.azimuth)
		point.azimuth = block.azimuth[block_idx];
	if (block.elevation)
		point.elevation = block.elevation[block_idx];
	if (block.invalid_spherical)
		point.invalid_spherical = block.invalid_spherical[block_idx] == 0 ? false : true;
	if (block.row_index)
		point.row_index = block.row_index[block_idx];
	if (block.column_index)
		point.column_index = block.column_index[block_idx];
	if (block.return_index)
		point.return_index = block.return_index[block_idx];
	if (block.return_count)
		point.return_count = block.return_count[block_idx];
	if (block.time_stamp)
		point.time_stamp = block.time_stamp[block_idx];
	if (block.invalid_time_stamp)
		point.invalid_time_stamp = block.invalid_time_stamp[block_idx] == 0 ? false : true;
	return point;
}

void E57_Scan::Block_from_E57_Point(E57_Block &block, __int64 block_idx, E57_Point point)
{
	if (block.x)
		block.x[block_idx] = point.x;
	if (block.y)
		block.y[block_idx] = point.y;
	if (block.z)
		block.z[block_idx] = point.z;
	if (block.invalid_xyz)
		block.invalid_xyz[block_idx] = point.invalid_xyz ? 1 : 0;
	if (block.intensity)
		block.intensity[block_idx] = point.intensity;
	if (block.invalid_intensity)
		block.invalid_intensity[block_idx] = point.invalid_intensity ? 1 : 0;
	if (block.red)
		block.red[block_idx] = point.red;
	if (block.green)
		block.green[block_idx] = point.green;
	if (block.blue)
		block.blue[block_idx] = point.blue;
	if (block.invalid_rgb)
		block.invalid_rgb[block_idx] = point.invalid_rgb ? 1 : 0;
	if (block.range)
		block.range[block_idx] = point.range;
	if (block.azimuth)
		block.azimuth[block_idx] = point.azimuth;
	if (block.elevation)
		block.elevation[block_idx] = point.elevation;
	if (block.invalid_spherical)
		block.invalid_spherical[block_idx] = point.invalid_spherical ? 1 : 0;
	if (block.row_index)
		block.row_index[block_idx] = (int32_t)point.row_index;
	if (block.column_index)
		block.column_index[block_idx] = (int32_t)point.column_index;
	if (block.return_index)
		block.return_index[block_idx] = (int8_t)point.return_index;
	if (block.return_count)
		block.return_count[block_idx] = (int8_t)point.return_count;
	if (block.time_stamp)
		block.time_stamp[block_idx] = point.time_stamp;
	if (block.invalid_time_stamp)
		block.invalid_time_stamp[block_idx] = point.invalid_time_stamp;
}

void E57_Scan::Initialize()
{
	if (e57_point)
	{
		Delete();
	}
	e57_point = new E57_Point[e57_header.point_count];
}
	
void E57_Scan::Initialize(E57_Header header)
{
	e57_header = header;
	Initialize();
}

void E57_Scan::Initialize(__int64 n_pts)
{
	e57_header.point_count = n_pts;
	Initialize();
}

void E57_Scan::Delete()
{
	if (e57_point)
	{
		delete[] e57_point;
		e57_point = NULL;
	}
}

void E57_Scan::Show_E57_Scan_Info()
{
	printf("\nE57 Scan Info:");
	printf("\n\tSource File: %s", source_file);
	printf("\n\tScan Index: %lli", scan_index);
	Print_E57_Header(e57_header);
}

void E57_Scan::Show_E57_Point(__int64 i_p)
{

}

// E57_io

void E57_io::Add_E57_Scan(char** file_list, int i_file_start, int i_file_end)
{
	int n_file = i_file_end - i_file_start + 1;
	for (int i_file = i_file_start; i_file <= i_file_end; i_file++)
	{
		__int64 n_scn = e57_scan_count(file_list[i_file]);
		E57_Scan* scn_tmp = new E57_Scan[n_scn];

		for (__int64 i_s = 0; i_s < n_scn; i_s++)
		{
			scn_tmp[i_s].Read_E57_Scan(file_list[i_file], i_s);
		}

		Add_E57_Scan(scn_tmp, n_scn);

#pragma omp parallel for schedule (guided)
		for (__int64 i_s = 0; i_s < n_scn; i_s++)
		{
			scn_tmp[i_s].Delete();
		}
		delete[] scn_tmp;
		scn_tmp = NULL;
	}
}

void E57_io::Add_E57_Scan(char* file_name)
{
	Add_E57_Scan(&file_name, 0, 0);
}

void E57_io::Add_E57_Scan(char* file_name, __int64 scan_index)
{
	e57::ustring fname(file_name);
	e57::Reader eReader(fname);

	__int64 n_scn = (__int64)eReader.GetData3DCount();

	eReader.Close();

	E57_Scan scn_tmp;
	scn_tmp.Read_E57_Scan(file_name, scan_index);

	Add_E57_Scan(scn_tmp);
	scn_tmp.Delete();

}

void E57_io::Add_E57_Scan(E57_Scan add_scans)
{
	Add_E57_Scan(&add_scans, 1);
}

void E57_io::Add_E57_Scan(E57_Scan* add_scans, __int64 n_add_scn)
{
	E57_Scan* scn_tmp = e57_scan;
	e57_scan = new E57_Scan[scan_count + n_add_scn];

#pragma omp parallel for schedule(guided)
	for (__int64 i_s = 0; i_s < scan_count; i_s++)
	{
		e57_scan[i_s].Copy_From_E57_Scan(scn_tmp[i_s]);
		scn_tmp[i_s].Delete();
	}

#pragma omp parallel for schedule(guided)
	for (__int64 i_s = scan_count; i_s < scan_count + n_add_scn; i_s++)
	{
		e57_scan[i_s].Copy_From_E57_Scan(add_scans[i_s - scan_count]);
	}

	scan_count += n_add_scn;
	if (scn_tmp)
	{
		delete[] scn_tmp;
		scn_tmp = NULL;
	}
}

void E57_io::Initialize()
{
	if (e57_scan)
	{
		Delete();
	}
	e57_scan = new E57_Scan[scan_count];
}

void E57_io::Initialize(int n_scn)
{
	scan_count = n_scn;
	Initialize();
}

E57_Header E57_io::e57_scan_header(__int64 scan_index)
{
	return e57_scan[scan_index].e57_header;
}

E57_Point E57_io::e57_scan_point(__int64 scan_index, __int64 point_index, bool local_or_global)
{
	return e57_scan[scan_index].e57_scan_point(point_index, local_or_global);
}

void E57_io::Write_E57_Scan(char* file_name, bool merge_or_seperate)
{
	if (merge_or_seperate)
	{
		e57::Writer eWriter(file_name, E57_IO_COORDINATE_METADATA);
		for (__int64 i_s = 0; i_s < scan_count; i_s++)
		{
			e57_scan[i_s].Write_E57_Scan(eWriter);
		}
		eWriter.Close();
	}
	else
	{
		for (__int64 i_s = 0; i_s < scan_count; i_s++)
		{
			char scan_name[E57_IO_MAX_STRING_SIZE];
			char* no_extension = remove_extension(file_name);
			if (scan_count < 1000)
			{
				sprintf(scan_name, "%s_%03lli.e57", no_extension, i_s + 1);
			}
			else if (scan_count < 1000000)
			{
				sprintf(scan_name, "%s_%06lli.e57", no_extension, i_s + 1);
			}
			else
			{
				sprintf(scan_name, "%s_%lli.e57", no_extension, i_s + 1);
			}

			e57::Writer eWriter(file_name, E57_IO_COORDINATE_METADATA);
			e57_scan[i_s].Write_E57_Scan(eWriter);
			eWriter.Close();

			delete[] no_extension;
		}
	}
	

}

void E57_io::Show_E57_Scan_Info(__int64 i_scan)
{
	printf("\nScan #%lli:", i_scan);
	e57_scan[i_scan].Show_E57_Scan_Info();
}

void E57_io::Delete()
{
	if (e57_scan)
	{
		for (__int64 i_s = 0; i_s < scan_count; i_s++)
		{
			e57_scan[i_s].Delete();
		}
		delete[] e57_scan;
		e57_scan = NULL;
	}
}


// Common

__int64 e57_scan_count(char* file_name)
{
	e57::Reader eReader(file_name);
	__int64 n_scn = (__int64)eReader.GetData3DCount();
	eReader.Close();
	return n_scn;
}

E57_Header* Get_E57_Header(char** file_list, int i_file_start, int i_file_end, __int64 &header_count)
{
	header_count = 0;

	for (int i_file = i_file_start; i_file <= i_file_end; i_file++)
	{
		e57::Reader eReader(file_list[i_file]);
		header_count += (__int64)eReader.GetData3DCount();
		eReader.Close();
	}

	if (header_count > 0)
	{
		__int64 i_header = 0;
		E57_Header* header = new E57_Header[header_count];
		for (int i_file = i_file_start; i_file <= i_file_end; i_file++)
		{
			__int64 n_scn = e57_scan_count(file_list[i_file]);
			for (__int64 i_s = 0; i_s < n_scn; i_s++)
			{
				E57_Scan scn_tmp;
				scn_tmp.Read_E57_Header(file_list[i_file], i_s);
				header[i_header] = scn_tmp.e57_header;
				i_header++;
			}
		}

		return header;
	}
	else
	{
		return NULL;
	}
}

E57_Header* Get_E57_Header(char* file_name, __int64 &header_count)
{
	return Get_E57_Header(&file_name, 0, 0, header_count);
}

void Print_E57_Header(E57_Header header)
{
	printf("\nE57 Header: ");
	printf("\n\tScan Info: ");
	printf("\n\t\tName: %s", header.name);
	printf("\n\t\tGUID: %s", header.guid);

	printf("\n\t\tOriginal GUIDs: ");
	for (int i = 0; i < header.original_guid_count; i++)
	{
		printf("\n\t\t\t#%i: %s", i, header.original_guid[i]);
	}
	printf("\n\t\tDescription: %s", header.description);

	printf("\n\tSensor Info: ");
	printf("\n\t\tVendor: %s", header.vendor);
	printf("\n\t\tModel: %s", header.model);
	printf("\n\t\tSerial Number: %s", header.serial_number);
	printf("\n\t\tHardware Version: %s", header.hardware_version);
	printf("\n\t\tSoftware Version: %s", header.software_version);
	printf("\n\t\tFirmware Version: %s", header.firmware_version);

	printf("\n\tAcquisition Info: ");
	printf("\n\t\tTemperautre: %f", header.temperature);
	printf("\n\t\tRelative Humidity: %f", header.relative_humidity);
	printf("\n\t\tAtmospheric Pressure: %f", header.atmospheric_pressure);

	e57::DateTime datetime;
	int year, month, day, hour, minute;
	float seconds;
	
	datetime.dateTimeValue = header.date_time_start;
	datetime.GetUTCDateTime(year, month, day, hour, minute, seconds);
	printf("\n\t\tDate Time Start: %04i%c%02i%c%02i, %02i%c%02i%c%02i", year, '-', month, '-', day, hour, ':', minute, ':', (int)seconds);
	printf("\n\t\tAtomic Clock Referenced: %s", header.data_time_start_isAtomicClockReferenced ? "Yes" : "No");

	datetime.dateTimeValue = header.date_time_end;
	datetime.GetUTCDateTime(year, month, day, hour, minute, seconds);
	printf("\n\t\tDate Time End: %04i%c%02i%c%02i, %02i%c%02i%c%02i", year, '-', month, '-', day, hour, ':', minute, ':', (int)seconds);
	printf("\n\t\tAtomic Clock Referenced: %s", header.data_time_end_isAtomicClockReferenced ? "Yes" : "No");	

	printf("\n\tScan Pose: ");
	printf("\n\t\tTranslation X: %lf", header.translation_x);
	printf("\n\t\tTranslation Y: %lf", header.translation_y);
	printf("\n\t\tTranslation Z: %lf", header.translation_z);

	printf("\n\t\tQuaternion x: %.12lf", header.quaternion_x);
	printf("\n\t\tQuaternion y: %.12lf", header.quaternion_y);
	printf("\n\t\tQuaternion z: %.12lf", header.quaternion_z);
	printf("\n\t\tQuaternion w: %.12lf", header.quaternion_w);

	printf("\n\tScale Factor: ");
	printf("\n\t\tDistance Minimum Limit: %lf", header.distance_minimum_limit);
	printf("\n\t\tDistance Maximum Limit: %lf", header.distance_maximum_limit);
	printf("\n\t\tDistance Scale: ");
	if (header.distance_scale_factor == -1)
	{
		printf("Integer");
	}
	else if (header.distance_scale_factor == 0)
	{
		printf("Float");
	}
	else
	{
		printf("\n%lf", header.distance_scale_factor);
	}

	printf("\n\t\tAngle Minimum Limit: %lf", header.angle_minimum_limit);
	printf("\n\t\tAngle Maximum Limit: %lf", header.angle_maximum_limit);
	printf("\n\t\tAngle Scale: ");
	if (header.angle_scale_factor == -1)
	{
		printf("Integer");
	}
	else if (header.angle_scale_factor == 0)
	{
		printf("Float");
	}
	else
	{
		printf("\n%lf", header.angle_scale_factor);
	}

	printf("\n\t\tIntsenity Minimum Limit: %lf", header.intensity_minimum_limit);
	printf("\n\t\tIntensity Maximum Limit: %lf", header.intensity_maximum_limit);
	printf("\n\t\tIntensity Scale: ");
	if (header.intensity_scale_factor == -1)
	{
		printf("Integer");
	}
	else if (header.intensity_scale_factor == 0)
	{
		printf("Float");
	}
	else
	{
		printf("\n%lf", header.intensity_scale_factor);
	}

	printf("\n\tXYZ Coorinate Info: ");
	printf("\n\t\tX: %s", header.x_field ? "Yes" : "No");
	printf("\n\t\tY: %s", header.y_field ? "Yes" : "No");
	printf("\n\t\tZ: %s", header.z_field ? "Yes" : "No");
	printf("\n\t\tInvalid Flag: %s", header.xyz_invalid_field ? "Yes" : "No");
	
	printf("\n\t\tX Minimum: %lf", header.x_minimum);
	printf("\n\t\tX Maximum: %lf", header.x_maximum);
	printf("\n\t\tY Minimum: %lf", header.y_minimum);
	printf("\n\t\tY Maximum: %lf", header.y_maximum);
	printf("\n\t\tZ Minimum: %lf", header.z_minimum);
	printf("\n\t\tZ Maximum: %lf", header.z_maximum);

	printf("\n\tSpherical Coordinate Info: ");
	printf("\n\t\tRange: %s", header.range_field ? "Yes" : "No");
	printf("\n\t\tAzimuth: %s", header.azimuth_field ? "Yes" : "No");
	printf("\n\t\tElevation: %s", header.elevation_field ? "Yes" : "No");
	printf("\n\t\tInvalid Flag: %s", header.spherical_invalid_field ? "Yes" : "No");

	printf("\n\t\tRange Minimum: %lf", header.range_minimum);
	printf("\n\t\tRange Maximum: %lf", header.range_maximum);
	printf("\n\t\tAzimuth Start: %lf", header.azimuth_start);
	printf("\n\t\tAzimuth End: %lf", header.azimuth_end);
	printf("\n\t\tElevation Minimum: %lf", header.elevation_minimum);
	printf("\n\t\tElevation Maximum: %lf", header.elevation_maximum);

	printf("\n\tRGB Info: ");
	printf("\n\t\tRed: %s", header.red_field ? "Yes" : "No");
	printf("\n\t\tGreen: %s", header.green_field ? "Yes" : "No");
	printf("\n\t\tBlue: %s", header.blue_field ? "Yes" : "No");
	printf("\n\t\tInvalid Flag: %s", header.rgb_invalid_field ? "Yes" : "No");

	printf("\n\t\tRed Minimum Limit: %lf", header.red_minimum_limit);
	printf("\n\t\tRed Maximum Limit: %lf", header.red_maximum_limit);
	printf("\n\t\tGreen Minimum Limit: %lf", header.green_minimum_limit);
	printf("\n\t\tGreen Maximum Limit: %lf", header.green_maximum_limit);
	printf("\n\t\tBlue Minimum Limit: %lf", header.blue_minimum_limit);
	printf("\n\t\tBlue Maximum Limit: %lf", header.blue_maximum_limit);

	printf("\n\tIntensity Info: ");
	printf("\n\t\tIntensity: %s", header.intensity_field ? "Yes" : "No");
	printf("\n\t\tInvalid Flag: %s", header.intensity_invalid_field ? "Yes" : "No");

	printf("\n\tIndex Info: ");
	printf("\n\t\tRow Index: %s", header.row_index_field ? "Yes" : "No");
	printf("\n\t\tColumn Index: %s", header.column_index_field ? "Yes" : "No");
	printf("\n\t\tReturn Index: %s", header.return_index_field ? "Yes" : "No");

	printf("\n\t\tRow Maximum Limit: %lli", header.row_maximum);
	printf("\n\t\tColumn Maximum Limit: %lli", header.column_maximum);
	printf("\n\t\tReturn Maximum Limit: %lli", header.return_maximum);
	printf("\n\t\tRow Minimum: %lli", header.row_minimum);
	printf("\n\t\tRow Maximum: %lli", header.row_maximum);
	printf("\n\t\tColumn Minimum: %lli", header.column_minimum);
	printf("\n\t\tColumn Maximum: %lli", header.column_maximum);
	printf("\n\t\tReturn Minimum: %lli", header.return_minimum);
	printf("\n\t\tReturn Maximum: %lli", header.return_maximum);

	printf("\n\tTime Info: ");
	printf("\n\t\tTime Stamp: %s", header.time_stamp_field ? "Yes" : "No");
	printf("\n\t\tInvalid Flag: %s", header.time_stamp_invalid_field ? "Yes" : "No");
	printf("\n\t\tTime Maximum Limit: %.12lf", header.time_stamp_maximum_limit);

	printf("\n\tPoint Count: %lli", header.point_count);
}
