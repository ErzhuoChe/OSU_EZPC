#pragma once
#ifdef _MSC_VER 
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS 
#endif // !_CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <algorithm>
#include <omp.h>
#include <string>
#include <cstring>
#include <ctime>

#define EZPC_LIC_RND_WRD_LEN 66
#define EZPC_PATH_BUFFER	512

namespace EZPC
{
	const char Extension_LAS[] = ".las";
	const char Extension_LAZ[] = ".laz";
	const char Extension_E57[] = ".e57";
	const char Extension_LST[] = ".lst";

	void Print_EZPC_Banner();
	void Check_EZ_License(char* dir, char* lic_key);
	void Write_EZ_License(char* lic_key, int year, int month, int date);
	void Activate_Multi_Threading();
	void Activate_Multi_Threading(int n_thread);
	
	char* replace_extension(char* path, char* new_ext, bool full_or_truncate = true);
	char* add_extension(char* path, char* new_ext, bool behind_name_or_ext = true);
	char* remove_extension(char* path);
	char* get_extension(char* path);
	char* replace_directory(char* path, char* new_dir);
	char* remove_directory(char* path);
	char* get_directory(char* path);

	bool is_standard_extension(char* path);
	bool is_extension_matched(char* path_1, char* path_2);
	char* match_extension(char* path, char* ref_path);

	int n_file_from_lst(char* input_lst);
	char** file_list_from_lst(char* input_lst, int &n_file);
	char** file_list_from_file(char** input_file_list, int i_file_first, int i_file_last, int& n_file);
	void Write_File_List(char** file_list, int n_file, char* output_file, bool overwrite_or_append);

	long long f_count(bool* flag_list, long long n_flag, bool flag);
	long long f_count(bool* flag_list, long long i_flag_first, long long i_flag_last, bool flag);
	bool* f_create(long long n_flag, bool flag);
	bool* f_create_cpy(bool* flag_list, long long n_flag);
	long long* f_idxlist(bool* flag_list, long long n_flag, bool flag, long long& nflag);
	bool* f_intersect(bool* flag_1, bool* flag_2, long long n_flag);
	bool* f_union(bool* flag_1, bool* flag_2, long long n_flag);

	char** c_create_list(int n_file);
	char** c_cpy_list(char** file_list, int i_first, int i_last, int& n_file);
	char* c_to_upper(char* str);
	char* c_to_lower(char* str);
	int c_str_cmpi(char* str1, char* str2);
	void c_delete_list(char** file_list, int n_file);
	void c_print_list(const char* before, char** file_list, const char* after, int n_file);
	void c_print_list(const char* before, char** file_list, int n_file);
	void c_print_list(const char** file_list, const char* after, int n_file);
	void c_sort_list(char** file_list, int i_first, int i_last);
	bool c_compare_string(char* str_1, char* str_2);
}
