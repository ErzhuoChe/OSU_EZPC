#include "EZPC_util.h"

namespace EZPC
{
	void Print_EZPC_Banner()
	{
		printf("\n-----------------------------------------------");
		printf("\n!!!!?!!?!!!!!!!!?!!!!!!!!!?!!!!!!!!!?!!!!!!!!!!");
		printf("\n!!?!!!!!         Buggy EZPC            !!!!!!!!");
		printf("\n!!!?!!!!         Buggilized by         !!!!?!!!");
		printf("\n!?!!!!!!         Buggy Ezra Che        !!!!!!!!");
		printf("\n!!!!!!!!!?!!!!!!!?!!!!!!!!!!!!!?!!!!!!!!?!!!!!!");
		printf("\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
		printf("\n");
	}

	void Check_EZ_License(char* dir, char* lic_key)
	{
		char* lic_path;
		char* lic_dir = get_directory(dir);
		lic_path = new char[strlen(lic_dir) + strlen("simple_license") + 1];
		memcpy(lic_path, lic_dir, strlen(lic_dir));
		memcpy(lic_path + strlen(lic_dir), "EZ_License", strlen("EZ_License"));
		lic_path[strlen(lic_dir) + strlen("EZ_License")] = '\0';

		FILE* fp = fopen(lic_path, "rb");
		if (!fp)
		{
			printf("\n\nCannot Find the License File!!");
			getchar();
			exit(1);
		}
		else
		{
			int len_key;
			char key[EZPC_LIC_RND_WRD_LEN];
			time_t rawtime;
			struct tm* timeinfo;

			time(&rawtime);
			timeinfo = gmtime(&rawtime);

			fseek(fp, EZPC_LIC_RND_WRD_LEN, SEEK_SET);

			fread(&len_key, sizeof(len_key), 1, fp);
			fread(key, len_key, 1, fp);

			if (strcmp(lic_key, key) != 0)
			{
				printf("\n\nWrong License Key!!");
				getchar();
				exit(1);
			}

			int yr;
			int mon;
			int day;

			fread(&yr, sizeof(yr), 1, fp);
			fread(&mon, sizeof(mon), 1, fp);
			fread(&day, sizeof(day), 1, fp);

			int val_lic = yr * 10000 + mon * 100 + day;
			int val_now = (timeinfo->tm_year + 1900) * 10000 + (timeinfo->tm_mon + 1) * 100 + timeinfo->tm_mday;

			if (val_lic - val_now < 0)
			{
				printf("\n\nThe License is Expired!!");
				getchar();
				exit(1);
			}
			else
			{
				printf("\n * The License is Valid until %02i-%02i-%04i * \n", mon, day, yr);
			}

			fclose(fp);
		}

		if (lic_dir)
		{
			delete[] lic_dir;
			lic_dir = NULL;
		}

		if (lic_path)
		{
			delete[] lic_path;
			lic_path = NULL;
		}
	}

	void Write_EZ_License(char* lic_key, int year, int month, int date)
	{
		FILE* fp = fopen("EZ_License", "wb");
		char random_words[EZPC_LIC_RND_WRD_LEN] = "EZPC_EZ_License!";

		char key[EZPC_LIC_RND_WRD_LEN];
		int len_key = (int)strlen(lic_key) + 1;
		memcpy(key, lic_key, strlen(lic_key));
		key[strlen(lic_key)] = '\0';

		fwrite(random_words, EZPC_LIC_RND_WRD_LEN, 1, fp);
		fwrite(&len_key, sizeof(len_key), 1, fp);
		fwrite(key, len_key, 1, fp);

		fwrite(&year, sizeof(year), 1, fp);
		fwrite(&month, sizeof(month), 1, fp);
		fwrite(&date, sizeof(date), 1, fp);

		fwrite(random_words, EZPC_LIC_RND_WRD_LEN, 1, fp);

		fclose(fp);
	}

	void Activate_Multi_Threading(int n_thread)
	{
		omp_set_num_threads(n_thread);
	}

	void Activate_Multi_Threading()
	{
		Activate_Multi_Threading(omp_get_max_threads());
	}

	char* replace_extension(char* path, char* new_ext, bool full_or_truncate)
	{
		char* no_extension = remove_extension(path);
		char* add_extension = new char[strlen(no_extension) + strlen(new_ext) + 1];
		memcpy(add_extension, no_extension, strlen(no_extension));
		memcpy(add_extension + strlen(no_extension), new_ext, strlen(new_ext));
		add_extension[strlen(no_extension) + strlen(new_ext)] = '\0';

		if (no_extension)
		{
			delete[] no_extension;
			no_extension = NULL;
		}
		
		if (full_or_truncate)
		{
			return add_extension;
		}
		else
		{
			char* no_directory = remove_directory(add_extension);
			if (add_extension)
			{
				delete[] add_extension;
				add_extension = NULL;
			}
			return no_directory;
		}
	}

	char* add_extension(char* path, char* new_ext, bool behind_name_or_ext)
	{
		char* nu_path = new char[strlen(path) + strlen(new_ext) + 1];

		if (behind_name_or_ext)
		{
			char* no_extension = remove_extension(path);
			char* ext = get_extension(path);
			memcpy(nu_path, no_extension, strlen(no_extension));
			memcpy(nu_path + strlen(no_extension), new_ext, strlen(new_ext));
			memcpy(nu_path + strlen(no_extension) + strlen(new_ext), ext, strlen(ext));

			if (ext)
			{
				delete[] ext;
				ext = NULL;
			}
			if (no_extension)
			{
				delete[] no_extension;
				no_extension = NULL;
			}
		}
		else
		{
			memcpy(nu_path, path, strlen(path));
			memcpy(nu_path + strlen(path), new_ext, strlen(new_ext));
		}

		nu_path[strlen(path) + strlen(new_ext)] = '\0';
		return nu_path;
	}

	char* remove_extension(char* path)
	{
		int len = (int)strlen(path);
		int pos = len;

		for (int i = 0; i < len; i++)
		{
			if (path[len - i - 1] == '.')
			{
				pos = len - i - 1;
				break;
			}
			if (path[len - i - 1] == '\\' || path[len - i - 1] == '/')
			{
				break;
			}
		}

		char* fname = new char[pos + 1];
		memcpy(fname, path, pos);
		fname[pos] = '\0';
		return fname;
	}

	char* get_extension(char* path)
	{
		int len = (int)strlen(path);
		int pos = len;

		for (int i = 0; i < len; i++)
		{
			if (path[len - i - 1] == '.')
			{
				pos = len - i - 1;
				break;
			}
			if (path[len - i - 1] == '\\' || path[len - i - 1] == '/')
			{
				break;
			}
		}

		char* ext = new char[len - pos + 1];
		memcpy(ext, path + pos, len - pos);
		ext[len - pos] = '\0';
		return ext;
	}

	char* replace_directory(char* path, char* new_dir)
	{
		if (new_dir)
		{
			size_t len_new_dir = strlen(new_dir);

			size_t dir_buff_size = len_new_dir + (new_dir[len_new_dir - 1] != '\\' && new_dir[len_new_dir - 1] != '/' ? 2 : 1);
			char* new_dir_buff = new char[dir_buff_size];

			memcpy(new_dir_buff, new_dir, len_new_dir);
			if (new_dir[len_new_dir - 1] != '\\' && new_dir[len_new_dir - 1] != '/')
			{
				new_dir[len_new_dir] = '\\';
			}
			new_dir_buff[dir_buff_size - 1] = '\0';

			char* no_directory = remove_directory(path);
			size_t new_path_size = strlen(no_directory) + strlen(new_dir_buff) + 1;
			char* new_path = new char[new_path_size];

			memcpy(new_path, new_dir_buff, strlen(new_dir_buff));
			memcpy(new_path + strlen(new_dir_buff), no_directory, strlen(no_directory));

			new_path[new_path_size - 1] = '\0';

			delete[] new_dir_buff;
			delete[] no_directory;

			return new_path;
		}
		else
		{
			return remove_directory(path);
		}
	}

	char* remove_directory(char* path)
	{
		int len = (int)strlen(path);
		int pos = -1;

		for (int i = 0; i < len; i++)
		{
			if (path[len - i - 1] == '\\' || path[len - i - 1] == '/')
			{
				pos = len - i - 1;
				break;
			}
		}

		char* fname = new char[len - pos];
		memcpy(fname, path + pos + 1, len - pos - 1);
		fname[len - pos - 1] = '\0';
		return fname;
	}

	char* get_directory(char* path)
	{
		int len = (int)strlen(path);
		int pos = -1;

		for (int i = 0; i < len; i++)
		{
			if (path[len - i - 1] == '\\' || path[len - i - 1] == '/')
			{
				pos = len - i - 1;
				break;
			}
		}

		char* dir = new char[pos + 2];
		memcpy(dir, path, pos + 1);
		dir[pos + 1] = '\0';
		return dir;
	}

	bool is_standard_extension(char* path)
	{
		char* ext = get_extension(path);
		size_t len = strlen(ext);
		if (ext)
		{
			delete[] ext;
			ext = NULL;
		}
		return len == 4;
	}

	bool is_extension_matched(char* path_1, char* path_2)
	{
		char* ext_1 = get_extension(path_1);
		char* ext_2 = get_extension(path_2);

		bool flag = (c_str_cmpi(ext_1, ext_2) == 0);
		
		if (ext_1)
		{
			delete[] ext_1;
			ext_1 = NULL;
		}
		if (ext_2)
		{
			delete[] ext_2;
			ext_2 = NULL;
		}
		return flag;
	}
	
	char* match_extension(char* path, char* ref_path)
	{
		char* nu_path = NULL;

		if (is_extension_matched(path, ref_path))
		{
			nu_path = new char[strlen(path) + 1];
			memcpy(nu_path, path, strlen(path));
			nu_path[strlen(path)] = '\0';
		}
		else
		{
			char* ref_ext = get_extension(ref_path);

			if (is_standard_extension(path))
			{
				nu_path = replace_extension(path, ref_ext, true);
			}
			else
			{
				nu_path = new char[strlen(path) + strlen(ref_ext) + 1];
				memcpy(nu_path, path, strlen(path));
				nu_path[strlen(path)] = '\0';
				strcat(nu_path, ref_ext);
			}

			if (ref_ext)
			{
				delete[] ref_ext;
				ref_ext = NULL;
			}
		}

		return nu_path;
	}

	int n_file_from_lst(char* input_lst)
	{
		int n_file = 0;
		FILE* fp = fopen(input_lst, "rt");

		char buff[EZPC_PATH_BUFFER];
		while (fgets(buff, EZPC_PATH_BUFFER, fp) != NULL)
		{
			n_file++;
		}

		fclose(fp);

		return n_file;
	}

	char** file_list_from_lst(char* input_lst, int& n_file)
	{
		n_file = n_file_from_lst(input_lst);
		char** list = c_create_list(n_file);

		FILE* fp = fopen(input_lst, "rt");
		for (int i = 0; i < n_file; i++)
		{
			fgets(list[i], EZPC_PATH_BUFFER, fp);
			strtok(list[i], "\n");
		}
		fclose(fp);

		return list;
	}

	char** file_list_from_file(char** input_file_list, int i_file_first, int i_file_last, int& n_file)
	{
		n_file = 0;
		for (int i_f = i_file_first; i_f <= i_file_last; i_f++)
		{
			if (is_extension_matched(input_file_list[i_f], (char*)Extension_LST))
			{
				n_file += n_file_from_lst(input_file_list[i_f]);
			}
			else if (is_extension_matched(input_file_list[i_f], (char*)Extension_LAS)
				|| is_extension_matched(input_file_list[i_f], (char*)Extension_LAZ)
				|| is_extension_matched(input_file_list[i_f], (char*)Extension_E57))
			{
				n_file++;
			}
		}

		char** list = c_create_list(n_file);

		int i_list = 0;
		for (int i_f = i_file_first; i_f <= i_file_last; i_f++)
		{
			if (is_extension_matched(input_file_list[i_f], (char*)Extension_LST))
			{
				int lst_size = 0;
				char** lst_list = file_list_from_lst(input_file_list[i_f], lst_size);
				for (int j = 0; j < lst_size; j++)
				{
					memcpy(list[i_list], lst_list[j], strlen(lst_list[j]) + 1);
					i_list++;
				}
				c_delete_list(lst_list, lst_size);
			}
			else if (is_extension_matched(input_file_list[i_f], (char*)Extension_LAS)
				|| is_extension_matched(input_file_list[i_f], (char*)Extension_LAZ)
				|| is_extension_matched(input_file_list[i_f], (char*)Extension_E57))
			{
				memcpy(list[i_list], input_file_list[i_f], strlen(input_file_list[i_f]) + 1);
				i_list++;
			}
		}

		return list;
	}

	void Write_File_List(char** file_list, int n_file, char* output_file, bool overwrite_or_append)
	{
		char* ofile = match_extension(output_file, (char*)Extension_LST);
		FILE* fp = fopen(ofile, overwrite_or_append ? "wt" : "a+");

		for (int i_file = 0; i_file < n_file; i_file++)
		{
			fprintf(fp, "%s\n", file_list[i_file]);
		}

		if (ofile)
		{
			delete[] ofile;
			ofile = NULL;
		}
		fclose(fp);
	}

	long long f_count(bool* flag_list, long long n_flag, bool flag)
	{
		return f_count(flag_list, 0, n_flag - 1, flag);
	}

	long long f_count(bool* flag_list, long long i_flag_first, long long i_flag_last, bool flag)
	{
		int n = 0;
#pragma omp parallel for schedule (guided)
		for (long long i = i_flag_first; i <= i_flag_last; i++)
		{
			if (flag_list[i] == flag)
			{
#pragma omp atomic
				n++;
			}
		}
		return n;
	}

	bool* f_create(long long n_flag, bool flag)
	{
		bool* flag_list = new bool[n_flag];

#pragma omp parallel for schedule (static)
		for (long long i = 0; i < n_flag; i++)
		{
			flag_list[i] = flag;
		}

		return flag_list;
	}

	bool* f_create_cpy(bool* flag_list, long long n_flag)
	{
		bool* flag_cpy = new bool[n_flag];
#pragma omp parallel for schedule (static)
		for (long long i = 0; i < n_flag; i++)
		{
			flag_cpy[i] = flag_list[i];
		}
		return flag_cpy;
	}

	long long* f_idxlist(bool* flag_list, long long n_flag, bool flag, long long& n_idx)
	{
		n_idx = f_count(flag_list, n_flag, flag);
		long long* list = new long long[n_idx];
		long long i_idx = 0;

		for (long long i = 0; i < n_flag; i++)
		{
			if (flag_list[i] == flag)
			{
				list[i_idx] = i;
				i_idx++;
			}
		}
		return list;
	}

	bool* f_intersect(bool* flag_1, bool* flag_2, long long n_flag)
	{
		bool* flag = new bool[n_flag];
#pragma omp parallel for schedule (static)
		for (long long i = 0; i < n_flag; i++)
		{
			flag[i] = flag_1[i] && flag_2[i] ? true : false;
		}
		return flag;
	}

	bool* f_union(bool* flag_1, bool* flag_2, long long n_flag)
	{
		bool* flag = new bool[n_flag];
#pragma omp parallel for schedule (static)
		for (long long i = 0; i < n_flag; i++)
		{
			flag[i] = flag_1[i] || flag_2[i] ? true : false;
		}
		return flag;
	}

	char** c_create_list(int n_file)
	{
		char** list = new char* [n_file];
		for (long long i = 0; i < n_file; i++)
		{
			list[i] = new char[EZPC_PATH_BUFFER];
		}
		return list;
	}

	char** c_cpy_list(char** file_list, int i_first, int i_last, int& n_file)
	{
		n_file = i_last - i_first + 1;
		char** list = c_create_list(n_file);
		for (int i = 0; i < n_file; i++)
		{
			sprintf(list[i], "%s%c", file_list[i + i_first], '\0');
		}

		return list;
	}

	void c_sort_list(char** file_list, int i_first, int i_last)
	{
		std::sort(file_list + i_first, file_list + i_last + 1, c_compare_string);
	}

	bool c_compare_string(char* str_1, char* str_2)
	{
		return c_str_cmpi(str_1, str_2) < 0;
	}

	char* c_to_upper(char* str)
	{
		size_t len = strlen(str);
		char* str_up = new char[len + 1];
		for (size_t i = 0; i < len; i++)
		{
			str_up[i] = toupper(str[i]);
		}
		str_up[len] = '\0';
		return str_up;
	}

	char* c_to_lower(char* str)
	{
		size_t len = strlen(str);
		char* str_lo = new char[len + 1];
		for (size_t i = 0; i < len; i++)
		{
			str_lo[i] = tolower(str[i]);
		}
		str_lo[len] = '\0';
		return str_lo;
	}

	int c_str_cmpi(char* str1, char* str2)
	{
		char* str1_up = c_to_upper(str1);
		char* str2_up = c_to_upper(str2);

		int cmp = strcmp(str1_up, str2_up);

		delete[] str1_up;
		delete[] str2_up;

		return cmp;
	}

	void c_delete_list(char** file_list, int n_file)
	{
		if (n_file > 0)
		{
			for (int i = 0; i < n_file; i++)
			{
				if (file_list[i])
				{
					delete[] file_list[i];
					file_list[i] = NULL;
				}
			}
			delete[] file_list;
			file_list = NULL;
		}
	}

	void c_print_list(const char* before, char** file_list, const char* after, int n_file)
	{
		for (int i = 0; i < n_file; i++)
		{
			printf("%s%s%s", before, file_list[i], after);
		}
	}

	void c_print_list(char** file_list, const char* after, int n_file)
	{
		c_print_list("\0", file_list, after, n_file);
	}

	void c_print_list(const char* before, char** file_list, int n_file)
	{
		c_print_list(before, file_list, "\0", n_file);
	}

	void c_print_list(char** file_list, int n_file)
	{
		c_print_list("\0", file_list, "\0", n_file);
	}
}