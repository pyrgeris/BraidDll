#ifndef _NATIVELIB_H_
#define _NATIVELIB_H_

#ifndef MYAPI
#define MYAPI
#endif

#ifdef __cplusplus
extern "C" {
#endif

	MYAPI __declspec(dllexport) void print_line(const char* str);
	MYAPI __declspec(dllexport) void braidcorrelation(const char* inputfile, int length, int sequences, int maxb, int interval);

#ifdef __cplusplus
}
#endif

#endif // _NATIVELIB_H_