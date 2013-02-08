Directly from vlfeat source except:
- comment VLEXPORT in host.h
	#if defined(VL_COMPILER_MSC)
	...
		#define VL_EXPORT //__declspec(dllimport)
	
