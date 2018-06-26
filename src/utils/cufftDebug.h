
#ifndef CUFFT_DEBUG_H
#define CUFFT_DEBUG_H
#include"debugTools.cuh"
#ifdef CUDA_ERROR_CHECK
#define CUFFT_ERROR_CHECK
#endif


#define CufftSafeCall(err) __cufftSafeCall(err, __FILE__, __LINE__)


  
const char* cufftGetErrorString(cufftResult_t err){
  switch (err) {
  case CUFFT_INVALID_PLAN:   return "CUFFT_INVALID_PLAN\n";
  case CUFFT_ALLOC_FAILED:   return "CUFFT_ALLOC_FAILED\n";
  case CUFFT_INVALID_TYPE:   return "CUFFT_INVALID_TYPE\n";
  case CUFFT_INVALID_VALUE:  return "CUFFT_INVALID_VALUE\n";
  case CUFFT_INTERNAL_ERROR: return "CUFFT_INTERNAL_ERROR\n";
  case CUFFT_EXEC_FAILED:    return "CUFFT_EXEC_FAILED\n";
  case CUFFT_SETUP_FAILED:   return "CUFFT_SETUP_FAILED\n";
  case CUFFT_INVALID_SIZE:   return "CUFFT_INVALID_SIZE\n";
  case CUFFT_UNALIGNED_DATA: return "CUFFT_UNALIGNED_DATA\n";
  case   CUFFT_INCOMPLETE_PARAMETER_LIST: return "CUFFT_INCOMPLETE_PARAMETER_LIST \n";
  case   CUFFT_INVALID_DEVICE:            return "CUFFT_INVALID_DEVICE \n";
  case   CUFFT_PARSE_ERROR:               return "CUFFT_PARSE_ERROR    \n";
  case   CUFFT_NO_WORKSPACE:              return "CUFFT_NO_WORKSPACE   \n";
  case   CUFFT_NOT_IMPLEMENTED:           return "CUFFT_NOT_IMPLEMENTED \n";
  case   CUFFT_LICENSE_ERROR:             return "CUFFT_LICENSE_ERROR  \n";
#ifdef CUFFT_NOT_SUPPORTED
  case   CUFFT_NOT_SUPPORTED:             return "CUFFT_NOT_SUPPORTED  \n";
#endif
  default: return "CUFFT Unknown error code\n";
  }
}

inline void __cufftSafeCall( cufftResult_t err, const char *file, const int line ){
  #ifdef CUFFT_ERROR_CHECK
  if ( CUFFT_SUCCESS != err )
    {
      fprintf( stderr, "cufftSafeCall() failed at %s:%i : %s - code: %i\n",
	       file, line, cufftGetErrorString( err ), err);
      exit( -1 );
    }
  #endif

  return;
}

#endif


