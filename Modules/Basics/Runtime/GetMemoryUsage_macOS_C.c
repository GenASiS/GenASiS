// https://stackoverflow.com/questions/18389581/memory-used-by-a-process-under-mac-os-x
// https://stackoverflow.com/questions/60751839/how-i-can-get-peak-memory-of-a-process-on-mac-os

#ifdef __APPLE__
#include <mach/mach.h>
#include <stdio.h>
#endif

void GetMemoryUsage_macOS_C ( double *hwm_kb, double *rss_kb )
{
  
#ifdef __APPLE__
  kern_return_t error;
  mach_msg_type_number_t outCount;
  mach_task_basic_info_data_t taskinfo;

  taskinfo.resident_size_max = 0;
  taskinfo.resident_size     = 0;
  outCount = MACH_TASK_BASIC_INFO_COUNT;
  error = task_info ( mach_task_self(), MACH_TASK_BASIC_INFO,
		      (task_info_t)&taskinfo, &outCount );
  if (error == KERN_SUCCESS) {
    // type is mach_vm_size_t
    *hwm_kb = ((double)((unsigned long long)taskinfo.resident_size_max)) / 1024;
    *rss_kb = ((double)((unsigned long long)taskinfo.resident_size)) / 1024;
  } else {
    printf("error %d\n", (int)error);
  }
#endif
  
}
