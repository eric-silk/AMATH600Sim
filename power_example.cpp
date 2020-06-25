#include <iostream>
#include "power_example.h"

int main(void)
{
  Load my_load;
  my_load.status = 1;
  
  PFData my_pfdata;
  my_pfdata.load = my_load;
  std::cout << "STATUS: " << my_pfdata.load.status << std::endl;

  return 0;
}
