// C++ program to demonstrate the use of std::max
// C++ program to demonstrate the use of std::max
#include <iostream>
#include <algorithm>
#include <array>
#include <string>

using namespace std;
int main()
{
  // Comparing ASCII values of a and b
  cout << std::max('a','b') << "\n";
  
  // Returns the first one if both the numbers
  // are same
  cout << std::max(7,7) << endl;

  cout << std::max(7, int(8.2)) << endl;
  
  std::array<std::pair<string, string>, 5> histTypes;
  histTypes[0] = {string("GenTau"), string("Generated Tau")};
  histTypes[1] = {string("VisTau"), string("Visible Tau")};
  histTypes[2] = {string("RecTau"), string("Reconstructed Tau")};
  histTypes[3] = {string("MchTau"), string("Correctly Matched Reconstructed Tau")};
  histTypes[4] = {string("WrgTau"), string("Wrongly Matched Reco Tau")};

  return 0;
} 
