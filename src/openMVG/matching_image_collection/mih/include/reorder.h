#include <vector>
using namespace std;

void greedyorder (int * order, UINT8 * B, size_t N, int d, int m);
int minmax (int *C, vector<int> cbits, int * freebits, int d);
int maxmax (int *C, vector<int> cbits, int * freebits, int d);

void reorder(UINT8 * out, UINT8 * in, size_t N, int d, int * order);
