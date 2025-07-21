#define HILBERT_HIGH_DIM_H

#include <cstdint>

namespace PSPT {
namespace LUT {
uint32_t deinterleave(uint32_t x) {
  x = x & 0x55555555;
  x = (x | (x >> 1)) & 0x33333333;
  x = (x | (x >> 2)) & 0x0F0F0F0F;
  x = (x | (x >> 4)) & 0x00FF00FF;
  x = (x | (x >> 8)) & 0x0000FFFF;
  return x;
}

uint32_t interleave(uint32_t x) {
  x = (x | (x << 8)) & 0x00FF00FF;
  x = (x | (x << 4)) & 0x0F0F0F0F;
  x = (x | (x << 2)) & 0x33333333;
  x = (x | (x << 1)) & 0x55555555;
  return x;
}

uint32_t prefixScan(uint32_t x) {
  x = (x >> 8) ^ x;
  x = (x >> 4) ^ x;
  x = (x >> 2) ^ x;
  x = (x >> 1) ^ x;
  return x;
}

uint32_t descan(uint32_t x) { return x ^ (x >> 1); }

void hilbertIndexToXY(uint32_t n, uint32_t i, uint32_t& x, uint32_t& y) {
  i = i << (32 - 2 * n);

  uint32_t i0 = deinterleave(i);
  uint32_t i1 = deinterleave(i >> 1);

  uint32_t t0 = (i0 | i1) ^ 0xFFFF;
  uint32_t t1 = i0 & i1;

  uint32_t prefixT0 = prefixScan(t0);
  uint32_t prefixT1 = prefixScan(t1);

  uint32_t a = (((i0 ^ 0xFFFF) & prefixT1) | (i0 & prefixT0));

  x = (a ^ i1) >> (16 - n);
  y = (a ^ i0 ^ i1) >> (16 - n);
}

uint32_t hilbertXYToIndex(uint32_t n, uint32_t x, uint32_t y) {
  x = x << (16 - n);
  y = y << (16 - n);

  uint32_t A, B, C, D;

  // Initial prefix scan round, prime with x and y
  {
    uint32_t a = x ^ y;
    uint32_t b = 0xFFFF ^ a;
    uint32_t c = 0xFFFF ^ (x | y);
    uint32_t d = x & (y ^ 0xFFFF);

    A = a | (b >> 1);
    B = (a >> 1) ^ a;

    C = ((c >> 1) ^ (b & (d >> 1))) ^ c;
    D = ((a & (c >> 1)) ^ (d >> 1)) ^ d;
  }

  {
    uint32_t a = A;
    uint32_t b = B;
    uint32_t c = C;
    uint32_t d = D;

    A = ((a & (a >> 2)) ^ (b & (b >> 2)));
    B = ((a & (b >> 2)) ^ (b & ((a ^ b) >> 2)));

    C ^= ((a & (c >> 2)) ^ (b & (d >> 2)));
    D ^= ((b & (c >> 2)) ^ ((a ^ b) & (d >> 2)));
  }

  {
    uint32_t a = A;
    uint32_t b = B;
    uint32_t c = C;
    uint32_t d = D;

    A = ((a & (a >> 4)) ^ (b & (b >> 4)));
    B = ((a & (b >> 4)) ^ (b & ((a ^ b) >> 4)));

    C ^= ((a & (c >> 4)) ^ (b & (d >> 4)));
    D ^= ((b & (c >> 4)) ^ ((a ^ b) & (d >> 4)));
  }

  // Final round and projection
  {
    uint32_t a = A;
    uint32_t b = B;
    uint32_t c = C;
    uint32_t d = D;

    C ^= ((a & (c >> 8)) ^ (b & (d >> 8)));
    D ^= ((b & (c >> 8)) ^ ((a ^ b) & (d >> 8)));
  }

  // Undo transformation prefix scan
  uint32_t a = C ^ (C >> 1);
  uint32_t b = D ^ (D >> 1);

  // Recover index bits
  uint32_t i0 = x ^ y;
  uint32_t i1 = b | (0xFFFF ^ (i0 | a));

  return ((interleave(i1) << 1) | interleave(i0)) >> (32 - 2 * n);
}

// These are multiplication tables of the alternating group A4,
// preconvolved with the mapping between Morton and Hilbert curves.
static uint8_t const mortonToHilbertTable[] = {
    48, 33, 27, 34, 47, 78, 28, 77, 66, 29, 51, 52, 65, 30, 72, 63,
    76, 95, 75, 24, 53, 54, 82, 81, 18, 3,  17, 80, 61, 4,  62, 15,
    0,  59, 71, 60, 49, 50, 86, 85, 84, 83, 5,  90, 79, 56, 6,  89,
    32, 23, 1,  94, 11, 12, 2,  93, 42, 41, 13, 14, 35, 88, 36, 31,
    92, 37, 87, 38, 91, 74, 8,  73, 46, 45, 9,  10, 7,  20, 64, 19,
    70, 25, 39, 16, 69, 26, 44, 43, 22, 55, 21, 68, 57, 40, 58, 67,
};

static uint8_t const hilbertToMortonTable[] = {
    48, 33, 35, 26, 30, 79, 77, 44, 78, 68, 64, 50, 51, 25, 29, 63,
    27, 87, 86, 74, 72, 52, 53, 89, 83, 18, 16, 1,  5,  60, 62, 15,
    0,  52, 53, 57, 59, 87, 86, 66, 61, 95, 91, 81, 80, 2,  6,  76,
    32, 2,  6,  12, 13, 95, 91, 17, 93, 41, 40, 36, 38, 10, 11, 31,
    14, 79, 77, 92, 88, 33, 35, 82, 70, 10, 11, 23, 21, 41, 40, 4,
    19, 25, 29, 47, 46, 68, 64, 34, 45, 60, 62, 71, 67, 18, 16, 49,
};

uint32_t transformCurve(uint32_t in, uint32_t bits,
                        uint8_t const* lookupTable) {
  uint32_t transform = 0;
  uint32_t out = 0;

  for (int32_t i = 3 * (bits - 1); i >= 0; i -= 3) {
    transform = lookupTable[transform | ((in >> i) & 7)];
    out = (out << 3) | (transform & 7);
    transform &= ~7;
  }

  return out;
}

uint32_t mortonToHilbert3D(uint32_t mortonIndex, uint32_t bits) {
  return transformCurve(mortonIndex, bits, mortonToHilbertTable);
}

uint32_t hilbertToMorton3D(uint32_t hilbertIndex, uint32_t bits) {
  return transformCurve(hilbertIndex, bits, hilbertToMortonTable);
}
}  // namespace LUT

namespace LUT_OPT {
static uint32_t const mortonToHilbertTable[] = {
    0x00000006, 0x20000005, 0x60000001, 0x40000007, 0xe0000001, 0xc000000c,
    0x80000005, 0xa000000e, 0x40000000, 0xa000000a, 0x6000000c, 0x8000000d,
    0x20000004, 0xc000000e, 0x00000007, 0xe0000008, 0x80000019, 0xe000001a,
    0x6000001b, 0x00000010, 0xa0000012, 0xc0000013, 0x4000001c, 0x2000001d,
    0x4000001a, 0x60000019, 0x20000018, 0x00000011, 0xa000001b, 0x8000001d,
    0xc0000019, 0xe000001e, 0x00000020, 0x60000026, 0xe000002a, 0x80000024,
    0x20000022, 0x40000023, 0xc000002c, 0xa000002d, 0x80000022, 0x60000023,
    0xa000002a, 0x40000020, 0xe0000025, 0x0000002a, 0xc000002e, 0x20000024,
    0x00000034, 0xe0000033, 0x20000032, 0xc0000038, 0x60000035, 0x80000034,
    0x40000036, 0xa000003c, 0x4000003d, 0x2000003c, 0xa000003b, 0xc000003a,
    0x60000038, 0x00000036, 0x8000003a, 0xe000003c, 0x8000004b, 0xa0000045,
    0xe0000048, 0xc0000047, 0x6000004f, 0x4000004c, 0x00000047, 0x2000004e,
    0xc000004d, 0xa000004c, 0x2000004b, 0x4000004a, 0xe000004c, 0x8000004f,
    0x00000046, 0x6000004d, 0xc0000058, 0x20000052, 0xe0000056, 0x00000051,
    0xa000005c, 0x40000056, 0x80000053, 0x60000052, 0xc000005a, 0xe000005f,
    0xa0000058, 0x80000053, 0x2000005b, 0x00000058, 0x40000059, 0x60000057,
};

static uint32_t const hilbertToMortonTable[] = {
    0x00000006, 0x20000005, 0x60000006, 0x40000000, 0xc0000007, 0xe000000c,
    0xa000000f, 0x80000002, 0xc0000001, 0x80000001, 0x00000002, 0x4000000d,
    0x6000000a, 0x2000000e, 0xa000000d, 0xe0000008, 0x60000013, 0xe000001b,
    0xc0000018, 0x4000001a, 0x0000001d, 0x80000013, 0xa0000010, 0x2000001c,
    0x60000012, 0x4000001b, 0x00000018, 0x2000001b, 0xa000001c, 0x8000001a,
    0xc0000019, 0xe000001e, 0x00000020, 0x80000027, 0xa0000024, 0x20000024,
    0x60000023, 0xe000002f, 0xc000002c, 0x4000002f, 0xa000002f, 0xe0000022,
    0x60000021, 0x20000021, 0x00000026, 0x4000002d, 0xc000002e, 0x80000026,
    0x00000034, 0x40000031, 0xc0000032, 0x80000032, 0xa0000035, 0xe000003e,
    0x6000003d, 0x20000035, 0xa0000033, 0x2000003c, 0x0000003f, 0x8000003f,
    0xc0000038, 0x4000003c, 0x6000003f, 0xe000003c, 0xc0000041, 0xe0000048,
    0xa000004b, 0x80000048, 0x0000004f, 0x20000041, 0x60000042, 0x4000004d,
    0xc0000040, 0x40000048, 0x6000004b, 0xe0000049, 0xa000004e, 0x20000048,
    0x0000004b, 0x8000004f, 0x60000052, 0x20000052, 0xa0000051, 0xe0000056,
    0xc0000051, 0x8000005d, 0x0000005e, 0x40000053, 0xa000005d, 0x8000005e,
    0xc000005d, 0xe0000053, 0x60000054, 0x4000005f, 0x0000005c, 0x20000059,
};

#if defined(_MSC_VER)
#include <intrin.h>
#define rotr64 _rotr64
#define rotl64 _rotl64
#else
inline uint64_t rotr64(uint64_t x, unsigned int n) {
  return (x >> n) | (x << (64 - n));
}
inline uint64_t rotl64(uint64_t x, unsigned int n) {
  return (x << n) | (x >> (64 - n));
}
#endif

uint32_t transformCurve(uint32_t in, uint32_t bits,
                        uint32_t const* lookupTable) {
  uint64_t x = rotr64(in, 3 * bits);

  for (uint32_t i = 0; i < bits; ++i) {
    x = rotl64(x, 3);
    x ^= lookupTable[uint32_t(x)];
  }

  return x >> 29;
}

uint32_t mortonToHilbert3D(uint32_t mortonIndex, uint32_t bits) {
  return transformCurve(mortonIndex, bits, mortonToHilbertTable);
}

uint32_t hilbertToMorton3D(uint32_t hilbertIndex, uint32_t bits) {
  return transformCurve(hilbertIndex, bits, hilbertToMortonTable);
}
}  // namespace LUT_OPT
}  // namespace PSPT
