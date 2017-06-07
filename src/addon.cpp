#include <node.h>

#include "MarchingCubes.h"
#include "LookUpTable.h"
#include "csg.h"
#include "PerlinNoise/PerlinNoise.hpp"

#define PI 3.14159265358979323846
#define SIZE 16
#define BUFFER 1
#define SIZE_BUFFERED (SIZE + (BUFFER * 2))
#define RADIUS 50

siv::PerlinNoise noise;
const float elevationNoiseFrequency = 0.05;
const unsigned int elevationNoiseOctaves = 8;

const float moistureNoiseFrequency = 0.04;
const unsigned int moistureNoiseOctaves =  6;

typedef struct {
  real x, y, z;
} Vector3;
typedef struct {
  int x, y, z;
} IntVector3;
typedef struct {
  real x, y;
} Vector2;

enum class Biome : unsigned int {
  // Features
  OCEAN = 0x44447a,
  // OCEAN = 0x000000,
  // COAST = 0x33335a,
  COAST = 0x333333,
  LAKESHORE = 0x225588,
  LAKE = 0x336699,
  RIVER = 0x225588,
  MARSH = 0x2f6666,
  // ICE = 0x99ffff,
  ICE = 0x99dddd,
  // BEACH = 0xa09077,
  BEACH = 0xa0b077,
  ROAD1 = 0x442211,
  ROAD2 = 0x553322,
  ROAD3 = 0x664433,
  BRIDGE = 0x686860,
  LAVA = 0xcc3333,

  // Terrain
  SNOW = 0xffffff,
  TUNDRA = 0xbbbbaa,
  BARE = 0x888888,
  SCORCHED = 0x555555,
  TAIGA = 0x99aa77,
  SHRUBLAND = 0x889977,
  TEMPERATE_DESERT = 0xc9d29b,
  TEMPERATE_RAIN_FOREST = 0x448855,
  TEMPERATE_DECIDUOUS_FOREST = 0x679459,
  GRASSLAND = 0x88aa55,
  SUBTROPICAL_DESERT = 0xd2b98b,
  TROPICAL_RAIN_FOREST = 0x337755,
  TROPICAL_SEASONAL_FOREST = 0x559944,
  MAGMA = 0xff3333,
};

float _sum(const Vector3 &v) {
  return v.x + v.y + v.z;
}
float _length(const Vector3 &v) {
  return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}
float _dot(const Vector3 &a, const Vector3 &b) {
  return a.x*b.x + a.y*b.y + a.z*b.z;
}
float _lenSq(const Vector3 &v) {
  return v.x*v.x + v.y*v.y + v.z*v.z;
}
template <typename T>
T _clamp(const T& n, const T& lower, const T& upper) {
  return std::max(lower, std::min(n, upper));
}
float _angleTo(const Vector3 &a, const Vector3 &b) {
  float dot = _dot(a, b);
  float lenSq1 = _lenSq(a);
  float lenSq2 = _lenSq(b);
  float angle = std::acos(_clamp(dot / std::sqrt(lenSq1 * lenSq2), (float)(-1), (float)(1)));
  return angle;
}
Vector3 _multiply(const Vector3 &a, const Vector3 &b) {
  return Vector3{
    a.x * b.x,
    a.y * b.y,
    a.z * b.z
  };
}
float _getValue(
  const Vector3 &vector,
  const IntVector3 &origin,
  const float noiseFrequency,
  const unsigned int noiseOctaves,
  const float minValue,
  const float factor,
  const float lengthPow,
  const float etherFactor
) {
  Vector3 localVector{
    vector.x - (SIZE_BUFFERED / 2),
    vector.y - (SIZE_BUFFERED / 2),
    vector.z - (SIZE_BUFFERED / 2)
  };
  Vector3 absoluteVector{
    localVector.x + (origin.x * SIZE),
    localVector.y + (origin.y * SIZE),
    localVector.z + (origin.z * SIZE)
  };

  float ether = std::abs(noise.octaveNoise(
    (absoluteVector.x * noiseFrequency),
    (absoluteVector.y * noiseFrequency),
    (absoluteVector.z * noiseFrequency),
    noiseOctaves
  ));
  float length = _length(absoluteVector);
  float lengthValue = -minValue + std::pow(length, lengthPow) / factor;
  float value = lengthValue - (ether * etherFactor);
  return value;
}

v8::Local<v8::Value> DoGenMarchCubes(
  v8::Isolate* isolate,
  unsigned int seedNumber,
  v8::Local<v8::Array> originValue,
  const float noiseFrequency,
  const unsigned int noiseOctaves,
  const float minValue,
  const float factor,
  const float lengthPow,
  const float etherFactor
) {
  v8::EscapableHandleScope scope(isolate);

  auto positionsKey = v8::String::NewFromUtf8(isolate, "positions");
  auto normalsKey = v8::String::NewFromUtf8(isolate, "normals");
  auto etherKey = v8::String::NewFromUtf8(isolate, "ether");

  IntVector3 originVector{
    originValue->Get(0)->Int32Value(),
    originValue->Get(1)->Int32Value(),
    originValue->Get(2)->Int32Value()
  };

  // begin mc
  // Init data
  noise.reseed(seedNumber);
  MarchingCubes mc(SIZE_BUFFERED, SIZE_BUFFERED, SIZE_BUFFERED);

  // Fills data structure
  for (unsigned int i = 0; i < SIZE_BUFFERED; i++) {
    for (unsigned int j = 0; j < SIZE_BUFFERED; j++) {
      for (unsigned int k = 0; k < SIZE_BUFFERED; k++) {
        float v = _getValue(
          Vector3{
            (float)i,
            (float)j,
            (float)k
          },
          originVector,
          noiseFrequency,
          noiseOctaves,
          minValue,
          factor,
          lengthPow,
          etherFactor
        );
        mc.set_data(v, i, j, k);
      }
    }
  }

  mc.run();
  // end mc

  // construct result
  unsigned int numTrigs = mc.ntrigs();
  unsigned int numVerts = numTrigs * 3;
  unsigned int numPositions = numVerts * 3;
  unsigned int numNormals = numVerts * 3;
  v8::Local<v8::ArrayBuffer> positionsBuffer = v8::ArrayBuffer::New(isolate, numPositions * 4);
  v8::Local<v8::Float32Array> positions = v8::Float32Array::New(positionsBuffer, 0, numPositions);
  v8::Local<v8::ArrayBuffer> normalsBuffer = v8::ArrayBuffer::New(isolate, numNormals * 4);
  v8::Local<v8::Float32Array> normals = v8::Float32Array::New(normalsBuffer, 0, numNormals);

  auto triangles = mc.triangles();
  auto vertices = mc.vertices();
  unsigned int positionIndex = 0;
  for (unsigned int i = 0; i < numTrigs; i++) {
    const Triangle &triangle = triangles[i];
    const Vertex &a = vertices[triangle.v1];
    const Vertex &b = vertices[triangle.v2];
    const Vertex &c = vertices[triangle.v3];

    const Vertex da{
      a.x - BUFFER,
      a.y - BUFFER,
      a.z - BUFFER
    };
    const Vertex db{
      b.x - BUFFER,
      b.y - BUFFER,
      b.z - BUFFER
    };
    const Vertex dc{
      c.x - BUFFER,
      c.y - BUFFER,
      c.z - BUFFER
    };
    if (
        da.x >= 0 && da.x <= SIZE && da.y >= 0 && da.y <= SIZE && da.z >= 0 && da.z <= SIZE &&
        db.x >= 0 && db.x <= SIZE && db.y >= 0 && db.y <= SIZE && db.z >= 0 && db.z <= SIZE &&
        dc.x >= 0 && dc.x <= SIZE && dc.y >= 0 && dc.y <= SIZE && dc.z >= 0 && dc.z <= SIZE
    ) {
      unsigned int baseIndex = positionIndex * 3 * 3;
      positions->Set(baseIndex + 0, v8::Number::New(isolate, a.x - (SIZE_BUFFERED / 2)));
      positions->Set(baseIndex + 1, v8::Number::New(isolate, a.y - (SIZE_BUFFERED / 2)));
      positions->Set(baseIndex + 2, v8::Number::New(isolate, a.z - (SIZE_BUFFERED / 2)));
      positions->Set(baseIndex + 3, v8::Number::New(isolate, b.x - (SIZE_BUFFERED / 2)));
      positions->Set(baseIndex + 4, v8::Number::New(isolate, b.y - (SIZE_BUFFERED / 2)));
      positions->Set(baseIndex + 5, v8::Number::New(isolate, b.z - (SIZE_BUFFERED / 2)));
      positions->Set(baseIndex + 6, v8::Number::New(isolate, c.x - (SIZE_BUFFERED / 2)));
      positions->Set(baseIndex + 7, v8::Number::New(isolate, c.y - (SIZE_BUFFERED / 2)));
      positions->Set(baseIndex + 8, v8::Number::New(isolate, c.z - (SIZE_BUFFERED / 2)));

      normals->Set(baseIndex + 0, v8::Number::New(isolate, a.nx));
      normals->Set(baseIndex + 1, v8::Number::New(isolate, a.ny));
      normals->Set(baseIndex + 2, v8::Number::New(isolate, a.nz));
      normals->Set(baseIndex + 3, v8::Number::New(isolate, b.nx));
      normals->Set(baseIndex + 4, v8::Number::New(isolate, b.ny));
      normals->Set(baseIndex + 5, v8::Number::New(isolate, b.nz));
      normals->Set(baseIndex + 6, v8::Number::New(isolate, c.nx));
      normals->Set(baseIndex + 7, v8::Number::New(isolate, c.ny));
      normals->Set(baseIndex + 8, v8::Number::New(isolate, c.nz));

      positionIndex++;
    }
  }

  const unsigned int numEthers = SIZE_BUFFERED * SIZE_BUFFERED * SIZE_BUFFERED;
  v8::Local<v8::ArrayBuffer> ethersBuffer = v8::ArrayBuffer::New(isolate, numEthers * 4);
  v8::Local<v8::Float32Array> ethers = v8::Float32Array::New(ethersBuffer, 0, numEthers);
  for (int i = 0; i < SIZE_BUFFERED; i++) {
    for (int j = 0; j < SIZE_BUFFERED; j++) {
      for (int k = 0; k < SIZE_BUFFERED; k++) {
        const unsigned int etherIndex = i + (j * SIZE_BUFFERED) + (k * SIZE_BUFFERED * SIZE_BUFFERED);
        const float v = mc.get_data(ivec3(i, j, k));
        ethers->Set(etherIndex, v8::Number::New(isolate, v));
      }
    }
  }

  v8::Local<v8::Object> result = v8::Object::New(isolate);
  result->Set(positionsKey, v8::Float32Array::New(positionsBuffer, 0, positionIndex * 3 * 3));
  result->Set(normalsKey, v8::Float32Array::New(normalsBuffer, 0, positionIndex * 3 * 3));
  result->Set(etherKey, ethers);

  return scope.Escape(result);
}

v8::Local<v8::Value> DoReMarchCubes(
  v8::Isolate* isolate,
  v8::Local<v8::Float32Array> etherValue
) {
  v8::EscapableHandleScope scope(isolate);

  auto positionsKey = v8::String::NewFromUtf8(isolate, "positions");
  auto normalsKey = v8::String::NewFromUtf8(isolate, "normals");
  auto etherKey = v8::String::NewFromUtf8(isolate, "ether");

  // begin mc
  // Init data
  MarchingCubes mc(SIZE_BUFFERED, SIZE_BUFFERED, SIZE_BUFFERED);

  // Fills data structure
  for (unsigned int i = 0; i < SIZE_BUFFERED; i++) {
    for (unsigned int j = 0; j < SIZE_BUFFERED; j++) {
      for (unsigned int k = 0; k < SIZE_BUFFERED; k++) {
        float v = (float)etherValue->Get(i + (j * SIZE_BUFFERED) + (k * SIZE_BUFFERED * SIZE_BUFFERED))->NumberValue();
        mc.set_data(v, i, j, k);
      }
    }
  }

  mc.run();
  // end mc

  // construct result
  unsigned int numTrigs = mc.ntrigs();
  unsigned int numVerts = numTrigs * 3;
  unsigned int numPositions = numVerts * 3;
  unsigned int numNormals = numVerts * 3;
  v8::Local<v8::ArrayBuffer> positionsBuffer = v8::ArrayBuffer::New(isolate, numPositions * 4);
  v8::Local<v8::Float32Array> positions = v8::Float32Array::New(positionsBuffer, 0, numPositions);
  v8::Local<v8::ArrayBuffer> normalsBuffer = v8::ArrayBuffer::New(isolate, numNormals * 4);
  v8::Local<v8::Float32Array> normals = v8::Float32Array::New(normalsBuffer, 0, numNormals);

  auto triangles = mc.triangles();
  auto vertices = mc.vertices();
  unsigned int positionIndex = 0;
  for (unsigned int i = 0; i < numTrigs; i++) {
    const Triangle &triangle = triangles[i];
    const Vertex &a = vertices[triangle.v1];
    const Vertex &b = vertices[triangle.v2];
    const Vertex &c = vertices[triangle.v3];

    const Vertex da{
      a.x - BUFFER,
      a.y - BUFFER,
      a.z - BUFFER
    };
    const Vertex db{
      b.x - BUFFER,
      b.y - BUFFER,
      b.z - BUFFER
    };
    const Vertex dc{
      c.x - BUFFER,
      c.y - BUFFER,
      c.z - BUFFER
    };
    if (
        da.x >= 0 && da.x <= SIZE && da.y >= 0 && da.y <= SIZE && da.z >= 0 && da.z <= SIZE &&
        db.x >= 0 && db.x <= SIZE && db.y >= 0 && db.y <= SIZE && db.z >= 0 && db.z <= SIZE &&
        dc.x >= 0 && dc.x <= SIZE && dc.y >= 0 && dc.y <= SIZE && dc.z >= 0 && dc.z <= SIZE
    ) {
      unsigned int baseIndex = positionIndex * 3 * 3;
      positions->Set(baseIndex + 0, v8::Number::New(isolate, a.x - (SIZE_BUFFERED / 2)));
      positions->Set(baseIndex + 1, v8::Number::New(isolate, a.y - (SIZE_BUFFERED / 2)));
      positions->Set(baseIndex + 2, v8::Number::New(isolate, a.z - (SIZE_BUFFERED / 2)));
      positions->Set(baseIndex + 3, v8::Number::New(isolate, b.x - (SIZE_BUFFERED / 2)));
      positions->Set(baseIndex + 4, v8::Number::New(isolate, b.y - (SIZE_BUFFERED / 2)));
      positions->Set(baseIndex + 5, v8::Number::New(isolate, b.z - (SIZE_BUFFERED / 2)));
      positions->Set(baseIndex + 6, v8::Number::New(isolate, c.x - (SIZE_BUFFERED / 2)));
      positions->Set(baseIndex + 7, v8::Number::New(isolate, c.y - (SIZE_BUFFERED / 2)));
      positions->Set(baseIndex + 8, v8::Number::New(isolate, c.z - (SIZE_BUFFERED / 2)));

      normals->Set(baseIndex + 0, v8::Number::New(isolate, a.nx));
      normals->Set(baseIndex + 1, v8::Number::New(isolate, a.ny));
      normals->Set(baseIndex + 2, v8::Number::New(isolate, a.nz));
      normals->Set(baseIndex + 3, v8::Number::New(isolate, b.nx));
      normals->Set(baseIndex + 4, v8::Number::New(isolate, b.ny));
      normals->Set(baseIndex + 5, v8::Number::New(isolate, b.nz));
      normals->Set(baseIndex + 6, v8::Number::New(isolate, c.nx));
      normals->Set(baseIndex + 7, v8::Number::New(isolate, c.ny));
      normals->Set(baseIndex + 8, v8::Number::New(isolate, c.nz));

      positionIndex++;
    }
  }

  const unsigned int numEthers = SIZE_BUFFERED * SIZE_BUFFERED * SIZE_BUFFERED;
  v8::Local<v8::ArrayBuffer> ethersBuffer = v8::ArrayBuffer::New(isolate, numEthers * 4);
  v8::Local<v8::Float32Array> ethers = v8::Float32Array::New(ethersBuffer, 0, numEthers);
  for (int i = 0; i < SIZE_BUFFERED; i++) {
    for (int j = 0; j < SIZE_BUFFERED; j++) {
      for (int k = 0; k < SIZE_BUFFERED; k++) {
        const unsigned int etherIndex = i + (j * SIZE_BUFFERED) + (k * SIZE_BUFFERED * SIZE_BUFFERED);
        float v = mc.get_data(ivec3(i, j, k));
        ethers->Set(etherIndex, v8::Number::New(isolate, v));
      }
    }
  }

  v8::Local<v8::Object> result = v8::Object::New(isolate);
  result->Set(positionsKey, v8::Float32Array::New(positionsBuffer, 0, positionIndex * 3 * 3));
  result->Set(normalsKey, v8::Float32Array::New(normalsBuffer, 0, positionIndex * 3 * 3));
  result->Set(etherKey, ethers);

  return scope.Escape(result);
}

v8::Local<v8::Value> DoHolesMarchCubes(
  v8::Isolate* isolate,
  v8::Local<v8::Array> originValue,
  v8::Local<v8::Float32Array> etherValue,
  v8::Local<v8::Int32Array> holesValue
) {
  v8::EscapableHandleScope scope(isolate);

  auto positionsKey = v8::String::NewFromUtf8(isolate, "positions");
  auto normalsKey = v8::String::NewFromUtf8(isolate, "normals");
  auto etherKey = v8::String::NewFromUtf8(isolate, "ether");

  // generate heightmap
  IntVector3 originVector{
    originValue->Get(0)->Int32Value(),
    originValue->Get(1)->Int32Value(),
    originValue->Get(2)->Int32Value()
  };

  // begin mc
  // Init data
  MarchingCubes mc(SIZE_BUFFERED, SIZE_BUFFERED, SIZE_BUFFERED);

  // Fills data structure
  for (unsigned int i = 0; i < SIZE_BUFFERED; i++) {
    for (unsigned int j = 0; j < SIZE_BUFFERED; j++) {
      for (unsigned int k = 0; k < SIZE_BUFFERED; k++) {
        float v = (float)etherValue->Get(i + (j * SIZE_BUFFERED) + (k * SIZE_BUFFERED * SIZE_BUFFERED))->NumberValue();
        mc.set_data(v, i, j, k);
      }
    }
  }

  // Adjust for holes
  unsigned int numHoles = holesValue->Length() / 3;
  for (unsigned int h = 0; h < numHoles; h++) {
    unsigned int holeIndexBase = h * 3;
    int x = holesValue->Get(holeIndexBase + 0)->Int32Value();
    int y = holesValue->Get(holeIndexBase + 1)->Int32Value();
    int z = holesValue->Get(holeIndexBase + 2)->Int32Value();

    for (int i = -2; i <= 2; i++) {
      int dx = BUFFER + (SIZE / 2) + x + i - (originVector.x * SIZE);

      if (dx >= 0 && dx < SIZE_BUFFERED) {
        for (int j = -2; j <= 2; j++) {
          int dy = BUFFER + (SIZE / 2) + y + j - (originVector.y * SIZE);

          if (dy >= 0 && dy < SIZE_BUFFERED) {
            for (int k = -2; k <= 2; k++) {
              int dz = BUFFER + (SIZE / 2) + z + k - (originVector.z * SIZE);

              if (dz >= 0 && dz < SIZE_BUFFERED) {
                float distance = std::sqrt((i * i) + (j * j) + (k * k));
                float distanceFactor = distance / std::sqrt(2*2 + 2*2 + 2*2);
                float valueFactor = std::pow(1 - distanceFactor, 3);
                float v = valueFactor * 500;

                mc.set_data(
                  mc.get_data(ivec3(dx, dy, dz)) + v,
                  dx,
                  dy,
                  dz
                );
              }
            }
          }
        }
      }
    }
  }

  mc.run();
  // end mc

  // construct result
  unsigned int numTrigs = mc.ntrigs();
  unsigned int numVerts = numTrigs * 3;
  unsigned int numPositions = numVerts * 3;
  unsigned int numNormals = numVerts * 3;
  v8::Local<v8::ArrayBuffer> positionsBuffer = v8::ArrayBuffer::New(isolate, numPositions * 4);
  v8::Local<v8::Float32Array> positions = v8::Float32Array::New(positionsBuffer, 0, numPositions);
  v8::Local<v8::ArrayBuffer> normalsBuffer = v8::ArrayBuffer::New(isolate, numNormals * 4);
  v8::Local<v8::Float32Array> normals = v8::Float32Array::New(normalsBuffer, 0, numNormals);

  auto triangles = mc.triangles();
  auto vertices = mc.vertices();
  unsigned int index = 0;
  for (unsigned int i = 0; i < numTrigs; i++) {
    const Triangle &triangle = triangles[i];
    const Vertex &a = vertices[triangle.v1];
    const Vertex &b = vertices[triangle.v2];
    const Vertex &c = vertices[triangle.v3];

    const Vertex da{
      a.x - BUFFER,
      a.y - BUFFER,
      a.z - BUFFER
    };
    const Vertex db{
      b.x - BUFFER,
      b.y - BUFFER,
      b.z - BUFFER
    };
    const Vertex dc{
      c.x - BUFFER,
      c.y - BUFFER,
      c.z - BUFFER
    };
    if (
        da.x >= 0 && da.x <= SIZE && da.y >= 0 && da.y <= SIZE && da.z >= 0 && da.z <= SIZE &&
        db.x >= 0 && db.x <= SIZE && db.y >= 0 && db.y <= SIZE && db.z >= 0 && db.z <= SIZE &&
        dc.x >= 0 && dc.x <= SIZE && dc.y >= 0 && dc.y <= SIZE && dc.z >= 0 && dc.z <= SIZE
    ) {
      unsigned int baseIndex = index * 3 * 3;
      positions->Set(baseIndex + 0, v8::Number::New(isolate, a.x - (SIZE_BUFFERED / 2)));
      positions->Set(baseIndex + 1, v8::Number::New(isolate, a.y - (SIZE_BUFFERED / 2)));
      positions->Set(baseIndex + 2, v8::Number::New(isolate, a.z - (SIZE_BUFFERED / 2)));
      positions->Set(baseIndex + 3, v8::Number::New(isolate, b.x - (SIZE_BUFFERED / 2)));
      positions->Set(baseIndex + 4, v8::Number::New(isolate, b.y - (SIZE_BUFFERED / 2)));
      positions->Set(baseIndex + 5, v8::Number::New(isolate, b.z - (SIZE_BUFFERED / 2)));
      positions->Set(baseIndex + 6, v8::Number::New(isolate, c.x - (SIZE_BUFFERED / 2)));
      positions->Set(baseIndex + 7, v8::Number::New(isolate, c.y - (SIZE_BUFFERED / 2)));
      positions->Set(baseIndex + 8, v8::Number::New(isolate, c.z - (SIZE_BUFFERED / 2)));

      normals->Set(baseIndex + 0, v8::Number::New(isolate, a.nx));
      normals->Set(baseIndex + 1, v8::Number::New(isolate, a.ny));
      normals->Set(baseIndex + 2, v8::Number::New(isolate, a.nz));
      normals->Set(baseIndex + 3, v8::Number::New(isolate, b.nx));
      normals->Set(baseIndex + 4, v8::Number::New(isolate, b.ny));
      normals->Set(baseIndex + 5, v8::Number::New(isolate, b.nz));
      normals->Set(baseIndex + 6, v8::Number::New(isolate, c.nx));
      normals->Set(baseIndex + 7, v8::Number::New(isolate, c.ny));
      normals->Set(baseIndex + 8, v8::Number::New(isolate, c.nz));

      index++;
    }
  }

  const unsigned int numEthers = SIZE_BUFFERED * SIZE_BUFFERED * SIZE_BUFFERED;
  v8::Local<v8::ArrayBuffer> ethersBuffer = v8::ArrayBuffer::New(isolate, numEthers * 4);
  v8::Local<v8::Float32Array> ethers = v8::Float32Array::New(ethersBuffer, 0, numEthers);
  for (int i = 0; i < SIZE_BUFFERED; i++) {
    for (int j = 0; j < SIZE_BUFFERED; j++) {
      for (int k = 0; k < SIZE_BUFFERED; k++) {
        const unsigned int etherIndex = i + (j * SIZE_BUFFERED) + (k * SIZE_BUFFERED * SIZE_BUFFERED);
        float v = mc.get_data(ivec3(i, j, k));
        ethers->Set(etherIndex, v8::Number::New(isolate, v));
      }
    }
  }

  v8::Local<v8::Object> result = v8::Object::New(isolate);
  result->Set(positionsKey, v8::Float32Array::New(positionsBuffer, 0, index * 3 * 3));
  result->Set(normalsKey, v8::Float32Array::New(normalsBuffer, 0, index * 3 * 3));
  result->Set(etherKey, ethers);

  return scope.Escape(result);
}

Biome _getBiome(float elevation, float moisture) {
  /* if (coast) {
    return Biome::BEACH;
  } else if (ocean) {
    return Biome::OCEAN;
  } else if (p.water) {
    if (elevation < (size * 0.1)) { return 'MARSH'; }
    if (elevation > (size * 0.25)) { return 'ICE'; }
    return Biome::LAKE;
  } else if (lava > 2) {
    return Biome::MAGMA;
  } else */if (elevation > (RADIUS * 0.5)) {
    if (moisture > 0.50) { return Biome::SNOW; }
    else if (moisture > 0.33) { return Biome::TUNDRA; }
    else if (moisture > 0.16) { return Biome::BARE; }
    else { return Biome::SCORCHED; }
  } else if (elevation > (RADIUS * 0.4)) {
    if (moisture > 0.66) { return Biome::TAIGA; }
    else if (moisture > 0.33) { return Biome::SHRUBLAND; }
    else { return Biome::TEMPERATE_DESERT; }
  } else if (elevation > (RADIUS * 0.15)) {
    if (moisture > 0.8) { return Biome::TEMPERATE_RAIN_FOREST; }
    else if (moisture > 0.5) { return Biome::TEMPERATE_DECIDUOUS_FOREST; }
    else if (moisture > 0.1) { return Biome::GRASSLAND; }
    else { return Biome::TEMPERATE_DESERT; }
  } else {
    if (moisture > 0.66) { return Biome::TROPICAL_RAIN_FOREST; }
    else if (moisture > 0.33) { return Biome::TROPICAL_SEASONAL_FOREST; }
    else if (moisture > 0.16) { return Biome::GRASSLAND; }
    else { return Biome::SUBTROPICAL_DESERT; }
  }
};
unsigned int _getBiomeColor(float elevation, float moisture) {
  auto biome = _getBiome(elevation, moisture);
  auto biomeColor = (unsigned int)biome;
  return biomeColor;
}
void GenMarchCubes(const v8::FunctionCallbackInfo<v8::Value>& args) {
  v8::Isolate* isolate = args.GetIsolate();
  v8::EscapableHandleScope scope(isolate);

  auto seedKey = v8::String::NewFromUtf8(isolate, "seed");
  auto originKey = v8::String::NewFromUtf8(isolate, "origin");
  auto landKey = v8::String::NewFromUtf8(isolate, "land");
  auto waterKey = v8::String::NewFromUtf8(isolate, "water");
  auto metadataKey = v8::String::NewFromUtf8(isolate, "metadata");
  auto positionsKey = v8::String::NewFromUtf8(isolate, "positions");
  auto colorsKey = v8::String::NewFromUtf8(isolate, "colors");
  auto moistureEtherKey = v8::String::NewFromUtf8(isolate, "moistureEther");

  // Check the number of arguments passed.
  if (args.Length() < 1) {
    // Throw an Error that is passed back to JavaScript
    isolate->ThrowException(v8::Exception::TypeError(
        v8::String::NewFromUtf8(isolate, "Wrong number of arguments")));
    return;
  }
  // Check the argument types
  if (!args[0]->IsObject()) {
    isolate->ThrowException(v8::Exception::TypeError(
        v8::String::NewFromUtf8(isolate, "Wrong arguments")));
    return;
  }

  v8::Local<v8::Object> opts = args[0]->ToObject();
  v8::Local<v8::Value> seed = opts->Get(seedKey);
  v8::Local<v8::Value> origin = opts->Get(originKey);
  if (!(seed->IsNumber() && origin->IsArray())) {
    isolate->ThrowException(v8::Exception::TypeError(
        v8::String::NewFromUtf8(isolate, "Invalid options")));
    return;
  }

  v8::Local<v8::Array> originValue = origin.As<v8::Array>();
  if (originValue->Length() != 3) {
    isolate->ThrowException(v8::Exception::TypeError(
        v8::String::NewFromUtf8(isolate, "Invalid options")));
    return;
  }

  IntVector3 originVector{
    originValue->Get(0)->Int32Value(),
    originValue->Get(1)->Int32Value(),
    originValue->Get(2)->Int32Value()
  };
  v8::Local<v8::Number> seedValue = seed.As<v8::Number>();
  unsigned int landSeedNumber = seedValue->Uint32Value() + 0;
  unsigned int waterSeedNumber = seedValue->Uint32Value() + 1;

  v8::Local<v8::Object> landResult = DoGenMarchCubes(
    isolate,
    landSeedNumber,
    originValue,
    elevationNoiseFrequency,
    elevationNoiseOctaves,
    0,
    1,
    2.25, // lengthPow
    2000 // etherFactor
  ).As<v8::Object>();

  std::default_random_engine generator(seedValue->Uint32Value() + 1);
  std::uniform_real_distribution<float> distribution(0.0f, 1.0f);
  const float moistureMultiplier = 0.2 + 0.8 * distribution(generator);
  const float moistureMinValue = 4000 * moistureMultiplier;
  const float moistureFactor = 0.2 * moistureMultiplier;
  const float mosistureLengthPow = 2.25;
  const float moistureEtherFactor = 2000;
  v8::Local<v8::Object> waterResult = DoGenMarchCubes(
    isolate,
    waterSeedNumber,
    originValue,
    moistureNoiseFrequency,
    moistureNoiseOctaves,
    moistureMinValue,
    moistureFactor,
    mosistureLengthPow,
    moistureEtherFactor
  ).As<v8::Object>();

  // begin latch moistures
  noise.reseed(waterSeedNumber);

  const unsigned int numMoistureEthers = SIZE_BUFFERED * SIZE_BUFFERED * SIZE_BUFFERED;
  v8::Local<v8::ArrayBuffer> moistureEthersBuffer = v8::ArrayBuffer::New(isolate, numMoistureEthers * 4);
  v8::Local<v8::Float32Array> moistureEthers = v8::Float32Array::New(moistureEthersBuffer, 0, numMoistureEthers);
  for (int i = 0; i < SIZE_BUFFERED; i++) {
    for (int j = 0; j < SIZE_BUFFERED; j++) {
      for (int k = 0; k < SIZE_BUFFERED; k++) {
        const unsigned int moistureEtherIndex = i + (j * SIZE_BUFFERED) + (k * SIZE_BUFFERED * SIZE_BUFFERED);
        IntVector3 localVector{
          i - (SIZE_BUFFERED / 2),
          j - (SIZE_BUFFERED / 2),
          k - (SIZE_BUFFERED / 2)
        };
        Vector3 absoluteVector{
          localVector.x + (originVector.x * SIZE),
          localVector.y + (originVector.y * SIZE),
          localVector.z + (originVector.z * SIZE)
        };
        const float v = noise.octaveNoise(
          (absoluteVector.x * moistureNoiseFrequency),
          (absoluteVector.y * moistureNoiseFrequency),
          (absoluteVector.z * moistureNoiseFrequency),
          moistureNoiseOctaves
        );
        moistureEthers->Set(moistureEtherIndex, v8::Number::New(isolate, v));
      }
    }
  }

  v8::Local<v8::Object> metadataResult = v8::Object::New(isolate);
  metadataResult->Set(moistureEtherKey, moistureEthers);
  // end latch moistures

  // begin add land colors
  v8::Local<v8::Float32Array> positions = landResult->Get(positionsKey).As<v8::Float32Array>();
  unsigned int numPositions = positions->Length();
  unsigned int numTriangles = numPositions / 3;
  v8::Local<v8::ArrayBuffer> colorsBuffer = v8::ArrayBuffer::New(isolate, numPositions * 4);
  v8::Local<v8::Float32Array> colors = v8::Float32Array::New(colorsBuffer, 0, numPositions);
  for (unsigned int i = 0; i < numTriangles; i++) {
    unsigned int triangleBaseIndex = i * 3 * 3;

    Vector3 pa{
      (float)positions->Get(triangleBaseIndex + 0)->NumberValue(),
      (float)positions->Get(triangleBaseIndex + 1)->NumberValue(),
      (float)positions->Get(triangleBaseIndex + 2)->NumberValue()
    };
    Vector3 pb{
      (float)positions->Get(triangleBaseIndex + 3)->NumberValue(),
      (float)positions->Get(triangleBaseIndex + 4)->NumberValue(),
      (float)positions->Get(triangleBaseIndex + 5)->NumberValue()
    };
    Vector3 pc{
      (float)positions->Get(triangleBaseIndex + 6)->NumberValue(),
      (float)positions->Get(triangleBaseIndex + 7)->NumberValue(),
      (float)positions->Get(triangleBaseIndex + 8)->NumberValue()
    };
    IntVector3 centerVector{
      (int)((pa.x + pb.x + pc.x) / 3),
      (int)((pa.y + pb.y + pc.y) / 3),
      (int)((pa.z + pb.z + pc.z) / 3)
    };
    IntVector3 localVector{
      centerVector.x + (SIZE_BUFFERED / 2),
      centerVector.y + (SIZE_BUFFERED / 2),
      centerVector.z + (SIZE_BUFFERED / 2)
    };
    Vector3 absoluteVector{
      centerVector.x + (originVector.x * SIZE),
      centerVector.y + (originVector.y * SIZE),
      centerVector.z + (originVector.z * SIZE)
    };
    float elevation = std::sqrt(
      absoluteVector.x * absoluteVector.x +
      absoluteVector.y * absoluteVector.y +
      absoluteVector.z * absoluteVector.z
    );
    float moisture = std::abs(
      (float)moistureEthers->Get(
        localVector.x +
        (localVector.y * SIZE_BUFFERED) +
        (localVector.z * SIZE_BUFFERED * SIZE_BUFFERED)
      )->NumberValue()
    ) * 10;
    unsigned int c = _getBiomeColor(elevation, moisture);
    float r = (float)((c >> (8 * 2)) & 0xFF) / 0xFF;
    float g = (float)((c >> (8 * 1)) & 0xFF) / 0xFF;
    float b = (float)((c >> (8 * 0)) & 0xFF) / 0xFF;
    for (unsigned int j = 0; j < 3; j++) {
      unsigned int positionBaseIndex = triangleBaseIndex + (j * 3);
      colors->Set(positionBaseIndex + 0, v8::Number::New(isolate, r));
      colors->Set(positionBaseIndex + 1, v8::Number::New(isolate, g));
      colors->Set(positionBaseIndex + 2, v8::Number::New(isolate, b));
    }
  }
  landResult->Set(colorsKey, colors);
  // end add land colors

  // construct result
  v8::Local<v8::Object> result = v8::Object::New(isolate);
  result->Set(landKey, landResult);
  result->Set(waterKey, waterResult);
  result->Set(metadataKey, metadataResult);

  args.GetReturnValue().Set(scope.Escape(result));
}
void ReMarchCubes(const v8::FunctionCallbackInfo<v8::Value>& args) {
  v8::Isolate* isolate = args.GetIsolate();
  v8::EscapableHandleScope scope(isolate);

  auto originKey = v8::String::NewFromUtf8(isolate, "origin");
  auto landKey = v8::String::NewFromUtf8(isolate, "land");
  auto waterKey = v8::String::NewFromUtf8(isolate, "water");
  auto metadataKey = v8::String::NewFromUtf8(isolate, "metadata");
  auto positionsKey = v8::String::NewFromUtf8(isolate, "positions");
  auto colorsKey = v8::String::NewFromUtf8(isolate, "colors");
  auto etherKey = v8::String::NewFromUtf8(isolate, "ether");
  auto moistureEtherKey = v8::String::NewFromUtf8(isolate, "moistureEther");

  // Check the number of arguments passed.
  if (args.Length() < 1) {
    // Throw an Error that is passed back to JavaScript
    isolate->ThrowException(v8::Exception::TypeError(
        v8::String::NewFromUtf8(isolate, "Wrong number of arguments")));
    return;
  }
  // Check the argument types
  if (!args[0]->IsObject()) {
    isolate->ThrowException(v8::Exception::TypeError(
        v8::String::NewFromUtf8(isolate, "Wrong arguments")));
    return;
  }

  v8::Local<v8::Object> opts = args[0]->ToObject();
  v8::Local<v8::Value> origin = opts->Get(originKey);
  v8::Local<v8::Value> landOption = opts->Get(landKey);
  v8::Local<v8::Value> waterOption = opts->Get(waterKey);
  v8::Local<v8::Value> metadataOption = opts->Get(metadataKey);
  if (!(origin->IsArray() && landOption->IsObject() && waterOption->IsObject() && metadataOption->IsObject())) {
    isolate->ThrowException(v8::Exception::TypeError(
        v8::String::NewFromUtf8(isolate, "Invalid options")));
    return;
  }

  v8::Local<v8::Array> originValue = origin.As<v8::Array>();

  v8::Local<v8::Object> landValue = landOption.As<v8::Object>();
  v8::Local<v8::Value> landEther = landValue->Get(etherKey);

  v8::Local<v8::Object> waterValue = waterOption.As<v8::Object>();
  v8::Local<v8::Value> waterEther = waterValue->Get(etherKey);

  v8::Local<v8::Object> metadataValue = metadataOption.As<v8::Object>();
  v8::Local<v8::Value> moistureEther = metadataValue->Get(moistureEtherKey);

  if (!(originValue->Length() == 3 && landEther->IsFloat32Array() && waterEther->IsFloat32Array() && moistureEther->IsFloat32Array())) {
    isolate->ThrowException(v8::Exception::TypeError(
        v8::String::NewFromUtf8(isolate, "Invalid options")));
    return;
  }

  IntVector3 originVector{
    originValue->Get(0)->Int32Value(),
    originValue->Get(1)->Int32Value(),
    originValue->Get(2)->Int32Value()
  };
  v8::Local<v8::Float32Array> landEtherValue = landEther.As<v8::Float32Array>();
  v8::Local<v8::Float32Array> waterEtherValue = waterEther.As<v8::Float32Array>();
  v8::Local<v8::Float32Array> moistureEtherValue = moistureEther.As<v8::Float32Array>();

  v8::Local<v8::Object> landResult = DoReMarchCubes(
    isolate,
    landEtherValue
  ).As<v8::Object>();

  v8::Local<v8::Object> waterResult = DoReMarchCubes(
    isolate,
    waterEtherValue
  ).As<v8::Object>();

  // begin add land colors
  v8::Local<v8::Float32Array> positions = landResult->Get(positionsKey).As<v8::Float32Array>();
  unsigned int numPositions = positions->Length();
  unsigned int numTriangles = numPositions / 3;
  v8::Local<v8::ArrayBuffer> colorsBuffer = v8::ArrayBuffer::New(isolate, numPositions * 4);
  v8::Local<v8::Float32Array> colors = v8::Float32Array::New(colorsBuffer, 0, numPositions);
  for (unsigned int i = 0; i < numTriangles; i++) {
    unsigned int triangleBaseIndex = i * 3 * 3;

    Vector3 pa{
      (float)positions->Get(triangleBaseIndex + 0)->NumberValue(),
      (float)positions->Get(triangleBaseIndex + 1)->NumberValue(),
      (float)positions->Get(triangleBaseIndex + 2)->NumberValue()
    };
    Vector3 pb{
      (float)positions->Get(triangleBaseIndex + 3)->NumberValue(),
      (float)positions->Get(triangleBaseIndex + 4)->NumberValue(),
      (float)positions->Get(triangleBaseIndex + 5)->NumberValue()
    };
    Vector3 pc{
      (float)positions->Get(triangleBaseIndex + 6)->NumberValue(),
      (float)positions->Get(triangleBaseIndex + 7)->NumberValue(),
      (float)positions->Get(triangleBaseIndex + 8)->NumberValue()
    };
    IntVector3 centerVector{
      (int)((pa.x + pb.x + pc.x) / 3),
      (int)((pa.y + pb.y + pc.y) / 3),
      (int)((pa.z + pb.z + pc.z) / 3)
    };
    IntVector3 localVector{
      centerVector.x + (SIZE_BUFFERED / 2),
      centerVector.y + (SIZE_BUFFERED / 2),
      centerVector.z + (SIZE_BUFFERED / 2)
    };
    Vector3 absoluteVector{
      centerVector.x + (originVector.x * SIZE),
      centerVector.y + (originVector.y * SIZE),
      centerVector.z + (originVector.z * SIZE)
    };
    float elevation = std::sqrt(
      absoluteVector.x * absoluteVector.x +
      absoluteVector.y * absoluteVector.y +
      absoluteVector.z * absoluteVector.z
    );
    float moisture = std::abs(
      (float)moistureEtherValue->Get(
        localVector.x +
        (localVector.y * SIZE_BUFFERED) +
        (localVector.z * SIZE_BUFFERED * SIZE_BUFFERED)
      )->NumberValue()
    ) * 10;
    unsigned int c = _getBiomeColor(elevation, moisture);
    float r = (float)((c >> (8 * 2)) & 0xFF) / 0xFF;
    float g = (float)((c >> (8 * 1)) & 0xFF) / 0xFF;
    float b = (float)((c >> (8 * 0)) & 0xFF) / 0xFF;
    for (unsigned int j = 0; j < 3; j++) {
      unsigned int positionBaseIndex = triangleBaseIndex + (j * 3);
      colors->Set(positionBaseIndex + 0, v8::Number::New(isolate, r));
      colors->Set(positionBaseIndex + 1, v8::Number::New(isolate, g));
      colors->Set(positionBaseIndex + 2, v8::Number::New(isolate, b));
    }
  }
  landResult->Set(colorsKey, colors);
  // end add land colors

  v8::Local<v8::Value> metadataResult = metadataOption;

  // construct result
  v8::Local<v8::Object> result = v8::Object::New(isolate);
  result->Set(landKey, landResult);
  result->Set(waterKey, waterResult);
  result->Set(metadataKey, metadataResult);

  args.GetReturnValue().Set(scope.Escape(result));
}
void HolesMarchCubes(const v8::FunctionCallbackInfo<v8::Value>& args) {
  v8::Isolate* isolate = args.GetIsolate();
  v8::EscapableHandleScope scope(isolate);

  auto originKey = v8::String::NewFromUtf8(isolate, "origin");
  auto holesKey = v8::String::NewFromUtf8(isolate, "holes");
  auto landKey = v8::String::NewFromUtf8(isolate, "land");
  auto waterKey = v8::String::NewFromUtf8(isolate, "water");
  auto metadataKey = v8::String::NewFromUtf8(isolate, "metadata");
  auto positionsKey = v8::String::NewFromUtf8(isolate, "positions");
  auto colorsKey = v8::String::NewFromUtf8(isolate, "colors");
  auto etherKey = v8::String::NewFromUtf8(isolate, "ether");
  auto moistureEtherKey = v8::String::NewFromUtf8(isolate, "moistureEther");

  // Check the number of arguments passed.
  if (args.Length() < 1) {
    // Throw an Error that is passed back to JavaScript
    isolate->ThrowException(v8::Exception::TypeError(
        v8::String::NewFromUtf8(isolate, "Wrong number of arguments")));
    return;
  }
  // Check the argument types
  if (!args[0]->IsObject()) {
    isolate->ThrowException(v8::Exception::TypeError(
        v8::String::NewFromUtf8(isolate, "Wrong arguments")));
    return;
  }

  v8::Local<v8::Object> opts = args[0]->ToObject();
  v8::Local<v8::Value> landOption = opts->Get(landKey);
  v8::Local<v8::Value> waterOption = opts->Get(waterKey);
  v8::Local<v8::Value> metadataOption = opts->Get(metadataKey);
  v8::Local<v8::Value> holes = opts->Get(holesKey);
  if (!(landOption->IsObject() && waterOption->IsObject() && metadataOption->IsObject() && holes->IsInt32Array())) {
    isolate->ThrowException(v8::Exception::TypeError(
        v8::String::NewFromUtf8(isolate, "Invalid options")));
    return;
  }

  v8::Local<v8::Array> originValue = opts->Get(originKey).As<v8::Array>();

  v8::Local<v8::Object> landValue = landOption.As<v8::Object>();
  v8::Local<v8::Value> landEther = landValue->Get(etherKey);

  v8::Local<v8::Object> waterValue = waterOption.As<v8::Object>();
  v8::Local<v8::Value> waterEther = waterValue->Get(etherKey);

  v8::Local<v8::Object> metadataValue = metadataOption.As<v8::Object>();
  v8::Local<v8::Value> moistureEther = metadataValue->Get(moistureEtherKey);

  if (!(originValue->Length() == 3 && landEther->IsFloat32Array() && waterEther->IsFloat32Array() && moistureEther->IsFloat32Array())) {
    isolate->ThrowException(v8::Exception::TypeError(
        v8::String::NewFromUtf8(isolate, "Invalid options")));
    return;
  }

  IntVector3 originVector{
    originValue->Get(0)->Int32Value(),
    originValue->Get(1)->Int32Value(),
    originValue->Get(2)->Int32Value()
  };
  v8::Local<v8::Int32Array> holesValue = holes.As<v8::Int32Array>();
  v8::Local<v8::Float32Array> landEtherValue = landEther.As<v8::Float32Array>();
  v8::Local<v8::Float32Array> waterEtherValue = waterEther.As<v8::Float32Array>();
  v8::Local<v8::Float32Array> moistureEtherValue = moistureEther.As<v8::Float32Array>();

  v8::Local<v8::Object> landResult = DoHolesMarchCubes(
    isolate,
    originValue,
    landEtherValue,
    holesValue
  ).As<v8::Object>();

  v8::Local<v8::Object> waterResult = DoHolesMarchCubes(
    isolate,
    originValue,
    waterEtherValue,
    holesValue
  ).As<v8::Object>();

  // begin add land colors
  v8::Local<v8::Float32Array> positions = landResult->Get(positionsKey).As<v8::Float32Array>();
  unsigned int numPositions = positions->Length();
  unsigned int numTriangles = numPositions / 3;
  v8::Local<v8::ArrayBuffer> colorsBuffer = v8::ArrayBuffer::New(isolate, numPositions * 4);
  v8::Local<v8::Float32Array> colors = v8::Float32Array::New(colorsBuffer, 0, numPositions);
  for (unsigned int i = 0; i < numTriangles; i++) {
    unsigned int triangleBaseIndex = i * 3 * 3;

    Vector3 pa{
      (float)positions->Get(triangleBaseIndex + 0)->NumberValue(),
      (float)positions->Get(triangleBaseIndex + 1)->NumberValue(),
      (float)positions->Get(triangleBaseIndex + 2)->NumberValue()
    };
    Vector3 pb{
      (float)positions->Get(triangleBaseIndex + 3)->NumberValue(),
      (float)positions->Get(triangleBaseIndex + 4)->NumberValue(),
      (float)positions->Get(triangleBaseIndex + 5)->NumberValue()
    };
    Vector3 pc{
      (float)positions->Get(triangleBaseIndex + 6)->NumberValue(),
      (float)positions->Get(triangleBaseIndex + 7)->NumberValue(),
      (float)positions->Get(triangleBaseIndex + 8)->NumberValue()
    };
    IntVector3 centerVector{
      (int)((pa.x + pb.x + pc.x) / 3),
      (int)((pa.y + pb.y + pc.y) / 3),
      (int)((pa.z + pb.z + pc.z) / 3)
    };
    IntVector3 localVector{
      centerVector.x + (SIZE_BUFFERED / 2),
      centerVector.y + (SIZE_BUFFERED / 2),
      centerVector.z + (SIZE_BUFFERED / 2)
    };
    Vector3 absoluteVector{
      centerVector.x + (originVector.x * SIZE),
      centerVector.y + (originVector.y * SIZE),
      centerVector.z + (originVector.z * SIZE)
    };
    float elevation = std::sqrt(
      absoluteVector.x * absoluteVector.x +
      absoluteVector.y * absoluteVector.y +
      absoluteVector.z * absoluteVector.z
    );
    float moisture = std::abs(
      (float)moistureEtherValue->Get(
        localVector.x +
        (localVector.y * SIZE_BUFFERED) +
        (localVector.z * SIZE_BUFFERED * SIZE_BUFFERED)
      )->NumberValue()
    ) * 10;
    unsigned int c = _getBiomeColor(elevation, moisture);
    float r = (float)((c >> (8 * 2)) & 0xFF) / 0xFF;
    float g = (float)((c >> (8 * 1)) & 0xFF) / 0xFF;
    float b = (float)((c >> (8 * 0)) & 0xFF) / 0xFF;
    for (unsigned int j = 0; j < 3; j++) {
      unsigned int positionBaseIndex = triangleBaseIndex + (j * 3);
      colors->Set(positionBaseIndex + 0, v8::Number::New(isolate, r));
      colors->Set(positionBaseIndex + 1, v8::Number::New(isolate, g));
      colors->Set(positionBaseIndex + 2, v8::Number::New(isolate, b));
    }
  }
  landResult->Set(colorsKey, colors);
  // end add land colors

  v8::Local<v8::Value> metadataResult = metadataOption;

  // construct result
  v8::Local<v8::Object> result = v8::Object::New(isolate);
  result->Set(landKey, landResult);
  result->Set(waterKey, waterResult);
  result->Set(metadataKey, metadataResult);

  args.GetReturnValue().Set(scope.Escape(result));
}

void Init(v8::Local<v8::Object> exports) {
  NODE_SET_METHOD(exports, "genMarchCubes", GenMarchCubes);
  NODE_SET_METHOD(exports, "reMarchCubes", ReMarchCubes);
  NODE_SET_METHOD(exports, "holesMarchCubes", HolesMarchCubes);
}

NODE_MODULE(addon, Init)
