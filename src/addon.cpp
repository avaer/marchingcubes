#include <node.h>

#include "MarchingCubes.h"
#include "LookUpTable.h"
#include "csg.h"
#include "PerlinNoise/PerlinNoise.hpp"

#define PI 3.14159265358979323846
#define SIZE 50

const siv::PerlinNoise elevationNoise(0);
const float elevationNoiseFrequency = 0.01;
const unsigned int elevationNoiseOctaves = 8;

const siv::PerlinNoise moistureNoise(0);
const float moistureNoiseFrequency = 0.04;
const unsigned int moistureNoiseOctaves = 6;

float heightmaps[SIZE * SIZE][6];

typedef struct {
  real x, y, z;
} Vector3;
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
template <size_t npos>
float _sideGenerator(const Vector3 &vector, const Vector3 &normal, const Vector3 &u, const Vector3 &v, float (&heightmap)[npos]) {
  float length = _length(vector);
  if (length > 0) {
    float angle = _angleTo(vector, normal);
    float angleFactor = 1 - (angle / PI);
    unsigned int uIndex = static_cast<unsigned int>(_sum(_multiply(u, vector))) + (SIZE / 2);
    unsigned int vIndex = static_cast<unsigned int>(_sum(_multiply(v, vector))) + (SIZE / 2);
    unsigned int index = uIndex + (vIndex * SIZE);
    float heightValue = heightmap[index];
    float insideOutsideValue = (length <= heightValue) ? -1 : 1;
    float etherValue = insideOutsideValue * angleFactor;
    return etherValue;
  } else {
    return -1;
  }
};
template <size_t npos, size_t nmap>
float _getValue(float x, float y, float z, float (&heightmaps)[npos][nmap]) {
  float dx = x - (SIZE / 2);
  float dy = y - (SIZE / 2);
  float dz = z - (SIZE / 2);

  float v =
    _sideGenerator( // front
      Vector3{dx, dy, dz},
      Vector3{0, 0, 1},
      Vector3{1, 0, 0},
      Vector3{0, 1, 0},
      heightmaps[0]
    ) +
    _sideGenerator( // top
      Vector3{dx, dy, dz},
      Vector3{0, 1, 0},
      Vector3{1, 0, 0},
      Vector3{0, 0, 1},
      heightmaps[1]
    ) +
    _sideGenerator( // bottom
      Vector3{dx, dy, dz},
      Vector3{0, 1, 0},
      Vector3{1, 0, 0},
      Vector3{0, 0, 1},
      heightmaps[2]
    ) +
    _sideGenerator( // left
      Vector3{dx, dy, dz},
      Vector3{1, 0, 0},
      Vector3{0, 0, 1},
      Vector3{0, -1, 0},
      heightmaps[3]
    ) +
    _sideGenerator( // right
      Vector3{dx, dy, dz},
      Vector3{1, 0, 0},
      Vector3{0, 0, 1},
      Vector3{0, 1, 0},
      heightmaps[4]
    ) +
    _sideGenerator( // back
      Vector3{dx, dy, dz},
      Vector3{0, 0, 1},
      Vector3{1, 0, 0},
      Vector3{0, 1, 0},
      heightmaps[5]
    );
  return v;
}
template <size_t npos>
void _setHeightmap(float (&heightmap)[npos], const Vector2 &uv) {
  for (unsigned int i = 0; i < SIZE; i++) {
    for (unsigned int j = 0; j < SIZE; j++) {
      unsigned int index = i + (j * SIZE);
      float v = 10 +
        elevationNoise.octaveNoise(
          (SIZE * 2) + (((uv.x * SIZE) + i) * elevationNoiseFrequency), // offset to avoid artifacts at the origin
          (SIZE * 2) + (((uv.y * SIZE) + j) * elevationNoiseFrequency),
          elevationNoiseOctaves
        ) * (SIZE - 10);
      heightmap[index] = v;
    }
  }
}

v8::Local<v8::Value> DoMarchCubes(const v8::FunctionCallbackInfo<v8::Value>& args) {
  v8::Isolate* isolate = args.GetIsolate();
  auto positionsKey = v8::String::NewFromUtf8(isolate, "positions");
  auto normalsKey = v8::String::NewFromUtf8(isolate, "normals");

  // Check the number of arguments passed.
  if (args.Length() < 1) {
    // Throw an Error that is passed back to JavaScript
    isolate->ThrowException(v8::Exception::TypeError(
        v8::String::NewFromUtf8(isolate, "Wrong number of arguments")));
    return v8::Null(isolate);
  }
  // Check the argument types
  if (!args[0]->IsObject()) {
    isolate->ThrowException(v8::Exception::TypeError(
        v8::String::NewFromUtf8(isolate, "Wrong arguments")));
    return v8::Null(isolate);
  }

  v8::Local<v8::Object> opts = args[0]->ToObject();
  v8::Local<v8::Value> holes = opts->Get(v8::String::NewFromUtf8(isolate, "holes"));
  if (!(holes->IsInt32Array())) {
    isolate->ThrowException(v8::Exception::TypeError(
        v8::String::NewFromUtf8(isolate, "Invalid options")));
    return v8::Null(isolate);
  }

  v8::Local<v8::Int32Array> holesArray = holes.As<v8::Int32Array>();

  // generate heightmap
  _setHeightmap(heightmaps[0], Vector2{0, 0}); // front
  _setHeightmap(heightmaps[1], Vector2{0, -1}); // top
  _setHeightmap(heightmaps[2], Vector2{0, 1}); // bottom
  _setHeightmap(heightmaps[3], Vector2{-1, 0}); // left
  _setHeightmap(heightmaps[4], Vector2{1, 0}); // right
  _setHeightmap(heightmaps[5], Vector2{2, 0}); // back

  // begin mc
  // Init data
  MarchingCubes mc;
  mc.set_resolution(SIZE, SIZE, SIZE);
  mc.init_all();

  // Fills data structure
  for (unsigned int i = 0; i < SIZE; i++) {
    for (unsigned int j = 0; j < SIZE; j++) {
      for (unsigned int k = 0; k < SIZE; k++) {
        float v = _getValue(i, j, k, heightmaps);
        mc.set_data(v, i, j, k);
      }
    }
  }

  // Adjust for holes
  unsigned int numHoles = holesArray->Length() / 3;
  for (unsigned int h = 0; h < numHoles; h++) {
    unsigned int holeIndexBase = h * 3;
    int x = holesArray->Get(holeIndexBase + 0)->Int32Value();
    int y = holesArray->Get(holeIndexBase + 1)->Int32Value();
    int z = holesArray->Get(holeIndexBase + 2)->Int32Value();

    for (int i = -1; i <= 1; i++) {
      int dx = x + i;

      if (dx >= 0 && dx < SIZE) {
        for (int j = -1; j <= 1; j++) {
          int dy = y + j;

          if (dx >= 0 && dx < SIZE) {
            for (int k = -1; k <= 1; k++) {
              int dz = z + k;

              if (dz >= 0 && dz < SIZE) {
                float distance = std::sqrt((i * i) + (j * j) + (k * k));
                float distanceFactor = distance / std::sqrt(3);
                float valueFactor = 1 - distanceFactor;
                float v = valueFactor * 1.5;

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

  // mc.set_method(true);
  mc.run();
  // end mc

  // construct result
  unsigned int numTrigs = mc.ntrigs();
  unsigned int numVerts = numTrigs * 3;
  unsigned int numPositions = numVerts * 3;
  unsigned int numNormals = numVerts * 3;
  v8::Local<v8::Float32Array> positions = v8::Float32Array::New(v8::ArrayBuffer::New(isolate, numPositions * 4), 0, numPositions);
  v8::Local<v8::Float32Array> normals = v8::Float32Array::New(v8::ArrayBuffer::New(isolate, numNormals * 4), 0, numNormals);

  auto triangles = mc.triangles();
  auto vertices = mc.vertices();
  for (unsigned int i = 0; i < numTrigs; i++) {
    const Triangle &triangle = triangles[i];
    const Vertex &a = vertices[triangle.v1];
    const Vertex &b = vertices[triangle.v2];
    const Vertex &c = vertices[triangle.v3];

    unsigned int baseIndex = i * 3 * 3;
    positions->Set(baseIndex + 0, v8::Number::New(isolate, a.x - (SIZE / 2)));
    positions->Set(baseIndex + 1, v8::Number::New(isolate, a.y - (SIZE / 2)));
    positions->Set(baseIndex + 2, v8::Number::New(isolate, a.z - (SIZE / 2)));
    positions->Set(baseIndex + 3, v8::Number::New(isolate, b.x - (SIZE / 2)));
    positions->Set(baseIndex + 4, v8::Number::New(isolate, b.y - (SIZE / 2)));
    positions->Set(baseIndex + 5, v8::Number::New(isolate, b.z - (SIZE / 2)));
    positions->Set(baseIndex + 6, v8::Number::New(isolate, c.x - (SIZE / 2)));
    positions->Set(baseIndex + 7, v8::Number::New(isolate, c.y - (SIZE / 2)));
    positions->Set(baseIndex + 8, v8::Number::New(isolate, c.z - (SIZE / 2)));

    normals->Set(baseIndex + 0, v8::Number::New(isolate, a.nx));
    normals->Set(baseIndex + 1, v8::Number::New(isolate, a.ny));
    normals->Set(baseIndex + 2, v8::Number::New(isolate, a.nz));
    normals->Set(baseIndex + 3, v8::Number::New(isolate, b.nx));
    normals->Set(baseIndex + 4, v8::Number::New(isolate, b.ny));
    normals->Set(baseIndex + 5, v8::Number::New(isolate, b.nz));
    normals->Set(baseIndex + 6, v8::Number::New(isolate, c.nx));
    normals->Set(baseIndex + 7, v8::Number::New(isolate, c.ny));
    normals->Set(baseIndex + 8, v8::Number::New(isolate, c.nz));
  }

  v8::Local<v8::Object> result = v8::Object::New(isolate);
  result->Set(positionsKey, positions);
  result->Set(normalsKey, normals);

  return result;
}

void MarchCubes(const v8::FunctionCallbackInfo<v8::Value>& args) {
  args.GetReturnValue().Set(DoMarchCubes(args));
}
Biome _getBiome(float elevation, float moisture, float size) {
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
  } else */if (elevation > (size * 0.3)) {
    if (moisture > 0.50) { return Biome::SNOW; }
    else if (moisture > 0.33) { return Biome::TUNDRA; }
    else if (moisture > 0.16) { return Biome::BARE; }
    else { return Biome::SCORCHED; }
  } else if (elevation > (size * 0.25)) {
    if (moisture > 0.66) { return Biome::TAIGA; }
    else if (moisture > 0.33) { return Biome::SHRUBLAND; }
    else { return Biome::TEMPERATE_DESERT; }
  } else if (elevation > (size * 0.1)) {
    if (moisture > 0.83) { return Biome::TEMPERATE_RAIN_FOREST; }
    else if (moisture > 0.50) { return Biome::TEMPERATE_DECIDUOUS_FOREST; }
    else if (moisture > 0.16) { return Biome::GRASSLAND; }
    else { return Biome::TEMPERATE_DESERT; }
  } else {
    if (moisture > 0.66) { return Biome::TROPICAL_RAIN_FOREST; }
    else if (moisture > 0.33) { return Biome::TROPICAL_SEASONAL_FOREST; }
    else if (moisture > 0.16) { return Biome::GRASSLAND; }
    else { return Biome::SUBTROPICAL_DESERT; }
  }
};
unsigned int _getBiomeColor(float elevation, float moisture, float size) {
  auto biome = _getBiome(elevation, moisture, size);
  auto biomeColor = (unsigned int)biome;
  return biomeColor;
}
void MarchCubesPlanet(const v8::FunctionCallbackInfo<v8::Value>& args) {
  v8::Isolate* isolate = args.GetIsolate();
  auto positionsKey = v8::String::NewFromUtf8(isolate, "positions");
  auto colorsKey = v8::String::NewFromUtf8(isolate, "colors");

  v8::Local<v8::Object> marchingCubes = DoMarchCubes(args).As<v8::Object>();
  if (!marchingCubes->IsNull()) {
    v8::Local<v8::Float32Array> positions = marchingCubes->Get(positionsKey).As<v8::Float32Array>();

    unsigned int numPositions = positions->Length();
    unsigned int numTriangles = numPositions / 3;
    v8::Local<v8::Float32Array> colors = v8::Float32Array::New(v8::ArrayBuffer::New(isolate, numPositions * 4), 0, numPositions);
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
      Vector3 center{
        (pa.x + pb.x + pc.x) / 3,
        (pa.y + pb.y + pc.y) / 3,
        (pa.z + pb.z + pc.z) / 3
      };
      float elevation = std::sqrt(center.x * center.x + center.y * center.y + center.z * center.z);
      float moisture = moistureNoise.octaveNoise(
        center.x * moistureNoiseFrequency,
        center.y * moistureNoiseFrequency,
        center.z * moistureNoiseFrequency,
        moistureNoiseOctaves
      );
      unsigned int c = _getBiomeColor(elevation, moisture, 50);
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

    marchingCubes->Set(colorsKey, colors);

    args.GetReturnValue().Set(marchingCubes);
  } else {
    args.GetReturnValue().Set(v8::Null(isolate));
  }
}

void Init(v8::Local<v8::Object> exports) {
  NODE_SET_METHOD(exports, "marchCubes", MarchCubes);
  NODE_SET_METHOD(exports, "marchCubesPlanet", MarchCubesPlanet);
}

NODE_MODULE(addon, Init)
