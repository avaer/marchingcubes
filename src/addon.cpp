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

v8::Local<v8::Value> DoMarchCubes(
  v8::Isolate* isolate,
  unsigned int seedNumber,
  v8::Local<v8::Array> originValue,
  v8::Local<v8::Int32Array> holesValue,
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

  // generate heightmap
  IntVector3 originVector{
    originValue->Get(0)->Int32Value(),
    originValue->Get(1)->Int32Value(),
    originValue->Get(2)->Int32Value()
  };

  noise.reseed(seedNumber);

  // begin mc
  // Init data
  MarchingCubes mc(SIZE_BUFFERED, SIZE_BUFFERED, SIZE_BUFFERED);
  // mc.set_resolution();
  // mc.init_all();

  // Fills data structure
  for (unsigned int i = 0; i < SIZE_BUFFERED; i++) {
    for (unsigned int j = 0; j < SIZE_BUFFERED; j++) {
      for (unsigned int k = 0; k < SIZE_BUFFERED; k++) {
        float v = _getValue(
          Vector3{
            i,
            j,
            k
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

  // Adjust for holes
  unsigned int numHoles = holesValue->Length() / 3;
  for (unsigned int h = 0; h < numHoles; h++) {
    unsigned int holeIndexBase = h * 3;
    int x = holesValue->Get(holeIndexBase + 0)->Int32Value();
    int y = holesValue->Get(holeIndexBase + 1)->Int32Value();
    int z = holesValue->Get(holeIndexBase + 2)->Int32Value();

    for (int i = -2; i <= 2; i++) {
      int dx = BUFFER + x + i - (originVector.x * SIZE);

      if (dx >= 0 && dx < SIZE_BUFFERED) {
        for (int j = -2; j <= 2; j++) {
          int dy = BUFFER + y + j - (originVector.y * SIZE);

          if (dy >= 0 && dy < SIZE_BUFFERED) {
            for (int k = -2; k <= 2; k++) {
              int dz = BUFFER + z + k - (originVector.z * SIZE);

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

  // mc.set_method(true);
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

  v8::Local<v8::Object> result = v8::Object::New(isolate);
  result->Set(positionsKey, v8::Float32Array::New(positionsBuffer, 0, index * 3 * 3));
  result->Set(normalsKey, v8::Float32Array::New(normalsBuffer, 0, index * 3 * 3));

  return scope.Escape(result);
}

void MarchCubes(const v8::FunctionCallbackInfo<v8::Value>& args) {
  v8::Isolate* isolate = args.GetIsolate();
  v8::EscapableHandleScope scope(isolate);

  auto seedKey = v8::String::NewFromUtf8(isolate, "seed");
  auto originKey = v8::String::NewFromUtf8(isolate, "origin");
  auto holesKey = v8::String::NewFromUtf8(isolate, "holes");

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
  v8::Local<v8::Value> holes = opts->Get(holesKey);
  if (!(seed->IsNumber() && origin->IsArray() && holes->IsInt32Array())) {
    isolate->ThrowException(v8::Exception::TypeError(
        v8::String::NewFromUtf8(isolate, "Invalid options")));
    return;
  }

  v8::Local<v8::Array> originValue = origin.As<v8::Array>();
  if (originValue->Length() != 3) {
    isolate->ThrowException(v8::Exception::TypeError(
        v8::String::NewFromUtf8(isolate, "Invalid origin array")));
    return;
  }

  v8::Local<v8::Number> seedValue = seed.As<v8::Number>();
  v8::Local<v8::Int32Array> holesValue = holes.As<v8::Int32Array>();

  v8::Local<v8::Value> result = DoMarchCubes(
    isolate,
    seedValue->Uint32Value(),
    originValue,
    holesValue,
    elevationNoiseFrequency,
    elevationNoiseOctaves,
    0,
    1,
    2.25, // lengthPow
    800 // etherFactor
  );

  args.GetReturnValue().Set(scope.Escape(result));
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
void MarchCubesPlanet(const v8::FunctionCallbackInfo<v8::Value>& args) {
  v8::Isolate* isolate = args.GetIsolate();
  v8::EscapableHandleScope scope(isolate);

  auto seedKey = v8::String::NewFromUtf8(isolate, "seed");
  auto originKey = v8::String::NewFromUtf8(isolate, "origin");
  auto holesKey = v8::String::NewFromUtf8(isolate, "holes");
  auto landKey = v8::String::NewFromUtf8(isolate, "land");
  auto waterKey = v8::String::NewFromUtf8(isolate, "water");
  auto positionsKey = v8::String::NewFromUtf8(isolate, "positions");
  auto colorsKey = v8::String::NewFromUtf8(isolate, "colors");

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
  v8::Local<v8::Value> holes = opts->Get(holesKey);
  if (!(seed->IsNumber() && origin->IsArray() && holes->IsInt32Array())) {
    isolate->ThrowException(v8::Exception::TypeError(
        v8::String::NewFromUtf8(isolate, "Invalid options")));
    return;
  }

  v8::Local<v8::Array> originValue = origin.As<v8::Array>();
  if (originValue->Length() != 3) {
    isolate->ThrowException(v8::Exception::TypeError(
        v8::String::NewFromUtf8(isolate, "Invalid origin array")));
    return;
  }

  v8::Local<v8::Number> seedValue = seed.As<v8::Number>();
  v8::Local<v8::Int32Array> holesValue = holes.As<v8::Int32Array>();

  v8::Local<v8::Object> land = DoMarchCubes(
    isolate,
    seedValue->Uint32Value() + 0,
    originValue,
    holesValue,
    elevationNoiseFrequency,
    elevationNoiseOctaves,
    0,
    1,
    2.25, // lengthPow
    2000 // etherFactor
  ).As<v8::Object>();

  std::default_random_engine generator(seedValue->Uint32Value() + 1);
  std::uniform_real_distribution<float> distribution(0.0f, 1.0f);
  const float moistureFactor = 0.2 + 0.8 * distribution(generator);
  v8::Local<v8::Object> water = DoMarchCubes(
    isolate,
    seedValue->Uint32Value() + 1,
    originValue,
    holesValue,
    moistureNoiseFrequency,
    moistureNoiseOctaves,
    4000 * moistureFactor,
    0.2 * moistureFactor,
    2.25, // lengthPow
    2000 // etherFactor
  ).As<v8::Object>();
  if (!land->IsNull() && !water->IsNull()) {
    v8::Local<v8::Object> opts = args[0]->ToObject();
    v8::Local<v8::Array> originValue = opts->Get(originKey).As<v8::Array>();
    v8::Local<v8::Number> seed = opts->Get(seedKey).As<v8::Number>();

    // begin add land colors
    noise.reseed(seedValue->Uint32Value() + 1);

    v8::Local<v8::Float32Array> positions = land->Get(positionsKey).As<v8::Float32Array>();
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
      Vector3 center{
        (pa.x + pb.x + pc.x) / 3,
        (pa.y + pb.y + pc.y) / 3,
        (pa.z + pb.z + pc.z) / 3
      };
      IntVector3 originVector{
        originValue->Get(0)->Int32Value(),
        originValue->Get(1)->Int32Value(),
        originValue->Get(2)->Int32Value()
      };
      Vector3 absoluteVector{
        center.x + (originVector.x * SIZE),
        center.y + (originVector.y * SIZE),
        center.z + (originVector.z * SIZE)
      };
      float elevation = std::sqrt(absoluteVector.x * absoluteVector.x + absoluteVector.y * absoluteVector.y + absoluteVector.z * absoluteVector.z);
      float moisture = std::abs(noise.octaveNoise(
        center.x * moistureNoiseFrequency,
        center.y * moistureNoiseFrequency,
        center.z * moistureNoiseFrequency,
        moistureNoiseOctaves
      )) * 10;
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
    land->Set(colorsKey, colors);
    // end add land colors

    // construct result
    v8::Local<v8::Object> result = v8::Object::New(isolate);
    result->Set(landKey, land);
    result->Set(waterKey, water);

    args.GetReturnValue().Set(scope.Escape(result));
  } else {
    args.GetReturnValue().Set(scope.Escape(v8::Null(isolate)));
  }
}

void Init(v8::Local<v8::Object> exports) {
  NODE_SET_METHOD(exports, "marchCubes", MarchCubes);
  NODE_SET_METHOD(exports, "marchCubesPlanet", MarchCubesPlanet);
}

NODE_MODULE(addon, Init)
