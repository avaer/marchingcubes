#include <node.h>

#include "MarchingCubes.h"
#include "LookUpTable.h"
#include "csg.h"
#include "PerlinNoise/PerlinNoise.hpp"

typedef struct
{
  real x, y, z;
} Vector;

// grid extension
const float xmin=-1.0f, xmax=1.0f,  ymin=-1.0f, ymax=1.0f,  zmin=-1.0f, zmax=1.0f;
const siv::PerlinNoise moistureNoise(0);
const float moistureNoiseFrequency = 0.04;
const unsigned int moistureNoiseOctaves = 6;

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
  v8::Local<v8::Value> width = opts->Get(v8::String::NewFromUtf8(isolate, "width"));
  v8::Local<v8::Value> height = opts->Get(v8::String::NewFromUtf8(isolate, "height"));
  v8::Local<v8::Value> depth = opts->Get(v8::String::NewFromUtf8(isolate, "depth"));
  v8::Local<v8::Value> data = opts->Get(v8::String::NewFromUtf8(isolate, "data"));
  if (!(width->IsNumber() && height->IsNumber() && depth->IsNumber() && data->IsFloat32Array())) {
    isolate->ThrowException(v8::Exception::TypeError(
        v8::String::NewFromUtf8(isolate, "Invalid options")));
    return v8::Null(isolate);
  }

  v8::Local<v8::Float32Array> dataArray = data.As<v8::Float32Array>();

  // begin mc
  // grid size control
  int size_x = width->ToInt32()->Value(), size_y = height->ToInt32()->Value(), size_z = depth->ToInt32()->Value();
  int half_size_x = size_x / 2, half_size_y = size_y / 2, half_size_z = size_z / 2;

  // Init data
  MarchingCubes mc;
  mc.set_resolution( size_x, size_y, size_z ) ;
  mc.init_all() ;

  // Fills data structure
  int i,j,k ;
  float w ;
  unsigned int index;
  for( i = 0 ; i < size_x ; i++ )
  {
    for( j = 0 ; j < size_y ; j++ )
    {
      for( k = 0 ; k < size_z ; k++ )
      {
        index = i + (j * size_x) + (k * size_x * size_y);
        w = (float)dataArray->Get(index)->NumberValue();
        mc.set_data( w, i, j, k ) ;
      }
    }
  }

  // mc.set_method(true);
  mc.run();
  // end mc

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
    positions->Set(baseIndex + 0, v8::Number::New(isolate, a.x - half_size_x));
    positions->Set(baseIndex + 1, v8::Number::New(isolate, a.y - half_size_y));
    positions->Set(baseIndex + 2, v8::Number::New(isolate, a.z - half_size_z));
    positions->Set(baseIndex + 3, v8::Number::New(isolate, b.x - half_size_x));
    positions->Set(baseIndex + 4, v8::Number::New(isolate, b.y - half_size_y));
    positions->Set(baseIndex + 5, v8::Number::New(isolate, b.z - half_size_z));
    positions->Set(baseIndex + 6, v8::Number::New(isolate, c.x - half_size_x));
    positions->Set(baseIndex + 7, v8::Number::New(isolate, c.y - half_size_y));
    positions->Set(baseIndex + 8, v8::Number::New(isolate, c.z - half_size_z));

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
unsigned int _getBiomeColor(float elevation, float moisture) { // XXX
  return 0x0;
}
void MarchCubesPlanet(const v8::FunctionCallbackInfo<v8::Value>& args) {
  v8::Isolate* isolate = args.GetIsolate();
  auto positionsKey = v8::String::NewFromUtf8(isolate, "positions");
  auto normalsKey = v8::String::NewFromUtf8(isolate, "normals");
  auto colorsKey = v8::String::NewFromUtf8(isolate, "colors");

  v8::Local<v8::Object> marchingCubes = DoMarchCubes(args).As<v8::Object>();
  v8::Local<v8::Float32Array> positions = marchingCubes->Get(positionsKey).As<v8::Float32Array>();
  v8::Local<v8::Float32Array> normals = marchingCubes->Get(normalsKey).As<v8::Float32Array>();

  unsigned int numPositions = positions->Length() / 3;
  unsigned int numTriangles = numPositions / 3;
  v8::Local<v8::Float32Array> colors = v8::Float32Array::New(v8::ArrayBuffer::New(isolate, numPositions * 4), 0, numPositions);
  for (unsigned int i = 0; i < numTriangles; i++) {
    unsigned int triangleBaseIndex = i * 3 * 3;

    Vector pa{
      positions->Get(triangleBaseIndex + 0)->NumberValue(),
      positions->Get(triangleBaseIndex + 1)->NumberValue(),
      positions->Get(triangleBaseIndex + 2)->NumberValue()
    };
    Vector pb{
      positions->Get(triangleBaseIndex + 3)->NumberValue(),
      positions->Get(triangleBaseIndex + 4)->NumberValue(),
      positions->Get(triangleBaseIndex + 5)->NumberValue()
    };
    Vector pc{
      positions->Get(triangleBaseIndex + 6)->NumberValue(),
      positions->Get(triangleBaseIndex + 7)->NumberValue(),
      positions->Get(triangleBaseIndex + 8)->NumberValue()
    };
    Vector center{
      (pa.x + pb.x + pb.x) / 3,
      (pa.y + pb.y + pb.y) / 3,
      (pa.z + pb.z + pb.z) / 3
    };
    float elevation = std::sqrt(center.x * center.x + center.y * center.y + center.z * center.z);
    float moisture = moistureNoise.octaveNoise(center.x * moistureNoiseFrequency, center.y * moistureNoiseFrequency, center.z, moistureNoiseOctaves);
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

  marchingCubes->Set(colorsKey, colors);

  args.GetReturnValue().Set(marchingCubes);
}

void Init(v8::Local<v8::Object> exports) {
  NODE_SET_METHOD(exports, "marchCubes", MarchCubes);
  NODE_SET_METHOD(exports, "marchCubesPlanet", MarchCubesPlanet);
}

NODE_MODULE(addon, Init)
