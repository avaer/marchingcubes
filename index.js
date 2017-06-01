const marchingcubes = require('./build/Debug/addon');

module.exports = marchingcubes;

if (require.main === module) {
  const ellipseFn = (x, y, z) => Math.sqrt(Math.pow(x * 2, 2) + Math.pow(y, 2) + Math.pow(z, 2)) - 0.5;

  const size = 50;
  const width = size;
  const height = size;
  const depth = size;
  const _getCoordIndex = (x, y, z) => x + (y * width) + (z * width * height);
  const data = (() => {
    const result = new Float32Array(width * height * depth);

    for (let x = 0; x < width; x++) {
      for (let y = 0; y < height; y++) {
        for (let z = 0; z < depth; z++) {
          const index = _getCoordIndex(x, y, z);
          result[index] = ellipseFn(
            ((x / width) - 0.5) * 2,
            ((y / height) - 0.5) * 2,
            ((z / depth) - 0.5) * 2
          );
        }
      }
    }

    return result;
  })();

  const result = marchingcubes.march({
    width: size,
    height: size,
    depth: size,
    data: data,
  });
  console.log('got marching cubes', result.positions.length);
};
