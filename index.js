const marchingcubes = require('./build/Release/addon');

const result = marchingcubes.march({
  width: 50,
  height: 50,
  depth: 50,
});
console.log('got marching cubes', result.positions.length);
