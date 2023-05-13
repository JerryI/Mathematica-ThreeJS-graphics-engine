# Graphics3D Parser based on THREE.js
*written with love in Javascript*

***Warning: Early Development Stage!***

![ocean](imgs/demo.gif)

```mathematica
Graphics3D[{Roughness[0], Reflectivity[1], IOR[2], Clearcoat[0.5], Table[{RGBColor[Normalize[i]], Sphere[i]}, {i, RandomReal[{-5,5}, {40,3}]}], SkyAndWater[]}, Background->None, Lighting->None]
```

See disscussion at [mathematica.stackexchange](https://mathematica.stackexchange.com/a/215025/53728).

__This is a core component of [Wolfram JS Frontend](https://github.com/JerryI/wolfram-js-frontend) project__
but one can use it independently as well

__This package depends on [WLJS Interpreter](https://github.com/JerryI/wljs-interpreter)__ (will be downloaded automatically on startup)

## Examples
To run examples gallery, you need to have `nodejs` and `wolframscript` installed
```bash
git clone https://github.com/JerryI/Mathematica-ThreeJS-graphics-engine
cd Mathematica-ThreeJS-graphics-engine
npm i
```
to generate gallery and install WL dependencies
```bash
npm run test
```
and run a dev server on your machine
```bash
npm run watch
```


## Contributing
------------

Please feel encouraged to contribute and expand features.

Issues
------
There a lot a functions which are not implemented such as ``Style[]``, ``Tube[]``, ``Ball[]``, ``Cone[]``, ``BezierCurve[]``...

Currently the minimum necessary set for the functioning of ``SphericalPlot3D``, ``Plot3D`` is already done
- ``Graphics3D`` - supported without styling, themes, custom lighting
- ``List`` - supported
- ``GraphicsGroup`` - supported
- ``RGBColor`` - supported
- ``Opacity`` - supported
- ``Tube`` - renders like arrows
- ``Sphere`` - supported
- ``Center`` - supported
- ``Tetrahedron`` - supported
- ``Cylinder`` - supported
- ``Polyhedron`` - supported
- ``GeometricTransformation`` - fully supported
- ``GraphicsComplex`` - supported
- ``Polygon`` - fully supported
- ``Line`` - supported


- ``THREE.Geometry`` is depricated in a new Three.JS version. Needed to be refactored.

### Extra features
- ``Emissive[]`` - property fro the object to emitt light
- ``IOR[]`` - specify the refractive index
- ``Reflectivity[]`` - reflectivity of the material
- ``SkyAndWater[]`` - apply shader to the scene with animated ocean and sun
- subsurface scattering
- bloom control from the menu

## Development

## Tests
you can easily add new scenes by adding files into `tests/src` dir. 
For the complex scenes use `LoadPage["templates/signlepage_nodom.wsp"]` instead of `LoadPage["templates/signlepage.wsp",{data = json}]`.

## License

Project is released under the GNU General Public License (GPL).
