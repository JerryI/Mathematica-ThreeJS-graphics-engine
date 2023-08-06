# Graphics3D Parser based on THREE.js
*written with love in Javascript*

***Warning: Early Development Stage!***

![](imgs/screenshot(17).png)

```mathematica
Graphics3D[{
    {Emissive[Red], Sphere[{0,0,2}]}, 
    {White, Sphere[]}
}, Lighting->None, RTX->True]
```

😋 Realtime pathtracing is now supported! See [dev.blog](https://jerryi.github.io/wljs-docs/blog/intro-transform-3d)

---

See disscussion at [mathematica.stackexchange](https://mathematica.stackexchange.com/a/215025/53728).

__This is a core component of [Wolfram JS Frontend](https://github.com/JerryI/wolfram-js-frontend) project__,
but you can __use it independently as well__ - [here is how](https://jerryi.github.io/wljs-docs/docs/interpreter/intro).

__Live demo @ [WLJS Interpreter](https://jerryi.github.io/wljs-interpreter/?example=boat.txt) sandbox__ and [here](https://jerryi.github.io/wljs-docs/docs/interpreter/intro) is DOCS for it

## Examples
Most Mathematica's functions for 3D plotting expands into a bunch of `Graphics3D` primitives

```mathematica
ContourPlot3D[x^3 + y^2 - z^2 == 0, {x, -2, 2}, {y, -2, 2}, {z, -2, 2}]
```

![](imgs/screenshot(20).png)

However, styling and labeling is not implemented

```mathematica
VectorPlot3D[{x, y, z}, {x, -1, 1}, {y, -1, 1}, {z, -1, 1}, VectorColorFunction -> Function[{x, y, z, vx, vy, vz, n}, ColorData["ThermometerColors"][x]]][[1]];
%/. {RGBColor[r_,g_,b_] :> Sequence[RGBColor[r/50,g/50,b/50], Emissive[RGBColor[r,g,b], 5]],};

Graphics3D[{%, Roughness[0], Sphere[{0,0,0}, 0.9]}, Lighting->None, RTX->True]
```

![](imgs/screenshot(19).png)

Custom lighting, mesh materials, shadows propeties are provided

![](imgs/screenshot(5).png)

![](imgs/screenshot(8).png)

## Docs?
Will be soon!

## Contributing

Please feel encouraged to contribute and expand features.

![](imgs/screenshot(16).png)

Issues
------
There a lot a functions which are not implemented such as ``Style[]``, ``Tube[]``, ``Ball[]``, ``Cone[]``, ``BezierCurve[]``...

Currently the minimum necessary set for the functioning of ``SphericalPlot3D``, ``Plot3D`` is already done
- ``Graphics3D`` - supported without styling, themes
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

### Extra features
- ``Emissive[]`` - property fro the object to emitt light
- ``Roughness[]`` - roughness of the material
- many other stuff, please, check `src/autocomplete.js`

## Development

```mathematica
wolframscript -f buildtests.wls
npm run watch
```

you can easily add new scenes by adding files into `tests/src` dir. 

> For the complex scenes use `LoadPage["templates/signlepage_nodom.wsp"]` instead of `LoadPage["templates/signlepage.wsp",{data = json}]`.

## License

Project is released under the GNU General Public License (GPL).
