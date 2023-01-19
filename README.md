3D Graphics drawer for Mathematica based on Three.js
===================
Written in JS parser and drawer allow to export or embed to web pages 3D graphics from Wolfram notebook. 
Unlike other build-in export functions it recreates pure Mathematica's functions like ``Sphere[]``, ``GraphicsComplex[]``, ``Polygon[]`` and etc. See disscussion at [mathematica.stackexchange](https://mathematica.stackexchange.com/a/215025/53728).

Some parts of the code which is responsible for rotation, zoom, dragging objects and lighting system were taken from Mathics project.
The home page of Mathics is http://mathics.github.io.

## Live example
----------
[See](https://jerryi.github.io/Mathematica-ThreeJS-graphics-engine/)

----------
## This project is a part of a bigger one [Wolfram Engine JS Frontend](https://github.com/JerryI/wolfram-js-frontend)
----------

## Usage
----------
1. Plot some graphics (used a low-poly mode for smaller code, see ``Example.nb``)

```Mathematica
Graphics3D[{
  SphericalPlot3D[
    2 SphericalHarmonicY[2, 0, t, p], {t, 0, Pi}, {p, 0, 2 Pi}, 
    PerformanceGoal -> "Speed"][[1]],
  Opacity[0.6], 
  Tetrahedron[{{1, 1, 1}, {-1, -1, 1}, {1, -1, -1}, {-1, 1, -1}}]
  }]
```
2. Export as a JSON string
```Mathematica
ExportString[%//N, "ExpressionJSON"]
```
```Mathematica
[
	"Graphics3D",
	[
		"List",
		[
			"GraphicsComplex",
			[
				"List",
				["List",
					0.0,
					0.0,
					1.2615662610100797
				]
				,
				["List",
					0.0,
					0.0,
					1.2615662610100797
				]
				,...
```

3. Run `index.html`

4. Paste the JSON code as a plain text into the textarea

5. Use

``Drag`` - to rotate;
``Ctrl+Drag`` - to zoom;
``Shift+Drag`` - to drag;

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
- ``RGBColor`` - supported
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

Possible to implement as well
------
There are a lot of features of Three.JS, which can bring extra control over the graphics repesentating and are not a part of FrontEnd of Mathematica. 

- ``FlatShading``, ``PhongShading``
- ``PathTracing``
- extended material properties

## Development

download the type file for three js.

```bash
npm i
```

to start dev server
```bash
npm run dev
```

and open you browser at `http://127.0.0.1:8090/dev.html`.
The index files relies on CDN (and release folder) and will not show any changes. 

the file watcher will automatically rebuild the all stuff on change in `src` dir.

## License
-------

Project is released under the GNU General Public License (GPL).
