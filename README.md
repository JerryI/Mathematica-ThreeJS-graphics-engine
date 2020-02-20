3D Graphics drawer for Mathematica based on Three.js
===================
Written in JS parser and drawer allows to export and embed to web pages 3D graphics from Wolfram Mathematica. 
Unlike build-in export functions it recreates pure Mathematica's functions like ``Sphere[]``, ``GraphicsComplex[]``, ``Polygon[]`` and etc.

Some parts of the code which is responsible for rotation, zoom, dragging objects and lighting system were taken from Mathics project.
The home page of Mathics is http://mathics.github.io.

Usage
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

3. Copy and paste it to data.js
```javascript
//\data.js

var JSONThree = [...
``` 

4. Run ``index.html``

``Drag`` - rotate;
``Ctrl+Drag`` - zoom;
``Shift+Drag`` - drag;


Single page export
----------
If you want to share your figure via e-mail, you will be able export it to autonomous .html page. Just pass a figure to a function ``Export2ThreeJS`` located in ``Export2ThreeJS.nb``. 

Contributing
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
- ``GeometricTransformation`` - fully supported
- ``GraphicsComplex`` - supported
- ``Polygon`` - fully supported
- ``Line`` - supported


License
-------

Project is released under the GNU General Public License (GPL).
