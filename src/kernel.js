const { MathUtils } = require('three');
{

  let g3d = {};
  g3d.name = "WebObjects/Graphics3D";
  interpretate.contextExpand(g3d); 

  ["AlignmentPoint", "AspectRatio", "AutomaticImageSize", "Axes", 
  "AxesEdge", "AxesLabel", "AxesOrigin", "AxesStyle", "Background", 
  "BaselinePosition", "BaseStyle", "Boxed", "BoxRatios", "BoxStyle", 
  "ClipPlanes", "ClipPlanesStyle", "ColorOutput", "ContentSelectable", 
  "ControllerLinking", "ControllerMethod", "ControllerPath", 
  "CoordinatesToolOptions", "DisplayFunction", "Epilog", "FaceGrids", 
  "FaceGridsStyle", "FormatType", "ImageMargins", "ImagePadding", 
  "ImageSize", "ImageSizeRaw", "LabelStyle", "Lighting", "Method", 
  "PlotLabel", "PlotRange", "PlotRangePadding", "PlotRegion", 
  "PreserveImageOptions", "Prolog", "RotationAction", 
  "SphericalRegion", "Ticks", "TicksStyle", "TouchscreenAutoZoom", 
  "ViewAngle", "ViewCenter", "ViewMatrix", "ViewPoint", 
  "ViewProjection", "ViewRange", "ViewVector", "ViewVertical", "Controls", "PointerLockControls"].map((e)=>{
    g3d[e] = () => e;
  });



  /**
  * @type {import('three')}
  */
  let THREE;

  function computeGroupCenter(group) {
    var center = new THREE.Vector3();
    var children = group.children;
    var count = children.length;
    for (var i = 0; i < count; i++) {
      center.add(children[i].position);
    }
    center.divideScalar(count);
    return center;
  }

  g3d.Style = core.List

  /**
   * @description https://threejs.org/docs/#api/en/materials/LineDashedMaterial
   */
  g3d.Dashing = (args, env) => {
    console.log("Dashing not implemented");
  }

  g3d.Annotation = core.List

  g3d.GraphicsGroup = async (args, env) => {
    const group = new THREE.Group();
    let copy = Object.assign({}, env);

    copy.mesh = group;

    for (const a of args) {
      await interpretate(a, copy);
    }

    env.mesh.add(group);
  };

  g3d.Metalness = (args, env) => {
    env.metalness = interpretate(args[0], env);
  }

  g3d.Emissive = async (args, env) => {
    const copy = {...env};
    await interpretate(args[0], copy);
    env.emissive = copy.color;
  }

  g3d.RGBColor = async (args, env) => {
    if (args.length !== 3 && args.length !== 1) {
      console.log("RGB format not implemented", args);
      console.error("RGB values should be triple!");
      return;
    }

    let a = [...args];

    if (args.length === 1) {
      a = await interpretate(args[0], env); // return [r, g, b] , 0<=r, g, b<=1
    }

    const r = await interpretate(a[0], env);
    const g = await interpretate(a[1], env);
    const b = await interpretate(a[2], env);

    env.color = new THREE.Color(r, g, b);
  };

  g3d.Roughness = (args, env) => {
    const o = interpretate(args[0], env);
    if (typeof o !== "number") console.error("Opacity must have number value!");
    console.log(o);
    env.roughness = o;  
  }

  g3d.Opacity = (args, env) => {
    var o = interpretate(args[0], env);
    if (typeof o !== "number") console.error("Opacity must have number value!");
    console.log(o);
    env.opacity = o;
  };

  g3d.ImageScaled = (args, env) => { };

  g3d.Thickness = (args, env) => { env.thickness = interpretate(args[0], env)};

  g3d.Arrowheads = async (args, env) => {
    if (args.length == 1) {
      env.arrowRadius = await interpretate(args[0], env);
    } else {
      env.arrowHeight = await interpretate(args[1], env);
      env.arrowRadius = await interpretate(args[0], env);
    }
  };

  g3d.TubeArrow = async (args, env) => {
    console.log('Context test');
    console.log(this);

    let radius = 1;
    if (args.length > 1) radius = await interpretate(args[1], env);
    /**
     * @type {THREE.Vector3}}
     */
    const coordinates = await interpretate(args[0], env);

    /**
     * @type {THREE.MeshPhysicalMaterial}}
     */  
    const material = new THREE.MeshPhysicalMaterial({
      color: env.color,
      transparent: false,
      roughness: env.roughness,
      opacity: env.opacity,
      metalness: env.metalness,
      emissive: env.emissive,
      reflectivity: env.reflectivity,
      clearcoat: env.clearcoat,
      ior: env.ior
    });

    //points 1, 2
    const p1 = new THREE.Vector3(...coordinates[0]);
    const p2 = new THREE.Vector3(...coordinates[1]);
    //direction
    const dp = p2.clone().addScaledVector(p1, -1);

    const geometry = new THREE.CylinderGeometry(radius, radius, dp.length(), 20, 1);

    //calculate the center (might be done better, i hope BoundingBox doest not envolve heavy computations)
    geometry.computeBoundingBox();
    let position = geometry.boundingBox;

    let center = position.max.addScaledVector(position.min, -1);

    //default geometry
    const cylinder = new THREE.Mesh(geometry, material);

    //cone
    const conegeometry = new THREE.ConeBufferGeometry(env.arrowRadius, env.arrowHeight, 32 );
    const cone = new THREE.Mesh(conegeometry, material);
    cone.position.y = dp.length()/2 + env.arrowHeight/2;

    let group = new THREE.Group();
    group.add(cylinder, cone);

    //the default axis of a Three.js cylinder is [010], then we rotate it to dp vector.
    //using https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    const v = new THREE.Vector3(0, 1, 0).cross(dp.normalize());
    const theta = Math.asin(v.length() / dp.length());
    const sc = Math.sin(theta);
    const mcs = 1.0 - Math.cos(theta);

    //Did not find how to write it using vectors
    const matrix = new THREE.Matrix4().set(
      1 - mcs * (v.y * v.y + v.z * v.z), mcs * v.x * v.y - sc * v.z,/*   */ sc * v.y + mcs * v.x * v.z,/*   */ 0,//
      mcs * v.x * v.y + sc * v.z,/*   */ 1 - mcs * (v.x * v.x + v.z * v.z), -(sc * v.x) + mcs * v.y * v.z,/**/ 0,//
      -(sc * v.y) + mcs * v.x * v.z,/**/ sc * v.x + mcs * v.y * v.z,/*   */ 1 - mcs * (v.x * v.x + v.y * v.y), 0,//
      0,/*                            */0,/*                            */ 0,/**                           */ 1
    );

    //middle target point
    const middle = p1.divideScalar(2.0).addScaledVector(p2, 0.5);

    //shift to the center and rotate
    //group.position = center;
    group.applyMatrix4(matrix);

    //translate its center to the middle target point
    group.position.addScaledVector(middle, -1);

    env.mesh.add(group);

    geometry.dispose();
    conegeometry.dispose();
    material.dispose();
  };

  g3d.Arrow = async (args, env) => {
    let arr = await interpretate(args[0], env);
    if (arr.length === 1) arr = arr[0];
    if (arr.length !== 2) {
      console.error("Tube must have 2 vectors!");
      return;
    }

    const points = [
      new THREE.Vector4(...arr[0], 1),
      new THREE.Vector4(...arr[1], 1),
    ];

    points.forEach((p) => {
      p = p.applyMatrix4(env.matrix);
    });

    const origin = points[0].clone();
    const dir = points[1].add(points[0].negate());

    const arrowHelper = new THREE.ArrowHelper(
      dir.normalize(),
      origin,
      dir.length(),
      env.color,
    );
    env.mesh.add(arrowHelper);
    arrowHelper.line.material.linewidth = env.thickness;
  };

  g3d.Sphere = async (args, env) => {
    var radius = 1;
    if (args.length > 1) radius = await interpretate(args[1], env);

    const material = new THREE.MeshPhysicalMaterial({
      color: env.color,
      roughness: env.roughness,
      opacity: env.opacity,
      metalness: env.metalness,
      emissive: env.emissive,
      reflectivity: env.reflectivity,
      clearcoat: env.clearcoat,
      ior: env.ior
    });

    function addSphere(cr) {
      const origin = new THREE.Vector4(...cr, 1);
      const geometry = new THREE.SphereGeometry(radius, 20, 20);
      const sphere = new THREE.Mesh(geometry, material);

      sphere.position.set(origin.x, origin.y, origin.z);

      env.mesh.add(sphere);
      geometry.dispose();
      return sphere;
    }

    let list = await interpretate(args[0], env);

    if (list.length === 1) list = list[0];
    if (list.length === 1) list = list[0];

    if (list.length === 3) {
      env.local.object = addSphere(list);
    } else if (list.length > 3) {

      env.local.multiple = true;
      env.local.object = [];

      list.forEach((el) => {
        env.local.object.push(addSphere(el));
      });
    } else {
      console.log(list);
      console.error("List of coords. for sphere object is less 1");
      return;
    }

    material.dispose();
  };

  g3d.Sphere.update = async (args, env) => {
    console.log('Sphere: updating the data!');
    console.log(args);
    console.log(env);

    if (env.local.multiple) {
      const data = await interpretate(args[0], env);
      let i = 0;
      data.forEach((c)=>{
        env.local.object[i].position.set(...c);
        ++i;
      });

      return;
    }

    const c = await interpretate(args[0], env);
    env.local.object.position.set(...c);

  }

  g3d.Sphere.destroy = (args, env) => {
    console.log('Sphere: destroy');
    console.log(args);
    console.log(env);

  }

  g3d.Sphere.virtual = true

  g3d.Sky = (args, env) => {
    const sky = new Sky();
  	sky.scale.setScalar( 10000 );
  	env.mesh.add( sky );
    env.sky = sky;
    env.sun = new THREE.Vector3();

  	const skyUniforms = sky.material.uniforms;

  	skyUniforms[ 'turbidity' ].value = 10;
  	skyUniforms[ 'rayleigh' ].value = 2;
  	skyUniforms[ 'mieCoefficient' ].value = 0.005;
  	skyUniforms[ 'mieDirectionalG' ].value = 0.8;
  }

  g3d.Water = (args, env) => {
    const waterGeometry = new THREE.PlaneGeometry( 10000, 10000 );

  	const water = new Water(
  		waterGeometry,
  		{
  			textureWidth: 512,
  			textureHeight: 512,
  			waterNormals: new THREE.TextureLoader().load( 'textures/waternormals.jpg', function ( texture ) {

          texture.wrapS = texture.wrapT = THREE.RepeatWrapping;
  			} ),

        sunDirection: new THREE.Vector3(),
  			sunColor: 0xffffff,
  			waterColor: 0x001e0f,
  			distortionScale: 3.7,
  			fog: true
  		}
  		);

  		water.rotation.x = - Math.PI / 2;

  		env.mesh.add( water );
      env.water = water;
  }

  g3d.Cuboid = async (args, env) => {
    //if (params.hasOwnProperty('geometry')) {
    //	var points = [new THREE.Vector4(...interpretate(func.args[0]), 1),
    //				new THREE.Vector4(...interpretate(func.args[1]), 1)];
    //}
    /**
     * @type {THREE.Vector4}
     */
    var diff;
    /**
     * @type {THREE.Vector4}
     */
    var origin;
    var p;

    if (args.length === 2) {
      var points = [
        new THREE.Vector4(...(await interpretate(args[0], env)), 1),
        new THREE.Vector4(...(await interpretate(args[1], env)), 1),
      ];

      origin = points[0]
        .clone()
        .add(points[1])
        .divideScalar(2);
      diff = points[0].clone().add(points[1].clone().negate());
    } else if (args.length === 1) {
      p = await interpretate(args[0], env);
      origin = new THREE.Vector4(...p, 1);
      diff = new THREE.Vector4(1, 1, 1, 1);

      //shift it
      origin.add(diff.clone().divideScalar(2));
    } else {
      console.error("Expected 2 or 1 arguments");
      return;
    }

    const geometry = new THREE.BoxGeometry(diff.x, diff.y, diff.z);
    const material = new THREE.MeshPhysicalMaterial({
      color: env.color,
      transparent: true,
      opacity: env.opacity,
      roughness: env.roughness,
      depthWrite: true,
      metalness: env.metalness,
      emissive: env.emissive,
      reflectivity: env.reflectivity,
      clearcoat: env.clearcoat,
      ior: env.ior
    });

    //material.side = THREE.DoubleSide;

    const cube = new THREE.Mesh(geometry, material);

    //var tr = new THREE.Matrix4();
    //	tr.makeTranslation(origin.x,origin.y,origin.z);

    //cube.applyMatrix(params.matrix.clone().multiply(tr));

    cube.position.set(origin.x, origin.y, origin.z);

    env.mesh.add(cube);

    geometry.dispose();
    material.dispose();
  };

  g3d.Center = (args, env) => {
    return "Center";
  };

  g3d.Cylinder = async (args, env) => {
    let radius = 1;
    if (args.length > 1) radius = await interpretate(args[1], env);
    /**
     * @type {THREE.Vector3}}
     */
    let coordinates = await interpretate(args[0], env);
    if (coordinates.length === 1) {
      coordinates = coordinates[0];
    }

    coordinates[0] = new THREE.Vector3(...coordinates[0]);
    coordinates[1] = new THREE.Vector3(...coordinates[1]);

    const material = new THREE.MeshPhysicalMaterial({
      color: env.color,
      transparent: false,
      roughness: env.roughness,
      opacity: env.opacity,
      metalness: env.metalness,
      emissive: env.emissive,
      reflectivity: env.reflectivity,
      clearcoat: env.clearcoat,
      ior: env.ior
    });

    console.log(coordinates);

    // edge from X to Y
    var direction = new THREE.Vector3().subVectors(coordinates[1], coordinates[0]);

    console.log(direction);
  
    // Make the geometry (of "direction" length)
    var geometry = new THREE.CylinderGeometry(radius, radius, direction.length(), 6, 4, false);
    // shift it so one end rests on the origin
    geometry.applyMatrix4(new THREE.Matrix4().makeTranslation(0, direction.length() / 2, 0));
    // rotate it the right way for lookAt to work
    geometry.applyMatrix4(new THREE.Matrix4().makeRotationX(THREE.MathUtils.degToRad(90)));
    // Make a mesh with the geometry
    var mesh = new THREE.Mesh(geometry, material);
    // Position it where we want

    mesh.position.copy(coordinates[0]);
    // And make it point to where we want
    mesh.lookAt(coordinates[1]); 

    env.mesh.add(mesh);

    //geometry.dispose();
    //material.dispose();
  };

  g3d.Tetrahedron = async (args, env) => {
    /**
     * @type {number[]}
     */
    var points = await interpretate(args[0], env);
    console.log("Points of tetrahedron:");
    console.log(points);
    var faces = [
      [points[0], points[1], points[2]],
      [points[0], points[1], points[3]],
      [points[1], points[2], points[3]],
      [points[0], points[3], points[2]],
    ];

    var fake = ["List"];

    var listVert = (cord) => ["List", cord[0], cord[1], cord[2]];

    faces.forEach((fs) => {
      fake.push([
        "Polygon",
        ["List", listVert(fs[0]), listVert(fs[1]), listVert(fs[2])],
      ]);
    });
    console.log(fake);
    await interpretate(fake, env);
  };

  g3d.Translate = async (args, env) => {
    let group = new THREE.Group();

    let p = await interpretate(args[1], env);

    //Backup of params
    let copy = Object.assign({}, env);
    copy.mesh = group;
    await interpretate(args[0], copy);

    group.translateX(p[0]);
    group.translateY(p[1]);
    group.translateZ(p[2]);

    env.local.mesh = group;

    env.mesh.add(group);
  };

  g3d.Translate.update = async (args, env) => {
    const p = await interpretate(args[1], env);
    const group = env.local.mesh;

    if (env.Lerp) {

      if (!env.local.lerp) {
        console.log('creating worker for lerp of movements..');
        const worker = {
          alpha: 0.05,
          target: new THREE.Vector3(...p),
          eval: () => {
            group.position.lerp(worker.target, 0.05)
          }
        };

        env.local.lerp = worker;  

        env.Handlers.push(worker);
      }

      env.local.lerp.target.fromArray(p);
      return;
    }

    group.translateX(p[0]);
    group.translateY(p[1]);
    group.translateZ(p[2]);
  }

  g3d.Translate.virtual = true  

  g3d.LookAt = async (args, env) => {
    const group = new THREE.Group();
    const dir = await interpretate(args[1], env);



    await interpretate(args[0], {...env, mesh:group});

    let bbox = new THREE.Box3().setFromObject(group);
    let center = bbox.max.clone().add(bbox.min).divideScalar(2);

    console.log('center: ');
    console.log(center);

    let translate = new THREE.Matrix4().makeTranslation(
      -center.x,
      -center.y,
      -center.z,
    );

    group.applyMatrix4(translate);

    group.lookAt(...dir);
    group.rotation.x = MathUtils.PI/2;

    translate = new THREE.Matrix4().makeTranslation(
      center.x,
      center.y,
      center.z,
    );

    group.applyMatrix4(translate);

    env.local.group = group;

    env.mesh.add(group);
  }

  g3d.LookAt.update = async (args, env) => {
    const dir = await interpretate(args[1], env);
    env.local.group.lookAt(...dir);
  }  

  g3d.LookAt.virtual = true

  g3d.GeometricTransformation = async (args, env) => {
    var group = new THREE.Group();
    //Если center, то наверное надо приметь matrix
    //к каждому объекту относительно родительской группы.
    var p = [...(await interpretate(args[1], {...env, hold:false}))];

    if (!p[0][0]) {
      //most likely this is Translate
      env.local.translationOnly = true;
    } else {
      //make it like Matrix4
      p.forEach((el) => {
        el.push(0);
      });
      p.push([0, 0, 0, 1]);
    }

    /*console.log(p);
    var centering = false;
    var centrans = [];

    if (p.length === 1) {
      p = p[0];
    }
    if (p.length === 1) {
      p = p[0];
    } else if (p.length === 2) {
      console.log(p);
      if (p[1] === "Center") {
        centering = true;
      } else {
        console.log("NON CENTERING ISSUE!!!");
        console.log(p);
        centrans = p[1];
        console.log("???");
      }
      //return;
      p = p[0];
    }

    if (p.length === 3) {
      if (typeof p[0] === "number") {
        var dir = p;
        var matrix = new THREE.Matrix4().makeTranslation(...dir, 1);
      } else {
        //make it like Matrix4
        p.forEach((el) => {
          el.push(0);
        });
        p.push([0, 0, 0, 1]);

        var matrix = new THREE.Matrix4();
        console.log("Apply matrix to group::");
        matrix.set(...aflatten(p));
      }
    } else {
      console.log(p);
      console.error("Unexpected length matrix: :: " + p);
    }

    //Backup of params
    var copy = Object.assign({}, env);
    copy.mesh = group;
    await interpretate(args[0], copy);
    console.log('MATRIX');
    console.log(matrix);

    if (centering || centrans.length > 0) {
      console.log("::CENTER::");
      var bbox = new THREE.Box3().setFromObject(group);
      console.log(bbox);
      var center = bbox.max.clone().add(bbox.min).divideScalar(2);
      if (centrans.length > 0) {
        console.log("CENTRANS");
        center = center.fromArray(centrans);
      }
      console.log(center);

      var translate = new THREE.Matrix4().makeTranslation(
        -center.x,
        -center.y,
        -center.z,
      );
      group.applyMatrix4(translate);
      group.applyMatrix4(matrix);
      translate = new THREE.Matrix4().makeTranslation(
        center.x,
        center.y,
        center.z
      );
      group.applyMatrix4(translate);
    } else {
      group.applyMatrix4(matrix);
    }*/

    await interpretate(args[0], {...env, mesh: group});

    let matrix;
    
    if (env.local.translationOnly) {
      matrix = new THREE.Matrix4().makeTranslation(...p);
    } else {
      matrix = new THREE.Matrix4().set(...aflatten(p));
    }

    group.matrixAutoUpdate = false;

    env.local.quaternion = new THREE.Quaternion();
    env.local.position = new THREE.Vector3();
    env.local.scale = new THREE.Vector3();    

    matrix.decompose(env.local.position, env.local.quaternion, env.local.scale);

    group.quaternion.copy( env.local.quaternion );
    group.position.copy( env.local.position );
    group.scale.copy( env.local.scale );

     // set initial values

    //group.quaternion.set(newQ);

    group.updateMatrix();

    env.local.group = group;

    env.mesh.add(group);
  };

  g3d.GeometricTransformation.update = async (args, env) => {
    let p = [...(await interpretate(args[1], {...env, hold:false}))];

    const group = env.local.group;
    
    if (!env.local.translationOnly) {
      p.forEach((el) => {
        el.push(0);
      });
      p.push([0, 0, 0, 1]);

      const matrix = new THREE.Matrix4().set(...aflatten(p));

      matrix.decompose(env.local.position, env.local.quaternion, env.local.scale); // set initial values
    }

    if (env.Lerp) {

      if (!env.local.lerp) {
        console.log('creating worker for lerp of matrix movements..');

        const worker = {
          alpha: 0.05,
          quaternion: env.local.quaternion.clone(), //target
          scale: env.local.scale.clone(), //target
          position: env.local.position.clone(), //target
          eval: () => {
            group.quaternion.slerp(worker.quaternion, worker.alpha); 
            group.position.lerp(worker.position, worker.alpha);
            group.updateMatrix();
          }
        };

        env.local.lerp = worker;  

        env.Handlers.push(worker);
      }

      if (!env.local.translationOnly)
        env.local.lerp.quaternion.copy(env.local.quaternion);
      else
        env.local.lerp.position.copy(env.local.position.set(...p));

      return;
    }

    if (!env.local.translationOnly) {
      group.quaternion.copy( env.local.quaternion );
      group.position.copy( env.local.position );
      group.scale.copy( env.local.scale );
    } else {
      group.position.copy( env.local.position.set(...p) );
    }

    group.updateMatrix();

    env.mesh.add(group);
  };  

  g3d.GeometricTransformation.virtual = true

  g3d.GraphicsComplex = async (args, env) => {
    var copy = Object.assign({}, env);

    copy.geometry = new THREE.Geometry();

    (await interpretate(args[0], copy)).forEach((el) => {
      if (typeof el[0] !== "number") console.error("not a triple of number" + el);
      copy.geometry.vertices.push(new THREE.Vector3(el[0], el[1], el[2]));
    });

    const group = new THREE.Group();

    await interpretate(args[1], copy);

    env.mesh.add(group);
    copy.geometry.dispose();
  };

  g3d.AbsoluteThickness = (args, env) => {}

  g3d.Polygon = async (args, env) => {
    if (env.hasOwnProperty("geometry")) {
      /**
       * @type {THREE.Geometry}
       */
      var geometry = env.geometry.clone();

      var createFace = (c) => {
        c = c.map((x) => x - 1);

        switch (c.length) {
          case 3:
            geometry.faces.push(new THREE.Face3(c[0], c[1], c[2]));
            break;

          case 4:
            geometry.faces.push(
              new THREE.Face3(c[0], c[1], c[2]),
              new THREE.Face3(c[0], c[2], c[3]),
            );
            break;

          case 5:
            geometry.faces.push(
              new THREE.Face3(c[0], c[1], c[4]),
              new THREE.Face3(c[1], c[2], c[3]),
              new THREE.Face3(c[1], c[3], c[4]),
            );
            break;
          /**
           * 0 1
           *5    2
           * 4  3
           */
          case 6:
            geometry.faces.push(
              new THREE.Face3(c[0], c[1], c[5]),
              new THREE.Face3(c[1], c[2], c[5]),
              new THREE.Face3(c[5], c[2], c[4]),
              new THREE.Face3(c[2], c[3], c[4])
            );
            break;
          default:
            console.log(c);
            console.log(c.length);
            console.error("Cant produce complex polygons! at", c);
        }
      };

      let a = await interpretate(args[0], env);
      if (a.length === 1) {
        a = a[0];
      }

      if (typeof a[0] === "number") {
        console.log("Create single face");
        createFace(a);
      } else {
        console.log("Create multiple face");
        a.forEach(createFace);
      }
    } else { 
      var geometry = new THREE.Geometry();
      let points = await interpretate(args[0], env);

      points.forEach((el) => {
        if (typeof el[0] !== "number") {
          console.error("not a triple of number", el);
          return;
        }
        geometry.vertices.push(new THREE.Vector3(el[0], el[1], el[2]));
      });

      console.log("points");
      console.log(points);

      switch (points.length) {
        case 3:
          geometry.faces.push(new THREE.Face3(0, 1, 2));
          break;

        case 4:
          geometry.faces.push(
            new THREE.Face3(0, 1, 2),
            new THREE.Face3(0, 2, 3));
          break;
        /**
         *  0 1
         * 4   2
         *   3
         */
        case 5:
          geometry.faces.push(
            new THREE.Face3(0, 1, 4),
            new THREE.Face3(1, 2, 3),
            new THREE.Face3(1, 3, 4));
          break;
        /**
         * 0  1
         *5     2
         * 4   3
         */
        case 6:
          geometry.faces.push(
            new THREE.Face3(0, 1, 5),
            new THREE.Face3(1, 2, 5),
            new THREE.Face3(5, 2, 4),
            new THREE.Face3(2, 3, 4)
          );
          break;
        default:
          console.log(points);
          console.error("Cant build complex polygon ::");
      }
    }

    const material = new THREE.MeshPhysicalMaterial({
      color: env.color,
      transparent: env.opacity < 0.9,
      opacity: env.opacity,
      roughness: env.roughness,
      metalness: env.metalness,
      emissive: env.emissive,
      reflectivity: env.reflectivity,
      clearcoat: env.clearcoat,
      ior: env.ior
      //depthTest: false
      //depthWrite: false
    });
    console.log(env.opacity);
    material.side = THREE.DoubleSide;

    geometry.computeFaceNormals();
    geometry.computeVertexNormals(true);
    const poly = new THREE.Mesh(geometry, material);

    //poly.frustumCulled = false;
    env.mesh.add(poly);
    material.dispose();
  };

  g3d.Polyhedron = async (args, env) => {
    if (args[1][1].length > 4) {
      //non-optimised variant to work with 4 vertex per face
      await interpretate(["GraphicsComplex", args[0], ["Polygon", args[1]]], env);
    } else {
      //reguar one. gpu-fiendly
      /**
       * @type {number[]}
       */
      const indices = await interpretate(args[1], env)
        .flat(4)
        .map((i) => i - 1);
      /**
       * @type {number[]}
       */
      const vertices = await interpretate(args[0], env).flat(4);

      const geometry = new THREE.PolyhedronGeometry(vertices, indices);

      var material = new THREE.MeshPhysicalMaterial({
        color: env.color,
        transparent: true,
        opacity: env.opacity,
        depthWrite: true,
        roughness: env.roughness,
        metalness: env.metalness,
        emissive: env.emissive,
        reflectivity: env.reflectivity,
        clearcoat: env.clearcoat,
        ior: env.ior
      });

      const mesh = new THREE.Mesh(geometry, material);
      env.mesh.add(mesh);
      geometry.dispose();
      material.dispose();
    }
  };

  g3d.GrayLevel = (args, env) => { };

  g3d.EdgeForm = (args, env) => { };

  g3d.Specularity = (args, env) => { };

  g3d.Text = (args, env) => { };

  g3d.Directive = (args, env) => { };

  g3d.PlaneGeometry = () => { new THREE.PlaneGeometry;  };

  g3d.Line = async (args, env) => {
    if (env.hasOwnProperty("geometry")) {
      const geometry = new THREE.Geometry();

      const points = await interpretate(args[0], env);
      points.forEach((el) => {
        geometry.vertices.push((env.geometry.vertices[el - 1]).clone(),);
      });

      const material = new THREE.LineBasicMaterial({
        linewidth: env.thickness,
        color: env.edgecolor,
      });
      const line = new THREE.Line(geometry, material);

      line.material.setValues({
        polygonOffset: true,
        polygonOffsetFactor: 1,
        polygonOffsetUnits: 1
      })

      env.mesh.add(line);

      geometry.dispose();
      material.dispose();
    } else {
      let arr = await interpretate(args[0], env);
      if (arr.length === 1) arr = arr[0];
      //if (arr.length !== 2) console.error( "Tube must have 2 vectors!");
      console.log("points: ", arr.length);

      const points = [];
      arr.forEach((p) => {
        points.push(new THREE.Vector4(...p, 1));
      });
      //new THREE.Vector4(...arr[0], 1)

      points.forEach((p) => {
        p = p.applyMatrix4(env.matrix);
      });

      const geometry = new THREE.Geometry().setFromPoints(points);
      const material = new THREE.LineBasicMaterial({
        color: env.edgecolor,
        linewidth: env.thickness,
      });

      env.mesh.add(new THREE.Line(geometry, material));
    }
  };

  let OrbitControls;
  let FirstPersonControls = false;

  let EffectComposer;
  let RenderPass; 
  let UnrealBloomPass;

  let GUI;

  g3d.ImageSize = () => "ImageSize"
  g3d.Background = () => "Background"
  g3d.AspectRatio = () => "AspectRatio"
  g3d.Lighting = () => "Lighting"
  g3d.Default = () => "Default"
  g3d.None = () => "None"
  g3d.Automatic = () => "Automatic"

  g3d.Graphics3D = async (args, env) => {
    /* lazy loading */


    THREE         = (await import('three'));
    OrbitControls = (await import("three/examples/jsm/controls/OrbitControls")).OrbitControls;
    EffectComposer= (await import('three/examples/jsm/postprocessing/EffectComposer')).EffectComposer;
    RenderPass    = (await import('three/examples/jsm/postprocessing/RenderPass')).RenderPass;
    UnrealBloomPass=(await import('three/examples/jsm/postprocessing/UnrealBloomPass')).UnrealBloomPass;
    GUI           = (await import('dat.gui')).GUI;


    /**
     * @type {Object}
     */   
    env.local.handlers = [];
    env.local.prolog   = [];

    const Handlers = [];

    /**
     * @type {Object}
     */  
    const options = await core._getRules(args, g3d);
    console.log(options);

    /**
     * @type {HTMLElement}
     */
    var container = env.element;

    /**
     * @type {[Number, Number]}
     */
    let ImageSize;
    
    if(options.ImageSize) {
      ImageSize = options.ImageSize;
    } else {
      ImageSize = [core.DefaultWidth, core.DefaultWidth*0.618034];
    } 

    let background = options.Background || new THREE.Color(0xffffff);

    const lighting = options.Lighting || "Default";

    const aspectratio = options.AspectRatio || 0.618034;

    //if only the width is specified
    if (!(ImageSize instanceof Array)) ImageSize = [ImageSize, ImageSize*aspectratio];
    console.log('Image size');
    console.log(ImageSize);

    //path tracing engine
    if (options.RTX) {
      //FullScreenQuad = (await import('THREE/examples/jsm/postprocessing/Pass.js')).Pass.FullScreenQuad;
      //RTX = await import('THREE-gpu-pathtracer/build/index.module');
      //PathTracingSceneGenerator   = RTX.PathTracingSceneGenerator;
      //PathTracingRenderer         = RTX.PathTracingRenderer;
      //PhysicalPathTracingMaterial = RTX.PhysicalPathTracingMaterial;
    }


    let controlObject = {
      init: (camera, dom) => {
        controlObject.o = new OrbitControls( camera, renderer.domElement );
        controlObject.o.target.set( 0, 1, 0 );
        controlObject.o.update();
      },

      dispose: () => {
        
      }
    };

    

    if (options.Controls) {

      if (options.Controls === 'PointerLockControls') {
        const o = (await import("three/examples/jsm/controls/PointerLockControls")).PointerLockControls;
        

        controlObject = {
          init: (camera, dom) => {
            controlObject.o = new o( camera, dom );
            env.local.scene.add( controlObject.o.getObject() );

            controlObject.onKeyDown = function ( event ) {
              switch ( event.code ) {
    
                case 'ArrowUp':
                case 'KeyW':
                  controlObject.moveForward = true;
                  break;
    
                case 'ArrowLeft':
                case 'KeyA':
                  controlObject.moveLeft = true;
                  break;
    
                case 'ArrowDown':
                case 'KeyS':
                  controlObject.moveBackward = true;
                  break;
    
                case 'ArrowRight':
                case 'KeyD':
                  controlObject.moveRight = true;
                  break;
    
                case 'Space':
                  if ( controlObject.canJump === true ) controlObject.velocity.y += 20;
                  controlObject.canJump = false;
                  break;
    
              }
              event.preventDefault();  
              event.returnValue = false;
              event.cancelBubble = true;
              return false;
              
    
            };
    
            controlObject.onKeyUp = function ( event ) {
    
              switch ( event.code ) {
    
                case 'ArrowUp':
                case 'KeyW':
                  controlObject.moveForward = false;
                  break;
    
                case 'ArrowLeft':
                case 'KeyA':
                  controlObject.moveLeft = false;
                  break;
    
                case 'ArrowDown':
                case 'KeyS':
                  controlObject.moveBackward = false;
                  break;
    
                case 'ArrowRight':
                case 'KeyD':
                  controlObject.moveRight = false;
                  break;
                
                
              }

              event.preventDefault();  
              event.returnValue = false;
              event.cancelBubble = true;
              return false;              
      
    
            };    
            

            env.local.handlers.push(controlObject.handler);

            const inst = document.createElement('div');
            inst.style.width="100%";
            inst.style.height="100%";
            inst.style.top = "0";
            inst.style.position = "absolute";


            
            env.element.appendChild(inst);

            renderer.domElement.addEventListener( 'keydown', controlObject.onKeyDown );
            renderer.domElement.addEventListener( 'keyup', controlObject.onKeyUp );

            renderer.domElement.tabIndex = 1;

            inst.addEventListener( 'click', function () {

              controlObject.o.lock();
    
            } );

            controlObject.o.addEventListener( 'lock', function () {

              inst.style.display = 'none';
              //blocker.style.display = 'none';
    
            } );
    
            controlObject.o.addEventListener( 'unlock', function () {
    
              //blocker.style.display = 'block';
              inst.style.display = '';
    
            } );

          },


          prevTime: performance.now(),

          handler: () => {
            const time = performance.now();
            const controls = controlObject.o;

            if ( controls.isLocked === true ) {
    
              //raycaster.ray.origin.copy( controls.getObject().position );
              //raycaster.ray.origin.y -= 10;
    
              //const intersections = raycaster.intersectObjects( objects, false );
    
              //const onObject = intersections.length > 0;
              const onObject = false;

              const delta = ( time - controlObject.prevTime ) / 1000;
    
              controlObject.velocity.x -= controlObject.velocity.x * 10.0 * delta;
              controlObject.velocity.z -= controlObject.velocity.z * 10.0 * delta;
    
              controlObject.velocity.y -= 9.8 * 4.0 * delta; // 100.0 = mass
    
              controlObject.direction.z = Number( controlObject.moveForward ) - Number( controlObject.moveBackward );
              controlObject.direction.x = Number( controlObject.moveRight ) - Number( controlObject.moveLeft );
              controlObject.direction.normalize(); // this ensures consistent movements in all directions
    
              if ( controlObject.moveForward || controlObject.moveBackward ) controlObject.velocity.z -= controlObject.direction.z * 40.0 * delta;
              if ( controlObject.moveLeft || controlObject.moveRight ) controlObject.velocity.x -= controlObject.direction.x * 40.0 * delta;
    
              if ( onObject === true ) {
    
                controlObject.velocity.y = Math.max( 0, controlObject.velocity.y );
                controlObject.canJump = true;
    
              }
    
              controls.moveRight( - controlObject.velocity.x * delta );
              controls.moveForward( - controlObject.velocity.z * delta );
    
              controls.getObject().position.y += ( controlObject.velocity.y * delta ); // new behavior
    
              if ( controls.getObject().position.y < 0.3 ) {
    
                controlObject.velocity.y = 0;
                controls.getObject().position.y = 0.3;
    
                controlObject.canJump = true;
    
              }
    
            }
    
            controlObject.prevTime = time;
          },

          moveBackward: false,
          moveForward: false,
          moveLeft: false,
          moveRight: false,
          canJump: false,
          velocity: new THREE.Vector3(),
          direction: new THREE.Vector3(),

          dispose: () =>{

            document.removeEventListener( 'keydown', controlObject.onKeyDown );
            document.removeEventListener( 'keyup', controlObject.onKeyUp );
            
          }  
        };

       

              

      } 
    }

    env.local.controlObject = controlObject;

    /**
    * @type {THREE.Mesh<THREE.Geometry>}
    */

    let camera, scene, renderer, composer;
    let controls, water, sun, mesh;

    const params = {
      exposure: 1,
      bloomStrength: 0.1,
      bloomThreshold: 0.5,
      bloomRadius: 0.11
    };

    await init();
    animate();
    

    function zoomExtents(scene, camera) {

    }

    async function init() {

      scene = new THREE.Scene();
      camera = new THREE.PerspectiveCamera( 55, ImageSize[0]/ImageSize[1], 0.01, 2000 );
      

      renderer = new THREE.WebGLRenderer();
      renderer.setPixelRatio( window.devicePixelRatio );
      renderer.setSize(ImageSize[0], ImageSize[1]);
      //renderer.toneMapping = THREE.ACESFilmicToneMapping;
      renderer.domElement.style = "margin:auto";
      container.appendChild( renderer.domElement );

      /* postprocess */
  		const renderScene = new RenderPass( scene, camera );

  		const bloomPass = new UnrealBloomPass( new THREE.Vector2( ImageSize[0], ImageSize[1] ), 1.5, 0.4, 0.85 );
  		bloomPass.threshold = params.bloomThreshold;
  		bloomPass.strength = params.bloomStrength;
  		bloomPass.radius = params.bloomRadius;

      composer = new EffectComposer( renderer );
  		composer.addPass( renderScene );
  		composer.addPass( bloomPass );

      composer.setSize(ImageSize[0], ImageSize[1]);

      function takeScheenshot() {
        //renderer.render( scene, camera );
        composer.render();
        renderer.domElement.toBlob(function(blob){
          var a = document.createElement('a');
          var url = URL.createObjectURL(blob);
          a.href = url;
          a.download = 'screenshot.png';
          a.click();
        }, 'image/png', 1.0);
      }

      const gui = new GUI({ autoPlace: false });
      const button = { Save:function(){ takeScheenshot() }};
      gui.add(button, 'Save');

      const bloomFolder = gui.addFolder('Bloom');

  		bloomFolder.add( params, 'exposure', 0.1, 2 ).onChange( function ( value ) {
  			renderer.toneMappingExposure = Math.pow( value, 4.0 );
  		} );

  		bloomFolder.add( params, 'bloomThreshold', 0.0, 1.0 ).onChange( function ( value ) {
  			bloomPass.threshold = Number( value );
  		} );

  		bloomFolder.add( params, 'bloomStrength', 0.0, 3.0 ).onChange( function ( value ) {
  			bloomPass.strength = Number( value );
  		} );

  		bloomFolder.add( params, 'bloomRadius', 0.0, 1.0 ).step( 0.01 ).onChange( function ( value ) {
  			bloomPass.radius = Number( value );
  		} );

      if (background instanceof THREE.Color) scene.background = background;


      const guiContainer = document.createElement('div');
      guiContainer.classList.add('graphics3d-controller');
      guiContainer.appendChild(gui.domElement);
      container.appendChild( guiContainer );    

      env.local.renderer = renderer;
      env.local.scene    = scene;
      env.local.camera = camera;

      if (lighting === "Default") g3d.DefaultLighting([], env);

      const group = new THREE.Group();

      const cameraMesh = {
        mesh: scene,
        pos: [3, 3, 10],
        set: (mesh, pos) => {
          cameraMesh.mesh = mesh;
          cameraMesh.pos = pos;
        }
      }

      const envcopy = {
        ...env,
        context: g3d,
        numerical: true,
        tostring: false,
        matrix: new THREE.Matrix4().set(
          1, 0, 0, 0,//
          0, 1, 0, 0,//
          0, 0, 1, 0,//
          0, 0, 0, 1),
        color: new THREE.Color(1, 1, 1),
        opacity: 1,
        thickness: 1,
        roughness: 0.5,
        edgecolor: new THREE.Color(0, 0, 0),
        mesh: group,
        metalness: 0,
        emissive: new THREE.Color(0, 0, 0),
        arrowHeight: 20,
        arrowRadius: 5,
        reflectivity: 0.5,
        clearcoat: 0,
        ior: 1.5,
        Lerp: options.Lerp,

        Handlers: Handlers,

        cameraMesh: cameraMesh
      }
    
      await interpretate(args[0], envcopy);

      envcopy.cameraMesh.mesh.add(camera);
      camera.position.set(...envcopy.cameraMesh.pos);

      group.applyMatrix4(new THREE.Matrix4().set(
        1, 0, 0, 0,
        0, 0, 1, 0,
        0,-1, 0, 0,
        0, 0, 0, 1));

      scene.add(group);

      
      
      env.local.controlObject.init( camera, renderer.domElement );
      zoomExtents(scene, camera);
    }


    function animate() {
      //added loop-handlers, void
      for (let i=0; i<Handlers.length; ++i) {
        Handlers[i].eval();
      }

      render();
      env.local.aid = requestAnimationFrame( animate );
    }

    function render() {
      //added loop-handlers, void
      env.local.handlers.forEach((f)=>{
        f();
      });

      //renderer.render( scene, camera );
      composer.render();
    }


  }; 

  let Water = false;
  let Sky   = false;

  g3d.Camera = (args, env) => {
    let pos = args;
    if (args.length === 0) pos = [3,3,1];
    env.cameraMesh.set(env.mesh, pos);
  }

  g3d.Reflectivity = (args, env) => {
    env.reflectivity = interpretate(args[0], env);
  }

  g3d.IOR = (args, env) => {
    env.ior = interpretate(args[0], env);
  }

  g3d.Clearcoat = (args, env) => {
    env.clearcoat = interpretate(args[0], env);
  }

  g3d.LightProbe = (args, env) => {
    //THREE.js light probe irradiance
  }

  g3d.DefaultLighting = (args, env) => {
    const lighting = [
      { type: "Ambient", color: [0.3, 0.2, 0.4] },
      {
        type: "Directional",
        color: [0.8, 0, 0],
        position: [2, 0, 2]
      },
      {
        type: "Directional",
        color: [0, 0.8, 0],
        position: [2, 2, 2]
      },
      {
        type: "Directional",
        color: [0, 0, 0.8],
        position: [0, 2, 2]
      }
    ];

    function addLight(l) {
      var color = new THREE.Color().setRGB(l.color[0], l.color[1], l.color[2]);
      var light;

      if (l.type === "Ambient") {
        light = new THREE.AmbientLight(color, 0.5);
      } else if (l.type === "Directional") {
        console.log('adding direction light');
        console.log(l);
        light = new THREE.DirectionalLight(color, 1);
        light.position.fromArray(l.position);

      } else if (l.type === "Spot") {
        light = new THREE.SpotLight(color);
        light.position.fromArray(l.position);
        light.target.position.fromArray(l.target);
        //light.target.updateMatrixWorld(); // This fixes bug in THREE.js
        light.angle = l.angle;
      } else if (l.type === "Point") {
        light = new THREE.PointLight(color);
        light.position.fromArray(l.position);

      } else {
        alert("Error: Internal Light Error", l.type);
        return;
      }
      return light;
    } 

    lighting.forEach((el) => env.local.camera.add(addLight(el)) );

  }

  g3d.SkyAndWater = async (args, env) => {
    if (!Water) {
      Water         = (await import('three/examples/jsm/objects/Water')).Water;
      Sky           = (await import('three/examples/jsm/objects/Sky')).Sky;  
    }

    let options = await core._getRules(args, env);
    console.log('options:');
    options.dims = options.Dims || [10000, 10000];
    options.skyscale = options.Skyscale || 10000;
    options.elevation = options.Elevation ||  8;
    options.azimuth = options.Azimuth || 180;

    options.turbidity = options.Turbidity || 10;
    options.rayleigh = options.Rayleigh || 2;
    options.mieCoefficient = options.MieCoefficient || 0.005;
    options.mieDirectionalG = options.MieDirectionalG || 0.8;

    console.log(options);

    let sun = new THREE.Vector3();
    let water;
    // Water

    const waterGeometry = new THREE.PlaneGeometry(...options.dims);

    water = new Water(
      waterGeometry,
      {
        textureWidth: 512,
        textureHeight: 512,
        waterNormals: new THREE.TextureLoader().load( 'https://cdn.statically.io/gh/JerryI/Mathematica-ThreeJS-graphics-engine/master/assets/waternormals.jpg', function ( texture ) {

          texture.wrapS = texture.wrapT = THREE.RepeatWrapping;

        } ),
        sunDirection: new THREE.Vector3(),
        sunColor: 0xffffff,
        waterColor: 0x001e0f,
        distortionScale: 3.7,
        fog: true
      }
    );

    water.rotation.x = - Math.PI / 2;
    
    env.local.water = water;

    // Skybox

    const sky = new Sky();
    sky.scale.setScalar( options.skyscale );

    env.local.sky = sky;  

    const skyUniforms = sky.material.uniforms;

    skyUniforms[ 'turbidity' ].value = options.turbidity;
    skyUniforms[ 'rayleigh' ].value = options.rayleigh;
    skyUniforms[ 'mieCoefficient' ].value = options.mieCoefficient;
    skyUniforms[ 'mieDirectionalG' ].value = options.mieDirectionalG;

    const parameters = {
      elevation: options.elevation,
      azimuth: options.azimuth
    };



    env.local.scene.add( water );
    env.local.scene.add( sky );

    const pmremGenerator = new THREE.PMREMGenerator( env.local.renderer );
    let renderTarget;

    const phi = THREE.MathUtils.degToRad( 90 - parameters.elevation );
    const theta = THREE.MathUtils.degToRad( parameters.azimuth );

    sun.setFromSphericalCoords( 1, phi, theta );

    sky.material.uniforms[ 'sunPosition' ].value.copy( sun );
    water.material.uniforms[ 'sunDirection' ].value.copy( sun ).normalize();

    if ( renderTarget !== undefined ) renderTarget.dispose();

    renderTarget = pmremGenerator.fromScene( sky );

    env.local.scene.environment = renderTarget.texture;  

    //every frame
    env.local.handlers.push(
      function() {
        env.local.water.material.uniforms[ 'time' ].value += 1.0 / 60.0;
      }
    );
  }

  g3d.Sky = async (args, env) => {
    if (!Sky) {
      Sky           = (await import('three/examples/jsm/objects/Sky')).Sky;  
    }

    let options = await core._getRules(args, env);
    console.log('options:');
    options.dims = options.Dims || [10000, 10000];
    options.skyscale = options.Skyscale || 10000;
    options.elevation = options.Elevation ||  8;
    options.azimuth = options.Azimuth || 180;

    options.turbidity = options.Turbidity || 10;
    options.rayleigh = options.Rayleigh || 2;
    options.mieCoefficient = options.MieCoefficient || 0.005;
    options.mieDirectionalG = options.MieDirectionalG || 0.8;

    console.log(options);

    let sun = new THREE.Vector3();

    // Skybox

    const sky = new Sky();
    sky.scale.setScalar( options.skyscale );

    env.local.sky = sky;  

    const skyUniforms = sky.material.uniforms;

    skyUniforms[ 'turbidity' ].value = options.turbidity;
    skyUniforms[ 'rayleigh' ].value = options.rayleigh;
    skyUniforms[ 'mieCoefficient' ].value = options.mieCoefficient;
    skyUniforms[ 'mieDirectionalG' ].value = options.mieDirectionalG;

    const parameters = {
      elevation: options.elevation,
      azimuth: options.azimuth
    };

    env.local.scene.add( sky );

    const pmremGenerator = new THREE.PMREMGenerator( env.local.renderer );
    let renderTarget;

    const phi = THREE.MathUtils.degToRad( 90 - parameters.elevation );
    const theta = THREE.MathUtils.degToRad( parameters.azimuth );

    sun.setFromSphericalCoords( 1, phi, theta );

    sky.material.uniforms[ 'sunPosition' ].value.copy( sun );
    env.local.sun = sun;
    //water.material.uniforms[ 'sunDirection' ].value.copy( sun ).normalize();

    if ( renderTarget !== undefined ) renderTarget.dispose();

    renderTarget = pmremGenerator.fromScene( sky );

    env.local.scene.environment = renderTarget.texture;  
  }
  
  g3d.Water = async (args, env) => {
    if (!Water) {
      Water         = (await import('three/examples/jsm/objects/Water')).Water;
    }

    let options = await core._getRules(args, env);
    console.log('options:');


    console.log(options);
    options.dims = options.Dims || [10000, 10000];

    let water;
    // Water

    const waterGeometry = new THREE.PlaneGeometry(...options.dims);

    water = new Water(
      waterGeometry,
      {
        textureWidth: 512,
        textureHeight: 512,
        waterNormals: new THREE.TextureLoader().load( 'https://cdn.statically.io/gh/JerryI/Mathematica-ThreeJS-graphics-engine/master/assets/waternormals.jpg', function ( texture ) {

          texture.wrapS = texture.wrapT = THREE.RepeatWrapping;

        } ),
        sunDirection: new THREE.Vector3(),
        sunColor: 0xffffff,
        waterColor: 0x001e0f,
        distortionScale: 3.7,
        fog: true
      }
    );

    water.rotation.x = - Math.PI / 2;
    
    env.local.water = water;

    env.local.scene.add( water );
    
    const sun = env.local.sun || (new THREE.Vector3(1,1,1));
    water.material.uniforms[ 'sunDirection' ].value.copy( sun ).normalize();

    //every frame
    env.local.handlers.push(
      function() {
        env.local.water.material.uniforms[ 'time' ].value += 1.0 / 60.0;
      }
    );
  }  

  g3d.Graphics3D.destroy = (args, env) => {
    console.log('Graphics3D was removed');
    console.log('env global:'); console.log(env.global);
    console.log('env local:'); console.log(env.local);
    env.local.controlObject.dispose();
    cancelAnimationFrame(env.local.aid);
  }

  g3d.Large = (args, env) => {
    return 1.0;
  }

}
