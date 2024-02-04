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
  "ViewProjection", "VertexTextureCoordinates", "RTX","ViewRange", "ViewVector", "ViewVertical", "Controls", "PointerLockControls", "VertexNormals", "VertexColors"].map((e)=>{
    g3d[e] = () => e;
  });

  g3d.Void = (args, env) => {console.warn(args); console.warn('went to the void...');};
  g3d.Void.update = () => {};
  g3d.Void.destroy = () => {};

  g3d.CapForm = g3d.Void;
  g3d.Appearance = g3d.Void;


  /**
  * @type {import('three')}
  */
  let THREE;
  let MathUtils;


  const lab2rgb = (lab) => {
    var y = (lab[0] + 16) / 116,
        x = lab[1] / 500 + y,
        z = y - lab[2] / 200,
        r, g, b;
  
    x = 0.95047 * ((x * x * x > 0.008856) ? x * x * x : (x - 16/116) / 7.787);
    y = 1.00000 * ((y * y * y > 0.008856) ? y * y * y : (y - 16/116) / 7.787);
    z = 1.08883 * ((z * z * z > 0.008856) ? z * z * z : (z - 16/116) / 7.787);
  
    r = x *  3.2406 + y * -1.5372 + z * -0.4986;
    g = x * -0.9689 + y *  1.8758 + z *  0.0415;
    b = x *  0.0557 + y * -0.2040 + z *  1.0570;
  
    r = (r > 0.0031308) ? (1.055 * Math.pow(r, 1/2.4) - 0.055) : 12.92 * r;
    g = (g > 0.0031308) ? (1.055 * Math.pow(g, 1/2.4) - 0.055) : 12.92 * g;
    b = (b > 0.0031308) ? (1.055 * Math.pow(b, 1/2.4) - 0.055) : 12.92 * b;
  
    return [Math.max(0, Math.min(1, r)) * 255, 
            Math.max(0, Math.min(1, g)) * 255, 
            Math.max(0, Math.min(1, b)) * 255]
  };

  g3d.LABColor =  async (args, env) => {
    let lab;
    if (args.length > 1)
      lab = [await interpretate(args[0], env), await interpretate(args[1], env), await interpretate(args[2], env)];
    else 
      lab = await interpretate(args[0], env);

    const color = lab2rgb(lab.map(e => e*255.0));
    console.log('LAB color');
    console.log(color);
    
    env.color = new THREE.Color(...(color.map(e => e/255.0)));
    if (args.length > 3) env.opacity = await interpretate(args[3], env);
    
    return env.color;   
  };

  g3d.LABColor.update = () => {};
  g3d.LABColor.destroy = () => {};

  g3d.Style = core.List;

  /**
   * @description https://threejs.org/docs/#api/en/materials/LineDashedMaterial
   */
  g3d.Dashing = (args, env) => {
    console.log("Dashing not implemented");
  };

  g3d.Annotation = core.List;

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
  };

  g3d.Emissive = async (args, env) => {
    const copy = {...env};
    await interpretate(args[0], copy);
    env.emissive = copy.color;
    if (args.length > 1) {
      env.emissiveIntensity = await interpretate(args[1], copy);
    }
  };


  g3d.Hue = async (args, env) => {
    const h = await interpretate(args[0], env);
    const l = await interpretate(args[1], env) * 100.0;
    const s = await interpretate(args[2], env) * 100.0;

    env.color = new THREE.Color(`hsl(${h}, ${l}%, ${s}%)`);
    return env.color;    
    
  };

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
    return env.color;
  };

  g3d.Roughness = (args, env) => {
    const o = interpretate(args[0], env);
    if (typeof o !== "number") console.error("Opacity must have number value!");
    console.log(o);
    env.roughness = o;  
  };

  g3d.Opacity = (args, env) => {
    var o = interpretate(args[0], env);
    if (typeof o !== "number") console.error("Opacity must have number value!");
    console.log(o);
    env.opacity = o;
  };

  g3d.ImageScaled = (args, env) => { };

  g3d.Thickness = async (args, env) => { env.thickness = await interpretate(args[0], env);
  };

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
    console.log(undefined);

    let radius = 1;
    if (args.length > 1) radius = await interpretate(args[1], env);
    /**
     * @type {THREE.Vector3}}
     */
    const coordinates = await interpretate(args[0], env);

    /**
     * @type {env.material}}
     */  
    const material = new env.material({
      color: env.color,
      transparent: false,
      roughness: env.roughness,
      opacity: env.opacity,
      metalness: env.metalness,
      emissive: env.emissive,
      emissiveIntensity: env.emissiveIntensity,
      
      
      
    });

    //points 1, 2
    const p2 = new THREE.Vector3(...coordinates[0]);
    const p1 = new THREE.Vector3(...coordinates[1]);
    //direction
    const dp = p2.clone().addScaledVector(p1, -1);

    const geometry = new THREE.CylinderGeometry(radius, radius, dp.length(), 32, 1);

    //calculate the center (might be done better, i hope BoundingBox doest not envolve heavy computations)
  

    //default geometry
    const cylinder = new THREE.Mesh(geometry, material);

    //cone
    const conegeometry = new THREE.ConeGeometry(env.arrowRadius, env.arrowHeight, 32 );
    const cone = new THREE.Mesh(conegeometry, material);
    cone.position.y = dp.length()/2 + env.arrowHeight/2;

    let group = new THREE.Group();
    group.add(cylinder, cone);


    var HALF_PI = Math.PI * .5;
    var position  = p1.clone().add(p2).divideScalar(2);

    var orientation = new THREE.Matrix4();//a new orientation matrix to offset pivot
    var offsetRotation = new THREE.Matrix4();//a matrix to fix pivot rotation
    new THREE.Matrix4();//a matrix to fix pivot position
    orientation.lookAt(p1,p2,new THREE.Vector3(0,1,0));//look at destination
    offsetRotation.makeRotationX(HALF_PI);//rotate 90 degs on X
    orientation.multiply(offsetRotation);//combine orientation with rotation transformations
    group.applyMatrix4(orientation);


    //group.position=position;    


    //translate its center to the middle target point
    group.position.addScaledVector(position, 1);

    env.mesh.add(group);

    geometry.dispose();
    conegeometry.dispose();
    material.dispose();

    return group;
  };

  g3d.Arrow = async (args, env) => {
    let arr;

    if (args.length === 1) {
      if (args[0][0] === 'Tube') {
        console.log('TUBE inside!');
        arr = await interpretate(args[0][1], env);
      } else {
        arr = await interpretate(args[0], env);
      }
    } else {
      arr = await interpretate(args[0], env);
    }
    
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
    arrowHelper.castShadow = env.shadows;
    arrowHelper.receiveShadow = env.shadows;
     
    if (env.PathRendering) {
      arrowHelper.line.material.emissive = env.emissive;
      arrowHelper.cone.material.emissive = env.emissive;
      arrowHelper.line.material.emissiveIntensity = env.emissiveIntensity;
      arrowHelper.cone.material.emissiveIntensity = env.emissiveIntensity;
    }

    env.mesh.add(arrowHelper);
    arrowHelper.line.material.linewidth = env.thickness;

    return arrowHelper;
  };

  g3d.Tube = g3d.Arrow;

  g3d.Sphere = async (args, env) => {
    var radius = 1;
    if (args.length > 1) radius = await interpretate(args[1], env);

    const material = new env.material({
      color: env.color,
      roughness: env.roughness,
      opacity: env.opacity,
      metalness: env.metalness,
      emissive: env.emissive,
emissiveIntensity: env.emissiveIntensity,
      
      
      
    });

    function addSphere(cr) {
      const origin = new THREE.Vector4(...cr, 1);
      const geometry = new THREE.SphereGeometry(radius, 40, 40);
      const sphere = new THREE.Mesh(geometry, material);

      sphere.position.set(origin.x, origin.y, origin.z);
      sphere.castShadow = env.shadows;
      sphere.receiveShadow = env.shadows;

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

    return env.local.object;
  };

  g3d.Sphere.update = async (args, env) => {
    console.log('Sphere: updating the data!');
    env.wake();

    const c = await interpretate(args[0], env);

    if (env.Lerp) {

      if (env.local.multiple) {

        if (!env.local.lerp) {
          console.log('creating worker for lerp of movements multiple..');
          const initial = c.map((e)=> new THREE.Vector3(...e));

          const worker = {
            alpha: 0.05,
            target: initial,
            eval: () => {
              for (let i=0; i<env.local.object.length; ++i)
                env.local.object[i].position.lerp(worker.target[i], 0.05);
            }
          };
  
          env.local.lerp = worker;  
  
          env.Handlers.push(worker);
        }
        
        for (let i=0; i<c.length; ++i)
          env.local.lerp.target[i].fromArray(c);

        return;

      } else {

        if (!env.local.lerp) {
          console.log('creating worker for lerp of movements..');
          const worker = {
            alpha: 0.05,
            target: new THREE.Vector3(...c),
            eval: () => {
              env.local.object.position.lerp(worker.target, 0.05);
            }
          };
  
          env.local.lerp = worker;  
  
          env.Handlers.push(worker);
        }
  
        env.local.lerp.target.fromArray(c);
        return;
      }
    }

    if (env.local.multiple) {
      let i = 0;
      c.forEach((cc)=>{
        env.local.object[i].position.set(...cc);
        ++i;
      });

      return;
    }

    
    env.local.object.position.set(...c);

  };

  g3d.Sphere.destroy = async (args, env) => {
    console.log('Sphere: destroy');
    console.log(args);
    console.log(env);

    for (const a of args)
      await interpretate(a, env);

  };

  g3d.Sphere.virtual = true;

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
  };

  g3d._Water = (args, env) => {
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
  };

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
        new THREE.Vector4(...(await interpretate(args[1], env)), 1),
        new THREE.Vector4(...(await interpretate(args[0], env)), 1),
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
    const material = new env.material({
      color: env.color,
      transparent: true,
      opacity: env.opacity,
      roughness: env.roughness,
      depthWrite: true,
      metalness: env.metalness,
      emissive: env.emissive,
emissiveIntensity: env.emissiveIntensity,
      
      
      
    });

    //material.side = THREE.DoubleSide;

    const cube = new THREE.Mesh(geometry, material);

    //var tr = new THREE.Matrix4();
    //	tr.makeTranslation(origin.x,origin.y,origin.z);

    //cube.applyMatrix(params.matrix.clone().multiply(tr));

    cube.position.set(origin.x, origin.y, origin.z);
    cube.receiveShadow = env.shadows;
    cube.castShadow = env.shadows;

    env.mesh.add(cube);

    geometry.dispose();
    material.dispose();

    return cube;
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

    const material = new env.material({
      color: env.color,
      transparent: false,
      roughness: env.roughness,
      opacity: env.opacity,
      metalness: env.metalness,
      emissive: env.emissive,
emissiveIntensity: env.emissiveIntensity,
      
      
      
    });

    console.log(coordinates);

    // edge from X to Y
    var direction = new THREE.Vector3().subVectors(coordinates[1], coordinates[0]);

    console.log(direction);
  
    // Make the geometry (of "direction" length)
    var geometry = new THREE.CylinderGeometry(radius, radius, direction.length(), 32, 4, false);
    // shift it so one end rests on the origin
    geometry.applyMatrix4(new THREE.Matrix4().makeTranslation(0, direction.length() / 2, 0));
    // rotate it the right way for lookAt to work
    geometry.applyMatrix4(new THREE.Matrix4().makeRotationX(THREE.MathUtils.degToRad(90)));
    // Make a mesh with the geometry
    var mesh = new THREE.Mesh(geometry, material);
    // Position it where we want
    mesh.receiveShadow = env.shadows;
    mesh.castShadow = env.shadows;

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
    
    return await interpretate(fake, env);
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
    env.wake();
    const p = await interpretate(args[1], env);
    const group = env.local.mesh;

    if (env.Lerp) {

      if (!env.local.lerp) {
        console.log('creating worker for lerp of movements..');
        const worker = {
          alpha: 0.05,
          target: new THREE.Vector3(...p),
          eval: () => {
            group.position.lerp(worker.target, 0.05);
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
  };

  g3d.Translate.virtual = true;  

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
  };

  g3d.LookAt.update = async (args, env) => {
    env.wake();
    const dir = await interpretate(args[1], env);
    env.local.group.lookAt(...dir);
  };  

  g3d.LookAt.virtual = true;

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
    env.wake();
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

  g3d.GeometricTransformation.virtual = true;

  g3d.GraphicsComplex = async (args, env) => {
    var copy = Object.assign({}, env);
    const options = await core._getRules(args, {...env, hold: true});
    

    const pts = (await interpretate(args[0], copy)).flat();
    copy.vertices = new Float32Array( pts );

    if ('VertexColors' in options) {
      copy.vertexcolors = await interpretate(options["VertexColors"], env);
    }

    const group = new THREE.Group();

    await interpretate(args[1], copy);

    env.mesh.add(group);
    //copy.geometry.dispose();
  };

  g3d.AbsoluteThickness = (args, env) => {};

  g3d.Polygon = async (args, env) => {
    var geometry = new THREE.BufferGeometry();
    let vertices;

    if (env.hasOwnProperty("vertices")) {
      vertices = env.vertices;

      let a = await interpretate(args[0], env);
      
      if (a[0].length === 3) {
        geometry.setIndex( a.flat().map((e)=>e-1) );
      } else {
        //more complicatec case, need to covert all polygons into triangles
        let extendedIndexes = [];

        console.log(a);

        if (Array.isArray(a[0])) {
       

        for (let i=0; i<a.length; ++i) {
          const b = a[i];
          switch (b.length) {
            case 3:
              extendedIndexes.push([b[0],b[1],b[2]]);
              break;
    
            case 4:
              extendedIndexes.push([b[0],b[1],b[2]]);
              extendedIndexes.push([b[0],b[2],b[3]]);
              break;
            /**
             *  0 1
             * 4   2
             *   3
             */
            case 5:
              extendedIndexes.push([b[0], b[1], b[4]]);
              extendedIndexes.push([b[1], b[2], b[3]]);
              extendedIndexes.push([b[1], b[3], b[4]]);
              break;
            /**
             * 0  1
             *5     2
             * 4   3
             */
            case 6:
              extendedIndexes.push([b[0], b[1], b[5]]);
              extendedIndexes.push([b[1], b[2], b[5]]);
              extendedIndexes.push([b[5], b[2], b[4]]);
              extendedIndexes.push([b[2], b[3], b[4]]);
              break;
            default:
             
              console.error("cannot build complex polygon");
              //FIXME
              break;
          }
        }   
      } else {
        extendedIndexes = a.flat();
      }

        geometry.setIndex( extendedIndexes.flat().map((e)=>e-1) );
        
        
      }


    } else { 
      
      let points = await interpretate(args[0], env);

      vertices = new Float32Array( points.flat() );

  

      console.log("points");
      console.log(points);

      switch (points.length) {
        case 3:
          geometry.setIndex( [0,1,2] );
          break;

        case 4:
          geometry.setIndex( [0,1,2,      0,2,3] );
          break;
        /**
         *  0 1
         * 4   2
         *   3
         */
        case 5:
          geometry.setIndex( [0, 1, 4,     1, 2, 3,     1, 3, 4] );
          break;
        /**
         * 0  1
         *5     2
         * 4   3
         */
        case 6:
          geometry.setIndex( [0, 1, 5,     1, 2, 5,     5, 2, 4,       2, 3, 4] );
          break;
        default:
          console.log(points);
          console.error("Cant build complex polygon ::");
      }
    }

    geometry.setAttribute( 'position', new THREE.BufferAttribute( vertices, 3 ) );
    geometry.computeVertexNormals();

    let material;
    
    if (env.vertexcolors) {
      console.log('vertex colors');
      geometry.setAttribute( 'color', new THREE.Float32BufferAttribute( new Float32Array( env.vertexcolors.flat() ), 3 ) );

      console.log(env.material);
      //!!! a bug. Cannot work with PhysicalBasedMaterial
      material = new THREE.MeshBasicMaterial({
        vertexColors: true,
        transparent: env.opacity < 1,
        opacity: env.opacity,
        roughness: env.roughness,
        metalness: env.metalness,
        emissive: env.emissive,
        emissiveIntensity: env.emissiveIntensity,        
      });
  } else {
    material = new env.material({
      color: env.color,
      transparent: env.opacity < 1,
      opacity: env.opacity,
      roughness: env.roughness,
      metalness: env.metalness,
      emissive: env.emissive,
      emissiveIntensity: env.emissiveIntensity,
      
      
      
      //depthTest: false
      //depthWrite: false
    });
  }

    

    console.log(env.opacity);
    material.side = THREE.DoubleSide;

    const poly = new THREE.Mesh(geometry, material);
    poly.receiveShadow = env.shadows;
    poly.castShadow = true;

    //poly.frustumCulled = false;
    env.mesh.add(poly);
    material.dispose();

    return poly;
  };

  g3d.Polyhedron = async (args, env) => {
    if (args[1][1].length > 4) {
      //non-optimised variant to work with 4 vertex per face
      return await interpretate(["GraphicsComplex", args[0], ["Polygon", args[1]]], env);
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

      var material = new env.material({
        color: env.color,
        transparent: true,
        opacity: env.opacity,
        depthWrite: true,
        roughness: env.roughness,
        metalness: env.metalness,
        emissive: env.emissive,
emissiveIntensity: env.emissiveIntensity,
        
        
        
      });

      const mesh = new THREE.Mesh(geometry, material);
      mesh.receiveShadow = env.shadows;
      mesh.castShadow = env.shadows;
      env.mesh.add(mesh);
      geometry.dispose();
      material.dispose();

      return mesh;
    }
  };

  g3d.GrayLevel = (args, env) => { };

  g3d.EdgeForm = (args, env) => { };

  g3d.Specularity = (args, env) => { };

  g3d.Text = (args, env) => { };

  g3d.Directive = (args, env) => { };

  g3d.PlaneGeometry = () => { new THREE.PlaneGeometry;  };


  g3d.Line = async (args, env) => {
    
    var geometry = new THREE.BufferGeometry();
    let vertices;

    if (env.hasOwnProperty("vertices")) {
      
      vertices = env.vertices;

      let a = await interpretate(args[0], env);
      

      geometry.setIndex( a.flat().map((e)=>e-1) );
     
      geometry.setAttribute( 'position', new THREE.BufferAttribute( vertices, 3 ) );


    } else { 
   
      const points = await interpretate(args[0], env);

      const pts = [];

      for (let i=0; i<points.length; ++i)
        pts.push(new THREE.Vector3().fromArray(points[i]));

        console.log(points);
   

      geometry.setFromPoints( pts );
    }

    
    

      const material = new THREE.LineBasicMaterial({
        linewidth: env.thickness,
        color: env.edgecolor,
        opacity: env.opacity,
        transparent: env.opacity < 1.0 ? true : false
      });
      const line = new THREE.Line(geometry, material);

 
      env.mesh.add(line);

      return line;

      //geometry.dispose();
      //material.dispose();
  };

  let GUI;

  g3d.ImageSize = () => "ImageSize";
  g3d.Background = () => "Background";
  g3d.AspectRatio = () => "AspectRatio";
  g3d.Lighting = () => "Lighting";
  g3d.Default = () => "Default";
  g3d.None = () => "None";
  g3d.Lightmap = () => "Lightmap";
  g3d.Automatic = () => "Automatic"; 

  let Water = false;
  let Sky   = false;

  g3d.Camera = (args, env) => {
    console.warn('temporary disabled');
    return;
  };

  g3d.Reflectivity = (args, env) => {
    env.reflectivity = interpretate(args[0], env);
  };

  g3d.IOR = (args, env) => {
    env.ior = interpretate(args[0], env);
  };

  g3d.Clearcoat = (args, env) => {
    env.clearcoat = interpretate(args[0], env);
  };

  g3d.LightProbe = (args, env) => {
    //THREE.js light probe irradiance
  };

  g3d.DefaultLighting = (args, env) => {
    console.warn('temporary disabled');
    return;

  };

  g3d.SkyAndWater = async (args, env) => {
    console.warn('temporary disabled');
    return;
  };

  g3d.Sky = async (args, env) => {
    console.warn('temporary disabled');
    return;
  };
  
  g3d.Water = async (args, env) => {
    
    
    if (!Water) {
      Water         = (await import('./Water-a573d7d4.js')).Water;
    }

    let options = await core._getRules(args, env);
    console.log('options:');


    console.log(options);
    options.dims = options.Size || [10000, 10000];

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
        sunDirection: new THREE.Vector3(1,1,1),
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
  };  


  g3d.Large = (args, env) => {
    return 1.0;
  };

const setImageSize = async (options, env) => {
  let ImageSize;

  if (options.ImageSize) {
    ImageSize = await interpretate(options.ImageSize, env);
    if (!(ImageSize instanceof Array)) ImageSize = [ImageSize, ImageSize*0.618034];
  } else {
    ImageSize = [core.DefaultWidth, core.DefaultWidth*0.618034];
  }

  return ImageSize;
};

let RTX = {};

const addDefaultLighting = (scene) => {
  const light = new THREE.PointLight(0xffffff, 2, 10);
  light.position.set(0, 10, 0);
  scene.add(light);
  var hemiLight = new THREE.HemisphereLight( 0xffffbb, 0x080820, 2 );
  scene.add( hemiLight );
};

g3d.PointLight = async (args, env) => {
  const copy = {...env};
  const options = await core._getRules(args, {...env, hold: true});

  console.log(options);

  let color = 0xffffff; if (args.length - options.length > 0) color = await interpretate(args[0], copy); 
  let intensity = 1; if (args.length - options.length > 1) intensity = await interpretate(args[1], env);
  let distance = 0; if (args.length - options.length > 2) distance = await interpretate(args[2], env);
  let decay = 2; if (args.length - options.length > 3) decay = await interpretate(args[3], env);

  let position = [0, 0, 10];
  if (options.Position) {
    position = await interpretate(options.Position, env);
  }

  console.log(position);
  console.log(color);

  const light = new THREE.PointLight(color, intensity, distance, decay);
  light.castShadow = env.shadows;
  light.position.set(...position);
  env.local.light = light;
  env.mesh.add(light);

  return light;
};

g3d.PointLight.update = async (args, env) => {
  env.wake();
  const options = await core._getRules(args, {...env, hold: true}); 

  if (options.Position) {
    const pos = await interpretate(options.Position, env);

    if (env.Lerp) {
        if (!env.local.lerp) {
          
          console.log('creating worker for lerp of movements..');
          const worker = {
            alpha: 0.05,
            target: new THREE.Vector3(...pos),
            eval: () => {
              env.local.light.position.lerp(worker.target, 0.05);
            }
          };
  
          env.local.lerp = worker;  
  
          env.Handlers.push(worker);
        }
  
        env.local.lerp.target.fromArray(pos);
        return;
    } 

    env.local.light.position.set(...pos);
  }  
};

g3d.PointLight.destroy = async (args, env) => {for (const i of args) await interpretate(i, env);};

g3d.PointLight.virtual = true;

g3d.SpotLight = async (args, env) => {
  const copy = {...env};
  const options = await core._getRules(args, {...env, hold: true});

  console.log(options);

  let color = 0xffffff; if (args.length - options.length > 0) color = await interpretate(args[0], copy); 
  let intensity = 1; if (args.length - options.length > 1) intensity = await interpretate(args[1], env);
  let distance = 0; if (args.length - options.length > 2) distance = await interpretate(args[2], env);
  let angle = Math.PI/3; if (args.length - options.length > 3) angle = await interpretate(args[3], env);
  let penumbra = 0; if (args.length - options.length > 4) penumbra = await interpretate(args[4], env);
  let decay = 2; if (args.length - options.length > 5) decay = await interpretate(decay[5], env);

  let position = [10, 10, 100];
  if (options.Position) {
    position = await interpretate(options.Position, env);
  }

  console.log(position);

  let target = [0,0,0];
  if (options.Target) {
    target = await interpretate(options.Target, env);
  }

  console.log(target);  

  const spotLight = new THREE.SpotLight( color, intensity, distance, angle, penumbra, decay );
  spotLight.position.set(...position);
  spotLight.target.position.set(...target);

  spotLight.castShadow = env.shadows;

  env.local.spotLight = spotLight;
  env.mesh.add(spotLight);
  env.mesh.add(spotLight.target);

  return spotLight;
};

g3d.SpotLight.update = async (args, env) => {
  env.wake();
  const options = await core._getRules(args, {...env, hold: true}); 

  if (options.Position) {
    const pos = await interpretate(options.Position, env);

    if (env.Lerp) {
        if (!env.local.lerp1) {
          
          console.log('creating worker for lerp of movements..');
          const worker = {
            alpha: 0.05,
            target: new THREE.Vector3(...pos),
            eval: () => {
              env.local.spotLight.position.lerp(worker.target, 0.05);
            }
          };
  
          env.local.lerp1 = worker;  
  
          env.Handlers.push(worker);
        }
  
        env.local.lerp1.target.fromArray(pos);
        
    } else {

      env.local.spotLight.position.set(...pos);
    }
  }
  
  if (options.Target) {
    const target = await interpretate(options.Target, env);

    if (env.Lerp) {
        if (!env.local.lerp2) {
          
          console.log('creating worker for lerp of movements..');
          const worker = {
            alpha: 0.05,
            target: new THREE.Vector3(...target),
            eval: () => {
              env.local.spotLight.target.position.lerp(worker.target, 0.05);
            }
          };
  
          env.local.lerp2 = worker;  
  
          env.Handlers.push(worker);
        }
  
        env.local.lerp2.target.fromArray(target);
        
    } else {

      env.local.spotLight.target.position.set(...target);
    }
  }  

};

g3d.SpotLight.destroy = async (args, env) => {for (const i of args) await interpretate(i, env);};

g3d.SpotLight.virtual = true;

g3d.Shadows = async (args, env) => {
  env.shadows = await interpretate(args[0], env);
};

g3d.Shadows.destroy = g3d.Shadows;

g3d.HemisphereLight = async (args, env) => {
  const copy = {...env};

  if (args.length > 0) await interpretate(args[0], copy); else copy.color = 0xffffbb;
  const skyColor = copy.color;

  if (args.length > 1) await interpretate(args[1], copy); else copy.color = 0x080820;
  const groundColor = copy.color;
 if (args.length > 2) await interpretate(args[1], env);

  const hemiLight = new THREE.HemisphereLight( skyColor, groundColor, 2 );
  env.local.scene.add( hemiLight );
};

g3d.MeshMaterial = async (args, env) => {
  const mat = await interpretate(args[0], env);
  env.material = mat;
};

g3d.MeshPhysicalMaterial = () => THREE.MeshPhysicalMaterial;
g3d.MeshLambertMaterial = () => THREE.MeshLambertMaterial;
g3d.MeshPhongMaterial = () => THREE.MeshPhongMaterial;
g3d.MeshToonMaterial = () => THREE.MeshToonMaterial;

let TransformControls = false;

g3d.EventListener = async (args, env) => {
    const rules = await interpretate(args[1], env);

    const copy = {...env};

    const object = await interpretate(args[0], env);

    if (!TransformControls) TransformControls = (await import('./TransformControls-2614e9c1.js')).TransformControls;
    rules.forEach((rule)=>{
      g3d.EventListener[rule.lhs](rule.rhs, object, copy);
    });

    return null;
};

  g3d.EventListener.destroy = (args, env) => {interpretate(args[0], env);};

  g3d.EventListener.transform = (uid, object, env) => {
    console.log(env.camera);
    const control = new TransformControls(env.camera, env.local.renderer.domElement);

    const orbit = env.controlObject.o;

    control.attach(object); 

    env.local.scene.add(control); 

    const updateData = throttle((x,y,z) => {
      server.kernel.emitt(uid, `<|"position"->{${x}, ${y}, ${z}}|>`, 'transform');
    });

    control.addEventListener( 'change', function(event) {
      updateData(object.position.x,object.position.y,object.position.z);
    } );

    control.addEventListener( 'dragging-changed', function ( event ) {
      console.log('changed');
      orbit.enabled = !event.value;
    } );
  };

let RGBELoader;
let OrbitControls;
let FullScreenQuad;

core.Graphics3D = async (args, env) => {  
  //Lazy loading

  THREE         = (await import('./three.module-e66a3903.js'));
  OrbitControls = (await import('./OrbitControls-07d4c9ef.js')).OrbitControls;
  GUI           = (await import('./dat.gui.module-6246c119.js')).GUI;  
  RGBELoader    = (await import('./RGBELoader-132ece2d.js')).RGBELoader; 
  FullScreenQuad = (await import('./Pass-bef7cacb.js')).FullScreenQuad; 
  MathUtils     = THREE.MathUtils;

  let sleeping = false;
  let timeStamp = performance.now();

  /**
   * @type {Object}
   */  
  let options = await core._getRules(args, {...env, context: g3d, hold:true});
  console.log(options);  

  if (Object.keys(options).length === 0 && args.length > 1) {
    options = await core._getRules(args[1], {...env, context: g3d, hold:true});
  }


  let PathRendering = false;
  if ('RTX' in options) {
    PathRendering = true;
    RTX = (await import('./index.module-c8bf6d87.js'));
  }


    /**
     * @type {Object}
     */   
    env.local.handlers = [];
    env.local.prolog   = [];

    const Handlers = [];

  /**
   * @type {HTMLElement}
   */
  const container = env.element;

  /**
   * @type {[Number, Number]}
   */
  const ImageSize = await setImageSize(options, env); 

  const params = 	{
    multipleImportanceSampling: true,
  	stableNoise: false,
  	denoiseEnabled: true,
  	denoiseSigma: 2.5,
  	denoiseThreshold: 0.1,
  	denoiseKSigma: 1.0,
  	environmentIntensity: 1,
  	environmentRotation: 0,
  	environmentBlur: 0.0,
  	backgroundBlur: 0.0,
  	bounces: 5,
  	samplesPerFrame: 1,
  	acesToneMapping: true,
  	resolutionScale: 0.5,
  	transparentTraversals: 20,
  	filterGlossyFactor: 0.5,
  	tiles: 1,
  	backgroundAlpha: 1,
  	checkerboardTransparency: true,
  	cameraProjection: 'Orthographic',
  };

  if (options.ViewProjection) { 
    params.cameraProjection = await interpretate(options.ViewProjection, env);
  }

  if (!PathRendering) params.resolutionScale = 1;

  //Setting GUI
  const gui = new GUI({ autoPlace: false, name: '...', closed:true });

  const guiContainer = document.createElement('div');
  guiContainer.classList.add('graphics3d-controller');
  guiContainer.appendChild(gui.domElement);
  container.appendChild( guiContainer );    

  function takeScheenshot() {
    animateOnce();
    renderer.domElement.toBlob(function(blob){
      var a = document.createElement('a');
      var url = URL.createObjectURL(blob);
      a.href = url;
      a.download = 'screenshot.png';
      a.click();
    }, 'image/png', 1.0);
  }

  const button = { Save:function(){ takeScheenshot(); }};
  gui.add(button, 'Save');

  //Setting up renderer
  let renderer, controls, ptRenderer, activeCamera, blitQuad, denoiseQuad;
  let perspectiveCamera, orthoCamera, equirectCamera;
  let envMap, envMapGenerator, scene;
  let PT_PROGRAM_ID;

  const orthoWidth = 5;

	renderer = new THREE.WebGLRenderer( { antialias: true } );
	renderer.toneMapping = THREE.ACESFilmicToneMapping;
	renderer.outputEncoding = THREE.sRGBEncoding;
	renderer.setClearColor( 0, 0 );
	container.appendChild( renderer.domElement );

	const aspect = ImageSize[0]/ImageSize[1];

  if (PathRendering) {
	  perspectiveCamera = new RTX.PhysicalCamera( 75, aspect, 0.025, 500 );
	  perspectiveCamera.position.set( - 4, 2, 3 );
  } else {
    perspectiveCamera = new THREE.PerspectiveCamera( 75, aspect, 0.025, 500 );

    renderer.shadowMap.enabled = true;
  }

  const wakeFunction = () => {
    timeStamp = performance.now();
    if (!sleeping) return;
    env.local.wakeThreadUp(); 
  };

  

	const orthoHeight = orthoWidth / aspect;
	orthoCamera = new THREE.OrthographicCamera( orthoWidth / - 2, orthoWidth / 2, orthoHeight / 2, orthoHeight / - 2, 0, 100 );
	orthoCamera.position.set( - 4, 2, 3 );

  activeCamera = orthoCamera;

  if (PathRendering) {
	  equirectCamera = new RTX.EquirectCamera();
	  equirectCamera.position.set( - 4, 2, 3 );

	  ptRenderer = new RTX.PathTracingRenderer( renderer );
	  ptRenderer.alpha = true;
	  ptRenderer.material = new RTX.PhysicalPathTracingMaterial();
	  ptRenderer.material.setDefine( 'TRANSPARENT_TRAVERSALS', params.transparentTraversals );
	  ptRenderer.material.setDefine( 'FEATURE_MIS', Number( params.multipleImportanceSampling ) );
	  ptRenderer.tiles.set( params.tiles, params.tiles );

	  blitQuad = new FullScreenQuad( new THREE.MeshBasicMaterial( {
		  map: ptRenderer.target.texture,
		  blending: THREE.CustomBlending,
		  premultipliedAlpha: renderer.getContextAttributes().premultipliedAlpha,
	  } ) );

	  denoiseQuad = new FullScreenQuad( new RTX.DenoiseMaterial( {
		  map: ptRenderer.target.texture,
		  blending: THREE.CustomBlending,
		  premultipliedAlpha: renderer.getContextAttributes().premultipliedAlpha,
	  } ) ); 
  } 

  let controlObject = {
    init: (camera, dom) => {
      controlObject.o = new OrbitControls( camera, renderer.domElement );
      controlObject.o.addEventListener('change', wakeFunction);
      controlObject.o.target.set( 0, 1, 0 );
      controlObject.o.update();
    },

    dispose: () => {
      
    }
  };

  

  if (options.Controls) {

    if ((await interpretate(options.Controls, env)) === 'PointerLockControls') {
      const o = (await import('./PointerLockControls-20643468.js')).PointerLockControls;
      

      controlObject = {
        init: (camera, dom) => {
          controlObject.o = new o( camera, dom );
          env.local.scene.add( controlObject.o.getObject() );
          controlObject.o.addEventListener('change', wakeFunction);

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

            const delta = ( time - controlObject.prevTime ) / 1000;
  
            controlObject.velocity.x -= controlObject.velocity.x * 10.0 * delta;
            controlObject.velocity.z -= controlObject.velocity.z * 10.0 * delta;
  
            controlObject.velocity.y -= 9.8 * 4.0 * delta; // 100.0 = mass
  
            controlObject.direction.z = Number( controlObject.moveForward ) - Number( controlObject.moveBackward );
            controlObject.direction.x = Number( controlObject.moveRight ) - Number( controlObject.moveLeft );
            controlObject.direction.normalize(); // this ensures consistent movements in all directions
  
            if ( controlObject.moveForward || controlObject.moveBackward ) controlObject.velocity.z -= controlObject.direction.z * 40.0 * delta;
            if ( controlObject.moveLeft || controlObject.moveRight ) controlObject.velocity.x -= controlObject.direction.x * 40.0 * delta;
  
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



  controlObject.init(activeCamera, renderer.domElement);
  controls = controlObject.o;

  env.local.controlObject = controlObject;
  env.local.renderer = renderer;

  if (PathRendering) {
	  controls.addEventListener( 'change', () => {
	  	ptRenderer.reset();
	  } ); 
  } 

	scene = new THREE.Scene();

  const group = new THREE.Group();



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
    material: THREE.MeshPhysicalMaterial,
    color: new THREE.Color(1, 1, 1),
    opacity: 1,
    thickness: 1,
    roughness: 0.5,
    edgecolor: new THREE.Color(0, 0, 0),
    mesh: group,
    metalness: 0,
    emissive: new THREE.Color(0, 0, 0),
    emissiveIntensity: 1,
    arrowHeight: 20,
    arrowRadius: 5,
    reflectivity: 0.5,
    clearcoat: 0,
    shadows: false,
    Lerp: options.Lerp || true,
    camera: activeCamera,
    controlObject: controlObject,

    Handlers: Handlers,
    wake: wakeFunction
  };  

  env.local.wakeThreadUp = () => {
    if (!sleeping) return;
    sleeping = false;
    console.warn("g3d >> waking up!");
    animate();
  };

  env.local.renderer = renderer;
  env.local.scene    = scene;
  envcopy.camera   = activeCamera;
  env.local.element  = container;

  if (PathRendering)
    envcopy.PathRendering = true;

  await interpretate(args[0], envcopy);

  var bbox = new THREE.Box3().setFromObject(group);

  group.position.set(-(bbox.min.x + bbox.max.x) / 2, -(bbox.min.y + bbox.max.y) / 2, 0);

  if (options.BoxRatios || options.Boxed) {
    const boxLine = [
      [[bbox.min.x, bbox.min.y, bbox.min.z], [bbox.max.x, bbox.min.y, bbox.min.z], [bbox.max.x, bbox.max.y, bbox.min.z], [bbox.min.x, bbox.max.y, bbox.min.z], [bbox.min.x, bbox.min.y, bbox.min.z]],
      [[bbox.min.x, bbox.min.y, bbox.max.z], [bbox.max.x, bbox.min.y, bbox.max.z], [bbox.max.x, bbox.max.y, bbox.max.z], [bbox.min.x, bbox.max.y, bbox.max.z], [bbox.min.x, bbox.min.y, bbox.max.z]],
      [[bbox.min.x, bbox.min.y, bbox.min.z], [bbox.min.x, bbox.min.y, bbox.max.z]],
      [[bbox.max.x, bbox.min.y, bbox.min.z], [bbox.max.x, bbox.min.y, bbox.max.z]],
      [[bbox.max.x, bbox.max.y, bbox.min.z], [bbox.max.x, bbox.max.y, bbox.max.z]],
      [[bbox.min.x, bbox.max.y, bbox.min.z], [bbox.min.x, bbox.max.y, bbox.max.z]]
    ];

    for (const l of boxLine) {
      await interpretate(['Line', ['JSObject', l]], {...envcopy, opacity:0.5});
    }  }
  

  group.applyMatrix4(new THREE.Matrix4().set( 
    1, 0, 0, 0,
    0, 0, 1, 0,
    0, -1, 0, 0,
    0, 0, 0, 1));



    

  scene.add(group);

  scene.updateMatrixWorld();

  //add some lighting
  if ('Lighting' in options) ; else {
    addDefaultLighting(scene);
  }
  
  if (PathRendering)
	  envMapGenerator = new RTX.BlurredEnvMapGenerator( renderer ); 

	let envMapPromise;
  
  if ('Lightmap' in options) {
    const url = await interpretate(options.Lightmap, env);

    envMapPromise = new RGBELoader().setDataType( THREE.FloatType )
		.loadAsync(url)
		.then( texture => {

      if (PathRendering) {
        envMap = texture;
        updateEnvBlur();
      }

      if (PathRendering) return;

      const localEnv = pmremGenerator.fromEquirectangular( texture ).texture;

      scene.environment = localEnv;

      scene.background = localEnv;

      texture.dispose();
      pmremGenerator.dispose();

		} );
  }
  

    if (!PathRendering) {
      var pmremGenerator = new THREE.PMREMGenerator( renderer );
      pmremGenerator.compileEquirectangularShader();
    }

    if (PathRendering) {
      var generator = new RTX.PathTracingSceneGenerator( scene );
      var sceneInfo = generator.generate( scene );
      var { bvh, textures, materials } = sceneInfo;

      var geometry = bvh.geometry;
      var material = ptRenderer.material;

      material.bvh.updateFrom( bvh );
      material.attributesArray.updateFrom(
        geometry.attributes.normal,
        geometry.attributes.tangent,
        geometry.attributes.uv,
        geometry.attributes.color,
      );

    material.materialIndexAttribute.updateFrom( geometry.attributes.materialIndex );
    material.textures.setTextures( renderer, 2048, 2048, textures );
    material.materials.updateFrom( materials, textures );
  }

  if ('Lightmap' in options)
    await Promise.all( [ envMapPromise ] );    



  function onResize() {

    const w = ImageSize[0];
    const h = ImageSize[1];
    const scale = params.resolutionScale;
    const dpr = window.devicePixelRatio;

    if (PathRendering) {
      ptRenderer.setSize( w * scale * dpr, h * scale * dpr );
      ptRenderer.reset();
    }
  
    renderer.setSize( w, h );
    renderer.setPixelRatio( window.devicePixelRatio * scale );
  
    const aspect = w / h;
    
    perspectiveCamera.aspect = aspect;
    perspectiveCamera.updateProjectionMatrix();
  
    const orthoHeight = orthoWidth / aspect;
    orthoCamera.top = orthoHeight / 2;
    orthoCamera.bottom = orthoHeight / - 2;
    orthoCamera.updateProjectionMatrix();
  
  }
  
  function reset() {
    if (PathRendering)
      ptRenderer.reset();
  }
  
  function updateEnvBlur() {
  
    const blurredTex = envMapGenerator.generate( envMap, params.environmentBlur );
    ptRenderer.material.envMapInfo.updateFrom( blurredTex );
    scene.environment = blurredTex;
    ptRenderer.reset();
  
  }
  
  function updateCamera( cameraProjection ) {
  
    if ( cameraProjection === 'Perspective' ) {
  
      if ( activeCamera ) {
  
        perspectiveCamera.position.copy( activeCamera.position );
  
      }
  
      activeCamera = perspectiveCamera;
  
    } else if ( cameraProjection === 'Orthographic' ) {
  
      if ( activeCamera ) {
  
        orthoCamera.position.copy( activeCamera.position );
  
      }
  
      activeCamera = orthoCamera;
  
    } 
  
    controls.object = activeCamera;
    if (PathRendering)
      ptRenderer.camera = activeCamera;
  
    controls.update();

    env.local.camera   = activeCamera;
    envcopy.camera   = activeCamera;


  
    reset();
  
  }

  let animate;

  if (PathRendering) {
    animate = () => {
      animateOnce();
      env.local.aid = requestAnimationFrame( animate );
    };
  } else {
    animate = () => {
      animateOnce();

      if (performance.now() - timeStamp > 1000) {
        sleeping = true;
        console.warn('g3d >> Sleeping...');
      } else {
        env.local.aid = requestAnimationFrame( animate );
      }

    };    
  }

  
  function updateSettings() {
    wakeFunction();

    if (PathRendering) {
      ptRenderer.material.materials.updateFrom( sceneInfo.materials, sceneInfo.textures );

  
      ptRenderer.material.filterGlossyFactor = params.filterGlossyFactor;
      ptRenderer.material.environmentIntensity = params.environmentIntensity;
      ptRenderer.material.backgroundBlur = params.backgroundBlur;
      ptRenderer.material.bounces = params.bounces;
      ptRenderer.material.backgroundAlpha = params.backgroundAlpha;
      ptRenderer.material.physicalCamera.updateFrom( activeCamera );
    }
  
    activeCamera.updateMatrixWorld();
  
    if ( params.backgroundAlpha < 1.0 ) {
  
      scene.background = null;
  
    } else {
  
      scene.background = scene.environment;
  
    }
  }

  function animateOnce() {
    
    if (PathRendering) {
      activeCamera.updateMatrixWorld();
      // Get the path tracing shader id. It will be the next program compiled and used here.
      if ( PT_PROGRAM_ID === undefined ) {
  
        PT_PROGRAM_ID = renderer.info.programs.length;
  
      }
  
      for ( let i = 0, l = params.samplesPerFrame; i < l; i ++ ) {
  
          ptRenderer.update();
  
      }
  

      denoiseQuad.material.sigma = params.denoiseSigma;
      denoiseQuad.material.threshold = params.denoiseThreshold;
      denoiseQuad.material.kSigma = params.denoiseKSigma;
  
      const quad = params.denoiseEnabled ? denoiseQuad : blitQuad;
  
      renderer.autoClear = false;
      quad.material.map = ptRenderer.target.texture;
      quad.render( renderer );
      renderer.autoClear = true;

    } else {
      renderer.render( scene, activeCamera );
    }

    for (let i=0; i<Handlers.length; ++i) {
      //if (Handlers[i].sleep) continue;
      Handlers[i].eval();
    }

    //added loop-handlers, void
    env.local.handlers.forEach((f)=>{
      f();
    });    
    /**/

    //env.wake();
  
    //samplesEl.innerText = `Samples: ${ Math.floor( ptRenderer.samples ) }`;
  
  }
  

	onResize();

	updateCamera( params.cameraProjection );

  if (PathRendering) {
	  const ptFolder = gui.addFolder( 'Path Tracing' );

	ptFolder.add( params, 'stableNoise' ).onChange( value => {

		ptRenderer.stableNoise = value;
    updateSettings();

	} );
	ptFolder.add( params, 'multipleImportanceSampling' ).onChange( value => {

		ptRenderer.material.setDefine( 'FEATURE_MIS', Number( value ) );
		ptRenderer.reset();
    updateSettings();

	} );
	ptFolder.add( params, 'tiles', 1, 4, 1 ).onChange( value => {

		ptRenderer.tiles.set( value, value );
    updateSettings();

	} );
	ptFolder.add( params, 'samplesPerFrame', 1, 10, 1 );
	ptFolder.add( params, 'filterGlossyFactor', 0, 1 ).onChange( () => {

		ptRenderer.reset();
    updateSettings();

	} );
	ptFolder.add( params, 'bounces', 1, 30, 1 ).onChange( () => {

		ptRenderer.reset();
    updateSettings();

	} );
	ptFolder.add( params, 'transparentTraversals', 0, 40, 1 ).onChange( value => {

		ptRenderer.material.setDefine( 'TRANSPARENT_TRAVERSALS', value );
		ptRenderer.reset();
    updateSettings();

	} );
	ptFolder.add( params, 'resolutionScale', 0.1, 1 ).onChange( () => {

		onResize();

	} );

	const denoiseFolder = gui.addFolder( 'Denoising' );
	denoiseFolder.add( params, 'denoiseEnabled' ).onChange(updateSettings);
	denoiseFolder.add( params, 'denoiseSigma', 0.01, 12.0 ).onChange(updateSettings);
	denoiseFolder.add( params, 'denoiseThreshold', 0.01, 1.0 ).onChange(updateSettings);
	denoiseFolder.add( params, 'denoiseKSigma', 0.0, 12.0 ).onChange(updateSettings);
	denoiseFolder.close();


	const envFolder = gui.addFolder( 'Environment' );
	envFolder.add( params, 'environmentIntensity', 0, 10 ).onChange( () => {

		ptRenderer.reset();
    updateSettings();

	} );
	envFolder.add( params, 'environmentRotation', 0, 2 * Math.PI ).onChange( v => {

		ptRenderer.material.environmentRotation.makeRotationY( v );
		ptRenderer.reset();
    updateSettings();

	} );
	envFolder.add( params, 'environmentBlur', 0, 1 ).onChange( () => {

		updateEnvBlur();
    updateSettings();

	} );
	envFolder.add( params, 'backgroundBlur', 0, 1 ).onChange( () => {

		ptRenderer.reset();
    updateSettings();

	} );
	envFolder.add( params, 'backgroundAlpha', 0, 1 ).onChange( () => {

		ptRenderer.reset();
    updateSettings();

	} );
	envFolder.close();
  }

	const cameraFolder = gui.addFolder( 'Camera' );
	cameraFolder.add( params, 'cameraProjection', [ 'Perspective', 'Orthographic' ] ).onChange( v => {

		updateCamera( v );
    updateSettings();

	} );

  cameraFolder.add( params, 'acesToneMapping' ).onChange( value => {

    renderer.toneMapping = value ? THREE.ACESFilmicToneMapping : THREE.NoToneMapping;
    if (PathRendering)
      blitQuad.material.needsUpdate = true;

    updateSettings();

  } );

  if (PathRendering) {
	cameraFolder.add( perspectiveCamera, 'focusDistance', 1, 100 ).onChange( reset );
	cameraFolder.add( perspectiveCamera, 'apertureBlades', 0, 10, 1 ).onChange( function ( v ) {

		perspectiveCamera.apertureBlades = v === 0 ? 0 : Math.max( v, 3 );
		this.updateDisplay();
		reset();
    updateSettings();

	} );
	cameraFolder.add( perspectiveCamera, 'apertureRotation', 0, 12.5 ).onChange( reset );
	cameraFolder.add( perspectiveCamera, 'anamorphicRatio', 0.1, 10.0 ).onChange( reset );
	cameraFolder.add( perspectiveCamera, 'bokehSize', 0, 50 ).onChange( reset ).listen();
	cameraFolder.add( perspectiveCamera, 'fStop', 0.3, 20 ).onChange( reset ).listen();
	cameraFolder.add( perspectiveCamera, 'fov', 25, 100 ).onChange( () => {

		perspectiveCamera.updateProjectionMatrix();
		reset();
    updateSettings();

	} ).listen();
  }

	cameraFolder.close();  

  animate();
};

core.Graphics3D.destroy = (args, env) => {
  console.log('Graphics3D was removed');
  console.log('env global:'); console.log(env.global);
  console.log('env local:'); console.log(env.local);
  env.local.wakeThreadUp = () => {};
  env.local.controlObject.dispose();
  cancelAnimationFrame(env.local.aid);
};
