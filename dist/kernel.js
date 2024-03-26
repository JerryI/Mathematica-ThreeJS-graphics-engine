var node = {};

Object.defineProperty(node, "__esModule", {
  value: true
});
var default_1 = node.default = void 0;
const t1 = 6 / 29;
const t2 = 3 * t1 * t1;

const lrgb2rgb = x => Math.round(255 * (x <= 0.0031308 ? 12.92 * x : 1.055 * Math.pow(x, 1 / 2.4) - 0.055)) || 0;

const lab2xyz = t => t > t1 ? t * t * t : t2 * (t - 4 / 29);

var _default = ({
  luminance,
  a,
  b
}) => {
  const baseY = (luminance + 16) / 116;
  const x = 0.96422 * lab2xyz(baseY + a / 500);
  const y = Number(lab2xyz(baseY));
  const z = 0.82521 * lab2xyz(baseY - b / 200);
  return {
    red: lrgb2rgb(3.1338561 * x - 1.6168667 * y - 0.4906146 * z),
    green: lrgb2rgb(-0.9787684 * x + 1.9161415 * y + 0.0334540 * z),
    blue: lrgb2rgb(0.0719453 * x - 0.2289914 * y + 1.4052427 * z)
  };
};

default_1 = node.default = _default;

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

  g3d.LABColor =  async (args, env) => {
    let lab;
    if (args.length > 1)
      lab = [await interpretate(args[0], env), await interpretate(args[1], env), await interpretate(args[2], env)];
    else 
      lab = await interpretate(args[0], env);

      const color = default_1({luminance: 100*lab[0], a: 100*lab[1], b: 100*lab[2]});
    console.log('LAB color');
    console.log(color);
    
    env.color = new THREE.Color(color.red / 255.0, color.green / 255.0, color.blue / 255.0);
    if (args.length > 3) env.opacity = await interpretate(args[3], env);
    
    return env.color;   
  };

  g3d.LABColor.update = () => {};


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
    const h = await interpretate(args[0], env) * 360.0;
    const l = await interpretate(args[1], env) * 100.0;
    const s = await interpretate(args[2], env) * 100.0;

    env.color = new THREE.Color(`hsl(${h}, ${l}%, ${s}%)`);
    return env.color;    
    
  };

  g3d.EdgeForm = async (args, env) => {
    env.edgecolor = await interpretate(args[0], {...env});
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

  g3d.GrayLevel = async (args, env) => { 
    const r = await interpretate(args[0], env);

    env.color = new THREE.Color(r, r, r);
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

  g3d.AbsoluteThickness = async (args, env) => { env.thickness = await interpretate(args[0], env);
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
    
    env.local.matrix = group.matrix.clone();
    group.applyMatrix4(orientation);


    //group.position=position;    


    //translate its center to the middle target point
    group.position.addScaledVector(position, 1);

    env.local.group = group;

    env.mesh.add(group);

    geometry.dispose();
    conegeometry.dispose();
    material.dispose();

    return group;
  };

  g3d.TubeArrow.update = async (args, env) => {
    /**
     * @type {THREE.Vector3}}
     */
    
    
    const coordinates = await interpretate(args[0], env);
    //points 1, 2
    const p2 = new THREE.Vector3(...coordinates[0]);
    const p1 = new THREE.Vector3(...coordinates[1]);
    //direction
    p2.clone().addScaledVector(p1, -1);

    //const geometry = new THREE.CylinderGeometry(radius, radius, dp.length(), 32, 1);

    //calculate the center (might be done better, i hope BoundingBox doest not envolve heavy computations)
  

    //default geometry
    //const cylinder = new THREE.Mesh(geometry, material);

    //cone
    //const conegeometry = new THREE.ConeGeometry(env.arrowRadius, env.arrowHeight, 32 );
    //const cone = new THREE.Mesh(conegeometry, material);
    //cone.position.y = dp.length()/2 + env.arrowHeight/2;

    ///let group = new THREE.Group();
    //group.add(cylinder, cone);


    var HALF_PI = Math.PI * .5;
    var position  = p1.clone().add(p2).divideScalar(2);

    var orientation = new THREE.Matrix4();//a new orientation matrix to offset pivot
    var offsetRotation = new THREE.Matrix4();//a matrix to fix pivot rotation
    new THREE.Matrix4();//a matrix to fix pivot position
    orientation.lookAt(p1,p2,new THREE.Vector3(0,1,0));//look at destination
    offsetRotation.makeRotationX(HALF_PI);//rotate 90 degs on X
    orientation.multiply(offsetRotation);//combine orientation with rotation transformations

    env.local.matrix.decompose( env.local.group.position, env.local.group.quaternion, env.local.group.scale );
    env.local.group.matrix.copy( env.local.matrix );

    env.local.group.applyMatrix4(orientation);


    //group.position=position;    


    //translate its center to the middle target point
    env.local.group.position.addScaledVector(position, 1);

    env.wake();
  };

  g3d.TubeArrow.virtual = true; 

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
 

    if (arr.length > 2) {
      var geometry = new THREE.BufferGeometry();
      const points = arr.slice(0, -1);

      geometry.setAttribute( 'position', new THREE.BufferAttribute( new Float32Array(points.flat()), 3 ) );

      const material = new THREE.LineBasicMaterial({
        linewidth: env.thickness,
        color: env.color,
        opacity: env.opacity,
        transparent: env.opacity < 1.0 ? true : false
      });
      const line = new THREE.Line(geometry, material);

      env.local.line = line;
 
      env.mesh.add(line);
    }

    const points = [
      new THREE.Vector4(...arr[arr.length-2], 1),
      new THREE.Vector4(...arr[arr.length-1], 1),
    ];

    points.forEach((p) => {
      p = p.applyMatrix4(env.matrix);
    });

    const origin = points[0].clone();
    const dir = points[1].add(points[0].negate());
    const len = dir.length();

    const arrowHelper = new THREE.ArrowHelper(
      dir.normalize(),
      origin,
      len,
      env.color
    );
    //arrowHelper.castShadow = env.shadows;
    //arrowHelper.receiveShadow = env.shadows;
     
    if (env.PathRendering) {
      arrowHelper.line.material.emissive = env.emissive;
      arrowHelper.cone.material.emissive = env.emissive;
      arrowHelper.line.material.emissiveIntensity = env.emissiveIntensity;
      arrowHelper.cone.material.emissiveIntensity = env.emissiveIntensity;
    }

    env.mesh.add(arrowHelper);
    arrowHelper.line.material.linewidth = env.thickness;

    env.local.arrow = arrowHelper;

    return arrowHelper;
  };

  g3d.Arrow.update = async (args, env) => {
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

    if (env.local.line) {
      //update line geometry
      const positionAttribute = env.local.line.geometry.getAttribute( 'position' );
      const points = arr.slice(0, -1);

      positionAttribute.needsUpdate = true;
  
      for ( let i = 0; i < positionAttribute.count; i ++ ) {
        positionAttribute.setXYZ( i, ...(points[i]));
      }
  
      env.local.line.geometry.computeBoundingBox();
      env.local.line.geometry.computeBoundingSphere();
    }

    const points = [
      new THREE.Vector4(...arr[arr.length-2], 1),
      new THREE.Vector4(...arr[arr.length-1], 1),
    ];

    points.forEach((p) => {
      p = p.applyMatrix4(env.matrix);
    });


    env.local.arrow.position.copy(points[0]);

    const dir = points[1].add(points[0].negate());

    const len = dir.length();

    env.local.arrow.setDirection(dir.normalize());
    env.local.arrow.setLength(len);

    env.wake();

  };

  g3d.Arrow.destroy = async (args, env) => {
    if (env.local.line) env.local.line.dispose();
    env.local.arrow.dispose();
  };

  g3d.Arrow.virtual = true;

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

    console.log(env.local);
    let list = await interpretate(args[0], env);
    console.log('DRAW A SPHERE');

    if (list.length === 3) {
      env.local.object = [addSphere(list)];
    } else {

      //env.local.multiple = true;
      env.local.object = [];

      list.forEach((el) => {
        env.local.object.push(addSphere(el));
      });
    } 

    material.dispose();

    return env.local.object;
  };

  g3d.Sphere.update = async (args, env) => {
    console.log('Sphere: updating the data!');
    env.wake();

    let c = await interpretate(args[0], env);

    if (env.local.object.length == 1) {
      c = [c];
    }

    if (env.Lerp) {

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
          env.local.lerp.target[i].fromArray(c[i]);

        return;

  
    }

    {
      let i = 0;
      c.forEach((cc)=>{
        env.local.object[i].position.set(...cc);
        ++i;
      });

      return;
    }

  };

  g3d.Sphere.destroy = async (args, env) => {
    console.log('Sphere: destroy');
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

    //env.local.prev = [diff.x, diff.y, diff.z];

    const geometry = new THREE.BoxGeometry(1, 1, 1);
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

    env.local.geometry = cube.geometry.clone();
    cube.geometry.applyMatrix4(new THREE.Matrix4().makeScale(diff.x, diff.y, diff.z));

    cube.receiveShadow = env.shadows;
    cube.castShadow = env.shadows;

    env.mesh.add(cube);

    env.local.box = cube;

    geometry.dispose();
    material.dispose();

    return cube;
  };

  g3d.Cuboid.update = async (args, env) => {
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
    } else {
      p = await interpretate(args[0], env);
      origin = new THREE.Vector4(...p, 1);
      diff = new THREE.Vector4(1, 1, 1, 1);
    
      //shift it
      origin.add(diff.clone().divideScalar(2));
    }


    console.log(diff.x, diff.y, diff.z);

    env.local.box.position.copy(origin);
    env.local.box.geometry.copy(env.local.geometry);
    env.local.box.geometry.applyMatrix4(new THREE.Matrix4().makeScale(diff.x, diff.y, diff.z));

    //env.local.box.updateMatrix();

 
    env.wake();

  }; 

  g3d.Cuboid.destroy = async (args, env) => {
    env.local.box.geometry.dispose();
  };

  g3d.Cuboid.virtual = true;


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
    var geometry = new THREE.CylinderGeometry(radius, radius, 1, 32, 4, false);
    // shift it so one end rests on the origin
    geometry.applyMatrix4(new THREE.Matrix4().makeTranslation(0, 1 / 2.0, 0));
    // rotate it the right way for lookAt to work
    geometry.applyMatrix4(new THREE.Matrix4().makeRotationX(THREE.MathUtils.degToRad(90)));
    // Make a mesh with the geometry

    //env.local.geometry = geometry.clone();

    //geometry.applyMatrix4(new THREE.Matrix4().makeScale(1, 1, direction.length()));


    var mesh = new THREE.Mesh(geometry, material);
    // Position it where we want
    mesh.receiveShadow = env.shadows;
    mesh.castShadow = env.shadows;

    //env.local.bmatrix = mesh.matrix.clone();

    mesh.position.copy(coordinates[0]);

    env.local.g = mesh.geometry.clone();
    mesh.geometry.applyMatrix4(new THREE.Matrix4().makeScale(1,1,direction.length()));
    //mesh.scale.set( 1,1,direction.length() );

    // And make it point to where we want
    mesh.geometry.lookAt(direction); 

    

    env.local.cylinder = mesh;
    //env.local.coordinates = coordinates;
    //mesh.matrixAutoUpdate = false;

    env.mesh.add(mesh);

    //geometry.dispose();
    //material.dispose();
  };

  g3d.Cylinder.update = async (args, env) => {
    let coordinates = await interpretate(args[0], env);
    if (coordinates.length === 1) {
      coordinates = coordinates[0];
    }

    coordinates[0] = new THREE.Vector3(...coordinates[0]);
    coordinates[1] = new THREE.Vector3(...coordinates[1]);   
    
    var direction = new THREE.Vector3().subVectors(coordinates[1], coordinates[0]);


    env.local.cylinder.position.copy(coordinates[0]);

    //env.local.cylinder.matrix.identity();
    env.local.cylinder.geometry.copy(env.local.g);
    env.local.cylinder.geometry.applyMatrix4(new THREE.Matrix4().makeScale(1,1,direction.length()));

    //env.local.cylinder.applyMatrix4(new THREE.Matrix4().makeScale(1, 1, direction.length()));
    //env.local.cylinder.scale.set( 1,1,direction.length() );
    // And make it point to where we want
    env.local.cylinder.geometry.lookAt(direction); 

    env.wake();

  };

  g3d.Cylinder.destroy = async (args, env) => {
    env.local.cylinder.geometry.dispose();
    env.local.g.dispose();
  };

  g3d.Cylinder.virtual = true;

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


  const decodeTransformation = (arrays, env) => {

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
    let matrix = [];

    if (!env.local.type) {
      if (arrays.length == 2) {
        console.warn('apply matrix3x3 + translation');
        //translation matrix + normal 3x3
        env.local.type = 'complex';
      } else {
        if (!Array.isArray(arrays[0])) {
          //most likely this is Translate
          console.warn('apply translation');
          env.local.type = 'translation';

        } else {
          env.local.type = 'normal';
          console.warn('apply matrix 3x3');
        }
      }
    }

    switch(env.local.type) {
      case 'normal':
        //make it like Matrix4

        arrays.forEach((el) => {
          el.push(0);
        });

        matrix = arrays;
        matrix.push([0, 0, 0, 1]);
        matrix = new THREE.Matrix4().set(...aflatten(matrix));
      break;

      case 'translation':
 
        matrix = new THREE.Matrix4().makeTranslation(...arrays);
      break;

      case 'complex':
        matrix = arrays[0];
        const v = arrays[1];
  
        matrix[0].push(v[0]);
        matrix[1].push(v[1]);
        matrix[2].push(v[2]);
  
        matrix.push([0, 0, 0, 1]);
        matrix = new THREE.Matrix4().set(...aflatten(matrix));
      break;

      default:
        throw 'undefined type of matrix or vector';
    }

    return matrix;
  };

  g3d.GeometricTransformation = async (args, env) => {  
    const data = await interpretate(args[1], env);

    if (data.length > 3) {
      //list of matrixes
      console.warn('multiple matrixes');
      env.local.entities = [];

      for (const m of data) {
        const group = new THREE.Group();
        const matrix = decodeTransformation(m, env);
  
        await interpretate(args[0], {...env, mesh: group});
  
        group.matrixAutoUpdate = false;
        
        const object = {};

        object.quaternion = new THREE.Quaternion();
        object.position = new THREE.Vector3();
        object.scale = new THREE.Vector3();    
    
        matrix.decompose(object.position, object.quaternion, object.scale);
    
        group.quaternion.copy( object.quaternion );
        group.position.copy( object.position );
        group.scale.copy( object.scale );
    
        group.updateMatrix();
    
        object.group = group;
    
        env.mesh.add(group);
        env.local.entities.push(object);
      }
  

      return env.local.entities[0];

    } else {
      console.warn('single matrix');

      const group = new THREE.Group();
      const matrix = decodeTransformation(data, env);

      await interpretate(args[0], {...env, mesh: group});

      group.matrixAutoUpdate = false;
  
      env.local.quaternion = new THREE.Quaternion();
      env.local.position = new THREE.Vector3();
      env.local.scale = new THREE.Vector3();    
  
      matrix.decompose(env.local.position, env.local.quaternion, env.local.scale);
  
      group.quaternion.copy( env.local.quaternion );
      group.position.copy( env.local.position );
      group.scale.copy( env.local.scale );
  
      group.updateMatrix();
  
      env.local.group = group;
  
      env.mesh.add(group);

      return group;
    }
    
  };

  g3d.GeometricTransformation.update = async (args, env) => {
    env.wake();
    const data = await interpretate(args[1], env);
  
    if (env.local.entities) {
      //list of matrixes
      console.log('multiple matrixes');

      for (let i =0; i<env.local.entities.length; ++i) {
        const group = env.local.entities[i].group;

        const matrix = decodeTransformation(data[i], env);
  
        //await interpretate(args[0], {...env, mesh: group});
  
        
 

        const quaternion = new THREE.Quaternion();
        const position = new THREE.Vector3();
        const scale = new THREE.Vector3();    
    
        matrix.decompose(position, quaternion, scale);
    
        group.quaternion.copy( quaternion );
        group.position.copy( position );
        group.scale.copy( scale );
    
        group.updateMatrix();
    
        //object.group = group;
    
        //env.mesh.add(group);
        //env.local.entities.push(object);
      }
  

      return env.local.entities[0];

    } else {
      console.log('single matrix');

      const group = env.local.group;
      const matrix = decodeTransformation(data, env);

 
  
      env.local.quaternion = new THREE.Quaternion();
      env.local.position = new THREE.Vector3();
      env.local.scale = new THREE.Vector3();    
  
      matrix.decompose(env.local.position, env.local.quaternion, env.local.scale);
  
      group.quaternion.copy( env.local.quaternion );
      group.position.copy( env.local.position );
      group.scale.copy( env.local.scale );
  
      group.updateMatrix();


      return group;
    }
    
  };  

  g3d.GeometricTransformation.destroy = (args, env) => {
    console.warn('Nothing to dispose!');
  };

  g3d.GeometricTransformation.virtual = true;

  g3d.GraphicsComplex = async (args, env) => {
    
    var copy = Object.assign({}, env);
    const options = await core._getRules(args, {...env, hold: true});
    

    const pts = (await interpretate(args[0], copy)).flat();
    const vertices = new Float32Array( pts );

    //local storage
    copy.vertices = {
      //geometry: new THREE.BufferGeometry(),
      position: new THREE.BufferAttribute( vertices, 3 ),
      colored: false,
      handlers: []
    };

    env.local.vertices = copy.vertices;

    //copy.vertices.geometry.setAttribute( 'position',  );

    if ('VertexColors' in options) {
      const colors = await interpretate(options["VertexColors"], env);
      copy.vertices.colored = true;
      copy.vertices.colors = new THREE.Float32BufferAttribute( new Float32Array( colors.flat() ), 3 );
    }

    const group = new THREE.Group();
    env.local.group = group;

    await interpretate(args[1], copy);

    env.mesh.add(group);
    //copy.geometry.dispose();
  };

  g3d.Reflectivity = () => {
    console.warn('not implemented');
  };

  g3d.GraphicsComplex.update = async (args, env) => {
    env.wake();

    const pts = (await interpretate(args[0], env)).flat();
    const vertices = new Float32Array( pts );

    env.local.vertices.position.set( vertices );
    env.local.vertices.position.needsUpdate = true;

    if (env.local.vertices.colored) {
      const options = await core._getRules(args, {...env, hold: true});
      const colors = await interpretate(options["VertexColors"], env);
      env.local.vertices.colors.set(new Float32Array( colors.flat() ));
      env.local.vertices.colors.needsUpdate = true;
    }

    for (let i=0; i<env.local.vertices.handlers.length; ++i) {
      env.local.vertices.handlers[i]();
    }
  };  

  g3d.GraphicsComplex.destroy = async (args, env) => {
    //env.local.vertices.position.dispose();
    //if (env.local.vertices.colored) env.local.vertices.colors.dispose();
  };  

  g3d.GraphicsComplex.virtual = true;


  g3d.Polygon = async (args, env) => {
    var geometry;
    let material;

    if (env.hasOwnProperty("vertices")) {
      geometry = new THREE.BufferGeometry();
      geometry.setAttribute('position', env.vertices.position);
      //env.vertices.geometry.clone();

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
        console.log('Set Index');
        geometry.setIndex( extendedIndexes.flat().map((e)=>e-1) );
        
        
      }

      //handler for future recomputations (in a case of update)
      env.vertices.handlers.push(() => {
        geometry.computeVertexNormals();
      });

      geometry.computeVertexNormals();

      //check if colored (Material BUG) !!!
      if (env.vertices.colored) {
        //geometry.setAttribute()
        geometry.setAttribute( 'color', env.vertices.colors );

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
          emissiveIntensity: env.emissiveIntensity
        });         
      }

    } else { 
      geometry = new THREE.BufferGeometry();
      let points = await interpretate(args[0], env);

      let vertices = new Float32Array( points.flat() );

  

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

      geometry.setAttribute( 'position', new THREE.BufferAttribute( vertices, 3 ) );
      geometry.computeVertexNormals();

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


    //console.log(env.opacity);
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



  g3d.Specularity = (args, env) => { };

  g3d.Text = (args, env) => { };

  g3d.Directive = async (args, env) => { 
    for (const i of args) {
      await interpretate(i, env);
    }
  };

  g3d.PlaneGeometry = () => { };


  g3d.Line = async (args, env) => {
    
    var geometry;
    //let vertices;

    if (env.hasOwnProperty("vertices")) {
      
      //vertices = env.vertices;
      geometry = new THREE.BufferGeometry();
      geometry.setAttribute('position', env.vertices.position);

      let a = await interpretate(args[0], env);
      

      geometry.setIndex( a.flat().map((e)=>e-1) );
     
      //geometry.setAttribute( 'position', new THREE.BufferAttribute( vertices, 3 ) );


    } else { 
      geometry = new THREE.BufferGeometry();
      const points = await interpretate(args[0], env);
      geometry.setAttribute( 'position', new THREE.BufferAttribute( new Float32Array(points.flat()), 3 ) );

    }


      const material = new THREE.LineBasicMaterial({
        linewidth: env.thickness,
        color: env.color,
        opacity: env.opacity,
        transparent: env.opacity < 1.0 ? true : false
      });
      const line = new THREE.Line(geometry, material);

      env.local.line = line;
 
      env.mesh.add(line);

      return line;

      //geometry.dispose();
      //material.dispose();
  };

  g3d.Line.update = async (args, env) => {
    
    const points = await interpretate(args[0], env);

    const positionAttribute = env.local.line.geometry.getAttribute( 'position' );

    positionAttribute.needsUpdate = true;

    for ( let i = 0; i < positionAttribute.count; i ++ ) {
      positionAttribute.setXYZ( i, ...(points[i]));
    }

    env.local.line.geometry.computeBoundingBox();
    env.local.line.geometry.computeBoundingSphere();

    env.wake();
  };

  g3d.Line.destroy = async (args, env) => {
    if (env.local.line) env.local.line.geometry.dispose();
  };

  g3d.Line.virtual = true;

  let GUI;

  g3d.ImageSize = () => "ImageSize";
  g3d.Background = () => "Background";
  g3d.AspectRatio = () => "AspectRatio";
  g3d.Lighting = () => "Lighting";
  g3d.Default = () => "Default";
  g3d.None = () => "None";
  g3d.Lightmap = () => "Lightmap";
  g3d.Automatic = () => "Automatic"; 

  g3d.AnimationFrameListener = async (args, env) => {
    await interpretate(args[0], env);

    const options = await core._getRules(args, {...env, hold:true});
    env.local.event = await interpretate(options.Event, env);
    
    const worker = {
      state: true,
      eval: () => {
        if (!env.local.worker.state) return;
        server.kernel.emitt(env.local.event, 'True', 'Frame');
        env.local.worker.state = false;
      }
    };

    env.local.worker = worker;  
    env.Handlers.push(worker);
  };

  g3d.AnimationFrameListener.update = async (args, env) => {
    env.local.worker.state = true;
  };

  g3d.AnimationFrameListener.destroy = async (args, env) => {
    console.warn('AnimationFrameListener does not exist anymore');
    env.local.worker.eval = () => {};
  };

  g3d.AnimationFrameListener.virtual = true;

  let Water = false;
  let Sky   = false;

  g3d.Camera = (args, env) => {
    console.warn('temporary disabled');
    return;
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

    env.global.scene.add( water );
    
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
  //const options = await core._getRules(args, {...env, hold: true});

  //console.log(options);
  //const keys = Object.keys(options);

  let position = [0, 0, 10];
  let color = 0xffffff; 
  
  if (args.length > 0) color = await interpretate(args[0], copy); 

  if (args.length > 1) {
    position = await interpretate(args[1], env);
    //position = [position[0], position[1], position[2]];
  }

  
  let intensity = 100; if (args.length > 2) intensity = await interpretate(args[2], env);
  let distance = 0; if (args.length > 3) distance = await interpretate(args[3], env);
  let decay = 2; if (args.length  > 4) decay = await interpretate(args[4], env);

  console.log(position);
  console.log(color);

  const light = new THREE.PointLight(color, intensity, distance, decay);
  light.castShadow = env.shadows;
  light.position.set(...position);
  light.shadow.bias = -0.01;
  env.local.light = light;
  env.mesh.add(light);

  return light;
};

g3d.PointLight.update = async (args, env) => {
  env.wake();
  //const options = await core._getRules(args, {...env, hold: true}); 

  if (args.length > 1) {
    let pos = await interpretate(args[1], env);
    //pos = [pos[0], pos[1], pos[2]];

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

g3d.PointLight.destroy = async (args, env) => {
  console.log('PointLight destroyed');
};

g3d.PointLight.virtual = true;

g3d.SpotLight = async (args, env) => {
  const copy = {...env};
  //const options = await core._getRules(args, {...env, hold: true});

  //console.log(options);
  //const keys = Object.keys(options);

  let color = 0xffffff; if (args.length > 0) color = await interpretate(args[0], copy);

  let position = [10, 100, 10];
  let target = [0,0,0];

  if (args.length > 1) {
    position = await interpretate(args[1], env);
    if (position.length == 2) {
      target = position[1];
      //target = [target[0], target[2], -target[1]];
      position = position[0];
    }
    //position = [position[0], position[2], -position[1]];
  }

  let angle = Math.PI/3; if (args.length > 2) angle = await interpretate(args[2], env);
  
  let intensity = 100; if (args.length > 3) intensity = await interpretate(args[3], env);
  let distance = 0; if (args.length > 4) distance = await interpretate(args[4], env);
  
  let penumbra = 0; if (args.length > 5) penumbra = await interpretate(args[5], env);
  let decay = 2; if (args.length > 6) decay = await interpretate(args[6], env);
 

  const spotLight = new THREE.SpotLight( color, intensity, distance, angle, penumbra, decay );
  spotLight.position.set(...position);
  spotLight.target.position.set(...target);

  spotLight.castShadow = env.shadows;
  spotLight.shadow.bias = -0.01;
  spotLight.shadow.mapSize.height = 1024;
  spotLight.shadow.mapSize.width = 1024;

  env.local.spotLight = spotLight;
  env.mesh.add(spotLight);
  env.mesh.add(spotLight.target);

  return spotLight;
};

g3d.SpotLight.update = async (args, env) => {
  env.wake();
  //const options = await core._getRules(args, {...env, hold: true}); 

  if (args.length > 1) {
    let position = await interpretate(args[1], env);
    if (position.length == 2) {
      let target = position[1];
      //target = [target[0], target[2], target[1]];
      position = position[0];
      //position = [position[0], position[2], -position[1]];

      if (env.Lerp) {
        if (!env.local.lerp1) {
          
          console.log('creating worker for lerp of movements..');
          const worker = {
            alpha: 0.05,
            target: new THREE.Vector3(...position),
            eval: () => {
              env.local.spotLight.position.lerp(worker.target, 0.05);
            }
          };
  
          env.local.lerp1 = worker;  
  
          env.Handlers.push(worker);
        }
  
        env.local.lerp1.target.fromArray(position);

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
        env.local.spotLight.position.set(...position);
        env.local.spotLight.target.position.set(...target);
      }
    } else {

      //position = [position[0], position[2], -position[1]];

      if (env.Lerp) {
        if (!env.local.lerp1) {
          
          console.log('creating worker for lerp of movements..');
          const worker = {
            alpha: 0.05,
            target: new THREE.Vector3(...position),
            eval: () => {
              env.local.spotLight.position.lerp(worker.target, 0.05);
            }
          };
  
          env.local.lerp1 = worker;  
  
          env.Handlers.push(worker);
        }
  
        env.local.lerp1.target.fromArray(position);

              
      } else {
        env.local.spotLight.position.set(...position);
      }
    }
    


  }

};

g3d.SpotLight.destroy = async (args, env) => {
  console.log('SpotLight destoyed');
};

g3d.SpotLight.virtual = true;

g3d.Shadows = async (args, env) => {
  env.shadows = await interpretate(args[0], env);
};



g3d.HemisphereLight = async (args, env) => {
  const copy = {...env};

  if (args.length > 0) await interpretate(args[0], copy); else copy.color = 0xffffbb;
  const skyColor = copy.color;

  if (args.length > 1) await interpretate(args[1], copy); else copy.color = 0x080820;
  const groundColor = copy.color;
 if (args.length > 2) await interpretate(args[1], env);

  const hemiLight = new THREE.HemisphereLight( skyColor, groundColor, 2 );
  env.global.scene.add( hemiLight );
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

    let object = await interpretate(args[0], env);
    if (Array.isArray(object)) object = object[0];

    if (!TransformControls) TransformControls = (await import('./TransformControls-2614e9c1.js')).TransformControls;
    rules.forEach((rule)=>{
      g3d.EventListener[rule.lhs](rule.rhs, object, copy);
    });

    return null;
};

  g3d.EventListener.transform = (uid, object, env) => {
    console.log(env);
    console.warn('Controls transform is enabled');
    const control = new TransformControls(env.camera, env.global.renderer.domElement);

    const orbit = env.controlObject.o;

    control.attach(object); 

    env.global.scene.add(control); 

    const updateData = throttle((x,y,z) => {
      server.kernel.emitt(uid, `<|"position"->{${x.toFixed(4)}, ${y.toFixed(4)}, ${z.toFixed(4)}}|>`, 'transform');
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
  OrbitControls = (await import('./OrbitControls-36d93d15.js')).OrbitControls;
  GUI           = (await import('./dat.gui.module-0f47b92e.js')).GUI;  
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

  const defaultMatrix = new THREE.Matrix4().set(
    1, 0, 0, 0,//
    0, 1, 0, 0,//
    0, 0, 1, 0,//
    0, 0, 0, 1);


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
    multipleImportanceSampling: false,
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

  if (!PathRendering) params.resolutionScale = 0.5;

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
  let perspectiveCamera, orthoCamera;
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
	orthoCamera = new THREE.OrthographicCamera( orthoWidth / - 2, orthoWidth / 2, orthoHeight / 2, orthoHeight / - 2, 0, 200 );
	orthoCamera.position.set( - 40, 20, 30 );

  activeCamera = orthoCamera;

  scene = new THREE.Scene();

  if (PathRendering) {
	  //equirectCamera = new RTX.EquirectCamera();
	  //equirectCamera.position.set( - 4, 2, 3 );

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
          scene.add( controlObject.o.getObject() );
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

	

  const group = new THREE.Group();

  const allowLerp = false;
  if (options.TransitionType) {
    const type = await interpretate(options.TransitionType, env);
    if (type === 'Linear') allowLerp = true;
  }

  const envcopy = {
    ...env,
    context: g3d,
    numerical: true,
    tostring: false,
    matrix: defaultMatrix,
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
    Lerp: allowLerp,
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

  env.global.renderer = renderer;
  env.global.scene    = scene;
  envcopy.camera   = activeCamera;
  env.local.element  = container;

  if (PathRendering)
    envcopy.PathRendering = true;

  await interpretate(args[0], envcopy);

  var bbox = new THREE.Box3().setFromObject(group);
  //console.error(bbox);
  group.position.set(-(bbox.min.x + bbox.max.x) / 2, -(bbox.min.y + bbox.max.y) / 2, -(bbox.min.z + bbox.max.z) / 2);
  //throw 'fuk';
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
      await interpretate(['Line', ['JSObject', l]], {...envcopy});
    }  }
  

  group.applyMatrix4(new THREE.Matrix4().set( 
    1, 0, 0, 0,
    0, 0, 1, 0,
    0, -1, 0, 0,
    0, 0, 0, 1));

  if ('BoxRatios' in options) {
    let size = [bbox.max.x - bbox.min.x, bbox.max.z - bbox.min.z, bbox.max.y - bbox.min.y];
    let max = Math.max(...size);
    const reciprocal = size.map((e) => 1.0/(e/max));

    console.warn('Rescaling....');

    let ratios = await interpretate(options.BoxRatios, env);
    ratios = [ratios[0], ratios[2], ratios[1]];

    max = Math.max(...ratios);
    ratios = ratios.map((e, index) => reciprocal[index] * e/max);

    group.applyMatrix4(new THREE.Matrix4().makeScale(...ratios));
  }  

  group.position.add(new THREE.Vector3(0,1,0));
  //recalculate
  bbox = new THREE.Box3().setFromObject(group);
  //console.log(bbox);

  if (envcopy.camera.isOrthographicCamera) {
    console.warn('fitting camera...');
    const camera = envcopy.camera;
    camera.zoom = Math.min(orthoWidth / (bbox.max.x - bbox.min.x),
    orthoHeight / (bbox.max.y - bbox.min.y)) * 0.8;
    camera.updateProjectionMatrix();
  }
    
 
  scene.add(group);
  //console.error(new THREE.Box3().setFromObject(scene));

  scene.updateMatrixWorld();

  //console.error(new THREE.Box3().setFromObject(scene));

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
